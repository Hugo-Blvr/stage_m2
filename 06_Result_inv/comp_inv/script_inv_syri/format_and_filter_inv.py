import pandas as pd
import warnings
import sys
warnings.simplefilter("ignore", UserWarning)


############################## CREATE DATAFRAME ##############################

def format_df_myinv(df_pathfile):
    df = pd.read_csv(df_pathfile, sep='\t')


    # ==================== PRORPRE A NOS DONNEES================================
    split_target = df['target'].str.split('_')
    df['target'] = split_target.str[0] + split_target.str[1]
    split_target = df['query'].str.split('_')
    df['query'] = split_target.str[0] + split_target.str[1]
    # ==========================================================================

    df_T = df[['target', 't_chr', 't_start', 't_end','query']].copy()
    df_Q = df[['query', 'q_chr', 'q_start', 'q_end','target']].copy()
    colnames = ['iso', 'chr', 'start', 'end','mapon']
    df_T.columns, df_Q.columns = colnames, colnames
    df_T['type'], df_Q['type'] = 'T', 'Q'   

    format_data = pd.concat([df_T, df_Q], ignore_index=True)
    format_data['size'] = format_data['end'] - format_data['start']    
    format_data = format_data.sort_values(by=['iso','mapon','chr','start']).reset_index(drop=True)

    return format_data
    

def filter_inversions(df_pathfile, outfile_name, filter = True, del_recover = True):
    df = format_df_myinv(df_pathfile)
    if not filter : 
        return df.to_csv(outfile_name, sep="\t", index=False)
    
    # Formatage et séparation des données
    df_t, df_q = df[df['type'] == 'T'].copy(), df[df['type'] == 'Q'].copy()
    df_t, df_q = df_t.sort_values(by=['iso','chr','mapon']), df_q.sort_values(by=['iso','chr','mapon'])

    def find_matching_inversions(t_row):
        """Trouve les inversions Q correspondant à une inversion T"""
        inversion_size = t_row['end'] - t_row['start']
        margin = 0.2 * inversion_size
        matching_mask = (
            (df_q['chr'] == t_row['chr']) &
            (df_q['iso'] == t_row['iso']) &
            (df_q['mapon'] == t_row['mapon']) &
            (df_q['start'].between(t_row['start'] - margin, t_row['start'] + margin)) &
            (df_q['end'].between(t_row['end'] - margin, t_row['end'] + margin))
        )
        return df_q[matching_mask]

    def create_combined_inversion(t_row, q_row):
        """Crée une nouvelle inversion combinée T-Q"""
        new_row = t_row.copy()
        new_row['type'] = 'T-Q'
        new_row['start'] = min(t_row['start'], q_row['start'])
        new_row['end'] = max(t_row['end'], q_row['end'])
        new_row['size'] = new_row['end'] - new_row['start']
        return new_row

    def filter_overlapping_inversions(group):
        """Filtre les inversions avec des chevauchements"""
        # Trier par taille décroissante
        group = group.sort_values('size', ascending=False)
        
        # Trouver les indices des inversions à supprimer
        to_remove = set()
        for i, row in group.iterrows():
            overlaps = group[ (group['start'] < row['end']) & (group['end'] > row['start'])]
            if len(overlaps) > 2: to_remove.add(i)
        
        # Supprimer les inversions redondantes
        group = group.drop(index=list(to_remove))
        
        # Conserver la plus grande inversion en cas de chevauchement
        final_group = []
        for row in group.itertuples(index=False):
            overlaps = [existing for existing in final_group if existing['start'] < row.end and existing['end'] > row.start]
        
            if not overlaps or row.size > max(overlap['size'] for overlap in overlaps):
                # Supprimer l'existant si la nouvelle est plus grande
                final_group = [existing for existing in final_group if existing['size'] >= row.size]
                final_group.append(row._asdict())
        
        return pd.DataFrame(final_group)

    # Processus principal
    result_rows = []
    
    for _, t_row in df_t.iterrows():
        matching_q_rows = find_matching_inversions(t_row)
        
        if not matching_q_rows.empty:
            # Sélectionner l'inversion T la plus proche
            best_q_row = matching_q_rows.iloc[(matching_q_rows['start'] - t_row['start']).abs().argmin()]
            combined_inversion = create_combined_inversion(t_row, best_q_row)
            result_rows.append(combined_inversion)

    # Gestion du cas où aucune inversion n'est trouvée
    if not result_rows: return pd.DataFrame()

    # Création du DataFrame final
    df_filter = pd.DataFrame(result_rows)
    
    # Filtrage par groupe
    if del_recover : 
        df_filter = df_filter.groupby(['iso', 'chr', 'mapon'], group_keys=False)[df_filter.columns].apply(filter_overlapping_inversions)
    
    df_filter = df_filter.drop(columns=['type'])
    df_filter = df_filter.drop(columns=['size'])
    df_filter = df_filter.sort_values(by=['iso','mapon','chr','start']).reset_index(drop=True)
    df_filter.to_csv(outfile_name, sep="\t", index=False)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 format_inv.py <in_pathname> <out_pathname>")
        sys.exit(1)

    df_pathfile = sys.argv[1]
    out_filname = sys.argv[2]
    if len(sys.argv) > 3:  filter = False
    else: filter = True
    
    filter_inversions(df_pathfile,out_filname, filter)