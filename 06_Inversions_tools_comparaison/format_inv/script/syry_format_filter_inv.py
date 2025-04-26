import pandas as pd
import numpy as np
import warnings
import sys
warnings.simplefilter("ignore", UserWarning)

def format_df_myinv(df_pathfile):
    df = pd.read_csv(df_pathfile, sep='\t')
    df = df.sort_values(by=['target','query','t_chr','t_start']).reset_index(drop=True)
    df['ID'] = df.index
    # ==================== PROPRE A NOS DONNEES================================
    split_target = df['target'].str.split('_')
    df['target'] = split_target.str[0] + '_' + split_target.str[1]
    split_target = df['query'].str.split('_')
    df['query'] = split_target.str[0] + '_' + split_target.str[1]
    # ==========================================================================
    df = df.sort_values(by=['target','query','t_chr','t_start']).reset_index(drop=True)

    #out_filname = "~/stage_m2/06_Inversions_tools_comparaison/format_inv/output/check_minimap_pain.tsv"
    #df.to_csv(out_filname, sep="\t", index=False)

    return df

def are_same_position(start1, end1, start2, end2, min_overlap):
    length1 = end1 - start1
    length2 = end2 - start2
    
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    
    if overlap_start < overlap_end:
        overlap_length = overlap_end - overlap_start
        overlap_percent1 = overlap_length / length1
        overlap_percent2 = overlap_length / length2
    
        if overlap_percent1 >= min_overlap and overlap_percent2 >= min_overlap:
            return True
    
    return False

def filter_englobing_inversions(df):
    """
    Filtre les inversions qui englobent exactement plusieurs autres inversions
    """
    # Créer une copie du DataFrame pour ne pas modifier l'original
    result_df = df.copy()
    
    # Ajouter une colonne pour faciliter le regroupement
    result_df['query_target_pair'] = result_df['query'] + '_' + result_df['target']
    
    # Listes pour stocker les résultats
    indices_to_keep = []
    indices_to_remove = []
    
    # Regrouper par paires query-target et chromosome
    for (qt_pair, chr_name), group in result_df.groupby(['query_target_pair', 't_chr']):
        if len(group) <= 1:  # Ignorer les groupes avec une seule inversion
            indices_to_keep.extend(group.index.tolist())
            continue
            
        # Trier le groupe par start et end
        group = group.sort_values(by=['t_start', 't_end', 'q_start', 'q_end'])
        
        # Pour chaque inversion, vérifier si elle englobe exactement d'autres inversions
        for i, row in group.iterrows():
            # Vérifier les inversions qui pourraient être englobées par celle-ci
            potential_subinversions = group[
                (group.index != i) & 
                (group['t_start'] >= row['t_start']) & 
                (group['t_end'] <= row['t_end']) &
                (group['q_start'] >= row['q_start']) & 
                (group['q_end'] <= row['q_end'])
            ]
            
            # Si cette inversion englobe d'autres inversions et couvre exactement 
            # le même intervalle que l'ensemble de ces inversions
            if len(potential_subinversions) > 0:
                # Vérifier si l'inversion englobe exactement les sous-inversions
                sub_t_start = potential_subinversions['t_start'].min()
                sub_t_end = potential_subinversions['t_end'].max()
                sub_q_start = potential_subinversions['q_start'].min()
                sub_q_end = potential_subinversions['q_end'].max()
                
                # Si l'inversion englobe exactement les sous-inversions (ou avec une marge très faible)
                marge = 0
                if (abs(row['t_start'] - sub_t_start) <= marge and 
                    abs(row['t_end'] - sub_t_end) <= marge and
                    abs(row['q_start'] - sub_q_start) <= marge and 
                    abs(row['q_end'] - sub_q_end) <= marge):
                    
                    # Cette inversion est artificielle, on la supprime
                    indices_to_remove.append(i)
                else:
                    # L'inversion ne correspond pas exactement à la combinaison de sous-inversions, on la garde
                    indices_to_keep.append(i)
            else:
                # Aucune sous-inversion trouvée, on garde l'inversion
                indices_to_keep.append(i)
    
    print(f"\nNombre d'inversions englobantes supprimées: {len(indices_to_remove)}")
    
    # Retourner le DataFrame filtré
    filtered_df = result_df.loc[indices_to_keep].drop(columns=['query_target_pair'])
    return filtered_df


def handle_duplicate_inversions(df, min_overlap):
    """
    Gère les doublons d'inversions (même paire d'individus avec query/target inversés)
    et garde une version fusionnée avec min-start et max-end pour chaque individu
    """
    data = df.sort_values(['target','query','t_chr','t_start']).copy()
    parsed = set()
    ID_fusion = []
    for (target, query, chr_name), group in data.groupby(['target', 'query', 't_chr']):
        
        id_group = f"{min(target, query)}_{max(target, query)}_{chr_name}"
        if id_group in parsed: continue
        parsed.add(id_group)

        group = group.sort_values('t_start')
        # groupe inverse : même chr, mais target <-> query inversés
        mask_inverse_group = ((data['target'] == query) & (data['query'] == target) & (data['t_chr'] == chr_name))
        inverse_group = data[mask_inverse_group].sort_values('q_start')

        j = 0  # index glissant pour sauter les j première ligne d'inverse_group
        for i, row1 in group.iterrows():
            # avancer j jusqu'à ce que row2['q_end'] soit plus grand que row1['t_start']
            while j < len(inverse_group) and row1['t_start'] >= inverse_group.iloc[j]['q_end']: j += 1

            max_start = row1['t_start'] + ((row1['t_end'] - row1['t_start']) * (1 - min_overlap))
            
            for k in range(j, len(inverse_group)):
                row2 = inverse_group.iloc[k]

                if row2['q_start'] > max_start: break  # plus aucune row2 ne peut matcher

                same_t = are_same_position(row1['t_start'], row1['t_end'], row2['q_start'], row2['q_end'], min_overlap)
                same_q = are_same_position(row1['q_start'], row1['q_end'], row2['t_start'], row2['t_end'], min_overlap)
                if same_t and same_q:                    
                    if row1['ID'] in ID_fusion : print(f'❌ ERREUR : Une inversion a deux id corespondant : {row1["ID"]}') 
                    ID_fusion.append(row1['ID'])
                    
                    mask_update = (data['ID'] == row1['ID'])
                    new_tstart, new_qstart = min(row1['t_start'],row2['q_start']), min(row1['q_start'],row2['t_start'])
                    new_tend, new_qend = max(row1['t_end'],row2['q_end']), max(row1['q_end'],row2['t_end'])
                    data.loc[mask_update, ['t_start', 't_end', 'q_start','q_end']] = [new_tstart, new_tend, new_qstart, new_qend]
                    
                    mask_delete = (data['ID'] == row2['ID'])
                    data.loc[mask_delete, 'ID'] = -1

                    if mask_update.sum() > 1 or mask_delete.sum() > 1: print(f"❌ ERREUR : ID doublé dans le data")

    print(f"\nSupression de {(data['ID'] == -1).sum()} inversions multiple")
    data = data[data['ID'] != -1].reset_index(drop=True)

    data = data.reset_index(drop=True)
    return data.reset_index(drop=True)



    
def final_format(df, min_overlap):

    df['ID'] = df.index
    df_T = df[['target', 't_chr', 't_start', 't_end','query', 'ID']].copy()
    df_Q = df[['query', 'q_chr', 'q_start', 'q_end','target', 'ID']].copy()
    colnames = ['iso', 'chr', 'start', 'end','mapon','ID']
    df_T.columns, df_Q.columns = colnames, colnames
    df_T['type'], df_Q['type'] = 'T', 'Q'   
    format_data = pd.concat([df_T, df_Q], ignore_index=True)
    format_data['size'] = format_data['end'] - format_data['start']   
    format_data = format_data.sort_values(by=['iso','chr','mapon','start'], ignore_index=True)

    for (iso, mapon, chr_name), group in format_data.groupby(['iso', 'mapon', 'chr']):
        if len(group) > 1:
            sorted_group = group.sort_values('start').reset_index(drop=True)
            clusters = []
            processed = set()
            for idx, row1 in sorted_group.iterrows():
                if idx in processed: continue

                cluster = [row1]
                processed.add(idx)
                max_start = row1['start'] + ((row1['end'] - row1['start']) * (1 - min_overlap))
                for next_idx, row2 in sorted_group[idx+1:].iterrows():
                    if row2['start'] > max_start: break
                    if next_idx not in processed and are_same_position(row1['start'], row1['end'], row2['start'], row2['end'], min_overlap):
                        cluster.append(row2)
                        are_same_position(row1['start'], row1['end'], row2['start'], row2['end'], min_overlap)
                        processed.add(next_idx)

                if len(cluster) > 1 : clusters.append(cluster)

            for cluster in clusters:
                reference_id = cluster[0]['ID']

                # Récupère les starts, ends, et types
                new_start = min([row['start'] for row in cluster])
                new_end = max([row['end'] for row in cluster])
                new_size = new_end - new_start
                types_in_cluster = set(row['type'] for row in cluster)
                if len(types_in_cluster) == 1: new_type = next(iter(types_in_cluster))
                else: new_type = types_in_cluster

                # Mettre à jour la ligne de référence
                mask_ref = ((format_data['ID'] == reference_id) & (format_data['iso'] == iso) &
                        (format_data['mapon'] == mapon) & (format_data['chr'] == chr_name))
                format_data.loc[mask_ref, ['start', 'end', 'size', 'type']] = [new_start, new_end, new_size, new_type]

                for row in cluster[1:]:
                    old_id = row['ID']
    
                    # Marquer pour suppression dans le groupe actuel
                    mask_delete = ((format_data['ID'] == old_id) & (format_data['iso'] == iso) &
                                (format_data['mapon'] == mapon) & (format_data['chr'] == chr_name))
                    format_data.loc[mask_delete, 'ID'] = -1

                    # Trouve les lignes complémentaires (iso <-> mapon inversés)
                    # Mettre à jour l'ID dans le groupe inverse
                    mask_update = ((format_data['ID'] == old_id) & (format_data['iso'] == mapon) &
                                (format_data['mapon'] == iso) & (format_data['chr'] == chr_name))
                    format_data.loc[mask_update, 'ID'] = reference_id


    inv_sup = (format_data['ID'] == -1).sum()
    format_data = format_data[format_data['ID'] != -1].reset_index(drop=True)

    print(f'\nSupression de {inv_sup} inversions multi-map')
    format_data = format_data.sort_values(by=['iso','chr','mapon','start']).reset_index(drop=True)
    return format_data.reset_index(drop=True)




def format_and_filter_inversions(df_pathfile, outfile_name, min_overlap):
    df = format_df_myinv(df_pathfile)
    
    # ÉTAPE 1: Filtrer les inversions englobantes
    englob_filtered_df = filter_englobing_inversions(df)
    
    # ÉTAPE 2: Gérer les doublons (même individus en query/target inversés)
    duplicate_filterd_df = handle_duplicate_inversions(englob_filtered_df, min_overlap)

    # ÉTAPE 3: Format inv
    final_df = final_format(duplicate_filterd_df,min_overlap)

    
    # Trier et sauvegarder le résultat
    final_df = final_df.sort_values(by=['iso','chr','mapon','start']).reset_index(drop=True)
    final_df.to_csv(outfile_name, sep="\t", index=False)
    return final_df




if __name__ == "__main__":

    if len(sys.argv) < 3:
        print(f"Usage: python3 {sys.argv[0]} <input_file> <output_file>")
        sys.exit(1)
    

    df_pathfile = sys.argv[1]
    out_filname = sys.argv[2]
    min_overlap = 0.85
    format_and_filter_inversions(df_pathfile, out_filname, min_overlap)


    """    
    df_pathfile = "~/stage_m2/05_Inversion_identification/output/inv_minimap.tsv"
    out_filname = "~/stage_m2/06_Inversions_tools_comparaison/format_inv/output/test_toto_minimap.tsv"
    min_overlap = 0.85
    format_and_filter_inversions(df_pathfile, out_filname, min_overlap)
    """







