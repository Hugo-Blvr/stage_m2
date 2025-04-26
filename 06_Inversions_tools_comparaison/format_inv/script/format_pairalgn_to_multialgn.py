import pandas as pd
import sys

def asso_lineage(df_pathfile,lineage_pathfile):
    df = pd.read_csv(df_pathfile, sep='\t')
    lineage_df = pd.read_csv(lineage_pathfile, sep='\t', usecols=[1, 3])
    lineage_df.columns = ['iso', 'lineage']
    lineage_df = lineage_df.drop_duplicates(subset='iso')
    # Merge pour récupérer le lineage de 'iso'
    df = df.merge(lineage_df, on='iso', how='left')
    df = df.rename(columns={'lineage': 'iso_lineage'})

    # Pour 'mapon', on renomme temporairement les colonnes pour éviter conflit
    lineage_df = lineage_df.rename(columns={'iso': 'mapon', 'lineage': 'mapon_lineage'})
    df = df.merge(lineage_df, on='mapon', how='left')

    return df


def pairwise_algn_to_multi(df, outfile_name):
    merged_intervals_all = []

    for iso in df['iso'].unique():
        df_iso = df[df['iso'] == iso]
        iso_lineage = df_iso['iso_lineage'].iloc[0]

        for chromosome in df_iso['chr'].unique():
            df_iso_chr = df_iso[df_iso['chr'] == chromosome].sort_values('start')

            for fusion_type in ['all', 'intra', 'inter']:
                if fusion_type == 'all':
                    subset = df_iso_chr
                elif fusion_type == 'intra':
                    subset = df_iso_chr[df_iso_chr['mapon_lineage'] == iso_lineage]
                elif fusion_type == 'inter':
                    subset = df_iso_chr[df_iso_chr['mapon_lineage'] != iso_lineage]

                merged_intervals_chr = []

                for _, row in subset.iterrows():
                    if not merged_intervals_chr or row['start'] > merged_intervals_chr[-1]['end']:
                        merged_intervals_chr.append({
                            'iso': iso,
                            'chr': chromosome,
                            'start': row['start'],
                            'end': row['end'],
                            'mapon': fusion_type  # ici on met la catégorie comme mapon
                        })
                    else:
                        merged_intervals_chr[-1]['end'] = max(merged_intervals_chr[-1]['end'], row['end'])

                merged_intervals_all.extend(merged_intervals_chr)

    result = pd.DataFrame(merged_intervals_all)
    result = result.sort_values(by=['iso', 'chr','mapon','start']).reset_index(drop=True)
    result.to_csv(outfile_name, sep="\t", index=False)



if __name__ == "__main__":
    
    
    if len(sys.argv) < 2:
        print(f"Usage: python3 {sys.argv[0]} <input_file> <output_file>")
        sys.exit(1)

    df_pathfile = sys.argv[1]
    out_filname = sys.argv[2]

    lineage_file = "/home/hugob/stage_m2/00_Data/Isolates.Master.Info.txt"
    df = asso_lineage(df_pathfile, lineage_file)
    pairwise_algn_to_multi(df,out_filname)