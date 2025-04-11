import pandas as pd
import sys

def syry_inv_to_bloc(df_pathfile, outfile_name):
    df = pd.read_csv(df_pathfile, sep='\t')
    
    # Liste pour stocker les résultats
    merged_intervals_all = []
    
    # Traiter chaque paire iso-chromosome séparément
    for iso in df['iso'].unique():
        df_iso = df[df['iso'] == iso]
        
        for chromosome in df_iso['chr'].unique():
            # Filtrer pour l'iso et le chromosome courants
            df_iso_chr = df_iso[df_iso['chr'] == chromosome].sort_values('start')
            merged_intervals_chr = []

            for _, row in df_iso_chr.iterrows():
                if not merged_intervals_chr or row['start'] > merged_intervals_chr[-1]['end']:
                    merged_intervals_chr.append({
                        'iso': iso,
                        'chr': chromosome,
                        'start': row['start'], 
                        'end': row['end']})
                else: 
                    merged_intervals_chr[-1]['end'] = max(merged_intervals_chr[-1]['end'], row['end'])
            
            # Ajouter les intervalles fusionnés de cette paire iso-chromosome à la liste globale
            merged_intervals_all.extend(merged_intervals_chr)
    
    result = pd.DataFrame(merged_intervals_all)

    result = result.sort_values(by=['iso', 'chr', 'start']).reset_index(drop=True)
    result.to_csv(outfile_name, sep="\t", index=False)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(f"Usage: python3 {sys.argv[0]} <input_file> <output_file>")
        sys.exit(1)

    df_pathfile = sys.argv[1]
    out_filname = sys.argv[2]
    syry_inv_to_bloc(df_pathfile,out_filname)