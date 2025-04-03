import os, re
from pathlib import Path
import pandas as pd
import warnings
warnings.simplefilter(action='ignore', category=pd.errors.ParserWarning)

def process_chromosome_name(chrom_name):
    try:
        chrom_num = int(chrom_name.replace('chr', ''))
        if chrom_num in [3, 16, 17]:
            return f'chr{chrom_num + 2000}'
        else:
            return f'chr{chrom_num + 1000}'
    except ValueError:
        return chrom_name


def create_df_all_TE(directory_path,corespondance_pathfile=None,outfile_name=None):
    all_dfs = []
    # Parcourir tous les fichiers .out dans le dossier
    for file_path in Path(directory_path).glob('*.out'):
        df = pd.read_csv(file_path, sep=r'\s+', engine='python', header=None, 
                        index_col=False, skip_blank_lines=True, skiprows=2)
        df.columns = ['score', 'div.', 'del.', 'ins.', 'sequence', 'begin', 'end', '(left)', 'strand', 
            'repeat', 'class/family', 'begin.1', 'end.1', '(left).1', 'ID']
        split = os.path.basename(file_path).split('_')
        isolat = split[0] + '_' + split[1]
        df['Isolat'] = isolat

        df['sequence'] = df['sequence'].apply(process_chromosome_name)
        all_dfs.append(df)



    df = pd.concat(all_dfs, ignore_index=True)

    if corespondance_pathfile:
        dico_corespondance = {}
        with open(corespondance_pathfile, "r") as infos_file:
                next(infos_file)
                for line in infos_file: dico_corespondance[line.split()[1]] = line.split()[3]
        df["Lineage"] = df[["Isolat"]].map(dico_corespondance.get)
    print(df)
    # Affichage des valeurs uniques de la colonne 'class/family'


    # POUR LAUREU
    df = df[(df['class/family'] != 'Low_complexity') & (df['class/family'] != 'Simple_repeat')]

    print("Valeurs uniques de class/family:")
    print(df['class/family'].unique())

    df['class/family'] = df['class/family'].replace(['Simple_repeat', 'Low_complexity', 'ClassI/Unclassified','RC/Helitron'], 'Others')

    df_filtre = df
    
    df_filtre['class/family'] = df_filtre['class/family'].replace({r'^LTR.*$': 'LTR', r'^DNA.*$': 'DNA', r'^LINE.*$': 'LINE'}, regex=True)

    # Calculer le pourcentage de pertes
    nb_lignes_initial, nb_lignes_apres = len(df), len(df_filtre)
    pourcentage_pertes = ((nb_lignes_initial - nb_lignes_apres) / nb_lignes_initial) * 100
    print(f"Pourcentage de pertes après filtrage sur TE: {pourcentage_pertes:.2f}%")

    # Exporter le DataFrame filtré dans un fichier CSV
    if outfile_name : df_filtre.to_csv(outfile_name, index=False, sep =';')

    return df_filtre



def create_df_TEbyIsolat(input_file, output_tsv):
    # Lire le fichier avec les bonnes colonnes
    print(f"Lecture du fichier {input_file}...")

    # Définir les noms de colonnes correspondant à votre format
    columns = ["score", "div", "del", "ins", "sequence", "begin", "end", 
               "left", "strand", "repeat", "class_family", "begin_1", "end_1", 
               "left_1", "ID", "Isolat", "Lineage"]

    df = pd.read_csv(input_file, sep=";", names=columns, low_memory=False)

    # Convertir 'begin' et 'end' en numérique, en forçant le type
    df['begin'] = pd.to_numeric(df['begin'], errors='coerce')
    df['end'] = pd.to_numeric(df['end'], errors='coerce')

    # Supprimer les lignes où begin ou end n'ont pas pu être convertis en nombres
    df = df.dropna(subset=['begin', 'end'])

    # Calculer la taille de chaque élément TE
    df['size'] = df['end'] - df['begin'] + 1

    # Extraire les informations de lignée pour chaque isolat
    print("Extraction des informations de lignée...")
    lineage_df = df[['Isolat', 'Lineage']].drop_duplicates()
    print(lineage_df)
    # S'assurer qu'un isolat n'a qu'une seule lignée
    if lineage_df.duplicated(subset=['Isolat']).any():
        print("Attention : Certains isolats ont plusieurs lignées. Vérifiez les données !")

    # Agrégation des données par isolat et classe/famille
    print("Agrégation des données par isolat et classe/famille...")
    te_sizes = df.groupby(['Isolat', 'class_family'])['size'].sum().reset_index()

    # Création du tableau pivot
    print("Création de la matrice pivot...")
    pivot_df = te_sizes.pivot(index='Isolat', columns='class_family', values='size').fillna(0).astype(int)

    # Ajouter la colonne Lineage
    pivot_df = pivot_df.merge(lineage_df, on='Isolat', how='left')

    # Enregistrement du fichier de sortie
    print(f"Enregistrement des résultats dans {output_tsv}...")
    pivot_df.to_csv(output_tsv, sep='\t', index=False)

    return pivot_df



if __name__ == "__main__":
    directory_path='/home/hugob/work/05_TE/data_annotation'
    corespondance_pathfile ='/home/hugob/work/00_assemblages/Isolates.Master.Info.txt'
    outfile_name = '/home/hugob/work/05_TE/all_TE2.csv'
    create_df_all_TE(directory_path,corespondance_pathfile,outfile_name)

    input_file = '/home/hugob/work/05_TE/all_TE2.csv'
    outfile_name = '/home/hugob/work/05_TE/TEbyIsolats2.tsv'
    result_df = create_df_TEbyIsolat(input_file, outfile_name)
    print(f"Résultats enregistrés dans {outfile_name}")
    print(f"Nombre d'isolats: {len(result_df.index)}")
    print(f"Nombre de classes/familles de TE: {len(result_df.columns) - 1}")  # -1 pour exclure Lineage
