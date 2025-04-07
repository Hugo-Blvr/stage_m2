#!/usr/bin/env python3
"""
Script pour traiter des fichiers funannotate .gff et les combiner en un seul fichier TSV.
Ce script compile les données, les filtre et les formate.
"""

import os
import sys
import argparse
from pathlib import Path
import pandas as pd
import warnings
# Supprimer les avertissements de pandas pour le parseur
warnings.simplefilter(action='ignore', category=pd.errors.ParserWarning)

def process_chromosome_name(chrom_name):
    """
    Traite les noms de chromosomes selon une règle spécifique.
    Args:
        chrom_name (str): Nom du chromosome à traiter (ex: 'chr3')
    Returns:
        str: Nom du chromosome transformé
    """
    try:
        chrom_num = int(chrom_name.replace('chr', ''))
        if chrom_num in [3, 16, 17]: return f'chr{chrom_num + 2000}'
        else: return f'chr{chrom_num + 1000}'
    except ValueError: return chrom_name


def create_df_all_gene(directory_path, outfile_path):
    """
    Crée un DataFrame unique à partir de tous les fichiers .gff dans un répertoire.
    Args:
        directory_path (str): Chemin vers le répertoire contenant les fichiers .gff
        outfile_path (str): Chemin pour le fichier de sortie TSV
        verbose (bool): Afficher des informations détaillées pendant le traitement
    Returns:
        bool: True si le traitement a réussi, False sinon
    """
    all_dfs = []
    file_count = 0
    
    # Vérifier que le répertoire existe
    if not os.path.isdir(directory_path):
        print(f"ERROR: Le répertoire '{directory_path}' n'existe pas", file=sys.stderr)
        return False
    
    # Créer le répertoire de sortie s'il n'existe pas
    outfile_dir = os.path.dirname(outfile_path)
    if outfile_dir and not os.path.exists(outfile_dir):
        try:
            os.makedirs(outfile_dir)
            print(f"Répertoire de sortie créé: {outfile_dir}")
        except OSError as e:
            print(f"Erreur lors de la création du répertoire de sortie: {e}", file=sys.stderr)
            return False
    
    # Parcourir tous les fichiers .gff dans le dossier
    out_files = list(Path(directory_path).glob('*.gff*'))
    if not out_files:
        print(f"Attention: Aucun fichier .gff trouvé dans '{directory_path}'", file=sys.stderr)
        return False
    
    print(f"Traitement de {len(out_files)} fichiers .gff...")
    
    for file_path in out_files:
        try:
            print(f"Traitement du fichier: {file_path.name}")   
            df = pd.read_csv(file_path, sep=r'\t',comment='##', header=None)
            
            # Vérifier que le fichier contient des données
            if df.empty:
                print(f"Attention: Le fichier {file_path.name} est vide ou mal formaté", file=sys.stderr)
                continue
            
            # Attribuer les noms de colonnes originaux
            df.columns = ['chr', 'tools', 'type', 'start', 'end', 'score', 'strand', 'phase', 'atr']
            df = df[df['type'] != 'CDS']  # supprimer les lignes CDS car == exon
            # Reformater le DataFrame
            df['iso'] = file_path.stem 
            df['chr'] = df['chr'].apply(process_chromosome_name)
            df = df[['iso', 'chr', 'start', 'end', 'type', 'atr']]

            # 1. Extraire `oid` de l'attribut
            df['oid'] = df['atr'].str.extract(r'ID=([^;]+)')[0].str.split('-').str[0]

            # 2. Extraire le `product` uniquement pour les RNA
            df['product'] = '.' 
            mask = df['type'].str.endswith('RNA', na=False)
            df.loc[mask, 'product'] = df.loc[mask, 'atr'].str.extract(r'product=([^;]+)')[0]
            df.drop(columns=['atr'], inplace=True)

            # 3. Récupérer les produits au niveau de chaque oid (si RNA présent)
            product_map = df[mask][['oid', 'product']].dropna().drop_duplicates()
            df = df[~mask] # supprimer les lignes RNA car == gene
            df = df.merge(product_map, on='oid', how='left', suffixes=('', '_from_rna'))
            df['product'] = df['product_from_rna'].combine_first(df['product'])
            df.drop(columns=['product_from_rna'], inplace=True)

            # 4. Générer les introns par isolement de régions non couvertes par les exons
            def generate_introns(subdf):
                gene = subdf[subdf['type'] == 'gene']     
                if len(gene) != 1:
                    print(f"ALERTE: {len(gene)} gènes trouvés pour l'oid {subdf['oid'].iloc[0]}")
                    return pd.DataFrame()
                gene = gene.iloc[0]

                exons = subdf[subdf['type'] == 'exon'].sort_values(by='start')
                introns = []
                if len(exons) == 0: return pd.DataFrame()
                
                # Vérifier s'il y a un espace entre le début du gène et le premier exon
                first_exon = exons.iloc[0]
                if gene['start'] < first_exon['start']:
                    introns.append({'iso': gene['iso'], 'chr': gene['chr'],
                        'start': gene['start'], 'end': first_exon['start'] - 1,
                        'type': 'intron', 'oid': gene['oid'], 'product': gene['product']})
                
                # Créer des introns entre les exons consécutifs
                for i in range(1, len(exons)):
                    prev_exon = exons.iloc[i-1]
                    curr_exon = exons.iloc[i]
                    
                    intron_start = prev_exon['end'] + 1
                    intron_end = curr_exon['start'] - 1
                    
                    if intron_start <= intron_end:
                        introns.append({'iso': gene['iso'], 'chr': gene['chr'],
                            'start': intron_start, 'end': intron_end,
                            'type': 'intron', 'oid': gene['oid'], 'product': gene['product']})
                
                # Vérifier s'il y a un espace entre le dernier exon et la fin du gène
                last_exon = exons.iloc[-1]
                if last_exon['end'] < gene['end']:
                    introns.append({'iso': gene['iso'], 'chr': gene['chr'],
                        'start': last_exon['end'] + 1, 'end': gene['end'],
                        'type': 'intron', 'oid': gene['oid'], 'product': gene['product']})
                
                return pd.DataFrame(introns)

            # Ajout introns
            intron_dfs = df.groupby('oid')[df.columns].apply(generate_introns).reset_index(drop=True)
            if not intron_dfs.empty: df = pd.concat([df, intron_dfs], ignore_index=True)
            
            df.reset_index(drop=True, inplace=True)
            all_dfs.append(df)
            file_count += 1
            
        except Exception as e:
            print(f"Erreur lors du traitement du fichier {file_path.name}: {e}", file=sys.stderr)
            continue
    
    # Vérifier si des fichiers ont été traités avec succès
    if not all_dfs:
        print("Erreur: Aucun fichier n'a pu être traité correctement", file=sys.stderr)
        return False
        
    try:
        # Fusionner tous les DataFrames
        df = pd.concat(all_dfs, ignore_index=True)
        # Exporter le DataFrame filtré dans un fichier TSV
        df = df.sort_values(by=['iso', 'chr', 'start'])
        df.to_csv(outfile_path, index=False, sep='\t')
        
        print(f"Traitement terminé: {file_count} fichiers traités")
        print(f"Nombre total de gene: {len(df['oid'].unique())}")
        print(f"Résultat enregistré dans: {outfile_path}")
            
        return True
        
    except Exception as e:
        print(f"Erreur lors de la fusion ou de l'exportation des données: {e}", file=sys.stderr)
        return False


def main():
    # Configurer le parseur d'arguments
    parser = argparse.ArgumentParser(
        description="Traite des fichiers funannotate .gff et les combine en un seul fichier TSV.",
        epilog="Exemple: %(prog)s -i /chemin/vers/fichiers -o resultat.tsv"
    )
    parser.add_argument('-i', '--input', required=True, help="Chemin vers le répertoire contenant les fichiers .gff")
    parser.add_argument('-o', '--output', required=True, help="Chemin pour le fichier de sortie TSV")

    # Analyser les arguments
    args = parser.parse_args()
    # Exécuter la fonction principale
    success = create_df_all_gene(args.input, args.output)
    # Définir le code de sortie
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()