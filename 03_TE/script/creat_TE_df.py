#!/usr/bin/env python3
"""
Script pour traiter des fichiers RepeatMasker .out et les combiner en un seul fichier TSV.
Ce script compile les données d'éléments transposables, les filtre et les formate.
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


def create_df_all_TE(directory_path, outfile_path, verbose=False):
    """
    Crée un DataFrame unique à partir de tous les fichiers .out dans un répertoire.
    Args:
        directory_path (str): Chemin vers le répertoire contenant les fichiers .out
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
            if verbose: print(f"Répertoire de sortie créé: {outfile_dir}")
        except OSError as e:
            print(f"Erreur lors de la création du répertoire de sortie: {e}", file=sys.stderr)
            return False
    
    # Parcourir tous les fichiers .out dans le dossier
    out_files = list(Path(directory_path).glob('*.out'))
    if not out_files:
        print(f"Attention: Aucun fichier .out trouvé dans '{directory_path}'", file=sys.stderr)
        return False
    
    if verbose: print(f"Traitement de {len(out_files)} fichiers .out...")
    
    for file_path in out_files:
        try:
            if verbose: print(f"Traitement du fichier: {file_path.name}")   
            df = pd.read_csv(file_path, sep=r'\s+', engine='python', header=None, 
                          index_col=False, skip_blank_lines=True, skiprows=2)
            
            # Vérifier que le fichier contient des données
            if df.empty:
                print(f"Attention: Le fichier {file_path.name} est vide ou mal formaté", file=sys.stderr)
                continue
            
            # Attribuer les noms de colonnes originaux
            df.columns = ['score', 'div.', 'del.', 'ins.', 'sequence', 'begin', 'end', '(left)', 'strand', 
                'repeat', 'class/family', 'begin.1', 'end.1', '(left).1', 'ID']
            
            # Reformater le DataFrame
            df['iso'] = file_path.stem 
            df['sequence'] = df['sequence'].apply(process_chromosome_name)
            df = df[['iso', 'sequence', 'begin', 'end', 'class/family', 'ID', 'div.']]
            df.columns = ['iso', 'chr', 'start', 'end', 'class', 'oid', 'div']
            
            all_dfs.append(df)
            file_count += 1
            
        except Exception as e:
            print(f"Erreur lors du traitement du fichier {file_path.name}: {e}", file=sys.stderr)
            continue
    
    # Vérifier si des fichiers ont été traités avec succès
    if not all_dfs:
        print("Erreur: Aucun fichier n'a pu être traité correctement", file=sys.stderr)
        return False
        
    # Fusionner tous les DataFrames
    try:
        df = pd.concat(all_dfs, ignore_index=True)
        
        # Filtrer et regrouper les classes d'éléments transposables
        df = df[(df['class'] != 'Low_complexity') & (df['class'] != 'Simple_repeat')]
        df['class'] = df['class'].replace(['Simple_repeat', 'Low_complexity', 'ClassI/Unclassified', 'RC/Helitron'], 'Others')
        df['class'] = df['class'].replace({r'^LTR.*$': 'LTR', r'^DNA.*$': 'DNA', r'^LINE.*$': 'LINE'}, regex=True)

        # Exporter le DataFrame filtré dans un fichier TSV
        df = df.sort_values(by=['iso', 'chr', 'start'])
        df.to_csv(outfile_path, index=False, sep='\t')
        
        if verbose:
            print(f"Traitement terminé: {file_count} fichiers traités")
            print(f"Nombre total d'éléments transposables: {len(df)}")
            print(f"Résultat enregistré dans: {outfile_path}")
            
        return True
        
    except Exception as e:
        print(f"Erreur lors de la fusion ou de l'exportation des données: {e}", file=sys.stderr)
        return False


def main():
    # Configurer le parseur d'arguments
    parser = argparse.ArgumentParser(
        description="Traite des fichiers RepeatMasker .out et les combine en un seul fichier TSV.",
        epilog="Exemple: %(prog)s -i /chemin/vers/fichiers -o resultat.tsv"
    )
    parser.add_argument('-i', '--input', required=True, help="Chemin vers le répertoire contenant les fichiers .out")
    parser.add_argument('-o', '--output', required=True, help="Chemin pour le fichier de sortie TSV")
    parser.add_argument('-v', '--verbose', action='store_true', help="Afficher des informations détaillées pendant le traitement")
    
    # Analyser les arguments
    args = parser.parse_args()
    # Exécuter la fonction principale
    success = create_df_all_TE(args.input, args.output, args.verbose)
    # Définir le code de sortie
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()