import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
from collections import defaultdict
import venn_tools as venn4
from itertools import combinations
import math

inv_files = [
    "/home/hugob/stage_m2/06_bench_tools_inv_calling/format_inv/output/inv_cactus_PtoMalgn.tsv",
    "/home/hugob/stage_m2/06_bench_tools_inv_calling/format_inv/output/inv_cactus_multialgn.tsv",
    "/home/hugob/stage_m2/06_bench_tools_inv_calling/format_inv/output/inv_syri_nucmer_PtoMalgn.tsv",
    "/home/hugob/stage_m2/06_bench_tools_inv_calling/format_inv/output/inv_syri_minimap_PtoMalgn.tsv"
]

chr_len_file = "/home/hugob/stage_m2/01_Processing/output/chr_len.tsv"
chr_len_df = pd.read_csv(chr_len_file, sep='\t')
grouped_sizes = chr_len_df.groupby('iso')['size'].sum()
genome_sizes = {iso: grouped_sizes.get(iso, None) for iso in chr_len_df['iso'].unique()}

lineage_file = "/home/hugob/stage_m2/00_Data/Isolates.Master.Info.txt"
lineage_df = pd.read_csv(lineage_file, sep='\t', usecols=[1, 3])
lineage_df.columns = ['iso', 'lineage']


def merge_inversions(df, max_gap):
    """
    Fusionne les inversions avec un écart maximum spécifié.
    
    Args:
        df: DataFrame avec les colonnes ['iso', 'chr', 'start', 'end']
        max_gap: Écart maximum autorisé entre les inversions à fusionner
        
    Returns:
        DataFrame des inversions fusionnées
    """
    # Vérifier si le dataframe est vide
    if df.empty:
        return pd.DataFrame(columns=['iso', 'chr', 'start', 'end'])
    
    # Tri efficace en utilisant uniquement les colonnes nécessaires
    colonnes_requises = ['iso', 'chr', 'start', 'end']
    df_trie = df[colonnes_requises].sort_values(by=['iso', 'chr', 'start'])
    
    # Pré-allocation de la liste de résultats
    inversions_fusionnees = []
    
    # Parcours des données groupées plus efficacement
    iso_courant = None
    chr_courant = None
    debut_courant = None
    fin_courant = None
    
    for _, row in df_trie.iterrows():
        iso, chr_nom, debut, fin = row['iso'], row['chr'], row['start'], row['end']
        
        # Gestion de la première ligne ou d'un nouveau groupe
        if iso_courant != iso or chr_courant != chr_nom:
            # Sauvegarder la dernière inversion fusionnée du groupe précédent si elle existe
            if iso_courant is not None:
                inversions_fusionnees.append({
                    'iso': iso_courant,
                    'chr': chr_courant,
                    'start': debut_courant,
                    'end': fin_courant
                })
            
            # Commencer un nouveau groupe
            iso_courant = iso
            chr_courant = chr_nom
            debut_courant = debut
            fin_courant = fin
            continue
        
        # Vérifier si l'inversion actuelle peut être fusionnée avec la précédente
        if debut - fin_courant <= max_gap:
            fin_courant = max(fin_courant, fin)
        else:
            # Sauvegarder l'inversion fusionnée complétée et en commencer une nouvelle
            inversions_fusionnees.append({
                'iso': iso_courant,
                'chr': chr_courant,
                'start': debut_courant,
                'end': fin_courant
            })
            debut_courant = debut
            fin_courant = fin
    
    # Ne pas oublier d'ajouter la dernière inversion
    if iso_courant is not None:
        inversions_fusionnees.append({
            'iso': iso_courant,
            'chr': chr_courant,
            'start': debut_courant,
            'end': fin_courant
        })
    
    # Créer un DataFrame à partir du résultat avec des colonnes prédéfinies pour de meilleures performances
    df_resultat = pd.DataFrame(inversions_fusionnees, columns=colonnes_requises)
    return df_resultat

# Lire tous les fichiers et organiser les données par iso et chr
def make_data(inv_files,threshold=0):
    methode = ['_'.join(os.path.splitext(os.path.basename(f))[0].split('_')[1:]) for f in inv_files]
    all_data = {}

    for i, file_path in enumerate(inv_files):
        try:
            df = pd.read_csv(file_path, sep='\t')
            chromosomes_à_exclure = ['chr2003', 'chr2016', 'chr2017']   
            df = df[~df['chr'].isin(chromosomes_à_exclure)]
            #df = df[(df['iso'] == 'Gd_00293-aad')]
            df = merge_inversions(df, threshold)
            # Parcourir chaque ligne
            for _, row in df.iterrows():
                iso = row['iso']
                chr_name = row['chr']
                start = row['start']
                end = row['end']
                
                # Initialiser la structure si nécessaire
                if iso not in all_data:
                    all_data[iso] = {}
                if chr_name not in all_data[iso]:
                    all_data[iso][chr_name] = {}
                if methode[i] not in all_data[iso][chr_name]:
                    all_data[iso][chr_name][methode[i]] = []
                        
                # Ajouter les coordonnées de l'inversion
                all_data[iso][chr_name][methode[i]].append((start, end))
        
        except Exception as e:
            print(f"Erreur lors de la lecture de {file_path}: {str(e)}")
    return all_data


def calculate_all_inversions_intersections(all_data, genome_sizes):
    """
    Calcule toutes les intersections de bases inversées pour tous les isolats
    en incluant le pourcentage par rapport à la taille du génome.
    
    Args:
        all_data: Dictionnaire structuré {iso -> {chr -> {method -> [(start, end)]}}}
        chr_len_df: DataFrame contenant les colonnes 'iso', 'chr', 'size'
    
    Returns:
        Un dictionnaire {iso -> {key_string: formatted_count}}
    """    
    results = {}
    
    
    # Pour chaque isolat
    for iso, chr_data_all in all_data.items():
        iso_results = {}
        all_methods = set()
        
        # Collecter toutes les méthodes uniques pour cet isolat
        for chr_name, methods in chr_data_all.items():
            all_methods.update(methods.keys())    
        method_to_index = {method: i for i, method in enumerate(sorted(all_methods))}

        all_methods = ["syri_minimap_PtoMalgn", "syri_nucmer_PtoMalgn", "cactus_PtoMalgn", "cactus_multialgn"]
        method_to_index = {method: i for i, method in enumerate(all_methods)}
    
        # Préparer les ensembles de bases pour chaque méthode
        method_bases = {i: set() for i in range(len(all_methods))}



        
        # Remplir les ensembles pour chaque chromosome et méthode
        for chr_name, methods in chr_data_all.items():
            for method, regions in methods.items():
                method_idx = method_to_index[method]
                for start, end in regions:
                    method_bases[method_idx].update(range(start, end + 1))
        
        # Calculer et stocker le nombre de bases pour chaque méthode et les intersections
        for r in range(1, len(all_methods) + 1):
            for combo in combinations(range(len(all_methods)), r):
                if r == 1:
                    # Méthode unique
                    count = len(method_bases[combo[0]])
                else:
                    # Intersection de plusieurs méthodes
                    intersection = set.intersection(*(method_bases[i] for i in combo))
                    count = len(intersection)
                
                # Calculer et formater le pourcentage
                if genome_sizes[iso]:
                    percentage = (count / genome_sizes[iso]) * 100
                    formatted_count = f"{percentage:.2f}%"
                else:
                    formatted_count = f"{count} (% inconnu)"
                
                iso_results[combo] = formatted_count
                
        results[iso] = iso_results
    return results


def plot_inv(all_data,chr_len_df,output_dir):
    methodes = set()
    for iso in all_data:
        for chr_name in all_data[iso]:
            for methode in all_data[iso][chr_name]:methodes.add(methode)
    methodes = sorted(methodes)
    methodes = ["syri_minimap_PtoMalgn", "syri_nucmer_PtoMalgn", "cactus_PtoMalgn", "cactus_multialgn"]


    # Générer un plot pour chaque isolat
    for iso in all_data:
        # Compter le nombre de chromosomes pour cet isolat pour définir la disposition
        chrs = list(all_data[iso].keys())
        n_chrs = len(chrs)
        
        # Définir une disposition de sous-graphiques
        n_cols = min(3, n_chrs)  # 3 colonnes maximum
        n_rows = int(np.ceil(n_chrs / 3))  # Calcul du nombre de lignes nécessaires

        # Créer une figure suffisamment grande
        plt.figure(figsize=(20, 3 * n_rows))
        
        # Couleurs pour les différentes méthodes
        colors = {
                "syri_minimap_PtoMalgn": "#3498db",  # Bleu
                "syri_nucmer_PtoMalgn": "#2ecc71",   # Vert
                "cactus_PtoMalgn": "#e74c3c",        # Rouge
                "cactus_multialgn": "#9b59b6"        # Violet
            }
        
        # Créer un sous-graphique pour chaque chromosome
        for c_idx, chr_name in enumerate(sorted(chrs)):
            ax = plt.subplot(n_rows, n_cols, c_idx + 1)
            
            # Trouver la longueur maximale du chromosome pour définir l'axe x
            chr_size = int(chr_len_df[(chr_len_df['iso'] == iso) & (chr_len_df['chr'] == chr_name)]['size'].iloc[0])
            # Créer des représentations graphiques pour chaque méthode
            y_offset = 0  # Décalage vertical pour empiler les méthodes
            for method in methodes[::-1]:
                if method in all_data[iso][chr_name]:
                    inv_list = all_data[iso][chr_name][method]
                    
                    # Tracer chaque inversion
                    for start, end in inv_list:
                        ax.plot([start, end], [y_offset, y_offset], color=colors[method], 
                                linewidth=2, solid_capstyle='butt')
                    
                    y_offset += 1  # Décaler vers le haut pour la prochaine méthode
            
            # Configurer le graphique
            ax.set_xlim(0, chr_size)
            ax.set_ylim(-0.5, len(methodes) - 0.5)
            ax.set_title(f"{chr_name}")
            ax.set_xlabel("Position (bp)")
            ax.grid(True, axis='x', linestyle='--', alpha=0.7)
            
            # Ajouter les étiquettes des méthodes sur l'axe y
            ax.set_yticks(range(len(methodes)))
            ax.set_yticklabels(methodes[::-1])
        
        # Titre global pour la figure
        plt.suptitle(f"Inversions détectées pour {iso}", fontsize=16)
        plt.tight_layout(rect=[0, 0, 1, 0.96])  # Ajuster pour laisser de la place au titre
        
        # Enregistrer la figure
        plt.savefig(os.path.join(output_dir, f"t0_inversions_{iso}.png"), dpi=300)
        plt.close()
        
        print(f"Image générée pour {iso}")

    print("Terminé! Toutes les plots ont été générées.")

output_dir = "/home/hugob/stage_m2/06_bench_tools_inv_calling/plot_inv/output/distrib_by_tools"
all_data = make_data(inv_files,50)
#plot_inv(all_data,chr_len_df,output_dir)

output_dir = "/home/hugob/stage_m2/06_bench_tools_inv_calling/plot_inv/output/venn_plot"
venn = True
if venn:
    methodes = set()
    for iso in all_data:
        for chr_name in all_data[iso]:
            for methode in all_data[iso][chr_name]:methodes.add(methode)
    methodes = sorted(methodes)
    methodes = ["syri_minimap_PtoMalgn", "syri_nucmer_PtoMalgn", "cactus_PtoMalgn", "cactus_multialgn"]
    labels = list(methodes)

    intersections = calculate_all_inversions_intersections(all_data, genome_sizes)


    isos = list(all_data.keys())
    n = len(isos)
    cols = 5
    rows = math.ceil(n / cols)

    # Crée la figure principale

    fig, axes = plt.subplots(rows, cols, figsize=(cols * 5, rows * 5))
    if rows == 1 and cols == 1:  # Un seul subplot (1x1)
        axes = np.array([axes])
    elif rows == 1:  # Une seule ligne
        axes = np.array(axes).reshape(1, cols).flatten()
    elif cols == 1:  # Une seule colonne
        axes = np.array(axes).reshape(rows, 1).flatten()
    else:  # Plusieurs lignes et colonnes
        axes = axes.flatten()


    for idx, iso in enumerate(isos):
        ax = axes[idx]
        genome_size = genome_sizes[iso]
        title = f"{iso} (genome size: {genome_size})"
        data = intersections[iso]

        venn4.venn_4(labels, data, title, ax=ax)
    # Supprime les axes en trop si n n'est pas un multiple de cols
    for j in range(idx + 1, len(axes)):
        fig.delaxes(axes[j])

    fig.legend(labels=labels, loc='lower center', ncol=4, frameon=False)
    plt.tight_layout(rect=[0, 0.05, 1, 1])  # laisse de la place en bas pour la légende

    output_path = f"{output_dir}/t50_all_venn_inv.png"
    plt.savefig(output_path, dpi=600, bbox_inches='tight')
    plt.close()
