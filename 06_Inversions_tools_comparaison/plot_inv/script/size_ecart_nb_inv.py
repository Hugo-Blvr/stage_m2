import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from collections import defaultdict

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Fonction pour regrouper les inversions proches selon un seuil
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


def plot_nb_size(df,thresholds, out_path):
    results = []

    for threshold in thresholds:
        # Fusionner les inversions selon le seuil actuel
        merged_df = merge_inversions(df, threshold)
        
        # Calculer la taille de chaque inversion
        merged_df['size'] = merged_df['end'] - merged_df['start'] + 1
        
        # Stocker les résultats
        results.append({
            'threshold': threshold,
            'nb_inversions': len(merged_df),
            'mean_size': merged_df['size'].mean(),
            'median_size': merged_df['size'].median(),
            'min_size': merged_df['size'].min(),
            'max_size': merged_df['size'].max(),
            'total_size': merged_df['size'].sum()
        })

        if threshold == 0 or threshold == 50: 
            mean_size = round(merged_df['size'].mean(),2)
            median_size = round(merged_df['size'].median(),2)
            print (f"Threshold {threshold}: nb = {len(merged_df)} ; mean_size = {mean_size}; median_size = {median_size}")

    # Créer un dataframe à partir des résultats
    results_df = pd.DataFrame(results)

    # Créer les graphiques
    plt.figure(figsize=(12, 10))

    # Graphique 1: Nombre d'inversions en fonction du seuil
    plt.subplot(2, 1, 1)
    plt.plot(results_df['threshold'], results_df['nb_inversions'], marker='o', linestyle='-')
    plt.title('Nombre d\'inversions en fonction du seuil d\'écart maximum')
    plt.xlabel('Seuil d\'écart maximum (bases)')
    plt.ylabel('Nombre d\'inversions')
    plt.grid(True)

    # Graphique 2: Taille moyenne et médiane des inversions en fonction du seuil
    plt.subplot(2, 1, 2)
    plt.plot(results_df['threshold'], results_df['mean_size'], marker='o', linestyle='-', label='Taille moyenne')
    plt.plot(results_df['threshold'], results_df['median_size'], marker='s', linestyle='--', label='Taille médiane')
    plt.title('Taille des inversions en fonction du seuil d\'écart maximum')
    plt.xlabel('Seuil d\'écart maximum (bases)')
    plt.ylabel('Taille des inversions (bases)')
    plt.legend()
    plt.grid(True)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300)


input_files = {"cactus_multialgn": "/home/hugob/stage_m2/06_bench_tools_inv_calling/format_inv/output/inv_cactus_multialgn.tsv", 
               "cactus_PtoMalgn": "/home/hugob/stage_m2/06_bench_tools_inv_calling/format_inv/output/inv_cactus_PtoMalgn.tsv", 
               "syri_minimap_PtoMalgn":"/home/hugob/stage_m2/06_bench_tools_inv_calling/format_inv/output/inv_syri_minimap_PtoMalgn.tsv",
               "syri_nucmer_PtoMalgn": "/home/hugob/stage_m2/06_bench_tools_inv_calling/format_inv/output/inv_syri_nucmer_PtoMalgn.tsv"}
output_dir = "/home/hugob/stage_m2/06_bench_tools_inv_calling/plot_inv/output/nb_size_ecart"


thresholds = list(range(0, 200, 5))
for input_file in input_files.keys(): 
    df = pd.read_csv(input_files[input_file], sep='\t')
    df = df[df['iso'] == 'Gd_00293-aad']

    out_path = f"{output_dir}/th_nb_size_{input_file}.png"
    print (f"Traitement {input_file}")
    #plot_nb_size(df,thresholds,out_path)
    print()


def analyze_inversions(input_files, output_dir, threshold = 0):
    """
    Analyse des fichiers TSV contenant des positions d'inversions par chromosomes d'isolats
    et génère deux visualisations différentes avec des échelles harmonisées par isolat.
    
    Args:
        input_files (dict): Dictionnaire contenant les chemins des fichiers TSV d'inversions
        output_dir (str): Chemin du répertoire de sortie pour les figures
        threshold (int): Seuil pour fusionner les inversions proches
    """
    
    # Charger les données depuis les fichiers TSV
    dfs = {}

    # Définir une palette de couleurs professionnelle pour chaque outil
    tools = ["syri_minimap_PtoMalgn", "syri_nucmer_PtoMalgn", "cactus_PtoMalgn", "cactus_multialgn"]
    colors = {
        "syri_minimap_PtoMalgn": "#3498db",  # Bleu
        "syri_nucmer_PtoMalgn": "#2ecc71",   # Vert
        "cactus_PtoMalgn": "#e74c3c",        # Rouge
        "cactus_multialgn": "#9b59b6"        # Violet
    }

    # Créer un dictionnaire pour les marqueurs de couleur à utiliser dans les plots
    color_dict = {}
    for tool in tools:
        color_dict[tool] = colors[tool]

    for tool in tools:
        file_path = input_files[tool]

        dfs[tool] = pd.read_csv(file_path, sep='\t')
        dfs[tool] = dfs[tool][dfs[tool]['iso'] == 'Gd_00293-aad']
        dfs[tool] = merge_inversions(dfs[tool], threshold)
        # Ajouter une colonne pour identifier l'outil
        dfs[tool]['tool'] = tool
 
    # Fusionner tous les dataframes
    all_data = pd.concat(dfs.values(), ignore_index=True)
    
    # Calculer la taille des inversions
    all_data['size'] = all_data['end'] - all_data['start']
    
    # Identifier les isolats communs à tous les jeux de données
    common_isos = set.intersection(*[set(df['iso']) for df in dfs.values()])

    # Filtrer les données pour ne garder que les isolats communs
    common_data = all_data[all_data['iso'].isin(common_isos)]
    
    # Calculer les écarts entre inversions consécutives pour chaque isolat et chromosome
    def calculate_gaps(df):
        gaps = []
        # Grouper par isolat, chromosome et outil
        for (iso, chrom, tool), group in df.groupby(['iso', 'chr', 'tool']):
            # Trier par position de début
            sorted_group = group.sort_values('start')
            if len(sorted_group) > 1:
                # Calculer l'écart entre la fin d'une inversion et le début de la suivante
                prev_end = sorted_group['end'].iloc[:-1].values
                next_start = sorted_group['start'].iloc[1:].values
                gap = next_start - prev_end
                
                # Créer un DataFrame avec les écarts
                gaps_df = pd.DataFrame({
                    'iso': iso,
                    'chr': chrom,
                    'tool': tool,
                    'gap': gap
                })
                gaps.append(gaps_df)
        
        if gaps:
            return pd.concat(gaps, ignore_index=True)
        else:
            return pd.DataFrame(columns=['iso', 'chr', 'tool', 'gap'])
    
    gaps_data = calculate_gaps(common_data)
    
    # 1. Violin plot des écarts entre les inversions avec les couleurs définies
    n_isos = len(common_isos)
    n_cols = min(3, n_isos)
    n_rows = (n_isos + n_cols - 1) // n_cols  # Arrondi supérieur pour le nombre de lignes
    
    fig1, axes1 = plt.subplots(n_rows, n_cols, figsize=(6*n_cols, 5*n_rows))
    axes1_flat = np.array(axes1).flatten() if isinstance(axes1, np.ndarray) else np.array([axes1])
    fig1.suptitle('Écarts entre inversions consécutives par isolat', fontsize=16)
    
    # Palette pour les violin plots avec les couleurs définies
    palette = color_dict

    for i, iso in enumerate(sorted(common_isos)):
        if i < len(axes1_flat):
            iso_gaps = gaps_data[gaps_data['iso'] == iso]
            if not iso_gaps.empty:
                ax = axes1_flat[i]
                sns.violinplot(data=iso_gaps, x='tool', y='gap', ax=ax, palette=palette)
                ax.set_title(f'Isolat: {iso}')
                ax.set_ylabel('Écart (pb)')
                ax.set_xlabel('tools')
                ax.set_yscale('log')  # Échelle logarithmique pour mieux visualiser les écarts
                ax.tick_params(axis='x', rotation=45)
            else:
                axes1_flat[i].set_visible(False)
    
    # Masquer les axes supplémentaires
    for j in range(i+1, len(axes1_flat)):
        axes1_flat[j].set_visible(False)
    
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    fig1.savefig(os.path.join(output_dir, f't{threshold}_ecarts_inversions.png'), dpi=300)
    plt.close(fig1)
    
    # 2. Histogrammes des tailles d'inversions - un histogramme par outil avec échelles harmonisées
    # Calculer le nombre de lignes et colonnes nécessaires (2 colonnes)
    n_tools = len(tools)
    n_cols_hist = 2
    n_rows_hist = (n_tools + n_cols_hist - 1) // n_cols_hist
    
    for iso in sorted(common_isos):
        iso_data = common_data[common_data['iso'] == iso]
        
        if not iso_data.empty:
            # Déterminer les limites globales pour l'axe x (taille) et y (fréquence) pour cet isolat
            all_sizes = []
            max_count = 0
            min_size_global = float('inf')
            max_size_global = 0
            
            for tool in tools:
                tool_data = iso_data[iso_data['tool'] == tool]
                if not tool_data.empty:
                    sizes = tool_data['size']
                    sizes = sizes[sizes > 0]
                    if not sizes.empty:
                        all_sizes.extend(sizes.tolist())
                        # Calculer l'histogramme pour connaître la hauteur maximale
                        if sizes.min() < sizes.max():
                            bins = np.logspace(np.log10(max(1, sizes.min())), np.log10(sizes.max()), 30)
                            counts, _ = np.histogram(sizes, bins=bins)
                            max_count = max(max_count, counts.max())
                            min_size_global = min(min_size_global, sizes.min())
                            max_size_global = max(max_size_global, sizes.max())
            
            # Vérifier si nous avons des données valides
            if all_sizes and min_size_global < max_size_global:
                fig3, axes3 = plt.subplots(n_rows_hist, n_cols_hist, figsize=(12, 4*n_rows_hist))
                axes3_flat = axes3.flatten() if isinstance(axes3, np.ndarray) else np.array([axes3])
                #fig3.suptitle(f'Distribution des tailles d\'inversions pour l\'isolat {iso}', fontsize=16)
                
                # Créer des bins communs pour tous les histogrammes
                bins = np.logspace(np.log10(max(1, min_size_global)), np.log10(max_size_global), 30)
                
                for i, tool in enumerate(tools):
                    if i < len(axes3_flat):
                        ax = axes3_flat[i]
                        tool_data = iso_data[iso_data['tool'] == tool]
                        
                        if not tool_data.empty:
                            sizes = tool_data['size']
                            sizes = sizes[sizes > 0]
                            
                            if not sizes.empty:
                                ax.hist(sizes, bins=bins, alpha=0.8, color=colors[tool])
                            
                            ax.set_title(f'{tool}')
                            ax.set_xlabel('Taille (pb)')
                            ax.set_ylabel('Fréquence')
                            ax.set_xscale('log')
                            
                            # Définir les mêmes limites pour tous les graphiques
                            ax.set_xlim(min_size_global * 0.9, max_size_global * 1.1)
                            ax.set_ylim(0, max_count * 1.1)  # Ajouter 10% de marge en haut
                        else:
                            ax.text(0.5, 0.5, f'Pas de données pour {tool}', 
                                    horizontalalignment='center', 
                                    verticalalignment='center',
                                    transform=ax.transAxes)
                            # Même avec absence de données, on définit les mêmes limites
                            ax.set_xscale('log')
                            ax.set_xlim(min_size_global * 0.9, max_size_global * 1.1)
                            ax.set_ylim(0, max_count * 1.1)
                
                # Masquer les axes supplémentaires si nécessaire
                for j in range(i+1, len(axes3_flat)):
                    axes3_flat[j].set_visible(False)
                    
                plt.tight_layout(rect=[0, 0, 1, 0.95])
                fig3.savefig(os.path.join(output_dir, f't{threshold}_distribution_tailles_inversions_{iso}.png'), dpi=300)
                plt.close(fig3)
    
    # 3. Distribution des écarts entre inversions - un histogramme par outil avec échelles harmonisées
    for iso in sorted(common_isos):
        iso_gaps = gaps_data[gaps_data['iso'] == iso]
        
        if not iso_gaps.empty:
            # Déterminer les limites globales pour l'axe x (écart) et y (fréquence) pour cet isolat
            all_gaps = []
            max_count = 0
            min_gap_global = float('inf')
            max_gap_global = 0
            
            for tool in tools:
                tool_gaps = iso_gaps[iso_gaps['tool'] == tool]
                if not tool_gaps.empty:
                    gaps = tool_gaps['gap']
                    gaps = gaps[gaps > 0]  # Filtrer les écarts négatifs ou nuls
                    if not gaps.empty:
                        all_gaps.extend(gaps.tolist())
                        # Calculer l'histogramme pour connaître la hauteur maximale
                        if gaps.min() < gaps.max():
                            bins = np.logspace(np.log10(max(1, gaps.min())), np.log10(gaps.max()), 30)
                            counts, _ = np.histogram(gaps, bins=bins)
                            max_count = max(max_count, counts.max())
                            min_gap_global = min(min_gap_global, gaps.min())
                            max_gap_global = max(max_gap_global, gaps.max())
            
            # Vérifier si nous avons des données valides
            if all_gaps and min_gap_global < max_gap_global:
                fig4, axes4 = plt.subplots(n_rows_hist, n_cols_hist, figsize=(12, 4*n_rows_hist))
                axes4_flat = axes4.flatten() if isinstance(axes4, np.ndarray) else np.array([axes4])
                #fig4.suptitle(f'Distribution des écarts entre inversions pour l\'isolat {iso}', fontsize=16)
                
                # Créer des bins communs pour tous les histogrammes
                bins = np.logspace(np.log10(max(1, min_gap_global)), np.log10(max_gap_global), 30)
                
                for i, tool in enumerate(tools):
                    if i < len(axes4_flat):
                        ax = axes4_flat[i]
                        tool_gaps = iso_gaps[iso_gaps['tool'] == tool]
                        
                        if not tool_gaps.empty:
                            gaps = tool_gaps['gap']
                            gaps = gaps[gaps > 0]  # Filtrer les écarts négatifs ou nuls
                            
                            if not gaps.empty:
                                ax.hist(gaps, bins=bins, alpha=0.8, color=colors[tool])
                            
                            ax.set_title(f'{tool}')
                            ax.set_xlabel('Écart (pb)')
                            ax.set_ylabel('Fréquence')
                            ax.set_xscale('log')
                            
                            # Définir les mêmes limites pour tous les graphiques
                            ax.set_xlim(min_gap_global * 0.9, max_gap_global * 1.1)
                            ax.set_ylim(0, max_count * 1.1)  # Ajouter 10% de marge en haut
                        else:
                            ax.text(0.5, 0.5, f'Pas de données pour {tool}', 
                                    horizontalalignment='center', 
                                    verticalalignment='center',
                                    transform=ax.transAxes)
                            # Même avec absence de données, on définit les mêmes limites
                            ax.set_xscale('log')
                            ax.set_xlim(min_gap_global * 0.9, max_gap_global * 1.1)
                            ax.set_ylim(0, max_count * 1.1)
                
                # Masquer les axes supplémentaires si nécessaire
                for j in range(i+1, len(axes4_flat)):
                    axes4_flat[j].set_visible(False)
                    
                plt.tight_layout(rect=[0, 0, 1, 0.95])
                fig4.savefig(os.path.join(output_dir, f't{threshold}_distribution_ecarts_inversions_{iso}.png'), dpi=300)
                plt.close(fig4)

    print(f"Les visualisations avec échelles harmonisées ont été enregistrées dans le répertoire {output_dir}")

analyze_inversions(input_files,output_dir,50)