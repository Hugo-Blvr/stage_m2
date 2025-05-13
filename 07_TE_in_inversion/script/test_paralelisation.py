import pandas as pd
from collections import defaultdict
import numpy as np
import seaborn as sns
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import os
import multiprocessing as mp
from functools import partial
import time

def del_invTE(df_inv, df_te, methode):
    inv_df = df_inv.copy()
    te_df = df_te.copy()
    
    # Précalcul des longueurs d'inversion pour éviter les calculs répétés
    inv_df['length'] = inv_df['end'] - inv_df['start']
    
    # Création des index pour la jointure spatiale
    inv_df.reset_index(inplace=True)  # Préserver l'index original
    inv_df.rename(columns={'index': 'original_idx'}, inplace=True)
    
    # Création des conditions de jointure spatiale
    # on regroupe par chromosome et iso pour accélérer
    to_remove = set()
    
    # Traiter chaque groupe (chr, iso) séparément
    print(f"Tri des inversions : {methode} ...")
    for (chr_val, iso_val), inv_group in inv_df.groupby(['chr', 'iso']):
        
        # Filtrer les TE pertinents en une seule fois
        matching_tes = te_df[(te_df['chr'] == chr_val) & (te_df['iso'] == iso_val)].copy()
        
        if matching_tes.empty:
            continue
        
        # Calculer les longueurs des TE
        matching_tes['te_length'] = matching_tes['end'] - matching_tes['start']
        
        # Vectorisation des calculs de chevauchement
        for _, inv in inv_group.iterrows():
            # Filtrage rapide: TE qui pourraient chevaucher l'inversion
            possible_overlaps = matching_tes[
                (matching_tes['start'] <= inv['end']) & 
                (matching_tes['end'] >= inv['start'])
            ]
            
            if possible_overlaps.empty:
                continue
                
            # Calcul vectorisé des chevauchements
            overlap_starts = np.maximum(inv['start'], possible_overlaps['start'])
            overlap_ends = np.minimum(inv['end'], possible_overlaps['end'])
            overlap_lengths = np.maximum(0, overlap_ends - overlap_starts)
            
            # Calcul des deux pourcentages de chevauchement (bidirectionnel)
            inv_overlap_percentages = (overlap_lengths / inv['length']) * 100
            te_overlap_percentages = (overlap_lengths / possible_overlaps['te_length']) * 100
            
            # Vérifier si au moins un TE a un chevauchement bidirectionnel ≥ 95%
            both_match = (inv_overlap_percentages >= 95) & (te_overlap_percentages >= 95)
            
            if both_match.any():
                to_remove.add(inv['original_idx'])
    
    nb_filtered = len(to_remove)
    filtered_inv_df = df_inv.drop(to_remove)  # Conserver la colonne 'iso'
    pourcentage = (nb_filtered / len(inv_df)) * 100
    print(f"{methode}: nombre d'inversions filtrées: {nb_filtered} ({pourcentage:.2f}%)")
    
    return filtered_inv_df

# Définition des fichiers d'entrée et de sortie
input_files = {"Cactus_PtoMalgn_maf": "/home/hugob/stage_m2/06_Inversions_tools_comparaison/format_inv/output/inv_cactus_PtoMalgn_maf.tsv", 
               "Syri_minimap_PtoMalgn":"/home/hugob/stage_m2/06_Inversions_tools_comparaison/format_inv/output/inv_syri_minimap_PtoMalgn.tsv"}
te_file = "/home/hugob/stage_m2/03_TE/output/all_TE.tsv"
#te_file = "/home/hugob/stage_m2/04_gene/output/format_gene.tsv"
chr_len_file = "/home/hugob/stage_m2/01_Processing/output/chr_len.tsv"
output_path = "/home/hugob/stage_m2/07_TE_in_inversion/output"

# Filtres
chromosomes_exclus = ['chr2003', 'chr2016', 'chr2017']
isos = ['Gd_00293-aad']
mapons = ['all']
inv_datas = {}


# Lecture des données TE
te_data = pd.read_csv(te_file, sep='\t')
te_data['iso'] = te_data['iso'].apply(lambda x: '_'.join(x.split('_')[:2]))
te_data = te_data.query("chr not in @chromosomes_exclus").sort_values(by='start')


# Lecture des fichiers d'inversions
for methode, patt_file in input_files.items():
    methode_name = methode.split('_')[0]
    df = pd.read_csv(patt_file, sep='\t')
    df = del_invTE(df,te_data,methode_name)
    inv_datas[methode_name] = df.query("chr not in @chromosomes_exclus and mapon in @mapons").sort_values(by='start')


# Lecture des longueurs de chromosomes
chr_len_df = pd.read_csv(chr_len_file, sep='\t')
chr_len_df = chr_len_df.query("chr not in @chromosomes_exclus")    
chr_size_dict = defaultdict(lambda: defaultdict(dict))
for iso, chr, size in zip(chr_len_df['iso'], chr_len_df['chr'], chr_len_df['size']):
    chr_size_dict[iso][chr] = size




def Te_density(df_inv, df_te, window_size_half, range_5to3, full_df, chr_size, inv_max_size):
    """
    Calcule la densité des éléments transposables autour des inversions.
    """
    te_tx_recover_5to3 = defaultdict(list)
    te_tx_recover_3to5 = defaultdict(list)
    range_3to5 = [-x for x in range_5to3[::-1]]

    if max(range_5to3) == inv_max_size // 2: 
        range_3to5 = range_3to5[1::]

    for _, inv in df_inv.iterrows():
        full_df_filetred = full_df[(full_df['start'] != inv['start']) & (full_df['end'] != inv['end'])].copy()
        mid_size_inv = (inv['end'] - inv['start']) // 2

        # Traitement du côté 5'→3'
        for range_val in range_5to3:
            if range_val not in te_tx_recover_5to3:te_tx_recover_5to3[range_val] = []
            pos = inv['start'] + range_val
            if range_val > mid_size_inv: continue 
            window_start = pos - window_size_half
            window_end = pos + window_size_half
            mask = (full_df_filetred['start'] < window_end) & (full_df_filetred['end'] > window_start) 
            a = True
            if (a or range_val >= 0) and window_start >= 0: 
                in_window = (df_te['start'] <= window_end) & (df_te['end'] >= window_start)
                te_in_window = df_te[in_window]

                # Calculer le chevauchement de chaque TE avec la fenêtre
                overlap_start = np.maximum(te_in_window['start'], window_start)
                overlap_end = np.minimum(te_in_window['end'], window_end)
                overlap_length = overlap_end - overlap_start + 1  # +1 si positions inclusives

                # Additionner toutes les longueurs de chevauchement
                total_bases_in_te = np.sum(overlap_length)
                density = total_bases_in_te / (2 * window_size_half)
                te_tx_recover_5to3[range_val].append(density)
        
        # Traitement du côté 3'→5'
        for range_val in range_3to5:
            if range_val not in te_tx_recover_3to5:  
                te_tx_recover_3to5[range_val] = []
            pos = inv['end'] + range_val 
            if -range_val > mid_size_inv: continue 
            window_start = pos - window_size_half
            window_end = pos + window_size_half
            mask = (full_df_filetred['start'] < window_end) & (full_df_filetred['end'] > window_start) 
            a = True
            if (a or range_val <= 0) and window_end <= chr_size: 
                in_window = (df_te['start'] <= window_end) & (df_te['end'] >= window_start)
                te_in_window = df_te[in_window]

                # Calculer le chevauchement de chaque TE avec la fenêtre
                overlap_start = np.maximum(te_in_window['start'], window_start)
                overlap_end = np.minimum(te_in_window['end'], window_end)
                overlap_length = overlap_end - overlap_start + 1  # +1 si positions inclusives

                # Additionner toutes les longueurs de chevauchement
                total_bases_in_te = np.sum(overlap_length)
                density = total_bases_in_te / (2 * window_size_half)
                te_tx_recover_3to5[range_val].append(density)

    # Calcul des moyennes et comptage
    te_tx_recover_5to3 = {key: [np.mean(values) if values else 0, len(values) if values else 0] 
                           for key, values in te_tx_recover_5to3.items()}
    te_tx_recover_3to5 = {key: [np.mean(values) if values else 0, len(values) if values else 0] 
                           for key, values in te_tx_recover_3to5.items()}
    
    pos = sorted(te_tx_recover_5to3) + sorted(te_tx_recover_3to5)
    value = [te_tx_recover_5to3[keys] for keys in sorted(te_tx_recover_5to3)] + \
            [te_tx_recover_3to5[keys] for keys in sorted(te_tx_recover_3to5)]

    # Formatage des étiquettes de position
    idx_zero = [i for i, val in enumerate(pos) if val == 0]
    pos[idx_zero[0]], pos[idx_zero[1]] = 'inv_start', 'inv_end'
    
    if max(range_5to3) == inv_max_size // 2: 
        idx_mid = [i for i, val in enumerate(pos) if val == inv_max_size // 2]
        pos[idx_mid[0]] = 'mid_inv'

    start_index = end_index = -1
    for i, item in enumerate(pos):
        if item == 'inv_start': start_index = i
        elif item == 'inv_end':
            end_index = i
            break 

    pos = [f'flank {x}' for x in pos[:start_index]] + \
          pos[start_index:end_index+1] + \
          [f'flank +{x}' for x in pos[end_index+1:]]

    return pos, value

def process_iso_chr(args):
    """
    Fonction pour traiter une combinaison isolat-chromosome spécifique.
    Cette fonction sera exécutée en parallèle.
    """
    methode, iso, chrom, df_iso_chr_inv, te_data, window_size_half, range_5to3, full_df, chr_size, inv_max_size = args
    
    print(f"Traitement de {methode} {iso} {chrom}...")
    
    # Appel à la fonction Te_density
    pos, value = Te_density(df_iso_chr_inv, te_data, window_size_half, 
        range_5to3, full_df, chr_size, inv_max_size)
    
    return methode, iso, chrom, pos, value

def make_visualization(all_result, output_path):
    """
    Fonction améliorée pour visualiser les résultats par isolat et créer un graphique global
    """
    # Structures de données pour stocker les résultats
    isolat_data = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    global_data = defaultdict(lambda: defaultdict(list))  # Pour le graphique global
    
    # Collecter les données pour chaque isolat et pour le graphique global
    for methode, iso_dict in all_result.items():
        for iso, chr_dict in iso_dict.items():
            for chr, data in chr_dict.items():
                # Si des données existent pour cette combinaison
                if data and len(data) == 2:
                    positions, values = data
                    for i, pos in enumerate(positions):
                        if i < len(values):
                            mean_density, count = values[i]
                            # Données pour l'isolat spécifique
                            isolat_data[iso][methode][pos].append((mean_density, count))
                            # Données pour le graphique global
                            global_data[methode][pos].append((mean_density, count))
    
    # Créer une visualisation pour chaque isolat
    for iso, methode_dict in isolat_data.items(): create_density_plot(methode_dict, f'Densité moyenne des TEs autour des inversions - Isolat {iso}', os.path.join(output_path, f'test_SF_te_density_isolat_{iso}.png'))
    
    # Créer le graphique global pour tous les isolats
    create_density_plot(global_data, 'Densité moyenne des TEs autour des inversions - Tous isolats', 
                        os.path.join(output_path, f'sort_SF_te_density_all_isolats.png'))

def create_density_plot(data_dict, title, output_file):
    """
    Fonction helper pour créer un graphique de densité
    """
    # Préparer les données pour le graphique
    positions_combined = None
    density_by_method = {}
    count_by_method = {}
    
    for methode, pos_data in data_dict.items():
        # Calculer la moyenne des densités pour chaque position
        positions = list(pos_data.keys())
        if not positions_combined: positions_combined = positions
        
        densities = []
        counts = []
        
        for pos in positions:
            values_at_pos = pos_data[pos]
            # Calculer la moyenne des moyennes pour chaque position
            mean_density = np.mean([v[0] for v in values_at_pos]) if values_at_pos else 0
            total_count = sum([v[1] for v in values_at_pos]) if values_at_pos else 0
            
            densities.append(mean_density)
            counts.append(total_count)
        
        density_by_method[methode] = densities
        count_by_method[methode] = counts
    
    # Créer la figure à deux sous-graphiques
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 12), gridspec_kw={'height_ratios': [3, 1]})
    
    # Graphique principal: densité des TE
    for methode, densities in density_by_method.items():
        ax1.plot(range(len(positions_combined)), densities, 
                label=methode, marker='o' if methode == 'Cactus' else 's',
                linestyle='-', linewidth=2, markersize=4)
    
    # Ajouter des lignes verticales pour les points d'inversion
    start_index = positions_combined.index('inv_start') if 'inv_start' in positions_combined else None
    end_index = positions_combined.index('inv_end') if 'inv_end' in positions_combined else None
    
    if start_index is not None: ax1.axvline(x=start_index, color='gray', linestyle='--', alpha=0.7)
    if end_index is not None: ax1.axvline(x=end_index, color='gray', linestyle='--', alpha=0.7)
    
    ax1.set_ylabel('Densité moyenne des TEs')
    ax1.set_ylim(0, 1)
    ax1.set_title(title)
    ax1.legend(title='Method')
    ax1.grid(True, linestyle='--', alpha=0.7)
    
    # Graphique inférieur: nombre d'inversions
    bar_width = 0.35
    index = np.arange(len(positions_combined))
    
    for i, (methode, counts) in enumerate(count_by_method.items()):
        ax2.bar(index + (i * bar_width - 0.5 * bar_width * (len(count_by_method) - 1)), 
               counts, 
               bar_width, 
               label=methode)
    
    ax2.set_xlabel('Position relative')
    ax2.set_ylabel('Nombre d\'inversions')
    ax2.legend(title='Method', loc='upper center')
    
    # Ajuster les étiquettes - uniquement sur le graphique du bas
    start_index = positions_combined.index('inv_start') if 'inv_start' in positions_combined else None
    end_index = positions_combined.index('inv_end') if 'inv_end' in positions_combined else None
    
    # Enlever les étiquettes du graphique du haut mais garder les ticks
    ax1.set_xticks(range(len(positions_combined)))
    ax1.set_xticklabels(['' for _ in range(len(positions_combined))])
    
    # Déterminer un espacement homogène basé sur les points d'ancrage
    max_ticks = 30  # Nombre maximum d'étiquettes souhaité
    
    if start_index is not None and end_index is not None:
        # Calculer les espaces avant, pendant et après l'inversion
        space_before = start_index
        space_inversion = end_index - start_index
        space_after = len(positions_combined) - end_index - 1
        
        # Calculer des pas homogènes pour chaque section
        step_before = max(1, space_before // (max_ticks // 3)) if space_before > 0 else 1
        step_inversion = max(1, space_inversion // (max_ticks // 3)) if space_inversion > 0 else 1
        step_after = max(1, space_after // (max_ticks // 3)) if space_after > 0 else 1
        
        # Créer les indices à afficher
        indices_before = list(range(0, start_index, step_before))
        if indices_before and indices_before[-1] != start_index - 1 and start_index - 1 >= 0:
            indices_before.append(start_index)  # Ajouter le point juste avant start_index
            
        indices_inversion = [start_index] + list(range(start_index + step_inversion, end_index, step_inversion)) + [end_index]
        
        indices_after = list(range(end_index + step_after, len(positions_combined), step_after))
        if indices_after and indices_after[0] != end_index + 1 and end_index + 1 < len(positions_combined):
            indices_after.insert(0, end_index)  # Ajouter le point juste après end_index
            
        indices_to_show = sorted(set(indices_before + indices_inversion + indices_after))
    else:
        # Fallback si start_index ou end_index n'est pas trouvé
        step = max(1, len(positions_combined) // max_ticks)
        indices_to_show = list(range(0, len(positions_combined), step))
    
    # Appliquer uniquement au graphique du bas (barplot)
    ax2.set_xticks([i for i in indices_to_show])
    ax2.set_xticklabels([positions_combined[i] for i in indices_to_show], rotation=90)
    
    plt.tight_layout()
    
    # Sauvegarder le graphique
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Figure sauvegardée dans {output_file}")
     
def visualize_te_inversion_relative_density_parallel(datas, te_data, chr_size_dict, output_dir, 
                                    inv_min_size=0, inv_max_size=100000, flanking_region_size=50000, step_size=500,
                                    num_processes=None):
    """
    Version parallélisée de la fonction d'analyse de densité des TE autour des inversions.
    
    Args:
        datas: Dictionnaire contenant les données d'inversion pour chaque méthode
        te_data: DataFrame contenant les données des éléments transposables
        chr_size_dict: Dictionnaire contenant les tailles de chromosomes
        output_dir: Répertoire de sortie pour les graphiques
        inv_min_size: Taille minimale des inversions à considérer
        inv_max_size: Taille maximale des inversions à considérer
        flanking_region_size: Taille de la région flanquante à analyser
        step_size: Taille du pas pour l'analyse
        num_processes: Nombre de processus à utiliser (None = nombre de cœurs disponibles)
    """
    # Définir le nombre de processus (utiliser tous les cœurs disponibles par défaut)
    if num_processes is None: num_processes = mp.cpu_count()
    
    print(f"Utilisation de {num_processes} processus parallèles")
    
    # Calculer les positions relatives à analyser
    range_5to3 = [-dist for dist in reversed(range(0, flanking_region_size + 1, step_size))] + \
                 [dist for dist in range(step_size, inv_max_size // 2 + 1, step_size)]    
    
    # Organiser les données TE par isolat et chromosome
    te_by_iso_chr = {(iso, chrom): chr_group
                     for iso, group in te_data.groupby('iso')
                     for chrom, chr_group in group.groupby('chr')}
    
    # Initialiser le dictionnaire pour stocker les résultats
    all_result = defaultdict(lambda: defaultdict(dict))
    window_size_half = step_size // 2 
    
    # Chronométrer l'exécution
    start_time = time.time()
    
    # Traiter chaque méthode
    for methode in sorted(datas.keys(), reverse=True):
        print(f"Traitement de {methode}...")
        full_df = datas[methode]
        df = full_df.copy()
        # Filtrage optionnel par taille d'inversion
        #df = full_df.query("(end - start) > @inv_min_size and (end - start) < @inv_max_size").copy()
        df['size'] = df['end'] - df['start']
        
        # Préparer les arguments pour le traitement parallèle
        args_list = []
        
        for iso, df_iso_inv in df.groupby('iso'):
            for chrom, df_iso_chr_inv in df_iso_inv.groupby('chr'):
                if (iso, chrom) in te_by_iso_chr:
                    args_list.append((
                        methode, 
                        iso, 
                        chrom, 
                        df_iso_chr_inv, 
                        te_by_iso_chr[(iso, chrom)], 
                        window_size_half, 
                        range_5to3, 
                        full_df.query("iso == @iso and chr == @chrom"), 
                        chr_size_dict[iso][chrom], 
                        inv_max_size
                    ))
        
        # Exécuter le traitement en parallèle
        with mp.Pool(processes=num_processes) as pool: results = pool.map(process_iso_chr, args_list)
        
        # Traiter les résultats
        for methode_res, iso, chrom, pos, value in results: all_result[methode_res][iso][chrom] = [pos, value]
    
    # Afficher le temps d'exécution
    end_time = time.time()
    print(f"Temps d'exécution total: {end_time - start_time:.2f} secondes")
    
    # Générer les visualisations
    make_visualization(all_result, output_dir)
    
    return all_result

# Point d'entrée principal
if __name__ == "__main__":
    # Exécuter l'analyse parallélisée
    results = visualize_te_inversion_relative_density_parallel(
        inv_datas, 
        te_data, 
        chr_size_dict, 
        output_path, 
        num_processes=mp.cpu_count() - 1  # Utiliser tous les cœurs sauf un
    )