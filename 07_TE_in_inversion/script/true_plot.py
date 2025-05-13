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
    """
    Version optimisée du filtrage des inversions par rapport aux TE.
    Utilise une approche vectorisée et un système d'indexation spatiale.
    """
    print(f"Tri des inversions : {methode} ...")
    inv_df = df_inv.copy()
    
    # Précalcul des longueurs d'inversion une seule fois
    inv_df['length'] = inv_df['end'] - inv_df['start']
    
    # Créer un dictionnaire pour un accès rapide aux TE par (chr, iso)
    te_by_chrom_iso = {}
    for (chr_val, iso_val), group in df_te.groupby(['chr', 'iso']):
        group_copy = group.copy()
        group_copy['te_length'] = group_copy['end'] - group_copy['start']
        te_by_chrom_iso[(chr_val, iso_val)] = group_copy
    
    # Identifier les inversions à filtrer plus efficacement
    to_remove = set()
    
    # Traiter par lots pour une meilleure efficacité mémoire
    for (chr_val, iso_val), inv_group in inv_df.groupby(['chr', 'iso']):
        # Récupérer les TE pour cette combinaison chr/iso
        if (chr_val, iso_val) not in te_by_chrom_iso: continue
            
        matching_tes = te_by_chrom_iso[(chr_val, iso_val)]
        
        # Prétraitement pour jointure spatiale efficace - créer des intervalles
        # Pour chaque inversion dans le groupe
        for idx, inv in inv_group.iterrows():
            # Filtrer les TE qui pourraient chevaucher cette inversion (filtre grossier)
            candidate_tes = matching_tes[
                (matching_tes['start'] <= inv['end']) & 
                (matching_tes['end'] >= inv['start'])
            ]
            
            if len(candidate_tes) == 0: continue
                
            # Calcul vectorisé des chevauchements
            overlap_starts = np.maximum(inv['start'], candidate_tes['start'].values)
            overlap_ends = np.minimum(inv['end'], candidate_tes['end'].values)
            overlap_lengths = np.maximum(0, overlap_ends - overlap_starts)
            
            # Calcul vectorisé des pourcentages de chevauchement
            inv_overlap_pcts = (overlap_lengths / inv['length']) * 100
            te_overlap_pcts = (overlap_lengths / candidate_tes['te_length'].values) * 100
            
            # Test de chevauchement bidirectionnel ≥ 95% 
            if ((inv_overlap_pcts >= 95) & (te_overlap_pcts >= 95)).any():
                to_remove.add(idx)
    
    # Filtrer les inversions à supprimer
    nb_filtered = len(to_remove)
    filtered_inv_df = inv_df.drop(to_remove)
    
    pourcentage = (nb_filtered / len(inv_df)) * 100 if len(inv_df) > 0 else 0
    print(f"{methode}: nombre d'inversions filtrées: {nb_filtered} ({pourcentage:.2f}%)")
    
    return filtered_inv_df


def Te_density(df_inv, df_te, window_size_half, range_5to3, full_df, chr_size, inv_max_size):
    """
    Version optimisée pour calculer la densité des éléments transposables autour des inversions.
    """
    # Précalculer les intervalles pour améliorer les performances
    te_tx_recover_5to3 = defaultdict(list)
    te_tx_recover_3to5 = defaultdict(list)
    range_3to5 = [-x for x in range_5to3[::-1]]

    if max(range_5to3) == inv_max_size // 2: range_3to5 = range_3to5[1::]
    
    # Prétraitement des TE - créer une structure de données efficace pour les requêtes spatiales
    # Créer deux arrays pour start et end pour un accès vectorisé
    te_starts = df_te['start'].values
    te_ends = df_te['end'].values
    
    
    for _, inv in df_inv.iterrows():
        mid_size_inv = (inv['end'] - inv['start']) // 2
        full_df_filetred = full_df[(full_df['start'] != inv['start']) & (full_df['end'] != inv['end'])].copy()
        
        # Traitement côté 5'→3'
        for range_val in range_5to3:

            # très important pour avoir toutes les range en clés ! 
            if range_val not in te_tx_recover_5to3:te_tx_recover_5to3[range_val] = []
            
            if range_val > mid_size_inv: continue
                
            pos = inv['start'] + range_val
            window_start = pos - window_size_half
            window_end = pos + window_size_half
        
            # Vérification des limites
            if window_start < 0: continue
            
            # Suppression des inversions chevauchantes
            mask = (full_df_filetred['start'] < window_end) & (full_df_filetred['end'] > window_start) 
            #if not mask.any() and range_val < 0 : continue


            # Sélection vectorisée des TE dans la fenêtre
            in_window = (te_starts <= window_end) & (te_ends >= window_start)
            if not any(in_window):
                te_tx_recover_5to3[range_val].append(0)  # Aucun TE dans cette fenêtre
                continue
                
            # Calculer les chevauchements vectorisés
            te_in_window_starts = te_starts[in_window]
            te_in_window_ends = te_ends[in_window]
            
            overlap_start = np.maximum(te_in_window_starts, window_start)
            overlap_end = np.minimum(te_in_window_ends, window_end)
            overlap_length = np.maximum(0, overlap_end - overlap_start)
            
            total_bases_in_te = np.sum(overlap_length)
            density = total_bases_in_te / (2 * window_size_half)
            te_tx_recover_5to3[range_val].append(density)
        
        # Traitement côté 3'→5' (similaire, optimisé de la même manière)
        for range_val in range_3to5:
            if range_val not in te_tx_recover_3to5:te_tx_recover_3to5[range_val] = []

            if -range_val > mid_size_inv: continue
                
            pos = inv['end'] + range_val
            window_start = pos - window_size_half
            window_end = pos + window_size_half
            
            # Vérification des limites
            if window_end > chr_size: continue
                
            # Suppression des inversions chevauchantes
            mask = (full_df_filetred['start'] < window_end) & (full_df_filetred['end'] > window_start) 
            #if not mask.any() and range_val > 0 : continue

            # Sélection vectorisée des TE dans la fenêtre
            in_window = (te_starts <= window_end) & (te_ends >= window_start)
            if not any(in_window):
                te_tx_recover_3to5[range_val].append(0)  # Aucun TE dans cette fenêtre
                continue
                
            # Calculer les chevauchements vectorisés
            te_in_window_starts = te_starts[in_window]
            te_in_window_ends = te_ends[in_window]
            
            overlap_start = np.maximum(te_in_window_starts, window_start)
            overlap_end = np.minimum(te_in_window_ends, window_end)
            overlap_length = np.maximum(0, overlap_end - overlap_start)
            
            total_bases_in_te = np.sum(overlap_length)
            density = total_bases_in_te / (2 * window_size_half)
            te_tx_recover_3to5[range_val].append(density)
    
    # Calcul des moyennes et comptage - utilisation de numpy pour plus d'efficacité
    te_tx_recover_5to3 = {key: [np.mean(values) if values else 0, len(values)] 
                           for key, values in te_tx_recover_5to3.items()}
    te_tx_recover_3to5 = {key: [np.mean(values) if values else 0, len(values)] 
                           for key, values in te_tx_recover_3to5.items()}
    
    # Formatage du résultat
    pos = sorted(te_tx_recover_5to3) + sorted(te_tx_recover_3to5)
    value = [te_tx_recover_5to3[keys] for keys in sorted(te_tx_recover_5to3)] + \
            [te_tx_recover_3to5[keys] for keys in sorted(te_tx_recover_3to5)]

    # Formatage des étiquettes de position
    idx_zero = [i for i, val in enumerate(pos) if val == 0]
    if len(idx_zero) >= 2:
        pos[idx_zero[0]], pos[idx_zero[1]] = 'inv_start', 'inv_end'
    
    if max(range_5to3) == inv_max_size // 2:
        idx_mid = [i for i, val in enumerate(pos) if val == inv_max_size // 2]
        if idx_mid: pos[idx_mid[0]] = 'mid_inv'

    # Identifier les positions de début et fin
    start_index = end_index = -1
    for i, item in enumerate(pos):
        if item == 'inv_start': start_index = i
        elif item == 'inv_end':
            end_index = i
            break
    
    # Renommer les positions flanquantes
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
    
    #print(f"Traitement de {methode} {iso} {chrom}...")
    
    # Appel à la fonction Te_density
    pos, value = Te_density(df_iso_chr_inv, te_data, window_size_half, 
        range_5to3, full_df, chr_size, inv_max_size)
    return methode, iso, chrom, pos, value


def make_visualization(all_result, output_path):
    """
    Fonction modifiée pour:
    1. Conserver le graphique global avec comptage d'inversions
    2. Regrouper les graphiques chromosomiques par méthode sur la même image
    3. Créer un graphique global par méthode montrant la densité par chromosome
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
                            isolat_data[iso][methode][pos].append((mean_density, count, chr))
                            # Données pour le graphique global
                            global_data[methode][pos].append((mean_density, count))
    
    # 1. Regrouper tous les graphiques par méthode pour tous les isolats
    for methode in global_data.keys():
        # Créer un dictionnaire pour stocker les données de tous les isolats pour cette méthode
        all_isolates_by_method = {}
        
        for iso in isolat_data.keys():
            if methode in isolat_data[iso]:
                # Ajouter les données de cet isolat au dictionnaire
                all_isolates_by_method[iso] = isolat_data[iso][methode]
        
        # Créer le graphique regroupé pour cette méthode
        create_grouped_isolates_plot(all_isolates_by_method, 
                                    f'Densité des TEs par chromosome - {methode} - Tous isolats',
                                    os.path.join(output_path, f'te_density_{methode}_all_isolats_by_chr.png'))
    
    # 2. Créer un graphique global par méthode montrant la densité par chromosome
    create_global_chr_comparison_plot(isolat_data, 
                                      'Densité des TEs par chromosome - Moyenne de tous les isolats',
                                      os.path.join(output_path, f'te_density_global_by_chr.png'))
    
    # Calculer le nombre total de TE pour chaque méthode
    te_count_by_method = {}
    for methode, pos_data in global_data.items():
        # Prendre une position de référence (par exemple, 'inv_start' si elle existe)
        ref_pos = next((p for p in pos_data.keys() if p == 'inv_start'), next(iter(pos_data.keys())))
        values_at_pos = pos_data[ref_pos]
        total_te_count = sum([v[1] for v in values_at_pos]) if values_at_pos else 0
        te_count_by_method[methode] = total_te_count
    
    # Créer le graphique global pour tous les isolats avec le comptage
    create_density_plot(global_data, 
                       'Densité moyenne des TEs autour des inversions - Tous isolats', 
                       os.path.join(output_path, f'sort_SF_te_density_all_isolats.png'),
                       te_count_by_method)


def create_grouped_isolates_plot(isolates_dict, title, output_file):
    """
    Nouvelle fonction pour créer un graphique regroupant tous les isolats 
    pour une méthode donnée, avec une sous-figure par isolat
    """
    # Calculer le nombre d'isolats et déterminer la disposition de la grille
    n_isolates = len(isolates_dict)
    if n_isolates == 0:
        print("Aucun isolat trouvé pour créer le graphique.")
        return
    
    # Calculer un layout optimal pour la grille
    n_cols = min(3, n_isolates)  # Maximum 3 colonnes
    n_rows = (n_isolates + n_cols - 1) // n_cols  # Arrondir vers le haut
    
    # Créer une figure avec des sous-graphiques
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6*n_cols, 4*n_rows), sharex=True, sharey=True)
    fig.suptitle(title, fontsize=16)
    
    # S'assurer que axes est toujours un tableau 2D même avec une seule ligne ou colonne
    if n_rows == 1 and n_cols == 1:
        axes = np.array([[axes]])
    elif n_rows == 1:
        axes = np.array([axes])
    elif n_cols == 1:
        axes = axes.reshape(-1, 1)
    
    # Variables pour stocker les informations sur les positions
    # Extraire les positions des données
    positions_ref = None
    for iso, data_dict in isolates_dict.items():
        positions_ref = list(data_dict.keys())
        if positions_ref:
            break
    
    if not positions_ref:
        print("Aucune position trouvée dans les données.")
        return
    
    # Chercher les indices de début et fin d'inversion
    start_index = positions_ref.index('inv_start') if 'inv_start' in positions_ref else None
    end_index = positions_ref.index('inv_end') if 'inv_end' in positions_ref else None
    
    # Tracer les données pour chaque isolat
    for idx, (iso, data_dict) in enumerate(isolates_dict.items()):
        row = idx // n_cols
        col = idx % n_cols
        ax = axes[row, col]
        
        # Organiser les données par chromosome
        chr_data = defaultdict(lambda: defaultdict(list))
        global_data = defaultdict(list)
        
        positions = positions_ref
        
        # Collecter les données par chromosome
        for pos, values_with_chr in data_dict.items():
            # Grouper par chromosome
            by_chr = defaultdict(list)
            for density, count, chr_name in values_with_chr:
                by_chr[chr_name].append((density, count))
                # Ajouter aussi aux données globales
                global_data[pos].append((density, count))
            
            # Calculer la moyenne pour chaque chromosome
            for chr_name, values in by_chr.items():
                mean_density = np.mean([v[0] for v in values]) if values else 0
                chr_data[chr_name][pos] = mean_density
        
        # Calculer la moyenne globale
        global_means = [np.mean([v[0] for v in global_data[pos]]) if global_data[pos] else 0 for pos in positions]
        
        # Tracer les lignes pour chaque chromosome
        for chr_name, pos_data in chr_data.items():
            densities = [pos_data[pos] if pos in pos_data else 0 for pos in positions]
            ax.plot(range(len(positions)), densities, alpha=0.3, linewidth=1, label=f"Chr {chr_name}")
        
        # Tracer la ligne globale (moyenne de tous les chromosomes) en gras
        ax.plot(range(len(positions)), global_means, color='black', linewidth=2, label="Moyenne")
        
        # Ajouter des lignes verticales pour les points d'inversion
        if start_index is not None: ax.axvline(x=start_index, color='red', linestyle='--', alpha=0.7)
        if end_index is not None: ax.axvline(x=end_index, color='red', linestyle='--', alpha=0.7)
        
        # Améliorer la lisibilité des étiquettes
        max_ticks = 20  # Nombre maximum d'étiquettes souhaité
        
        if start_index is not None and end_index is not None:
            # Calculer les espaces avant, pendant et après l'inversion
            space_before = start_index
            space_inversion = end_index - start_index
            space_after = len(positions) - end_index - 1
            
            # Calculer des pas homogènes pour chaque section
            step_before = max(1, space_before // (max_ticks // 3)) if space_before > 0 else 1
            step_inversion = max(1, space_inversion // (max_ticks // 3)) if space_inversion > 0 else 1
            step_after = max(1, space_after // (max_ticks // 3)) if space_after > 0 else 1
            
            # Créer les indices à afficher
            indices_before = list(range(0, start_index, step_before))
            if indices_before and indices_before[-1] != start_index - 1 and start_index - 1 >= 0:
                indices_before.append(start_index)  # Ajouter le point juste avant start_index
                
            indices_inversion = [start_index] + list(range(start_index + step_inversion, end_index, step_inversion)) + [end_index]
            
            indices_after = list(range(end_index + step_after, len(positions), step_after))
            if indices_after and indices_after[0] != end_index + 1 and end_index + 1 < len(positions):
                indices_after.insert(0, end_index)  # Ajouter le point juste après end_index
                
            indices_to_show = sorted(set(indices_before + indices_inversion + indices_after))
        else:
            # Fallback si start_index ou end_index n'est pas trouvé
            step = max(1, len(positions) // max_ticks)
            indices_to_show = list(range(0, len(positions), step))
        
        # Ajouter les étiquettes de l'axe X avec les noms des positions
        ax.set_xticks([i for i in indices_to_show])
        ax.set_xticklabels([positions[i] for i in indices_to_show], rotation=90)
        
        # Configurer le sous-graphique
        ax.set_title(f"Isolat {iso}")
        ax.set_ylim(0, 1)
        ax.grid(True, linestyle='--', alpha=0.5)
        
        # Ne pas ajouter de légende dans les sous-graphiques individuels
    # La légende commune sera ajoutée en dehors des sous-graphiques
    
    # Masquer les axes inutilisés
    for idx in range(n_isolates, n_rows * n_cols):
        row = idx // n_cols
        col = idx % n_cols
        axes[row, col].set_visible(False)
    
    # Ajouter une légende commune en bas de la figure
    # Récupérer tous les chromosomes uniques et ajouter la moyenne globale
    all_chromosomes = set()
    for iso_data in isolates_dict.values():
        for pos, values_with_chr in iso_data.items():
            for _, _, chr_name in values_with_chr:
                all_chromosomes.add(chr_name)
    
    # Créer des handles et labels pour la légende
    legend_handles = []
    legend_labels = []
    
    # Ajouter les chromosomes
    for chr_name in sorted(all_chromosomes):
        legend_handles.append(plt.Line2D([0], [0], color='C' + str(list(sorted(all_chromosomes)).index(chr_name) % 10), lw=1, alpha=0.3))
        legend_labels.append(f"Chr {chr_name}")
    
    # Ajouter la moyenne globale
    legend_handles.append(plt.Line2D([0], [0], color='black', lw=2))
    legend_labels.append("Moyenne globale")
    
    # Ajouter la légende en bas avec 5 éléments par ligne
    fig.legend(legend_handles, legend_labels, loc='lower center', 
               bbox_to_anchor=(0.5, 0), ncol=min(5, len(legend_handles)), 
               fontsize='small', frameon=True)
    
    # Ajouter des étiquettes communes pour les axes x et y
    fig.text(0.5, 0.08, 'Position relative', ha='center', va='center', fontsize=14)
    fig.text(0.06, 0.5, 'Densité moyenne des TEs', ha='center', va='center', rotation='vertical', fontsize=14)
    
    # Ajuster la mise en page pour laisser de l'espace pour la légende en bas
    plt.tight_layout(rect=[0.07, 0.15, 1, 0.95])
    
    # Sauvegarder le graphique
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Figure groupée sauvegardée dans {output_file}")


def create_global_chr_comparison_plot(all_isolates_data, title, output_file, positions_ref=None):
    """
    Nouvelle fonction pour créer un graphique global par méthode 
    montrant la densité par chromosome pour tous les isolats combinés
    """
    # Extraire toutes les méthodes disponibles
    all_methods = set()
    for iso_data in all_isolates_data.values():
        all_methods.update(iso_data.keys())
    
    # Calculer le nombre de méthodes et déterminer la disposition de la grille
    n_methods = len(all_methods)
    if n_methods == 0:
        print("Aucune méthode trouvée pour créer le graphique.")
        return
    
    # Calculer un layout optimal pour la grille
    n_cols = min(2, n_methods)  # Maximum 2 colonnes
    n_rows = (n_methods + n_cols - 1) // n_cols  # Arrondir vers le haut
    
    # Créer une figure avec des sous-graphiques
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(8*n_cols, 6*n_rows), sharex=True, sharey=True)
    fig.suptitle(title, fontsize=16)
    
    # S'assurer que axes est toujours un tableau 2D même avec une seule ligne ou colonne
    if n_rows == 1 and n_cols == 1:
        axes = np.array([[axes]])
    elif n_rows == 1:
        axes = np.array([axes])
    elif n_cols == 1:
        axes = axes.reshape(-1, 1)
    
    # Variables pour stocker les informations sur les positions
    start_index, end_index = None, None
    
    # Si positions_ref n'est pas fourni, on le détermine à partir des données
    if positions_ref is None:
        for iso, iso_data in all_isolates_data.items():
            for methode in iso_data:
                if iso_data[methode]:
                    positions_ref = list(iso_data[methode].keys())
                    break
            if positions_ref:
                break
    
    # Déterminer les indices de début et fin d'inversion
    if positions_ref:
        start_index = positions_ref.index('inv_start') if 'inv_start' in positions_ref else None
        end_index = positions_ref.index('inv_end') if 'inv_end' in positions_ref else None
    
    # Tracer les données pour chaque méthode
    for idx, methode in enumerate(sorted(all_methods)):
        row = idx // n_cols
        col = idx % n_cols
        ax = axes[row, col]
        
        # Dictionnaire pour stocker les données agrégées par chromosome pour cette méthode
        chr_data_combined = defaultdict(lambda: defaultdict(list))
        
        # Collecter et combiner les données de tous les isolats pour cette méthode
        for iso, iso_data in all_isolates_data.items():
            if methode in iso_data:
                for pos, values_with_chr in iso_data[methode].items():
                    # Regrouper par chromosome
                    for density, count, chr_name in values_with_chr:
                        chr_data_combined[chr_name][pos].append(density)
        
        # Si aucune donnée trouvée pour cette méthode, continuer
        if not chr_data_combined:
            ax.text(0.5, 0.5, f"Pas de données pour {methode}", ha='center', va='center', transform=ax.transAxes)
            continue
        
        positions = positions_ref if positions_ref else list(next(iter(chr_data_combined.values())).keys())
        
        # Calculer les moyennes par chromosome et globale
        chr_means = {}
        all_densities_by_pos = defaultdict(list)
        
        for chr_name, pos_data in chr_data_combined.items():
            chr_means[chr_name] = {}
            for pos in positions:
                if pos in pos_data and pos_data[pos]:
                    mean_density = np.mean(pos_data[pos])
                    chr_means[chr_name][pos] = mean_density
                    all_densities_by_pos[pos].append(mean_density)
                else:
                    chr_means[chr_name][pos] = 0
        
        # Calculer la moyenne globale (tous chromosomes confondus)
        global_means = [np.mean(all_densities_by_pos[pos]) if pos in all_densities_by_pos and all_densities_by_pos[pos] else 0 for pos in positions]
        
        # Tracer les lignes pour chaque chromosome
        for chr_name, pos_data in chr_means.items():
            densities = [pos_data[pos] if pos in pos_data else 0 for pos in positions]
            ax.plot(range(len(positions)), densities, alpha=0.3, linewidth=1, label=f"Chr {chr_name}")
        
        # Tracer la ligne globale (moyenne de tous les chromosomes) en gras
        ax.plot(range(len(positions)), global_means, color='black', linewidth=3, label="Moyenne globale")
        
        # Ajouter des lignes verticales pour les points d'inversion
        if start_index is not None: ax.axvline(x=start_index, color='red', linestyle='--', alpha=0.7)
        if end_index is not None: ax.axvline(x=end_index, color='red', linestyle='--', alpha=0.7)
        
        # Améliorer la lisibilité des étiquettes
        max_ticks = 20  # Nombre maximum d'étiquettes souhaité
        
        if start_index is not None and end_index is not None:
            # Calculer les espaces avant, pendant et après l'inversion
            space_before = start_index
            space_inversion = end_index - start_index
            space_after = len(positions) - end_index - 1
            
            # Calculer des pas homogènes pour chaque section
            step_before = max(1, space_before // (max_ticks // 3)) if space_before > 0 else 1
            step_inversion = max(1, space_inversion // (max_ticks // 3)) if space_inversion > 0 else 1
            step_after = max(1, space_after // (max_ticks // 3)) if space_after > 0 else 1
            
            # Créer les indices à afficher
            indices_before = list(range(0, start_index, step_before))
            if indices_before and indices_before[-1] != start_index - 1 and start_index - 1 >= 0:
                indices_before.append(start_index)  # Ajouter le point juste avant start_index
                
            indices_inversion = [start_index] + list(range(start_index + step_inversion, end_index, step_inversion)) + [end_index]
            
            indices_after = list(range(end_index + step_after, len(positions), step_after))
            if indices_after and indices_after[0] != end_index + 1 and end_index + 1 < len(positions):
                indices_after.insert(0, end_index)  # Ajouter le point juste après end_index
                
            indices_to_show = sorted(set(indices_before + indices_inversion + indices_after))
        else:
            # Fallback si start_index ou end_index n'est pas trouvé
            step = max(1, len(positions) // max_ticks)
            indices_to_show = list(range(0, len(positions), step))
            
        # Ajouter les étiquettes de l'axe X avec les noms des positions
        ax.set_xticks([i for i in indices_to_show])
        ax.set_xticklabels([positions[i] for i in indices_to_show], rotation=90)
        
        # Configurer le sous-graphique
        ax.set_title(f"Méthode: {methode}")
        ax.set_ylim(0, 1)
        ax.grid(True, linestyle='--', alpha=0.5)
        
        # Ajouter légende uniquement au premier graphique pour économiser de l'espace
        if row == 0 and col == 0:
            ax.legend(loc='upper right', fontsize='small')
    
    # Masquer les axes inutilisés
    for idx in range(n_methods, n_rows * n_cols):
        row = idx // n_cols
        col = idx % n_cols
        if idx < len(axes.flat):  # Vérifier que l'indice est valide
            axes[row, col].set_visible(False)
    
    # Ajouter une légende commune en bas de la figure
    # Récupérer tous les chromosomes uniques et ajouter la moyenne globale
    all_chromosomes = set()
    for chr_data in chr_data_combined.keys():
        all_chromosomes.add(chr_data)
    
    # Créer des handles et labels pour la légende
    legend_handles = []
    legend_labels = []
    
    # Ajouter les chromosomes
    for chr_name in sorted(all_chromosomes):
        legend_handles.append(plt.Line2D([0], [0], color='C' + str(list(sorted(all_chromosomes)).index(chr_name) % 10), lw=1, alpha=0.3))
        legend_labels.append(f"Chr {chr_name}")
    
    # Ajouter la moyenne globale
    legend_handles.append(plt.Line2D([0], [0], color='black', lw=3))
    legend_labels.append("Moyenne globale")
    
    # Ajouter la légende en bas avec 5 éléments par ligne
    fig.legend(legend_handles, legend_labels, loc='lower center', 
               bbox_to_anchor=(0.5, 0), ncol=min(5, len(legend_handles)), 
               fontsize='small', frameon=True)
    
    # Ajouter des étiquettes communes pour les axes
    fig.text(0.06, 0.5, 'Densité moyenne des TEs', ha='center', va='center', rotation='vertical', fontsize=14)
    
    # Ajuster la mise en page pour laisser de l'espace pour la légende en bas
    plt.tight_layout(rect=[0.07, 0.15, 1, 0.95])
    
    # Sauvegarder le graphique
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Figure globale par méthode sauvegardée dans {output_file}")


def create_chr_comparison_plot(data_dict, title, output_file):
    """
    Fonction pour créer un graphique comparant tous les chromosomes
    et la moyenne globale pour une méthode et un isolat
    """
    # Organiser les données par chromosome et calculer aussi les données globales
    chr_data = defaultdict(lambda: defaultdict(list))
    global_data = defaultdict(list)
    
    # Extraire les positions une seule fois
    positions = list(data_dict.keys())
    
    # Collecter les données par chromosome
    for pos, values_with_chr in data_dict.items():
        # Grouper par chromosome
        by_chr = defaultdict(list)
        for density, count, chr_name in values_with_chr:
            by_chr[chr_name].append((density, count))
            # Ajouter aussi aux données globales
            global_data[pos].append((density, count))
        
        # Calculer la moyenne pour chaque chromosome
        for chr_name, values in by_chr.items():
            mean_density = np.mean([v[0] for v in values]) if values else 0
            chr_data[chr_name][pos] = mean_density
    
    # Calculer la moyenne globale
    global_means = [np.mean([v[0] for v in global_data[pos]]) if global_data[pos] else 0 for pos in positions]
    
    # Créer la figure
    plt.figure(figsize=(16, 9))
    
    # Tracer les lignes pour chaque chromosome
    for chr_name, pos_data in chr_data.items():
        densities = [pos_data[pos] if pos in pos_data else 0 for pos in positions]
        plt.plot(range(len(positions)), densities, alpha=0.3, linewidth=1, label=f"Chr {chr_name}")
    
    # Tracer la ligne globale (moyenne de tous les chromosomes) en gras
    plt.plot(range(len(positions)), global_means, color='black', linewidth=3, label="Moyenne globale")
    
    # Ajouter des lignes verticales pour les points d'inversion
    start_index = positions.index('inv_start') if 'inv_start' in positions else None
    end_index = positions.index('inv_end') if 'inv_end' in positions else None
    
    if start_index is not None: plt.axvline(x=start_index, color='gray', linestyle='--', alpha=0.7)
    if end_index is not None: plt.axvline(x=end_index, color='gray', linestyle='--', alpha=0.7)
    
    # Améliorer la lisibilité des étiquettes
    max_ticks = 20  # Nombre maximum d'étiquettes souhaité
    
    if start_index is not None and end_index is not None:
        # Calculer les espaces avant, pendant et après l'inversion
        space_before = start_index
        space_inversion = end_index - start_index
        space_after = len(positions) - end_index - 1
        
        # Calculer des pas homogènes pour chaque section
        step_before = max(1, space_before // (max_ticks // 3)) if space_before > 0 else 1
        step_inversion = max(1, space_inversion // (max_ticks // 3)) if space_inversion > 0 else 1
        step_after = max(1, space_after // (max_ticks // 3)) if space_after > 0 else 1
        
        # Créer les indices à afficher
        indices_before = list(range(0, start_index, step_before))
        if indices_before and indices_before[-1] != start_index - 1 and start_index - 1 >= 0:
            indices_before.append(start_index)  # Ajouter le point juste avant start_index
            
        indices_inversion = [start_index] + list(range(start_index + step_inversion, end_index, step_inversion)) + [end_index]
        
        indices_after = list(range(end_index + step_after, len(positions), step_after))
        if indices_after and indices_after[0] != end_index + 1 and end_index + 1 < len(positions):
            indices_after.insert(0, end_index)  # Ajouter le point juste après end_index
            
        indices_to_show = sorted(set(indices_before + indices_inversion + indices_after))
    else:
        # Fallback si start_index ou end_index n'est pas trouvé
        step = max(1, len(positions) // max_ticks)
        indices_to_show = list(range(0, len(positions), step))
    
    plt.xticks([i for i in indices_to_show], [positions[i] for i in indices_to_show], rotation=90)
    
    plt.ylabel('Densité moyenne des TEs')
    plt.ylim(0, 1)
    plt.title(title)
    plt.grid(True, linestyle='--', alpha=0.7)
    
    # Ajouter légende avec taille réduite pour la lisibilité
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=5, fontsize='small')
    
    plt.tight_layout()
    
    # Sauvegarder le graphique
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Figure sauvegardée dans {output_file}")


def create_density_plot(data_dict, title, output_file, te_count_by_method=None):
    """
    Fonction pour créer un graphique de densité avec nombre de TE dans le titre
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
        # Ajouter le nombre de TE dans la légende si disponible
        label = methode
        if te_count_by_method and methode in te_count_by_method:
            label = f"{methode} (n={te_count_by_method[methode]})"
            
        ax1.plot(range(len(positions_combined)), densities, 
                label=label, marker='o' if methode == 'Cactus' else 's',
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
        # Ajouter le nombre de TE dans la légende si disponible
        label = methode
        if te_count_by_method and methode in te_count_by_method:
            label = f"{methode} (n={te_count_by_method[methode]})"
            
        ax2.bar(index + (i * bar_width - 0.5 * bar_width * (len(count_by_method) - 1)), 
               counts, 
               bar_width, 
               label=label)
    
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
                                    inv_min_size=0, inv_max_size=50000, flanking_region_size=20000, step_size=500,
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
    # Définition des fichiers d'entrée et de sortie
    input_files = {"Cactus_PtoMalgn_maf": "/home/hugob/stage_m2/06_Inversions_tools_comparaison/format_inv/output/inv_cactus_PtoMalgn_maf.tsv", 
                "Syri_minimap_PtoMalgn":"/home/hugob/stage_m2/06_Inversions_tools_comparaison/format_inv/output/inv_syri_minimap_PtoMalgn.tsv"}
    te_file = "/home/hugob/stage_m2/03_TE/output/all_TE.tsv"
    #te_file = "/home/hugob/stage_m2/04_gene/output/format_gene.tsv"
    chr_len_file = "/home/hugob/stage_m2/01_Processing/output/chr_len.tsv"
    output_path = "/home/hugob/stage_m2/07_TE_in_inversion/output2"

    # Filtres
    chromosomes_exclus = ['chr2003', 'chr2016', 'chr2017']
    isos = ['Gd_00293-aad']
    isos = ['Gd_00614-ba']
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

    # Exécuter l'analyse parallélisée
    results = visualize_te_inversion_relative_density_parallel(
        inv_datas, 
        te_data, 
        chr_size_dict, 
        output_path, 
        num_processes=mp.cpu_count() - 1  # Utiliser tous les cœurs sauf un
    )