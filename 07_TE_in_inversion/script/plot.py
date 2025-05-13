import pandas as pd
from collections import defaultdict
import numpy as np
import seaborn as sns
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import os
import multiprocessing as mp
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor
import psutil



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


def to_mean(data):
    for methode, dic_iso in data.items():
        for iso, dic_chr in dic_iso.items():
            for chr, dic_te_class in dic_chr.items():
                for te_class, strand_dic in dic_te_class.items():
                    for strand, pos_dic in strand_dic.items():
                        for pos, values in pos_dic.items():
                            if len(values) > 0: 
                                data[methode][iso][chr][te_class][strand][pos] = [np.mean(values), len(values)]
                            else: data[methode][iso][chr][te_class][strand][pos] = [0, 0]


    for methode, dic_iso in data.items():
        for iso, dic_chr in dic_iso.items():
            data[methode][iso]['globale'] = {}
            for chr, dic_te_class in dic_chr.items():
                if chr == 'globale' : continue
                df_chr_global = data[methode][iso][chr]['globale'] 
                for strand, pos_dic in df_chr_global.items():
                    if strand not in data[methode][iso]['globale'] : data[methode][iso]['globale'][strand] = {}
                    for pos, values in pos_dic.items():
                        if pos not in data[methode][iso]['globale'][strand] : data[methode][iso]['globale'][strand][pos] = [[], 0]
                        data[methode][iso]['globale'][strand][pos][0].append(values[0])
                        data[methode][iso]['globale'][strand][pos][1] += values[1]
            for strand, pos_dic in data[methode][iso]['globale'].items():   
                for pos, values in pos_dic.items():
                    if len(values[0]) > 0: 
                        data[methode][iso]['globale'][strand][pos] = [np.mean(values[0]), values[1]]

    for methode, dic_iso in data.items():
        data[methode]['globale'] = {}
        for iso, dic_chr in dic_iso.items():
            if iso == 'globale' : continue
            df_iso_global = data[methode][iso]['globale']
            for strand, pos_dic in df_iso_global.items():
                if strand not in data[methode]['globale'] : data[methode]['globale'][strand] = {}
                for pos, values in pos_dic.items():
                    if pos not in data[methode]['globale'][strand] : data[methode]['globale'][strand][pos] = [[], 0]
                    data[methode]['globale'][strand][pos][0].append(values[0])
                    data[methode]['globale'][strand][pos][1] += values[1]
        for strand, pos_dic in data[methode]['globale'].items():   
            for pos, values in pos_dic.items():
                if len(values[0]) > 0: 
                    data[methode]['globale'][strand][pos] = [np.mean(values[0]), values[1]]

    return data


def process_inversion(inv, full_df_filtered, te_starts, te_ends, range_5to3, range_3to5, window_size_half, mid_size_inv, chr_size):
    """Traite une seule inversion - fonction pour la parallélisation"""
    results = {'5to3': {}, '3to5': {}}
    
    # Traitement côté 5'→3'
    for range_val in range_5to3:
        if range_val > mid_size_inv: continue
            
        pos = inv['start'] + range_val
        window_start = pos - window_size_half
        window_end = pos + window_size_half
        
        # Vérification des limites
        if window_start < 0: continue
        
        # Vérification des inversions chevauchantes
        if range_val < 0:
            mask = (full_df_filtered['start'] < window_end) & (full_df_filtered['end'] > window_start)
            #if not mask.any(): continue
        
        # Sélection vectorisée des TE dans la fenêtre
        in_window = (te_starts <= window_end) & (te_ends >= window_start)
        if not np.any(in_window):
            results['5to3'][range_val] = 0  # Aucun TE dans cette fenêtre
            continue
        
        # Calcul optimisé des chevauchements
        te_in_window_starts = te_starts[in_window]
        te_in_window_ends = te_ends[in_window]
        
        overlap_start = np.maximum(te_in_window_starts, window_start)
        overlap_end = np.minimum(te_in_window_ends, window_end)
        overlap_length = np.maximum(0, overlap_end - overlap_start)
        
        total_bases_in_te = np.sum(overlap_length)
        density = total_bases_in_te / (2 * window_size_half)
        results['5to3'][range_val] = density
    
    # Traitement côté 3'→5'
    for range_val in range_3to5:
        if -range_val > mid_size_inv: continue
            
        pos = inv['end'] + range_val
        window_start = pos - window_size_half
        window_end = pos + window_size_half
        
        # Vérification des limites
        if window_end > chr_size: continue
        
        # Vérification des inversions chevauchantes
        if range_val > 0:
            mask = (full_df_filtered['start'] < window_end) & (full_df_filtered['end'] > window_start)
            #if not mask.any(): continue
        
        # Sélection vectorisée des TE dans la fenêtre
        in_window = (te_starts <= window_end) & (te_ends >= window_start)
        if not np.any(in_window):
            results['3to5'][range_val] = 0  # Aucun TE dans cette fenêtre
            continue
        
        # Calcul optimisé des chevauchements
        te_in_window_starts = te_starts[in_window]
        te_in_window_ends = te_ends[in_window]
        
        overlap_start = np.maximum(te_in_window_starts, window_start)
        overlap_end = np.minimum(te_in_window_ends, window_end)
        overlap_length = np.maximum(0, overlap_end - overlap_start)
        
        total_bases_in_te = np.sum(overlap_length)
        density = total_bases_in_te / (2 * window_size_half)
        results['3to5'][range_val] = density
    
    return results


def process_chromosome_te_class(args):
    """Traite un chromosome et une classe de TE"""
    methode, iso, chr, te_class, df_chr, df_te_class, range_5to3, range_3to5, window_size_half, chr_size = args
    
    result_dict = {'5to3': {}, '3to5': {}}
    
    # Initialisation des dictionnaires pour tous les ranges
    for range_val in range_5to3: result_dict['5to3'][range_val] = []
    for range_val in range_3to5: result_dict['3to5'][range_val] = []
    
    # Pré-calcul des valeurs pour éviter les recalculs répétés
    te_starts = df_te_class['start'].values
    te_ends = df_te_class['end'].values
    
    # Traitement de chaque inversion
    for _, inv in df_chr.iterrows():
        full_df_filtered = df_chr[(df_chr['start'] != inv['start']) & (df_chr['end'] != inv['end'])].copy()
        mid_size_inv = (inv['end'] - inv['start']) // 2
        
        # Traitement d'une inversion
        inv_results = process_inversion(inv, full_df_filtered, te_starts, te_ends, 
                                       range_5to3, range_3to5, window_size_half, mid_size_inv, chr_size)
        
        # Fusionner les résultats
        for direction in ['5to3', '3to5']:
            for range_val, density in inv_results[direction].items():
                result_dict[direction][range_val].append(density)
    

    return methode, iso, chr, te_class, result_dict


def make_data(inv_datas, te_data, chr_size_dict, inv_size=10000, flanking_region_size=5000, step_size=500, n_jobs=None):
    """
    Analyse optimisée et parallélisée de la densité de TE autour des inversions
    
    Parameters:
    -----------
    inv_datas : dict
        Dictionnaire contenant les données d'inversion par méthode
    te_data : pd.DataFrame
        DataFrame contenant les données d'éléments transposables
    chr_size_dict : dict
        Dictionnaire des tailles de chromosomes
    inv_size : int, default=10000
        Taille standard des inversions
    flanking_region_size : int, default=4000
        Taille des régions flanquantes à analyser
    step_size : int, default=1000
        Pas d'échantillonnage
    n_jobs : int, default=None
        Nombre de processus parallèles à utiliser. Si None, utilise (CPU count - 1)
        
    Returns:
    --------
    dict
        Dictionnaire structuré contenant les densités de TE
    """
    # Import de tqdm pour la barre de progression

    
    # Déterminer le nombre optimal de processus
    if n_jobs is None:
        n_jobs = max(1, psutil.cpu_count(logical=True) - 1)  # Physique - 1
    
    print(f"Utilisation de {n_jobs} processus parallèles")
    
    # Précalcul des intervalles une seule fois
    range_5to3 = [-dist for dist in reversed(range(0, flanking_region_size + 1, step_size))] + \
                 [dist for dist in range(step_size, inv_size // 2 + 1, step_size)]
    range_3to5 = [-x for x in range_5to3[::-1]]
    
    # Demi-taille de la fenêtre d'analyse
    window_size_half = step_size
        
    all_te_classes = te_data.copy()
    all_te_classes['class'] = "globale"
    te_data = pd.concat([te_data, all_te_classes], axis=0)
    
    # Créer une liste de toutes les tâches à traiter
    tasks = []
    full_data = {}
    for methode, df_inv in inv_datas.items():
        full_data[methode] = {}
        for iso, df in df_inv.groupby('iso'):
            full_data[methode][iso] = {}
            for chr, df_chr in df.groupby('chr'):
                full_data[methode][iso][chr] = {}
                df_te = te_data.query("iso == @iso and chr == @chr")
                chr_size = chr_size_dict[iso][chr]
                if df_te.empty: continue
                for te_class, df_te_class in df_te.groupby('class'):
                    # Ajouter la tâche à la liste
                    tasks.append((
                        methode, iso, chr, te_class, df_chr, df_te_class,
                        range_5to3, range_3to5, window_size_half, chr_size))
    
  
    # Initialiser le dictionnaire complet pour tous les niveaux
    for methode, df_inv in inv_datas.items():
        full_data[methode] = {}
        for iso, df in df_inv.groupby('iso'):
            full_data[methode][iso] = {}
            for chr, df_chr in df.groupby('chr'):
                full_data[methode][iso][chr] = {}
                df_te = te_data.query("iso == @iso and chr == @chr")
                if not df_te.empty:
                    for te_class in df_te['class'].unique():
                        full_data[methode][iso][chr][te_class] = {'5to3': {}, '3to5': {}}
    
    print(f"Nombre total de tâches à traiter : {len(tasks)}")
    
    # Exécution parallèle des tâches avec barre de progression
    results = []
    with ProcessPoolExecutor(max_workers=n_jobs) as executor:
        # Utilisation de tqdm pour afficher la progression
        for result in tqdm(executor.map(process_chromosome_te_class, tasks), 
                          total=len(tasks), 
                          desc="Traitement des données", 
                          unit="tâche"):
            results.append(result)
            
    # Traitement des résultats
    for methode, iso, chr, te_class, result_dict in results:
        for direction in ['5to3', '3to5']:
            for range_val, densities in result_dict[direction].items():
                full_data[methode][iso][chr][te_class][direction][range_val] = densities
                
    full_data = to_mean(full_data)
    
    return full_data


def reduire_intervalles(liste, x):
    if x < 3: return liste

    liste = sorted(liste)
    pas_de_base = liste[1] - liste[0]
    min_val = liste[0]
    max_val = liste[-1]
    meilleur_resultat = None
    meilleur_total = 0

    for facteur in range(1, (max_val - 0) // pas_de_base + 1):
        pas = facteur * pas_de_base
        # On ne garde le pas que si min_val et max_val sont alignés avec ce pas
        if (0 - min_val) % pas != 0 or (max_val - 0) % pas != 0: continue

        gauche = list(range(min_val + pas, 0, pas))
        droite = list(range(pas, max_val, pas))
        total = 1 + len(gauche) + 1 + len(droite) + 1  # min + gauche + 0 + droite + max
        if total <= x and total > meilleur_total:
            meilleur_resultat = [min_val] + gauche + [0] + droite + [max_val]
            meilleur_total = total

    if meilleur_resultat is None: return [min_val, 0, max_val]

    return meilleur_resultat


def plot_density(data, title, show_histogram=True):
    """
    Plot density data with optional histograms showing second values.
    
    Parameters:
    -----------
    data : dict
        Input data dictionary with nested structure
    title : str
        Title for the plot
    show_histogram : bool, default=False
        Whether to show histograms with second values below the plots
    """
    _5to3 = {}
    _3to5 = {}


    for main_key, subdict in data.items():
        for subkey, directions in subdict.items():
            for pos, value in directions.get('5to3', {}).items():
                _5to3.setdefault(pos, {})[main_key] = value
            for pos, value in directions.get('3to5', {}).items():
                _3to5.setdefault(pos, {})[main_key] = value

    df_5to3 = pd.DataFrame.from_dict(_5to3, orient='index').sort_index()
    df_3to5 = pd.DataFrame.from_dict(_3to5, orient='index').sort_index()



    def extract_value(df,i):
        df_values = df.copy()
        for col in df_values.columns:
            df_values[col] = df_values[col].apply(lambda x: x[i])
        return df_values
    

    # Extraction des valeurs
    df_5to3_values = extract_value(df_5to3.copy(),0)
    df_3to5_values = extract_value(df_3to5.copy(),0)
    
    x_ticks_5to3 = reduire_intervalles(sorted(df_5to3_values.index.tolist()), 30)
    x_tick_labels_5to3 = [str(int(x)) if x != 0 else 'start_inv' for x in x_ticks_5to3]
    x_ticks_3to5 = reduire_intervalles(sorted(df_3to5_values.index.tolist()), 30)
    x_tick_labels_3to5 = [str(int(x)) if x != 0 else 'start_inv' for x in x_ticks_3to5]

    # Création des graphiques
    if show_histogram:
        df_5to3_histo = extract_value(df_5to3.copy(),1)
        df_3to5_histo = extract_value(df_3to5.copy(),1)
        # Avec histogrammes (2 rangées)
        fig = plt.figure(figsize=(16, 12))
        gs = fig.add_gridspec(2, 2, height_ratios=[3, 1], wspace=0.1, hspace=0.1)
        ax1 = fig.add_subplot(gs[0, 0])  # Plot principal gauche
        ax2 = fig.add_subplot(gs[0, 1])  # Plot principal droit
        ax3 = fig.add_subplot(gs[1, 0])  # Histogramme gauche
        ax4 = fig.add_subplot(gs[1, 1])  # Histogramme droit
    else:
        # Sans histogrammes (1 rangée)
        fig = plt.figure(figsize=(16, 8))
        gs = fig.add_gridspec(1, 2, wspace=0.1)
        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[0, 1])

    # Ajout d'un titre global
    fig.suptitle("TROUVER TITRE", fontsize=16, y=0.98)

    
    
    
    
    
    
    
    # Première graphique - avant et début de l'inversion
    for column in df_5to3_values.columns:
        ax1.plot(df_5to3_values.index, df_5to3_values[column], marker='o', label=column, markersize=4)
    if not show_histogram : ax1.set_xlabel("Distance")
    ax1.set_ylabel("Valeur")
    ax1.set_ylim(0, 1)  
    ax1.spines['right'].set_visible(False)
    
    if show_histogram: x_tick_labels = ["" for x in x_tick_labels_5to3]
    else : x_tick_labels = x_tick_labels_5to3
    # Appliquer les ticks sélectionnés
    ax1.set_xticks(x_ticks_5to3)
    ax1.set_xticklabels(x_tick_labels, rotation=90)
    ax1.grid(True, color='gray', alpha=0.4)
    ax1.axvline(x=0, color='black', linestyle='--', alpha=0.8)
    ax1.set_title("Flanc gauche et début de l'inversion")

    # Deuxième graphique - fin de l'inversion et après
    for column in df_3to5_values.columns:
        ax2.plot(df_3to5_values.index, df_3to5_values[column], marker='o', label=column, markersize=4)
    ax2.set_ylim(0, 1)  
    if not show_histogram : ax2.set_xlabel("Distance")
    ax2.spines['left'].set_visible(False)
    ax2.tick_params(axis='y', which='both', length=0,labelleft=False)

    
    if show_histogram: x_tick_labels = ["" for x in x_tick_labels_3to5]
    else : x_tick_labels = x_tick_labels_3to5
    # Appliquer les ticks sélectionnés
    ax2.set_xticks(x_ticks_3to5)
    ax2.set_xticklabels(x_tick_labels, rotation=90)
    ax2.grid(True, color='gray', alpha=0.4)
    ax2.axvline(x=0, color='black', linestyle='--', alpha=0.8)
    ax2.set_title("Fin de l'inversion puis flanc droit")
    
    if show_histogram:
        bar_width = 0.8


        # Histogramme gauche (5'→3')
        x_5to3 = np.arange(len(df_5to3_histo.index))  # Positions numériques pour les barres
        x_mapping = dict(zip(df_5to3_histo.index, x_5to3))  # Crée un mapping entre index réels et positions numériques

        for i, column in enumerate(df_5to3_histo.columns):
            width = bar_width / len(df_5to3_histo.columns)
            offset = (i - (len(df_5to3_histo.columns) - 1) / 2) * width
            values = df_5to3_histo[column].values
            ax3.bar(x_5to3 + offset, values, width=width, label=column, alpha=0.7)

        ax3.set_xlabel("Distance")
        ax3.set_ylabel("Nombre d'inversions")
        ax3.grid(True,axis='y', color='gray', alpha=0.4)
        ax3.spines['right'].set_visible(False)
        if 0 in x_mapping: ax3.axvline(x=x_mapping[0], color='black', linestyle='--', alpha=0.8)
        # Convertissez vos ticks personnalisés en positions numériques pour qu'ils s'alignent avec les barres
        mapped_ticks = [x_mapping.get(tick, tick) for tick in x_ticks_5to3 if tick in x_mapping]
        ax3.set_xticks(mapped_ticks)
        ax3.set_xticklabels([str(int(x)) if x != 0 else 'start_inv' for x in x_ticks_5to3 if x in df_5to3_histo.index], rotation=90) 
            
            
    
    


        
        # Histogramme droit (3'→5')
        x_3to5 = np.arange(len(df_3to5_histo.index))  # Positions numériques pour les barres
        x_mapping = dict(zip(df_3to5_histo.index, x_3to5))  # Crée un mapping entre index réels et positions numériques

        for i, column in enumerate(df_3to5_histo.columns):
            width = bar_width / len(df_3to5_histo.columns)
            offset = (i - (len(df_3to5_histo.columns) - 1) / 2) * width
            values = df_3to5_histo[column].values
            ax4.bar(x_3to5 + offset, values, width=width, label=column, alpha=0.7)

        ax4.set_xlabel("Distance")
        ax4.tick_params(axis='y', which='both', length=0,labelleft=False)
        ax4.grid(True, axis='y', color='gray', alpha=0.4)
        ax4.spines['left'].set_visible(False)
        if 0 in x_mapping: ax4.axvline(x=x_mapping[0], color='black', linestyle='--', alpha=0.8)
        # Convertissez vos ticks personnalisés en positions numériques pour qu'ils s'alignent avec les barres
        mapped_ticks = [x_mapping.get(tick, tick) for tick in x_ticks_3to5 if tick in x_mapping]
        ax4.set_xticks(mapped_ticks)
        ax4.set_xticklabels([str(int(x)) if x != 0 else 'end_inv' for x in x_ticks_3to5 if x in df_3to5_histo.index], rotation=90)









    # Création d'une légende commune
    lines, labels = [], []
    for ax in [ax1, ax2]:
        line, label = ax.get_legend_handles_labels()
        lines.extend(line)
        labels.extend(label)
    
    # Supprimer les légendes individuelles
    ax1.get_legend().remove() if ax1.get_legend() else None
    ax2.get_legend().remove() if ax2.get_legend() else None
    
    # La légende est placée en bas de la figure
    by_label = dict(zip(labels, lines))  # Éliminer les doublons
    fig.legend(by_label.values(), by_label.keys(), loc='lower center', ncol=len(by_label), bbox_to_anchor=(0.5, 0.02))
    plt.subplots_adjust(bottom=0.15, left=0.05, right=0.95, top=0.9)  
    
    # Enregistrement de l'image
    plt.savefig(f"/home/hugob/stage_m2/07_TE_in_inversion/out_plot/{title}.png", dpi=300)
    print(f"L'image a été enregistrée sous le nom {title}.png")
    
    return fig


def test(data):
    pass

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
    isos = ['Gd_00293-aad', 'Gd_00045-a2ab']
    #isos = ['Gd_00614-ba']
    mapons = ['all']
    inv_datas = {}


    chrs = ['chr1001', 'chr1002']
    # Lecture des données TE
    te_data = pd.read_csv(te_file, sep='\t')
    te_data['iso'] = te_data['iso'].apply(lambda x: '_'.join(x.split('_')[:2]))
    te_data = te_data.query("chr not in @chromosomes_exclus").sort_values(by='start')
    te_data = te_data.query("chr in @chrs and iso in @isos").sort_values(by='start')



    # Lecture des fichiers d'inversions
    for methode, patt_file in input_files.items():
        methode_name = methode.split('_')[0]
        df = pd.read_csv(patt_file, sep='\t')
        df = df.query("chr not in @chromosomes_exclus and mapon in @mapons").sort_values(by='start')
        df = df.query("chr in @chrs and iso in @isos").sort_values(by='start')
        #inv_datas[methode_name] = del_invTE(df,te_data,methode_name) 
        inv_datas[methode_name] = df


    # Lecture des longueurs de chromosomes
    chr_len_df = pd.read_csv(chr_len_file, sep='\t')
    chr_len_df = chr_len_df.query("chr not in @chromosomes_exclus")    
    chr_size_dict = defaultdict(lambda: defaultdict(dict))
    for iso, chr, size in zip(chr_len_df['iso'], chr_len_df['chr'], chr_len_df['size']):
        chr_size_dict[iso][chr] = size

    full_data = make_data(inv_datas, te_data, chr_size_dict)

    plot_density(full_data, "inversion_analysis")
    test(full_data)