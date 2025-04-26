import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import math
from collections import defaultdict
import re
import seaborn as sns
from matplotlib.ticker import FuncFormatter
from matplotlib.gridspec import GridSpec
from itertools import combinations
import venn_tools


# ================== Foncions de formatage

def print_data(d, indent=0):
    for key, value in d.items():
        if isinstance(value, dict):
            print("  " * indent + str(key) + "/")
            print_data(value, indent + 1)
        else:
            print("  " * indent + f"{key} -> {len(value)} inversions")

def make_data(inv_files: dict):
    all_data = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(list))))

    for methode, file_path in inv_files.items():
        try:
            df = pd.read_csv(file_path, sep='\t')
            chromosomes_à_exclure = ['chr2003', 'chr2016', 'chr2017']
            df = df[~df['chr'].isin(chromosomes_à_exclure)]
            if 'mapon' not in df.columns: df['mapon'] = 'all'

            for _, row in df.iterrows():
                all_data[row['iso']][methode][row['chr']][row['mapon']].append((row['start'], row['end']))

        except Exception as e:
            print(f"Erreur lors de la lecture de {file_path} ({methode}): {str(e)}")

    return all_data

def merge_inversions(data, method, mapon, max_gap):
    """
    Fusionne les inversions directement à partir de la structure de données imbriquée
    pour une méthode et un mapon donnés.
    
    Args:
        data: dictionnaire imbriqué {iso -> {method -> {chr -> {mapon -> list[(start, end)]}}}}
        method: la méthode d'analyse à considérer
        mapon: le mapon à considérer
        max_gap: seuil d'écart maximum pour la fusion des inversions
        
    Returns:
        list: Liste des inversions fusionnées sous forme de tuples (start, end)
    """
    # Collecter toutes les inversions avec leurs métadonnées
    all_inversions = []
    
    for iso in data:
        if method in data[iso]:
            for chr_ in data[iso][method]:
                if mapon in data[iso][method][chr_]:
                    inversions = data[iso][method][chr_][mapon]
                    for start, end in inversions:
                        all_inversions.append((iso, chr_, start, end))
    
    # Si aucune inversion, retourner une liste vide
    if not all_inversions:
        return []
    
    # Trier les inversions par iso, chr, puis start
    all_inversions.sort(key=lambda x: (x[0], x[1], x[2]))
    
    # Fusionner les inversions par groupes (iso, chr)
    merged_inversions = []
    current_group = None
    current_start, current_end = None, None
    
    for inv in all_inversions:
        iso, chr_, start, end = inv
        group = (iso, chr_)
        
        # Premier élément ou nouveau groupe
        if current_group is None or current_group != group:
            if current_group is not None:
                merged_inversions.append((current_start, current_end))
            
            current_group = group
            current_start, current_end = start, end
            continue
        
        # Vérifier si on peut fusionner avec l'inversion précédente
        if start - current_end <= max_gap:
            current_end = max(current_end, end)
        else:
            merged_inversions.append((current_start, current_end))
            current_start, current_end = start, end
    
    # Ajouter la dernière inversion
    if current_group is not None:
        merged_inversions.append((current_start, current_end))
    
    return merged_inversions

def calculate_all_inversions_intersections(all_data, labels, genome_sizes, mapon="all"):
    """
    Calcule toutes les intersections de bases inversées pour tous les isolats
    en incluant le pourcentage par rapport à la taille du génome.
    
    Args:
        all_data: Dictionnaire structuré {iso -> {method -> {chr -> {mapon -> [(start, end)]}}}}
        labels: Liste des méthodes à considérer
        genome_sizes: Dictionnaire avec les tailles des génomes par isolat
        mapon: Le mapon à considérer (par défaut "default")
    
    Returns:
        Un dictionnaire {iso -> {combo: formatted_count}}
    """    
    results = {}
    
    method_to_index = {method: i for i, method in enumerate(sorted(labels))}
    
    # Pour chaque isolat
    for iso, methods_data in all_data.items():
        iso_results = {}
        
        # Préparer les ensembles de bases pour chaque méthode
        method_bases = {i: set() for i in range(len(labels))}
        
        # Remplir les ensembles pour chaque méthode, chromosome et mapon
        for method, chr_data in methods_data.items():
            if method not in method_to_index:
                continue
                
            method_idx = method_to_index[method]
            
            for chr_name, mapon_data in chr_data.items():
                if mapon in mapon_data:
                    regions = mapon_data[mapon]
                    for start, end in regions:
                        method_bases[method_idx].update(range(start, end + 1))
        
        # Calculer et stocker le nombre de bases pour chaque méthode et les intersections
        for r in range(1, len(labels) + 1):
            for combo in combinations(range(len(labels)), r):
                if r == 1:
                    # Méthode unique
                    count = len(method_bases[combo[0]])
                else:
                    # Intersection de plusieurs méthodes
                    intersection = set.intersection(*(method_bases[i] for i in combo))
                    count = len(intersection)
                
                # Calculer et formater le pourcentage
                if iso in genome_sizes and genome_sizes[iso]:
                    percentage = (count / genome_sizes[iso]) * 100
                    formatted_count = f"{percentage:.2f}%"
                else:
                    formatted_count = f"{count} (% inconnu)"
                
                # Garder combo comme clé (tuple d'indices)
                iso_results[combo] = formatted_count
                
        results[iso] = iso_results
    
    return results


# ================== PLOT 1
def plot_inv_by_mapon(df:dict,iso,chr_size_dict,lineage_dict, output_dir):
    

    for methode, chr_dic in df.items():
        # Couleur unique pour les inversions mappées
        inversion_color = 'green'

        n_chrs = len(chr_dic)
        n_rows = math.ceil(n_chrs / 3)  # 3 colonnes, donc diviser par 2 pour obtenir le nombre de lignes
        
        fig, axes = plt.subplots(n_rows, 3, figsize=(27, 4*n_rows))
        
        # Si un seul chromosome, ajuster les axes
        if n_chrs == 1:
            axes = np.array([[axes]])
        elif n_rows == 1:
            axes = axes.reshape(1, -1)
        
        # Parcourir les chromosomes
        i = 0
        for chr_name, mapon_dic in chr_dic.items():
            row_idx = i // 3
            col_idx = i % 3
            ax = axes[row_idx, col_idx]
            i+=1

            # Utiliser la taille du chromosome depuis le dic si disponible, sinon calculer
            if chr_size_dict[iso][chr_name]: max_pos =chr_size_dict[iso][chr_name]
            else: max_pos = max([end for mapon_invs in mapon_dic.values() for start, end in mapon_invs])

            # Trier les mapons par lignée
            sorted_mapons = []
            mapon_lineages = {}
            
            # Récupérer la lignée pour chaque mapon
            for mapon in mapon_dic.keys():
                lineage = lineage_dict.get(mapon, mapon)
                mapon_lineages[mapon] = lineage
                
            # Trier les mapons d'abord par lignée, puis par nom
            sorted_mapons = sorted(mapon_dic.keys(), key=lambda m: (mapon_lineages.get(m, mapon), m))
            sorted_mapons = sorted_mapons [::-1]

                # Obtenir palette de couleurs pour lignées
            unique_lineage = sorted(list(set(mapon_lineages.values())))
            lineage_palette = dict(zip(unique_lineage, sns.color_palette("tab20", len(unique_lineage))))
            
            # Une ligne par mapon
            for j, mapon in enumerate(sorted_mapons):
                y_pos = j
                inversions = mapon_dic[mapon]
                
                # Tracer les inversions (toutes de la même couleur)
                for start, end in inversions:
                    ax.plot([start, end], [y_pos, y_pos], linewidth=2, color=inversion_color)
            
            # Configuration de l'axe
            ax.set_title(f"{iso} - {chr_name}")
            ax.set_xlim(0, max_pos)
            ax.set_ylim(-0.5, len(sorted_mapons) - 0.5)
            ax.set_yticks(range(len(sorted_mapons)))
            ax.set_yticklabels(sorted_mapons)
            
            # Colorier chaque label y individuellement par lignée
            for j, mapon in enumerate(sorted_mapons):
                lineage = mapon_lineages.get(mapon, mapon)
                color = lineage_palette.get(lineage, 'black')
                ax.get_yticklabels()[j].set_color(color)
            
            ax.grid(True, linestyle='--', alpha=0.7)
            ax.set_xlabel('Position (bp)')
        
        # Cacher les axes non utilisés dans la dernière ligne si nombre impair de chromosomes
        if n_chrs % 2 != 0:
            axes[n_rows-1, 1].axis('off')
        
        plt.tight_layout()
        plt.savefig(f"{output_dir}/{iso}_{methode}.png", dpi=600, bbox_inches='tight')
        plt.close(fig)

# ================== PLOT 2
def plot_inversion_by_mapon(data, chr_size_dict, mapon_values, output_dir):
    """
    Crée une visualisation des inversions chromosomiques avec tous les chromosomes par isolat,
    et deux isolats par ligne. L'organisation s'adapte aux chromosomes disponibles pour chaque isolat.
    
    Args:
        data (dict): Dictionnaire structuré comme data[iso][methode][chr][mapon]
        chr_size_dict (dict): Dictionnaire contenant les tailles des chromosomes
        mapon_values (list): Liste des valeurs de mapon à traiter
        output_dir (str): Répertoire de sortie pour les figures
    """
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import numpy as np
    import os
    import math
    from matplotlib.ticker import FuncFormatter
    from matplotlib.gridspec import GridSpec
    
    plt.style.use('seaborn-v0_8-whitegrid')
    mpl.rcParams['font.family'] = 'sans-serif'
    mpl.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
    mpl.rcParams['axes.labelsize'] = 12
    mpl.rcParams['axes.titlesize'] = 14
    mpl.rcParams['xtick.labelsize'] = 10
    mpl.rcParams['ytick.labelsize'] = 10
    
    # Pour chaque mapon, créer une figure avec tous les isolats
    for mapon in mapon_values:
        # Collecte de toutes les méthodes et tous les chromosomes disponibles par isolat
        all_methodes = set()
        isolat_chromosomes = {}  # Chromosomes par isolat
        isolat_methodes = {}     # Méthodes disponibles par isolat
        all_isolats = list(data.keys())
        
        # Collecter les chromosomes et méthodes disponibles pour chaque isolat
        for iso in all_isolats:
            isolat_chromosomes[iso] = set()
            isolat_methodes[iso] = set()
            if iso in data:
                for methode, chr_dic in data[iso].items():
                    all_methodes.add(methode)
                    isolat_methodes[iso].add(methode)
                    for chr_name in chr_dic.keys():
                        isolat_chromosomes[iso].add(chr_name)
        
        # Tri des méthodes et isolats
        all_methodes = sorted(list(all_methodes))
        all_isolats = sorted(all_isolats)
        
        # Disposition avec 2 isolats par ligne
        n_isolats = len(all_isolats)
        n_cols = min(2, n_isolats)
        n_rows = math.ceil(n_isolats / n_cols)
                
        # Créer une grande figure
        fig = plt.figure(figsize=(20, 5 * n_rows), facecolor='white')
        
        
        # Générer les couleurs pour les méthodes
        cmap = plt.get_cmap("tab10")
        method_colors = {
            methode: cmap(i % 10) 
            for i, methode in enumerate(all_methodes)
        }
        
        # Définir un trieur naturel pour les chromosomes
        def natural_sort_key(s):
            import re
            return [int(text) if text.isdigit() else text.lower() 
                    for text in re.split(r'(\d+)', s)]
        
        # Pour chaque isolat, créer un sous-plot avec tous ses chromosomes
        for iso_idx, iso in enumerate(all_isolats):
            row = iso_idx // n_cols
            col = iso_idx % n_cols


            
            # Trier les chromosomes de cet isolat
            iso_chromosomes = sorted(list(isolat_chromosomes[iso]), key=natural_sort_key)
            if not iso_chromosomes:
                continue  # Passer si pas de chromosomes pour cet isolat
            
            # Méthodes disponibles pour cet isolat
            iso_methodes = sorted(list(isolat_methodes[iso]))
            # Calcul des positions pour la grille
            vertical_spacing = 0.15  # Espace vertical entre les rangées
            left_margin = 0.1 + col * 0.5
            right_margin = left_margin + 0.45
            bottom_margin = 0.1 + (n_rows - row - 1) * (1/n_rows) * (1 + vertical_spacing)
            top_margin = bottom_margin + (1/n_rows) * 0.85
            
            # Créer un sous-subplot pour cet isolat avec ses chromosomes et légende
            n_rows_subplot = len(iso_chromosomes)
            height_ratios = [1] * n_rows_subplot

            # Réserver un espace pour les chromosomes et la légende
            gs = GridSpec(n_rows_subplot, 1,
                        left=left_margin, right=right_margin,
                        bottom=bottom_margin, top=top_margin,
                        height_ratios=height_ratios,
                        hspace=0.05)
            
            # Trouver la taille maximale des chromosomes pour cet isolat
            max_chr_size = 0
            for chr_name in iso_chromosomes:
                if chr_name in chr_size_dict.get(iso, {}):
                    max_chr_size = max(max_chr_size, chr_size_dict[iso][chr_name])
            
            # Ajouter un titre pour cet isolat
            plt.figtext(0.1 + col * 0.5 + 0.25, 
                      0.11 + (n_rows - row) * (1/n_rows) - 0.03, 
                      f"{iso}", 
                      ha="center", fontsize=14, fontweight='bold')
            
            # Pour chaque chromosome de cet isolat
            for chr_idx, chr_name in enumerate(iso_chromosomes):
                ax = fig.add_subplot(gs[chr_idx, 0])
                
                # Configuration de base de l'axe
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.spines['left'].set_visible(False)
                
                # Ajouter le nom du chromosome à gauche
                ax.text(-0.01 * max_chr_size, 0.5, chr_name, 
                       transform=ax.transData, 
                       ha='right', va='center', 
                       fontweight='bold',
                       fontsize=11,
                       color='#333333')
                
                # Configuration de l'axe
                ax.set_ylim(0, 1)
                ax.set_yticks([])
                ax.grid(True, axis='x', linestyle='-', alpha=0.2, color='#cccccc')
                ax.axhline(y=0, color='#aaaaaa', linestyle='-', alpha=0.3, linewidth=0.8)
                ax.axhline(y=1, color='#aaaaaa', linestyle='-', alpha=0.3, linewidth=0.8)
                
                # Adapter dynamiquement le placement des méthodes selon leur nombre
                # Pour un positionnement plus centré
                n_methods = len(iso_methodes)
                
                if n_methods > 0:
                    # Centrer les méthodes dans l'axe Y
                    # Utiliser 70% de l'espace vertical pour les méthodes
                    method_height = min(0.7 / n_methods, 0.15)  # Limiter la hauteur maximale
                    spacing = min((0.7 - method_height * n_methods) / (n_methods + 1), 0.05)
                    
                    y_start = 0.15  # Commencer à 15% du bas
                    
                    # Calculer les positions y des méthodes de manière centrée
                    y_heights = {}
                    for idx, methode in enumerate(iso_methodes):
                        y_pos = y_start + idx * (method_height + spacing)
                        y_heights[methode] = y_pos + method_height / 2
                    
                    # Tracer les inversions
                    if iso in data:
                        for methode in iso_methodes:
                            if (methode in data[iso] and 
                                chr_name in data[iso][methode] and 
                                mapon in data[iso][methode][chr_name]):
                                
                                y_pos = y_heights[methode]
                                height = method_height * 0.8  # Légèrement plus petit que l'espace alloué
                                
                                # Tracer chaque inversion
                                for pos in data[iso][methode][chr_name][mapon]:
                                    start, end = pos[0], pos[1]
                                    rect = plt.Rectangle((start, y_pos - height/2), 
                                                      end - start, height, 
                                                      facecolor=method_colors[methode],
                                                      edgecolor='none',
                                                      alpha=0.8,
                                                      zorder=10)
                                    ax.add_patch(rect)
                
                # Formatter pour Mpb
                def format_x_tick(x, pos):
                    return '{:.1f}'.format(x / 1e6)
                
                # Configurer l'axe X seulement pour le dernier chromosome
                if chr_idx == len(iso_chromosomes) - 1:
                    ax.xaxis.set_major_formatter(FuncFormatter(format_x_tick))
                    ax.set_xlabel('Position (Mpb)', fontweight='bold', color='#333333', labelpad=10)
                else:
                    ax.set_xticklabels([])
                
                # Grille et ticks
                major_ticks = np.arange(0, max_chr_size, 1e6)  # Tous les 1 Mpb
                ax.set_xticks(major_ticks)
                ax.grid(True, which='major', axis='x', linestyle='-', alpha=0.2, color='#dddddd')
                
                # Limites des axes X
                ax.set_xlim(-0.02 * max_chr_size, 1.05 * max_chr_size)
            
            # Ajouter une légende spécifique à cet isolat en bas à droite
            if iso_methodes:
                legend_ax = fig.add_subplot(gs[-1, 0])
                legend_ax.axis('off')  # Cacher les axes
                
                # Créer des handles uniquement pour les méthodes présentes dans cet isolat
                handles = [plt.Rectangle((0, 0), 1, 1, color=method_colors[methode], alpha=0.8) 
                          for methode in iso_methodes]
                
                # Positionner la légende en bas à droite
                legend = fig.legend(handles, iso_methodes, 
                                  bbox_to_anchor=(right_margin, bottom_margin + 0.06),
                                  loc='lower right',
                                  frameon=True, 
                                  framealpha=1,
                                  fontsize=9)
        
        # Enregistrer le graphique
        output_file = os.path.join(output_dir, f"mapon_{mapon}_all_isolats.png")
        plt.savefig(output_file, dpi=200, bbox_inches='tight', facecolor='white')
        plt.close(fig)

# ================== PLOT 3 et 4
def plot_nb_size_by_mapon(data, thresholds, out_dir):
    """
    Crée des graphiques pour analyser les inversions par mapon et méthode.
    
    Args:
        data: dictionnaire imbriqué {iso -> {method -> {chr -> {mapon -> list[(start, end)]}}}}
        thresholds: liste des seuils à tester
        out_dir: répertoire de sortie pour les figures
    """
    # Pour chaque mapon, on va créer une image
    all_mapons = set()
    all_methods = set()
    
    # Identifier tous les mapons et méthodes disponibles
    for iso in data:
        for method in data[iso]:
            all_methods.add(method)
            for chr_ in data[iso][method]:
                for mapon in data[iso][method][chr_]:
                    all_mapons.add(mapon)
    
    all_methods = sorted(list(all_methods))
    all_mapons = sorted(list(all_mapons))
    
    # Pour chaque mapon, créer une figure avec des subplots pour chaque méthode
    for mapon in all_mapons:
        # Créer une figure avec 2 rangées (nb inversions, taille) et assez de colonnes pour chaque méthode
        fig, axes = plt.subplots(2, len(all_methods), figsize=(6*len(all_methods), 10))
        
        # Gérer le cas où il n'y a qu'une seule méthode
        if len(all_methods) == 1:
            axes = axes.reshape(2, 1)
        
        # Pour chaque méthode
        for m_idx, method in enumerate(all_methods):
            # Collecter et analyser les résultats pour différents seuils
            results = []
            
            for threshold in thresholds:
                # Fusionner les inversions directement
                merged_inversions = merge_inversions(data, method, mapon, threshold)
                
                # Si des inversions sont trouvées
                if merged_inversions:
                    # Calculer les statistiques
                    sizes = [end - start + 1 for start, end in merged_inversions]
                    nb_inversions = len(sizes)
                    mean_size = sum(sizes) / nb_inversions if nb_inversions > 0 else 0
                    
                    # Calculer la médiane manuellement
                    sizes.sort()
                    if nb_inversions % 2 == 0:
                        median_size = (sizes[nb_inversions//2 - 1] + sizes[nb_inversions//2]) / 2
                    else:
                        median_size = sizes[nb_inversions//2]
                    
                    # Stocker les résultats
                    results.append({
                        'threshold': threshold,
                        'nb_inversions': nb_inversions,
                        'mean_size': mean_size,
                        'median_size': median_size,
                        'min_size': min(sizes) if sizes else 0,
                        'max_size': max(sizes) if sizes else 0,
                        'total_size': sum(sizes)
                    })
                    
                    if threshold == 0 or threshold == 50:
                        mean_size_rounded = round(mean_size, 2)
                        median_size_rounded = round(median_size, 2)
                        print(f"Mapon {mapon}, Method {method}, Threshold {threshold}: nb = {nb_inversions}; "
                              f"mean_size = {mean_size_rounded}; median_size = {median_size_rounded}")
                else:
                    # Aucune inversion pour cette combinaison
                    results.append({
                        'threshold': threshold,
                        'nb_inversions': 0,
                        'mean_size': 0,
                        'median_size': 0,
                        'min_size': 0,
                        'max_size': 0,
                        'total_size': 0
                    })
            
            # Créer un dataframe à partir des résultats
            results_df = pd.DataFrame(results)
            
            if not results_df.empty and results_df['nb_inversions'].sum() > 0:
                # Graphique 1: Nombre d'inversions en fonction du seuil
                axes[0, m_idx].plot(results_df['threshold'], results_df['nb_inversions'], marker='o', linestyle='-')
                axes[0, m_idx].set_title(f'Nombre d\'inversions - {method}')
                axes[0, m_idx].set_xlabel('Seuil d\'écart minimum (bases)')
                axes[0, m_idx].set_ylabel('Nombre d\'inversions')
                axes[0, m_idx].grid(True)
                
                # Graphique 2: Taille moyenne et médiane des inversions en fonction du seuil
                axes[1, m_idx].plot(results_df['threshold'], results_df['mean_size'], marker='o', linestyle='-', label='Moyenne')
                axes[1, m_idx].plot(results_df['threshold'], results_df['median_size'], marker='s', linestyle='--', label='Médiane')
                axes[1, m_idx].set_title(f'Taille des inversions - {method}')
                axes[1, m_idx].set_xlabel('Seuil d\'écart minimum (bases)')
                axes[1, m_idx].set_ylabel('Taille (bases)')
                axes[1, m_idx].legend()
                axes[1, m_idx].grid(True)
            else:
                # Cas où aucune inversion n'est trouvée pour cette méthode/mapon
                axes[0, m_idx].text(0.5, 0.5, 'Aucune inversion', ha='center', va='center')
                axes[1, m_idx].text(0.5, 0.5, 'Aucune inversion', ha='center', va='center')
                
        # Titre global pour cette figure (mapon)
        fig.suptitle(f'Analyse des inversions pour le mapon: {mapon}', fontsize=16)
        plt.tight_layout(rect=[0, 0, 1, 0.96])  # Ajuster pour le titre
        
        # Sauvegarder la figure
        out_path = os.path.join(out_dir, f'inversions_analysis_{mapon}.png')
        plt.savefig(out_path, dpi=300)
        plt.close(fig)
        print(f"Figure sauvegardée pour mapon {mapon}: {out_path}")

def plot_nb_distrib_by_mapon(data, thresholds, out_dir):
    """
    Crée des graphiques pour analyser les distributions d'écarts et de tailles d'inversions 
    par mapon et méthode en utilisant des boxplots.
    
    Args:
        data: dictionnaire imbriqué {iso -> {method -> {chr -> {mapon -> list[(start, end)]}}}}
        thresholds: liste des seuils à tester
        out_dir: répertoire de sortie pour les figures
    """
    # Pour chaque mapon, on va créer une image
    all_mapons = set()
    all_methods = set()
    
    # Identifier tous les mapons et méthodes disponibles
    for iso in data:
        for method in data[iso]:
            all_methods.add(method)
            for chr_ in data[iso][method]:
                for mapon in data[iso][method][chr_]:
                    all_mapons.add(mapon)
    
    all_methods = sorted(list(all_methods))
    all_mapons = sorted(list(all_mapons))
    
    # Vérifier l'existence du répertoire de sortie
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    # Sélectionner les seuils spécifiques (0, max, et 3 intermédiaires multiples de 5)
    max_threshold = max(thresholds)
    selected_thresholds = [0]
    
    # Sélectionner 3 seuils intermédiaires équidistants et multiples de 5
    step = max_threshold / 4  # Pour avoir 3 points intermédiaires entre 0 et max
    for i in range(1, 4):
        # Trouver le multiple de 5 le plus proche du point équidistant
        intermediate = round(i * step / 5) * 5
        if intermediate > 0 and intermediate < max_threshold:
            selected_thresholds.append(intermediate)
    
    selected_thresholds.append(max_threshold)
    selected_thresholds = sorted(list(set(selected_thresholds)))  # Éliminer les doublons éventuels
        
    # Définir une palette de couleurs agréable
    palette = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    
    # Pour chaque mapon, créer une figure avec des subplots pour chaque méthode
    for mapon in all_mapons:
        # Créer une figure avec 2 rangées (distribution des écarts, distribution des tailles) 
        # et assez de colonnes pour chaque méthode
        fig, axes = plt.subplots(2, len(all_methods), figsize=(6*len(all_methods), 10))
        
        # Gérer le cas où il n'y a qu'une seule méthode
        if len(all_methods) == 1:
            axes = axes.reshape(2, 1)
        
        # Pour chaque méthode
        for m_idx, method in enumerate(all_methods):
            # Collecter les données pour chaque seuil sélectionné
            all_gaps_data = {}
            all_sizes_data = {}
            
            for threshold in selected_thresholds:
                # Récupérer les inversions et leurs données
                merged_inversions = merge_inversions(data, method, mapon, threshold)
                
                if merged_inversions:
                    # Collecter les écarts entre inversions consécutives
                    gaps = []
                    sorted_inversions = sorted(merged_inversions)
                    for i in range(1, len(sorted_inversions)):
                        prev_end = sorted_inversions[i-1][1]
                        curr_start = sorted_inversions[i][0]
                        gap = curr_start - prev_end - 1
                        if gap > 0:  # Ne considérer que les écarts positifs
                            gaps.append(gap)
                    
                    # Collecter les tailles des inversions
                    sizes = [end - start + 1 for start, end in merged_inversions]
                    
                    all_gaps_data[threshold] = gaps
                    all_sizes_data[threshold] = sizes
                else:
                    all_gaps_data[threshold] = []
                    all_sizes_data[threshold] = []
            
            # Vérifier si des données ont été trouvées pour cette méthode/mapon
            has_data = any(len(gaps) > 0 for gaps in all_gaps_data.values()) and \
                     any(len(sizes) > 0 for sizes in all_sizes_data.values())
            
            if has_data:
                # Graphique 1: Boxplot des écarts entre inversions
                boxplot_data_gaps = []
                boxplot_labels_gaps = []
                
                for threshold in selected_thresholds:
                    gaps = all_gaps_data[threshold]
                    if gaps and len(gaps) > 0:
                        boxplot_data_gaps.append(gaps)
                        boxplot_labels_gaps.append(f'Seuil={threshold}')
                
                if boxplot_data_gaps:
                    # Créer le boxplot horizontal
                    bplot1 = axes[0, m_idx].boxplot(boxplot_data_gaps, 
                                                  vert=False,  # Orientation horizontale
                                                  patch_artist=True,  # Remplir les boîtes
                                                  showfliers=False)  # Ne pas montrer les valeurs aberrantes
                    
                    # Utiliser des couleurs plus agréables
                    for i, patch in enumerate(bplot1['boxes']):
                        patch.set_facecolor(palette[i % len(palette)])
                    
                    # Configurer les étiquettes
                    axes[0, m_idx].set_yticklabels(boxplot_labels_gaps)
                    axes[0, m_idx].set_title(f'Distribution des écarts - {method}')
                    axes[0, m_idx].set_xlabel('Écart entre inversions (bases)')
                    
                    # Grille horizontale uniquement, avec style affiné
                    axes[0, m_idx].grid(True, axis='x', linestyle='--', alpha=0.7)
                    axes[0, m_idx].grid(False, axis='y')  # Désactiver la grille verticale
                else:
                    axes[0, m_idx].text(0.5, 0.5, 'Données insuffisantes', ha='center', va='center')
                
                # Graphique 2: Boxplot des tailles d'inversions
                boxplot_data_sizes = []
                boxplot_labels_sizes = []
                
                for threshold in selected_thresholds:
                    sizes = all_sizes_data[threshold]
                    if sizes and len(sizes) > 0:
                        boxplot_data_sizes.append(sizes)
                        boxplot_labels_sizes.append(f'Seuil={threshold}')
                
                if boxplot_data_sizes:
                    # Créer le boxplot horizontal
                    bplot2 = axes[1, m_idx].boxplot(boxplot_data_sizes, 
                                                  vert=False,  # Orientation horizontale
                                                  patch_artist=True,  # Remplir les boîtes
                                                  showfliers=False)  # Ne pas montrer les valeurs aberrantes
                    
                    # Utiliser des couleurs plus agréables
                    for i, patch in enumerate(bplot2['boxes']):
                        patch.set_facecolor(palette[i % len(palette)])
                    
                    # Configurer les étiquettes
                    axes[1, m_idx].set_yticklabels(boxplot_labels_sizes)
                    axes[1, m_idx].set_title(f'Distribution des tailles - {method}')
                    axes[1, m_idx].set_xlabel('Taille des inversions (bases)')
                    
                    # Grille horizontale uniquement, avec style affiné
                    axes[1, m_idx].grid(True, axis='x', linestyle='--', alpha=0.7)
                    axes[1, m_idx].grid(False, axis='y')  # Désactiver la grille verticale
                else:
                    axes[1, m_idx].text(0.5, 0.5, 'Données insuffisantes', ha='center', va='center')
                
            else:
                # Cas où aucune inversion n'est trouvée pour cette méthode/mapon
                axes[0, m_idx].text(0.5, 0.5, 'Aucune inversion', ha='center', va='center')
                axes[1, m_idx].text(0.5, 0.5, 'Aucune inversion', ha='center', va='center')
        
        # Titre global pour cette figure (mapon)
        fig.suptitle(f'Analyse des distributions pour le mapon: {mapon}', fontsize=16)
        plt.tight_layout(rect=[0, 0, 1, 0.96])  # Ajuster pour le titre
        
        # Sauvegarder la figure
        out_path = os.path.join(out_dir, f'distributions_boxplot_{mapon}.png')
        plt.savefig(out_path, dpi=300)
        plt.close(fig)
        print(f"Figure sauvegardée pour mapon {mapon}: {out_path}")

# ================== PLOT 5

def plot_venn_inv(data,genome_sizes,out_dir):
    
    mapons = set()
    methodes = set()
    for iso, methode_dict in data.items():
        methodes.update(methode_dict.keys())
        for methode, chr_dict in methode_dict.items():
            for chr_name, mapon_dict in chr_dict.items():
                mapons.update(mapon_dict.keys())
    mapons = sorted(mapons)
    methodes = sorted(methodes)

    for mapon in mapons:
        intersections = calculate_all_inversions_intersections(data, methodes,genome_sizes,mapon)


        isos = list(data.keys())
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
            df = intersections[iso]

            venn_tools.venn_4(methodes, df, title, ax=ax)
        # Supprime les axes en trop si n n'est pas un multiple de cols
        for j in range(idx + 1, len(axes)):
            fig.delaxes(axes[j])

        fig.legend(labels=methodes, loc='lower center', ncol=4, frameon=False)
        plt.tight_layout(rect=[0, 0.05, 1, 1])  # laisse de la place en bas pour la légende

        output_path = f"{out_dir}/{mapon}_venn_inv.png"
        plt.savefig(output_path, dpi=600, bbox_inches='tight')
        plt.close()