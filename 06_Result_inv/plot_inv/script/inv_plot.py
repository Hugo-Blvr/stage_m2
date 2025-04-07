import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import numpy as np
import sys
import yaml

def plot_inversions(fig, ax, df, chr_length, lineage_pathfile = None):

    if lineage_pathfile : 
        with open(lineage_pathfile, "r") as file:
            next(file)  # Skip header
            dico_corespondance = {parts[1]: parts[3] for parts in map(str.split, file)}
        df['t_lineage'] = df['target'].map(dico_corespondance)
        df['q_lineage'] = df['mapon'].map(dico_corespondance)
    else: 
        df['t_lineage'] = 'x'
        df['q_lineage'] = 'x'

    

    print(df)
    # Trier les données par lignée puis mapon
    df = df.sort_values(by='start',ascending = True)
    df = df.sort_values(by=['q_lineage', 'mapon'],ascending = False)
    ordered_mapons = df['mapon'].unique()

    # Palette de couleurs pour les lignées
    unique_lineages = df["q_lineage"].unique()
    palette = sns.color_palette("tab20", n_colors=len(unique_lineages))
    lineage_colors = dict(zip(unique_lineages[::-1], palette))

    # Couleurs pour les types d'inversions
    if 'type' in df.columns : colors = {'T': '#2CA02C', 'Q': '#D62728'}

    # Fonction pour minimiser les niveaux en regroupant les chevauchements
    def assign_levels(group_df):
        events = sorted([(row['start'], row['end'], idx) for idx, row in group_df.iterrows()])
        active_intervals = []  # [(end, level)]
        level_map = {}

        for start, end, idx in events:
            # Libérer les niveaux des inversions terminées
            active_intervals = [iv for iv in active_intervals if iv[0] > start]
            # Trouver le premier niveau disponible
            used_levels = {lvl for _, lvl in active_intervals}
            level = 0
            while level in used_levels: level += 1
            level_map[idx] = level
            active_intervals.append((end, level))

        return level_map

    # Appliquer l'optimisation des niveaux pour chaque mapon
    for mapon in ordered_mapons:
        group_df = df[df['mapon'] == mapon].copy()
        level_map = assign_levels(group_df)
        df.loc[group_df.index, 'optimized_level'] = df.loc[group_df.index].index.map(level_map)

    # Initialisation des tracés
    mapon_data = {}
    current_y_base = 0
    y_ticks, y_labels, y_colors = [], [], []

    # Organisation des données pour chaque mapon
    for mapon in ordered_mapons:
        group_df = df[df['mapon'] == mapon]
        lineage = group_df.iloc[0]['q_lineage']
        max_level = group_df['optimized_level'].max()

        mapon_data[mapon] = {
            'group_df': group_df,
            'lineage': lineage,
            'color': lineage_colors.get(lineage, 'black'),
            'max_level': max_level,
            'y_base': current_y_base
        }

        # Calcul du centre pour l'axe Y
        middle_y = current_y_base + max_level / 2
        y_ticks.append(middle_y)
        y_labels.append(str(mapon))
        y_colors.append(lineage_colors.get(lineage, 'black'))

        # Mise à jour du y_base avec un espacement minimal
        current_y_base += max_level + 1

    # Configuration de l'axe X

    ticks = np.linspace(0, chr_length, 5)
    ax.set_xticks(ticks)
    ax.set_xlim(0, chr_length)

    # Tracé des inversions
    for mapon in ordered_mapons:
        data = mapon_data[mapon]
        group_df = data['group_df']
        y_base = data['y_base']

        for _, row in group_df.iterrows():
            if 'type' in df.columns : 
                type_val = row['type'].strip("[]'").split("'")[0]
                color = colors.get(type_val, 'green')  
            else: color = 'green'
            y_level = y_base + row['optimized_level']
            ax.plot([row['start'], row['end']], [y_level, y_level], '-', 
                    linewidth=3, color=color, solid_capstyle='butt')

        # Ligne de séparation entre les mapons
        if mapon != ordered_mapons[-1]:
            ax.axhline(y=y_base + data['max_level'] + 0.5, color='gray', linestyle='--', alpha=0.7)

    # Création des légendes
    if 'type' in df.columns : custom_lines = [plt.Line2D([0], [0], color=colors['T'], lw=2), plt.Line2D([0], [0], color=colors['Q'], lw=2)]
    
    if 'type' in df.columns and len(df['type'].unique()) > 1 : 
        legend1 = ax.legend(custom_lines, ['Target', 'mapon'], loc='upper right')
        ax.add_artist(legend1)

    lineage_handles = [plt.Line2D([0], [0], marker='o', color='w', 
                                   markerfacecolor=lineage_colors[lineage], 
                                   markersize=8, label=lineage) 
                       for lineage in unique_lineages[::-1]]
    legend2 = ax.legend(handles=lineage_handles, title="Lineages",
                        loc='center right', bbox_to_anchor=(1.15, 0.5))

    

    # Configuration finale de l'axe des inversions
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_labels)
    ax.set_ylim(-1, current_y_base)
    ax.set_ylabel('Isolats')

    # Coloration des étiquettes selon les lignées
    for tick, color in zip(ax.get_yticklabels(), y_colors):
        tick.set_color(color)
        tick.set_fontsize(10)

    return ax, current_y_base

"""

def calculate_inversion_distribution(df, chr_length, window_size=25000):
    # Calculer le nombre de fenêtres
    n_windows = int(chr_length / window_size) + 1
    inversion_counts = np.zeros(n_windows)
    
    # Pour chaque inversion, compter le nombre de bases dans chaque fenêtre
    for _, row in df.iterrows():
        start = int(row['start'])
        end = int(row['end'])
        # Déterminer les fenêtres affectées
        start_window = start // window_size
        end_window = min(end // window_size, n_windows - 1)
        
        # Ajouter les contributions à chaque fenêtre
        for i in range(start_window, end_window + 1):
            window_start = i * window_size
            window_end = min((i + 1) * window_size, chr_length)
            # Calculer l'intersection entre l'inversion et la fenêtre
            overlap_start = max(start, window_start)
            overlap_end = min(end, window_end)
            if overlap_end > overlap_start:
                inversion_counts[i] += (overlap_end - overlap_start)
    
    # Normaliser par la longueur du chromosome
    normalized_counts = inversion_counts / window_size * 100
    # Positions centrales des fenêtres - CORRECTION: s'assurer que positions a la même taille que normalized_counts
    positions = np.arange(window_size/2, window_size/2 + n_windows*window_size, window_size)[:n_windows]
    
    # S'assurer que les deux tableaux ont la même taille
    if len(positions) > len(normalized_counts): positions = positions[:len(normalized_counts)]
    elif len(normalized_counts) > len(positions): normalized_counts = normalized_counts[:len(positions)]
    
    return positions, normalized_counts


def plot_inversion_distribution(fig, ax, df, chr_length, window_size=25000):
    # Calculer la distribution
    positions, normalized_counts = calculate_inversion_distribution(df, chr_length, window_size)
    
    # Tracer le barplot
    ax.bar(positions, normalized_counts, width=window_size*0.9, color='steelblue', alpha=0.7)
    ax.set_xlim(0, chr_length)
    ax.set_ylabel("Sc d'inv cumulée", fontsize=9)
    
    ax.set_xticks(np.linspace(0, chr_length, 5))
    ax.set_xticklabels([f"{int(t):,}" for t in ax.get_xticks()])
    
    return ax

"""

def plot_genome_features(df_inversions,chr_length,output_file):
    iso = df_inversions['iso'].unique()[0]
    chr = df_inversions['chr'].unique()[0]
    fig = plt.figure(figsize=(20, 5))  # Augmentation de la taille verticale pour ajouter le barplot
    
    # Créer une grille avec 2 sous-graphiques (inversions et barplot)
    gs = fig.add_gridspec(2, 1, height_ratios=[5, 1], hspace=0.2)
    ax_inv = fig.add_subplot(gs[0])
    ax_distrib = fig.add_subplot(gs[1], sharex=ax_inv)
    
    # Tracer les inversions
    plot_inversions(fig, ax_inv, df_inversions, chr_length)
    ax_inv.set_title(f'Inversions du {chr} de {iso}')
    ax_inv.set_xticklabels([])  # Ne pas afficher les étiquettes x pour l'axe des inversions
    
    # Tracer le barplot de distribution des inversions
    # -------- plot_inversion_distribution(fig, ax_distrib, df_inversions, chr_length)
    ax_distrib.set_xlabel(f"Position (bp)")
    ax_distrib.set_xticklabels([f"{int(t):,}" for t in ax_distrib.get_xticks()])
    
    plt.tight_layout(rect=[0, 0, 0.85, 1])
        
    
    # Enregistrer la figure
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Figure enregistrée dans {output_file}")



if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 inv_plot.py <yaml_pathname>")
        sys.exit(1)

    config_yaml = sys.argv[1]

    with open(config_yaml, 'r') as file: config = yaml.safe_load(file)
    # Parametre
    isolat = config['required']['isolat']
    chromosome = config['required']['chr']
    path_inv_file = config['required']['path_inv_file']
    # Parametre facultative 
    chr_len = config['optional']['chr_len']
    lineage_file = config['optional']['lineage_file']
    # Chemin de sortie
    path_out_plot = config['output']['path_out_plot'].format(isolat=isolat, chr=chr)

    df_inversions = pd.read_csv(path_inv_file, sep='\t')
    if chr_len: 
        try : 
            pd.read_csv(chr_len, sep='\t')
        except : chr_len = int(chr_len)
    else : 
        chr_len = int(df_inversions['end'].max())

    df_inversions = df_inversions[(df_inversions['iso'] == isolat) & (df_inversions['chr'] == chromosome)]
    plot_genome_features(df_inversions,chr_len,path_out_plot)
