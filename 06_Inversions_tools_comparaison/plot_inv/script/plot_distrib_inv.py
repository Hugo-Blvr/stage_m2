import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import math
import numpy as np
from collections import defaultdict



# Créer une visualisation par isolat
def plot_inv(df,chr_len_df,lineage_df, output_dir,out_finame):
    # Créer dict de taille des chromosomes pour accès rapide
    chr_size_dict = dict(zip(zip(chr_len_df['iso'], chr_len_df['chr']), chr_len_df['size']))
    # Obtenir palette de couleurs pour lignées
    unique_lineages = lineage_df['lineage'].unique()
    lineage_palette = dict(zip(unique_lineages, sns.color_palette("tab20", len(unique_lineages))))
    # Créer un mapping iso -> lineage
    iso_lineage = dict(zip(lineage_df['iso'], lineage_df['lineage']))
    # Organiser par iso, chr, mapon
    iso_chr_mapon_data = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for _, row in df.iterrows(): iso_chr_mapon_data[row['iso']][row['chr']][row['mapon']].append((row['start'], row['end']))
    # Couleur unique pour les inversions mappées
    inversion_color = 'green'


    for iso, chr_data in iso_chr_mapon_data.items():
        n_chrs = len(chr_data)
        n_rows = math.ceil(n_chrs / 3)  # 3 colonnes, donc diviser par 2 pour obtenir le nombre de lignes
        
        fig, axes = plt.subplots(n_rows, 3, figsize=(27, 4*n_rows))
        
        # Si un seul chromosome, ajuster les axes
        if n_chrs == 1:
            axes = np.array([[axes]])
        elif n_rows == 1:
            axes = axes.reshape(1, -1)
        
        # Parcourir les chromosomes
        for i, (chr_name, mapon_data) in enumerate(sorted(chr_data.items())):
            row_idx = i // 3
            col_idx = i % 3
            
            ax = axes[row_idx, col_idx]
            
            # Utiliser la taille du chromosome depuis le df si disponible, sinon calculer
            if (iso, chr_name) in chr_size_dict:
                max_pos = chr_size_dict[(iso, chr_name)]
            else:
                max_pos = max([end for mapon_invs in mapon_data.values() for _, end in mapon_invs])
            
            # Trier les mapons par lignée
            sorted_mapons = []
            mapon_lineages = {}
            
            # Récupérer la lignée pour chaque mapon
            for mapon in mapon_data.keys():
                lineage = iso_lineage.get(mapon, "Unknown")
                mapon_lineages[mapon] = lineage
                
            # Trier les mapons d'abord par lignée, puis par nom
            sorted_mapons = sorted(mapon_data.keys(), key=lambda m: (mapon_lineages.get(m, "Unknown"), m))
            sorted_mapons = sorted_mapons [::-1]
            # Une ligne par mapon
            for j, mapon in enumerate(sorted_mapons):
                y_pos = j
                inversions = mapon_data[mapon]
                
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
                lineage = mapon_lineages.get(mapon, "Unknown")
                color = lineage_palette.get(lineage, 'black')
                ax.get_yticklabels()[j].set_color(color)
            
            ax.grid(True, linestyle='--', alpha=0.7)
            ax.set_xlabel('Position (bp)')
        
        # Cacher les axes non utilisés dans la dernière ligne si nombre impair de chromosomes
        if n_chrs % 2 != 0:
            axes[n_rows-1, 1].axis('off')
        
        plt.tight_layout()
        plt.savefig(f"{output_dir}/{iso}_{out_finame}.png", dpi=600, bbox_inches='tight')
        plt.close(fig)


# ======================================== PLOTS
input_files = {"cactus_multialgn": "/home/hugob/stage_m2/06_bench_tools_inv_calling/format_inv/output/inv_cactus_multialgn.tsv",
               "cactus_parialgn": "/home/hugob/stage_m2/06_bench_tools_inv_calling/format_inv/output/inv_cactus_pairalgn.tsv",
               "cactus_PtoMalgn": "/home/hugob/stage_m2/06_bench_tools_inv_calling/format_inv/output/inv_cactus_PtoMalgn.tsv", 
               "syri_minimap_pairalgn":"/home/hugob/stage_m2/06_bench_tools_inv_calling/format_inv/output/inv_syri_minimap_pairalgn.tsv",
               "syri_minimap_PtoMalgn":"/home/hugob/stage_m2/06_bench_tools_inv_calling/format_inv/output/inv_syri_minimap_PtoMalgn.tsv",
               "syri_nucmer_pairalgn":"/home/hugob/stage_m2/06_bench_tools_inv_calling/format_inv/output/inv_syri_nucmer_pairalgn.tsv",
               "syri_nucmer_PtoMalgn": "/home/hugob/stage_m2/06_bench_tools_inv_calling/format_inv/output/inv_syri_nucmer_PtoMalgn.tsv"}

chr_len_file = "/home/hugob/stage_m2/01_Processing/output/chr_len.tsv"
lineage_file = "/home/hugob/stage_m2/00_Data/Isolates.Master.Info.txt"
output_dir = "/home/hugob/stage_m2/06_bench_tools_inv_calling/plot_inv/output/distrib_inv_plot"

# Charger les données
chr_len_df = pd.read_csv(chr_len_file, sep='\t')
lineage_df = pd.read_csv(lineage_file, sep='\t', usecols=[1, 3])
lineage_df.columns = ['iso', 'lineage']

for input_file in input_files.keys(): 
    df = pd.read_csv(input_files[input_file], sep='\t')
    if not 'mapon' in df.columns : df['mapon'] = ''
    df = df[df['iso'] == 'Gd_00293-aad']
    out_finame = input_file
    plot_inv(df,chr_len_df,lineage_df,output_dir,out_finame)
    print(f'Plot {input_file} terminée.')

