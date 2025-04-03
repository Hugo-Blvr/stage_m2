import os, re
from pathlib import Path
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import warnings
warnings.simplefilter(action='ignore', category=pd.errors.ParserWarning)
#pd.set_option('display.max_rows', None)



def get_chr_number(chr_name):
    num = re.findall(r'\d+', chr_name)
    return int(num[0]) if num else float('inf')  # Pour gérer les cas comme chrX, chrY

def plot_TE_size_distribution(df, lineage, out_dir='./te_size_boxplot', by_family=True):
    """
    Création d'histogrammes de distribution de taille des TE par chromosome
    
    Parameters:
    - df : DataFrame contenant les données de transposons
    - lineage : identifiant de lignée pour le nom du fichier de sortie
    - out_dir : répertoire de sortie pour les graphiques
    - by_family : si True, crée un graphique par famille de TE; sinon, crée un graphique global
    """
    os.makedirs(out_dir, exist_ok=True)
    
    # Récupérer les isolats uniques de cette lignée
    isolats = sorted(df['Isolat'].unique())
    
    # Vérifier s'il y a des données
    if len(isolats) == 0:
        print(f"Pas de données pour la lignée {lineage}")
        return
    
    # Calculer la longueur des TE
    df = df.copy()
    df['length'] = df['end'] - df['begin'] + 1
    
    # Familles de TE
    all_classes = sorted(df['class/family'].unique())
    colors = sns.color_palette("colorblind", len(all_classes))
    color_dict = dict(zip(all_classes, colors))
    
    # Pour chaque isolat
    for isolat in isolats:
        # Filtrer les TE pour cet isolat
        df_iso = df[df['Isolat'] == isolat]
        
        # Récupérer les chromosomes pour cet isolat
        chromosomes = sorted(df_iso['sequence'].unique(), key=lambda x: int(''.join(filter(str.isdigit, x))) if any(c.isdigit() for c in x) else x)
        
        if by_family:
            # Mode par famille: un graphique séparé pour chaque famille de TE
            for te_class in all_classes:
                # Filtrer pour cette famille de TE
                df_class = df_iso[df_iso['class/family'] == te_class]
                
                # S'il n'y a pas de données pour cette famille, passer à la suivante
                if len(df_class) == 0:
                    continue
                
                # Déterminer la disposition des sous-graphiques (un par chromosome)
                n_chr = len(chromosomes)
                n_cols = min(3, n_chr)  # Maximum 3 colonnes
                n_rows = (n_chr + n_cols - 1) // n_cols
                
                # Créer la figure pour cette famille de TE
                fig = plt.figure(figsize=(15, 6 * n_rows))
                fig.suptitle(f"Distribution de taille des TE - Famille {te_class} - Isolat {isolat} (Lignée {lineage})", fontsize=18)
                
                # Traiter chaque chromosome
                for chr_idx, chromosome in enumerate(chromosomes, 1):
                    ax = plt.subplot(n_rows, n_cols, chr_idx)
                    
                    # Données pour ce chromosome et cette famille
                    df_chr_class = df_class[df_class['sequence'] == chromosome]
                    
                    if len(df_chr_class) > 0:
                        # Créer l'histogramme
                        ax.hist(
                            df_chr_class['length'], 
                            bins=30,  # Nombre de bins pour l'histogramme
                            color=color_dict[te_class],
                            edgecolor='black',
                            linewidth=0.5
                        )
                        
                        # Ajouter une ligne pour la moyenne
                        mean_length = df_chr_class['length'].mean()
                        ax.axvline(x=mean_length, color='red', linestyle='--', 
                                   label=f'Moyenne: {mean_length:.0f} pb')
                        
                        # Ajouter la médiane
                        median_length = df_chr_class['length'].median()
                        ax.axvline(x=median_length, color='green', linestyle=':', 
                                   label=f'Médiane: {median_length:.0f} pb')
                        
                        # Ajouter d'abord l'annotation avec n = X à droite
                        n_ann = ax.annotate(f'n = {len(df_chr_class)}', 
                                           xy=(0.55, 0.95),  # Position à droite en haut
                                           xycoords='axes fraction',
                                           ha='right', 
                                           va='top', 
                                           bbox=dict(boxstyle='round', fc='white', alpha=0.8))
                        
                        # Positionner la légende à droite de l'annotation
                        ax.legend(loc='upper right', fontsize=8, bbox_to_anchor=(0.98, 0.98))

                    ax.set_title(f"Chromosome {chromosome}", fontsize=14)
                    ax.set_xlabel("Longueur (pb)", fontsize=12)
                    ax.set_ylabel("Nombre d'éléments", fontsize=12)
                
                plt.tight_layout(rect=[0, 0.03, 1, 0.95])
                plt.savefig(f'{out_dir}/TE_size_distribution_{isolat}_{lineage}_{te_class.replace("/", "_")}.svg', 
                           dpi=300, bbox_inches='tight', format='svg')
                plt.close()
                
        else:
            # Mode global: un seul graphique pour toutes les familles avec barres cumulatives
            # Déterminer la disposition des sous-graphiques (un par chromosome)
            n_chr = len(chromosomes)
            n_cols = min(3, n_chr)  # Maximum 3 colonnes
            n_rows = (n_chr + n_cols - 1) // n_cols
            
            # Créer la figure pour l'ensemble des TE
            fig = plt.figure(figsize=(18, 6 * n_rows))
            
            fig.suptitle(f"Distribution de taille des TE - Toutes familles - Isolat {isolat} (Lignée {lineage})", fontsize=18)
            
            # Traiter chaque chromosome
            for chr_idx, chromosome in enumerate(chromosomes, 1):
                ax = plt.subplot(n_rows, n_cols, chr_idx)
                
                # Données pour ce chromosome
                df_chr = df_iso[df_iso['sequence'] == chromosome]
                
                if len(df_chr) > 0:
                    # Déterminer les bins pour garantir la cohérence
                    bins = np.linspace(df_chr['length'].min(), df_chr['length'].max(), 31)
                    
                    # Histogramme empilé avec différentes couleurs par famille
                    hist_data = []
                    labels = []
                    
                    for te_class in all_classes:
                        df_chr_class = df_chr[df_chr['class/family'] == te_class]
                        if len(df_chr_class) > 0:
                            hist_data.append(df_chr_class['length'].values)
                            labels.append(te_class)
                    
                    if hist_data:  # S'assurer qu'il y a des données
                        # Créer l'histogramme cumulatif
                        n, bins, patches = ax.hist(
                            hist_data,
                            bins=bins,
                            stacked=True,
                            label=labels,
                            edgecolor='black',
                            linewidth=0.5,
                            color=[color_dict[label] for label in labels]
                        )
                        
                        # Ajouter une ligne pour la moyenne globale
                        mean_length = df_chr['length'].mean()
                        ax.axvline(x=mean_length, color='red', linestyle='--', 
                                  label=f'Moyenne: {mean_length:.0f} pb')
                        
                        # Ajouter la médiane globale
                        median_length = df_chr['length'].median()
                        ax.axvline(x=median_length, color='green', linestyle=':', 
                                  label=f'Médiane: {median_length:.0f} pb')
                        
                        # Ajouter d'abord l'annotation avec n = X à droite
                        n_ann = ax.annotate(f'n = {len(df_chr)}', 
                                           xy=(0.65, 0.95),  # Position à droite en haut
                                           xycoords='axes fraction',
                                           ha='right', 
                                           va='top', 
                                           bbox=dict(boxstyle='round', fc='white', alpha=0.8))
                        
                        # Positionner la légende à droite de l'annotation
                        ax.legend(loc='upper right', fontsize=8, bbox_to_anchor=(0.98, 0.98))
                    
                    ax.set_title(f"Chromosome {chromosome}", fontsize=14)
                    ax.set_xlabel("Longueur (pb)", fontsize=12)
                    ax.set_ylabel("Nombre d'éléments", fontsize=12)
            
            # Ajuster la mise en page
            plt.tight_layout(rect=[0, 0.03, 1, 0.95])
            
            # Sauvegarder la figure
            plt.savefig(f'{out_dir}/TE_size_distribution_{isolat}_{lineage}_all.svg', 
                       dpi=300, bbox_inches='tight', format='svg')
            plt.close()
                
        print(f"Graphiques créés pour l'isolat {isolat}")

def distrib_TE_byIsolat(df, isolat, out_dir='.', window_size=10000):
    # Create output directory
    os.makedirs(out_dir, exist_ok=True)
    
    # Sort chromosomes 
    chromosomes = sorted(df['sequence'].unique(), key=get_chr_number)
    # Colorblind-friendly color palette
    all_classes = sorted(df['class/family'].unique())
    colors = sns.color_palette("colorblind", len(all_classes))
    color_dict = dict(zip(all_classes, colors))

    # Calculate subplot grid
    n_chr, n_cols = len(chromosomes), 5
    n_rows = (n_chr + n_cols - 1) // n_cols

    # Create figure with subplots
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(25, 5*n_rows), squeeze=False)
    axes = axes.flatten()
        
    # Process each chromosome
    for idx, chromosome in enumerate(chromosomes):
        df_chr = df[df['sequence'] == chromosome].copy()
        
        # Get max chromosome position
        max_position = df_chr['end'].max()
        
        # Create the plot for this chromosome
        ax = axes[idx]
        # Plot KDE for each TE class/family
        for te_class in all_classes:
            # Filter for specific TE class
            df_class = df_chr[df_chr['class/family'] == te_class]
            
            # If no elements of this class, skip
            if len(df_class) == 0:
                continue
            
            # Use kernel density estimation
            kde_positions = df_class['begin']
            kde = sns.kdeplot(
                x=kde_positions, 
                ax=ax, 
                color=color_dict[te_class], 
                alpha=0.8,  # transparence 
                linewidth=3,  # épaisseur de la ligne
                label=te_class,
                warn_singular=False
            )
        
        # Customize the plot
        ax.set_title(f"TE Distribution on {chromosome}")
        ax.set_xlabel("Genomic Position (Mb)")
        ax.set_ylabel("Density")
        # Set x-ticks to Mb
        ax.set_xlim(0, max_position)
        xticks = np.linspace(0, max_position, 6)  # 6 positions pour plus de lisibilité
        ax.set_xticks(xticks)
        ax.set_xticklabels([f"{x/1e6:.1f}" for x in xticks], rotation=45)
        # Remove individual legends
        if ax.get_legend(): ax.get_legend().remove()

    # Remove any unused subplots
    for idx in range(len(chromosomes), len(axes)):
        fig.delaxes(axes[idx])

    # Create a common legend
    handles = [plt.Rectangle((0,0),1,1, color=color_dict[label], alpha=0.5) for label in all_classes]
    fig.legend(handles, all_classes, 
               bbox_to_anchor=(1.05, 0.5), 
               loc='center left',
               title="TE Families", 
               title_fontsize=28, 
               fontsize=25)

    # Adjust layout and save
    plt.tight_layout()
    plt.savefig(f'{out_dir}/TE_distribution_{isolat}.svg', dpi=300, bbox_inches='tight')
    plt.close()

def pie_plot_TEfam_by_isolat(df, lineage, genome_size_path, out_dir='.'):
    """
    Création de pie charts pour les familles de transposons par lignée
    
    Parameters:
    - df : DataFrame contenant les données de transposons
    - lineage : identifiant de lignée pour le nom du fichier de sortie
    - genome_size_path : chemin vers le fichier de tailles de génome
    - out_dir : répertoire de sortie pour les graphiques
    """
    os.makedirs(out_dir, exist_ok=True)
    
    # Charger les tailles de génome
    genome_sizes_df = pd.read_csv(genome_size_path, sep='\s+')
    
    def autopct_format(pct):
        """Formater le pourcentage pour l'affichage"""
        return f'{pct:.1f}%' if pct > 2 else ''
    
    # Calculer la longueur des TE
    df = df.copy()
    df['length'] = df['end'] - df['begin'] + 1
    
    # Récupérer les isolats uniques de cette lignée
    isolats = sorted(df['Isolat'].unique())
    
    # Vérifier s'il y a des données
    if len(isolats) == 0:
        print(f"Pas de données pour la lignée {lineage}")
        return
    
    # Déterminer la disposition des sous-graphiques
    n_isolats = len(isolats)
    n_cols = min(4, n_isolats)  # Maximum 5 colonnes
    n_rows = (n_isolats + n_cols - 1) // n_cols 
    
    # Créer la figure
    fig = plt.figure(figsize=(8 * n_cols, 8 * n_rows))
    fig.suptitle(f"Proportion of TE families per genome - Lineage {lineage}", fontsize=30)
    
    # Couleurs pour toutes les familles
    all_classes = sorted(df['class/family'].unique())
    all_classes.append('Norepeat')  # Ajouter la catégorie non-répétitive
    colors = sns.color_palette("colorblind", len(all_classes)).as_hex()
    color_dict = dict(zip(all_classes, colors))
    
    # Traiter chaque isolat
    for idx, isolat in enumerate(isolats, 1):
        plt.subplot(n_rows, n_cols, idx)
        
        # Filtrer le génome pour cet isolat
        try:
            genome_row = genome_sizes_df[genome_sizes_df['Isolat'] == isolat]
            if len(genome_row) == 0:
                print(f"Pas d'information de taille de génome pour {isolat}")
                continue
            
            genome_size = genome_row.iloc[0, 1:].sum()
        except Exception as e:
            print(f"Erreur lors du traitement de {isolat}: {e}")
            continue
        
        # Filtrer les TE pour cet isolat
        df_iso = df[df['Isolat'] == isolat]
        
        # Calculer la taille totale des TE par famille
        te_sizes_by_family = df_iso.groupby('class/family')['length'].sum()
        
        # Calculer la taille totale des TE
        total_te_size = te_sizes_by_family.sum()
        
        # Créer un dictionnaire avec les tailles des TE et la région non-répétitive
        sizes_dict = te_sizes_by_family.to_dict()
        sizes_dict['Norepeat'] = max(0, genome_size - total_te_size)
        
        # Créer le pie chart
        wedges, texts, autotexts = plt.pie(
            list(sizes_dict.values()), 
            labels=None, 
            autopct=autopct_format,
            startangle=90, 
            colors=[color_dict[cat] for cat in sizes_dict.keys()]
        )
        
        plt.setp(autotexts, size=15, weight="bold")
        plt.title(f"{isolat}\n({genome_size:,} bp)", fontsize=16)
    
    # Légende commune
    fig.legend(
        [f"{cat}" for cat in all_classes], 
        title="TE families", 
        loc="center left", 
        bbox_to_anchor=(1, 0.5),
        title_fontsize=20, 
        fontsize=15
    )
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(f'{out_dir}/TE_genome_proportion_{lineage}.svg', 
                dpi=300, bbox_inches='tight', format='svg')
    plt.close()

def plot_divergence_distribution_by_isolat(df, isolat, out_dir):
    """
    Génère un histogramme + KDE pour la distribution de 'div.' par isolat.
    
    :param df: DataFrame contenant les données des TE pour un isolat donné.
    :param isolat: Nom de l'isolat en cours d'analyse.
    :param out_dir: Dossier de sortie pour enregistrer les graphiques.
    """
    if 'div.' not in df.columns:
        print(f"⚠️  La colonne 'div.' est absente pour {isolat}.")
        return

    plt.figure(figsize=(8, 6))
    sns.histplot(df['div.'], kde=True, bins=30, color='royalblue', edgecolor="black", alpha=0.7)
    
    plt.xlabel("Divergence (%)")
    plt.ylabel("Fréquence")
    plt.title(f"Distribution de la divergence des TE ({isolat})")
    
    os.makedirs(out_dir, exist_ok=True)
    plt.savefig(os.path.join(out_dir, f"{isolat}_divergence_distribution.svg"), format='svg', dpi=300)
    plt.close()


if __name__ == "__main__":
    chrLen_path = "/home/hugob/work/len_chr.tsv"

    input_file = '/home/hugob/work/05_TE/all_TE.csv'
    df = pd.read_csv(input_file, sep=";",low_memory=False)
    df['Lineage'] = df['Lineage'].replace({'Psp1': 'Psp'})
    df_filtre = df

    isolats=[]
    #isolats = sorted(df_filtre['Isolat'].unique())
    for isolat in isolats:
        print(isolat)
        df_iso = df_filtre[df_filtre['Isolat'] == isolat]
        # ----- PLOT 1 --> TE le long des chr 
        distrib_TE_byIsolat(df_iso, isolat, './distribTe')

    lineages = []
    #lineages = sorted(df_filtre['Lineage'].unique())
    for lineage in lineages:
        print(lineage)
        df_line = df_filtre[df_filtre['Lineage'] == lineage]
        # ------PLOT 2 --> camembert répartition class TE par génome
        pie_plot_TEfam_by_isolat(df_line, lineage,chrLen_path, './piechart')

    lineages = []
    #lineages = sorted(df_filtre['Lineage'].unique())
    for lineage in lineages:
        print(lineage)

        df_line = df_filtre[df_filtre['Lineage'] == lineage]        
        # Plot de distribution de taille avec des boîtes à moustaches
        plot_TE_size_distribution(df_line, lineage, by_family=True)

        
    lineages = sorted(df_filtre['Lineage'].unique())
    for lineage in lineages:
        print(lineage)
        df_line = df_filtre[df_filtre['Lineage'] == lineage]        
        isolats = sorted(df_line['Isolat'].unique())
        for isolat in isolats:
            df_iso = df_line[df_line['Isolat'] == isolat]
            plot_divergence_distribution_by_isolat(df_iso, isolat, out_dir='./distribDIV')