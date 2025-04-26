import pandas as pd
import utils as uts
from collections import defaultdict


input_files = {"cactus_pairalgn_maf": "/home/hugob/stage_m2/06_Inversions_tools_comparaison/format_inv/output/inv_cactus_pairalgn_maf.tsv",  
               "cactus_PtoMalgn_maf": "/home/hugob/stage_m2/06_Inversions_tools_comparaison/format_inv/output/inv_cactus_PtoMalgn_maf.tsv", 
               "cactus_pairalgn_perbase": "/home/hugob/stage_m2/06_Inversions_tools_comparaison/format_inv/output/inv_cactus_pairalgn_perbase.tsv",
               "cactus_PtoMalgn_perbase": "/home/hugob/stage_m2/06_Inversions_tools_comparaison/format_inv/output/inv_cactus_PtoMalgn_perbase.tsv", 
               "syri_minimap_PtoMalgn":"/home/hugob/stage_m2/06_Inversions_tools_comparaison/format_inv/output/inv_syri_minimap_PtoMalgn.tsv",
               "syri_minimap_pairalgn":"/home/hugob/stage_m2/06_Inversions_tools_comparaison/format_inv/output/inv_syri_minimap_pairalgn.tsv",
               "syri_nucmer_pairalgn":"/home/hugob/stage_m2/06_Inversions_tools_comparaison/format_inv/output/inv_syri_nucmer_pairalgn.tsv",
               "syri_nucmer_PtoMalgn": "/home/hugob/stage_m2/06_Inversions_tools_comparaison/format_inv/output/inv_syri_nucmer_PtoMalgn.tsv"}
chr_len_file = "/home/hugob/stage_m2/01_Processing/output/chr_len.tsv"
lineage_file = "/home/hugob/stage_m2/00_Data/Isolates.Master.Info.txt"
output_path = "/home/hugob/stage_m2/06_Inversions_tools_comparaison/plot_inv/output"


all_data = uts.make_data(input_files)
#uts.print_data(all_data,2)
isos_in_all_data = all_data.keys()

chr_len_df = pd.read_csv(chr_len_file, sep='\t')
chr_len_df_filtered = chr_len_df[chr_len_df['iso'].isin(isos_in_all_data)]
chr_size_dict = defaultdict(lambda: defaultdict(dict))
for iso, chr, size in zip(chr_len_df_filtered['iso'], chr_len_df_filtered['chr'], chr_len_df_filtered['size']):
    chr_size_dict[iso][chr] = size
genome_size_dic = {iso: sum(chr_sizes.values()) for iso, chr_sizes in chr_size_dict.items()}

lineage_df = pd.read_csv(lineage_file, sep='\t', usecols=[1, 3], names=['iso', 'lineage'], header=0)
lineage_df_filtered = lineage_df[lineage_df['iso'].isin(isos_in_all_data)]
iso_lineage_dict = dict(zip(lineage_df_filtered['iso'], lineage_df_filtered['lineage']))


# ===================== PLOT 1 : DISTRIB ALL INV BY MAPON FOR ONE ISO ON ALL CHR
output_dir = f"{output_path}/distrib_inv_plot"
isolats = ["Gd_00293-aad","Gd_01882-ad",'Gd_00994-aaa']

for iso in isolats:
    data = all_data[iso]
    uts.plot_inv_by_mapon(data, iso, chr_size_dict,iso_lineage_dict,output_dir)


# ================== PLOT 2 : DISTRIB MULTIALGN INV BY METHODE FOR ALL ISO ON ALL CHR
output_dir = f"{output_path}/distrib_by_tools"
isolats = sorted(all_data)
methodes = ['cactus_PtoMalgn_maf','cactus_PtoMalgn_perbase','syri_minimap_PtoMalgn','syri_nucmer_PtoMalgn']
data = {iso: {methode: all_data[iso][methode] for methode in all_data[iso] if methode in methodes}
        for iso in isolats}

mapons = ['all','intra','inter']
uts.plot_inversion_by_mapon(data,chr_size_dict,mapons,output_dir)


# ================== PLOT 3 et 4 : 

output_dir = f"{output_path}/nb_size_ecart"
isolats = sorted(all_data)
methodes = ["cactus_PtoMalgn_maf", "cactus_PtoMalgn_perbase", "syri_minimap_PtoMalgn", "syri_nucmer_PtoMalgn"]
mapons = ['all','intra','inter']

data = {iso: {method: {chr_: {mapon: all_data[iso][method][chr_][mapon]
                for mapon in mapons if mapon in all_data[iso][method][chr_]}
            for chr_ in all_data[iso][method]}
        for method in  methodes if method in all_data[iso]}
    for iso in isolats if iso in all_data}

thresholds = list(range(0, 200, 1))
uts.plot_nb_size_by_mapon(data, thresholds, output_dir)

# ===

output_dir = f"{output_path}/nb_size_ecart"
isolats = sorted(all_data)
methodes = ["cactus_PtoMalgn_maf", "cactus_PtoMalgn_perbase", "syri_minimap_PtoMalgn", "syri_nucmer_PtoMalgn"]
mapons = ['all','intra','inter']

data = {iso: {method: {chr_: {mapon: all_data[iso][method][chr_][mapon]
                for mapon in mapons if mapon in all_data[iso][method][chr_]}
            for chr_ in all_data[iso][method]}
        for method in  methodes if method in all_data[iso]}
    for iso in isolats if iso in all_data}

thresholds = list(range(0, 2000, 50))
uts.plot_nb_distrib_by_mapon(data, thresholds, output_dir)


# ================== PLOT 5 : 
output_dir = f"{output_path}/venn_plot"
isolats = sorted(all_data)
methodes = ["cactus_PtoMalgn_maf", "cactus_PtoMalgn_perbase", "syri_minimap_PtoMalgn", "syri_nucmer_PtoMalgn"]
mapons = ['all','inter','intra']

data = {iso: {method: {chr_: {mapon: all_data[iso][method][chr_][mapon]
                for mapon in mapons if mapon in all_data[iso][method][chr_]}
            for chr_ in all_data[iso][method]}
        for method in  methodes if method in all_data[iso]}
    for iso in isolats if iso in all_data}

uts.plot_venn_inv(data,genome_size_dic,output_dir)

