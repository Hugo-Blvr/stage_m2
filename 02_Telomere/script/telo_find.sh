#!/bin/bash

path_script='.'

# Affichage de l'aide
usage() {
    echo "Usage: $0 -g <genome.fasta> -l <longreads.fastq> [-o <output_dir>]"
    echo "  -g    Path to the genome FASTA file (required)"
    echo "  -l    Path to the long reads FASTQ file (required)"
    echo "  -o    Output directory (default: current directory)"
    exit 1
}

# Initialisation des variables
outdir="."

# Lecture des arguments
while getopts ":g:l:o:h" opt; do
    case ${opt} in
        g ) genome="$OPTARG" ;;
        l ) longreads="$OPTARG" ;;
        o ) outdir="$OPTARG" ;;
        h ) usage ;;
        * ) echo "Invalid option: -$OPTARG" >&2; usage ;;
    esac
done

# Vérifier les paramètres obligatoires
[[ -z "$genome" || -z "$longreads" ]] && { echo "Error: Missing required arguments."; usage; }

# Création des répertoires de sortie en une seule commande
bam_output_dir="$outdir/bam"
soft_clip_output_dir="$outdir/soft_clip"
trf_output_dir="$outdir/trf"
mkdir -p "$bam_output_dir" "$soft_clip_output_dir" "$trf_output_dir"

# Définition des fichiers
genome_base="${genome%.fasta}"
mmi_file="$bam_output_dir/${genome_base}.mmi"
sorted_bam="$bam_output_dir/${genome_base}_mapping.sorted.bam"

# Étape 1 : Alignement des long reads
echo "Indexation du génome..."
minimap2 -d "$mmi_file" "$genome"
echo "Alignement des long reads..."
minimap2 -t 8 -ax map-ont "$mmi_file" "$longreads" | samtools view -bS | samtools sort -@ 8 -o "$sorted_bam"
echo "Alignement terminé : $sorted_bam"

# Étape 2 : Extraction des soft-clips d’intérêt
soft_clip_start_tsv="$soft_clip_output_dir/${genome_base}_softclip_start.tsv"
soft_clip_end_tsv="$soft_clip_output_dir/${genome_base}_softclip_end.tsv"
soft_clip_start_fasta="$soft_clip_output_dir/${genome_base}_softclip_start.fasta"
soft_clip_end_fasta="$soft_clip_output_dir/${genome_base}_softclip_end.fasta"

# Extraction des longueurs des chromosomes et écriture dans chr_lengths.tsv
echo "Extraction des longueurs des chromosomes..."
samtools idxstats "$sorted_bam" > "$outdir/chr_lengths.tsv"

# Lecture des longueurs des chromosomes
declare -A chr_lengths
while read -r chr len _; do
    chr_lengths["$chr"]=$len
done < "$outdir/chr_lengths.tsv"

echo "Extraction des soft-clips..."
samtools view "$sorted_bam" | awk -v OFS="\t" '
BEGIN {
    while ((getline < "chr_lengths.tsv") > 0) { chr_lengths[$1] = $2 }
    close("chr_lengths.tsv")
}
{
    if ($10 == "*") next  # Ignorer séquences vides
    chr = $3; start = $4
    if (!(chr in chr_lengths)) next
    ref_length = chr_lengths[chr]

    if (match($6, /^([0-9]+)S/, m)) { soft_clip = m[1] } else { soft_clip = 0 }

    cigar = $6; mapped_length = 0
    while (match(cigar, /([0-9]+)([MID])/)) {
        len = substr(cigar, RSTART, RLENGTH - 1)
        op = substr(cigar, RLENGTH, 1)
        if (op == "M" || op == "D") mapped_length += len
        cigar = substr(cigar, RSTART + RLENGTH)
    }

    end = start + mapped_length - 1
    min_start = ref_length - 8000
    min_end = ref_length - 1000

    if (start <= 1000 && end >= 8000 && soft_clip >= 50) {
        print $1, chr, start, end, $10 >> "'"$soft_clip_start_tsv"'"
    }

    if (start <= min_start && end >= min_end && soft_clip >= 50) {
        print $1, chr, start, end, $10 >> "'"$soft_clip_end_tsv"'"
    }
}'

echo "Formatages soft-clips..."
awk '{print ">"$2":"$1":"$3"-"$4 "\n" $5}' "$soft_clip_start_tsv" > "$soft_clip_start_fasta"
awk '{print ">"$2":"$1":"$3"-"$4 "\n" $5}' "$soft_clip_end_tsv" > "$soft_clip_end_fasta"

# Étape 3 : Détection des motifs avec TRF
echo "Recherche de motifs..."

cd $trf_output_dir
trf "../${soft_clip_start_fasta}" 2 7 7 80 4 10 500 -f -h -d
trf "../${soft_clip_end_fasta}" 2 7 7 80 4 10 500 -f -h -d
cd ..

dat_file_start=$(realpath "$trf_output_dir"/*start*.dat)
dat_file_end=$(realpath "$trf_output_dir"/*end*.dat)
filtered_repeats_start="$trf_output_dir/${genome_base}_filtered_repeats_start.fasta"
filtered_repeats_end="$trf_output_dir/${genome_base}_filtered_repeats_end.fasta"

# Filtrage des répétitions
filter_repeats() {
    awk 'BEGIN {seq_name=""}
        /^Sequence:/ {seq_name=$2}
        $3 >= 10 && $4 >= 6 {
            print ">" seq_name " repeat:" $1 "-" $2 " motif=" $14 " copies=" $4
            print $15
        }' "$1" > "$2"
}

filter_repeats "$dat_file_start" "$filtered_repeats_start"
filter_repeats "$dat_file_end" "$filtered_repeats_end"

# Étape 4 : Génération de la liste des télomères
python_script="$path_script/creat_list_motif.py"
out_name_start="$outdir/${genome_base}_list_telomere_start.tsv"
out_name_end="$outdir/${genome_base}_list_telomere_end.tsv"

python3 "$python_script" "$filtered_repeats_start" "$out_name_start"
python3 "$python_script" "$filtered_repeats_end" "$out_name_end"

echo "Analyse terminée. Résultats dans $outdir"
