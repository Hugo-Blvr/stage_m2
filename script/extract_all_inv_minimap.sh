#!/bin/bash

# Fonction help
show_help() {
    echo "Usage: $0 [options]"
    echo
    echo "Analyse des inversions génomiques entre les fichiers FASTA dans un dossier."
    echo
    echo "Options:"
    echo "  -h, --help               Affiche ce message d'aide"
    echo "  -d, --directory DIR      Dossier contenant les fichiers FASTA à analyser"
    echo "  -o, --output FILE        Spécifie le fichier de sortie (défaut: ./all_inv_extract_minimap_05.tsv)"
    echo "  -i, --identity VALUE     Spécifie l'identité de séquence minimale de séquence (0-1, défaut: 0.99)"
    echo "  -t, --threads NUM        Nombre de threads pour minimap2 (défaut: 8)"
    echo
    echo "Exemple:"
    echo "  $0 -d data/fasta_folder -o results.tsv -i 0.95"
    echo
}

# Valeurs par défaut
OUTPUT_FILE="./inv_calling.tsv"
MIN_IDENTITY=0
THREADS=8
INPUT_DIR=""

# Traitement des arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_help
            exit 0
            ;;
        -d|--directory)
            INPUT_DIR="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_FILE="$2"
            shift 2
            ;;
        -i|--identity)
            MIN_IDENTITY="$2"
            if ! [[ "$MIN_IDENTITY" =~ ^[0-9]*\.?[0-9]+$ ]] || (( $(echo "$MIN_IDENTITY < 0 || $MIN_IDENTITY > 1" | bc -l) )); then
                echo "Erreur: L'identité de séquence doit être un nombre entre 0 et 1."
                exit 1
            fi
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            if ! [[ "$THREADS" =~ ^[0-9]+$ ]]; then
                echo "Erreur: Le nombre de threads doit être un entier positif."
                exit 1
            fi
            shift 2
            ;;
        -*)
            echo "Option inconnue: $1"
            show_help
            exit 1
            ;;
        *)
            echo "Erreur: Arguments non reconnus."
            show_help
            exit 1
            ;;
    esac
done

# Vérification si le dossier d'entrée est spécifié
if [ -z "$INPUT_DIR" ]; then
    echo "Erreur: Vous devez spécifier un dossier d'entrée avec l'option -d ou --directory."
    show_help
    exit 1
fi

# Vérification si l'argument est un dossier
if [ ! -d "$INPUT_DIR" ]; then
    echo "Erreur: $INPUT_DIR n'est pas un dossier valide."
    exit 1
fi

# Affichage des paramètres
echo
echo "=== Paramètres de l'analyse ==="
echo "Dossier d'entrée: $INPUT_DIR"
echo "Fichier de sortie: $OUTPUT_FILE"
echo "Identité de séquence minimale: $MIN_IDENTITY"
echo "Nombre de threads: $THREADS"
echo "============================="
echo

# Récupération de la liste des fichiers FASTA dans le dossier
fasta_files=("$INPUT_DIR"/*.fasta)

# Vérification de la présence de fichiers FASTA
if [ ${#fasta_files[@]} -eq 0 ]; then
    echo "Erreur: Aucun fichier FASTA trouvé dans le dossier $INPUT_DIR."
    exit 1
fi

# Extraction du répertoire de sortie à partir du chemin du fichier de sortie
OUTPUT_DIR=$(dirname "$OUTPUT_FILE")
if [ "$OUTPUT_DIR" = "." ]; then
    OUTPUT_DIR=$(pwd)
fi

# Vérification si le répertoire de sortie existe
if [ ! -d "$OUTPUT_DIR" ]; then
    echo "Erreur: Le répertoire de sortie $OUTPUT_DIR n'existe pas."
    exit 1
fi

# Définition des fichiers/dossiers de sortie dans le même répertoire
error_log="$OUTPUT_DIR/extraction_errors.log"
output_dir="$OUTPUT_DIR/tmp_dir_output"

# Nettoyage et préparation des fichiers/dossiers de sortie
[ -f "$OUTPUT_FILE" ] && rm "$OUTPUT_FILE"
[ -f "$error_log" ] && rm "$error_log"
mkdir -p "$output_dir"
[ "$(ls -A "$output_dir" 2>/dev/null)" ] && rm "$output_dir"/*

# Affichage du nombre de fichiers trouvés
echo "Nombre de fichiers FASTA trouvés: ${#fasta_files[@]}"

# Calcul du nombre total de comparaisons
total_comparisons=$(( ${#fasta_files[@]} * (${#fasta_files[@]} - 1) ))
current_comparison=0

for ((i=0; i<${#fasta_files[@]}; i++)); do
    for ((j=0; j<${#fasta_files[@]}; j++)); do
        if [[ $i -ne $j ]]; then 
            current_comparison=$((current_comparison + 1))
            target_file="${fasta_files[i]}"
            query_file="${fasta_files[j]}"
            target_name=$(basename "${target_file}" | cut -d. -f1)
            query_name=$(basename "${query_file}" | cut -d. -f1)
            
            # Calcul et affichage du pourcentage de progression
            percent=$((current_comparison * 100 / total_comparisons))
            echo "[$percent%] Extraction des inversions: $target_name - $query_name"
            
            # Extraire les chromosomes pour chaque fichier et garder les communs
            chrs_query=$(grep -o '^>[^ ]*' "$query_file" | sed 's/^>//' | sort | uniq)
            chrs_target=$(grep -o '^>[^ ]*' "$target_file" | sed 's/^>//' | sort | uniq)
            common_chrs=$(comm -12 <(echo "$chrs_target") <(echo "$chrs_query"))

            # Boucler sur les chromosomes communs
            for chr in $common_chrs; do
                echo -ne "\rTraitement du chromosome: $chr"
                
                # Noms de sortie
                chr_target_file="$output_dir/target.fasta"
                chr_query_file="$output_dir/query.fasta"
                output_bam="$output_dir/algnt.bam"
                syri_out="$output_dir/syri.out"

                # Extraire les séquences du chromosome commun
                awk -v chr=">$chr" '
                    $0 ~ chr"([ \t]|$)" {print_flag=1} 
                    $0 ~ "^>" && $0 !~ chr"([ \t]|$)" {print_flag=0} 
                    print_flag
                ' "$target_file" > "$chr_target_file"

                awk -v chr=">$chr" '
                    $0 ~ chr"([ \t]|$)" {print_flag=1} 
                    $0 ~ "^>" && $0 !~ chr"([ \t]|$)" {print_flag=0} 
                    print_flag
                ' "$query_file" > "$chr_query_file"

                # Exécuter minimap2 et filtrer les résultats
                minimap2 -ax asm5 -t $THREADS --cs --eqx "$chr_target_file" "$chr_query_file" 2>/dev/null | \
                samtools view -h | \
                awk -v min_identity="$MIN_IDENTITY" 'BEGIN {OFS="\t"} 
                    /^@/ {print; next} 
                    {                        
                        matches = mismatches = indels = 0;
                        cigar = $6;
                        while (match(cigar, /([0-9]+)([=XID])/)) {
                            val = substr(cigar, RSTART, RLENGTH-1);
                            type = substr(cigar, RSTART+RLENGTH-1, 1);
                            if (type == "=") matches += val;
                            else if (type == "X") mismatches += val;
                            else if (type == "I" || type == "D") indels += val;
                            cigar = substr(cigar, RSTART+RLENGTH);
                        }
                        
                        aligned_bases = matches + mismatches + indels;
                        if (aligned_bases > 0) {
                            identity = matches / aligned_bases;
                            if (identity >= min_identity) print;
                        }
                    }' | \
                samtools sort -O BAM -o "$output_bam" 2>>"$error_log"
                samtools index "$output_bam" 2>>"$error_log"

                # Exécuter syri et vérifier si syri.out existe
                syri -c "$output_bam" -r "$chr_target_file" -q "$chr_query_file" -F B \
                    --dir "$output_dir" \
                    --nc $THREADS > /dev/null 2>&1

                if [ ! -f "$syri_out" ]; then
                    echo -e "\nErreur pour $target_name - $query_name : $chr" | tee -a "$error_log"
                    continue
                fi

                # Extraction et formatage des inversions génomiques: filtre les lignes INV, calcule les coordonnées et tailles, 
                # déduplique les résultats, puis écrit au format TSV dans le fichier de sortie
                awk -F'\t' -v OFS='\t' -v tgt_name="$target_name" -v qry_name="$query_name" '
                ($9 ~ /INV/) {
                    tstart = ($2 < $3) ? $2 : $3;
                    tend = ($2 > $3) ? $2 : $3;
                    qstart = ($7 < $8) ? $7 : $8;
                    qend = ($7 > $8) ? $7 : $8;
                    pair = tstart "-" tend "-" qstart "-" qend;
                    sep = "-"
                    if (!(pair in seen)) {
                        seen[pair] = 1;
                        tgt_size = tend - tstart;
                        qry_size = qend - qstart;
                        print tgt_name, $1, tstart, tend, tgt_size, sep, qry_name, $6, qstart, qend, qry_size, $9, $10
                    }
                }' "$syri_out" >> "$OUTPUT_FILE"
            done
            echo
        fi
    done
done

# Nettoyage et rapport final
rm -r "$output_dir"

echo
if [ -f "$OUTPUT_FILE" ]; then
    inversions_count=$(wc -l < "$OUTPUT_FILE")
    echo "Traitement terminé avec succès !"
    echo "Nombre d'inversions détectées: $inversions_count"
    echo "Résultats enregistrés dans: $OUTPUT_FILE"
    
    if [ -f "$error_log" ]; then
        if [ -s "$error_log" ]; then
            error_count=$(wc -l < "$error_log")
            echo "Attention: $error_count erreurs ont été enregistrées dans $error_log"
        else
            # Supprime le fichier log s'il existe mais est vide
            rm "$error_log"
        fi
    fi
else
    echo "Aucune inversion n'a été détectée."
fi