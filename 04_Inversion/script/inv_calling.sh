#!/bin/bash
#
# ============================== DESC =====================================
#
#
#
# ========================================================================

RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[0;33m';
BLUE='\033[0;34m'; NC='\033[0m'

help() {
    echo -e "\n${BLUE}Usage:${NC} $0 -d <in_directory> -o <out_directory/file.tsv> [options]\n"
    echo -e "${BLUE}Options obligatoires:${NC}"
    echo "  -d, --directory <dossier>   Dossier contenant les fichiers FASTA à analyser"
    echo -e "${BLUE}Options facultatives:${NC}"
    echo "  -h, --help                  Affiche ce message d'aide"
    echo "  -o, --output <paht_file>    Chemin du fichier de sortie des résultats (défaut: ./inv_calling.tsv)"
    echo "  -i, --identity <float>     Identité de séquence minimale (0-1, défaut: 0)"
    echo "  -t, --threads <int>      Nombre de threads à utiliser (défaut: 8)\n"    
    echo -e "\n${YELLOW}IMPORTANT: ${NC}Les chromosomes homologues doivent avoir le même identifiant et être sur le même brin.\n"
    exit 0
}

verifier_prerequis() {
    local missing_tools=()
    for tool in minimap2 samtools awk syri; do
        command -v "$tool" &>/dev/null || missing_tools+=("$tool")
    done

    if ((${#missing_tools[@]})); then
        echo -e "\n${RED}Erreur: Les outils suivants ne sont pas installés ou ne sont pas dans le PATH:${NC}"
        printf "${RED}  - %s\n${NC}" "${missing_tools[@]}"
        echo -e "Veuillez installer ces outils avant d'exécuter ce script.\n"
        exit 1
    fi
}

verifier_dossier() {
    local dossier="$1"
    local description="$2"
    # Vérification de l'existence et de la lisibilité du dossier
    if [ ! -d "$dossier" ] || [ ! -r "$dossier" ]; then
        echo -e "${RED}Erreur: Le dossier $description '$dossier' n'existe pas ou n'est pas lisible.${NC}"
        journaliser "ERREUR: Dossier '$dossier' introuvable ou non lisible"
        exit 1
    fi
}

creer_dossier() {
    local dossier="$1"
    # Si le dossier n'existe pas, tente de le créer
    if [ ! -d "$dossier" ]; then
        journaliser "Création du dossier '$dossier'"
        mkdir -p "$dossier" || { 
            echo -e "${RED}Erreur: Impossible de créer le dossier '$dossier'.${NC}"
            journaliser "ERREUR: Création du dossier '$dossier' échouée"
            exit 1
        }
    fi
    # Vérifie que le dossier est accessible en écriture
    if [ ! -w "$dossier" ]; then
        echo -e "${RED}Erreur: Le dossier '$dossier' n'est pas accessible en écriture.${NC}"
        journaliser "ERREUR: Dossier '$dossier' non accessible en écriture"
        exit 1
    fi
}

journaliser() { echo -e "[$(date '+%Y-%m-%d %H:%M:%S')] $1" >> "$LOG_FILE"; }

executer_commande() {
    local description="$1"           
    local commande="$2"              
    local ignore_error="${3:-false}" # Un flag qui indique si les erreurs doivent être ignorées (par défaut "false")
    
    journaliser "EXÉCUTION: $description" 
    
    OUTPUT=$(eval "$commande" 2>&1) # Exécution de la commande dans un sous-shell et récupération de la sortie
    local STATUS=$?  # Récupération du code de retour de la commande
    
    if [ $STATUS -eq 0 ] || [ "$ignore_error" == "true" ]; then
        return 0  # si la commande s'est exécutée avec succès ou si on ignore les erreurs
    else
        # Si une erreur s'est produite, afficher un message d'erreur et journaliser
        echo -e "  ${RED}✗${NC} Erreur lors de $description"
        journaliser "ERREUR: $description: $OUTPUT"
        return 1
    fi
}

# Initialisation des variables par défaut
OUTFILE="./inv_calling.tsv"
MIN_IDENTITY=0
THREADS=8
INDIR=""

# Traitement des options
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -h|--help)
            help
            ;;
        -d|--directory)
            INDIR="$2"
            shift 2
            ;;
        -o|--output)
            OUTFILE="$2"
            shift 2
            ;;
        -i|--identity)
            MIN_IDENTITY="$2"
            if ! [[ "$MIN_IDENTITY" =~ ^[0-9]*\.?[0-9]+$ ]] || (( $(echo "$MIN_IDENTITY < 0 || $MIN_IDENTITY > 1" | bc -l) )); then
                echo -e "${RED}Erreur: L'identité de séquence doit être un nombre entre 0 et 1.${NC}"
                exit 1
            fi
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            if ! [[ "$THREADS" =~ ^[0-9]+$ ]]; then
                echo -e "${RED}Erreur: Le nombre de threads doit être un entier positif.${NC}"
                exit 1
            fi
            shift 2
            ;;
        *)
            echo -e "${RED}Erreur: Option inconnue '$1'${NC}"
            echo "Utilisez '$0 --help' pour voir les options disponibles."
            exit 1
            ;;
    esac
done

# =========== Validation des paramètre et préparation de l’environnement de sortie ========
if [ -z "$INDIR" ]; then
    echo -e "${RED}Erreur: L'option -d/--directory est obligatoire.${NC}"
    echo "Utilisez '$0 --help' pour plus d'informations."
    exit 1
fi

if [[ "$OUTFILE" != *.tsv ]]; then
    echo -e "${RED}Erreur: Le chemin $OUTFILE n'est pas un chemin de fichier .tsv valide.${NC}"
    exit 1 
fi

verifier_prerequis

OUTDIR=$(dirname "$OUTFILE")
if [ "$OUTDIR" = "." ]; then
    OUTDIR=$(pwd)
fi
OUTDIR=$(realpath "$OUTDIR")

LOG_FILE="$OUTDIR/inversions_calling.log"
[ -f "$LOG_FILE" ] && rm "$LOG_FILE"

verifier_dossier "$INDIR" "d'entrée"
INDIR=$(realpath "$INDIR")

creer_dossier "$OUTDIR"

OUTFILE=$(realpath "$OUTFILE")
[ -f "$OUTFILE" ] && rm "$OUTFILE"

TMP_DIR="$OUTDIR/tmp_dir_output"
creer_dossier "$TMP_DIR"
[ "$(ls -A "$TMP_DIR" 2>/dev/null)" ] && rm "$TMP_DIR"/*
# ======================================================================================



# ===== Récupération et vérification de la liste des fichiers FASTA dans le dossier =====
fasta_files=($(find "$INDIR" -maxdepth 1 -type f -iname "*.fasta"))
if [ ${#fasta_files[@]} -eq 0 ]; then
    echo -e "${RED}Erreur: Aucun fichier .fasta trouvé dans le dossier $INDIR.${NC}"
    journaliser "ERREUR: Aucun fichier .fasta trouvé dans le dossier $INDIR"
    exit 1
fi
# ======================================================================================

# Affiche un résumé des paramètres
echo -e "${BLUE}====== Paramètres de l'analyse ======${NC}"
echo -e "Dossier d'entrée:            ${GREEN}$INDIR${NC}"
echo -e "Fichier de sortie:           ${GREEN}$OUTFILE${NC}"
echo -e "Identité de séquence minimale:           ${GREEN}$MIN_IDENTITY${NC}"
echo -e "Threads:                     ${GREEN}$THREADS${NC}"
echo -e "Fichiers FASTA:              ${GREEN}${#fasta_files[@]}${NC}\n"

# Calcul du nombre total de comparaisons
total_comparisons=$(( ${#fasta_files[@]} * (${#fasta_files[@]} - 1) ))
current_comparison=0

# Boucle principale pour comparer toutes les paires de génomes
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
            echo -e "\n${BLUE}[$percent%] Analyse:${NC} ${GREEN}$target_name${NC} - ${GREEN}$query_name${NC}"
            journaliser "\n\t======= Analyse de la paire: $target_name - $query_name ($percent%) ======="
            
            # Extraire les chromosomes pour chaque fichier et garder les communs
            chrs_query=$(grep -o '^>[^ ]*' "$query_file" | sed 's/^>//' | sort | uniq)
            chrs_target=$(grep -o '^>[^ ]*' "$target_file" | sed 's/^>//' | sort | uniq)
            common_chrs=$(comm -12 <(echo "$chrs_target") <(echo "$chrs_query"))
            
            # Vérifier s'il y a des chromosomes communs
            if [ -z "$common_chrs" ]; then
                echo -e "  ${YELLOW}⚠ Aucun chromosome commun trouvé${NC}"
                journaliser "AVERTISSEMENT: Aucun chromosome commun entre $target_name et $query_name"
                continue
            fi

            # Affichage de la progression pour les chromosomes
            num_chrs=$(echo "$common_chrs" | wc -l)
            chr_idx=0            
            # Boucler sur les chromosomes communs
            for chr in $common_chrs; do
                chr_idx=$((chr_idx + 1))
                echo -ne "\r\033[K[${chr_idx}/${num_chrs}] Traitement ${GREEN}$chr${NC}..."
                journaliser " --> Traitement du chromosome: $chr ($target_name - $query_name)"
                
                # Noms de fichiers temporaires
                chr_target_file="$TMP_DIR/target.fasta"
                chr_query_file="$TMP_DIR/query.fasta"
                output_bam="$TMP_DIR/algnt.bam"
                syri_out="$TMP_DIR/syri.out"

                # ==================== Extraire les séquences du chromosome commun ====================
                AWK_EXTRACT="awk -v chr=\">$chr\" '
                    \$0 ~ chr\"([ \\t]|$)\" {print_flag=1} 
                    \$0 ~ \"^>\" && \$0 !~ chr\"([ \\t]|$)\" {print_flag=0} 
                    print_flag
                ' \"$target_file\" > \"$chr_target_file\""
                if ! executer_commande "Extraction du chromosome $chr de $target_name" "$AWK_EXTRACT"; then
                    echo -e "${YELLOW}✗${NC}"
                    continue
                fi
                
                AWK_EXTRACT="awk -v chr=\">$chr\" '
                    \$0 ~ chr\"([ \\t]|$)\" {print_flag=1} 
                    \$0 ~ \"^>\" && \$0 !~ chr\"([ \\t]|$)\" {print_flag=0} 
                    print_flag
                ' \"$query_file\" > \"$chr_query_file\""
                if ! executer_commande "Extraction du chromosome $chr de $query_name" "$AWK_EXTRACT"; then
                    echo -e "${YELLOW}✗${NC}"
                    continue
                fi

                # Vérifier que les fichiers extraits ne sont pas vides
                if [ ! -s "$chr_target_file" ] || [ ! -s "$chr_query_file" ]; then
                    echo -e "${YELLOW}✗${NC} (fichier vide)"
                    journaliser "AVERTISSEMENT: Fichier d'extraction vide pour le chromosome $chr"
                    continue
                fi
                # =====================================================================================


                # ========= Exécuter minimap2 pour mapper les chr et filtrer les résultats =============
                MINIMAP_CMD="minimap2 -ax asm5 -t $THREADS --cs --eqx \"$chr_target_file\" \"$chr_query_file\" 2>/dev/null | 
                samtools view -h | 
                awk -v min_identity=\"$MIN_IDENTITY\" 'BEGIN {OFS=\"\\t\"} 
                    /^@/ {print; next} 
                    {                        
                        matches = mismatches = indels = 0;
                        cigar = \$6;
                        while (match(cigar, /([0-9]+)([=XID])/)) {
                            val = substr(cigar, RSTART, RLENGTH-1);
                            type = substr(cigar, RSTART+RLENGTH-1, 1);
                            if (type == \"=\") matches += val;
                            else if (type == \"X\") mismatches += val;
                            else if (type == \"I\" || type == \"D\") indels += val;
                            cigar = substr(cigar, RSTART+RLENGTH);
                        }
                        
                        aligned_bases = matches + mismatches + indels;
                        if (aligned_bases > 0) {
                            identity = matches / aligned_bases;
                            if (identity >= min_identity) print;
                        }
                    }' | 
                samtools sort -O BAM -o \"$output_bam\"  && 
                samtools index \"$output_bam\""
                if ! executer_commande "Alignement des $chr avec minimap2" "$MINIMAP_CMD"; then
                    echo -e "${YELLOW}✗${NC}"
                    continue
                fi
                # =====================================================================================


                # ========== Exécuter syri pour la détection des variants structuraux =================
                SYRI_CMD="syri -c \"$output_bam\" -r \"$chr_target_file\" -q \"$chr_query_file\" -F B --dir \"$TMP_DIR\" --nc $THREADS > /dev/null 2>&1"
                if ! executer_commande "Détection des variants structuraux avec SyRI" "$SYRI_CMD" true; then
                    echo -e "${YELLOW}✗${NC}"
                    echo "Erreur pour $target_name - $query_name : $chr" >> "$LOG_FILE"
                    continue
                fi
                # Vérifier si le fichier de sortie de SyRI existe
                if [ ! -f "$syri_out" ]; then
                    echo -e "${YELLOW}✗${NC} (fichier manquant)"
                    echo "Fichier de sortie SyRI manquant pour $target_name - $query_name : $chr" >> "$LOG_FILE"
                    continue
                fi
                # =====================================================================================


                # ==================== Extraction et formatage des inversions ===========================
                AWK_EXTRACT="awk -F'\\t' -v OFS='\\t' -v tgt_name=\"$target_name\" -v qry_name=\"$query_name\" '
                (\$9 ~ /INV/) {
                    tstart = (\$2 < \$3) ? \$2 : \$3;
                    tend = (\$2 > \$3) ? \$2 : \$3;
                    qstart = (\$7 < \$8) ? \$7 : \$8;
                    qend = (\$7 > \$8) ? \$7 : \$8;
                    pair = tstart \"-\" tend \"-\" qstart \"-\" qend;
                    sep = \"-\"
                    if (!(pair in seen)) {
                        seen[pair] = 1;
                        tgt_size = tend - tstart;
                        qry_size = qend - qstart;
                        print tgt_name, \$1, tstart, tend, tgt_size, sep, qry_name, \$6, qstart, qend, qry_size
                    }
                }' \"$syri_out\" >> \"$OUTFILE\""
                if ! executer_commande "Extraction des inversions pour $chr" "$AWK_EXTRACT" true; then
                    echo -e "${YELLOW}✗${NC}"
                    continue
                fi
                # =====================================================================================

            done
        fi
    done
done

rm -r $TMP_DIR

journaliser "Analyse terminée. Résultats disponibles dans '$OUTFILE'"
echo -e "\n\n${BLUE}Analyse terminée.${NC}\nRésultats disponibles dans : ${GREEN}'$OUTFILE'${NC}"
exit 0