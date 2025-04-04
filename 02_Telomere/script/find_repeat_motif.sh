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
    echo
    echo -e "${BLUE}Usage:${NC} $0 -g <genome.fasta> -l <longreads.fastq> [options]\n"
    echo -e "${BLUE}Options obligatoires:${NC}"
    echo "  -g, --genome <fichier>    Fichier FASTA du génome de référence"
    echo "  -l, --longreads <fichier> Fichier FASTQ des lectures longues (ONT)"
    echo -e "${BLUE}Options facultatives:${NC}"
    echo "  -h, --help                Affiche ce message d'aide"
    echo "  -o, --output <dossier>    Dossier de sortie des résultats (par défaut: dossier courant)"
    echo "  -t, --threads <nombre>    Nombre de threads à utiliser (par défaut: 8)"
    echo -e "\n${YELLOW}IMPORTANT: ${NC}Le script python 'creat_list_motif.py' doit être prèsent dans le même dossier que $0\n"
    exit 0
}

verifier_prerequis() {
    local missing_tools=()
    for tool in minimap2 samtools trf python3; do
        command -v "$tool" &>/dev/null || missing_tools+=("$tool")
    done

    if ((${#missing_tools[@]})); then
        echo -e "\n${RED}Erreur: Les outils suivants ne sont pas installés ou ne sont pas dans le PATH:${NC}"
        printf "${RED}  - %s\n${NC}" "${missing_tools[@]}"
        echo -e "Veuillez installer ces outils avant d'exécuter ce script.\n"
        exit 1
    fi
}

verifier_fichier() {
    local fichier="$1"
    local description="$2"

    if [ ! -f "$fichier" ] || [ ! -r "$fichier" ]; then
        echo -e "${RED}Erreur: Le fichier $description '$fichier' est introuvable ou non lisible.${NC}"
        journaliser "ERREUR: Fichier '$fichier' introuvable ou non lisible"
        exit 1
    elif [ ! -s "$fichier" ]; then
        echo -e "${YELLOW}Avertissement: Le fichier $description '$fichier' est vide.${NC}"
        journaliser "AVERTISSEMENT: Fichier '$fichier' vide"
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
    
    echo -e "${BLUE}$description${NC}"
    journaliser "EXÉCUTION: $description" 
    
    OUTPUT=$(eval "$commande" 2>&1) # Exécution de la commande dans un sous-shell et récupération de la sortie
    local STATUS=$?  # Récupération du code de retour de la commande
    
    if [ $STATUS -eq 0 ] || [ "$ignore_error" == "true" ] || { [[ "$description" == *"TRF"* ]] && [ -f "$DAT_FILE_START" ]; }; then
        return 0  # si la commande s'est exécutée avec succès ou si on ignore les erreurs
    else
        # Si une erreur s'est produite, afficher un message d'erreur et journaliser
        echo -e "  ${RED}✗${NC} Erreur lors de $description"
        journaliser "ERREUR: $description: $OUTPUT"
        return 1
    fi
}

# Initialisation des variables par défaut
SCRIPT_DIR=$(dirname "$(realpath "$0")")
GENOME=""
LONGREADS=""
OUTDIR="."
THREADS=8

# Traitement des options
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -h|--help)
            help
            ;;
        -g|--genome)
            GENOME="$2"
            shift 2
            ;;
        -l|--longreads)
            LONGREADS="$2"
            shift 2
            ;;
        -o|--output)
            OUTDIR="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
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
if [ -z "$GENOME" ] || [ -z "$LONGREADS" ]; then
    echo -e "${RED}Erreur: Les options -g/--genome et -l/--longreads sont obligatoires.${NC}"
    echo "Utilisez '$0 --help' pour plus d'informations."
    exit 1
fi

verifier_prerequis
verifier_fichier "$GENOME" "du génome"
verifier_fichier "$LONGREADS" "des lectures longues"
GENOME=$(realpath "$GENOME")
LONGREADS=$(realpath "$LONGREADS")

# Définition du fichier de log dans le dossier de sortie
LOG_FILE="$OUTDIR/telomere_analysis.log"
[ -f "$LOG_FILE" ] && rm "$LOG_FILE"

# Création des répertoires de sortie
creer_dossier "$OUTDIR"
OUTDIR=$(realpath "$OUTDIR")
BAM_OUTPUT_DIR="$OUTDIR/bam"
SOFT_CLIP_OUTPUT_DIR="$OUTDIR/soft_clip"
TRF_OUTPUT_DIR="$OUTDIR/trf"
creer_dossier "$BAM_OUTPUT_DIR"
creer_dossier "$SOFT_CLIP_OUTPUT_DIR" 
creer_dossier "$TRF_OUTPUT_DIR"

# Définition des fichiers
GENOME_BASE=$(basename "${GENOME}" | cut -d. -f1)
MMI_FILE="$BAM_OUTPUT_DIR/${GENOME_BASE}.mmi"
SORTED_BAM="$BAM_OUTPUT_DIR/${GENOME_BASE}_mapping.sorted.bam"
CHR_LENGTHS_FILE="$OUTDIR/${GENOME_BASE}_chr_lengths.tsv"

# Fichiers soft-clips
SOFT_CLIP_START_TSV="$SOFT_CLIP_OUTPUT_DIR/${GENOME_BASE}_softclip_start.tsv"
SOFT_CLIP_END_TSV="$SOFT_CLIP_OUTPUT_DIR/${GENOME_BASE}_softclip_end.tsv"
SOFT_CLIP_START_FASTA="$SOFT_CLIP_OUTPUT_DIR/${GENOME_BASE}_softclip_start.fasta"
SOFT_CLIP_END_FASTA="$SOFT_CLIP_OUTPUT_DIR/${GENOME_BASE}_softclip_end.fasta"

# Fichiers TRF
DAT_FILE_START="$TRF_OUTPUT_DIR/${GENOME_BASE}_softclip_start.fasta.2.7.7.80.4.10.500.dat"
DAT_FILE_END="$TRF_OUTPUT_DIR/${GENOME_BASE}_softclip_end.fasta.2.7.7.80.4.10.500.dat"
FILTERED_REPEATS_START="$TRF_OUTPUT_DIR/${GENOME_BASE}_filtered_repeats_start.fasta"
FILTERED_REPEATS_END="$TRF_OUTPUT_DIR/${GENOME_BASE}_filtered_repeats_end.fasta"

# Fichiers de sortie
OUT_NAME_START="$OUTDIR/${GENOME_BASE}_list_telomere_start.tsv"
OUT_NAME_END="$OUTDIR/${GENOME_BASE}_list_telomere_end.tsv"
# ======================================================================================


# Affiche un résumé des paramètres
echo -e "${BLUE}====== Paramètres de l'analyse ======${NC}"
echo -e "Génome:                 ${GREEN}$GENOME${NC}"
echo -e "Lectures longues:       ${GREEN}$LONGREADS${NC}"
echo -e "Dossier de sortie:      ${GREEN}$OUTDIR${NC}"
echo -e "Threads:                ${GREEN}$THREADS${NC}\n"


############################# Étape 1 : Alignement des long reads ############################
if [ ! -f "$MMI_FILE" ]; then
    executer_commande "Indexation du génome" "minimap2 -d \"$MMI_FILE\" \"$GENOME\""
fi
if [ ! -f "$SORTED_BAM" ]; then
    executer_commande "Alignement des lectures longues" \
        "minimap2 -t $THREADS -ax map-ont \"$MMI_FILE\" \"$LONGREADS\" | \
        samtools view -bS | samtools sort -@ $THREADS -o \"$SORTED_BAM\""
fi
if [ ! -f "${SORTED_BAM}.bai" ]; then
    executer_commande "Indexation du BAM" "samtools index \"$SORTED_BAM\""
fi
# Extraction des longueurs des chromosomes
if [ ! -f "$CHR_LENGTHS_FILE" ]; then
    executer_commande "Extraction des longueurs des chromosomes" \
        "samtools idxstats \"$SORTED_BAM\" > \"$CHR_LENGTHS_FILE\""
fi
################################################################################################################


############################# Étape 2 : Extraction des soft-clips d'intérêt ############################
echo
# Vérification si les fichiers de sortie existent
if [ ! -f "$SOFT_CLIP_START_TSV" ] || [ ! -f "$SOFT_CLIP_END_TSV" ]; then
    # Construction de la commande AWK avec le script en ligne
    AWK_CMD="samtools view \"$SORTED_BAM\" | awk -v OFS='\t' -v chr_lengths_file=\"$CHR_LENGTHS_FILE\" -v soft_clip_start_tsv=\"$SOFT_CLIP_START_TSV\" -v soft_clip_end_tsv=\"$SOFT_CLIP_END_TSV\" '
    BEGIN {
        while ((getline < chr_lengths_file) > 0) { chr_lengths[\$1] = \$2 }
        close(chr_lengths_file)
    }
    {
        if (\$10 == \"*\") next  # Ignorer séquences vides
        chr = \$3; start = \$4
        if (!(chr in chr_lengths)) next
        ref_length = chr_lengths[chr]
        
        if (match(\$6, /^([0-9]+)S/, m)) { soft_clip = m[1] } else { soft_clip = 0 }
        
        cigar = \$6; mapped_length = 0
        while (match(cigar, /([0-9]+)([MID])/)) {
            len = substr(cigar, RSTART, RLENGTH - 1)
            op = substr(cigar, RSTART + RLENGTH - 1, 1)
            if (op == \"M\" || op == \"D\") mapped_length += len
            cigar = substr(cigar, RSTART + RLENGTH)
        }
        
        end = start + mapped_length - 1
        min_start = ref_length - 8000
        min_end = ref_length - 1000
        
        if (start <= 1000 && end >= 8000 && soft_clip >= 50) {
            print \$1, chr, start, end, \$10 >> soft_clip_start_tsv
        }
        
        if (start <= min_start && end >= min_end && soft_clip >= 50) {
            print \$1, chr, start, end, \$10 >> soft_clip_end_tsv
        }
    }'"
    
    # Exécution avec la fonction executer_commande
    if ! executer_commande "Extraction des soft-clips" "$AWK_CMD" false; then
        exit 1
    fi
    
    # Vérification supplémentaire
    if [ -f "$SOFT_CLIP_START_TSV" ] && [ -f "$SOFT_CLIP_END_TSV" ]; then
        journaliser "SUCCÈS: Extraction des soft-clips"
    else
        journaliser "ERREUR: Extraction des soft-clips échouée"
        exit 1
    fi
fi

# Formatage en fichiers FASTA
if [ ! -f "$SOFT_CLIP_START_FASTA" ]; then
    executer_commande "Formatage des soft-clips (début)" "awk -F '\t' '{print \">\"\$2\":\"\$1\":\"\$3\"-\"\$4\"\\n\"\$5}' \"$SOFT_CLIP_START_TSV\" > \"$SOFT_CLIP_START_FASTA\""
fi
if [ ! -f "$SOFT_CLIP_END_FASTA" ]; then
    executer_commande "Formatage des soft-clips (fin)" "awk -F '\t' '{print \">\"\$2\":\"\$1\":\"\$3\"-\"\$4\"\\n\"\$5}' \"$SOFT_CLIP_END_TSV\" > \"$SOFT_CLIP_END_FASTA\""
fi

# Vérification de la présence de séquences
for fasta_file in "$SOFT_CLIP_START_FASTA" "$SOFT_CLIP_END_FASTA"; do
    if [ ! -s "$fasta_file" ]; then
        echo -e "${YELLOW}Avertissement: Le fichier '$fasta_file' est vide. Aucune soft-clip n'a été trouvée.${NC}"
        journaliser "AVERTISSEMENT: Fichier '$fasta_file' vide"
    fi
done

rm $CHR_LENGTHS_FILE
################################################################################################################


############################ Étape 3 : Détection des motifs avec TRF ############################
echo
# Exécution de TRF pour les séquences de début
if [ -s "$SOFT_CLIP_START_FASTA" ] && [ ! -f "$DAT_FILE_START" ]; then
    executer_commande "Analyse TRF des soft-clips de début" "(cd \"$TRF_OUTPUT_DIR\" && trf \"$SOFT_CLIP_START_FASTA\" 2 7 7 80 4 10 500 -f -h -d)"
fi
# Exécution de TRF pour les séquences de fin
if [ -s "$SOFT_CLIP_END_FASTA" ] && [ ! -f "$DAT_FILE_END" ]; then
    executer_commande "Analyse TRF des soft-clips de fin" "(cd \"$TRF_OUTPUT_DIR\" && trf \"$SOFT_CLIP_END_FASTA\" 2 7 7 80 4 10 500 -f -h -d)"
fi

# Fonction pour filtrer les répétitions
filter_repeats() {
    local input_file="$1"
    local output_file="$2"
    
    if [ -f "$input_file" ] && [ -s "$input_file" ] && [ ! -f "$output_file" ]; then
        awk 'BEGIN {seq_name=""}
            /^Sequence:/ {seq_name=$2}
            $4 >= 3 {
                print ">" seq_name " repeat:" $1 "-" $2 " motif=" $14 " copies=" $4
                print $15
            }' "$input_file" > "$output_file"
            
        if [ -s "$output_file" ]; then
            journaliser "SUCCÈS: Filtrage des répétitions de '$(basename "$input_file")'"
            return 0
        else
            echo -e "${YELLOW}Avertissement: Aucun motif répétitif trouvé dans $(basename "$input_file").${NC}"
            journaliser "AVERTISSEMENT: Aucun motif répétitif trouvé dans '$(basename "$input_file")'"
            touch "$output_file"  # Création d'un fichier vide pour éviter les erreurs
            return 0  # On considère ce cas comme un succès
        fi
    else
        echo -e "${YELLOW}Avertissement: Le fichier '$(basename "$input_file")' n'existe pas ou est vide.${NC}"
        journaliser "AVERTISSEMENT: Fichier '$(basename "$input_file")' non disponible pour filtrage"
        touch "$output_file"  # Création d'un fichier vide pour éviter les erreurs
        return 0  # On considère ce cas comme un succès
    fi
}

# Filtrage des répétitions
journaliser "EXÉCUTION: Filtrage des motifs répétitifs"
if [ -s "$DAT_FILE_START" ] && [ ! -f "$FILTERED_REPEATS_START" ]; then
    executer_commande "Filtrage des motifs répétitifs (début)" "filter_repeats '$DAT_FILE_START' '$FILTERED_REPEATS_START'"
fi
if [ -s "$DAT_FILE_END" ] && [ ! -f "$FILTERED_REPEATS_END" ]; then
    executer_commande "Filtrage des motifs répétitifs (fin)" "filter_repeats '$DAT_FILE_END' '$FILTERED_REPEATS_END'"
fi
################################################################################################################


############################# Étape 4 : Génération de la liste des télomères ############################
echo
# Vérification de l'existence du script Python
PYTHON_SCRIPT="$SCRIPT_DIR/creat_list_motif.py"
if [ ! -f "$PYTHON_SCRIPT" ]; then
    echo -e "${RED}Erreur: Le script Python '$PYTHON_SCRIPT' n'existe pas.${NC}"
    journaliser "ERREUR: Script Python '$PYTHON_SCRIPT' introuvable"
    exit 1
fi

# Génération des listes de télomères
echo -e "${BLUE}Génération des listes de télomères...${NC}"
journaliser "EXÉCUTION: Génération des listes de télomères"

# Traitement pour les motifs de début
if [ -s "$FILTERED_REPEATS_START" ] && [ ! -f "$OUT_NAME_START" ]; then
    executer_commande "Génération de la liste des télomères de début" "python3 \"$PYTHON_SCRIPT\" \"$FILTERED_REPEATS_START\" \"$OUT_NAME_START\""
elif [ ! -s "$FILTERED_REPEATS_START" ]; then
    echo -e "${YELLOW}Avertissement: Pas de motifs répétitifs de début à traiter.${NC}"
    journaliser "AVERTISSEMENT: Pas de motifs répétitifs de début à traiter"
    touch "$OUT_NAME_START"  # Création d'un fichier vide 
fi

# Traitement pour les motifs de fin
if [ -s "$FILTERED_REPEATS_END" ] && [ ! -f "$OUT_NAME_END" ]; then
    executer_commande "Génération de la liste des télomères de fin" "python3 \"$PYTHON_SCRIPT\" \"$FILTERED_REPEATS_END\" \"$OUT_NAME_END\""
elif [ ! -s "$FILTERED_REPEATS_END" ]; then
    echo -e "${YELLOW}Avertissement: Pas de motifs répétitifs de fin à traiter.${NC}"
    journaliser "AVERTISSEMENT: Pas de motifs répétitifs de fin à traiter"
    touch "$OUT_NAME_END"  # Création d'un fichier vide
fi
################################################################################################################


journaliser "\nAnalyse terminée. Résultats disponibles dans '$OUTDIR'"
echo -e "\n${BLUE}Analyse terminée.${NC}\nRésultats disponibles dans : ${GREEN}'$OUTDIR'${NC}"
exit 0
