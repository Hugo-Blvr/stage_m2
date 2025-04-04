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
    echo -e "${BLUE}Usage:${NC} $0 -d <dossier> -m <motif> [options]\n"
    echo -e "${BLUE}Options obligatoires:${NC}"
    echo "  -d, --dossier <dossier>  Chemin vers le dossier contenant les fichiers FASTA (.fasta ou .fa)"
    echo "  -m, --motif <motif>      Motif ADN à rechercher"
    echo -e "${BLUE}Options facultatives:${NC}"
    echo "  -h, --help               Affiche ce message d'aide"
    echo "  -o, --output <dossier>   Spécifie le dossier de sortie (par défaut: ./out_telomere_region)"
    echo "  -v, --visu <visu_file>   Nom du fichier de sortie pour la visualisation"
    echo -e "\n${YELLOW}IMPORTANT: ${NC}Le script python 'visu_telo.py' doit être prèsent dans le même dossier que $0"
    echo -e " si l'option -v/--visu est utilisée.\n"
    exit 0
}

verifier_prerequis() {
    local tools=(tidk find)
    [[ -n "$VISU_FILE_NAME" ]] && tools+=(python)  # Ajout de python si VISU_FILE_NAME est défini

    local missing_tools=()
    for tool in "${tools[@]}"; do
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
    
    echo -e "${BLUE}$description...${NC}"
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
SCRIPT_DIR=$(dirname "$(realpath "$0")")
INDIR=""
MOTIF=""
OUTDIR="./out_telomere_region"
VISU_FILE_NAME=""

# Traitement des options
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -h|--help)
            help
            ;;
        -d|--dossier)
            INDIR="$2"
            shift 2
            ;;
        -m|--motif)
            MOTIF="$2"
            shift 2
            ;;
        -o|--output)
            OUTDIR="$2"
            shift 2
            ;;
        -v|--visu)
            VISU_FILE_NAME="$(realpath "$2")"
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
if [ -z "$INDIR" ] || [ -z "$MOTIF" ]; then
    echo -e "${RED}Erreur: Les options -d/--dossier et -m/--motif sont obligatoires.${NC}"
    echo "Utilisez '$0 --help' pour plus d'informations."
    exit 1
fi

verifier_prerequis
LOG_FILE="$OUTDIR/telo_search.log"
[ -f "$LOG_FILE" ] && rm "$LOG_FILE"

verifier_dossier "$INDIR" "d'entrée"
INDIR=$(realpath "$INDIR")

OUTDIR=$(realpath "$OUTDIR")
creer_dossier "$OUTDIR"

MOTIF_TRIPLE="${MOTIF}${MOTIF}${MOTIF}"

# Vérifie l'existence du script Python si la visualisation est demandée + vérifie l'extention de la visualisation
PYTHON_SCRIPT="$SCRIPT_DIR/visu_telo.py"
if [ -n "$VISU_FILE_NAME" ]; then
    if [ ! -f "$PYTHON_SCRIPT" ]; then
        echo -e "${RED}Erreur: Le script Python '$PYTHON_SCRIPT' n'existe pas.${NC}"
        echo "Ce script est nécessaire pour l'option de visualisation."
        journaliser "ERREUR: Script Python '$PYTHON_SCRIPT' introuvable"
        exit 1
    fi
    if [ ! -r "$PYTHON_SCRIPT" ]; then
        echo -e "${RED}Erreur: Le script Python '$PYTHON_SCRIPT' n'est pas lisible.${NC}"
        journaliser "ERREUR: Script Python '$PYTHON_SCRIPT' non lisible"
        exit 1
    fi

    if [[ "$VISU_FILE_NAME" != *.html ]]; then
        VISU_FILE_NAME="${VISU_FILE_NAME}.html"
    fi

    if ls "$OUTDIR"/*.tsv 1> /dev/null 2>&1; then
        echo -e "${RED}ERREUR : Des fichiers .tsv existent dans $OUTDIR.${NC}"
        journaliser "ERREUR : Des fichiers .tsv existent dans $OUTDIR."
        echo "Veuillez les supprimer ou choisir un autre dossier pour ne pas compromettre la visualisation."
        exit 1
    fi
fi
# ======================================================================================



# ===== Récupération et vérification de la liste des fichiers FASTA dans le dossier =====
fasta_files=($(find "$INDIR" -maxdepth 1 -type f -iname "*.fasta"))
if [ ${#fasta_files[@]} -eq 0 ]; then
    echo -e "${RED}Erreur: Aucun fichier .fasta trouvé dans le dossier $INDIR.${NC}"
    journaliser "ERREUR: Aucun fichier .fasta trouvé dans le dossier $INDIR"
    exit 1
fi
# ======================================================================================

# Affiche un résumé des opérations à effectuer
echo -e "${BLUE}=== Résumé des opérations ===${NC}"
echo -e "Dossier d'entrée:    ${GREEN}$INDIR${NC}"
echo -e "Dossier de sortie:   ${GREEN}$OUTDIR${NC}"
echo -e "Motif rechercher:    ${GREEN}$MOTIF${NC}"
echo -e "Fichiers FASTA:      ${GREEN}${#fasta_files[@]}${NC}"
if [ -n "$VISU_FILE_NAME" ]; then
    echo -e "Fichier html:        ${GREEN}$VISU_FILE_NAME${NC}\n"
else
    echo
fi


# Initialise le compteur
FICHIERS_TRAITES=0

# Boucle sur chaque fichier FASTA
echo -e "${BLUE}Démarrage du traitement des fichiers...${NC}"

for FASTA_FILE in "${fasta_files[@]}"; do
    echo
    FILENAME=$(basename "$FASTA_FILE" | cut -d. -f1)  # Récupère le nom du fichier
    echo -e "${NC}Traitement du fichier: ${GREEN}$FILENAME${NC} (${YELLOW}$((++FICHIERS_TRAITES))/${#fasta_files[@]}${NC})"
    journaliser "Traitement de '$FILENAME'"

    # Vérifie que le fichier est lisible
    if [ ! -r "$FASTA_FILE" ]; then
        echo -e "${RED}Avertissement: Le fichier '$FASTA_FILE' n'est pas lisible. Ignoré.${NC}"
        journaliser "ERREUR: Fichier non lisible: '$FASTA_FILE'"
        continue
    fi
    # Vérifie que le fichier n'est pas vide
    if [ ! -s "$FASTA_FILE" ]; then
        echo -e "${YELLOW}Avertissement: Le fichier '$FASTA_FILE' est vide. Ignoré.${NC}"
        journaliser "AVERTISSEMENT: Fichier vide: '$FASTA_FILE'"
        continue
    fi

    # Exécute la commande tidk search avec le motif triplé
    TIDK_CMD="tidk search -s \"$MOTIF_TRIPLE\" -o \"${FILENAME}_${MOTIF}\" -d \"$OUTDIR\" \"$FASTA_FILE\""
    if executer_commande "Analyse TIDK de '$FILENAME'" "$TIDK_CMD" true; then
        RESULT_FILE="$OUTDIR/${FILENAME}_${MOTIF}_telomeric_repeat_windows.tsv"
    else
        journaliser "ERREUR: Échec de l'analyse de '$FILENAME'"
    fi

done

journaliser "Analyse terminée. Résultats disponibles dans '$OUTDIR'"
echo -e "\n${BLUE}Analyse terminée.${NC}\nRésultats disponibles dans : ${GREEN}'$OUTDIR'${NC}"

# Si un fichier de visualisation est spécifié, lance le script de visualisation
if [ -n "$VISU_FILE_NAME" ]; then
    echo
    VISU_CMD="python \"$PYTHON_SCRIPT\" --input_folder \"$OUTDIR\" --output_file \"$VISU_FILE_NAME\""
    if ! executer_commande "Génération de la visualisation" "$VISU_CMD"; then
        echo -e "${RED}Erreur lors de la génération de la visualisation${NC}"
    else
        journaliser "Visualisation générée avec succès dans: '$VISU_FILE_NAME'"
        echo -e "Visualisation générée avec succès dans: ${GREEN}'$VISU_FILE_NAME'${NC}"
    fi
fi

exit 0