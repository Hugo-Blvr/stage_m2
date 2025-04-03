#!/bin/bash
#
# Ce script recherche un motif spécifié triplé dans tous les fichiers FASTA 
# d'un répertoire donné en utilisant l'outil tidk.

# Définition des couleurs pour une meilleure lisibilité
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Fonction d'affichage de l'aide
help() {
    echo
    echo -e "${BLUE}Usage:${NC} $0 -d <dossier> -m <motif> [options]"
    echo
    echo "Description:"
    echo "  Ce script recherche un motif ADN triplé dans tous les fichiers FASTA d'un dossier."
    echo
    echo -e "${BLUE}Options obligatoires:${NC}"
    echo "  -d, --dossier <dossier>  Chemin vers le dossier contenant les fichiers FASTA (.fasta ou .fa)"
    echo "  -m, --motif <motif>      Motif ADN à rechercher (sera triplé automatiquement)"
    echo
    echo -e "${BLUE}Options facultatives:${NC}"
    echo "  -h, --help               Affiche ce message d'aide"
    echo "  -o, --output <dossier>   Spécifie le dossier de sortie (par défaut: ./out_telomere_region)"
    echo "  -v, --visu <visu_file>   Nom du fichier de sortie pour la visualisation"
    echo
    echo -e "${YELLOW}Note:${NC}"
    echo " Le script python 'visu_telo.py' doit être présent dans le même dossier que $0"
    echo " si l'option -v/--visu est utilisée."
    echo
    exit 0
}

# Fonction pour vérifier les prérequis
verifier_prerequis() {
    local missing_tools=()
    
    for tool in tidk find; do
        if ! command -v $tool &> /dev/null; then
            missing_tools+=("$tool")
        fi
    done

    # Vérifie si python est installé uniquement si visu est demandée
    if [ -n "$VISU_FILE_NAME" ] && ! command -v python &> /dev/null; then
        missing_tools+=("python")
    fi
    
    if [ ${#missing_tools[@]} -gt 0 ]; then
        echo -e "${RED}Erreur: Les outils suivants ne sont pas installés ou ne sont pas dans le PATH:${NC}"
        for tool in "${missing_tools[@]}"; do
            echo -e "${RED}  - $tool${NC}"
        done
        echo "Veuillez installer ces outils avant d'exécuter ce script."
        journaliser "ERREUR: Outils manquants: ${missing_tools[*]}"
        exit 1
    fi
}

# Fonction pour journaliser les actions
journaliser() {
    echo -e "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_FILE"
}

# Fonction pour vérifier l'existence et la lisibilité d'un dossier
verifier_dossier() {
    local dossier="$1"
    local description="$2"
    
    if [ ! -d "$dossier" ]; then
        echo -e "${RED}Erreur: Le dossier $description '$dossier' n'existe pas.${NC}"
        journaliser "ERREUR: Dossier '$dossier' introuvable"
        exit 1
    fi
    
    if [ ! -r "$dossier" ]; then
        echo -e "${RED}Erreur: Le dossier $description '$dossier' n'est pas lisible.${NC}"
        journaliser "ERREUR: Dossier '$dossier' non lisible"
        exit 1
    fi
}

# Fonction pour créer et vérifier un dossier
creer_dossier() {
    local dossier="$1"
    
    if [ ! -d "$dossier" ]; then
        journaliser "Création du dossier '$dossier'"
        mkdir -p "$dossier"
        
        if [ $? -ne 0 ]; then
            echo -e "${RED}Erreur: Impossible de créer le dossier '$dossier'.${NC}"
            journaliser "ERREUR: Création du dossier '$dossier' échouée"
            exit 1
        fi
    fi
    
    if [ ! -w "$dossier" ]; then
        echo -e "${RED}Erreur: Le dossier '$dossier' n'est pas accessible en écriture.${NC}"
        journaliser "ERREUR: Dossier '$dossier' non accessible en écriture"
        exit 1
    fi
}

# Fonction pour exécuter une commande avec gestion d'erreur
executer_commande() {
    local description="$1"
    local commande="$2"
    local ignore_error="${3:-false}"
    
    echo -e "${BLUE}$description...${NC}"
    journaliser "EXÉCUTION: $description"
    
    OUTPUT=$(eval "$commande" 2>&1)
    local STATUS=$?
    
    if [ $STATUS -eq 0 ] || [ "$ignore_error" == "true" ]; then
        journaliser "SUCCÈS: $description"
        return 0
    else
        echo -e "  ${RED}✗${NC} Erreur lors de $description"
        journaliser "ERREUR: $description: $OUTPUT"
        if [ "$ignore_error" != "true" ]; then
            return 1
        fi
    fi
}

# Initialisation des variables par défaut
SCRIPT_DIR=$(dirname "$(realpath "$0")")
DOSSIER=""
MOTIF=""
DOSSIER_SORTIE="./out_telomere_region"
VISU_FILE_NAME=""

# Traitement des options
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -h|--help)
            help
            ;;
        -d|--dossier)
            DOSSIER=$(realpath "$2")
            shift 2
            ;;
        -m|--motif)
            MOTIF="$2"
            shift 2
            ;;
        -o|--output)
            DOSSIER_SORTIE=$(realpath "$2")
            shift 2
            ;;
        -v|--visu)
            VISU_FILE_NAME="$2"
            shift 2
            ;;
        *)
            echo -e "${RED}Erreur: Option inconnue '$1'${NC}"
            echo "Utilisez '$0 --help' pour voir les options disponibles."
            exit 1
            ;;
    esac
done

# Vérifie si les arguments obligatoires sont fournis
if [ -z "$DOSSIER" ] || [ -z "$MOTIF" ]; then
    echo -e "${RED}Erreur: Les options -d/--dossier et -m/--motif sont obligatoires.${NC}"
    echo "Utilisez '$0 --help' pour plus d'informations."
    exit 1
fi

# Définition du fichier de log dans le dossier de sortie
creer_dossier "$DOSSIER_SORTIE"
LOG_FILE="$DOSSIER_SORTIE/tidk_search.log"

# Vérifie les prérequis
verifier_prerequis

# Vérifie si le dossier d'entrée existe et est lisible
verifier_dossier "$DOSSIER" "d'entrée"

# Triple le motif ADN
MOTIF_TRIPLE="${MOTIF}${MOTIF}${MOTIF}"
journaliser "Motif original: $MOTIF, Motif triplé: $MOTIF_TRIPLE"

# Stocke la liste des fichiers FASTA dans une variable
FASTA_FILES=$(find "$DOSSIER" -type f \( -iname "*.fasta" -o -iname "*.fa" \))

# Vérifie si des fichiers FASTA ont été trouvés
if [ -z "$FASTA_FILES" ]; then
    echo -e "${RED}Erreur: Aucun fichier FASTA (.fasta ou .fa) trouvé dans '$DOSSIER'.${NC}"
    journaliser "ERREUR: Aucun fichier FASTA trouvé dans '$DOSSIER'"
    exit 1
fi

# Compte le nombre de fichiers FASTA trouvés
NOMBRE_FICHIERS=$(echo "$FASTA_FILES" | wc -l)
journaliser "Nombre de fichiers FASTA trouvés: $NOMBRE_FICHIERS"

# Vérifie l'existence du script Python si la visualisation est demandée
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
fi

# Affiche un résumé des opérations à effectuer
echo
echo -e "${BLUE}=== Résumé des opérations ===${NC}"
echo -e "Dossier d'entrée:      ${GREEN}$DOSSIER${NC}"
echo -e "Motif à rechercher:    ${GREEN}$MOTIF${NC} (sera triplé en ${GREEN}$MOTIF_TRIPLE${NC})"
echo -e "Nombre de fichiers:    ${GREEN}$NOMBRE_FICHIERS${NC}"
echo -e "Dossier de sortie:     ${GREEN}$DOSSIER_SORTIE${NC}"
if [ -n "$VISU_FILE_NAME" ]; then
    echo -e "Fichier visualisation: ${GREEN}$VISU_FILE_NAME${NC}"
fi
echo -e "${BLUE}===========================${NC}"
echo

# Initialise le compteur
FICHIERS_TRAITES=0

# Boucle sur chaque fichier FASTA
echo -e "${BLUE}Démarrage du traitement des fichiers...${NC}"

for FASTA_FILE in $FASTA_FILES; do
    echo
    FILENAME=$(basename "$FASTA_FILE" | cut -d. -f1)  # Récupère juste le nom du fichier

    echo -e "${BLUE}Traitement du fichier ${GREEN}$FILENAME${NC} (${YELLOW}$((++FICHIERS_TRAITES))/$NOMBRE_FICHIERS${NC})"
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
    TIDK_CMD="tidk search -s \"$MOTIF_TRIPLE\" -o \"${FILENAME}_${MOTIF}\" -d \"$DOSSIER_SORTIE\" \"$FASTA_FILE\""
    
    if executer_commande "Analyse TIDK de '$FILENAME'" "$TIDK_CMD" true; then
        # Vérifie si des résultats ont été trouvés
        RESULT_FILE="$DOSSIER_SORTIE/${FILENAME}_${MOTIF}_telomeric_repeat_windows.tsv"
    else
        journaliser "ERREUR: Échec de l'analyse de '$FILENAME'"
    fi
done

# Résumé final
echo
journaliser "Traitement des fichiers terminé."
echo

# Si un fichier de visualisation est spécifié, lance le script de visualisation
if [ -n "$VISU_FILE_NAME" ]; then
    VISU_CMD="python \"$PYTHON_SCRIPT\" --input_folder \"$DOSSIER_SORTIE\" --output_file \"$VISU_FILE_NAME\""
    if ! executer_commande "Génération de la visualisation" "$VISU_CMD"; then
        echo -e "${RED}Erreur lors de la génération de la visualisation${NC}"
    else
        echo -e "${GREEN}Visualisation générée avec succès dans '$VISU_FILE_NAME'${NC}"
    fi
fi

echo 
journaliser "${GREEN}Toutes les opérations sont terminées. Résultats disponibles dans '$DOSSIER_SORTIE'${NC}"

exit 0