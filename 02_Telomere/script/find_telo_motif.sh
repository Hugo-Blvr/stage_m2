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
    echo "  -o, --output <dossier>   Spécifie le dossier de sortie (par défaut: ./out_tidk)"
    echo
    echo -e "${BLUE}Exemple:${NC}"
    echo "  $0 -d ./data -m ACGT"
    echo "  $0 --dossier ./data --motif ACGT --output ./resultats"
    echo
    exit 0
}

# Fonction pour vérifier les prérequis
verifier_prerequis() {
    # Vérifie si tidk est installé
    if ! command -v tidk &> /dev/null; then
        echo -e "${RED}Erreur: L'outil 'tidk' n'est pas installé ou n'est pas dans le PATH.${NC}"
        echo "Veuillez installer tidk avant d'exécuter ce script."
        exit 1
    fi
}

# Fonction pour journaliser les actions
journaliser() {
    echo -e "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "./tidk_search.log"
}

# Initialisation des variables par défaut
DOSSIER=""
MOTIF=""
DOSSIER_SORTIE="./out_tidk"

# Traitement des options
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -h|--help)
            help
            ;;
        -d|--dossier)
            DOSSIER="$2"
            shift 2
            ;;
        -m|--motif)
            MOTIF="$2"
            shift 2
            ;;
        -o|--output)
            DOSSIER_SORTIE="$2"
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

# Vérifie les prérequis
verifier_prerequis

# Vérifie si le dossier d'entrée existe et est lisible
if [ ! -d "$DOSSIER" ]; then
    echo -e "${RED}Erreur: Le dossier '$DOSSIER' n'existe pas.${NC}"
    exit 1
fi

if [ ! -r "$DOSSIER" ]; then
    echo -e "${RED}Erreur: Le dossier '$DOSSIER' n'est pas lisible.${NC}"
    exit 1
fi

# Triple le motif ADN
MOTIF_TRIPLE="${MOTIF}${MOTIF}${MOTIF}"
journaliser "Motif original: $MOTIF, Motif triplé: $MOTIF_TRIPLE"

# Crée le dossier de sortie s'il n'existe pas
if [ ! -d "$DOSSIER_SORTIE" ]; then
    journaliser "Création du dossier de sortie '$DOSSIER_SORTIE'"
    mkdir -p "$DOSSIER_SORTIE"
fi

# Vérifie que le dossier de sortie est accessible en écriture
if [ ! -w "$DOSSIER_SORTIE" ]; then
    echo -e "${RED}Erreur: Le dossier de sortie '$DOSSIER_SORTIE' n'est pas accessible en écriture.${NC}"
    exit 1
fi

# Stocke la liste des fichiers FASTA dans une variable
# Utilisation de find avec -iname pour être insensible à la casse
FASTA_FILES=$(find "$DOSSIER" -type f \( -iname "*.fasta" -o -iname "*.fa" \))

# Vérifie si des fichiers FASTA ont été trouvés
if [ -z "$FASTA_FILES" ]; then
    echo -e "${RED}Erreur: Aucun fichier FASTA (.fasta ou .fa) trouvé dans '$DOSSIER'.${NC}"
    exit 1
fi

# Compte le nombre de fichiers FASTA trouvés
NOMBRE_FICHIERS=$(echo "$FASTA_FILES" | wc -l)
journaliser "Nombre de fichiers FASTA trouvés: $NOMBRE_FICHIERS"

# Affiche un résumé des opérations à effectuer
echo
echo -e "${BLUE}=== Résumé des opérations ===${NC}"
echo -e "Dossier d'entrée:    ${GREEN}$DOSSIER${NC}"
echo -e "Motif à rechercher:  ${GREEN}$MOTIF${NC} (sera triplé en ${GREEN}$MOTIF_TRIPLE${NC})"
echo -e "Nombre de fichiers:  ${GREEN}$NOMBRE_FICHIERS${NC}"
echo -e "Dossier de sortie:   ${GREEN}$DOSSIER_SORTIE${NC}"
echo -e "${BLUE}===========================${NC}"
echo

# Initialise les compteurs
FICHIERS_TRAITES=0
FICHIERS_AVEC_RESULTATS=0

# Boucle sur chaque fichier FASTA
echo -e "${BLUE}Démarrage du traitement des fichiers...${NC}"
for FASTA_FILE in $FASTA_FILES; do
    FILENAME=$(basename "$FASTA_FILE")  # Récupère juste le nom du fichier
    ISOLAT="${FILENAME%%-*}"            # Extrait la partie avant le premier '-'
    
    echo -e "Traitement du fichier ${GREEN}$FILENAME${NC} (${YELLOW}$((++FICHIERS_TRAITES))/$NOMBRE_FICHIERS${NC})"
    journaliser "Traitement de '$FASTA_FILE' avec motif '$MOTIF_TRIPLE'"

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
    # La commande est exécutée dans un sous-shell pour capturer sa sortie et son code de retour
    OUTPUT=$(tidk search -s "$MOTIF_TRIPLE" -o "${ISOLAT}_${MOTIF}" -d "$DOSSIER_SORTIE" "$FASTA_FILE" 2>&1)
    TIDK_STATUS=$?
    
    if [ $TIDK_STATUS -eq 0 ]; then
        echo -e "  ${GREEN}✓${NC} Analyse terminée"
        
        # Vérifie si des résultats ont été trouvés
        RESULT_FILE="$DOSSIER_SORTIE/${ISOLAT}_${MOTIF}.csv"
        if [ -f "$RESULT_FILE" ] && [ -s "$RESULT_FILE" ]; then
            NOMBRE_RESULTATS=$(grep -c . "$RESULT_FILE" 2>/dev/null || echo "0")
            ((NOMBRE_RESULTATS--))  # Soustrait 1 pour l'en-tête
            if [ "$NOMBRE_RESULTATS" -gt 0 ]; then
                echo -e "  ${GREEN}✓${NC} $NOMBRE_RESULTATS occurrence(s) du motif triplé trouvée(s)"
                ((FICHIERS_AVEC_RESULTATS++))
            else
                echo -e "  ${YELLOW}!${NC} Aucune occurrence du motif triplé trouvée"
            fi
        else
            echo -e "  ${YELLOW}!${NC} Aucun résultat généré"
        fi
    else
        echo -e "  ${RED}✗${NC} Erreur lors de l'analyse"
        echo -e "  ${RED}Détails: $OUTPUT${NC}"
        journaliser "ERREUR: tidk search a échoué pour '$FASTA_FILE': $OUTPUT"
    fi
    
    echo
done

# Affiche un récapitulatif
echo -e "${BLUE}=== Récapitulatif ===${NC}"
echo -e "Fichiers traités: ${GREEN}$FICHIERS_TRAITES/${NOMBRE_FICHIERS}${NC}"
echo -e "Fichiers avec occurrences du motif: ${GREEN}$FICHIERS_AVEC_RESULTATS/${FICHIERS_TRAITES}${NC}"
echo -e "Résultats disponibles dans: ${GREEN}$DOSSIER_SORTIE/${NC}"
echo -e "${BLUE}===================${NC}"

journaliser "Traitement terminé. $FICHIERS_TRAITES fichiers traités, $FICHIERS_AVEC_RESULTATS fichiers avec occurrences du motif."

exit 0
