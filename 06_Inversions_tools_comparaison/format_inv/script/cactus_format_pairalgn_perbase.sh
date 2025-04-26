#!/bin/bash

# Vérification des arguments
if [ -z "$2" ]; then
    echo "Usage: $0 <input_file> <output_file>"
    exit 1
fi

ARCHIVE_FILE=$1
OUTPUT_FILE=$2
TEMP_DIR=$(mktemp -d)
MERGED_TEMP=$(mktemp)
SORTED_TEMP=$(mktemp)
FINAL_TEMP=$(mktemp)

# Vérification de l'existence du fichier d'archive
if [ ! -f "$ARCHIVE_FILE" ]; then
    echo "Erreur: Le fichier d'archive '$ARCHIVE_FILE' n'existe pas."
    exit 1
fi

echo "Extraction de l'archive..."
tar -xzf "$ARCHIVE_FILE" -C "$TEMP_DIR"

# Créer l'en-tête dans le fichier temporaire
echo -e "iso\tchr\tstart\tend\tmapon" > "$MERGED_TEMP"

echo "Traitement des fichiers CSV..."
# Trouve tous les fichiers CSV dans le répertoire temporaire et les traite
find "$TEMP_DIR" -name "*.csv" -type f | while read csv_file; do
    echo "Traitement de $csv_file"
    
    # Utilise AWK pour transformer les données en se basant sur le 5ème champ comme identifiant de bloc
    awk -F, '
    {
        query = $1;
        target = $2;
        pos = $3;
        chr = $4;
        bloc = $5;
        sv = $6;
        
        if (sv == 1) {
            # Clé unique pour chaque combinaison de query, chr, target et bloc
            key = query "," chr "," target "," bloc;
            
            # Stocker les positions min et max pour chaque bloc
            if (!(key in min_pos) || pos < min_pos[key]) {
                min_pos[key] = pos;
            }
            if (!(key in max_pos) || pos > max_pos[key]) {
                max_pos[key] = pos;
            }
        }
    }
    
    END {
        # Sortir les résultats pour chaque bloc
        for (key in min_pos) {
            split(key, comp, ",");
            query = comp[1];
            chr = comp[2];
            target = comp[3];
            bloc = comp[4];
            
            printf "%s\t%s\t%d\t%d\t%s\n", query, chr, min_pos[key], max_pos[key], target;
        }
    }
    ' "$csv_file" >> "$MERGED_TEMP"
done

echo "Tri du fichier de sortie..."
# Trier le fichier temporaire par iso, chr, mapon, start
(head -n 1 "$MERGED_TEMP"; tail -n +2 "$MERGED_TEMP" | sort -k1,1 -k2,2 -k5,5 -k3,3n) > "$SORTED_TEMP"

echo "Fusion des inversions adjacentes..."
# Conserver l'en-tête
head -n 1 "$SORTED_TEMP" > "$FINAL_TEMP"

# Fusionner les blocs adjacents ou qui se chevauchent avec un écart de 0 ou 1 pb
awk ' BEGIN { FS=OFS="\t" } 
NR == 1 { next } # Skip header 
{
    if (NR == 2) {
        # First data record
        last_iso = $1;
        last_chr = $2;
        last_start = $3;
        last_end = $4;
        last_mapon = $5;
    } else {
        # If same iso, chr and mapon, and positions are adjacent or overlapping
        if ($1 == last_iso && $2 == last_chr && $5 == last_mapon &&
            ($3 <= last_end + 1)) {
            # Merge keeping start of first and end of last
            if ($4 > last_end) last_end = $4;
        } else {
            # Write previous record and start a new one
            print last_iso, last_chr, last_start, last_end, last_mapon;
            last_iso = $1;
            last_chr = $2;
            last_start = $3;
            last_end = $4;
            last_mapon = $5;
        }
    }
} 
END {
    # Write the last record
    print last_iso, last_chr, last_start, last_end, last_mapon;
}' "$SORTED_TEMP" >> "$FINAL_TEMP"


# Copier le résultat final dans le fichier de sortie
cp "$FINAL_TEMP" "$OUTPUT_FILE"


# filtre les inversion de - de 100pb
echo "Filtrage des inversions de taille inférieure à 100pb..."
temp_file=$(mktemp)
awk '
BEGIN {FS = OFS = "\t";}
NR == 1 || ($4 - $3) >= 100 {print;}
' "$OUTPUT_FILE" > "$temp_file"
mv "$temp_file" "$OUTPUT_FILE"

# Nettoyage
echo "Nettoyage des fichiers temporaires..."
rm -rf "$TEMP_DIR" "$MERGED_TEMP" "$SORTED_TEMP" "$FINAL_TEMP"

echo "Conversion terminée. Résultat disponible dans $OUTPUT_FILE"