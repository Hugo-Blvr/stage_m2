#!/bin/bash

# Vérification des arguments
if [ $# -ne 2 ]; then
    echo "Usage: $0 03_all_data_per_nt.tar.gz output.tsv"
    exit 1
fi

ARCHIVE_FILE=$1
OUTPUT_FILE=$2
TEMP_DIR=$(mktemp -d)
MERGED_TEMP=$(mktemp)

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
    
    # Utilise AWK pour transformer les données
    awk -F, '
    NR == 1 {
        for (i=1; i<=NF; i++) {
            if ($i=="query") query_idx=i;
            else if ($i=="target") target_idx=i;
            else if ($i=="pos") pos_idx=i;
            else if ($i=="chr") chr_idx=i;
            else if ($i=="SV") sv_idx=i;
        }
        next;
    }
    $sv_idx == 1 {
        key = $query_idx "," $chr_idx "," $target_idx;
        pos[key][length(pos[key])+1] = $pos_idx;
    }
    
    END {
        for (key in pos) {
            split(key, comp, ",");
            query = comp[1];
            chr = comp[2];
            target = comp[3];
            
            n = asort(pos[key]);
            
            if (n > 0) {
                start = pos[key][1];
                prev = start;
                
                for (i=2; i<=n; i++) {
                    if (pos[key][i] != prev + 1) {
                        printf "%s\t%s\t%d\t%d\t%s\n", query, chr, start, prev, target;
                        start = pos[key][i];
                    }
                    prev = pos[key][i];
                }

                printf "%s\t%s\t%d\t%d\t%s\n", query, chr, start, prev, target;
            }
        }
    }
    ' "$csv_file" >> "$MERGED_TEMP"
done

echo "Tri du fichier de sortie..."
# Trier le fichier temporaire par iso, chr, start (en gardant l'en-tête en première ligne)
(head -n 1 "$MERGED_TEMP"; tail -n +2 "$MERGED_TEMP" | sort -k1,1 -k5,5 -k2,2 -k3,3n) > "$OUTPUT_FILE"

# Nettoyage
echo "Nettoyage des fichiers temporaires..."
rm -rf "$TEMP_DIR" "$MERGED_TEMP"

echo "Conversion terminée. Résultat disponible dans $OUTPUT_FILE"