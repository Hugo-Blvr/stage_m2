#!/bin/bash

# Usage: ./detect_inversions.sh input.maf output.tsv

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input.maf output.tsv"
    exit 1
fi

INPUT_FILE="$1"
OUTPUT_FILE="$2"

# Créer un fichier temporaire pour les résultats avant tri
TMP_FILE=$(mktemp)

# Écrire l'en-tête dans le fichier temporaire
echo -e "iso\tchr\tstart\tend\tmapon" > "$OUTPUT_FILE"
echo -e "iso\tchr\tstart\tend\tmapon" > "$TMP_FILE"

# Traitement principal avec AWK
awk '
BEGIN {
    in_block = 0;
    block_num = 0;
    # Initialiser explicitement le tableau des inversions
    delete inversions;
}

/^a/ {
    if (in_block) {
        process_block();
    }
    in_block = 1;
    block_num++;
    if (block_num % 100 == 0) {
        printf "\rTraitement du bloc numéro %d", block_num > "/dev/stderr";
    }
    delete lines_by_chr;
    delete chrs_in_block;
    chr_count = 0;
    next;
}

/^s/ {
    if (in_block) {
        split($2, name_parts, "\\.");
        isolat = name_parts[1];
        chr = name_parts[2];
        start = $3;
        size = $4;
        strand = $5;

        if (!(chr in chrs_in_block)) {
            chrs_in_block[chr] = ++chr_count;
        }
        
        lines_by_chr[chr][length(lines_by_chr[chr]) + 1] = isolat "," start "," strand "," size;
    }
    next;
}

/^$/ || /^#/ {
    next;
}

END {
    if (in_block) {
        process_block();
    }
    print "" > "/dev/stderr";
    for (inv_key in inversions) {
        print_inversion(inv_key);
    }
}

function print_inversion(key) {
    split(key, parts, ":");
    iso = parts[1];
    chr = parts[2];
    mapon = parts[3];
    start = parts[4];
    end = parts[5];
    
    print iso "\t" chr "\t" start "\t" end "\t" mapon;
}

function process_block() {
    for (chr in chrs_in_block) {
        has_plus = 0;
        has_minus = 0;

        # Vérifier si le bloc contient des brins + et -
        for (i = 1; i <= length(lines_by_chr[chr]); i++) {
            split(lines_by_chr[chr][i], data, ",");
            strand = data[3];
            if (strand == "+") has_plus = 1;
            if (strand == "-") has_minus = 1;
        }
        
        # Si le bloc contient les deux orientations, c"est potentiellement une inversion
        if (has_plus && has_minus) {
            for (i = 1; i <= length(lines_by_chr[chr]); i++) {
                split(lines_by_chr[chr][i], data_i, ",");
                iso_i = data_i[1];
                start_i = data_i[2];
                strand_i = data_i[3];
                size = data_i[4];
                
                for (j = 1; j <= length(lines_by_chr[chr]); j++) {
                    if (i == j) continue;
                    
                    split(lines_by_chr[chr][j], data_j, ",");
                    iso_j = data_j[1];
                    start_j = data_j[2];
                    strand_j = data_j[3];
                    
                    if (strand_i != strand_j) {
                        start_pos = start_i;
                        end_pos = start_i + size - 1;
                        
                        # Créer une inversion et l"ajouter à la liste
                        add_inversion(iso_i, chr, iso_j, start_pos, end_pos);
                    }
                }
            }
        }
    }
}

function add_inversion(iso, chr, mapon, start, end) {
    merged = 0;
    
    # Tenter de fusionner avec des inversions existantes
    for (inv_key in inversions) {
        split(inv_key, parts, ":");
        inv_iso = parts[1];
        inv_chr = parts[2];
        inv_mapon = parts[3];
        inv_start = parts[4];
        inv_end = parts[5];
        
        if (inv_iso == iso && inv_chr == chr && inv_mapon == mapon) {
            if ((start <= inv_end + 1) && (end >= inv_start - 1)) {
                delete inversions[inv_key];
                
                new_start = (start < inv_start) ? start : inv_start;
                new_end = (end > inv_end) ? end : inv_end;
                new_key = iso ":" chr ":" mapon ":" new_start ":" new_end;
                inversions[new_key] = 1;
                
                merged = 1;
                break;
            }
        }
    }

    # Si pas de fusion, ajouter comme nouvelle inversion
    if (!merged) {
        new_key = iso ":" chr ":" mapon ":" start ":" end;
        inversions[new_key] = 1;
    }
}
' "$INPUT_FILE" >> "$TMP_FILE"

# Trier le fichier résultat
sort -k1,1 -k2,2 -k5,5 -k3,3n "$TMP_FILE" >> "$OUTPUT_FILE"

# Supprimer le fichier temporaire
rm "$TMP_FILE"

# Filtrer les inversions de moins de 100pb
temp_file=$(mktemp)
awk '
BEGIN {FS = OFS = "\t";}
NR == 1 || ($4 - $3) >= 100 {print;}
' "$OUTPUT_FILE" > "$temp_file"
mv "$temp_file" "$OUTPUT_FILE"

echo "Traitement terminé. Résultats dans $OUTPUT_FILE"