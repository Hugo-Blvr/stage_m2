#!/bin/bash

# ============ IMPORTANT ==============
# ECART MIN ENTRE LES INVERSION : 10pb 
# =====================================

if [ -z "$2" ]; then
  echo "Usage: $0 <input_file> <output_file>"
  exit 1
fi

input_file="$1"
output_file="$2"
temp_file=$(mktemp)



# Écrire directement l'en-tête dans le fichier de sortie final
echo -e "iso\tchr\tstart\tend" > "$temp_file"
echo -e "iso\tchr\tstart\tend" > "$output_file"
# Traitement avec awk pour extraire les inversions
awk '
BEGIN {
    block = "";
}

/^a/ {
    if (block != "") {
        has_inversion = 0;
        for (chr_key in chr_plus) {
            if (chr_key in chr_minus) {
                has_inversion = 1;
                n = split(block, lines, "\n");
                for (i = 1; i <= n; i++) {
                    line = lines[i];
                    if (substr(line, 1, 1) == "s") {
                        split(line, fields);
                        seq = fields[2];
                        start = fields[3];
                        size = fields[4];
                        
                        split(seq, parts, ".");
                        iso = parts[1];
                        chr_name = parts[2];
                        for (j = 3; j <= length(parts); j++) {
                            chr_name = chr_name "." parts[j];
                        }
                        
                        if (chr_name == chr_key) {
                            end = start + size -1;
                            print iso "\t" chr_name "\t" start "\t" end;
                        }
                    }
                }
            }
        }
    }

    block = $0;
    delete chr_plus;
    delete chr_minus;
    next;
}

/^s/ {
    block = block "\n" $0;

    split($2, parts, ".");
    chr_name = parts[2];
    for (j = 3; j <= length(parts); j++) {
        chr_name = chr_name "." parts[j];
    }
    
    orientation = $5;

    if (orientation == "+") {
        chr_plus[chr_name] = 1;
    } else if (orientation == "-") {
        chr_minus[chr_name] = 1;
    }
    
    next;
}

{
    block = block "\n" $0;
}

END {
    if (block != "") {
        has_inversion = 0;
        for (chr_key in chr_plus) {
            if (chr_key in chr_minus) {
                has_inversion = 1;
                n = split(block, lines, "\n");
                for (i = 1; i <= n; i++) {
                    line = lines[i];
                    if (substr(line, 1, 1) == "s") {
                        split(line, fields);
                        seq = fields[2];
                        start = fields[3];
                        size = fields[4];
                        
                        split(seq, parts, ".");
                        iso = parts[1];
                        chr_name = parts[2];
                        for (j = 3; j <= length(parts); j++) {
                            chr_name = chr_name "." parts[j];
                        }
                        
                        if (chr_name == chr_key) {
                            end = start + size - 1;
                            print iso "\t" chr_name "\t" start "\t" end;
                        }
                    }
                }
            }
        }
    }
}
' "$input_file" | sort -k1,1 -k2,2 -k3,3n > "$temp_file"

# Fusion des alignements consécutifs pour la même paire iso-chr
awk '
BEGIN {
    FS = OFS = "\t";
    last_iso = "";
    last_chr = "";
    last_start = "";
    last_end = "";
}
{
    if (NR == 1) {
        print;
    } else {
        if ($1 == last_iso && $2 == last_chr && $3 <= last_end + 10) {
            if ($4 > last_end) {
                last_end = $4;
            }
        } else {
            if (last_iso != "") {
                print last_iso, last_chr, last_start, last_end;
            }
            last_iso = $1;
            last_chr = $2;
            last_start = $3;
            last_end = $4;
        }
    }
}
END {
    if (last_iso != "") {
        print last_iso, last_chr, last_start, last_end;
    }
}
' "$temp_file" >> "$output_file"

# Supprime le fichier temporaire
rm "$temp_file"

# filtre les inversion de - de 100pb
temp_file=$(mktemp)
awk '
BEGIN {FS = OFS = "\t";}
NR == 1 || ($4 - $3) >= 100 {print;}
' "$output_file" > "$temp_file"
# Remplacez le fichier de sortie par la version filtrée
mv "$temp_file" "$output_file"