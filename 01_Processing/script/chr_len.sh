#!/bin/bash

if [ -z "$2" ]; then
  echo "Usage: $0 <fasta_folder> <output_file>"
  exit 1
fi

folder="$1"
output="$2"

# Vide le fichier de sortie s'il existe déjà
> "$output"

# Parcourt tous les fichiers .fasta dans le dossier
for fasta in "$folder"/*.fasta; do
  filename=$(basename "$fasta" .fasta)
  IFS='_' read -r part1 part2 _ <<< "$filename"
  iso="${part1}_${part2}"

  awk -v iso="$iso" '
    /^>/ {
      if (seqname != "") {
        printf "%s\t%s\t%d\n", iso, seqname, length(seq) >> "'$output'"
      }
      seqname = substr($0, 2)
      seq = ""
      next
    }
    {
      seq = seq $0
    }
    END {
      if (seqname != "") {
        printf "%s\t%s\t%d\n", iso, seqname, length(seq) >> "'$output'"
      }
    }
  ' "$fasta"
done

# Ajoute l'en-tête avec des tabulations
sed -i '1i iso\tchr\tsize' "$output"
