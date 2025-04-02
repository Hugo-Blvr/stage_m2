#!/bin/bash

# Vérification des arguments
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 dossier_fasta"
    exit 1
fi

DIR="$1"

# Vérification que le dossier existe
if [ ! -d "$DIR" ]; then
    echo "Erreur : le dossier $DIR n'existe pas."
    exit 1
fi

# Traitement des fichiers FASTA
for file in "$DIR"/*.fasta "$DIR"/*.fa; do
    [ -e "$file" ] || continue  # Vérifie si des fichiers existent
    awk '{if($0 ~ /^>/) print $0; else print toupper($0)}' "$file" > "$file.tmp" && mv "$file.tmp" "$file"
done

echo "Mise en majuscules terminée pour tous les fichiers dans $DIR."
