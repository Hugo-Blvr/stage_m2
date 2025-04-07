# InversionCaller

## Description

Outil conçu pour détecter des inversions génomiques entre plusieurs génomes.  
Le script fonctionne en extrayant les séquences des chromosomes communs entre chaque paire de génomes, en les alignant avec minimap2, puis en utilisant syri pour détecter les inversions. Les résultats sont ensuite filtrés en fonction d'un seuil d'identité de séquence minimale et compilés dans un fichier de sortie au format TSV.  

Ce processus permet de comparer exhaustivement tous les génomes présents dans le répertoire spécifié, facilitant ainsi l'identification des inversions chromosomiques potentielles entre eux.  

## Prérequis

Le script nécessite l'installation préalable des outils suivants :  
- **minimap2** : pour l'alignement des séquences génomiques  
- **samtools** : pour la manipulation des fichiers SAM/BAM  
- **syri** : pour la détection des variants structuraux  

## Installation

1. Assurez-vous que tous les prérequis sont installés et disponibles dans votre `PATH`.  
2. Téléchargez le script et rendez-le exécutable :  
   ```bash
   chmod +x inversioncaller.sh
   ```

## Utilisation

### Syntaxe de base

```bash
./inversioncaller.sh -d <dossier_entree> -o <fichier_sortie.tsv> [options]
```

### Options

| Option        | Description                                           | Valeur par défaut |
|--------------|------------------------------------------------------|------------------|
| `-d`, `--directory` | **[Obligatoire]** Dossier contenant les fichiers FASTA à analyser | - |
| `-o`, `--output`    | Chemin du fichier TSV de sortie               | `./inv_calling.tsv` |
| `-i`, `--identity`  | Identité de séquence minimale (valeur entre 0 et 1) | `0` |
| `-t`, `--threads`   | Nombre de threads à utiliser                   | `8` |
| `-h`, `--help`      | Affiche l'aide                                 | - |

### Exemple

```bash
./inversioncaller.sh -d ./genomes -o ./resultats/inversions.tsv -i 0.8 -t 16
```

## Format de sortie

Le fichier de sortie est au format TSV avec les colonnes suivantes :

1. Nom du génome cible  
2. Chromosome cible  
3. Position de début dans le génome cible  
4. Position de fin dans le génome cible  
5. Nom du génome requête  
6. Chromosome requête  
7. Position de début dans le génome requête  
8. Position de fin dans le génome requête  

## Fonctionnement détaillé

### Vue d'ensemble du pipeline

1. **Préparation** : Vérification des prérequis et des paramètres d'entrée  
2. **Comparaison par paires** : Pour chaque paire de génomes FASTA dans le dossier d'entrée  
3. **Analyse par chromosome** : Pour chaque chromosome commun entre deux génomes  
   - Extraction des séquences chromosomiques  
   - Alignement avec minimap2  
   - Détection des variants structuraux avec SyRI  
   - Extraction et formatage des inversions détectées  

### Choix techniques

#### Utilisation de minimap2
Minimap2 a été choisi pour sa capacité à aligner efficacement de longs fragments d'ADN, ce qui est crucial pour la détection précise d'inversions chromosomiques. Le mode `asm5` est optimisé pour comparer des assemblages génomiques avec ~5% de divergence. L'utilisation de ce script avec des génomes trop divergent le rendra donc moins prècis.
→ https://github.com/lh3/minimap2

#### Filtrage par identité de séquence
Le script permet de filtrer les alignements selon un seuil d'identité de séquence, ce qui aide à éliminer les alignements de faible qualité tout en conservant les signaux biologiquement pertinents.

L'identité est définie comme : \( \text{identity} = \frac{\text{matches}}{\text{sequence length}} \)


#### SyRI pour la détection des variants
SyRI (**Synteny and Rearrangement Identifier**) est spécialisé dans la détection de réarrangements génomiques, y compris les inversions, à partir d'alignements génomiques.
→ https://github.com/schneebergerlab/syri


## Journal d'exécution
Un fichier journal (`inversions_calling.log`) est généré dans le dossier de sortie, contenant des informations détaillées sur le déroulement de l'analyse et d'éventuelles erreurs.

## Limitations
Le script traite uniquement les chromosomes portant le même identifiant entre génomes, présupposant qu'ils représentent des chromosomes homologues. De plus pour un maximum de prècision les génomes ne doivent pas dépasser les 5% de divergence. 

+ filter with filter_inv.py ... 
