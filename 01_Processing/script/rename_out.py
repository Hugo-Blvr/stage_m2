import pandas as pd
from pathlib import Path
import re
import sys
import warnings
warnings.simplefilter(action='ignore', category=pd.errors.ParserWarning)

# ============= ⚠ HARDCODED VALUES TO ADAPT IF INPUT DATA FORMAT CHANGES =============
DELETE_CHR = ['chr17_1', 'chr18', 'chr19']
ACCESSORY_CHR = ['chr3', 'chr16', 'chr17']
# =====================================================================================

def rename_chr(chr_name, filname, acces_chr = ACCESSORY_CHR):
    """Rename chromosome according to predefined rules."""
    match = re.match(r'^chr(\d+)$', chr_name)
    if match:
        num = int(match.group(1))
        return f'chr{num + 2000}' if chr_name in acces_chr else f'chr{num + 1000}'
    else:
        print(f"\n\t ⚠ Warning: chromosome {chr_name} has an unexpected format in {filname} \n")
        return chr_name

def process_files(input_dir, output_dir, delete_chr = DELETE_CHR):
    """Process all .out files in input_dir, apply chromosome renaming, and save to output_dir."""
    input_path, output_path = Path(input_dir), Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    for file_path in input_path.glob('*.out'):
        df = pd.read_csv(file_path, sep=r'\s+', engine='python', header=None, 
                        index_col=False, skip_blank_lines=True, skiprows=2)
        df.columns = ['score', 'div.', 'del.', 'ins.', 'sequence', 'begin', 'end', '(left)', 'strand', 
                    'repeat', 'class/family', 'begin', 'end', '(left)', 'ID']
        if df.empty or df.shape[1] < 5:
            print(f"Skipping file {file_path.name} due to unexpected format.")
            continue

        df = df[~df['sequence'].isin(delete_chr)] # Filter out unwanted chromosomes
        df['sequence'] = df['sequence'].apply(lambda chr_name: rename_chr(chr_name, filname=file_path)) # Rename chromosomes
        
        output_file = output_path / f'{file_path.name.split('.')[0]}_rename_chr.out'
        df.to_csv(output_file, sep='\t', index=False)
        print(f"Processed: {file_path.name} → {output_file.name}")

if __name__ == "__main__":
    if len(sys.argv) != 3: print(f"Usage: python3 {sys.argv[0]} <input_dir> <output_dir>")
    else: process_files(sys.argv[1], sys.argv[2])