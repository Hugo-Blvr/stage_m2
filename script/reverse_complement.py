from Bio import SeqIO
from Bio.Seq import Seq
import sys
    
def reverse_complement_contig(fasta_file, contig_name):
    records = list(SeqIO.parse(fasta_file, "fasta"))
    
    found = False
    for record in records:
        if record.id == contig_name:
            record.seq = record.seq.reverse_complement()
            found = True
            break
    
    if not found:
        print(f"Erreur: Contig '{contig_name}' non trouvé dans le fichier.")
        return
    
    with open(fasta_file, 'w') as output_handle: SeqIO.write(records, output_handle, "fasta")    
    print(f"Le contig {contig_name} a été remplacé par son complément inverse.")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 reverse_complement.py <fichier_fasta> <nom_contig>")
        sys.exit(1)
    
    reverse_complement_contig(sys.argv[1], sys.argv[2])