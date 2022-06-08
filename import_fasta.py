#file to import sequence data from fasta files and output it to a text file with each line seperated by a newline
from Bio import SeqIO

def import_fasta(fasta_file):
    """
    Import fasta file and return a list of sequences
    """
    sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(record.seq)
        print('.', end='')
    return sequences

def load_from_file(file_name, output_file):
    """
    Load a list of sequences from a file
    """
    with open(file_name) as f:
        import_fasta(f)
    
    with open(output_file, "w") as f:
        for sequence in import_fasta(file_name):
            f.write(str(sequence) + "\n")

#Change these lines to match your file names
FilePath = "./Testing/mouse.1.rna.fna"
OutputFile = "./Testing/MouseTranscriptome.txt"
load_from_file(FilePath, OutputFile)