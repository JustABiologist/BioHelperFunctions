import os
import glob
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import molecular_weight

def calculate_mw_with_biopython(sequence):
    return molecular_weight(sequence, seq_type='protein') / 1000  # Convert to kDa

def process_fasta_files(directory, sequence, output):
    results = []
    for fasta_file in glob.glob(os.path.join(directory, '*.fasta')):
        for record in SeqIO.parse(fasta_file, "fasta"):
            # Append the sequence to the C-terminal end if provided
            new_seq = str(record.seq) + sequence if sequence else str(record.seq)
            # Calculate the molecular weight
            mw = calculate_mw_with_biopython(new_seq)
            # Add to the results list
            results.append({'Name': record.id, 'MW': mw})
    df = pd.DataFrame(results)
    # Check the file extension to decide whether to save as CSV or Excel
    if output.endswith('.xlsx'):
        with pd.ExcelWriter(output) as writer:
            df.to_excel(writer, index=False)
    else:
        df.to_csv(output, index=False)
    
    print(f"Results saved to {output}")

def main():
    # Command-line argument parser
    parser = argparse.ArgumentParser(description="Append sequence and calculate molecular weight of FASTA files.")
    parser.add_argument('directory', help="Directory containing FASTA files")
    parser.add_argument('sequence', nargs='?', default=None, help="Optional: Sequence to append to the C-terminal end of the sequences")
    parser.add_argument('--output', help="Output CSV file name", default="output.csv")
    args = parser.parse_args()

    # Process the FASTA files
    process_fasta_files(args.directory, args.sequence, args.output)

if __name__ == "__main__":
    main()
