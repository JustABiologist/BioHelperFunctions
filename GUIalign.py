import tkinter as tk
from tkinter import filedialog
from tkinter import messagebox
import threading
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Align

# Sequence transformation functions
def reverse_seq(seq):
    return seq[::-1]

def reverse_complement_seq(seq):
    return str(Seq(seq).reverse_complement())

def complement_seq(seq):
    return str(Seq(seq).complement())

# Alignment and translation function
def align_and_translate(seq, pac_seq):
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'  # Change to local alignment
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.match_score = 2
    aligner.mismatch_score = -1

    alignments = aligner.align(seq, pac_seq)
    best_alignment = next(alignments)

    target_aligned_regions = best_alignment.aligned[0]
    if len(target_aligned_regions) > 0:
        start_index_pac = target_aligned_regions[0][0]
        end_index_pac = target_aligned_regions[-1][1]

        konsensus_seq = pac_seq[:start_index_pac] + seq + pac_seq[end_index_pac:]

        # Find the start codon (ATG) and translate
        start_codon_index = konsensus_seq.find("ATG")
        if start_codon_index != -1:
            translated_seq = Seq(konsensus_seq[start_codon_index:]).translate(to_stop=True)
            return translated_seq, best_alignment.score

    return "", 0

# Process sequences function
def process_sequences(pac_file, fasta_file, output_dir):
    with open(pac_file, 'r') as file:
        pac_seq = file.read().strip()

    for record in SeqIO.parse(fasta_file, "fasta"):
        best_score = None
        best_consensus = None
        for transform_func in [lambda x: x, reverse_seq, reverse_complement_seq, complement_seq]:
            transformed_seq = transform_func(str(record.seq))
            consensus_seq, score = align_and_translate(transformed_seq, pac_seq)

            if best_score is None or score > best_score:
                best_score = score
                best_consensus = consensus_seq

        if best_consensus:
            output_filename = output_dir + "/" + record.id + "_consensus.fasta"
            with open(output_filename, "w") as output_file:
                output_file.write(f">{record.id}\n{best_consensus}\n")
            gui_log(f"Processed: {record.id}")
        else:
            gui_log(f"No suitable alignment found for: {record.id}")

    gui_log("Processing completed.")

def gui_log(message):
    output_text.insert(tk.END, message + "\n")
    output_text.see(tk.END)

def browse_file(entry):
    filename = filedialog.askopenfilename()
    entry.delete(0, tk.END)
    entry.insert(0, filename)

def browse_directory(entry):
    directory = filedialog.askdirectory()
    entry.delete(0, tk.END)
    entry.insert(0, directory)

def start_thread():
    pac_file = pac_entry.get()
    fasta_file = fasta_entry.get()
    output_dir = output_entry.get()
    if not pac_file or not fasta_file or not output_dir:
        messagebox.showerror("Error", "Please select all files and output directory")
        return
    threading.Thread(target=process_sequences, args=(pac_file, fasta_file, output_dir)).start()

root = tk.Tk()
root.title("Sequence Alignment and Consensus Generation Tool")

tk.Label(root, text="PAC Sequence File:").grid(row=0, column=0, sticky="w")
pac_entry = tk.Entry(root, width=50)
pac_entry.grid(row=0, column=1)
tk.Button(root, text="Browse", command=lambda: browse_file(pac_entry)).grid(row=0, column=2)

tk.Label(root, text="FASTA File:").grid(row=1, column=0, sticky="w")
fasta_entry = tk.Entry(root, width=50)
fasta_entry.grid(row=1, column=1)
tk.Button(root, text="Browse", command=lambda: browse_file(fasta_entry)).grid(row=1, column=2)

tk.Label(root, text="Output Directory:").grid(row=2, column=0, sticky="w")
output_entry = tk.Entry(root, width=50)
output_entry.grid(row=2, column=1)
tk.Button(root, text="Browse", command=lambda: browse_directory(output_entry)).grid(row=2, column=2)

tk.Button(root, text="Start Processing", command=start_thread).grid(row=3, column=1, sticky="ew")

output_text = tk.Text(root, height=10, width=60)
output_text.grid(row=4, column=0, columnspan=3)
scroll = tk.Scrollbar(root, command=output_text.yview)
scroll.grid(row=4, column=3, sticky='nsew')
output_text['yscrollcommand'] = scroll.set

root.mainloop()