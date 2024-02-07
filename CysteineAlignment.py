import numpy as np
from Bio import SeqIO, pairwise2
import matplotlib.pyplot as plt
import mplcursors
import pickle

def find_cysteines(sequence):
    return [i for i, amino_acid in enumerate(sequence) if amino_acid == 'C']

def compute_pairwise_distances(cysteine_positions):
    cysteine_positions = np.array(cysteine_positions)
    distance_matrix = np.abs(cysteine_positions - cysteine_positions[:, np.newaxis])
    return distance_matrix

def compare_matrices(mat1, mat2):
    if mat1.shape != mat2.shape:
        if mat1.size < mat2.size:
            mat1 = np.pad(mat1, ((0, mat2.shape[0] - mat1.shape[0]), (0, mat2.shape[1] - mat1.shape[1])), 'constant')
        else:
            mat2 = np.pad(mat2, ((0, mat1.shape[0] - mat2.shape[0]), (0, mat1.shape[1] - mat2.shape[1])), 'constant')
    return np.sum(np.abs(mat1 - mat2))

def calculate_sequence_identity(seq1, seq2):
    alignments = pairwise2.align.globalxx(seq1, seq2)
    alignment = alignments[0]
    aligned_seq1, aligned_seq2, _, _, _ = alignment

    matches = sum(res1 == res2 for res1, res2 in zip(aligned_seq1, aligned_seq2))
    total_length = max(len(aligned_seq1), len(aligned_seq2))
    return (matches / total_length) * 100

def process_fasta(file_path):
    pairwise_distances = {}
    sequences = SeqIO.parse(file_path, "fasta")
    for record in sequences:
        cysteine_positions = find_cysteines(str(record.seq))
        if cysteine_positions:  # Only include sequences with cysteines
            pairwise_distances[record.id] = compute_pairwise_distances(cysteine_positions)
    return pairwise_distances

def plot_cysteine_vs_identity(results):
    cysteine_scores = []
    identities = []
    labels = []  # To store the pair names for annotations

    zoom_cysteine_scores = []
    zoom_identities = []
    zoom_labels = []

    for pair, values in results.items():
        cysteine_score, identity = values
        cysteine_scores.append(cysteine_score)
        identities.append(identity)
        labels.append(pair)  # Append the pair name

        # Check if the pair falls within the specified range for zoomed plot
        if identity < 30 and cysteine_score < 100:
            zoom_cysteine_scores.append(cysteine_score)
            zoom_identities.append(identity)
            zoom_labels.append(pair)

    # Creating the full scatter plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

    scatter1 = ax1.scatter(cysteine_scores, identities)
    ax1.set_xlabel('Cysteine Score')
    ax1.set_ylabel('Sequence Identity (%)')
    ax1.set_title('Cysteine Score vs Sequence Identity')
    ax1.grid(True)

    # Creating the zoomed-in scatter plot
    scatter2 = ax2.scatter(zoom_cysteine_scores, zoom_identities)
    ax2.set_xlabel('Cysteine Score (Zoomed)')
    ax2.set_ylabel('Sequence Identity (%) (Zoomed)')
    ax2.set_title('Zoomed In: Cysteine Score < 100 & Identity < 30%')
    ax2.set_xlim(0, 100)
    ax2.set_ylim(0, 30)
    ax2.grid(True)

    # Adding interactivity with mplcursors
    tooltip1 = mplcursors.cursor(scatter1, hover=True)
    tooltip1.connect("add", lambda sel: sel.annotation.set_text(labels[sel.target.index]))

    tooltip2 = mplcursors.cursor(scatter2, hover=True)
    tooltip2.connect("add", lambda sel: sel.annotation.set_text(zoom_labels[sel.target.index]))

    plt.tight_layout()
    plt.show()

def score_and_compare_sequences(distances_dict):
    results = {}
    for id1, mat1 in distances_dict.items():
        for id2, mat2 in distances_dict.items():
            if id1 != id2:
                cysteine_score = compare_matrices(mat1, mat2)
                identity = calculate_sequence_identity(str(distances_dict[id1]), str(distances_dict[id2]))
                results[id1 + "_" + id2] = [cysteine_score, identity]
    
    # Sort and find the first instance with the lowest identity and cysteine score of 0
    best_match = min((pair for pair in results.items() if pair[1][0] == 0), key=lambda x: x[1][1], default=None)
    
    if best_match:
        print(f"Best Match: {best_match[0]}, Cysteine Score: {best_match[1][0]}, Sequence Identity: {best_match[1][1]:.2f}%")
    else:
        print("No matches found with a cysteine score of 0.")

    return results

# Example usage
fasta_file = "/Users/floriangrun/Downloads/uniparc_taxonomy_id_5455_AND_protein_2024_02_02.fasta"  # Replace with your FASTA file path
export_data = "/Users/floriangrun/Desktop/data_cysteine_dmaps.pkl"
distances = process_fasta(fasta_file)
#results = score_and_compare_sequences(distances)

#plot_cysteine_vs_identity(results)
with open(export_data, "wb") as path:
    pickle.dump(distances, path)

#for i,j in results.items():
#    print(f"Pair: {i}, Cysteine Score: {j[0]}, Sequence Identity: {j[1]:.2f}%")

