import MDAnalysis as mda
from MDAnalysis.analysis import distances
import glob

# Load the MD trajectory
traject = sorted(glob.glob("/Users/floriangrun/Downloads/all-2/trajectories/frame_*.pdb"))

u = mda.Universe('/Users/floriangrun/Downloads/all-2/topology.pdb', traject)

# Select sulfur atoms in Cysteine residues
cysteines = u.select_atoms('resname CYS and name SG')

def meets_criteria(dist_matrix, threshold=2.2):
    # Check if all pairs of Cysteine sulfur atoms are within a certain distance threshold.
    n_cysteines = len(cysteines)
    for i in range(n_cysteines):
        for j in range(i + 1, n_cysteines):
            if dist_matrix[i, j] > threshold:
                return False
    return True



# Iterate over each frame of the trajectory
for ts in u.trajectory:
    # Calculate the distance matrix for the sulfur atoms
    dist_matrix = distances.distance_array(cysteines.positions, cysteines.positions)
    
    # Check if the frame meets the criteria for potential disulfide bonds
    if meets_criteria(dist_matrix):
        print(f"Frame {ts.frame} meets the criteria for disulfide bond formation")
        # Export this frame as a PDB file
        filename = f"frame_{ts.frame}.pdb"
        with mda.Writer(filename, u.atoms.n_atoms) as W:
            W.write(u.atoms)
        print(f"Exported frame {ts.frame} to {filename}")
