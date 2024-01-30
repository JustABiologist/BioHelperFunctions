import numpy as np
import matplotlib.pyplot as plt

# Define your amino acid positions of interest
selected_positions = [46, 47, 48 ,49, 50, 51, 52]  # Example positions, replace with your actual positions

# Create a color map for the amino acids
colors = {selected_positions[i]: plt.cm.tab10(i) for i in range(len(selected_positions))}

# Initialize a dictionary to hold our data
data = {pos: [] for pos in selected_positions}

# Read the file and parse the data
with open('../Strukturen_116/EC116Dimer_cdcab/rama_all_003.xvg', 'r') as file:
    for line in file:
        if not line.startswith('@') and not line.startswith('#'):
            parts = line.split()
            position = int(parts[2].split('-')[-1])  # Extract the residue position
            if position in selected_positions:
                phi, psi = float(parts[0]), float(parts[1])
                data[position].append((phi, psi))

# Determine the number of frames (assuming all selected_positions have the same number of frames)
num_frames = len(data[selected_positions[0]])

# Prepare the plot
plt.figure(figsize=(8, 6))

# Loop over each selected position and plot
for pos, angles in data.items():
    phi, psi = zip(*angles)  # This unzips into two lists: phi angles and psi angles
    # Decrease opacity linearly from 1 to 0.1 over the frames
    alphas = np.linspace(.1, 1, num_frames)
    for i in range(num_frames):
        plt.scatter(phi[i], psi[i], color=colors[pos], alpha=alphas[i])

# Adding a custom legend
legend_elements = [plt.Line2D([0], [0], marker='o', color='w', label=f'Residue {pos}',
                              markerfacecolor=colors[pos], markersize=5) for pos in selected_positions]
plt.legend(handles=legend_elements)

plt.title('Ramachandran Plot with Time-Dependent Opacity')
plt.xlabel('Phi angles (degrees)')
plt.ylabel('Psi angles (degrees)')
plt.xlim(-180, 180)
plt.ylim(-180, 180)
plt.grid(True)

# Show plot
plt.show()
