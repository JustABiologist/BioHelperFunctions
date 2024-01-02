import Bio.PDB
import argparse

def get_cysteines(pdb_file):
    """
    Extracts cysteine residues from the PDB file.
    """
    parser = Bio.PDB.PDBParser()
    structure = parser.get_structure('PDB', pdb_file)
    cysteines = []

    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname() == 'CYS':
                    cysteines.append((chain.id, residue.id[1], residue))

    return cysteines

def choose_disulfide_bonds(cysteines):
    """
    Allows user to select pairs of cysteines to form disulfide bonds.
    """
    available_cysteines = cysteines.copy()
    selected_bonds = []

    while available_cysteines:
        print("\nAvailable Cysteine Residues:")
        for i, (chain, res_id, _) in enumerate(available_cysteines):
            print(f"{i+1}. Chain {chain}, Residue {res_id}")

        choice = input("\nSelect a cysteine for disulfide bond (or type 'skip' to leave it unbonded, 'exit' to finish): ")
        if choice.lower() == 'exit':
            break
        if choice.lower() == 'skip':
            continue

        try:
            choice = int(choice) - 1
            if choice < 0 or choice >= len(available_cysteines):
                raise ValueError
        except ValueError:
            print("Invalid selection. Please try again.")
            continue

        chosen_cys = available_cysteines.pop(choice)
        print(f"Selected: Chain {chosen_cys[0]}, Residue {chosen_cys[1]}")

        print("\nChoose a partner for this cysteine:")
        for i, (chain, res_id, _) in enumerate(available_cysteines):
            print(f"{i+1}. Chain {chain}, Residue {res_id}")

        partner_choice = input("\nSelect a partner cysteine (or type 'skip' to leave it unbonded): ")
        if partner_choice.lower() == 'skip':
            continue

        try:
            partner_choice = int(partner_choice) - 1
            if partner_choice < 0 or partner_choice >= len(available_cysteines):
                raise ValueError
        except ValueError:
            print("Invalid selection. Please try again.")
            available_cysteines.append(chosen_cys)
            continue

        partner_cys = available_cysteines.pop(partner_choice)
        selected_bonds.append((chosen_cys, partner_cys))
        print(f"Bond formed between Chain {chosen_cys[0]}, Residue {chosen_cys[1]} and Chain {partner_cys[0]}, Residue {partner_cys[1]}")

    return selected_bonds

def print_disulfide_distances(selected_bonds):
    """
    Prints the distances between selected disulfide bond pairs.
    """
    for (chain1, res_id1, residue1), (chain2, res_id2, residue2) in selected_bonds:
        sg1 = residue1['SG'].get_vector()
        sg2 = residue2['SG'].get_vector()
        distance = (sg1 - sg2).norm()

        print(f"Distance between CYS {res_id1} (Chain {chain1}) and CYS {res_id2} (Chain {chain2}): {distance:.2f} Å")


def append_ssbond_records(pdb_file, selected_bonds):
    """
    Appends SSBOND records to the PDB file.
    """
    with open(pdb_file, 'a') as file:
        for i, ((chain1, res_id1, _), (chain2, res_id2, _)) in enumerate(selected_bonds, start=1):
            file.write(f"SSBOND  {i:>2} CYS {chain1} {res_id1:>4}    CYS {chain2} {res_id2:>4}\n")


def main(pdb_file):
    cysteines = get_cysteines(pdb_file)
    selected_bonds = choose_disulfide_bonds(cysteines)
    print_disulfide_distances(selected_bonds)
    append_ssbond_records(pdb_file, selected_bonds)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process a PDB file to select and append disulfide bonds.")
    parser.add_argument("pdb_file", help="Path to the PDB file")
    args = parser.parse_args()

    main(args.pdb_file)