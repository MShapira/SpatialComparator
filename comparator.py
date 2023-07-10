from math import sqrt, exp
from Bio.PDB import Chain, PDBParser, Residue, PDBIO, Superimposer
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure
import sys
import os
import string

class ResidueDifference:
    def __init__(self, string: str):
        components = string.split('|')
        self.first_residue_id = components[0]
        self.first_residue_index = int(components[1])
        self.second_residue_id = components[2]
        self.second_residue_index = int(components[3])
        self.distance = components[4]
        self.similarity = components[5]
        self.powers_score = components[6]
        self.pam_similarity = components[7]


def load_grantham_table(path: str) -> dict:
    table = dict()

    with open(path) as file:
        header = file.readline()
        header_components = header[:-1].split(',')

        for line in file.readlines():
            line_components = line[:-1].split(',')
            first_aminoacid = line_components[0]

            similarities = dict()

            for second_aminoacid_index in range(1, len(line_components)):
                value = float(line_components[second_aminoacid_index])
                second_aminoacid = header_components[second_aminoacid_index]
                similarities[second_aminoacid] = value

            table[first_aminoacid] = similarities

    return table


def load_pam_table(path: str) -> dict:
    table = dict()

    with open(path) as file:
        header = file.readline()
        header_components = header.rstrip('\n').split(',')

        for line in file.readlines():
            line_components = line.rstrip('\n').split(',')
            first_aminoacid = line_components[0]

            similarities = dict()

            for second_aminoacid_index in range(1, len(line_components)):
                value = int(line_components[second_aminoacid_index])
                second_aminoacid = header_components[second_aminoacid_index]
                similarities[second_aminoacid] = value

            table[first_aminoacid] = similarities

    return table


def load_aminoacid_letters_table() -> dict:
    return {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H',
            'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
            'TYR': 'Y', 'VAL': 'V', 'MSE': 'M', 'CYM': 'C', 'HIE': 'H', 'HIP': 'H'}


grantham_table = load_grantham_table('Database/grantham-74.csv')
pam_table = load_pam_table('Database/pam-250.csv')
aminoacid_letters_table = load_aminoacid_letters_table()


def is_aminoacid(name: str) -> bool:
    return name in grantham_table.keys()


def get_c_alpha_position(residue: Residue) -> list:
    for atom in residue.get_atoms():
        if atom.id == 'CA':
            return atom.coord
    return None


def get_distance(first_position: list, second_position: list) -> float:
    return sqrt((first_position[0] - second_position[0])**2 + (first_position[1] - second_position[1])**2 + (first_position[2] - second_position[2])**2)


def get_aminoacid_similarity(first_aminoacid: str, second_aminoacid: str) -> float:
    return grantham_table[first_aminoacid][second_aminoacid]


def get_residue_similarity(first_residue: Residue, second_residue: Residue) -> float:
    if is_aminoacid(first_residue.get_resname()) and is_aminoacid(second_residue.get_resname()):
        return get_aminoacid_similarity(first_residue.get_resname(), second_residue.get_resname())
    else:
        return 0.0


def get_pam_similarity(first_residue: Residue, second_residue: Residue) -> int:
    if is_aminoacid(first_residue.get_resname()) and is_aminoacid(second_residue.get_resname()):
        first_letter = aminoacid_letters_table[first_residue.get_resname()]
        second_letter = aminoacid_letters_table[second_residue.get_resname()]
        return pam_table[first_letter][second_letter]
    else:
        return -100


def calculate_powers_score(first_residue: Residue, second_residue: Residue, distance: float) -> float:
    similarity = get_residue_similarity(first_residue, second_residue)

    if distance <= 1.0:
        return similarity
    else:
        return similarity * exp(1-distance)**2


def format_difference_string(first_residue: Residue, second_residue: Residue, distance: float, powers_score: float) -> str:
    return f'{first_residue.get_resname()}|{first_residue.id[1]}|{second_residue.get_resname()}|{second_residue.id[1]}|{distance:.6f}|{get_residue_similarity(first_residue, second_residue):.6f}|{powers_score}|{get_pam_similarity(first_residue, second_residue)};'


def get_protein_difference(first_protein: Chain, second_protein: Chain) -> str:
    difference = ''

    total_powers_score = 0.0

    for first_residue in first_protein.get_residues():
        first_c_alpha_position = get_c_alpha_position(first_residue)
        if first_c_alpha_position is None:
            continue

        # look for nearest residue from second protein
        nearest_residue = None
        nearest_residue_distance = 100000.0
        for second_residue in second_protein.get_residues():
            second_c_alpha_position = get_c_alpha_position(second_residue)
            if second_c_alpha_position is None:
                continue

            distance = get_distance(first_c_alpha_position, second_c_alpha_position)
            if distance < nearest_residue_distance:
                nearest_residue = second_residue
                nearest_residue_distance = distance
        if nearest_residue is None:
            print(f'something went wrong: failed to find nearest residue for {first_residue.id} of {first_protein.id} in {second_protein.id}')
            continue

        # check similarity
        powers_score = calculate_powers_score(first_residue, nearest_residue, nearest_residue_distance)
        if True:  # threshold for identical aminoacids spaced apart 1.5 angstroms
            difference += format_difference_string(first_residue, nearest_residue, nearest_residue_distance, powers_score)
        total_powers_score += powers_score

    difference += f'TPS={total_powers_score}'

    return difference


def construct_superimposed_pdb() -> Structure:
    # collect PDBs for future chains
    pdb_filenames = dict()
    while True:
        # PDB filename (or empty to break the loop)
        pdb_filename = input('Add new chain PDB path (or just hit ENTER to finish collecting PDBs): ')
        if pdb_filename == '':
            break

        # chain name (valid non-empty string)
        chain_name = input('Set chain name (or just hit ENTER to use file name): ')
        if chain_name == '':
            chain_name = os.path.basename(pdb_filename).split('.pdb')[0]
        pdb_filenames[chain_name] = pdb_filename

    print(f'PDBs to combine: {len(pdb_filenames)}')
    for index, chain_name in enumerate(pdb_filenames.keys()):
        print(f'{index+1}. {chain_name}: {pdb_filenames[chain_name]}')

    if len(pdb_filenames) == 0:
        print('Warning: no chains provided, exiting')
        exit()

    chain_names = list()
    chains = list()

    # add PDBs to superimposed structure, one by one
    parser = PDBParser(QUIET=True)
    for chain_index, chain_name in enumerate(pdb_filenames):
        pdb_filename = pdb_filenames[chain_name]
        structure = parser.get_structure(chain_name, pdb_filename)
        # pick first chain from PDB structure
        chain = next(structure.get_chains())
        # use consecutive alphabet letters (uppercase) as chain IDs
        chain.id = string.ascii_uppercase[chain_index]
        chain_names.append(chain_name)
        chains.append(chain)

    # superimpose chains
    reference_chain = chains[0]
    reference_atoms = [atom for atom in reference_chain.get_atoms() if atom.name == 'CA']
    for chain in chains[1:]:
        superimposer = Superimposer()
        chain_atoms = [atom for atom in chain.get_atoms() if atom.name == 'CA']
        min_atom_count = min(len(reference_atoms), len(chain_atoms))
        superimposer.set_atoms(reference_atoms[:min_atom_count-1], chain_atoms[:min_atom_count-1])
        rotation, transform = superimposer.rotran
        rotation = rotation.astype('f')
        transform = transform.astype('f')
        for atom in chain.get_atoms():
            atom.transform(rotation, transform)

    # save superimposed PDB
    superimposed_name = input('Enter superimposed PDB name: ')
    if superimposed_name == '':
        superimposed_name = 'result'
    superimposed_structure = Structure(superimposed_name)
    for chain_name, chain in zip(chain_names, chains):
        model = Model(chain_name)
        model.add(chain)
        superimposed_structure.add(model)
    io = PDBIO()
    io.set_structure(superimposed_structure)
    io.save(f'{superimposed_name}.pdb')

    return superimposed_structure


if __name__ == '__main__':
    structure = construct_superimposed_pdb()

    pdb_name = structure.id
    output_path = pdb_name

    # make folder for output files
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    parser = PDBParser(QUIET=True)

    # get proteins
    protein_names = [model.id for model in structure.get_models()]
    proteins = [chain for chain in structure.get_chains()]

    selected_protein_name = protein_names[0]

    # fill protein difference table
    differences = dict()
    for first_protein_index in range(len(proteins)):
        first_protein = proteins[first_protein_index]
        first_protein_name = protein_names[first_protein_index]
        differences[first_protein_name] = dict()
        for second_protein_index in range(len(proteins)):
            second_protein = proteins[second_protein_index]
            second_protein_name = protein_names[second_protein_index]
            if first_protein_index == second_protein_index:
                differences[first_protein_name][second_protein_name] = ''
            else:
                difference = get_protein_difference(first_protein, second_protein)
                differences[first_protein_name][second_protein_name] = difference
        print(f'{first_protein_index+1} / {len(proteins)}')

    # get detailed difference tables for selected protein
    selected_protein_index = protein_names.index(selected_protein_name)
    selected_protein = proteins[selected_protein_index]
    rmsds = dict()
    similarities = dict()
    powers_scores = dict()
    residues = dict()
    residue_differences = dict()
    pam_similarities = dict()
    for protein_index in range(len(proteins)):
        if protein_index == selected_protein_index:
            continue
        protein = proteins[protein_index]
        protein_name = protein_names[protein_index]
        current_protein_rmsds = list()
        current_protein_similarities = list()
        current_protein_powers_scores = list()
        current_protein_residues = list()
        current_protein_pam_similarities = list()

        difference = differences[selected_protein_name][protein_name]

        difference_elements = [ResidueDifference(string) for string in difference.split(';')[:-1]]  # skip TPS component
        for residue in selected_protein.get_residues():
            difference_found = [difference for difference in difference_elements if difference.first_residue_index == residue.id[1]]
            if len(difference_found) != 0:
                current_protein_rmsds.append(difference_found[0].distance)
                current_protein_similarities.append(difference_found[0].similarity)
                current_protein_powers_scores.append(difference_found[0].powers_score)
                current_protein_residues.append(difference_found[0].second_residue_id)
                current_protein_pam_similarities.append(difference_found[0].pam_similarity)
            else:
                current_protein_rmsds.append('')
                current_protein_similarities.append('')
                current_protein_powers_scores.append('')
                current_protein_residues.append('')
                current_protein_pam_similarities.append(' ')
        rmsds[protein_name] = current_protein_rmsds
        similarities[protein_name] = current_protein_similarities
        powers_scores[protein_name] = current_protein_powers_scores
        residues[protein_name] = current_protein_residues
        pam_similarities[protein_name] = current_protein_pam_similarities
        residue_differences[protein_name] = difference_elements
    print('difference data reorganized')

    other_protein_names = [name for name in protein_names if name != selected_protein_name]

    with open(f'{output_path}/{pdb_name}.{selected_protein_name}.RMSD.csv', 'w') as file:
        file.write('x,y,z,residue,index,' + ','.join(other_protein_names) + '\n')
        residue_index = 0
        for residue in selected_protein.get_residues():
            c_alpha_position = get_c_alpha_position(residue)
            if c_alpha_position is None:
                residue_index += 1
                continue
            file.write(f'{c_alpha_position[0]:.6f},{c_alpha_position[1]:.6f},{c_alpha_position[2]:.6f},')
            file.write(f'{residue.get_resname()},{residue.id[1]},')
            file.write(','.join([rmsds[protein_name][residue_index] for protein_name in other_protein_names]))
            file.write('\n')
            residue_index += 1
    print('RMSDs saved')
    with open(f'{output_path}/{pdb_name}.{selected_protein_name}.Matrix.csv', 'w') as file:
        file.write('x,y,z,residue,index,' + ','.join(other_protein_names) + '\n')
        residue_index = 0
        for residue in selected_protein.get_residues():
            c_alpha_position = get_c_alpha_position(residue)
            if c_alpha_position is None:
                residue_index += 1
                continue
            file.write(f'{c_alpha_position[0]:.6f},{c_alpha_position[1]:.6f},{c_alpha_position[2]:.6f},')
            file.write(f'{residue.get_resname()},{residue.id[1]},')
            file.write(','.join([similarities[protein_name][residue_index] for protein_name in other_protein_names]))
            file.write('\n')
            residue_index += 1
    print('similarities saved')
    with open(f'{output_path}/{pdb_name}.{selected_protein_name}.BF.csv', 'w') as file:
        file.write('x,y,z,residue,index,' + ','.join(other_protein_names) + '\n')
        residue_index = 0
        for residue in selected_protein.get_residues():
            c_alpha_position = get_c_alpha_position(residue)
            if c_alpha_position is None:
                residue_index += 1
                continue
            file.write(f'{c_alpha_position[0]:.6f},{c_alpha_position[1]:.6f},{c_alpha_position[2]:.6f},')
            file.write(f'{residue.get_resname()},{residue.id[1]},')
            file.write(','.join([powers_scores[protein_name][residue_index] for protein_name in other_protein_names]))
            file.write('\n')
            residue_index += 1
    print('Powers scores saved')
    with open(f'{output_path}/{pdb_name}.{selected_protein_name}.Residues.csv', 'w') as file:
        file.write('x,y,z,residue,index,' + ','.join(other_protein_names) + '\n')
        residue_index = 0
        for residue in selected_protein.get_residues():
            c_alpha_position = get_c_alpha_position(residue)
            if c_alpha_position is None:
                residue_index += 1
                continue
            file.write(f'{c_alpha_position[0]:.6f},{c_alpha_position[1]:.6f},{c_alpha_position[2]:.6f},')
            file.write(f'{residue.get_resname()},{residue.id[1]},')
            file.write(','.join([residues[protein_name][residue_index] for protein_name in other_protein_names]))
            file.write('\n')
            residue_index += 1
    print('residues saved')
    with open(f'{output_path}/{pdb_name}.{selected_protein_name}.PAM.csv', 'w') as file:
        file.write('x,y,z,residue,index,' + ','.join(other_protein_names) + '\n')
        residue_index = 0
        for residue in selected_protein.get_residues():
            c_alpha_position = get_c_alpha_position(residue)
            if c_alpha_position is None:
                residue_index += 1
                continue
            file.write(f'{c_alpha_position[0]:.6f},{c_alpha_position[1]:.6f},{c_alpha_position[2]:.6f},')
            file.write(f'{residue.get_resname()},{residue.id[1]},')
            file.write(','.join([pam_similarities[protein_name][residue_index] for protein_name in other_protein_names]))
            file.write('\n')
            residue_index += 1
    print('PAM similarities saved')

    with open(f'{output_path}/{pdb_name}.{selected_protein_name}.Differences.csv', 'w') as file:
        file.write('-,' + ','.join([name for name in protein_names]) + '\n')
        for first_protein in differences:
            file.write(f'{first_protein},' + ','.join([differences[first_protein][second_protein_name] for name in protein_names]) + '\n')
    print('raw difference data saved')

    # save parts to separate PDBs with values stored as B factor
    for protein_name in other_protein_names:
        protein = selected_protein
        path = f'{output_path}/{selected_protein_name}-{protein_name}.RMSD.pdb'
        if True:
            for residue in protein.get_residues():
                b_factor = 0.0  # default value for RMSD
                # search for this residue in differences
                for difference in residue_differences[protein_name]:
                    if difference.first_residue_index == residue.id[1]:
                        b_factor = float(difference.distance)
                        break
                for atom in residue.get_atoms():
                    atom.set_bfactor(b_factor)
            io = PDBIO()
            new_structure = Structure(f'{pdb_name}.{protein_name}.RMSD')
            new_model = Model(f'{pdb_name}.{protein_name}.RMSD')
            new_model.add(protein)
            new_structure.add(new_model)
            io.set_structure(new_structure)
            io.save(path)
        path = f'{output_path}/{selected_protein_name}-{protein_name}.Matrix.pdb'
        if True:
            for residue in protein.get_residues():
                b_factor = 1.0  # default value for similarity
                # search for this residue in differences
                for difference in residue_differences[protein_name]:
                    if difference.first_residue_index == residue.id[1]:
                        b_factor = float(difference.similarity)
                        break
                for atom in residue.get_atoms():
                    atom.set_bfactor(b_factor)
            io = PDBIO()
            new_structure = Structure(f'{pdb_name}.{protein_name}.Matrix')
            new_model = Model(f'{pdb_name}.{protein_name}.Matrix')
            new_model.add(protein)
            new_structure.add(new_model)
            io.set_structure(new_structure)
            io.save(path)
        path = f'{output_path}/{selected_protein_name}-{protein_name}.BF.pdb'
        if True:
            for residue in protein.get_residues():
                b_factor = 1.0  # default value for Powers score
                # search for this residue in differences
                for difference in residue_differences[protein_name]:
                    if difference.first_residue_index == residue.id[1]:
                        b_factor = float(difference.powers_score)
                        break
                for atom in residue.get_atoms():
                    atom.set_bfactor(b_factor)
            io = PDBIO()
            new_structure = Structure(f'{pdb_name}.{protein_name}.BF')
            new_model = Model(f'{pdb_name}.{protein_name}.BF')
            new_model.add(protein)
            new_structure.add(new_model)
            io.set_structure(new_structure)
            io.save(path)
        path = f'{output_path}/{selected_protein_name}-{protein_name}.PAM.pdb'
        if True:
            for residue in protein.get_residues():
                b_factor = 0.0  # default value for PAM similarity
                # search for this residue in differences
                for difference in residue_differences[protein_name]:
                    if difference.first_residue_index == residue.id[1]:
                        b_factor = float(difference.pam_similarity)
                        break
                for atom in residue.get_atoms():
                    atom.set_bfactor(b_factor)
            io = PDBIO()
            new_structure = Structure(f'{pdb_name}.{protein_name}.BF')
            new_model = Model(f'{pdb_name}.{protein_name}.BF')
            new_model.add(protein)
            new_structure.add(new_model)
            io.set_structure(new_structure)
            io.save(path)
        print(f'PDBs created for {protein_name}')
