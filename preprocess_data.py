from collections import defaultdict
import os
import pickle
import sys

import numpy as np

from rdkit import Chem
#import torch

def create_atoms(mol):
    """Create a list of atom (e.g., hydrogen and oxygen) IDs
    considering the aromaticity."""
    atoms = [a.GetSymbol() for a in mol.GetAtoms()]
    for a in mol.GetAromaticAtoms():
        i = a.GetIdx()
        atoms[i] = (atoms[i], 'aromatic')
    atoms = [atom_dict[a] for a in atoms]
    return np.array(atoms)


def create_ijbonddict(mol):
    """Create a dictionary, which each key is a node ID
    and each value is the tuples of its neighboring node
    and bond (e.g., single and double) IDs."""
    i_jbond_dict = defaultdict(lambda: [])
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        bond = bond_dict[str(b.GetBondType())]
        i_jbond_dict[i].append((j, bond))
        i_jbond_dict[j].append((i, bond))
    return i_jbond_dict


def extract_fingerprints(atoms, i_jbond_dict, radius):
    """Extract the r-radius subgraphs (i.e., fingerprints)
    from a molecular graph using Weisfeiler-Lehman algorithm."""

    if (len(atoms) == 1) or (radius == 0):
        fingerprints = [fingerprint_dict[a] for a in atoms]

    else:
        nodes = atoms
        i_jedge_dict = i_jbond_dict

        for _ in range(radius):

            """Update each node ID considering its neighboring nodes and edges
            (i.e., r-radius subgraphs or fingerprints)."""
            fingerprints = []
            for i, j_edge in i_jedge_dict.items():
                neighbors = [(nodes[j], edge) for j, edge in j_edge]
                fingerprint = (nodes[i], tuple(sorted(neighbors)))
                fingerprints.append(fingerprint_dict[fingerprint])
            nodes = fingerprints

            """Also update each edge ID considering two nodes
            on its both sides."""
            _i_jedge_dict = defaultdict(lambda: [])
            for i, j_edge in i_jedge_dict.items():
                for j, edge in j_edge:
                    both_side = tuple(sorted((nodes[i], nodes[j])))
                    edge = edge_dict[(both_side, edge)]
                    _i_jedge_dict[i].append((j, edge))
            i_jedge_dict = _i_jedge_dict

    return np.array(fingerprints)


def create_adjacency(mol):
    adjacency = Chem.GetAdjacencyMatrix(mol)
    return np.array(adjacency)


def split_sequence(sequence, ngram):
    sequence = '-' + sequence + '='
    words = [word_dict[sequence[i:i+ngram]]
             for i in range(len(sequence)-ngram+1)]
    return np.array(words)


def dump_dictionary(dictionary, filename):
    with open(filename, 'wb') as f:
        pickle.dump(dict(dictionary), f)


if __name__ == "__main__":

    DATASET, radius, ngram = ['dude', 2, 3]
    radius, ngram = map(int, [radius, ngram])

#    with open('../dataset/' + DATASET + '/original/train_data.txt', 'r') as f:
#        data_list = f.read().strip().split('\n')
#        
    ## 得到训练数据分批的index
#    N_index=[32908,44088, 62040, 89077, 94811, 132417, 148833, 164535, 181534, 188676, 198032, 200996, 216022,
#             223121, 241827, 252176, 284721, 295893, 325019, 329031, 336788, 349091, 361810, 365346, 371245, 391387, 433838,
#             468903, 486853, 523127, 544572, 565480, 571967, 593176, 596088, 601604, 602340, 608445, 662186, 671414, 687161,
#             691311, 703669, 710438, 721042, 731790, 738757, 776430, 796203, 805386, 810453, 815383, 825016, 827405, 836328,
#             843071, 850180, 861041, 864039, 873131, 901670, 911391, 917841, 923274, 934951, 939718, 946618, 983965, 1023011,
#             1031494, 1039801, 1046250, 1051593, 1082764, 1111296, 1122489, 1146425, 1153459, 1160708, 1181083, 1194603, 1221193,
#             1237451, 1245109]
#    for i in range(0,84):
#        with open('E:/DUDE_dataset/'+str(i)+'train_data.txt', 'r') as f:
#            data_train_list = f.read().strip().split('\n')
#            #data_train_list = [d for d in data_train_list if '.' not in d.strip().split()[0]]
#            MM=len(data_train_list)
#            N_index.append(MM)
#        data_train_list=[]

    with open('D:/SERT/SERT_01_data.txt', 'r') as f:
        data_list = f.read().strip().split('\n')
    """Exclude data contains '.' in the SMILES format."""
    data_list = [d for d in data_list if '.' not in d.strip().split()[0]]
    N = len(data_list)
    print(N)

    atom_dict = defaultdict(lambda: len(atom_dict))
    bond_dict = defaultdict(lambda: len(bond_dict))
    fingerprint_dict = defaultdict(lambda: len(fingerprint_dict))
    edge_dict = defaultdict(lambda: len(edge_dict))
    word_dict = defaultdict(lambda: len(word_dict))

#    dir_input=('../dataset/' + DATASET + '/input/for_train/')
#    os.makedirs(dir_input, exist_ok=True)

    Smiles, compounds, adjacencies, proteins, interactions = '', [], [], [], []
    

    for no, data in enumerate(data_list):

        #print('/'.join(map(str, [no+1, N])))

            smiles, sequence, interaction = data.strip().split()
            Smiles += smiles + '\n'
    
            mol = Chem.AddHs(Chem.MolFromSmiles(smiles))  # Consider hydrogens.
            atoms = create_atoms(mol)
            i_jbond_dict = create_ijbonddict(mol)
    
            fingerprints = extract_fingerprints(atoms, i_jbond_dict, radius)
            compounds.append(fingerprints)
    
            adjacency = create_adjacency(mol)
            adjacencies.append(adjacency)
            #print(len(adjacencies))
            words = split_sequence(sequence, ngram)
            proteins.append(words)
            #print(len())
            interactions.append(np.array([float(interaction)]))
            
#            if no+1 in N_index :
#                print('/'.join(map(str, [no+1, N])))
#                np.save(dir_input + 'adjacencies_'+str(N_index.index(no+1)), adjacencies)
#                np.save(dir_input + 'proteins_'+str(N_index.index(no+1)), proteins)
#                #torch.cuda.empty_cache()
#                adjacencies=[]
#                proteins=[]
    
    #            words = split_sequence(sequence, ngram)
    #            proteins.append(words)

#        interactions.append(np.array([float(interaction)]))

#    dir_input = ('../dataset/' + DATASET + '/input/'
#                 'radius' + str(radius) + '_ngram' + str(ngram) + '/')
    dir_input=('D:/SERT/SERT_dataset/')
    os.makedirs(dir_input, exist_ok=True)

    with open(dir_input + 'Smiles.txt', 'w') as f:
        f.write(Smiles)
    np.save(dir_input + 'compounds', compounds)
    np.save(dir_input + 'adjacencies', adjacencies)
    np.save(dir_input + 'proteins', proteins)
    np.save(dir_input + 'interactions', interactions)
    dump_dictionary(fingerprint_dict, dir_input + 'fingerprint_dict.pickle')
    dump_dictionary(word_dict, dir_input + 'word_dict.pickle')

    print('The preprocess of ' + DATASET + ' dataset has finished!')
