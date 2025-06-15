import torch
from torch import Tensor
from rdkit import Chem
from torch_geometric.data import Data

def get_node_features(mol: Chem.Mol) -> Tensor:
    """
    Extract node features from the molecule.
    Each atom is represented by its atomic number.
    """
    mol_features = []
    for atom in mol.GetAtoms():
        atom_features = []

        atom_features.append(atom.GetAtomicNum())
        atom_features.append(atom.GetDegree())
        atom_features.append(atom.GetFormalCharge())
        atom_features.append(atom.GetIsAromatic())
        atom_features.append(atom.GetHybridization().real)
        



    return 


def get_edge_index(mol: Chem.Mol) -> Tensor:
    """
    Extract edge indices from the molecule.
    Each bond is represented by a pair of atom indices.
    """
    return None


def get_edge_attr(mol: Chem.Mol) -> Tensor:
    return None


def mol_to_graph(mol: Chem.Mol) -> Data:
    """
    Convert a molecule to a graph representation.
    """
    node_features = get_node_features(mol)
    edge_index = get_edge_index(mol)
    edge_attr = get_edge_attr(mol)

    return Data(x=node_features, edge_index=edge_index, edge_attr=edge_attr)