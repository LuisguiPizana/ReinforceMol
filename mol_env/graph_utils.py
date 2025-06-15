import torch
from torch import Tensor
from rdkit import Chem
from torch_geometric.data import Data
from config import ATOM_TYPES, BOND_TYPES
from typing import TypeVar, Sequence

# Type variable for generic type hinting
T = TypeVar('T')

def one_hot(value: T, category_list: Sequence[T]) -> list:
    """
    Convert a string value to a one-hot encoded list. It uses the global ATOM_TYPES
    to create the one-hot encoding.
    """
    if value not in category_list:
        raise ValueError(f"Value '{value}' not found in category list: {category_list}")

    return [1 if value == category else 0 for category in category_list]


def get_node_features(mol: Chem.Mol) -> Tensor:
    """
    Extract each node/atom features from the molecule.
    """
    mol_features = []
    for atom in mol.GetAtoms():
        atom_features = []

        atom_features += one_hot(atom.GetSymbol(), category_list=ATOM_TYPES)
        atom_features.append(atom.GetDegree())
        atom_features.append(atom.GetFormalCharge())
        atom_features.append(int(atom.GetIsAromatic()))

        mol_features.append(atom_features)
    x = torch.tensor(mol_features, dtype = torch.float)

    return x


def get_edge_index(mol: Chem.Mol) -> Tensor:
    """
    Extract edge indices from the molecule.
    Each bond is represented by a pair of atom indices.
    """
    edge_index_list = []
    for bond in mol.GetBonds():
        start_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        edge_index_list.append([start_idx, end_idx])
        edge_index_list.append([end_idx, start_idx])
    
    edge_index = torch.tensor(edge_index_list, dtype=torch.long)
    
    if edge_index.numel() > 0:
        return edge_index.t().contiguous()
    else:
        # If there are no edges, return an empty tensor
        return torch.empty((2, 0), dtype=torch.int32)



def get_edge_attr(mol: Chem.Mol) -> Tensor:
    mol_bond_attr = []
    for bond in mol.GetBonds():
        bond_attr = []
        bond_attr += one_hot(bond.GetBondType(), category_list=BOND_TYPES)
        bond_attr.append(int(bond.GetIsAromatic()))

        # Add bond features for both directions (start to end and end to start)
        mol_bond_attr.append(bond_attr)
        mol_bond_attr.append(bond_attr)

    bond_attr_tensor = torch.tensor(mol_bond_attr, dtype=torch.float)

    if bond_attr_tensor.numel() > 0:
        return bond_attr_tensor
    else:
        # Important to adjust for the number of features
        return torch.empty((0, len(BOND_TYPES) + 1), dtype=torch.float)


def mol_to_graph(mol: Chem.Mol) -> Data:
    """
    Convert a molecule to a graph representation.
    """
    node_features = get_node_features(mol)
    edge_index = get_edge_index(mol)
    edge_attr = get_edge_attr(mol)

    return Data(x=node_features, edge_index=edge_index, edge_attr=edge_attr)