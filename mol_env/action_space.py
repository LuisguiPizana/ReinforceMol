from .config import ATOM_TYPES, BOND_TYPES, MAX_ATOMS, ACTION_TYPES
from typing import List, Tuple, Dict, Any


def create_action_list() -> List[Tuple]:
    action_list = []


    # 1. Add atom
    for atom_type in ATOM_TYPES:
        for target_idx in range(MAX_ATOMS):
            for bond_type in BOND_TYPES:
                action_list.append(('add_atom', atom_type, target_idx, bond_type))

    # 2. Add bond
    for atom1_idx in range(MAX_ATOMS):
        for atom2_idx in range(MAX_ATOMS):
            for bond_type in BOND_TYPES:
                if atom1_idx != atom2_idx:
                    action_list.append(('add_bond', atom1_idx, atom2_idx, bond_type))

    # 3. Terminate action
    action_list.append(('terminate',))

    return action_list

def decode_action(action: int) -> Tuple:
    action_list = create_action_list()
    return action_list[action]

def encode_action(action: Tuple) -> int:
    action_list = create_action_list()
    return action_list.index(action)

def get_action_space_size() -> int:
    """Return total number of possible actions"""
    return len(create_action_list())