from typing import NamedTuple, List
from rdkit.Chem import Mol, EditableMol, Atom
from rdkit.Chem.rdchem import BondType

from state import copy_molecule, sanitize_molecule

# This really goes on the config file, but for now we keep it here
ALLOWED_ATOMS: List[str] = ["C", "N", "O", "F"]
ALLOWED_BONDS: List[BondType] = [BondType.SINGLE, BondType.DOUBLE]

class Action(NamedTuple):
    # Lets say for now we only allow adding atoms
    target_index: int
    added_atom: str
    bond_type: BondType = None

def get_valid_actions(mol: Mol) -> List[Action]:
    """
    Generate a list of valid actions for the given molecule.
    """
    valid_actions = []
    for atom_indx, _ in enumerate(mol.GetAtoms()):
        for atom in ALLOWED_ATOMS:
            for bond in ALLOWED_BONDS:
                candidate_action = Action(target_index=atom_indx, added_atom = atom, bond_type=bond)
                if is_action_valid(mol, candidate_action):
                    valid_actions.append(candidate_action)

    # Add terminate action
    valid_actions.append(Action(target_index=-1, added_atom="", bond_type=None))
    return valid_actions


def apply_action(mol: Mol, action: Action) -> Mol:
    """
    Apply the given action to the molecule and return the resulting molecule.
    """
    if action.target_index < 0:
        #Terminate
        return copy_molecule(mol)
    else:
        new_mol = copy_molecule(mol)
        editable_mol = EditableMol(new_mol)
        new_index = editable_mol.AddAtom(Atom(action.added_atom))
        editable_mol.AddBond(action.target_index, new_index, action.bond_type)
        return editable_mol.GetMol()


def is_action_valid(mol: Mol, action: Action) -> bool:
    """
    Check if the action is valid for the given molecule.
    """
    test_mol = copy_molecule(mol)
    test_mol = apply_action(test_mol, action)
    
    return sanitize_molecule(test_mol)


# Optionally I could create a function for debugging/logging the actions