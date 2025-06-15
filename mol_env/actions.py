from .config import Action
from .action_space import create_action_list
from .actions import apply_action, copy_molecule  # Use Person A's functions
from rdkit import Chem
from typing import List

def get_valid_action_ids(mol) -> List[int]:
    """Return list of valid action IDs for current molecule using sanitization"""
    valid_ids = []
    action_list = create_action_list()
    
    for action_id, action in enumerate(action_list):
        if is_action_valid_by_sanitization(mol, action):
            valid_ids.append(action_id)
    
    return valid_ids

def is_action_valid_by_sanitization(mol, action: Action) -> bool:
    """Test if action is valid by trying it and checking sanitization"""
    try:
        # Try applying the action
        test_mol = copy_molecule(mol)
        result_mol = apply_action(test_mol, action)
        
        # Test sanitization
        Chem.SanitizeMol(result_mol)
        return True
    except:
        return False