from rdkit import Chem

"""
This module provides functions to create molecules based on SMILES strings, copy them and validate
their chemical structure.

Its supposed to be used in:
env.py: 
    -On reset() -> calls initialize_molecule()
    -After step() -> calls copy_molecule() and sanitize_molecule()

actions.py:
    -Before applying graph edits (to copy and keep the original molecule intact)
    
reward.py:
    -To validate the molecule before calculating the reward
"""

def initialize_molecule(seed: str = None) -> Chem.Mol:
    if seed is None:
        mol = Chem.MolFromSmiles("C") # Default to single carbon atom
    else:
        mol = Chem.MolFromSmiles(seed)

    return mol

def copy_molecule(mol: Chem.Mol) -> Chem.Mol:
    if mol is None:
        raise ValueError("Input molecule is None")
    return Chem.Mol(mol)

def sanitize_molecule(mol: Chem.Mol) -> bool:
    try:
        Chem.SanitizeMol(mol)
        return True
    except Exception as e:
        print(f"Sanitization failed: {e}")
        return False