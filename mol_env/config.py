"""
Configuration file for the Molecular RL Environment
Contains all constants and parameters for environment setup
"""

from typing import NamedTuple
from rdkit.Chem.rdchem import BondType

# ============================================================================
# ACTION CLASS DEFINITION
# ============================================================================
class Action(NamedTuple):
    action_type: str                    # 'add_atom', 'add_bond', 'terminate'
    atom_type: str = None              # 'C', 'N', 'O', 'F' (for add_atom)
    target_idx: int = None             # Target atom index (for add_atom)
    atom1_idx: int = None              # First atom index (for add_bond)
    atom2_idx: int = None              # Second atom index (for add_bond)
    bond_type: BondType = None         # Bond type (RDKit enum)

# ============================================================================
# MOLECULE CONSTRAINTS
# ============================================================================
MAX_ATOMS = 20              # Maximum atoms per molecule
MAX_STEPS = 15              # Maximum steps per episode
MIN_ATOMS = 1               # Minimum atoms (prevent empty molecules)

# ============================================================================
# CHEMISTRY PARAMETERS  
# ============================================================================
ATOM_TYPES = ['C', 'N', 'O', 'F']  # Available atom types for addition
BOND_TYPES = [BondType.SINGLE, BondType.DOUBLE]  # Available bond types (RDKit enums)

# Bond type mapping for string-based logic
BOND_TYPE_MAP = {
    'SINGLE': BondType.SINGLE,
    'DOUBLE': BondType.DOUBLE
}

# Atom valence limits (for validation)
ATOM_VALENCES = {
    'C': 4,  # Carbon max 4 bonds
    'N': 3,  # Nitrogen max 3 bonds  
    'O': 2,  # Oxygen max 2 bonds
    'F': 1,  # Fluorine max 1 bond
}

# ============================================================================
# STARTING MOLECULES
# ============================================================================
STARTING_MOLECULES = [
    'C',        # Methane
    'CC',       # Ethane  
    'c1ccccc1'  # Benzene 
]

# ============================================================================
# ACTION SPACE CONFIGURATION
# ============================================================================
ACTION_TYPES = [
    'add_atom',     # Add new atom with bond
    'add_bond',     # Add bond between existing atoms
    'terminate'     # End episode
]

# ============================================================================
# REWARD WEIGHTS (for later integration)
# ============================================================================
REWARD_WEIGHTS = {
    'qed': 0.4,         # Drug-likeness weight
    'logp': 0.3,        # Lipophilicity weight  
    'sa_score': 0.3,    # Synthetic accessibility weight
    'validity_bonus': 0.1  # Bonus for valid molecules
}

# ============================================================================
# ENVIRONMENT SETTINGS
# ============================================================================
INVALID_ACTION_PENALTY = -1.0    # Penalty for invalid actions
EPISODE_END_REWARD = 0.0         # Reward when episode terminates normally
DEFAULT_STARTING_MOL = 'C'       # Default if random choice fails