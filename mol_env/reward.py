from rdtkit.Chem import QED, Descriptor
from rdtkit.Chem import Mol
from config import REWARD_WEIGHTS

def qed_reward(mol: Mol) -> float:
    """
    Calculate the QED (Quantitative Estimation of Drug-likeness) reward for a molecule.
    The QED value is a measure of drug-likeness, with higher values indicating better drug-like properties.
    """
    qed_value = QED.qed(mol)
    return qed_value

def logp_reward(mol: Mol) -> float:
    """
    Calculate the LogP (Partition Coefficient) reward for a molecule.
    The LogP value is a measure of hydrophobicity, with higher values indicating more hydrophobic molecules.
    """
    logp_value = Descriptor.MolLogP(mol)
    return logp_value

def sa_reward(mol: Mol) -> float:
    """
    Calculate the Synthetic Accessibility (SA) reward for a molecule.
    The SA score indicates how easy it is to synthesize the molecule, with lower values indicating easier synthesis.
    """
    sa_value = Descriptor.SAIndex(mol)
    return sa_value



def compute_reward(mol: Mol) -> float:
    """
    Compute the total reward for a molecule.
    The total reward is a weighted sum of QED, LogP, and SA scores.
    """
    qed_reward_value = qed_reward(mol)
    logp_reward_value = logp_reward(mol)
    sa_reward_value = sa_reward(mol)

    total_reward = (
        REWARD_WEIGHTS['qed'] * qed_reward_value +
        REWARD_WEIGHTS['logp'] * logp_reward_value +
        REWARD_WEIGHTS['sa_score'] * sa_reward_value
    )

    return total_reward