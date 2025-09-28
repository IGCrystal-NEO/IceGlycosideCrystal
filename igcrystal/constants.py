"""Constants and default data for igcrystal."""
from __future__ import annotations
from typing import Dict

# Physical constants
AVOGADRO: float = 6.02214076e23  # mol^-1

# Safety limit for supercell generation
DEFAULT_MAX_ATOMS: int = 1_000_000

# Default element mass table (g/mol) - extend or override as needed
ELEMENT_MASSES: Dict[str, float] = {
    "H": 1.00794, "He": 4.002602, "Li": 6.941, "Be": 9.012182,
    "B": 10.811, "C": 12.0107, "N": 14.0067, "O": 15.999,
    "F": 18.9984032, "Ne": 20.1797, "Na": 22.98976928, "Mg": 24.3050,
    "Al": 26.9815386, "Si": 28.0855, "P": 30.973762, "S": 32.065,
    "Cl": 35.453, "Ar": 39.948, "K": 39.0983, "Ca": 40.078,
    "Fe": 55.845, "Cu": 63.546, "Zn": 65.38, "Ag": 107.8682,
    "Au": 196.966569, "Pb": 207.2
}
