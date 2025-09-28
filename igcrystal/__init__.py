"""igcrystal package

Lightweight and robust crystal structure analysis & manipulation toolkit.

Public API re-exports the most commonly used classes and helpers for convenience:
- Crystal
- Atom
- create_simple_cubic
- create_diamond_structure
- extend_element_masses
- ELEMENT_MASSES

This package is refactored from the original single-file IGCrystal.py into decoupled modules.
"""
from .constants import AVOGADRO, DEFAULT_MAX_ATOMS, ELEMENT_MASSES
from .core import Atom, Crystal
from .factories import create_diamond_structure, create_simple_cubic
from .utils import extend_element_masses

__all__ = [
    "Crystal",
    "Atom",
    "create_simple_cubic",
    "create_diamond_structure",
    "extend_element_masses",
    "ELEMENT_MASSES",
    "AVOGADRO",
    "DEFAULT_MAX_ATOMS",
]
