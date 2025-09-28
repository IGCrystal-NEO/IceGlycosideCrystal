"""Utility helpers for igcrystal."""
from __future__ import annotations

from typing import Dict

from .constants import ELEMENT_MASSES


def _normalize_element_symbol(el: str) -> str:
    """Normalize element symbol: 'fe'|'FE' -> 'Fe'."""
    el = str(el).strip()
    if not el:
        raise ValueError("Empty element symbol")
    return el[0].upper() + el[1:].lower() if len(el) > 1 else el.upper()


def extend_element_masses(mapping: Dict[str, float], override: bool = True) -> None:
    """
    Extend the global ELEMENT_MASSES table. If override is True, existing keys will be overwritten.
    Example:
        extend_element_masses({'U': 238.02891})
    """
    for k, v in mapping.items():
        key = _normalize_element_symbol(k)
        if override or key not in ELEMENT_MASSES:
            ELEMENT_MASSES[key] = float(v)
