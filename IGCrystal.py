#!/usr/bin/env python3
"""Compatibility wrapper for the decoupled igcrystal package.

This script keeps the previous entry-point behavior while internally importing
from the new package modules. Users can continue running `python IGCrystal.py`
for the demo CLI, or import `igcrystal` as a package.
"""
from __future__ import annotations
import logging
from igcrystal import (
    Crystal,
    create_diamond_structure,
)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def main() -> None:
    logger.info("IGCrystal - Crystal Structure Toolkit (package wrapper)")
    logger.info("=" * 60)

    diamond = create_diamond_structure()
    logger.info("\n%s", diamond)

    # Supercell example
    try:
        supercell = diamond.generate_supercell(2, 2, 2)
        logger.info("Supercell has %d atoms and volume %.3f Å³", len(supercell.atoms), supercell.get_volume())
    except Exception as e:
        logger.error("Supercell generation failed: %s", e)

    # Density examples
    try:
        density_auto = diamond.get_density()
        logger.info("Diamond density (from atomic masses): %.3f g/cm³", density_auto)
    except Exception as e:
        logger.error("Density (auto) failed: %s", e)

    carbon_mw = 12.01
    try:
        density_mz = diamond.get_density(molecular_weight=carbon_mw, Z=len(diamond.atoms))
        logger.info("Diamond density (molecular_weight & Z=len(atoms)): %.3f g/cm³", density_mz)
    except Exception as e:
        logger.error("Density (molecular_weight & Z) failed: %s", e)

    # Save / Load demo
    diamond.save_to_file("diamond_example.json")
    loaded = Crystal.load_from_file("diamond_example.json")
    logger.info("Loaded: %s, atoms: %d", loaded.name, len(loaded.atoms))


if __name__ == "__main__":
    main()

