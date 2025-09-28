"""Factory helpers to create common crystal prototypes."""
from __future__ import annotations
from .core import Crystal


def create_simple_cubic(element: str = "Si", lattice_constant: float = 5.43) -> Crystal:
    crystal = Crystal(f"Simple cubic {element}")
    crystal.set_lattice_parameters(lattice_constant, lattice_constant, lattice_constant)
    crystal.add_atom(element, 0.0, 0.0, 0.0)
    return crystal


def create_diamond_structure(lattice_constant: float = 5.43) -> Crystal:
    crystal = Crystal("Diamond")
    crystal.set_lattice_parameters(lattice_constant, lattice_constant, lattice_constant)
    crystal.space_group = "Fd-3m"
    positions = [
        (0.0, 0.0, 0.0), (0.25, 0.25, 0.25),
        (0.5, 0.5, 0.0), (0.75, 0.75, 0.25),
        (0.5, 0.0, 0.5), (0.75, 0.25, 0.75),
        (0.0, 0.5, 0.5), (0.25, 0.75, 0.75)
    ]
    for x, y, z in positions:
        crystal.add_atom("C", x, y, z)
    return crystal
