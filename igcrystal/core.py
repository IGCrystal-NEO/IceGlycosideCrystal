"""Core data structures and operations for igcrystal."""
from __future__ import annotations
from dataclasses import dataclass
from typing import List, Tuple, Optional, Dict, Any
import numpy as np
import json
import math
import logging

from .constants import AVOGADRO, DEFAULT_MAX_ATOMS, ELEMENT_MASSES
from .utils import _normalize_element_symbol

# Attempt to import fast neighbor search (optional)
try:
    from scipy.spatial import cKDTree  # type: ignore
    _HAS_CKDTree = True
except Exception:
    _HAS_CKDTree = False

logger = logging.getLogger(__name__)


@dataclass
class Atom:
    element: str
    x: float  # fractional coordinate
    y: float
    z: float
    occupancy: float = 1.0

    def to_list(self) -> List[Any]:
        # New format includes occupancy for round-trip compatibility.
        return [self.element, float(self.x), float(self.y), float(self.z), float(self.occupancy)]


class Crystal:
    """
    Represents a crystal cell with lattice parameters and atoms (fractional coordinates).

    Parameters
    ----------
    name : str
        Human-friendly name of the structure.
    mass_table : Optional[Dict[str, float]]
        Per-instance element mass table (g/mol). If provided, used in preference to module-level ELEMENT_MASSES.
    """

    def __init__(self, name: str = "Unknown Crystal", mass_table: Optional[Dict[str, float]] = None):
        self.name: str = name
        self.lattice_parameters = {
            "a": 1.0, "b": 1.0, "c": 1.0,
            "alpha": 90.0, "beta": 90.0, "gamma": 90.0
        }
        self.atoms: List[Atom] = []
        self.space_group: str = "P1"
        # instance mass table (optional). Normalize keys if provided.
        if mass_table:
            self.mass_table = { _normalize_element_symbol(k): float(v) for k,v in mass_table.items() }
        else:
            self.mass_table = None

    # -------------------------
    # Basic setters / validation
    # -------------------------
    def set_lattice_parameters(self, a: float, b: float, c: float,
                               alpha: float = 90.0, beta: float = 90.0, gamma: float = 90.0) -> None:
        if a <= 0 or b <= 0 or c <= 0:
            raise ValueError("Lattice constants must be positive.")
        for ang in (alpha, beta, gamma):
            if ang <= 0 or ang >= 180:
                raise ValueError("Lattice angles must be in (0, 180) degrees.")
        self.lattice_parameters = {
            "a": float(a), "b": float(b), "c": float(c),
            "alpha": float(alpha), "beta": float(beta), "gamma": float(gamma)
        }

    def add_atom(self, element: str, x: float, y: float, z: float, occupancy: float = 1.0) -> None:
        """Add an atom with fractional coordinates. occupancy in [0,1]."""
        el = _normalize_element_symbol(element)
        if not (0.0 <= occupancy <= 1.0):
            raise ValueError("occupancy must be between 0.0 and 1.0")
        self.atoms.append(Atom(el, float(x), float(y), float(z), float(occupancy)))

    # -------------------------
    # Lattice vector construction
    # -------------------------
    def _lattice_vectors(self) -> np.ndarray:
        """
        Returns a 3x3 matrix whose columns are the lattice vectors (in Å).
        Convention:
         - alpha is angle between b and c
         - beta is angle between a and c
         - gamma is angle between a and b
        """
        a = self.lattice_parameters["a"]
        b = self.lattice_parameters["b"]
        c = self.lattice_parameters["c"]
        alpha = math.radians(self.lattice_parameters["alpha"])
        beta = math.radians(self.lattice_parameters["beta"])
        gamma = math.radians(self.lattice_parameters["gamma"])

        cos_alpha = math.cos(alpha)
        cos_beta = math.cos(beta)
        cos_gamma = math.cos(gamma)
        sin_gamma = math.sin(gamma)
        if abs(sin_gamma) < 1e-12:
            raise ValueError("gamma angle too close to 0 or 180 degrees; sin(gamma) ≈ 0.")

        a_vec = np.array([a, 0.0, 0.0])
        b_vec = np.array([b * cos_gamma, b * sin_gamma, 0.0])

        # c vector components (robust)
        c_x = c * cos_beta
        term = (cos_alpha - cos_beta * cos_gamma) / sin_gamma
        c_y = c * term
        c_z_sq = 1.0 - cos_beta**2 - term**2
        c_z_sq = max(c_z_sq, 0.0)
        c_z = c * math.sqrt(c_z_sq)
        c_vec = np.array([c_x, c_y, c_z])

        mat = np.column_stack((a_vec, b_vec, c_vec))
        return mat

    # -------------------------
    # Volume and density
    # -------------------------
    def get_volume(self) -> float:
        """Compute unit cell volume in cubic Angstroms (Å^3)."""
        mat = self._lattice_vectors()
        volume = abs(np.linalg.det(mat))
        if not np.isfinite(volume) or volume <= 0.0:
            raise ValueError("Calculated unit cell volume is non-positive or invalid. Check lattice parameters.")
        return float(volume)

    def get_density(self, molecular_weight: Optional[float] = None, Z: Optional[int] = None) -> float:
        """
        Compute density in g/cm^3.

        Two modes:
        1) molecular_weight (g/mol of a formula unit) + Z (number of formula units per cell).
           density = Z * M / (N_A * V)
        2) otherwise, sum atomic masses in the unit cell (from mass table) weighted by occupancy.
        """
        V_A3 = self.get_volume()  # Å^3
        volume_cm3 = V_A3 * 1e-24  # Å^3 -> cm^3

        if molecular_weight is not None and Z is not None:
            if molecular_weight <= 0 or Z <= 0:
                raise ValueError("molecular_weight and Z must be positive.")
            mass_per_cell_g = (float(molecular_weight) * float(Z)) / AVOGADRO
        else:
            if len(self.atoms) == 0:
                raise ValueError("No atoms in cell and molecular_weight/Z not provided; cannot compute density.")
            total_mass = 0.0
            missing = set()
            # choose instance table if given else module-level
            mass_src = self.mass_table if self.mass_table is not None else ELEMENT_MASSES
            for atom in self.atoms:
                el = _normalize_element_symbol(atom.element)
                if el in mass_src:
                    total_mass += mass_src[el] * getattr(atom, "occupancy", 1.0)
                else:
                    missing.add(el)
            if missing:
                raise ValueError(
                    f"Missing atomic masses for elements: {sorted(list(missing))}. "
                    "Provide molecular_weight & Z, extend mass table, or pass mass_table to Crystal."
                )
            mass_per_cell_g = total_mass / AVOGADRO

        density = mass_per_cell_g / volume_cm3
        return float(density)

    # -------------------------
    # Coordinate transforms
    # -------------------------
    def fractional_to_cartesian(self, frac: Tuple[float, float, float]) -> np.ndarray:
        """Convert fractional (u,v,w) to Cartesian coordinates in Å."""
        mat = self._lattice_vectors()
        frac_vec = np.array(frac, dtype=float)
        return mat @ frac_vec

    def cartesian_to_fractional(self, cart: Tuple[float, float, float]) -> np.ndarray:
        """Convert Cartesian coordinates (Å) to fractional (u,v,w).

        Uses np.linalg.solve for numerical stability instead of computing inverse explicitly.
        """
        mat = self._lattice_vectors()
        vec = np.array(cart, dtype=float)
        try:
            frac = np.linalg.solve(mat, vec)
        except np.linalg.LinAlgError:
            raise ValueError("Lattice vectors matrix is singular or ill-conditioned; cannot solve. Check lattice parameters.")
        return frac

    # -------------------------
    # Supercell and normalization
    # -------------------------
    def normalize_fractional(self) -> None:
        """Normalize fractional coordinates into [0,1) by modulo 1."""
        for atom in self.atoms:
            atom.x = atom.x % 1.0
            atom.y = atom.y % 1.0
            atom.z = atom.z % 1.0

    def generate_supercell(self, nx: int, ny: int, nz: int, max_atoms: int = DEFAULT_MAX_ATOMS) -> "Crystal":
        """
        Generate a supercell with nx x ny x nz repetition.

        Raises ValueError if the estimated number of atoms exceeds max_atoms.
        """
        if not (isinstance(nx, int) and isinstance(ny, int) and isinstance(nz, int)):
            raise ValueError("nx, ny, nz must be integers.")
        if nx < 1 or ny < 1 or nz < 1:
            raise ValueError("nx, ny, nz must be >= 1.")

        estimated = nx * ny * nz * max(1, len(self.atoms))
        if estimated > max_atoms:
            raise ValueError(f"Requested supercell would create ~{estimated} atoms, exceeding safe limit {max_atoms}.")

        supercell = Crystal(f"{self.name}_supercell_{nx}x{ny}x{nz}", mass_table=self.mass_table)
        lp = self.lattice_parameters.copy()
        lp["a"] *= nx
        lp["b"] *= ny
        lp["c"] *= nz
        supercell.lattice_parameters = lp
        supercell.space_group = self.space_group

        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    for atom in self.atoms:
                        new_x = (atom.x + i) / nx
                        new_y = (atom.y + j) / ny
                        new_z = (atom.z + k) / nz
                        supercell.add_atom(atom.element, new_x, new_y, new_z, getattr(atom, "occupancy", 1.0))

        supercell.normalize_fractional()
        return supercell

    # -------------------------
    # Plotting & bonds
    # -------------------------
    def plot_structure(self, ax=None, show_bonds: bool = False, show_unit_cell: bool = True,
                       bond_cutoff: float = 1.6, save_as: Optional[str] = None, show: bool = True,
                       element_colors: Optional[Dict[str, str]] = None):
        """
        Plot the crystal in Cartesian coordinates (Å) using matplotlib 3D.

        - save_as: optional filename to save the figure (png, svg, ...)
        - show: whether to call plt.show(). If False, the figure will be closed (safe for headless env).
        - element_colors: optional mapping element->color to override defaults.
        """
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

        cart_coords = []
        elements = []
        for atom in self.atoms:
            cart = self.fractional_to_cartesian((atom.x, atom.y, atom.z))
            cart_coords.append(cart)
            elements.append(atom.element)
        cart_coords = np.array(cart_coords) if cart_coords else np.zeros((0, 3))

        created_fig = False
        if ax is None:
            fig = plt.figure(figsize=(8, 6))
            ax = fig.add_subplot(111, projection='3d')
            created_fig = True

        # sensible default colors (avoid pure white for visibility)
        default_colors = {
            'H': 'lightgray', 'C': 'black', 'N': 'blue', 'O': 'red',
            'Si': 'gray', 'Al': 'silver', 'Fe': 'orange', 'Ca': 'green'
        }
        if element_colors:
            default_colors.update(element_colors)

        unique_elements = sorted(set(elements), key=lambda e: elements.index(e))
        for el in unique_elements:
            idx = [i for i, e in enumerate(elements) if e == el]
            pts = cart_coords[idx]
            if pts.size == 0:
                continue
            color = default_colors.get(el, 'purple')
            ax.scatter(pts[:, 0], pts[:, 1], pts[:, 2], s=60, alpha=0.9, label=el, color=color)

        ax.set_xlabel('X (Å)')
        ax.set_ylabel('Y (Å)')
        ax.set_zlabel('Z (Å)')
        ax.set_title(f'Crystal Structure: {self.name}')

        # bonds: prefer cKDTree for performance if available
        if show_bonds and cart_coords.shape[0] > 1:
            if _HAS_CKDTree:
                try:
                    tree = cKDTree(cart_coords)
                    pairs = tree.query_pairs(r=bond_cutoff)
                except Exception:
                    pairs = set()
                    # fallback to naive
                    for i in range(len(cart_coords)):
                        for j in range(i + 1, len(cart_coords)):
                            if np.linalg.norm(cart_coords[i] - cart_coords[j]) <= bond_cutoff:
                                pairs.add((i, j))
            else:
                # naive O(N^2) fallback
                pairs = set()
                for i in range(len(cart_coords)):
                    for j in range(i + 1, len(cart_coords)):
                        if np.linalg.norm(cart_coords[i] - cart_coords[j]) <= bond_cutoff:
                            pairs.add((i, j))
            for i, j in pairs:
                xs = [cart_coords[i, 0], cart_coords[j, 0]]
                ys = [cart_coords[i, 1], cart_coords[j, 1]]
                zs = [cart_coords[i, 2], cart_coords[j, 2]]
                ax.plot(xs, ys, zs, linewidth=1)

        # unit cell edges
        if show_unit_cell:
            mat = self._lattice_vectors()
            origin = np.zeros(3)
            corners = [
                origin,
                mat[:, 0],
                mat[:, 1],
                mat[:, 2],
                mat[:, 0] + mat[:, 1],
                mat[:, 0] + mat[:, 2],
                mat[:, 1] + mat[:, 2],
                mat[:, 0] + mat[:, 1] + mat[:, 2]
            ]
            edges = [
                (0, 1), (0, 2), (0, 3),
                (1, 4), (1, 5),
                (2, 4), (2, 6),
                (3, 5), (3, 6),
                (4, 7), (5, 7), (6, 7)
            ]
            for e0, e1 in edges:
                p0 = corners[e0]
                p1 = corners[e1]
                ax.plot([p0[0], p1[0]], [p0[1], p1[1]], [p0[2], p1[2]], linestyle='--', linewidth=0.8)

        ax.legend()
        import matplotlib.pyplot as plt
        plt.tight_layout()
        if save_as:
            plt.savefig(save_as)
            logger.info("Saved plot to %s", save_as)
        if show:
            plt.show()
        else:
            if created_fig:
                plt.close()

    # -------------------------
    # I/O
    # -------------------------
    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "lattice_parameters": self.lattice_parameters,
            "space_group": self.space_group,
            "atoms": [atom.to_list() for atom in self.atoms],
            "format_version": 1
        }

    def save_to_file(self, filename: str, verbose: bool = True) -> str:
        """Save structure to JSON. Returns filename."""
        with open(filename, "w", encoding="utf-8") as f:
            json.dump(self.to_dict(), f, indent=2, ensure_ascii=False)
        if verbose:
            logger.info("Crystal structure saved to %s", filename)
        return filename

    @classmethod
    def load_from_file(cls, filename: str, mass_table: Optional[Dict[str, float]] = None) -> "Crystal":
        """Load from JSON file. Accepts both old ([el,x,y,z]) and new ([el,x,y,z,occupancy]) atom formats."""
        with open(filename, "r", encoding="utf-8") as f:
            data = json.load(f)
        name = data.get("name", "Loaded Crystal")
        crystal = cls(name, mass_table=mass_table)
        lp = data.get("lattice_parameters")
        if lp:
            try:
                crystal.set_lattice_parameters(lp["a"], lp["b"], lp["c"],
                                               lp.get("alpha", 90.0), lp.get("beta", 90.0), lp.get("gamma", 90.0))
            except Exception as e:
                raise ValueError(f"Invalid lattice_parameters in file '{filename}': {e}")
        crystal.space_group = data.get("space_group", crystal.space_group)
        atoms = data.get("atoms", [])
        for atom_data in atoms:
            # support [el,x,y,z] and [el,x,y,z,occupancy]
            if len(atom_data) >= 4:
                el = atom_data[0]
                x = atom_data[1]
                y = atom_data[2]
                z = atom_data[3]
                occ = float(atom_data[4]) if len(atom_data) > 4 else 1.0
                crystal.add_atom(el, x, y, z, occupancy=occ)
        return crystal
