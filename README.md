# IceGlycosideCrystal

**Lightweight and Robust Crystal Structure Analysis & Manipulation Toolkit (Python)** — for representing unit cells, coordinate transformations, volume/density calculations, supercell generation, visualization, and simple I/O.

Note: The codebase has been refactored from a single-file script into a small, decoupled package at `igcrystal/`, while keeping `IGCrystal.py` as a backward-compatible CLI/demo wrapper.



## Key Features

* Fractional ↔ Cartesian coordinate conversion for arbitrary triclinic lattices (numerically stable).
* Unit-cell volume and density calculations: either via provided (`molecular_weight, Z`) or computed automatically from atomic masses inside the cell (respecting `occupancy`).
* Supercell generation (`generate_supercell`) with a safety cap to prevent OOM.
* Plotting (matplotlib): supports saving to file; points/bonds/unit-cell edges optional. When available, `scipy.spatial.cKDTree` is used to accelerate bond searching.
* JSON serialization/deserialization (backward compatible with old atom format `[el, x, y, z]` and new format `[el, x, y, z, occupancy]`).
* Automatic normalization of element symbols (case-tolerant) and support for extending element mass tables either per-instance or globally.



## Requirements

* Python 3.8+ (3.9 / 3.10 recommended)
* Required: `numpy`, `matplotlib`
* Optional (for accelerated bond search): `scipy` (for `scipy.spatial.cKDTree`)

Installation (example):

```bash
pip install numpy matplotlib
# optional
pip install scipy
```



## Quick Start

Use the package directly in Python:

```python
from igcrystal import create_diamond_structure, Crystal

# create a diamond structure and inspect
diamond = create_diamond_structure()
print(diamond)

# volume
print("Volume (Å^3):", diamond.get_volume())

# density computed from atomic masses (using built-in ELEMENT_MASSES)
print("Density (auto):", diamond.get_density(), "g/cm^3")

# density from molecular_weight and Z (example)
# if formula unit = C (12.01 g/mol) and the cell contains 8 C atoms:
print("Density (molecular_weight & Z):", diamond.get_density(molecular_weight=12.01, Z=8))

# generate a 2x2x2 supercell (note safety cap)
supercell = diamond.generate_supercell(2, 2, 2)
print("Supercell atoms:", len(supercell.atoms))

# plot (save to file — avoids blocking on headless servers)
diamond.plot_structure(save_as="diamond.png", show=False)

# save / load
diamond.save_to_file("diamond_example.json")
loaded = Crystal.load_from_file("diamond_example.json")
```

Or run the legacy wrapper (same demo output):

```bash
python IGCrystal.py
```



## Main API (Summary)

### Classes & Constructors

* `Crystal(name: str = "Unknown Crystal", mass_table: Optional[Dict[str, float]] = None)`

  * `mass_table`: optional instance-level element mass table (keys are auto-normalized).

### Lattice & Atoms

* `set_lattice_parameters(a, b, c, alpha=90.0, beta=90.0, gamma=90.0)`

  * Parameters must be positive; angles in (0, 180) degrees.
* `add_atom(element, x, y, z, occupancy=1.0)`

  * `element` is auto-formatted with initial capital letter. `occupancy` in `[0, 1]`.

### Calculations

* `get_volume()` → returns unit cell volume in Å³ (float).
* `get_density(molecular_weight: Optional[float], Z: Optional[int])` → returns density in g/cm³.

  * If `molecular_weight` and `Z` are both provided: use `density = Z * M / (N_A * V)`. Otherwise sum atomic masses in the cell (weighted by occupancy); mass table comes from instance `mass_table` (if provided) or module default `ELEMENT_MASSES`.

### Coordinate Transforms

* `fractional_to_cartesian((u, v, w))` → returns Cartesian coordinates in Å.
* `cartesian_to_fractional((x, y, z))` → returns fractional coordinates (uses `np.linalg.solve` for numeric stability).

### Supercell & Normalization

* `generate_supercell(nx, ny, nz, max_atoms=DEFAULT_MAX_ATOMS)` → returns a `Crystal` supercell object.

  * If the estimated atom count exceeds `max_atoms` (default `1_000_000`), a `ValueError` is raised.
* `normalize_fractional()` → reduces fractional coordinates modulo 1 into `[0, 1)`.

### Plotting & Bonds

* `plot_structure(ax=None, show_bonds=False, show_unit_cell=True, bond_cutoff=1.6, save_as=None, show=True, element_colors=None)`

  * `save_as`: save image file if provided.
  * `show=False`: does not call `plt.show()` (suitable for headless environments).
  * If `scipy` is installed, bond detection will use `cKDTree` for much faster neighbor searches.

### I/O

* `to_dict()` / `save_to_file(filename, verbose=True)`
* `Crystal.load_from_file(filename, mass_table=None)` — loads files and is compatible with both old/new atom formats.

### Global / Utilities

* `extend_element_masses(mapping: Dict[str, float], override: bool = True)` — extend/override module-level `ELEMENT_MASSES`.



## JSON Format Example

Current `to_dict()` output (example):

```json
{
  "name": "Diamond",
  "lattice_parameters": {"a": 5.43, "b": 5.43, "c": 5.43, "alpha": 90, "beta": 90, "gamma": 90},
  "space_group": "Fd-3m",
  "atoms": [
    ["C", 0.0, 0.0, 0.0, 1.0],
    ["C", 0.25, 0.25, 0.25, 1.0]
    // ... each atom includes occupancy (old format [el,x,y,z] is still accepted)
  ],
  "format_version": 1
}
```



## Troubleshooting & Common Issues

* **`ValueError: gamma angle too close to 0 or 180 degrees`**

  * Means `gamma` is near 0°/180°, making `sin(gamma)` extremely small or zero — the lattice vector matrix degenerates. Check input angles.

* **`ValueError: Lattice vectors matrix is singular` (cartesian→fractional)**

  * Lattice vectors are linearly dependent or numerically ill-conditioned. Verify `a,b,c` and angle combinations.

* **`ValueError: Missing atomic masses for elements`**

  * Some element masses required for automatic density computation are missing from the mass table. Add masses with `extend_element_masses`, pass an instance `mass_table` when creating `Crystal`, or compute density with `get_density(molecular_weight, Z)`.

* **Supercell generation errors exceeding `max_atoms`**

  * Reduce `nx, ny, nz` or increase `max_atoms` (careful), or use more advanced libraries (ASE, pymatgen) for very large systems.

* **Plotting blocks/errors on headless servers**

  * Use `plot_structure(save_as="out.png", show=False)` to save without displaying, or set matplotlib backend to `Agg` in the environment.



## Advanced Suggestions

* To interoperate with the materials-science ecosystem (pymatgen, ASE, spglib), implement `to_cif()` / `from_cif()` and space-group standardization (hook into `spglib`).
* To support more elements and isotopic precision, expand `ELEMENT_MASSES` or depend on packages like `mendeleev` / `periodictable`.
* For neighbor/bond detection on large systems, install `scipy` and use `cKDTree` for significant performance gains.

---

## Testing

* Add `pytest` unit tests (recommended): volume, density, coordinate inverses, supercell counts, JSON compatibility, etc.

---
### Project Layout (after refactor)

```
igcrystal/
  __init__.py      # public API exports (Crystal, Atom, factories, constants, utils)
  constants.py     # AVOGADRO, DEFAULT_MAX_ATOMS, ELEMENT_MASSES
  utils.py         # helpers like _normalize_element_symbol, extend_element_masses
  core.py          # Atom, Crystal core logic (volume, density, transforms, supercell, I/O, plotting)
  factories.py     # convenience constructors (e.g., create_diamond_structure)
IGCrystal.py       # backward-compatible CLI/demo wrapper
```


