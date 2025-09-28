import os
import json
import tempfile
import unittest
import numpy as np

from igcrystal import (
    Crystal,
    create_diamond_structure,
)


class TestIGCrystal(unittest.TestCase):
    def test_volume_cubic(self):
        diamond = create_diamond_structure()  # a=b=c=5.43 by default in this repo
        vol = diamond.get_volume()
        self.assertAlmostEqual(vol, 5.43 ** 3, places=3)

    def test_density_auto_diamond(self):
        diamond = create_diamond_structure()
        rho = diamond.get_density()
        # with a=5.43 (Si lattice constant), density ~ 1.0 g/cm^3 for 8 C atoms
        self.assertAlmostEqual(rho, 0.997, places=3)

    def test_frac_cart_invertibility_triclinic(self):
        c = Crystal("triclinic")
        c.set_lattice_parameters(a=5.0, b=6.0, c=7.0, alpha=80.0, beta=100.0, gamma=110.0)
        c.add_atom("C", 0.1, 0.2, 0.3)
        for vec in [(0.1, 0.2, 0.3), (0.9, 0.8, 0.7), (0.0, 0.0, 0.0)]:
            cart = c.fractional_to_cartesian(vec)
            back = c.cartesian_to_fractional(tuple(cart))
            np.testing.assert_allclose(back, np.array(vec), rtol=1e-10, atol=1e-10)

    def test_supercell_counts_and_volume(self):
        base = create_diamond_structure()
        v0 = base.get_volume()
        sc = base.generate_supercell(2, 3, 1)
        self.assertEqual(len(sc.atoms), len(base.atoms) * 2 * 3 * 1)
        self.assertAlmostEqual(sc.get_volume(), v0 * 2 * 3 * 1, places=6)

    def test_save_load_roundtrip(self):
        c = Crystal("save_load")
        c.set_lattice_parameters(3.0, 4.0, 5.0, 90.0, 100.0, 70.0)
        c.add_atom("C", 0.1, 0.2, 0.3, occupancy=0.8)
        c.add_atom("O", 0.4, 0.5, 0.6)
        fd, path = tempfile.mkstemp(suffix=".json")
        os.close(fd)
        try:
            c.save_to_file(path)
            loaded = Crystal.load_from_file(path)
            self.assertEqual(loaded.name, c.name)
            self.assertEqual(len(loaded.atoms), 2)
            self.assertEqual(loaded.lattice_parameters["a"], 3.0)
        finally:
            os.remove(path)

    def test_load_old_format_atoms(self):
        # old format atoms: [el, x, y, z] without occupancy
        data = {
            "name": "old_format",
            "lattice_parameters": {"a": 3, "b": 3, "c": 3, "alpha": 90, "beta": 90, "gamma": 90},
            "space_group": "P1",
            "atoms": [["C", 0.0, 0.0, 0.0], ["O", 0.5, 0.5, 0.5]],
        }
        fd, path = tempfile.mkstemp(suffix=".json")
        os.close(fd)
        try:
            with open(path, "w", encoding="utf-8") as f:
                json.dump(data, f)
            loaded = Crystal.load_from_file(path)
            self.assertEqual(len(loaded.atoms), 2)
            # occupancy defaults to 1.0
            self.assertTrue(all(getattr(a, "occupancy", 1.0) == 1.0 for a in loaded.atoms))
        finally:
            os.remove(path)

    def test_normalize_fractional(self):
        c = Crystal("norm")
        c.set_lattice_parameters(3, 3, 3)
        c.add_atom("C", -0.1, 1.2, 2.3)
        c.normalize_fractional()
        for a in c.atoms:
            self.assertTrue(0.0 <= a.x < 1.0)
            self.assertTrue(0.0 <= a.y < 1.0)
            self.assertTrue(0.0 <= a.z < 1.0)

    def test_gamma_near_singular_raises(self):
        c = Crystal("singular")
        # gamma extremely small but > 0 to pass input validation; should fail when using lattice vectors
        # Implementation raises when abs(sin(gamma)) < 1e-12; set gamma << 1e-12 (in degrees)
        c.set_lattice_parameters(3, 3, 3, 90, 90, 1e-13)
        with self.assertRaises(ValueError):
            _ = c.get_volume()

    def test_missing_mass_error(self):
        c = Crystal("missing_mass")
        c.set_lattice_parameters(3, 3, 3)
        c.add_atom("Xx", 0.0, 0.0, 0.0)  # element not in default table
        with self.assertRaises(ValueError):
            _ = c.get_density()


if __name__ == "__main__":
    unittest.main()
