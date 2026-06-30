"""
tests/test_periodic.py
─────────────────────
V3.6 periodic topology tests — LQG / crystal / MOF support.

Validates:
  1. Periodic bond @tx,ty,tz parse round-trip
  2. Lattice parameter encoding/decoding in [[context;a,b,c,α,β,γ]]
  3. Canonical ranking distinguishes atoms with different translation environments
  4. Serialization: canonical string contains @tx,ty,tz on cross-cell bonds
  5. Real-world-style cases: NaCl (rock salt), simple MOF pillar, rutile TiO2
"""

import sys, os, math
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import pytest
from script.parser import SCRIPTParser
from script.canonical import SCRIPTCanonicalizer

parser    = SCRIPTParser()
canonicalizer = SCRIPTCanonicalizer()

# ─── helpers ──────────────────────────────────────────────────────────────────

def parse_ok(s):
    r = parser.parse(s)
    assert r is not None and r.get('success'), \
        f"Parse failed for {s!r}: {r.get('error') if r else 'None returned'}"
    return r['molecule'] if r else None

def parse_fail(s):
    r = parser.parse(s)
    assert r is None or not r.get('success'), \
        f"Expected parse failure for {s!r} but got success"


# ─── 1. Periodic bond parsing ─────────────────────────────────────────────────

class TestPeriodicBondParsing:

    def test_single_bond_pos_translation(self):
        """Fe-@1,0,0Fe: single bond with +a lattice shift."""
        mol = parse_ok("Fe-@1,0,0Fe")
        periodic_bonds = [b for b in mol.bonds
                          if getattr(b, 'translation', (0,0,0)) != (0,0,0)]
        assert len(periodic_bonds) == 1
        assert periodic_bonds[0].translation == (1, 0, 0)

    def test_single_bond_neg_translation(self):
        """Fe-@-1,0,0Fe: single bond with -a lattice shift."""
        mol = parse_ok("Fe-@-1,0,0Fe")
        periodic_bonds = [b for b in mol.bonds
                          if getattr(b, 'translation', (0,0,0)) != (0,0,0)]
        assert len(periodic_bonds) == 1
        assert periodic_bonds[0].translation == (-1, 0, 0)

    def test_dative_bond_c_translation(self):
        """Ti->@0,0,1O: dative bond with +c shift (MOF pillar)."""
        mol = parse_ok("Ti->@0,0,1O")
        assert mol.is_periodic
        pb = [b for b in mol.bonds if getattr(b, 'translation', (0,0,0)) != (0,0,0)]
        assert len(pb) == 1
        assert pb[0].translation == (0, 0, 1)

    def test_intracell_bond_zero_translation(self):
        """C-@0,0,0C: explicit zero translation — same as plain C-C."""
        mol = parse_ok("C-@0,0,0C")
        # is_periodic should be False (zero translation doesn't set the flag)
        assert not mol.is_periodic

    def test_multiple_periodic_bonds(self):
        """Ti->@0,0,1O->@1,0,0Ti: two cross-cell bonds in one chain."""
        mol = parse_ok("Ti->@0,0,1O->@1,0,0Ti")
        pb = [b for b in mol.bonds if getattr(b, 'translation', (0,0,0)) != (0,0,0)]
        assert len(pb) == 2
        translations = {b.translation for b in pb}
        assert (0, 0, 1) in translations
        assert (1, 0, 0) in translations

    def test_is_periodic_flag(self):
        """mol.is_periodic must be True when any bond has non-zero translation."""
        mol = parse_ok("[Na+]-@1,0,0[Cl-]")
        assert mol.is_periodic

    def test_mixed_periodic_and_intracell(self):
        """C=C-@0,1,0N: one normal bond + one periodic bond."""
        mol = parse_ok("C=C-@0,1,0N")
        normal = [b for b in mol.bonds if getattr(b, 'translation', (0,0,0)) == (0,0,0)]
        periodic = [b for b in mol.bonds if getattr(b, 'translation', (0,0,0)) != (0,0,0)]
        assert len(normal) == 1
        assert len(periodic) == 1


# ─── 2. Lattice parameter encoding ───────────────────────────────────────────

class TestLatticeParameters:

    def test_cubic_lattice(self):
        """[[NaCl;5.6402,5.6402,5.6402,90.0,90.0,90.0]] Na-@1,0,0Cl"""
        s = "[[NaCl;5.6402,5.6402,5.6402,90.0,90.0,90.0]] Na-@1,0,0Cl"
        mol = parse_ok(s)
        assert mol.macroscopic_context == "NaCl"
        assert mol.lattice_vectors is not None
        # Check a ≈ 5.6402
        a = math.sqrt(sum(x**2 for x in mol.lattice_vectors[0]))
        assert abs(a - 5.6402) < 0.01

    def test_tetragonal_rutile(self):
        """[[Rutile;4.593,4.593,2.959,90.0,90.0,90.0]]"""
        s = "[[Rutile;4.593,4.593,2.959,90.0,90.0,90.0]] Ti->@0,0,1O"
        mol = parse_ok(s)
        assert mol.macroscopic_context == "Rutile"
        assert mol.lattice_vectors is not None
        a = math.sqrt(sum(x**2 for x in mol.lattice_vectors[0]))
        c = math.sqrt(sum(x**2 for x in mol.lattice_vectors[2]))
        assert abs(a - 4.593) < 0.01
        assert abs(c - 2.959) < 0.01

    def test_context_without_lattice(self):
        """[[BCC]] Fe — context label only, no lattice params."""
        mol = parse_ok("[[BCC]] Fe")
        assert mol.macroscopic_context == "BCC"
        assert mol.lattice_vectors is None

    def test_lattice_roundtrip(self):
        """Lattice params encoded in context must survive canon round-trip."""
        s = "[[Rutile;4.5930,4.5930,2.9590,90.00,90.00,90.00]] Ti->@0,0,1O"
        mol = parse_ok(s)
        canon = canonicalizer.canonicalize_core(mol)
        assert canon is not None
        assert "[[Rutile;" in canon
        # Re-parse the canonical string
        mol2 = parse_ok(canon)
        assert mol2.macroscopic_context == "Rutile"
        assert mol2.lattice_vectors is not None


# ─── 3. Canonical ranking with translations ───────────────────────────────────

class TestPeriodicCanonicalRanking:

    def test_different_translations_different_strings(self):
        """
        The same atom pair with different translation directions must
        produce different canonical strings — they represent different
        crystal topologies.
        """
        mol_a = parse_ok("Fe-@1,0,0Fe")
        mol_b = parse_ok("Fe-@0,1,0Fe")
        s_a = canonicalizer.canonicalize_core(mol_a)
        s_b = canonicalizer.canonicalize_core(mol_b)
        assert s_a != s_b, (
            "Fe-@1,0,0Fe and Fe-@0,1,0Fe are different nets "
            "and must have different canonical strings"
        )

    def test_same_topology_same_string(self):
        """
        Two strings representing the same periodic topology must canonicalize
        to the same string regardless of atom traversal order.
        """
        # NaCl: Na bonded to Cl with +a translation
        s1 = "Na-@1,0,0Cl"
        s2 = "Cl-@-1,0,0Na"   # same bond, opposite direction
        mol1 = parse_ok(s1)
        mol2 = parse_ok(s2)
        c1 = canonicalizer.canonicalize_core(mol1)
        c2 = canonicalizer.canonicalize_core(mol2)
        assert c1 == c2, (
            f"Same topology should give same canonical string.\n"
            f"  {s1!r} → {c1!r}\n"
            f"  {s2!r} → {c2!r}"
        )

    def test_periodic_string_contains_translation(self):
        """Canonical string of a periodic mol must contain @tx,ty,tz."""
        mol = parse_ok("Ti->@0,0,1O")
        canon = canonicalizer.canonicalize_core(mol)
        assert canon is not None
        assert "@" in canon, f"Expected @ in canonical string, got: {canon!r}"

    def test_no_sentinel_prefix(self):
        """~P~ sentinel must NOT appear in canonical output (it was removed)."""
        mol = parse_ok("Fe-@1,0,0Fe")
        canon = canonicalizer.canonicalize_core(mol)
        assert canon is not None
        assert "~P~" not in canon, f"Unexpected ~P~ sentinel in: {canon!r}"


# ─── 4. Real-world crystal cases ─────────────────────────────────────────────

class TestRealWorldCrystals:

    def test_nacl_rock_salt(self):
        """NaCl: Na⁺ and Cl⁻ connected along all three lattice vectors."""
        s = "[[NaCl;5.6402,5.6402,5.6402,90.0,90.0,90.0]] [Na+]-@1,0,0[Cl-]-@0,1,0[Na+]-@0,0,1[Cl-]"
        mol = parse_ok(s)
        assert mol.is_periodic
        pb = [b for b in mol.bonds if getattr(b, 'translation', (0,0,0)) != (0,0,0)]
        assert len(pb) >= 2

    def test_rutile_tio2(self):
        """TiO2 rutile: Ti in octahedral coordination with O along c-axis."""
        s = "[[Rutile;4.5930,4.5930,2.9590,90.00,90.00,90.00]] Ti(->@0,0,1O)(->@0,0,-1O)"
        mol = parse_ok(s)
        assert mol.macroscopic_context == "Rutile"
        assert mol.is_periodic

    def test_mof_pillar_bond(self):
        """Simplified MOF: Zn node connected to oxygen pillar across cell."""
        s = "[[MOF-5;25.8320,25.8320,25.8320,90.0,90.0,90.0]] Zn->@1,0,0O"
        mol = parse_ok(s)
        assert mol.macroscopic_context == "MOF-5"
        assert mol.lattice_vectors is not None
        assert mol.is_periodic

    def test_periodic_plus_organic(self):
        """Periodic bond followed by organic fragment — mixed system."""
        s = "Fe-@0,0,1Fe.CC(=O)O"
        r = parser.parse(s)
        assert r is not None and r.get('success'), f"Parse failed: {r.get('error') if r else 'None'}"
        # Multi-component: result may be list or single mol with fragments
        # The periodic fragment must have is_periodic True
        mol = r['molecule']
        # For multi-component, check the first periodic fragment
        if isinstance(mol, list):
            periodic_mols = [m for m in mol if getattr(m, 'is_periodic', False)]
            assert len(periodic_mols) >= 1
        else:
            assert mol.is_periodic


# ─── 5. Serialization round-trips ─────────────────────────────────────────────

class TestPeriodicSerializationRoundtrip:

    def test_translation_survives_roundtrip(self):
        """Parse → canonicalize → reparse: translation vectors must survive."""
        original = "Fe-@1,0,0Fe"
        mol1 = parse_ok(original)
        canon = canonicalizer.canonicalize_core(mol1)
        assert canon is not None
        mol2 = parse_ok(canon)
        pb2 = [b for b in mol2.bonds if getattr(b, 'translation', (0,0,0)) != (0,0,0)]
        assert len(pb2) >= 1

    def test_lattice_and_translation_roundtrip(self):
        """Full periodic SCRIPT with lattice params round-trips cleanly."""
        s = "[[NaCl;5.6402,5.6402,5.6402,90.00,90.00,90.00]] Na-@1,0,0Cl"
        mol1 = parse_ok(s)
        canon = canonicalizer.canonicalize_core(mol1)
        assert canon is not None
        mol2 = parse_ok(canon)
        assert mol2.macroscopic_context == "NaCl"
        assert mol2.is_periodic

    def test_canonical_is_idempotent(self):
        """canonicalize(canonicalize(s)) == canonicalize(s) for periodic strings."""
        s = "Ti->@0,0,1O"
        mol1  = parse_ok(s)
        c1    = canonicalizer.canonicalize_core(mol1)
        mol2  = parse_ok(c1)
        c2    = canonicalizer.canonicalize_core(mol2)
        assert c1 == c2, f"Canonicalization not idempotent:\n  c1={c1!r}\n  c2={c2!r}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
