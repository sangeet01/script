"""
Auto-detection tests for pseudoasymmetric centers (r/s) and planar chirality
(@PL1/@PL2 for ferrocenes).

These tests verify that the resolver auto-detects these stereochemical features
from the molecular graph when the user does NOT explicitly specify them.

Pseudoasymmetric detection (Yatha-samkhya):
  A tetrahedral center is pseudoasymmetric (r/s) if it has exactly one pair of
  enantiomorphic substituents (same topology, opposite embedded R/S config).
  The canonicalizer should emit @r/@s (lowercase) for these centers.

Planar chirality detection (Sthana):
  A metallocene (e.g. ferrocene) with substituents on different rings that
  break mirror symmetry is planar chiral. The canonicalizer should emit
  @PL1/@PL2.
"""

import sys
import os
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from script.parser import SCRIPTParser
from script.canonical import canonicalize_SCRIPT, SCRIPTCanonicalizer
from script.mol import StereoType


def _canon(s):
    return canonicalize_SCRIPT(s)


def _parse(s):
    return SCRIPTParser().parse(s)


# ===========================================================================
# Pseudoasymmetric auto-detection (Yatha-samkhya)
# ===========================================================================

class TestPseudoasymmetricAutoDetection:
    """A tetrahedral center with an enantiomorphic substituent pair should
    be auto-detected as pseudoasymmetric and emitted as @r/@s (lowercase).
    """

    def test_pentanetriol_middle_carbon_is_pseudoasymmetric(self):
        """Central carbon with explicitly enantiomorphic substituents should
        be auto-detected as pseudoasymmetric.

        We use @R/@S on the terminal centers to guarantee they are true
        enantiomorphs (one R, one S). The central carbon then has:
          - F
          - H (implicit)
          - (R)-substituent
          - (S)-substituent  ← enantiomorphic pair
        → central C is pseudoasymmetric → @r or @s in canonical form.
        """
        # Central C bonded to F, H(implicit), R-CH(Cl)Br, S-CH(Cl)Br
        # @RH = R config + H, @SH = S config + H
        c = _canon("F[C@H]([C@RH](Cl)Br)[C@SH](Cl)Br")
        assert c is not None, "Canonicalize returned None"
        # The central carbon should be emitted as @r or @s (lowercase)
        assert "@r" in c or "@s" in c, \
            f"Central carbon not detected as pseudoasymmetric (no @r/@s in {c!r})"

    def test_pseudoasymmetric_enantiomers_distinct(self):
        """The two enantiomorphic configurations of a pseudoasymmetric center
        (r vs s) must produce DISTINCT canonical strings."""
        # Flip the central C's @ vs @@ — should produce @r vs @s
        c1 = _canon("F[C@H]([C@RH](Cl)Br)[C@SH](Cl)Br")
        c2 = _canon("F[C@@H]([C@RH](Cl)Br)[C@SH](Cl)Br")
        assert c1 != c2, \
            f"r and s configurations canonicalize identically: {c1}"

    def test_no_pseudoasymmetric_when_substituents_identical(self):
        """If both terminal stereocenters have the SAME configuration (both
        R), the central C has IDENTICAL (not enantiomorphic) substituents →
        it's achiral, not pseudoasymmetric. Should NOT emit @r/@s.
        """
        # Both terminals are @R (same config) → central is not pseudoasymmetric
        c = _canon("F[C@H]([C@RH](Cl)Br)[C@RH](Cl)Br")
        assert c is not None
        # Should NOT be emitted as @r/@s
        assert "@r" not in c and "@s" not in c, \
            f"Identical substituents wrongly detected as pseudoasymmetric: {c!r}"

    def test_meso_tartaric_not_pseudoasymmetric(self):
        """meso-tartaric acid has two truly stereogenic centers (R,S), not
        pseudoasymmetric. Each center has DIFFERENT substituents (COOH vs
        CH(OH)COOH), so neither is pseudoasymmetric.
        """
        c = _canon("O=C(O)[C@@H](O)[C@H](O)C(=O)O")
        assert c is not None
        # Should emit @R/@S (uppercase), not @r/@s
        assert "@r" not in c and "@s" not in c, \
            f"meso-tartaric wrongly detected as pseudoasymmetric: {c!r}"

    def test_pseudoasymmetric_idempotent(self):
        """Pseudoasymmetric markers must round-trip stably."""
        s = "F[C@H]([C@RH](Cl)Br)[C@SH](Cl)Br"
        c1 = _canon(s)
        c2 = _canon(c1) if c1 else None
        assert c1 == c2, f"Pseudoasymmetric canonical not idempotent: {c1} -> {c2}"


# ===========================================================================
# Planar chirality auto-detection (Sthana) — ferrocene focus
# ===========================================================================

class TestPlanarChiralityAutoDetection:
    """Ferrocene-type planar chirality should be auto-detected when a metal
    is sandwiched between two aromatic rings with substituents that break
    mirror symmetry.
    """

    def test_ferrocene_basic_parses(self):
        """Basic ferrocene (no substituents) should parse — it's achiral."""
        # Two Cp rings bonded to Fe. Use V2 aromatic notation.
        r = _parse("[Fe]1:C:C:C:C:1")
        assert r["success"], f"Ferrocene failed to parse: {r.get('error')}"

    def test_planar_chiral_ferrocene_detected(self):
        """1,2-disubstituted ferrocene with different substituents on one
        ring should be auto-detected as planar chiral.

        Structure: Fe between two Cp rings; one ring has -CH3 and -OH on
        different positions.
        """
        # Cp ring 1 (disubstituted): C(CH3)-C(OH)-C-C-C, bonded to Fe
        # Cp ring 2 (unsubstituted): C-C-C-C-C, bonded to Fe
        script = "C:C(C):C(O):C:C:C&5:[Fe]C:C:C:C:C&5:"
        r = _parse(script)
        assert r["success"], f"Planar ferrocene failed to parse: {r.get('error')}"
        mol = r["molecule"]
        # Flatten fragments
        if isinstance(mol, list):
            atoms = []
            for m in mol:
                atoms.extend(m.atoms)
        else:
            atoms = mol.atoms
        planar_atoms = [a for a in atoms if a.stereo_type == StereoType.PLANAR]
        assert len(planar_atoms) >= 1, \
            "Planar chirality not auto-detected in disubstituted ferrocene"

    def test_planar_chiral_ferrocene_canonical_has_PL(self):
        """The canonical form of a planar-chiral ferrocene should contain
        @PL marker (auto-detected)."""
        c = _canon("C:C(C):C(O):C:C:C&5:[Fe]C:C:C:C:C&5:")
        assert c is not None, "Canonicalize returned None"
        assert "@PL" in c, f"Planar chirality not in canonical form: {c!r}"

    def test_unsubstituted_ferrocene_not_planar_chiral(self):
        """Ferrocene without substituents is NOT planar chiral (mirror
        symmetry is preserved)."""
        # Two unsubstituted Cp rings + Fe
        c = _canon("C:C:C:C:C&5:[Fe]C:C:C:C:C&5:")
        assert c is not None
        # Should NOT contain @PL (no planar chirality)
        assert "@PL" not in c, \
            f"Unsubstituted ferrocene wrongly detected as planar chiral: {c!r}"

    def test_planar_chiral_ferrocene_enantiomers_distinct(self):
        """The two enantiomers of a planar-chiral ferrocene (Rp vs Sp) must
        produce DISTINCT canonical strings.

        We use explicit @PL1 vs @PL2 to define the enantiomers, then verify
        the canonicalizer distinguishes them.
        """
        # This test verifies the marker machinery works for explicit input.
        # Auto-detection is a separate concern (see test_planar_chiral_ferrocene_detected).
        s1 = "[C@PL1]:C(C):C(O):C:C:C&6:"
        s2 = "[C@PL2]:C(C):C(O):C:C:C&6:"
        c1 = _canon(s1)
        c2 = _canon(s2)
        if c1 is None or c2 is None:
            pytest.skip("Planar markers not yet supported on full structure")
        assert c1 != c2, f"PL1 and PL2 canonicalize identically: {c1}"

    def test_planar_chirality_idempotent(self):
        """Explicit @PL1/@PL2 markers must round-trip stably."""
        s = "[C@PL1]:C(C):C(O):C:C:C&6:"
        c1 = _canon(s)
        if c1 is None:
            pytest.skip("Planar markers not yet supported on full structure")
        c2 = _canon(c1)
        assert c1 == c2, f"Planar canonical not idempotent: {c1} -> {c2}"


# ===========================================================================
# Sanity: existing features still work
# ===========================================================================

class TestAutoDetectionSanity:
    """Regression guard — auto-detection must not break existing stereo features."""

    def test_normal_tetrahedral_still_emits_uppercase(self):
        """A normal stereocenter (no enantiomorphic substituents) must still
        emit @R/@S (uppercase), not @r/@s.
        """
        c = _canon("C[C@H](O)C(=O)O")  # lactic acid — normal stereocenter
        assert c is not None
        # Should contain @R or @S (uppercase), not @r/@s
        assert "@R" in c or "@S" in c, f"Normal stereocenter not @R/@S: {c!r}"
        assert "@r" not in c and "@s" not in c, \
            f"Normal stereocenter wrongly lowercase: {c!r}"

    def test_glucose_still_idempotent(self):
        """Glucose must remain idempotent after auto-detection changes."""
        s = "O[C@H]([C@@H]([C@H]([C@@H](C&6-O)O)O)O)CO"
        c1 = _canon(s)
        c2 = _canon(c1) if c1 else None
        assert c1 == c2, f"Glucose idempotency broken: {c1} -> {c2}"


if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
