"""
SCRIPT Stereochemistry Correctness Test Suite
=============================================

Comprehensive tests for the SCRIPT stereo system, covering all 7 known bugs
plus the existing behaviour that already works.

Designed to drive the fixes: each test class targets one bug from the audit.
Tests are expected to FAIL until the corresponding fix is applied.

Paninian framework mapping (for reference):
  Pratyaya = @/@@ marker (the suffix that signals stereo)
  Vak      = DFS neighbor order (the order atoms are "spoken")
  Sandhi   = CIP priority reconciliation (how two reference frames combine)
  Lopa     = elided/ghost neighbors (implicit H, implicit lone pair, ring-closure partner)
"""

import sys
import os
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from script.parser import SCRIPTParser
from script.canonical import SCRIPTCanonicalizer, canonicalize_SCRIPT
from script.mol import StereoType, BondType


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _canon(s):
    """Convenience: canonicalize a SCRIPT string."""
    return canonicalize_SCRIPT(s)


def _distinct(s1, s2):
    """Two strings must canonicalize to DIFFERENT outputs."""
    c1 = _canon(s1)
    c2 = _canon(s2)
    return c1, c2, c1 != c2


def _same(s1, s2):
    """Two strings must canonicalize to the SAME output."""
    c1 = _canon(s1)
    c2 = _canon(s2)
    return c1, c2, c1 == c2


def _idempotent(s):
    """Canonical form must be stable across round-trips."""
    c1 = _canon(s)
    if c1 is None:
        return False
    c2 = _canon(c1)
    return c1 == c2


# ===========================================================================
# Bug 1: Enantiomers must produce distinct canonical strings
# ===========================================================================

class TestEnantiomerDistinction:
    """Critical: (R,R) and (S,S) enantiomers must NOT canonicalize identically.

    The CIP resolver assigns correct chirality bits, but the canonicalizer's
    parity transform was collapsing enantiomers to the same string. This is
    the most fundamental violation of the "one molecule, one string" promise.
    """

    def test_butanediol_RR_vs_SS_distinct(self):
        c1, c2, distinct = _distinct(
            "C[C@@H](O)[C@@H](O)C",   # (R,R)
            "C[C@H](O)[C@H](O)C",     # (S,S)
        )
        assert distinct, f"(R,R) and (S,S)-2,3-butanediol canonicalize identically: {c1}"

    def test_dichlorobutane_RR_vs_SS_distinct(self):
        c1, c2, distinct = _distinct(
            "C[C@@H](Cl)[C@@H](Cl)C",  # (R,R)
            "C[C@H](Cl)[C@H](Cl)C",    # (S,S)
        )
        assert distinct, f"(R,R) and (S,S)-2,3-dichlorobutane canonicalize identically: {c1}"

    def test_tartaric_RR_vs_SS_distinct(self):
        c1, c2, distinct = _distinct(
            "O=C(O)[C@@H](O)[C@@H](O)C(=O)O",  # (R,R)
            "O=C(O)[C@H](O)[C@H](O)C(=O)O",    # (S,S)
        )
        assert distinct, f"(R,R) and (S,S)-tartaric acid canonicalize identically: {c1}"

    def test_threonine_L_vs_D_distinct(self):
        c1, c2, distinct = _distinct(
            "N[C@@H](C)[C@H](O)C(=O)O",  # L-Threonine
            "N[C@H](C)[C@@H](O)C(=O)O",  # D-Threonine
        )
        assert distinct, f"L- and D-Threonine canonicalize identically: {c1}"

    def test_threonine_L_vs_Lallo_distinct(self):
        # L-Threonine vs L-allo-Threonine (diastereomers, should also be distinct)
        c1, c2, distinct = _distinct(
            "N[C@@H](C)[C@H](O)C(=O)O",      # L-Threonine
            "N[C@@H](C)[C@@H](O)C(=O)O",     # L-allo-Threonine
        )
        assert distinct, f"L-Threonine and L-allo-Threonine canonicalize identically: {c1}"

    def test_meso_distinct_from_both_enantiomers(self):
        # meso-tartaric must be distinct from both (R,R) and (S,S)
        meso = "O=C(O)[C@@H](O)[C@H](O)C(=O)O"
        rr   = "O=C(O)[C@@H](O)[C@@H](O)C(=O)O"
        ss   = "O=C(O)[C@H](O)[C@H](O)C(=O)O"
        c_meso, c_rr, distinct1 = _distinct(meso, rr)
        c_meso2, c_ss, distinct2 = _distinct(meso, ss)
        assert distinct1, f"meso and (R,R) canonicalize identically: {c_meso}"
        assert distinct2, f"meso and (S,S) canonicalize identically: {c_meso}"


# ===========================================================================
# Bug 2: In-ring chiral centers must preserve @/@@
# ===========================================================================

class TestInRingStereo:
    """Chiral atoms that are part of a ring (proline, menthol, etc.) must
    preserve their @/@@ marker through canonicalization.

    Root cause: ring-closure bonds (added via _handle_v2_ring_closure / add_ring)
    were not updating the chiral center's _initial_nbrs list, so the resolver
    saw only 2 of 4 neighbors and silently dropped the chirality.
    """

    def test_proline_L_vs_D_distinct(self):
        c1, c2, distinct = _distinct(
            "O=C(O)[C@@H]1CCCN1",  # L-Proline
            "O=C(O)[C@H]1CCCN1",   # D-Proline
        )
        assert distinct, f"L- and D-Proline canonicalize identically: {c1}"

    def test_proline_canonical_contains_chiral_marker(self):
        c = _canon("O=C(O)[C@@H]1CCCN1")
        assert c is not None
        assert "@" in c, f"L-Proline canonical lost chiral marker: {c}"

    def test_menthol_canonical_contains_chiral_marker(self):
        c = _canon("CC(C)[C@@H]1CC[C@@H](C(C)C)[C@H](O)C1")
        assert c is not None
        # Menthol has 3 stereocenters; at least one @ must survive
        assert "@" in c, f"L-Menthol canonical lost all chiral markers: {c}"


# ===========================================================================
# Bug 3: E/Z alkene canonicalization must be idempotent
# ===========================================================================

class TestAlkeneEZStereo:
    """Cis/trans alkene stereo must round-trip stably. Previously the
    canonicalizer flipped / and \\ on every pass, never converging.
    """

    def test_2butene_E_idempotent(self):
        assert _idempotent("C/C=C/C"), f"2-butene E not idempotent"

    def test_2butene_Z_idempotent(self):
        assert _idempotent("C/C=C\\C"), f"2-butene Z not idempotent"

    def test_2butene_E_vs_Z_distinct(self):
        c1, c2, distinct = _distinct("C/C=C/C", "C/C=C\\C")
        assert distinct, f"2-butene E and Z canonicalize identically: {c1}"

    def test_dichloroethene_cis_vs_trans_distinct(self):
        c1, c2, distinct = _distinct("Cl/C=C/Cl", "Cl/C=C\\Cl")
        assert distinct, f"DCE cis and trans canonicalize identically: {c1}"

    def test_dichloroethene_cis_idempotent(self):
        assert _idempotent("Cl/C=C\\Cl"), f"DCE cis not idempotent"

    def test_dichloroethene_trans_idempotent(self):
        assert _idempotent("Cl/C=C/Cl"), f"DCE trans not idempotent"


# ===========================================================================
# Bug 4: @AX1 vs @AX2 allene distinction must be preserved
# ===========================================================================

class TestAlleneStereo:
    """Allene axial stereo @AX1 (CW) vs @AX2 (CCW) must produce distinct
    canonical strings. The grammar tokenised them but the parser dispatch
    didn't handle them — they fell through silently.
    """

    def test_allene_AX1_vs_AX2_distinct(self):
        c1, c2, distinct = _distinct(
            "C[C@AX1]=C=CC",
            "C[C@AX2]=C=CC",
        )
        assert distinct, f"@AX1 and @AX2 canonicalize identically: {c1}"

    def test_allene_AX1_canonical_contains_marker(self):
        c = _canon("C[C@AX1]=C=CC")
        assert c is not None
        assert "AX" in c, f"@AX1 marker lost in canonical: {c}"

    def test_allene_AX2_canonical_contains_marker(self):
        c = _canon("C[C@AX2]=C=CC")
        assert c is not None
        assert "AX" in c, f"@AX2 marker lost in canonical: {c}"

    def test_allene_AX1_idempotent(self):
        assert _idempotent("C[C@AX1]=C=CC"), f"@AX1 not idempotent"

    def test_allene_AX2_idempotent(self):
        assert _idempotent("C[C@AX2]=C=CC"), f"@AX2 not idempotent"


# ===========================================================================
# Bug 5: Non-tetrahedral stereo must distinguish geometric isomers
# ===========================================================================

class TestPolyhedralStereo:
    """@SP/@OH/@TB/@PY must distinguish cis vs trans, fac vs mer, etc.

    We extend the Paninian pratyaya (marker) scheme: @SP1/@SP2 for cis/trans
    square planar; @OH1-5 for octahedral; @TB1/@TB2 for trigonal bipyramidal
    apical/equatorial; @PY1/@PY2 for pyramidal enantiomers.
    """

    def test_square_planar_cis_vs_trans_distinct(self):
        c1, c2, distinct = _distinct(
            "[Pt@SP1](Cl)([NH3])(Cl)([NH3])",   # cis
            "[Pt@SP2](Cl)([NH3])([NH3])(Cl)",   # trans
        )
        assert distinct, f"cis and trans Pt@SP canonicalize identically: {c1}"

    def test_octahedral_fac_vs_mer_distinct(self):
        # fac = 3 identical ligands on one face; mer = 3 in a meridian
        c1, c2, distinct = _distinct(
            "[Co@OH1](Cl)(Cl)(Cl)(N)(N)N",   # fac
            "[Co@OH2](Cl)(Cl)(Cl)(N)(N)N",   # mer
        )
        assert distinct, f"fac and mer Co@OH canonicalize identically: {c1}"


# ===========================================================================
# Bug 6: Sulfoxide and other 3-coordinate chiral centers must preserve stereo
# ===========================================================================

class TestLonePairStereo:
    """Sulfoxide (R-S(=O)-R'), chiral phosphine oxide, and other 3-coordinate
    chiral centers must preserve their @/@@ marker. The lone pair counts as
    a 4th "ghost" neighbor (Lopa in Paninian terms).
    """

    def test_sulfoxide_R_vs_S_distinct(self):
        c1, c2, distinct = _distinct(
            "C[S@](=O)C",
            "C[S@@](=O)C",
        )
        assert distinct, f"Sulfoxide R and S canonicalize identically: {c1}"

    def test_sulfoxide_canonical_contains_marker(self):
        c = _canon("C[S@](=O)C")
        assert c is not None
        assert "@" in c, f"Sulfoxide canonical lost chiral marker: {c}"

    def test_sulfoxide_idempotent(self):
        assert _idempotent("C[S@](=O)C"), f"Sulfoxide not idempotent"


# ===========================================================================
# Bug 7: Quaternary ammonium — explicit + required (documented behaviour)
# ===========================================================================

class TestQuaternaryNitrogen:
    """Bare 4-coordinate N is rejected by Sandhi (correct chemistry —
    a 4-coordinate N really IS a cation). The fix is documentation, not
    a code change: SCRIPT requires explicit + on quaternary N.
    """

    def test_bare_quaternary_N_rejected(self):
        # This is the EXPECTED behaviour — Sandhi correctly rejects over-valent bare N
        parser = SCRIPTParser()
        r = parser.parse("[N@@](C)(C)(C)C")
        assert not r["success"], "Bare 4-coordinate N should be rejected by Sandhi"

    def test_explicit_plus_quaternary_N_accepted(self):
        parser = SCRIPTParser()
        r = parser.parse("[N@@+](C)(C)(C)C")
        assert r["success"], f"Explicit + on quaternary N should be accepted: {r.get('error')}"

    def test_quaternary_N_enantiomers_distinct(self):
        c1, c2, distinct = _distinct(
            "[N@@+](C)(C)(C)C",
            "[N@+](C)(C)(C)C",
        )
        assert distinct, f"Quaternary N enantiomers canonicalize identically: {c1}"


# ===========================================================================
# Bug 8: Transition metals (Ir, Hf, Ta, Re, Os, etc.) must parse
# ===========================================================================

class TestTransitionMetalParsing:
    """The ATOM lexer token was missing many transition metals. They are in
    the `element` non-terminal rule but not the ATOM terminal, so they lexed
    as the first letter (e.g. Ir -> I + r).
    """

    @pytest.mark.parametrize("symbol", [
        "Ir", "Hf", "Ta", "Re", "Os", "Tc", "Ru", "Rh", "Nb", "Mo",
        "Yb", "Lu", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
        "Tb", "Dy", "Ho", "Er", "Tm",
        "Th", "Pa", "Np", "Pu", "Am", "Cm",
    ])
    def test_element_parses_in_bracket(self, symbol):
        parser = SCRIPTParser()
        # Build a minimal molecule: [El]
        r = parser.parse(f"[{symbol}]")
        assert r["success"], f"Failed to parse [{symbol}]: {r.get('error')}"
        mol = r["molecule"]
        # Should be 1 atom with the right symbol
        if isinstance(mol, list):
            mol = mol[0]
        assert mol.atoms[0].symbol == symbol, \
            f"[{symbol}] parsed as symbol={mol.atoms[0].symbol!r}"

    def test_iridium_complex_parses(self):
        # The original failing case from the audit
        parser = SCRIPTParser()
        r = parser.parse("[Ir@SP](Cl)(Cl)([PH3])([PH3])")
        assert r["success"], f"Ir complex failed to parse: {r.get('error')}"
        mol = r["molecule"]
        assert mol.atoms[0].symbol == "Ir"
        assert mol.atoms[0].stereo_type == StereoType.SQUARE_PLANAR


# ===========================================================================
# Sanity checks — features that already work and must not regress
# ===========================================================================

class TestStereoSanity:
    """Things that already work. These must continue to pass after the fixes."""

    def test_single_center_R_vs_S_distinct(self):
        c1, c2, distinct = _distinct(
            "[C@H](F)(Cl)Br",
            "[C@@H](F)(Cl)Br",
        )
        assert distinct, f"CHFClBr R and S canonicalize identically: {c1}"

    def test_lactic_acid_R_vs_S_distinct(self):
        c1, c2, distinct = _distinct(
            "C[C@H](O)C(=O)O",
            "C[C@@H](O)C(=O)O",
        )
        assert distinct, f"Lactic acid R and S canonicalize identically: {c1}"

    def test_lactic_acid_two_ways_same_enantiomer(self):
        c1, c2, same = _same(
            "C[C@H](O)C(=O)O",
            "O=C(O)[C@H](C)O",
        )
        assert same, f"Same enantiomer written two ways canonicalize differently: {c1} vs {c2}"

    def test_chiral_phosphine_R_vs_S_distinct(self):
        # Use 4 DIFFERENT substituents — [P@](C)(C)(C)C with 4 identical
        # methyls is achiral (no stereochemistry to distinguish).
        c1, c2, distinct = _distinct(
            "[P@](C)(CC)(CCC)(F)",
            "[P@@](C)(CC)(CCC)(F)",
        )
        assert distinct, f"Chiral phosphine enantiomers canonicalize identically: {c1}"

    def test_chiral_silicon_R_vs_S_distinct(self):
        # 4 different substituents required for chirality
        c1, c2, distinct = _distinct(
            "[Si@](C)(CC)(CCC)(F)",
            "[Si@@](C)(CC)(CCC)(F)",
        )
        assert distinct, f"Chiral Si enantiomers canonicalize identically: {c1}"

    def test_glucose_idempotent(self):
        assert _idempotent("O[C@H]([C@@H]([C@H]([C@@H](C&6-O)O)O)O)CO"), \
            "Glucose canonical not idempotent"

    def test_multi_fragment_stereo_preserved(self):
        # Lactic acid + Na+ salt — stereo must survive
        c = _canon("C[C@H](O)C(=O)O.[Na+]")
        assert c is not None
        assert "@" in c, f"Lactic+Na salt lost stereo: {c}"
        assert "." in c, f"Multi-fragment separator lost: {c}"


if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
