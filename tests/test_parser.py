"""
Test suite for SCRIPT parser - comprehensive coverage
Includes: basics, stereochemistry, reactions, query atoms, organometallics,
polymers, salts, and boss-fight molecules (Taxol, Strychnine, Glucose).
"""

import pytest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from script.parser import SCRIPTParser, parse_script
from script.canonical import SCRIPTCanonicalizer
from script.mol import BondType, StereoType


class TestSCRIPTParser:
    """Test cases for SCRIPT parser functionality"""

    def setup_method(self):
        """Setup for each test"""
        self.parser = SCRIPTParser()
        self.canon = SCRIPTCanonicalizer()

    # Tier 0: Simple molecules

    def test_simple_molecules(self):
        """Test parsing of simple molecules"""
        test_cases = ["C", "CC", "CCO", "CCN"]
        for script in test_cases:
            result = self.parser.parse(script)
            assert result["success"], f"Failed to parse {script}: {result['error']}"
            assert len(result["atoms"]) > 0, f"No atoms found for {script}"

    def test_branched_molecules(self):
        """Test parsing of branched structures"""
        test_cases = ["CC(C)C", "CC(O)C", "CC(C)(C)C"]
        for script in test_cases:
            result = self.parser.parse(script)
            assert result["success"], f"Failed to parse {script}: {result['error']}"

    def test_ring_molecules(self):
        """Test parsing of ring structures with topological closures"""
        test_cases = [
            "CCCCCC&6-",
            "C:C:C:C:C:C&6:",
            "CCCC(C&5-)CC&7-",
        ]
        for script in test_cases:
            result = self.parser.parse(script)
            assert result["success"], f"Failed to parse {script}: {result['error']}"

    def test_bracket_atoms(self):
        """Test parsing of bracket notation"""
        test_cases = [
            "[C]", "[C+]", "[C-]", "[13C]", "[C@H]"
        ]
        for script in test_cases:
            result = self.parser.parse(script)
            assert result["success"], f"Failed to parse {script}: {result['error']}"

    def test_multi_component(self):
        """Test parsing of multi-component systems"""
        test_cases = [
            "CCO.CCO",
            "[Na+]~[Cl-]",
        ]
        for script in test_cases:
            result = self.parser.parse(script)
            assert result["success"], f"Failed to parse {script}: {result['error']}"

    # Tier 1: Typed IR features

    def test_dative_bonds(self):
        """Test dative bond parsing and BondType enum"""
        result = self.parser.parse("N->[B](F)(F)F")
        assert result["success"]
        mol = result["molecule"]
        assert mol.bonds[0].bond_type == BondType.DATIVE

    def test_coordinate_bonds(self):
        """Test coordinate/haptic bond parsing"""
        result = self.parser.parse("Fe>1C")
        assert result["success"]
        mol = result["molecule"]
        assert any(b.bond_type == BondType.COORDINATE for b in mol.bonds)

    def test_tautomeric_bonds(self):
        """Test mobile tautomeric bond"""
        for script in ["C=:C-C", "C~:C-C"]:
            result = self.parser.parse(script)
            assert result["success"], f"Failed to parse tautomeric bond: {script}"
            mol = result["molecule"]
            assert any(b.bond_type == BondType.TAUTOMERIC for b in mol.bonds), f"No tautomeric bond found in {script}"

    def test_non_tetrahedral_stereo(self):
        """Test extended stereochemistry types"""
        cases = [
            ("[Pt@SP](Cl)(Cl)([NH3])[NH3]", StereoType.SQUARE_PLANAR),
            ("[Co@OH](Cl)(Cl)(Cl)(N)(N)N", StereoType.OCTAHEDRAL),
        ]
        for script, expected in cases:
            result = self.parser.parse(script)
            assert result["success"], f"Failed: {script}"
            atom = result["molecule"].atoms[0]
            assert atom.stereo_type == expected, f"Wrong stereo type for {script}"

    # Tier 2: Semantic metadata

    def test_reactions(self):
        """Test reaction parsing"""
        result = self.parser.parse("CC>>CCO")
        assert result["success"]
        rxn = result["reaction"]
        assert rxn is not None
        assert len(rxn.reactants) == 1
        assert len(rxn.products) == 1

    def test_reactions_with_agents(self):
        """Test reaction with catalyst/solvent"""
        result = self.parser.parse("CC>[Pd]>CCO")
        assert result["success"]
        rxn = result["reaction"]
        assert rxn is not None
        assert len(rxn.agents) == 1

    def test_atom_mapping(self):
        """Test atom mapping across reaction"""
        result = self.parser.parse("[C:1]=O>>[C:1]O")
        assert result["success"]

    def test_salts_and_solvates(self):
        """Test fragment separator semantics"""
        result_dot = self.parser.parse("CC.O")
        assert result_dot["success"]
        result_tilde = self.parser.parse("[Na+]~[Cl-]")
        assert result_tilde["success"]

    def test_polymers(self):
        """Test polymer notation"""
        cases = [
            "{[CC]}n",
            "{[CC]}<n:50-100>",
            "{[CC]}<n:10>",
        ]
        for script in cases:
            result = self.parser.parse(script)
            assert result["success"], f"Failed: {script}"

    def test_materials_alloys(self):
        """Test alloy occupancy notation"""
        result = self.parser.parse("Ti<~0.9>N<~0.1>")
        assert result["success"]

    def test_materials_crystal(self):
        """Test crystallographic context"""
        result = self.parser.parse("[[Rutile]] Ti(O)2")
        assert result["success"]

    def test_materials_surface(self):
        """Test surface phase boundary"""
        result = self.parser.parse("[[Pt_111]] | >C=O")
        assert result["success"]

    def test_quantum_states(self):
        """Test spin and excited states"""
        result = self.parser.parse("O=O<s:3>")
        assert result["success"]

    # Tier 3: Query atoms, biopolymers, allenes

    def test_query_atoms(self):
        """Test SMARTS-style query atoms"""
        cases = [
            "[#6]CC",
            "[R]CC",
            "[r5]CC",
            "[v3]CC",
            "[!N]CC",
            "[a]CC",
            "[#6,#7]C",
        ]
        for script in cases:
            result = self.parser.parse(script)
            assert result["success"], f"Failed: {script}"
            atom = result["molecule"].atoms[0]
            assert atom.is_query, f"Not marked as query: {script}"
            if script.startswith("[#6]"):
                assert atom.query_atomic_nums == [6]
                assert atom.symbol == "C"

    def test_peptides(self):
        """Test peptide notation parsing"""
        cases = ["{A}", "{A.G.S}", "{pS.G.acK.V}"]
        for script in cases:
            result = self.parser.parse(script)
            assert result["success"], f"Failed: {script}"

    def test_nucleic_acids(self):
        """Test nucleic acid notation"""
        cases = [
            "{dA.dG.dC.dT}",
            "{rA.rG.rC.rU}",
            "{m5C.dG.dA.dT}",
        ]
        for script in cases:
            result = self.parser.parse(script)
            assert result["success"], f"Failed: {script}"

    # Boss Fights: Canonicalization stress tests

    def test_glucose(self):
        """Glucose: 4 stereocenters, pyranose ring"""
        script = "O[C@H]([C@@H]([C@H]([C@@H](C&6-O)O)O)O)CO"
        result = self.parser.parse(script)
        assert result["success"]
        mol = result["molecule"]
        chiral = [a for a in mol.atoms if a.chirality != 0]
        assert len(chiral) >= 4, f"Expected 4 stereocenters, got {len(chiral)}"
        out = self.canon.canonicalize_core(mol)
        assert out is not None
        result2 = self.parser.parse(out)
        assert result2["success"]
        out2 = self.canon.canonicalize_core(result2["molecule"])
        assert out2 == out, f"Glucose canonical roundtrip failed: {out} != {out2}"

    def test_aspirin(self):
        """Aspirin: fused ring, ester, carboxylic acid"""
        scripts = [
            "C(C(=O)O):C(OC(=O)C):C:C:C:C&6:",
            "C(OC(=O)C):C(C(=O)O):C:C:C:C&6:"
        ]
        canonical_results = []
        for script in scripts:
            result = self.parser.parse(script)
            assert result["success"]
            mol = result["molecule"]
            out = self.canon.canonicalize_core(mol)
            assert out is not None
            canonical_results.append(out)
        assert canonical_results[0] == canonical_results[1], \
            f"Aspirin variants produced different canonical strings: {canonical_results}"

    def test_strychnine(self):
        """Strychnine: 7 fused rings, 6 stereocenters"""
        script = "C1CC2C3C4CC5C6C(C34)C2C1C56"
        result = self.parser.parse(script)
        assert result is not None

    def test_taxol(self):
        """Taxol (Paclitaxel): 11 stereocenters, fused bridged system"""
        script = "CC(=O)O[C@H]1C[C@H]2C(C)(C)[C@@H](O)C[C@H]2C1(C)C"
        result = self.parser.parse(script)
        assert result["success"]
        mol = result["molecule"]
        out = self.canon.canonicalize_core(mol)
        assert out is not None

    # Validation and Sandhi

    def test_invalid_syntax(self):
        """Test that invalid syntax is rejected"""
        invalid_cases = ["C(", "C[", "", "C(C)(C)(C)(C)(C)"]
        for script in invalid_cases:
            result = self.parser.parse(script)
            assert not result["success"], f"Should have failed: {script}"

    def test_sandhi_valence_guard(self):
        """Test Paninian Sandhi valence rejection"""
        result = self.parser.parse("C(C)(C)(C)(C)(C)")
        assert not result["success"], "6-valent carbon should be rejected"

    def test_canonicalization_roundtrip(self):
        """Test that canonicalization produces stable output"""
        scripts = ["CCO", "C(C)C", "C:C:C:C:C:C&6:"]
        for s in scripts:
            r1 = self.parser.parse(s)
            assert r1["success"]
            c1 = self.canon.canonicalize_core(r1["molecule"])
            r2 = self.parser.parse(c1)
            assert r2["success"]
            c2 = self.canon.canonicalize_core(r2["molecule"])
            assert c1 == c2, f"Canonical unstable: {s} -> {c1} -> {c2}"


def test_convenience_function():
    """Test the convenience parse_script function"""
    result = parse_script("CCO")
    assert result["success"]
    assert len(result["atoms"]) == 3


class TestMultiFragmentCanonicalization:
    """Regression tests for multi-fragment (salt/solvate/ionic pair) handling.

    The parser intentionally returns ``list[CoreMolecule]`` when the input
    SCRIPT contains ``.`` (solvate) or ``~`` (ionic pair) separators.
    ``canonicalize_core`` only accepts a single CoreMolecule and would crash
    on the list shape; ``canonicalize_mols`` was added to handle both.
    These tests pin that contract so the bug cannot silently return.
    """

    def setup_method(self):
        self.parser = SCRIPTParser()
        self.canon = SCRIPTCanonicalizer()

    def test_parser_returns_list_for_dot_separator(self):
        """``.`` between components must produce a list of fragments."""
        result = self.parser.parse("CN(C)C(=N)N=C(N)N.[Cl-]")
        assert result["success"]
        mol = result["molecule"]
        assert isinstance(mol, list), "Multi-fragment input must return a list"
        assert len(mol) == 2
        assert len(mol[0].atoms) == 9   # metformin cation
        assert len(mol[1].atoms) == 1   # chloride anion

    def test_parser_returns_list_for_tilde_separator(self):
        """``~`` between components must produce a list with ionic-pair sep."""
        result = self.parser.parse("[Na+]~[Cl-]")
        assert result["success"]
        mol = result["molecule"]
        assert isinstance(mol, list)
        assert len(mol) == 2
        # The second fragment carries the separator that preceded it.
        assert mol[1].fragment_separator == "~"

    def test_canonicalize_mols_handles_single_mol(self):
        """Passing a single CoreMolecule must delegate to canonicalize_core."""
        result = self.parser.parse("CCO")
        assert result["success"]
        mol = result["molecule"]
        # Single-fragment input -> CoreMolecule, not list
        from script.mol import CoreMolecule
        assert isinstance(mol, CoreMolecule)
        out = self.canon.canonicalize_mols(mol)
        assert out is not None
        # And the unwrapped-list path must give the same string.
        out_list = self.canon.canonicalize_mols([mol])
        assert out == out_list

    def test_canonicalize_mols_handles_multi_fragment(self):
        """canonicalize_mols must accept a list and join with separators."""
        result = self.parser.parse("CN(C)C(=N)N=C(N)N.[Cl-]")
        mol = result["molecule"]
        out = self.canon.canonicalize_mols(mol)
        assert out is not None
        assert "." in out, "Multi-fragment canonical must contain `.` separator"
        # Round-trip: re-parse the canonical string, must succeed and stay stable.
        r2 = self.parser.parse(out)
        assert r2["success"]
        out2 = self.canon.canonicalize_mols(r2["molecule"])
        assert out == out2, f"Multi-fragment canonical not stable: {out} -> {out2}"

    def test_canonicalize_script_handles_multi_fragment(self):
        """The top-level canonicalize_SCRIPT helper must handle multi-fragment."""
        from script.canonical import canonicalize_SCRIPT
        out = canonicalize_SCRIPT("CN(C)C(=N)N=C(N)N.[Cl-]")
        assert out is not None
        assert "." in out
        # Idempotent.
        assert canonicalize_SCRIPT(out) == out

    def test_canonicalize_script_preserves_tilde_separator(self):
        """``~`` (ionic pair) separator must survive canonicalization."""
        from script.canonical import canonicalize_SCRIPT
        out = canonicalize_SCRIPT("[Na+]~[Cl-]")
        assert out is not None
        assert "~" in out, "Ionic-pair `~` separator must survive canonicalization"
        assert canonicalize_SCRIPT(out) == out


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
