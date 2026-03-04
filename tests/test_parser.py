"""
Test suite for SCRIPT parser
"""

import pytest
import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from script.parser import SCRIPTParser, parse_script


class TestSCRIPTParser:
    """Test cases for SCRIPT parser functionality"""
    
    def setup_method(self):
        """Setup for each test"""
        self.parser = SCRIPTParser()
    
    def test_simple_molecules(self):
        """Test parsing of simple molecules"""
        test_cases = [
            "C",      # methane
            "CC",     # ethane  
            "CCO",    # ethanol
            "CCN",    # ethylamine
        ]
        
        for script in test_cases:
            result = self.parser.parse(script)
            assert result["success"], f"Failed to parse {script}: {result['error']}"
            assert len(result["atoms"]) > 0, f"No atoms found for {script}"
    
    def test_branched_molecules(self):
        """Test parsing of branched structures"""
        test_cases = [
            "CC(C)C",     # isobutane
            "CC(O)C",     # isopropanol
            "CC(C)(C)C",  # tert-butane
        ]
        
        for script in test_cases:
            result = self.parser.parse(script)
            assert result["success"], f"Failed to parse {script}: {result['error']}"
    
    def test_ring_molecules(self):
        """Test parsing of ring structures"""
        test_cases = [
            "C1CCCCC6",        # cyclohexane (local ring)
            "C=CC=CC=C6",      # benzene (Kekulé)
        ]
        
        for script in test_cases:
            result = self.parser.parse(script)
            assert result["success"], f"Failed to parse {script}: {result['error']}"
    
    def test_bracket_atoms(self):
        """Test parsing of bracket notation"""
        test_cases = [
            "[C]",        # carbon
            "[C+]",       # carbocation
            "[C-]",       # carbanion
            "[13C]",      # carbon-13
            "[C@H]",      # chiral carbon
        ]
        
        for script in test_cases:
            result = self.parser.parse(script)
            assert result["success"], f"Failed to parse {script}: {result['error']}"
    
    def test_multi_component(self):
        """Test parsing of multi-component systems"""
        test_cases = [
            "CCO.CCO",        # ethanol dimer
            "[Na+].[Cl-]",    # sodium chloride
        ]
        
        for script in test_cases:
            result = self.parser.parse(script)
            assert result["success"], f"Failed to parse {script}: {result['error']}"
    
    def test_invalid_syntax(self):
        """Test that invalid syntax is rejected"""
        invalid_cases = [
            "C(",         # unclosed branch
            "C[",         # unclosed bracket
            "C1CC",       # unmatched ring (should fail in validation)
            "",           # empty string
        ]
        
        for script in invalid_cases:
            result = self.parser.parse(script)
            assert not result["success"], f"Should have failed to parse {script}"
    
    def test_peptide_mode(self):
        """Test peptide notation parsing"""
        test_cases = [
            "{A}",           # alanine
            "{A.G.S}",       # Ala-Gly-Ser
            "{C[A].G.C[A]}", # disulfide bridge
        ]
        
        for script in test_cases:
            result = self.parser.parse(script)
            assert result["success"], f"Failed to parse peptide {script}: {result['error']}"
    
    def test_polymer_mode(self):
        """Test polymer notation parsing"""
        test_cases = [
            "{CC}n",      # polyethylene
            "{CC}100",    # 100-mer
        ]
        
        for script in test_cases:
            result = self.parser.parse(script)
            assert result["success"], f"Failed to parse polymer {script}: {result['error']}"


def test_convenience_function():
    """Test the convenience parse_script function"""
    result = parse_script("CCO")
    assert result["success"]
    assert len(result["atoms"]) == 3


if __name__ == "__main__":
    # Run tests directly
    pytest.main([__file__, "-v"])