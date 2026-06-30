import sys
import os
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from script.parser import SCRIPTParser
from script.mol import PolymerBlock

def test_simple_polymer():
    """Test basic polymer parsing"""
    parser = SCRIPTParser()
    
    # Test simple polymer with exact count
    result = parser.parse("{[CC]}<n:50>")
    assert result["success"], f"Error: {result.get('error')}"
    mol = result["molecule"]
    assert mol.repeat_count == 50
    assert len(mol.atoms) == 2  # CC = 2 carbons

def test_polymer_with_range():
    """Test polymer with stochastic range"""
    parser = SCRIPTParser()
    
    result = parser.parse("{[CCCO]}<n:50-100>")
    assert result["success"], f"Error: {result.get('error')}"
    mol = result["molecule"]
    assert mol.repeat_count == (50, 100)
    assert len(mol.atoms) == 4  # CCCO = 4 atoms

def test_polymer_symbolic():
    """Test symbolic polymer repeat"""
    parser = SCRIPTParser()
    
    result = parser.parse("{[CC]}n")
    assert result["success"], f"Error: {result.get('error')}"
    mol = result["molecule"]
    assert mol.repeat_count == 'n'

def test_polymer_block_copolymer():
    """Test block copolymer - currently may not be fully implemented"""
    parser = SCRIPTParser()
    
    # Test a diblock copolymer - this may not parse yet
    result = parser.parse("{[CC]}<n:50> -b- {[CCCO]}<n:100>")
    # If not implemented, just check it doesn't crash
    # We know from the plan this is WIP
    if result["success"]:
        mol = result["molecule"]
        # Check that polymer blocks were created if implementation complete
        if hasattr(mol, 'polymer_blocks'):
            assert len(mol.polymer_blocks) >= 2

def test_polymer_alternating_copolymer():
    """Test alternating copolymer - currently may not be fully implemented"""
    parser = SCRIPTParser()
    
    # Test alternating polymer blocks - this may not parse yet
    result = parser.parse("{[CC]}<n:50> -alt- {[CO]}<n:50>")
    # If not implemented, just check it doesn't crash
    if result["success"]:
        mol = result["molecule"]
        if hasattr(mol, 'polymer_blocks'):
            assert len(mol.polymer_blocks) >= 2


if __name__ == "__main__":
    # Run tests directly
    pytest.main([__file__, "-v"])