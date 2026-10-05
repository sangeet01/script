"""
SCRIPT Grammar - Optional Lark/PLY parser from EBNF
"""

import os
from pathlib import Path

# SCRIPT Grammar in Lark EBNF format
SCRIPT_GRAMMAR = '''
// SCRIPT.lark - Stub grammar for grammar.py module
// WARNING: this is NOT the production grammar.
// The production grammar is grammar.lark (33KB).
// Based on SCRIPT specification v1.0 (organic subset only)

%import common (WS, NUMBER, SIGNED_NUMBER)
%ignore WS

// === Start rule ===
start: script

// === Top-level: multi-component with . separator ===
script: component ("." component)*

component: molecular_chain
         | peptide_block
         | polymer_block

// === Core molecular chain ===
molecular_chain: chain_element+

chain_element: atom_expr
             | bond
             | local_ring
             | branch

atom_expr: organic_atom
         | bracket_atom

// Organic subset (no brackets needed)
organic_atom: "B" "r"?
            | "C" "l"?
            | "N"
            | "O" 
            | "P"
            | "S"
            | "F"
            | "I"

// Bracket atoms [element@charge]
bracket_atom: "[" isotope? element chiral? hcount? charge? ring_class? "]"

element: /[A-Z][a-z]?/
isotope: NUMBER
chiral: "@" | "@@" | "@SP" | "@OH" | "@AX"
hcount: "H" NUMBER?
charge: "+" NUMBER? | "-" NUMBER? | "++" | "--"
ring_class: ":" /[A-Z]/

// Chemical bonds
bond: "-" | "=" | "#" | ":" | "/" | "\\\\" | "->"

// Branches in parentheses
branch: "(" molecular_chain ")"

// Local ring closures (digit after atom)
local_ring: DIGIT | "%" DIGIT DIGIT

// === Peptide blocks ===
peptide_block: "{" amino_sequence "}" local_ring?

amino_sequence: amino_acid ("." amino_acid)*

amino_acid: AA_CODE | PTM_CODE | bridge_ref

AA_CODE: "A"|"R"|"N"|"D"|"C"|"E"|"Q"|"G"|"H"|"I"|"L"|"K"|"M"|"F"|"P"|"S"|"T"|"W"|"Y"|"V"
PTM_CODE: "pS" | "pT" | "pY" | "mK"
bridge_ref: "[" /[A-Z]/ "]"

// === Polymer blocks ===
polymer_block: "{" molecular_chain "}" repeat_spec?

repeat_spec: "n" | NUMBER | stochastic
stochastic: NUMBER "." NUMBER

// === Terminals ===
DIGIT: "0".."9"
NUMBER: DIGIT+

// Comments (optional)
COMMENT: /#.*/
%ignore COMMENT
'''

def get_grammar() -> str:
    """Get the SCRIPT grammar string"""
    return SCRIPT_GRAMMAR

def save_grammar_file(filepath: str = None) -> str:
    """Save the stub grammar to a file and return the path.

    WARNING: never call without an explicit filepath.  The default (None)
    would write next to this module -- potentially overwriting grammar.lark,
    the 33 KB production grammar.  Pass an explicit path such as a temp dir.
    """
    if filepath is None:
        raise ValueError(
            "save_grammar_file() requires an explicit filepath.  "
            "Calling without one risks overwriting grammar.lark, the "
            "production LALR grammar.  Pass a safe destination, e.g. "
            "a temporary directory."
        )

    with open(filepath, 'w') as f:
        f.write(SCRIPT_GRAMMAR)

    return str(filepath)

def load_grammar_from_file(filepath: str = None) -> str:
    """Load grammar from file"""
    if filepath is None:
        filepath = Path(__file__).parent / "grammar.lark"
    
    if os.path.exists(filepath):
        with open(filepath, 'r') as f:
            return f.read()
    else:
        # Return built-in grammar
        return SCRIPT_GRAMMAR

# Grammar validation and testing
def test_grammar():
    """Test grammar with sample inputs"""
    try:
        from lark import Lark
        
        parser = Lark(SCRIPT_GRAMMAR, start='start', parser='lalr')
        
        test_cases = [
            "C",           # Methane
            "CCO",         # Ethanol
            "C=O",         # Formaldehyde
            "C1CCCCC6",    # Cyclohexane (local ring)
            "{A}",         # Alanine
            "{A.G.S}",     # Tripeptide
            "CC{A}O",      # Mixed structure
        ]
        
        results = {}
        for case in test_cases:
            try:
                tree = parser.parse(case)
                results[case] = "OK"
            except Exception as e:
                results[case] = "FAIL: " + str(e)
        
        return results
    
    except ImportError:
        return {"error": "Lark not available for testing"}

if __name__ == "__main__":
    # Test the grammar
    print("Testing SCRIPT grammar...")
    results = test_grammar()
    
    for case, result in results.items():
        print("  {:12} -> {}".format(case, result))
    
    # Save grammar file
    grammar_path = save_grammar_file()
    print(f"\nGrammar saved to: {grammar_path}")