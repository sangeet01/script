"""
SCRIPT Parser - Parse SCRIPT strings to RDKit Mol objects
"""

import os
import re
from pathlib import Path
from lark import Lark, Transformer, v_args, Token, Tree
from typing import List, Dict, Any, Optional, Tuple

from .peptide import PeptideHandler
from .local_rings import LocalRingHandler
from .validator import SCRIPTValidator


from .state_machine import GenerativeStateMachine
from lark.visitors import Interpreter

class SCRIPTInterpreter(Interpreter):
    """Interpret parse tree into molecular graph via GenerativeStateMachine"""
    
    def __init__(self):
        super().__init__()
        self.state = GenerativeStateMachine()
        self._next_bond_order = -1
        self._next_bond_dir = 0
        
    def start(self, tree):
        self.visit_children(tree)
        return self.state.mol

    def molecular_chain(self, tree):
        for child in tree.children:
            if not isinstance(child, Tree): continue
            
            # Use raw rule names as we patched them with !
            data = child.data.lstrip('!')
            if data == 'atom_expr':
                self.visit(child)
                self._next_bond_order = -1 
                self._next_bond_dir = 0
            elif data == 'bond':
                self.visit(child)
            elif data == 'local_ring':
                self.visit(child)
                self._next_bond_order = -1 
                self._next_bond_dir = 0
            elif data == 'branch':
                self.state.open_branch()
                self.visit(child)
                self.state.close_branch()

    def branch_content(self, tree):
        self.visit_children(tree)

    def bond(self, tree):
        self._next_bond_order, self._next_bond_dir = self._get_bond_info(tree)

    def _get_bond_info(self, bond_node):
        # Extract the bond type from children (Tokens)
        tokens = [t for t in bond_node.scan_values(lambda x: isinstance(x, Token))]
        if not tokens: return 1, 0
        
        t_type = tokens[0].type
        order = 1
        direction = 0
        if t_type == 'DOUBLE_BOND': order = 2
        elif t_type == 'TRIPLE_BOND': order = 3
        elif t_type == 'AROMATIC_BOND': order = 4
        elif t_type == 'EXPLICIT_MOBILE': order = 4  # Treat explicit mobile =: as aromatic/delocalized internally
        elif t_type == 'UP_BOND': direction = 3
        elif t_type == 'DOWN_BOND': direction = 4
        return order, direction

    def atom_expr(self, tree):
        # Find the actual atom child
        for child in tree.children:
            if isinstance(child, Tree):
                data = child.data.lstrip('!')
                if data == 'bracket_atom':
                    self._handle_bracket_atom(child)
                    return
            else:
                # Organic atom or wildcard
                symbol = str(child)
                self.state.add_atom(symbol, bond_order=self._next_bond_order, bond_dir=self._next_bond_dir)
                return
            
    def _handle_bracket_atom(self, node):
        element = "C"
        isotope = 0
        chiral = None
        hcount = None
        charge = 0
        
        for child in node.children:
            if not isinstance(child, Tree): continue
            t = child.data.lstrip('!')
            
            val = "".join([str(leaf) for leaf in child.scan_values(lambda x: not isinstance(x, Tree))])
            
            if t == 'element': element = val
            elif t == 'isotope': isotope = int(val) if val else 0
            elif t == 'chiral': chiral = val
            elif t == 'hcount': hcount = self._parse_hcount(val)
            elif t == 'charge': charge = self._parse_charge(val)
            
        self.state.add_atom(element, charge=charge, isotope=isotope, 
                            hcount=hcount, chiral=chiral,
                            bond_order=self._next_bond_order,
                            bond_dir=self._next_bond_dir,
                            is_bracket=True)

    def _parse_hcount(self, h_str: str) -> int:
        if h_str == "H": return 1
        if h_str.startswith("H"):
            try: return int(h_str[1:])
            except: return 1
        return 0

    def _parse_charge(self, charge_str: str) -> int:
        if charge_str == "+": return 1
        if charge_str == "-": return -1
        if charge_str == "++": return 2
        if charge_str == "--": return -2
        if charge_str.startswith("+"): 
            try: return int(charge_str[1:])
            except: return 1
        if charge_str.startswith("-"):
            try: return -int(charge_str[1:])
            except: return -1
        return 0

    def local_ring(self, tree):
        # Check for V2 ring closures: &INT: or &INT.
        ring_closure_nodes = list(tree.find_data('ring_closure'))
        if ring_closure_nodes:
            # It's an Anubandha ring closure
            node = ring_closure_nodes[0]
            
            # The children of ring_closure are the parsed tokens/trees.
            # However, because of the '!' patch, we scan all values.
            tokens = [str(t) for t in node.scan_values(lambda x: isinstance(x, Token))]
            
            ring_size_str = "".join([t for t in tokens if t.isdigit()])
            ring_size = int(ring_size_str) if ring_size_str else 0
            
            is_resonant = ":" in tokens
            
            self.state.add_v2_ring(ring_size, is_resonant, bond_order=self._next_bond_order)
            return

        # Look for named ring or digits (Legacy)
        named = list(tree.find_data('named_ring'))
        if named:
            letter = "".join([str(t) for t in named[0].scan_values(lambda x: isinstance(x, Token)) if str(t).isalpha()])
            self.state.add_ring(letter, bond_order=self._next_bond_order)
        else:
            ring_num = self._parse_ring_num(tree)
            self.state.add_ring(ring_num, bond_order=self._next_bond_order)

    def _parse_ring_num(self, tree) -> int:
        # Collect all tokens in order
        tokens = [str(t) for t in tree.scan_values(lambda x: isinstance(x, Token))]
        if not tokens: return 0
        digits = "".join([t for t in tokens if t.isdigit()])
        return int(digits) if digits else 0

class SCRIPTParser:
    """Main SCRIPT parser class - converts SCRIPT strings to CoreMolecule graphs"""
    
    def __init__(self):
        grammar_path = Path(__file__).parent / "grammar.lark"
        with open(grammar_path, 'r') as f:
            grammar_content = f.read()
            
        # Patch grammar: enforce preservation of ALL anonymous terminals
        # Prepend '!' to all lowercase rules (non-terminals)
        grammar_content = re.sub(rf"(?m)^([a-z][a-z0-9_]*)\s*:", r"!\1:", grammar_content)
        
        self.parser = Lark(grammar_content, start='start', parser='lalr')
        self.interpreter = SCRIPTInterpreter()
    
    def parse(self, script_string: str) -> Dict[str, Any]:
        try:
            tree = self.parser.parse(script_string)
            self.interpreter.state = GenerativeStateMachine()
            self.interpreter._next_bond_order = -1
            self.interpreter._next_bond_dir = 0
            mol = self.interpreter.visit(tree)
            
            return {
                "success": True,
                "molecule": mol,
                "error": None
            }
        except Exception as e:
            return {
                "success": False,
                "molecule": None,
                "error": str(e)
            }