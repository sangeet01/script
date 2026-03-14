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
from .chiral import ChiralResolver
from .mol import CoreMolecule


from .state_machine import GenerativeStateMachine
from lark.visitors import Interpreter

class SCRIPTInterpreter(Interpreter):
    """Interpret parse tree into molecular graph via GenerativeStateMachine"""
    
    def __init__(self):
        super().__init__()
        self.state = GenerativeStateMachine()
        self._next_bond_order = -1
        self._next_bond_dir = 0
        self._next_hapticity = 0
        self._next_bond_class = ""
        
    def entry(self, tree):
        res = self.visit_children(tree)
        # res can be a MacroscopicSystem if visit_children hit macroscopic_structure
        return res[0] if res else None

    def macroscopic_structure(self, tree):
        context = None
        entities = []
        
        for child in tree.children:
            if isinstance(child, Token):
                # CONTEXT_LABEL token or VBAR token - skip VBAR
                continue
            if not isinstance(child, Tree): continue
            t = child.data.lstrip('!')
            if t == 'macroscopic_context':
                # Scan all token values inside the context node
                for tok in child.children:
                    if isinstance(tok, Token) and tok.type == 'CONTEXT_LABEL':
                        context = str(tok)
                        break
                    elif isinstance(tok, Token):
                        # grammar patched with ! renames it - just grab first non-bracket token
                        val = str(tok)
                        if val not in ('[[', ']]'):
                            context = val
                            break
            elif t in ('reaction', 'script'):
                mols = self.visit(child)
                def apply_context(obj):
                    if isinstance(obj, list):
                        for item in obj: apply_context(item)
                    elif isinstance(obj, CoreMolecule):
                        obj.macroscopic_context = context
                
                apply_context(mols)
                entities.append(mols)
        
        if len(entities) == 1:
            return entities[0]
        return entities if entities else None

    def reaction(self, tree):
        # Reactions return a list of "sides", each side is a list of molecules (scripts)
        sides = []
        for child in tree.children:
            if isinstance(child, Tree) and child.data.lstrip('!') == 'script':
                # We DON'T reset state here because script() itself will 
                # manage states for its components.
                res = self.visit(child)
                if isinstance(res, CoreMolecule):
                    sides.append([res])
                else:
                    sides.append(res)
        return sides

    def script(self, tree):
        # Visit children and collect components
        molecules = []
        for child in tree.children:
            if isinstance(child, Tree) and child.data.lstrip('!') == 'component':
                # Each component starts fresh if it's separated by "." or "~"
                # but wait, the first component also needs a fresh state 
                # relative to whatever was in self.state before.
                self.state = GenerativeStateMachine()
                self.visit(child)
                self.state.finalize_valences()
                molecules.append(self.state.mol)
        
        if len(molecules) == 1:
            return molecules[0]
        
        return molecules

    def fragment_separator(self, tree):
        # We don't do anything here, handled in script() loop above
        # for proper state management between components.
        pass

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
                count = self._get_multiplier(child)
                for _ in range(count):
                    self.state.open_branch()
                    self.visit(child)
                    self.state.close_branch()
            elif data == 'peptide_chain':
                ph = PeptideHandler(self.state)
                ph.handle(child)
            elif data == 'polymer':
                self.visit_children(child)
            elif data == 'polymer_block':
                # Polymer blocks restart components or just continue?
                # User example: [ {PEG} ]<n:50-60> -b- [ {Styrene} ]<n:100>
                # We can treat each block as a component with a special multiplier
                mol_data = self.visit(child.children[0]) # (molecular_chain | peptide_chain | script)
                # multiplier is in child.children[1] or missing
                # For now just collect it.
                pass

    def branch_content(self, tree):
        self.visit_children(tree)

    def bond(self, tree):
        self._next_bond_order, self._next_bond_dir, self._next_hapticity, self._next_bond_class = self._get_bond_info(tree)

    def hbond(self, tree):
        # hbond = STAR_BOND INT  -> haptic bond with explicit hapticity number
        tokens = [t for t in tree.scan_values(lambda x: isinstance(x, Token))]
        hapticity = 0
        for tok in tokens:
            if tok.type == 'INT':
                try: hapticity = int(str(tok))
                except: pass
        self._next_bond_order = 4
        self._next_bond_class = "star"
        self._next_hapticity = hapticity

    def _get_bond_info(self, bond_node):
        # Extract bond type, direction, and bond_class from a bond node.
        # Note: hapticity (eta-n) is now handled by hbond() directly.
        tokens = [t for t in bond_node.scan_values(lambda x: isinstance(x, Token))]
        if not tokens: return 1, 0, 0, ""
        
        t_type = tokens[0].type
        hapticity = 0
        
        order = 1
        direction = 0
        bond_class = ""
        if t_type == 'DOUBLE_BOND': order = 2
        elif t_type == 'TRIPLE_BOND': order = 3
        elif t_type == 'AROMATIC_BOND': order = 4
        elif t_type == 'EXPLICIT_MOBILE': order = 4
        elif t_type == 'UP_BOND': direction = 3
        elif t_type == 'DOWN_BOND': direction = 4
        elif t_type == 'COORDINATE_BOND': order = 1; bond_class = "coordinate"
        elif t_type == 'STAR_BOND': order = 4; bond_class = "star"
        elif t_type == 'DATIVE': order = 1; bond_class = "dative"
        elif t_type == 'REV_DATIVE': order = 1; bond_class = "rev_dative"
        
        return order, direction, hapticity, bond_class

    def atom_expr(self, tree):
        count = self._get_multiplier(tree)
        
        # Find the actual atom child
        atom_node = None
        for child in tree.children:
            if isinstance(child, Tree):
                data = child.data.lstrip('!')
                if data == 'bracket_atom' or data == 'dhatu':
                    atom_node = child
                    break
            elif isinstance(child, Token) and child.type in ('ORGANIC_ATOM', 'WILDCARD'):
                atom_node = child
                break
        
        if atom_node is None: return

        for _ in range(count):
            if isinstance(atom_node, Tree):
                data = atom_node.data.lstrip('!')
                if data == 'bracket_atom':
                    self._handle_bracket_atom(atom_node)
                elif data == 'dhatu':
                    self._handle_dhatu(atom_node)
            else:
                symbol = str(atom_node)
                self.state.add_atom(symbol, bond_order=self._next_bond_order, bond_dir=self._next_bond_dir)
            
            self._next_bond_order = -1
            self._next_bond_dir = 0
            self._next_hapticity = 0
            self._next_bond_class = ""
            
    def _handle_bracket_atom(self, node):
        element = "C"
        isotope = 0
        chiral = None
        hcount = None
        charge = 0
        mapping = 0
        
        for child in node.children:
            if not isinstance(child, Tree): continue
            t = child.data.lstrip('!')
            
            val = "".join([str(leaf) for leaf in child.scan_values(lambda x: not isinstance(x, Tree))])
            
            if t == 'element': element = val
            elif t == 'isotope': isotope = int(val) if val else 0
            elif t == 'chiral': chiral = val
            elif t == 'hcount': hcount = self._parse_hcount(val)
            elif t == 'charge': charge = self._parse_charge(val)
            elif t == 'ring_class': 
                # :INT, extract digits
                mapping = int("".join(filter(str.isdigit, val)))
            
        self.state.add_atom(element, charge=charge, isotope=isotope, 
                            hcount=hcount, chiral=chiral,
                            bond_order=self._next_bond_order,
                            bond_dir=self._next_bond_dir,
                            is_bracket=True, mapping=mapping)

    def _handle_dhatu(self, node):
        element = "C"
        isotope = 0
        charge = 0
        hcount = None
        chiral = None
        mapping = 0
        occupancy = 1.0
        spin = 0
        is_excited = False
        
        has_bracket_attr = False
        for child in node.children:
            if isinstance(child, Token):
                if child.type == 'ATOM': element = str(child)
                continue
                
            if isinstance(child, Tree):
                t = child.data.lstrip('!')
                val = "".join([str(leaf) for leaf in child.scan_values(lambda x: not isinstance(x, Tree))])
                if t == 'state_block':
                    for sub in child.children:
                        if not isinstance(sub, Tree): continue
                        st = sub.data.lstrip('!')
                        sval = "".join([str(leaf) for leaf in sub.scan_values(lambda x: not isinstance(x, Tree))])
                        if st == 'isotope': isotope = int(sval) if sval else 0
                        elif st == 'charge': charge = self._parse_charge(sval)
                        elif st == 'hcount': hcount = self._parse_hcount(sval)
                        elif sval.startswith('~'):
                            try: occupancy = float(sval[1:])
                            except: occupancy = 1.0
                        elif sval.startswith('s:'):
                            try: spin = int(sval[2:])
                            except: spin = 0
                        elif sval == '*':
                            is_excited = True
                    has_bracket_attr = True
                elif t == 'chiral':
                    has_bracket_attr = True
                    chiral = "".join([str(leaf) for leaf in child.scan_values(lambda x: not isinstance(x, Tree))])
                elif t == 'ring_class':
                    has_bracket_attr = True
                    mapping = int("".join(filter(str.isdigit, val)))

        self.state.add_atom(element, charge=charge, isotope=isotope,
                            hcount=hcount, chiral=chiral,
                            bond_order=self._next_bond_order,
                            bond_dir=self._next_bond_dir,
                            is_bracket=has_bracket_attr,
                            mapping=mapping,
                            occupancy=occupancy, spin=spin, is_excited=is_excited)

    def _get_multiplier(self, tree) -> int:
        for child in tree.children:
            if isinstance(child, Tree) and child.data.lstrip('!') == 'multiplier':
                for node in child.children:
                    if isinstance(node, Token) and node.type == 'INT':
                        return int(str(node))
        return 1

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
        with open(grammar_path, 'r', encoding='utf-8') as f:
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
            # SCRIPTInterpreter.visit(tree) will return either a mol or a list of mols
            res = self.interpreter.visit(tree)
            
            # Resolve chirality as a post-pass (Paninian Sandhi)
            # res can be CoreMolecule, List[CoreMolecule] (script), or List[List[CoreMolecule]] (reaction)
            def resolve_all(obj):
                if isinstance(obj, list):
                    for item in obj: resolve_all(item)
                elif isinstance(obj, CoreMolecule):
                    # Valences are finalized during visit() now, but we'll be safe
                    # Only resolved if not already handled
                    ChiralResolver(obj).resolve()

            resolve_all(res)

            return {
                "success": True,
                "molecule": res, # Can be CoreMolecule or List[CoreMolecule]
                "error": None
            }
        except Exception as e:
            return {
                "success": False,
                "molecule": None,
                "error": str(e)
            }