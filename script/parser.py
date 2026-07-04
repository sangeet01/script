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
from .mol import CoreMolecule, Reaction


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
        self._next_translation = (0, 0, 0)   # V3.6 periodic topology
        
    def entry(self, tree):
        res = self.visit_children(tree)
        # res can be a MacroscopicSystem if visit_children hit macroscopic_structure
        return res[0] if res else None

    def macroscopic_structure(self, tree):
        context = None
        lattice = None   # V3.6: lattice vectors from context block
        entities = []
        phase_labels = []  # labels between VBAR tokens

        for child in tree.children:
            if isinstance(child, Token):
                if child.type == 'VBAR':
                    phase_labels.append('|')
                continue
            if not isinstance(child, Tree): continue
            t = child.data.lstrip('!')
            if t == 'macroscopic_context':
                lattice = None
                for tok in child.children:
                    if isinstance(tok, Token) and tok.type == 'CONTEXT_LABEL':
                        context = str(tok)
                    elif isinstance(tok, Tree) and tok.data.lstrip('!') == 'lattice_params':
                        # Parse a,b,c,alpha,beta,gamma from 6 FLOAT tokens
                        floats = [float(str(f)) for f in tok.scan_values(
                            lambda x: isinstance(x, Token) and x.type == 'FLOAT')]
                        if len(floats) == 6:
                            a, b, c, alpha, beta, gamma = floats
                            import math
                            # Compute 3x3 lattice vectors from cell parameters
                            ar = math.radians(alpha)
                            br = math.radians(beta)
                            gr = math.radians(gamma)
                            cos_a = math.cos(ar); cos_b = math.cos(br)
                            cos_g = math.cos(gr); sin_g = math.sin(gr)
                            vol_fac = math.sqrt(max(0, 1 - cos_a**2 - cos_b**2 - cos_g**2
                                                    + 2*cos_a*cos_b*cos_g))
                            lattice = (
                                (a,           0,                                       0),
                                (b*cos_g,     b*sin_g,                                 0),
                                (c*cos_b,     c*(cos_a - cos_b*cos_g)/sin_g,
                                              c*vol_fac/sin_g),
                            )
                    elif isinstance(tok, Token) and tok.type not in ('[[', ']]', ';'):
                        if context is None:
                            context = str(tok)
            elif t in ('reaction', 'script'):
                mols = self.visit(child)
                def apply_context(obj, ctx=context, lat=lattice):
                    if isinstance(obj, list):
                        for item in obj: apply_context(item, ctx, lat)
                    elif isinstance(obj, CoreMolecule):
                        obj.macroscopic_context = ctx
                        if lat is not None:
                            obj.lattice_vectors = lat
                apply_context(mols)

                # Tag phase boundary: each entity after the first gets a label
                phase_idx = len(entities)
                if phase_idx > 0 and phase_labels:
                    def apply_phase(obj, label='|'):
                        if isinstance(obj, list):
                            for item in obj: apply_phase(item, label)
                        elif isinstance(obj, CoreMolecule):
                            obj.phase_boundary = label
                    apply_phase(mols)

                entities.append(mols)

        if len(entities) == 1:
            return entities[0]
        return entities if entities else None

    def script(self, tree):
        # Visit children and collect components, preserving separator semantics
        molecules = []
        separators = []  # separator token between each pair of components

        for child in tree.children:
            if not isinstance(child, Tree):
                continue
            data = child.data.lstrip('!')
            if data == 'component':
                self.state = GenerativeStateMachine()
                self.visit(child)
                self.state.finalize_valences()
                # Sandhi gate: propagate any valence violations recorded during
                # add_atom -> add_bond as a hard parse failure.
                if self.state.valence_violations:
                    u, v, order = self.state.valence_violations[0]
                    atoms = self.state.mol.atoms
                    sym_u = atoms[u].symbol if u < len(atoms) else '?'
                    sym_v = atoms[v].symbol if v < len(atoms) else '?'
                    raise ValueError(
                        f"Sandhi violation: bond {sym_u}({u})-{sym_v}({v}) "
                        f"order={order} exceeds maximum valence"
                    )
                mol = self.state.mol
                # Apply pending separator from the preceding fragment_separator
                if separators:
                    mol.fragment_separator = separators[-1]
                molecules.append(mol)
            elif data == 'fragment_separator':
                # Capture the actual token: '.' = solvate, '~' = ionic pair
                sep_tokens = list(child.scan_values(lambda t: True))
                sep = str(sep_tokens[0]) if sep_tokens else '.'
                separators.append(sep)

        if len(molecules) == 1:
            return molecules[0]
        return molecules

    def reaction(self, tree):
        """Return a Reaction object instead of a raw list-of-lists."""
        from .mol import Reaction
        sides = []
        for child in tree.children:
            if isinstance(child, Tree) and child.data.lstrip('!') == 'script':
                res = self.visit(child)
                if isinstance(res, list):
                    sides.append(res)
                elif res is not None:
                    sides.append([res])
                else:
                    sides.append([])
        # SCRIPT supports: reactants >> products  (2 sides)
        #              or: reactants > agents > products  (3 sides, middle = agents)
        if len(sides) == 2:
            return Reaction(reactants=sides[0], products=sides[1])
        elif len(sides) == 3:
            return Reaction(reactants=sides[0], agents=sides[1], products=sides[2])
        elif len(sides) == 1:
            return Reaction(reactants=sides[0], products=[])
        return Reaction(reactants=[], products=[])

    def molecular_chain(self, tree):
        for child in tree.children:
            if not isinstance(child, Tree): continue
            
            # Use raw rule names as we patched them with !
            data = child.data.lstrip('!')
            if data == 'atom_expr':
                self.visit(child)
                self._next_bond_order = -1 
                self._next_bond_dir = 0
                self._next_translation = (0, 0, 0)
            elif data == 'bond':
                self.visit(child)
            elif data == 'maybe_periodic_bond':
                # V3.6: bond with optional @tx,ty,tz translation vector
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
            elif data == 'polymer':
                # Visit polymer_unit to build the molecular graph,
                # then extract repeat_count from polymer_suffix if present
                for sub in child.children:
                    if not isinstance(sub, Tree): continue
                    sub_data = sub.data.lstrip('!')
                    if sub_data == 'polymer_unit':
                        self.visit_children(sub)
                    elif sub_data == 'polymer_suffix':
                        self._parse_polymer_suffix(sub)
            elif data == 'polymer_block':
                self._parse_polymer_block(child)

    def peptide_chain(self, tree):
        """Handle peptide_chain whether as direct component or inside molecular_chain."""
        ph = PeptideHandler(self.state)
        ph.handle(tree)

    def polymer(self, tree):
        """Handle top-level polymer node (direct child of component)."""
        from lark import Tree as LTree, Token as LToken
        for child in tree.children:
            if not isinstance(child, LTree): continue
            sub = child.data.lstrip('!')
            if sub == 'polymer_unit':
                self.visit_children(child)
            elif sub == 'polymer_suffix':
                self._parse_polymer_suffix(child)

    def polymer_block(self, tree):
        """Handle top-level polymer_block node (block copolymer at component level).

        The grammar rule is:
            polymer_block: polymer polymer_junction polymer (polymer_junction polymer)*

        Each polymer child is parsed as its own unit; junction tokens are extracted
        to set block_kind on each resulting PolymerBlock object.
        """
        from lark import Tree as LTree, Token as LToken
        from .mol import PolymerBlock
        from .state_machine import GenerativeStateMachine

        _KIND_MAP = {
            '-b-': 'diblock',
            '-a-': 'alternating', '-alt-': 'alternating',
            '-r-': 'random', '-ran-': 'random', '-stat-': 'random',
            '-g-': 'graft',
        }

        current_kind = ''
        for child in tree.children:
            if isinstance(child, LToken):
                current_kind = _KIND_MAP.get(str(child), str(child).strip('-'))
                continue
            if not isinstance(child, LTree):
                continue
            node_data = child.data.lstrip('!')
            if node_data == 'polymer':
                # Parse this block unit in its own state machine
                saved_state = self.state
                self.state = GenerativeStateMachine()
                self.visit(child)
                self.state.finalize_valences()
                unit_mol = self.state.mol
                self.state = saved_state

                block = PolymerBlock(
                    unit=unit_mol,
                    repeat_count=unit_mol.repeat_count,
                    block_kind=current_kind,
                )
                self.state.mol.polymer_blocks.append(block)
                if self.state.mol.block_topology is None and current_kind:
                    self.state.mol.block_topology = current_kind


        self.visit_children(tree)

    def bond(self, tree):
        self._next_bond_order, self._next_bond_dir, self._next_hapticity, self._next_bond_class = self._get_bond_info(tree)
        self._next_translation = (0, 0, 0)   # reset on plain bond visit

    def maybe_periodic_bond(self, tree):
        """
        V3.6: bond optionally followed by @tx,ty,tz translation vector.
        Children: [bond_tree] or [bond_tree, signed_int, signed_int, signed_int]
        """
        from lark import Tree as LarkTree
        children = [c for c in tree.children if c is not None]
        # First child is always the bond subtree
        bond_tree = children[0]
        self._next_bond_order, self._next_bond_dir, self._next_hapticity, self._next_bond_class = self._get_bond_info(bond_tree)
        # If there are signed_int children, extract translation
        sint_nodes = [c for c in children[1:] if isinstance(c, LarkTree)]
        if len(sint_nodes) == 3:
            tx = self._parse_signed_int(sint_nodes[0])
            ty = self._parse_signed_int(sint_nodes[1])
            tz = self._parse_signed_int(sint_nodes[2])
            self._next_translation = (tx, ty, tz)
        else:
            self._next_translation = (0, 0, 0)

    def neg_int(self, tree):
        return -int(str(tree.children[0]))

    def pos_int(self, tree):
        return int(str(tree.children[0]))

    def _parse_signed_int(self, node) -> int:
        """Parse a signed_int tree node to Python int.
        neg_int children: [Token(SINGLE_BOND,'-'), Token(INT,'N')] — take last child.
        pos_int children: [Token(INT,'N')] — take first (only) child.
        """
        from lark import Tree, Token
        if isinstance(node, Tree):
            if node.data == 'neg_int':
                # Last child is always the INT token; first child is '-' (SINGLE_BOND)
                return -int(str(node.children[-1]))
            elif node.data == 'pos_int':
                return int(str(node.children[0]))
        if isinstance(node, Token):
            return int(str(node))
        return 0

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
        elif t_type == 'EXPLICIT_MOBILE': order = 1; bond_class = 'tautomeric'
        elif t_type == 'UP_BOND': direction = 3
        elif t_type == 'DOWN_BOND': direction = 4
        elif t_type == 'COORDINATE_BOND': order = 1; bond_class = "coordinate"
        elif t_type == 'STAR_BOND': order = 4; bond_class = "star"
        elif t_type == 'DATIVE': order = 1; bond_class = "dative"
        elif t_type == 'REV_DATIVE': order = 1; bond_class = "rev_dative"
        
        return order, direction, hapticity, bond_class

    def atom_expr(self, tree):
        count = self._get_multiplier(tree)
        legacy_ring_count = count if self._is_legacy_ring_closure_candidate(tree, count) else None
        if legacy_ring_count is not None:
            count = 1
        
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
                is_wildcard = (atom_node.type == 'WILDCARD' or symbol == '*')
                self.state.add_atom(
                    symbol,
                    bond_order=self._next_bond_order,
                    bond_dir=self._next_bond_dir,
                    bond_class=self._next_bond_class,
                    translation=self._next_translation,
                )
                if is_wildcard:
                    self.state.mol.atoms[-1].is_wildcard = True
                    self.state.mol.atoms[-1].atomic_num = 0
            
            self._next_bond_order = -1
            self._next_bond_dir = 0
            self._next_hapticity = 0
            self._next_bond_class = ""
            self._next_translation = (0, 0, 0)

        if legacy_ring_count is not None:
            self.state.add_ring(legacy_ring_count, bond_order=1)
            
    def _handle_bracket_atom(self, node):
        element = "C"
        isotope = 0
        chiral = None
        hcount = None
        charge = 0
        mapping = 0
        radical = 0
        is_wildcard = False
        is_query = False
        query_data = {}

        # Detect [*] form
        for child in node.children:
            if isinstance(child, Token) and str(child) == '*':
                is_wildcard = True
                element = '*'
                break

        # Detect [query_atom] form
        if not is_wildcard:
            for child in node.children:
                if isinstance(child, Tree) and child.data.lstrip('!') == 'query_atom':
                    is_query = True
                    query_data = self._parse_query_atom(child)
                    element = query_data.get('symbol', '*')
                    break

        if not is_wildcard and not is_query:
            for child in node.children:
                if not isinstance(child, Tree): continue
                t = child.data.lstrip('!')
                val = "".join([str(leaf) for leaf in child.scan_values(lambda x: not isinstance(x, Tree))])
                if t == 'element': element = val
                elif t == 'isotope': isotope = int(val) if val else 0
                elif t == 'chiral': chiral = val
                elif t == 'hcount': hcount = self._parse_hcount(val)
                elif t == 'charge': charge = self._parse_charge(val)
                elif t == 'radical':
                    radical = len([c for c in val if c == '.'])
                elif t == 'ring_class':
                    mapping = int("".join(filter(str.isdigit, val)))

        self.state.add_atom(element, charge=charge, isotope=isotope,
                            hcount=hcount, chiral=chiral,
                            bond_order=self._next_bond_order,
                            bond_dir=self._next_bond_dir,
                            bond_class=self._next_bond_class,
                            is_bracket=True, mapping=mapping,
                            radical=radical,
                            translation=self._next_translation)
        atom = self.state.mol.atoms[-1]
        if is_wildcard:
            atom.is_wildcard = True
            atom.atomic_num = 0
        self._next_translation = (0, 0, 0)
        if is_query:
            atom.is_query = True
            atom.is_wildcard = True  # query atoms are wildcards in canonical sense
            atom.atomic_num = 0
            for k, v in query_data.items():
                if hasattr(atom, k):
                    setattr(atom, k, v)

    def _parse_query_atom(self, query_tree) -> dict:
        """Parse a query_atom tree into a dict of query fields for CoreAtom."""
        from .mol import ATOMIC_NUM_TO_SYMBOL

        result = {'symbol': '*', 'query_atomic_nums': [], 'query_not': False,
                  'query_ring': None, 'query_valence': None,
                  'query_hcount': None, 'query_aromatic': None,
                  'query_primitives': []}
        self._parse_query_primitives(query_tree, result)
        # If only one atomic number, use it as the symbol for query readability
        if len(result['query_atomic_nums']) == 1:
            result['symbol'] = ATOMIC_NUM_TO_SYMBOL.get(result['query_atomic_nums'][0], '*')
        return result

    def _parse_query_primitives(self, tree, result: dict):
        """Recursively extract query primitives."""
        for child in tree.children:
            if not isinstance(child, Tree): continue
            t = child.data.lstrip('!')
            if t == 'query_primitive':
                tokens = list(child.scan_values(lambda x: True))
                if not tokens: continue
                first = str(tokens[0])
                vals = [str(t) for t in tokens]
                if first == '!':
                    result['query_not'] = True
                    # Recurse into negated primitive
                    self._parse_query_primitives(child, result)
                elif first == '#':
                    if len(vals) > 1:
                        try: result['query_atomic_nums'].append(int(vals[1]))
                        except ValueError: pass
                elif first == 'R':
                    result['query_ring'] = int(vals[1]) if len(vals) > 1 else 0
                elif first == 'r':
                    result['query_ring'] = int(vals[1]) if len(vals) > 1 else 0
                elif first == 'v':
                    if len(vals) > 1:
                        try: result['query_valence'] = int(vals[1])
                        except ValueError: pass
                elif first == 'H':
                    result['query_hcount'] = int(vals[1]) if len(vals) > 1 else 1
                elif first == 'a':
                    result['query_aromatic'] = True
                elif first == 'A':
                    result['query_aromatic'] = False
                else:
                    # Could be an element symbol
                    result['symbol'] = first
                result['query_primitives'].append({'type': first, 'vals': vals})
            elif t in ('query_atom', 'query_primitive'):
                self._parse_query_primitives(child, result)

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
                            bond_class=self._next_bond_class,
                            is_bracket=has_bracket_attr,
                            mapping=mapping,
                            occupancy=occupancy, spin=spin, is_excited=is_excited,
                            translation=self._next_translation)
        self._next_translation = (0, 0, 0)

    def _parse_polymer_suffix(self, suffix_tree) -> None:
        """Extract repeat count from polymer_suffix and store on current mol."""
        from lark import Token as LarkToken
        for child in suffix_tree.children:
            if not isinstance(child, Tree): continue
            sub = child.data.lstrip('!')
            if sub == 'polymer_multiplier':
                ints = [int(str(t)) for t in child.scan_values(
                    lambda t: isinstance(t, LarkToken) and t.type == 'INT')]
                if len(ints) == 1:
                    self.state.mol.repeat_count = ints[0]
                elif len(ints) == 2:
                    self.state.mol.repeat_count = (ints[0], ints[1])
            elif sub == 'repeat_spec':
                # bare 'n' = unspecified repeat, or INT for exact
                tokens = list(child.scan_values(lambda t: True))
                val = str(tokens[0]) if tokens else 'n'
                try:
                    self.state.mol.repeat_count = int(val)
                except ValueError:
                    self.state.mol.repeat_count = 'n'  # symbolic

    def _parse_polymer_block(self, block_tree) -> None:
        """
        Handle a polymer_block junction node in a block copolymer.

        The grammar produces polymer_block when consecutive polymer units are
        joined by a block-junction token (e.g. '-b-', '-alt-', '-ran-', '-g-').
        Each block is parsed into its own CoreMolecule unit and stored as a
        PolymerBlock in the parent molecule's polymer_blocks list.

        Recognised junction keywords:
            diblock / b     -> 'diblock'
            triblock        -> 'triblock'
            alternating/alt -> 'alternating'
            random/ran      -> 'random'
            graft/g         -> 'graft'

        The parent CoreMolecule (self.state.mol) receives:
            .polymer_blocks  - appended with the new PolymerBlock
            .block_topology  - set to the junction kind string
        """
        from .mol import PolymerBlock, CoreMolecule
        from .state_machine import GenerativeStateMachine
        from lark import Token as LarkToken

        # Map grammar tokens / literal strings to canonical kind names
        _KIND_MAP = {
            'b': 'diblock', 'diblock': 'diblock',
            'triblock': 'triblock',
            'alt': 'alternating', 'alternating': 'alternating',
            'ran': 'random', 'random': 'random',
            'g': 'graft', 'graft': 'graft',
        }

        parent_mol = self.state.mol
        block_kind = ''
        unit_mol = None
        block_repeat: Any = None

        for child in block_tree.children:
            if isinstance(child, LarkToken):
                raw = str(child).strip('-').lower()
                block_kind = _KIND_MAP.get(raw, raw)
                continue
            if not isinstance(child, Tree):
                continue
            sub_data = child.data.lstrip('!')

            if sub_data == 'polymer_unit':
                # Parse the block's repeat unit into a fresh graph
                saved_state = self.state
                self.state = GenerativeStateMachine()
                self.visit_children(child)
                self.state.finalize_valences()
                unit_mol = self.state.mol
                self.state = saved_state

            elif sub_data == 'polymer_suffix':
                # Extract per-block repeat count from the suffix node
                from lark import Token as LT
                for sc in child.children:
                    if not isinstance(sc, Tree): continue
                    ss = sc.data.lstrip('!')
                    if ss == 'polymer_multiplier':
                        ints = [int(str(t)) for t in sc.scan_values(
                            lambda t: isinstance(t, LT) and t.type == 'INT')]
                        if len(ints) == 1:
                            block_repeat = ints[0]
                        elif len(ints) == 2:
                            block_repeat = (ints[0], ints[1])
                    elif ss == 'repeat_spec':
                        tokens = list(sc.scan_values(lambda t: True))
                        val = str(tokens[0]) if tokens else 'n'
                        try:
                            block_repeat = int(val)
                        except ValueError:
                            block_repeat = 'n'

            elif sub_data in ('polymer', 'polymer_block'):
                # Nested or chained block — recurse so triblock etc. work
                if sub_data == 'polymer_block':
                    self._parse_polymer_block(child)
                else:
                    saved_state = self.state
                    self.state = GenerativeStateMachine()
                    self.visit(child)
                    self.state.finalize_valences()
                    unit_mol = self.state.mol
                    self.state = saved_state

        if unit_mol is not None:
            block = PolymerBlock(unit=unit_mol,
                                 repeat_count=block_repeat,
                                 block_kind=block_kind)
            parent_mol.polymer_blocks.append(block)
            # Update topology on parent: first block wins if already set
            if parent_mol.block_topology is None and block_kind:
                parent_mol.block_topology = block_kind


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

    def _attach_parents(self, tree, parent=None):
        if isinstance(tree, Tree):
            tree.parent = parent
            for child in tree.children:
                self._attach_parents(child, tree)

    def _is_legacy_ring_closure_candidate(self, tree, count: int) -> bool:
        if count < 3:
            return False
        parent = getattr(tree, 'parent', None)
        if parent is None or parent.data.lstrip('!') != 'molecular_chain':
            return False
        child_trees = [c for c in parent.children if isinstance(c, Tree)]
        return child_trees and child_trees[-1] is tree

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
        import re as _re
        from .mol import Reaction, CoreMolecule
        from .validator import SCRIPTValidator

        # 1. Pre-validation (grounds the probabilistic generator)
        if not SCRIPTValidator().is_valid(script_string):
            return {
                "success": False,
                "molecule": None,
                "atoms": [],
                "bonds": [],
                "reaction": None,
                "error": "Invalid SCRIPT syntax (Validator rejection)"
            }

        # Pre-process 3-part reactions: "R>A>P" (lone > not part of >>)
        _agent_parts = _re.split(r'(?<!>)>(?!>)', script_string)
        if len(_agent_parts) == 3:
            r_str, a_str, p_str = [s.strip() for s in _agent_parts]
            r_res = self.parse(r_str)
            a_res = self.parse(a_str) if a_str else {"success": True, "molecule": None, "reaction": None, "error": None}
            p_res = self.parse(p_str)
            if r_res["success"] and p_res["success"]:
                def _to_list(m):
                    if m is None: return []
                    return m if isinstance(m, list) else [m]
                rxn = Reaction(
                    reactants=_to_list(r_res["molecule"]),
                    agents=_to_list(a_res.get("molecule")),
                    products=_to_list(p_res["molecule"]),
                )
                return {
                    "success": True,
                    "molecule": rxn.reactants + rxn.agents + rxn.products,
                    "atoms": [], # Reactions don't have a single atom list
                    "bonds": [],
                    "reaction": rxn,
                    "error": None
                }

        try:
            from .mol import Reaction, CoreMolecule
            tree = self.parser.parse(script_string)
            self.interpreter._attach_parents(tree)
            self.interpreter.state = GenerativeStateMachine()
            self.interpreter._next_bond_order = -1
            self.interpreter._next_bond_dir = 0
            res = self.interpreter.visit(tree)

            # Resolve chirality post-pass (Paninian Sandhi)
            def resolve_all(obj):
                if isinstance(obj, Reaction):
                    for mol in obj.reactants + obj.agents + obj.products:
                        ChiralResolver(mol).resolve()
                elif isinstance(obj, list):
                    for item in obj: resolve_all(item)
                elif isinstance(obj, CoreMolecule):
                    ChiralResolver(obj).resolve()

            resolve_all(res)

            # If result is a Reaction, expose it in both 'molecule' and 'reaction' keys
            if isinstance(res, Reaction):
                return {
                    "success": True,
                    "molecule": res.reactants + res.agents + res.products,
                    "atoms": [],
                    "bonds": [],
                    "reaction": res,
                    "error": None
                }

            return {
                "success": True,
                "molecule": res,  # CoreMolecule or List[CoreMolecule]
                "atoms": res.atoms if isinstance(res, CoreMolecule) else [],
                "bonds": res.bonds if isinstance(res, CoreMolecule) else [],
                "reaction": None,
                "error": None
            }
        except Exception as e:
            return {
                "success": False,
                "molecule": None,
                "reaction": None,
                "error": str(e)
            }


def parse_script(script_string: str):
    """Parse SCRIPT string to internal representation"""
    parser = SCRIPTParser()
    return parser.parse(script_string)
