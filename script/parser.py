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
        # V4 Lattice Extension pending state
        self._pending_lattice_type = None     # str, set by [Lattice:T]
        self._pending_thickness = None        # (cls, args_tuple), set by [Thickness:..]
        self._pending_post_process = []       # list of (op_name, args_tuple)
        self._pending_typed_tags = []         # [V4.1] list of (ns, val, args)
        self._next_control_points = None      # [V4.2 Q5] list of (x,y,z) tuples
        
    def entry(self, tree):
        res = self.visit_children(tree)
        # res can be a MacroscopicSystem if visit_children hit macroscopic_structure
        return res[0] if res else None

    def macroscopic_structure(self, tree):
        context = None
        namespace = None   # [V4 L6] "geom" / "xtal" / "chem" / None
        lattice = None   # V3.6: lattice vectors from context block
        entities = []
        phase_labels = []  # labels between VBAR tokens

        # Reset pending post-process ops at the start of a structure
        self._pending_post_process = []

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
                    if isinstance(tok, Token) and tok.type == 'NAMESPACE':
                        # [V4 L6] "geom:" / "xtal:" / "chem:" prefix
                        ns = str(tok)
                        namespace = ns.rstrip(':')
                    elif isinstance(tok, Token) and tok.type == 'CONTEXT_LABEL':
                        context = str(tok)
                    elif isinstance(tok, Tree) and tok.data.lstrip('!') == 'lattice_params':
                        # Parse a,b,c,alpha,beta,gamma from 6 numeric tokens.
                        # lattice_val accepts FLOAT or INT (bare integers
                        # like "90" for right angles are common and should
                        # not require "90.0"), so both token types must be
                        # collected here or the extraction silently yields
                        # fewer than 6 values and lattice_vectors stays None.
                        floats = [float(str(f)) for f in tok.scan_values(
                            lambda x: isinstance(x, Token) and x.type in ('FLOAT', 'INT'))]
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
                def apply_context(obj, ctx=context, lat=lattice, ns=namespace):
                    if isinstance(obj, list):
                        for item in obj: apply_context(item, ctx, lat, ns)
                    elif isinstance(obj, CoreMolecule):
                        obj.macroscopic_context = ctx
                        if ns: obj.context_namespace = ns
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
            elif t == 'post_process':
                # [V4 L5] Visit and consume one |> op(args) operation.
                self.visit(child)

        # [V4 L5] Attach pending post-process ops to the final entity.
        if self._pending_post_process:
            if entities:
                target = entities[-1]
                def apply_pp(obj, ops=self._pending_post_process):
                    if isinstance(obj, list):
                        for item in obj: apply_pp(item, ops)
                    elif isinstance(obj, CoreMolecule):
                        obj.post_process_ops.extend(ops)
                apply_pp(target)
            self._pending_post_process = []

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
                self._next_control_points = None  # [V4.2 Q5]
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
                count = self._get_multiplier(child) or 1
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

        V3.7 fix: each 'polymer' child is visited via self.visit(child),
        which auto-dispatches to the polymer() method below. That method
        visits polymer_unit DIRECTLY into the shared self.state (not an
        isolated sub-state), so atoms across blocks share one contiguous
        index space in the same CoreMolecule. The junction bond between
        blocks needs no special-case code: as long as
        self.state.current_atom_idx is left pointing at the previous
        block's last atom when the next block's first atom_expr is
        visited, the existing add_atom/add_bond implicit-bond resolution
        (bond_order=-1 -> SINGLE, unless both atoms are aromatic) creates
        the junction bond automatically.

        (The previous version of this method built PolymerBlock metadata
        from a throwaway isolated GenerativeStateMachine per block — those
        objects never contributed atoms to the real graph; the actual
        merged graph came only from an unconditional visit_children(tree)
        call at the end. That worked by accident for connectivity but left
        block_kind and atom_start/atom_end unpopulated, since
        polymer_junction arrives here as a Tree node, not a bare Token.
        This version tracks atom ranges and junction kind correctly around
        the same working shared-state mechanism, and drops the redundant
        isolated-state loop entirely.)

        polymer_junction: "-b-" | "-alt-" | "-a-" | "-ran-" | "-stat-" | "-g-"
        """
        from lark import Tree as LTree, Token as LToken
        from .mol import PolymerBlock

        _KIND_MAP = {
            'b': 'diblock', 'diblock': 'diblock',
            'alt': 'alternating', 'a': 'alternating', 'alternating': 'alternating',
            'ran': 'random', 'r': 'random', 'stat': 'random', 'random': 'random',
            'g': 'graft', 'graft': 'graft',
        }

        parent_mol = self.state.mol
        pending_kind = ''

        for child in tree.children:
            if not isinstance(child, LTree):
                continue
            node_data = child.data.lstrip('!')

            if node_data == 'polymer_junction':
                raw_tokens = [str(t) for t in child.scan_values(lambda x: True)]
                raw = raw_tokens[0].strip('-').lower() if raw_tokens else ''
                pending_kind = _KIND_MAP.get(raw, raw)

            elif node_data == 'polymer':
                atom_start = len(parent_mol.atoms)

                # Reset bond state so this block's first atom doesn't
                # inherit stale bond-order/direction from tail tokens of
                # the previous block; current_atom_idx is deliberately
                # left untouched so the junction bond forms via the
                # normal implicit-bond mechanism in add_atom/add_bond.
                self._next_bond_order = -1
                self._next_bond_dir = 0
                self._next_translation = (0, 0, 0)

                # Dispatch to polymer() below — visits polymer_unit
                # directly into the shared self.state, and extracts this
                # block's own repeat_count via _parse_polymer_suffix.
                self.visit(child)

                atom_end = len(parent_mol.atoms) - 1
                block_repeat = parent_mol.repeat_count
                parent_mol.repeat_count = None  # consumed; avoid leaking
                                                  # into the next block or
                                                  # the whole-molecule level

                block = PolymerBlock(unit=None,
                                     repeat_count=block_repeat,
                                     block_kind=pending_kind,
                                     atom_start=atom_start,
                                     atom_end=atom_end)
                parent_mol.polymer_blocks.append(block)
                if parent_mol.block_topology is None and pending_kind:
                    parent_mol.block_topology = pending_kind
                pending_kind = ''  # consumed

    def bond(self, tree):
        self._next_bond_order, self._next_bond_dir, self._next_hapticity, self._next_bond_class = self._get_bond_info(tree)
        self._next_translation = (0, 0, 0)   # reset on plain bond visit

    # =============================================================
    # V4 Lattice Extension interpreter methods
    # =============================================================
    def lattice_cell(self, tree):
        """[V4 L3] Parse [Lattice:T] and/or [Thickness:Cls(args)] then visit the chain.

        The tags are sticky — they set self._pending_lattice_type /
        self._pending_thickness so that when the inner molecular_chain
        or polymer is visited (which calls self.state = GenerativeStateMachine()),
        the new mol gets tagged before it leaves script().

        [V4.1] Also handles generalized typed_tag children like
        [Mesh:Icosphere(2)] or [Material:Steel]. These are stashed on
        self._pending_typed_tags and applied to the mol as a list.
        """
        for child in tree.children:
            if not isinstance(child, Tree): continue
            data = child.data.lstrip('!')
            if data == 'lattice_tag':
                self._visit_lattice_tag(child)
            elif data == 'thickness_tag':
                self._visit_thickness_tag(child)
            elif data == 'typed_tag':
                # [V4.1] Generalized typed tag
                self._visit_typed_tag(child)
            elif data == 'lattice_subject':
                # [V4 L3] lattice_subject wraps molecular_chain / polymer /
                # polymer_block / peptide_chain. Visit its single child.
                for sub in child.children:
                    if not isinstance(sub, Tree): continue
                    sub_data = sub.data.lstrip('!')
                    # Set state machine with pending tags applied to its mol
                    self.state = GenerativeStateMachine()
                    if self._pending_lattice_type:
                        self.state.mol.lattice_type = self._pending_lattice_type
                        self._pending_lattice_type = None
                    if self._pending_thickness is not None:
                        cls, args = self._pending_thickness
                        self.state.mol.thickness_class = cls
                        self.state.mol.thickness_args = args
                        self._pending_thickness = None
                    if self._pending_typed_tags:
                        self.state.mol.typed_tags.extend(self._pending_typed_tags)
                        self._pending_typed_tags = []
                    # Visit the chain into this state
                    self.visit(sub)
                    self.state.finalize_valences()

    def _visit_typed_tag(self, tree):
        """[V4.1] Parse a generalized typed tag [Namespace:Value(args)].

        Stash as (namespace, value, args_tuple) on _pending_typed_tags.
        [V4.3] Now handles STRING, signed_num, and ; separators.
        """
        ns = None; val = None; args = []
        for child in tree.children:
            if isinstance(child, Token):
                if child.type == 'TAG_NAMESPACE': ns = str(child)
                elif child.type == 'TAG_VALUE':   val = str(child)
                elif child.type == 'STRING' and val is None:
                    val = str(child).strip('"')   # [V4.3] STRING as tag value
            elif isinstance(child, Tree) and child.data.lstrip('!') == 'tag_args':
                # Collect args — positional (FLOAT/INT/IDENTIFIER/STRING/signed_num) and kv (k=v)
                for sub in child.children:
                    if not isinstance(sub, Tree): continue
                    sd = sub.data.lstrip('!')
                    if sd == 'tag_arg':
                        for leaf in sub.children:
                            if isinstance(leaf, Tree) and leaf.data.lstrip('!') == 'tag_kv':
                                k = None; v = None
                                for t in leaf.scan_values(lambda x: isinstance(x, Token)):
                                    if t.type == 'IDENTIFIER' and k is None:
                                        k = str(t)
                                    elif t.type in ('FLOAT', 'INT'):
                                        v = float(str(t)) if t.type == 'FLOAT' else int(str(t))
                                    elif t.type == 'IDENTIFIER':
                                        v = str(t)
                                    elif t.type == 'STRING':
                                        v = str(t).strip('"')
                                if k: args.append((k, v))
                            elif isinstance(leaf, Tree) and leaf.data.lstrip('!') in ('pos_num', 'neg_num', 'num'):
                                # signed_num in tag_arg
                                val_num = self._parse_signed_num(leaf)
                                args.append(val_num)
                            elif isinstance(leaf, Token):
                                if leaf.type in ('FLOAT', 'INT'):
                                    args.append(float(str(leaf)) if leaf.type == 'FLOAT' else int(str(leaf)))
                                elif leaf.type == 'IDENTIFIER':
                                    args.append(str(leaf))
                                elif leaf.type == 'STRING':
                                    args.append(str(leaf).strip('"'))
        if ns and val:
            self._pending_typed_tags.append((ns, val, tuple(args)))

    def _visit_lattice_tag(self, tree):
        """Extract LATTICE_TYPE token and stash on _pending_lattice_type."""
        for t in tree.scan_values(lambda x: isinstance(x, Token)):
            if t.type == 'LATTICE_TYPE':
                self._pending_lattice_type = str(t)
                return
        self._pending_lattice_type = 'Custom'

    def _visit_thickness_tag(self, tree):
        """Extract THICKNESS_CLASS + optional args. Stash as (cls, args_tuple)."""
        cls = None
        args = []
        for child in tree.children:
            if isinstance(child, Token) and child.type == 'THICKNESS_CLASS':
                cls = str(child)
            elif isinstance(child, Tree) and child.data.lstrip('!') == 'thickness_args':
                # Collect FLOAT / INT / IDENTIFIER tokens in order
                for sub in child.iter_subtrees():
                    for t in sub.children if hasattr(sub, 'children') else []:
                        pass
                # Simpler: scan all tokens in thickness_args
                for t in child.scan_values(lambda x: isinstance(x, Token)):
                    s = str(t)
                    if t.type in ('FLOAT', 'INT'):
                        try:
                            args.append(float(s) if t.type == 'FLOAT' else int(s))
                        except: pass
                    elif t.type == 'IDENTIFIER':
                        args.append(s)
        self._pending_thickness = (cls or 'Constant', tuple(args))

    def post_process(self, tree):
        """[V4 L5] Parse one |> op(args) post-process operation.

        Stash on self._pending_post_process so macroscopic_structure can
        attach the ordered list to the final CoreMolecule.
        """
        op_name = None
        args = []
        for child in tree.children:
            if isinstance(child, Token) and child.type == 'PIPE_OP':
                op_name = str(child)
            elif isinstance(child, Tree) and child.data.lstrip('!') == 'pipe_args':
                # Collect positional + kv args
                for sub in child.children:
                    if not isinstance(sub, Tree): continue
                    sd = sub.data.lstrip('!')
                    if sd == 'pipe_arg':
                        # pipe_arg: pipe_kv | FLOAT | INT
                        for leaf in sub.children:
                            if isinstance(leaf, Tree) and leaf.data.lstrip('!') == 'pipe_kv':
                                k = None; v = None
                                for t in leaf.scan_values(lambda x: isinstance(x, Token)):
                                    if t.type == 'IDENTIFIER' and k is None:
                                        k = str(t)
                                    elif t.type in ('FLOAT', 'INT'):
                                        v = float(str(t)) if t.type == 'FLOAT' else int(str(t))
                                if k: args.append((k, v))
                            elif isinstance(leaf, Token) and leaf.type in ('FLOAT', 'INT'):
                                args.append(float(str(leaf)) if leaf.type == 'FLOAT' else int(str(leaf)))
        if op_name:
            self._pending_post_process.append((op_name, tuple(args)))

    def maybe_periodic_bond(self, tree):
        """
        V3.6: bond optionally followed by @tx,ty,tz translation vector.
        Children: [bond_tree] or [bond_tree, signed_num, signed_num, signed_num]

        [V4.2 Q2] signed_num now accepts FLOAT, so fractional translations
        like @0.5,0.5,0.5 work for crystallographic coordinates.
        """
        from lark import Tree as LarkTree
        children = [c for c in tree.children if c is not None]
        # First child is always the bond subtree
        bond_tree = children[0]
        self._next_bond_order, self._next_bond_dir, self._next_hapticity, self._next_bond_class = self._get_bond_info(bond_tree)
        # If there are signed_num children, extract translation
        sint_nodes = [c for c in children[1:] if isinstance(c, LarkTree)]
        if len(sint_nodes) == 3:
            tx = self._parse_signed_num(sint_nodes[0])
            ty = self._parse_signed_num(sint_nodes[1])
            tz = self._parse_signed_num(sint_nodes[2])
            self._next_translation = (tx, ty, tz)
        else:
            self._next_translation = (0, 0, 0)

    def neg_int(self, tree):
        return -int(str(tree.children[0]))

    def pos_int(self, tree):
        return int(str(tree.children[0]))

    # [V4.2 Q2] num / neg_num / pos_num handlers for float translations
    def neg_num(self, tree):
        return -self._parse_num(tree.children[-1])

    def pos_num(self, tree):
        return self._parse_num(tree.children[0])

    def num(self, tree):
        return self._parse_num(tree.children[0])

    def _parse_num(self, node) -> float:
        """Parse a num token (INT or FLOAT) to Python float.
        Handles both Token (direct) and Tree (num rule wrapping a token).
        """
        from lark import Tree, Token
        if isinstance(node, Tree):
            # num rule wraps a FLOAT or INT token — extract it
            if node.data == 'num':
                for child in node.children:
                    if isinstance(child, Token):
                        return self._parse_num(child)
            return 0.0
        if isinstance(node, Token):
            s = str(node)
            if node.type == 'FLOAT':
                return float(s)
            elif node.type == 'INT':
                return float(s)
        try:
            return float(str(node))
        except (ValueError, TypeError):
            return 0.0

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

    def _parse_signed_num(self, node) -> float:
        """[V4.2 Q2] Parse a signed_num tree node to Python float.
        Accepts both int and float translations.
        """
        from lark import Tree, Token
        if isinstance(node, Tree):
            if node.data == 'neg_num':
                return -self._parse_num(node.children[-1])
            elif node.data == 'pos_num':
                return self._parse_num(node.children[0])
        if isinstance(node, Token):
            return self._parse_num(node)
        return 0.0

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

    def spline_explicit(self, tree):
        """[V4.2 Q5] Parse ~[x,y,z;x,y,z;...] explicit spline control points.

        Sets bond_class='spline' and stashes control points on
        self._next_control_points for the next bond construction.
        """
        from lark import Tree as LarkTree
        ctrl_points = []
        for child in tree.children:
            if isinstance(child, LarkTree) and child.data.lstrip('!') == 'ctrl_point':
                # ctrl_point: signed_num "," signed_num "," signed_num
                nums = []
                for sub in child.children:
                    if isinstance(sub, LarkTree):
                        val = self._parse_signed_num(sub)
                        nums.append(val)
                if len(nums) == 3:
                    ctrl_points.append(tuple(nums))
        self._next_bond_order = 1
        self._next_bond_class = "spline"
        self._next_control_points = ctrl_points if ctrl_points else None

    def _get_bond_info(self, bond_node):
        # Extract bond type, direction, and bond_class from a bond node.
        # Note: hapticity (eta-n) is now handled by hbond() directly.
        # [V4.2 Q5] Handle spline_explicit as a child Tree of bond
        from lark import Tree as LarkTree
        for child in bond_node.children:
            if isinstance(child, LarkTree) and child.data.lstrip('!') == 'spline_explicit':
                # Spline with explicit control points — dispatch to handler
                self.spline_explicit(child)
                return self._next_bond_order, 0, 0, self._next_bond_class

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
        elif t_type == 'SPLINE_BOND': order = 1; bond_class = "spline"   # [V4 L1]
        elif t_type == 'BRIDGE_BOND': order = 1; bond_class = "bridge"   # [V4.3] 3c2e

        return order, direction, hapticity, bond_class

    def atom_expr(self, tree):
        count = self._get_multiplier(tree)
        # count is None when there's no digit after the atom (the common case).
        # A non-None count means there IS a digit — it's either a multiplier
        # (C3 = three carbons) or a SMILES-style ring closure (C1...C1).
        # _is_legacy_ring_closure_candidate decides which interpretation applies.
        if count is None:
            count = 1
            legacy_ring_count = None
        else:
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
                    control_points=self._next_control_points,
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
            # SMILES-style ring closures use digit pairs (1...1, 2...2).
            # The digit is a register LABEL, not a ring size. Pass as string
            # so add_ring uses the named-register path (first call stores
            # the register, second call closes the ring). Passing as int
            # would hit the V2 back-count path and compute the wrong target.
            self.state.add_ring(str(legacy_ring_count), bond_order=1)
            
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
            beam_radius = None  # [V4.2] extract from state_block
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
                elif t == 'state_block':
                    # [V4.2] Extract beam_radius (and spin/occupancy/is_excited)
                    # from state_block inside bracket atoms
                    for sub in child.children:
                        if not isinstance(sub, Tree): continue
                        st = sub.data.lstrip('!')
                        sval = "".join([str(leaf) for leaf in sub.scan_values(lambda x: not isinstance(x, Tree))])
                        if sval.startswith('r:'):
                            try: beam_radius = float(sval[2:])
                            except: beam_radius = None

        self.state.add_atom(element, charge=charge, isotope=isotope,
                            hcount=hcount, chiral=chiral,
                            bond_order=self._next_bond_order,
                            bond_dir=self._next_bond_dir,
                            bond_class=self._next_bond_class,
                            is_bracket=True, mapping=mapping,
                            radical=radical,
                            translation=self._next_translation,
                            control_points=self._next_control_points,
                            beam_radius=beam_radius if 'beam_radius' in dir() else None)
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
        beam_radius = None   # [V4 L2]
        
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
                        elif sval.startswith('r:'):
                            # [V4 L2] beam radius ratio <r:0.5>
                            try: beam_radius = float(sval[2:])
                            except: beam_radius = None
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
                            translation=self._next_translation,
                            beam_radius=beam_radius,
                            control_points=self._next_control_points)
        self._next_translation = (0, 0, 0)
        self._next_control_points = None  # [V4.2 Q5] reset after consumption

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


    def _get_multiplier(self, tree):
        """Return the explicit multiplier INT, or None if no multiplier is present.

        Returning None (instead of the previous default 1) lets callers
        distinguish "no digit after this atom" from "digit 1 after this
        atom".  This is critical for SMILES-style ring closures: ``C1``
        means "C with ring-closure digit 1", while ``C`` (no digit) means
        just "C".  The previous default of 1 caused every last-in-chain
        atom to be misidentified as a ring-closure candidate.
        """
        for child in tree.children:
            if isinstance(child, Tree) and child.data.lstrip('!') == 'multiplier':
                for node in child.children:
                    if isinstance(node, Token) and node.type == 'INT':
                        return int(str(node))
        return None

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

    def _is_legacy_ring_closure_candidate(self, tree, count) -> bool:
        # Any atom with an explicit digit multiplier is a SMILES-style ring
        # closure (e.g. C1...C1, [C@@H]1CCCN1).  The digit is a register
        # label, not a ring size or atom count.
        #
        # The previous implementation required the atom to be the LAST child
        # of molecular_chain, which is wrong — ring closures appear in the
        # MIDDLE of chains (e.g. [C@@H]1 in proline is followed by CCCN1).
        # That caused the ring bond to be silently dropped, which in turn
        # caused the chiral center to see only 2 of its 4 neighbors.
        if count is None or count < 1:
            return False
        parent = getattr(tree, 'parent', None)
        if parent is None or parent.data.lstrip('!') != 'molecular_chain':
            return False
        return True

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
            # Legacy SMILES-style ring closures use digit pairs (1...1, 2...2).
            # The digit is a register LABEL, not a ring size. Pass as string
            # so add_ring uses the named-register path.
            self.state.add_ring(str(ring_num), bond_order=self._next_bond_order)

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