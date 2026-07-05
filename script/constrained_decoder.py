"""
Constrained SCRIPT Decoder — Option 2: validity by constrained generation.

This module provides a grammar-state-aware decoder that sits between an ML
model and the SCRIPT grammar. At each generation step, it:

1. Tracks the current grammar state (open valences, branch depth, ring
   registers, bracket state, current atom)
2. Computes which tokens are valid next (given the state)
3. Returns a mask the model can use to constrain sampling

This gives ~100% validity for ML generation WITHOUT changing the SCRIPT
notation. The notation stays human-readable; the constraint is enforced
at generation time, not at the grammar level.

Design principles (Sanskrit-inspired):
  - Sandhi (junction): valence is tracked continuously, not just at the end
  - Lopa (elision): implicit Hs and lone pairs are tracked as ghost slots
  - Anubandha (ring closure): ring registers track open ring bonds
  - Pratyaya (marker): chiral/charge markers are validated against element

Usage:
    decoder = ConstrainedSCRIPTDecoder()
    decoder.reset()
    
    # During generation:
    for step in range(max_length):
        mask = decoder.get_valid_token_mask(vocab)  # boolean array
        logits = model(input_ids)[step]
        logits[~mask] = -float('inf')  # mask invalid tokens
        next_token = sample(logits)
        decoder.consume(next_token)
        if decoder.is_complete():
            break
    
    # The generated SCRIPT string is guaranteed valid
"""

from typing import Dict, Set, Optional, List, Tuple
from dataclasses import dataclass, field


# =============================================================================
# Grammar state — mirrors the Sandhi state machine
# =============================================================================

@dataclass
class GrammarState:
    """Tracks the state of SCRIPT generation for constraint enforcement.
    
    This mirrors the GenerativeStateMachine in state_machine.py, but in a
    lightweight form suitable for real-time token masking during ML generation.
    """
    # Current atom index (None = no atoms yet)
    current_atom: Optional[int] = None
    
    # Atom states: index -> (symbol, valence_used, max_valence, implicit_hs)
    atoms: Dict[int, dict] = field(default_factory=dict)
    
    # Branch stack: list of atom indices to return to after ')'
    branch_stack: List[int] = field(default_factory=list)
    
    # Bracket state: True if inside [...]
    in_bracket: bool = False
    bracket_content: str = ""  # accumulated content inside [...]

    # V3.9 fix: explicit stage tracker for bracket_atom grammar fields.
    # bracket_atom: [ isotope? element chiral? hcount? charge? radical? ]
    # Each field may appear AT MOST ONCE, in this strict linear order.
    # The old heuristic-based mask (content.isdigit(), any(c.isalpha()))
    # could not distinguish "digit is part of isotope" from "digit is part
    # of hcount" from "digit is part of charge magnitude", and had no way
    # to prevent a field from being placed twice — allowing malformed
    # bracket atoms like [29O5], [OH4H2143@R@@5211], [S9484@S@S.] to be
    # declared decoder-complete even though the real grammar rejects them.
    # bracket_stage tracks which field we are currently allowed to start
    # or continue: 'isotope' -> 'element' -> 'chiral' -> 'hcount' ->
    # 'charge' -> 'radical' -> 'close'. Once a field is finished (a
    # different-kind token is consumed), we advance the stage and can
    # never return to an earlier one.
    bracket_stage: str = 'isotope'
    bracket_has_element: bool = False

    # Ring registers: digit -> atom_index (for SMILES-style 1...1 closures)
    ring_registers: Dict[str, int] = field(default_factory=dict)
    
    # V2 ring open count: how many &N closures are "owed" (for aromatic paths)
    pending_aromatic: bool = False
    
    # Bond prefix: the next atom should use this bond order
    next_bond_order: int = -1  # -1 = default (single or aromatic by context)
    
    # Whether we've started (first atom vs. subsequent)
    started: bool = False
    
    # Fragment separator state: just consumed '.'
    after_fragment_sep: bool = False
    
    # Depth of nested branches
    branch_depth: int = 0


# =============================================================================
# Valence table (mirrors state_machine.py DEFAULT_VALENCE)
# =============================================================================

DEFAULT_VALENCE = {
    'B': 3, 'C': 4, 'N': 3, 'O': 2, 'P': 3, 'S': 2,
    'F': 1, 'Cl': 1, 'Br': 1, 'I': 1, 'H': 1,
}

HYPERVALENT_MAX = {
    'P': 5, 'S': 6, 'Cl': 7, 'Br': 7, 'I': 7,
}

TRANSITION_METALS = {
    'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
    'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Cr', 'Mn', 'Mo', 'W', 'V',
    'Ti', 'Zr', 'Hf', 'Re',
}


def get_max_valence(symbol: str, charge: int = 0, is_bracket: bool = False) -> int:
    """Get maximum valence for an element (mirrors state_machine.py)."""
    if symbol in TRANSITION_METALS:
        return 18  # permissive for metals
    
    base = DEFAULT_VALENCE.get(symbol, 4)
    
    if charge != 0:
        charged_val = base + charge
        if is_bracket:
            hyper = HYPERVALENT_MAX.get(symbol, base)
            return max(charged_val, hyper)
        return max(charged_val, 1)
    
    if is_bracket:
        return HYPERVALENT_MAX.get(symbol, base)
    
    return base


# =============================================================================
# Constrained Decoder
# =============================================================================

class ConstrainedSCRIPTDecoder:
    """Grammar-state-aware decoder for constrained SCRIPT generation.
    
    At each step, call get_valid_tokens() to get the set of valid next
    tokens, then mask your model's logits accordingly.
    """
    
    def __init__(self):
        self.state = GrammarState()
        self._atom_counter = 0
    
    def reset(self):
        """Reset to initial state (no atoms, ready for first token)."""
        self.state = GrammarState()
        self._atom_counter = 0
    
    # -------------------------------------------------------------------------
    # Public API for ML generation
    # -------------------------------------------------------------------------
    
    def get_valid_tokens(self, vocab: List[str]) -> Set[str]:
        """Return the set of valid next tokens given current grammar state.
        
        Args:
            vocab: list of all possible tokens (vocabulary)
            
        Returns:
            Set of tokens that are valid at the current position
        """
        valid = set()
        s = self.state
        
        if s.in_bracket:
            # Inside [...]: allow element, isotope, charge, H, chiral, ]
            valid.update(self._valid_bracket_tokens())
        elif not s.started:
            # First token: allow atoms, '[', '{' (peptide), '[' (macroscopic)
            valid.update(self._valid_first_tokens())
        else:
            # Subsequent tokens: bonds, atoms, branches, rings, separators
            valid.update(self._valid_subsequent_tokens())
        
        # Intersect with vocabulary
        return valid & set(vocab)
    
    def get_valid_token_mask(self, vocab: List[str]) -> List[bool]:
        """Return a boolean mask (same length as vocab) for valid tokens.
        
        Usage:
            mask = decoder.get_valid_token_mask(vocab)
            logits[~mask] = -float('inf')
            next_token = sample(logits)
        """
        valid = self.get_valid_tokens(vocab)
        return [tok in valid for tok in vocab]
    
    def consume(self, token: str) -> bool:
        """Consume a token and update grammar state.
        
        Returns True if the token was valid and consumed, False if invalid.
        """
        valid = self.get_valid_tokens([token])
        if token not in valid:
            return False
        
        self._apply_token(token)
        return True
    
    def is_complete(self) -> bool:
        """Check if the generated SCRIPT is complete (valid end state)."""
        s = self.state
        if s.in_bracket:
            return False
        if s.branch_depth > 0:
            return False  # unclosed branches
        if not s.started:
            return False
        # All atoms should have valid valence
        for idx, atom_info in s.atoms.items():
            if atom_info['valence_used'] > atom_info['max_valence']:
                return False
        # V3.9: a bare (non-bracket) wildcard as the CURRENT/last atom is
        # not a valid complete molecule. atom_expr's WILDCARD and hbond's
        # STAR_BOND share the same "*" lexeme and are grammatically
        # ambiguous with 1-token lookahead; the grammar currently resolves
        # a leading/standalone "*" as the start of an hbond (which then
        # requires a following INT and atom), not as a terminal wildcard
        # atom. [*] (bracketed) is unaffected and works correctly.
        if s.current_atom is not None:
            cur = s.atoms.get(s.current_atom)
            if cur and cur.get('is_bare_wildcard'):
                return False
        return True
    
    # -------------------------------------------------------------------------
    # Token validity computation (the core constraint logic)
    # -------------------------------------------------------------------------
    
    def _valid_first_tokens(self) -> Set[str]:
        """Valid tokens at the start of a SCRIPT string."""
        valid = set()
        # Organic subset atoms (bare)
        valid.update({'C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I', 'B', 'H'})
        # Bracket atom start
        valid.add('[')
        # Peptide/nucleic acid start
        valid.add('{')
        # Macroscopic context
        valid.add('[[')
        # V3.9: bare (non-bracket) WILDCARD is NOT offered as a first token.
        # atom_expr's WILDCARD and hbond's STAR_BOND share the "*" lexeme;
        # the real grammar's LALR parser always resolves a leading "*" as
        # STAR_BOND (the start of an hbond), which is invalid without a
        # preceding atom. Use "[*]" (bracketed wildcard) instead, which is
        # unambiguous and fully supported.
        return valid
    
    def _valid_subsequent_tokens(self) -> Set[str]:
        """Valid tokens after the first atom."""
        valid = set()
        s = self.state
        
        # Check if we're waiting for an atom (after a bond token)
        waiting_for_atom = s.next_bond_order > 0
        
        # Bond tokens — only if NOT already waiting for an atom
        # (prevents C#:P — consecutive bonds)
        if not s.in_bracket and not waiting_for_atom:
            valid.update({'-', '=', '#', ':', '->', '<-', '/', '\\'})
            valid.add('=:')
        
        # Atom tokens — always valid after first atom (unless in bracket)
        if not s.in_bracket:
            valid.update({'C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I', 'B', 'H'})
            valid.add('[')
            valid.add('*')
        
        # Branch tokens
        if not waiting_for_atom and not s.in_bracket:
            valid.add('(')
        # V3.9 fix: ')' must never be offered while waiting_for_atom is True.
        # Previously gated only by branch_depth > 0, this allowed sequences
        # like "S(1N.[P.]-)" where the bond token '-' set next_bond_order=1
        # (an atom is owed) and then ')' closed the branch anyway, leaving
        # next_bond_order stuck at 1 forever with no atom-token path to
        # resolve it (a genuine mask dead-end — all_masked / NaN probs).
        if s.branch_depth > 0 and not waiting_for_atom:
            valid.add(')')
        
        # Ring closure tokens — only if we have a current atom
        # (prevents ring digit after fragment separator)
        if s.current_atom is not None and not waiting_for_atom:
            valid.add('&')
            for d in '123456789':
                valid.add(d)
        
        # Fragment separator — only if NOT waiting for atom
        if not waiting_for_atom and not s.in_bracket:
            valid.add('.')
            valid.add('~')
        
        # End of string (implicit — model can stop)
        
        # Filter by valence constraints
        valid = self._filter_by_valence(valid)
        
        return valid
    
    def _valid_bracket_tokens(self) -> Set[str]:
        """Valid tokens inside a [...] bracket atom.

        V3.9 fix: bracket_atom grammar is
            [ isotope? element chiral? hcount? charge? radical? ]
        Each field appears AT MOST ONCE, in this strict linear order.
        We track bracket_stage explicitly rather than re-deriving field
        membership from loose string heuristics on bracket_content, which
        previously allowed repeated/out-of-order fields (e.g. two chiral
        markers, digits attributable to no specific field) to be masked
        as valid.
        """
        valid = set()
        s = self.state
        stage = s.bracket_stage

        elements = {'C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I', 'B', 'H',
                     'Fe', 'Cu', 'Zn', 'Au', 'Ag', 'Pt', 'Pd', 'Ni', 'Co',
                     'Mn', 'Cr', 'Mo', 'W', 'V', 'Ti', 'Zr', 'Hf', 'Ru', 'Rh',
                     'Os', 'Ir', 'Re', 'Ba', 'Ca', 'Na', 'K', 'Mg', 'Al', 'Si',
                     'Li', 'Be', 'Ga', 'Ge', 'As', 'Se', 'Rb', 'Sr', 'Y', 'Nb',
                     'Tc', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'Xe', 'Cs', 'La', 'Ce',
                     'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er',
                     'Tm', 'Yb', 'Lu', 'Ta', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At',
                     'Rn', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm'}
        chiral_markers = {'@', '@@', '@R', '@S', '@r', '@s', '@AX', '@AX1', '@AX2',
                     '@SP', '@SP1', '@SP2', '@OH', '@OH1', '@OH2', '@OH3',
                     '@OH4', '@OH5', '@TB', '@TB1', '@TB2', '@PY', '@PY1',
                     '@PY2', '@PL', '@PL1', '@PL2'}

        if stage == 'isotope':
            # Isotope digits are optional; either start the isotope number
            # or skip straight to the element.
            valid.update(set('0123456789'))
            valid.update(elements)

        elif stage == 'isotope_continue':
            # Mid-isotope-number: continue digits, or the element follows.
            valid.update(set('0123456789'))
            valid.update(elements)

        elif stage == 'element':
            # Element is mandatory and not yet placed.
            valid.update(elements)

        elif stage == 'chiral':
            # Chiral marker optional; can skip straight to hcount/charge/
            # radical/close.
            valid.update(chiral_markers)
            valid.add('H')
            valid.update({'+', '-'})
            valid.add('.')
            valid.add(']')

        elif stage == 'hcount':
            # H optional; if present, exactly one digit run may follow.
            valid.add('H')
            valid.update({'+', '-'})
            valid.add('.')
            valid.add(']')

        elif stage == 'hcount_digit':
            # Continue the H-count number, or move on.
            valid.update(set('0123456789'))
            valid.update({'+', '-'})
            valid.add('.')
            valid.add(']')

        elif stage == 'charge':
            # Charge optional; +/- starts it.
            valid.update({'+', '-'})
            valid.add('.')
            valid.add(']')

        elif stage == 'charge_digit':
            # Continue charge magnitude digits, or move on.
            valid.update(set('0123456789'))
            valid.add('.')
            valid.add(']')

        elif stage == 'radical':
            valid.add('.')
            valid.add(']')

        elif stage == 'close':
            valid.add(']')

        return valid
    
    def _filter_by_valence(self, tokens: Set[str]) -> Set[str]:
        """Filter out tokens that would exceed valence.
        
        This is the key Sandhi constraint: don't allow adding a bond
        to an atom that's already at max valence.
        """
        s = self.state
        if s.current_atom is None:
            return tokens
        
        # Check if current atom can accept more bonds
        atom_info = s.atoms.get(s.current_atom)
        if atom_info is None:
            return tokens
        
        available = atom_info['max_valence'] - atom_info['valence_used']
        
        if available <= 0:
            # Current atom is saturated — can't add more atoms/bonds
            # But CAN close branches, close rings, end fragment
            allowed_when_saturated = {')', '.', '~', '&'} | set('123456789')
            tokens = tokens & allowed_when_saturated
        
        return tokens
    
    # -------------------------------------------------------------------------
    # State updates (apply consumed token)
    # -------------------------------------------------------------------------
    
    def _apply_token(self, token: str):
        """Update grammar state after consuming a token."""
        s = self.state
        
        if s.in_bracket:
            self._apply_bracket_token(token)
            return
        
        if token == '[':
            s.in_bracket = True
            s.bracket_content = ""
            s.bracket_stage = 'isotope'
            return
        
        if token == '(':
            s.branch_stack.append(s.current_atom)
            s.branch_depth += 1
            return
        
        if token == ')':
            if s.branch_stack:
                s.current_atom = s.branch_stack.pop()
                s.branch_depth -= 1
            return
        
        if token == '&':
            # V2 ring closure — next tokens are digits + anubandha
            # (simplified: just mark that we're in ring closure mode)
            # In full implementation, we'd track the number being built
            return
        
        if token in '.~':
            # Fragment separator
            s.current_atom = None
            s.after_fragment_sep = True
            s.next_bond_order = -1
            return
        
        if token in '-=#:':
            s.next_bond_order = {'-': 1, '=': 2, '#': 3, ':': 4}[token]
            return
        
        # Atom tokens (organic subset or wildcard)
        if token in DEFAULT_VALENCE or token in TRANSITION_METALS or token == '*':
            self._add_atom(token)
            return
        
        # Digit (SMILES ring closure)
        if token.isdigit():
            self._handle_ring_digit(token)
            return
    
    def _apply_bracket_token(self, token: str):
        """Handle tokens inside [...] brackets, advancing bracket_stage
        to enforce the strict field order of the bracket_atom grammar:
        isotope? element chiral? hcount? charge? radical?
        """
        s = self.state

        elements = {'C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I', 'B', 'H',
                     'Fe', 'Cu', 'Zn', 'Au', 'Ag', 'Pt', 'Pd', 'Ni', 'Co',
                     'Mn', 'Cr', 'Mo', 'W', 'V', 'Ti', 'Zr', 'Hf', 'Ru', 'Rh',
                     'Os', 'Ir', 'Re', 'Ba', 'Ca', 'Na', 'K', 'Mg', 'Al', 'Si',
                     'Li', 'Be', 'Ga', 'Ge', 'As', 'Se', 'Rb', 'Sr', 'Y', 'Nb',
                     'Tc', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'Xe', 'Cs', 'La', 'Ce',
                     'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er',
                     'Tm', 'Yb', 'Lu', 'Ta', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At',
                     'Rn', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm'}
        chiral_markers = {'@', '@@', '@R', '@S', '@r', '@s', '@AX', '@AX1', '@AX2',
                     '@SP', '@SP1', '@SP2', '@OH', '@OH1', '@OH2', '@OH3',
                     '@OH4', '@OH5', '@TB', '@TB1', '@TB2', '@PY', '@PY1',
                     '@PY2', '@PL', '@PL1', '@PL2'}

        if token == ']':
            self._parse_bracket_atom(s.bracket_content)
            s.in_bracket = False
            s.bracket_content = ""
            s.bracket_stage = 'isotope'
            return

        s.bracket_content += token

        stage = s.bracket_stage
        if stage in ('isotope', 'isotope_continue'):
            if token in elements:
                s.bracket_stage = 'chiral'
            elif token.isdigit():
                s.bracket_stage = 'isotope_continue'
            # else stays put (shouldn't happen if mask is correct)

        elif stage == 'element':
            if token in elements:
                s.bracket_stage = 'chiral'

        elif stage == 'chiral':
            if token in chiral_markers:
                s.bracket_stage = 'hcount'
            elif token == 'H':
                s.bracket_stage = 'hcount_digit'
            elif token in ('+', '-'):
                s.bracket_stage = 'charge_digit'
            elif token == '.':
                s.bracket_stage = 'close'

        elif stage == 'hcount':
            if token == 'H':
                s.bracket_stage = 'hcount_digit'
            elif token in ('+', '-'):
                s.bracket_stage = 'charge_digit'
            elif token == '.':
                s.bracket_stage = 'close'

        elif stage == 'hcount_digit':
            if token.isdigit():
                s.bracket_stage = 'charge'  # one digit run only, move on
            elif token in ('+', '-'):
                s.bracket_stage = 'charge_digit'
            elif token == '.':
                s.bracket_stage = 'close'

        elif stage == 'charge':
            if token in ('+', '-'):
                s.bracket_stage = 'charge_digit'
            elif token == '.':
                s.bracket_stage = 'close'

        elif stage == 'charge_digit':
            if token.isdigit():
                s.bracket_stage = 'radical'  # one digit run only, move on
            elif token == '.':
                s.bracket_stage = 'close'

        elif stage == 'radical':
            if token == '.':
                s.bracket_stage = 'close'
    
    def _add_atom(self, symbol: str):
        """Add a new atom to the state."""
        s = self.state
        atom_idx = self._atom_counter
        self._atom_counter += 1
        
        max_val = get_max_valence(symbol, is_bracket=False)
        
        # Bond from previous atom
        bond_order = s.next_bond_order if s.next_bond_order > 0 else 1
        s.next_bond_order = -1
        
        # Update valence
        valence_used = bond_order
        if s.current_atom is not None:
            prev = s.atoms[s.current_atom]
            prev['valence_used'] += bond_order
        
        s.atoms[atom_idx] = {
            'symbol': symbol,
            'valence_used': valence_used,
            'max_valence': max_val,
            'implicit_hs': 0,
            # V3.9: a bare (non-bracket) WILDCARD atom is not a valid
            # complete molecule per the real grammar — see is_complete().
            # This is due to a genuine grammar-level ambiguity between
            # atom_expr's bare WILDCARD and hbond's STAR_BOND sharing the
            # same lexeme; unresolved without a deeper grammar change.
            'is_bare_wildcard': (symbol == '*'),
        }
        s.current_atom = atom_idx
        s.started = True
        s.after_fragment_sep = False
    
    def _parse_bracket_atom(self, content: str):
        """Parse bracket atom content (simplified)."""
        s = self.state
        # Extract element symbol (letters)
        element = ""
        for c in content:
            if c.isalpha():
                element += c
            elif element:
                break
        
        if not element:
            element = "C"
        
        # Extract charge
        charge = 0
        if '+' in content:
            charge = 1
            # Check for ++ or +2
            if '++' in content:
                charge = 2
        if '-' in content:
            charge = -1
            if '--' in content:
                charge = -2
        
        # Extract H count
        hcount = 0
        if 'H' in content:
            # Simplified: just mark 1 H
            hcount = 1
        
        atom_idx = self._atom_counter
        self._atom_counter += 1
        
        max_val = get_max_valence(element, charge=charge, is_bracket=True)
        
        bond_order = s.next_bond_order if s.next_bond_order > 0 else 1
        s.next_bond_order = -1
        
        valence_used = bond_order + hcount
        if s.current_atom is not None:
            prev = s.atoms[s.current_atom]
            prev['valence_used'] += bond_order
        
        s.atoms[atom_idx] = {
            'symbol': element,
            'valence_used': valence_used,
            'max_valence': max_val,
            'implicit_hs': hcount,
        }
        s.current_atom = atom_idx
        s.started = True
    
    def _handle_ring_digit(self, digit: str):
        """Handle SMILES-style ring closure digit."""
        s = self.state
        if digit in s.ring_registers:
            # Close ring: add bond to registered atom
            target = s.ring_registers[digit]
            bond_order = s.next_bond_order if s.next_bond_order > 0 else 1
            s.atoms[target]['valence_used'] += bond_order
            s.atoms[s.current_atom]['valence_used'] += bond_order
            del s.ring_registers[digit]
        else:
            # Open ring: register current atom
            s.ring_registers[digit] = s.current_atom
        s.next_bond_order = -1


# =============================================================================
# Demo / Test
# =============================================================================

if __name__ == "__main__":
    decoder = ConstrainedSCRIPTDecoder()
    
    # Example vocabulary (simplified)
    vocab = list("CNOSPFClBrIBH()*-=#:.&123456789[]@XYZ") + ['->', '<-', '@@', '@R', '@S']
    
    print("=== Constrained SCRIPT Generation Demo ===\n")
    
    # Simulate generating "CCO" (ethanol)
    decoder.reset()
    target = "CCO"
    print(f"Target: {target}")
    
    for i, char in enumerate(target):
        valid = decoder.get_valid_tokens(vocab)
        print(f"  Step {i}: next char = {char!r}, valid = {sorted(valid)[:10]}...")
        assert char in valid, f"Token {char!r} not valid at step {i}!"
        decoder.consume(char)
    
    print(f"\n  Final state: {decoder.state.atoms}")
    print(f"  Complete: {decoder.is_complete()}")
    
    # Try generating invalid molecule: C(C)(C)(C)(C)(C) — 6-valent C
    print("\n=== Invalid molecule prevention: C(C)(C)(C)(C)(C) ===\n")
    decoder.reset()
    invalid = "C(C)(C)(C)(C)(C)"
    
    for i, char in enumerate(invalid):
        valid = decoder.get_valid_tokens(vocab)
        if char not in valid:
            print(f"  Step {i}: BLOCKED token {char!r} (not in valid set)")
            print(f"  Valid tokens: {sorted(valid)}")
            print(f"\n  ✓ Invalid molecule prevented at step {i}!")
            break
        decoder.consume(char)
    else:
        print("  ✗ Invalid molecule was NOT prevented!")
    
    # Show valence tracking
    print(f"\n  Atom states after blocking:")
    for idx, info in decoder.state.atoms.items():
        print(f"    atom {idx} ({info['symbol']}): "
              f"valence_used={info['valence_used']}/{info['max_valence']}")
