"""
SCRIPT Validator
Delegates to the production Lark parser for structural correctness.
The lightweight pre-checks (bracket balance, X-outside-brackets) serve
as a fast path that avoids spinning up the parser for obviously-broken inputs.
"""

import re
from typing import Optional, Set


class SCRIPTValidator:
    """Validator for SCRIPT strings.

    Strategy:
    1. Fast pre-checks (O(n) character scan) to reject obviously-broken strings.
    2. Full parse via SCRIPTParser as the authoritative structural check.
       The parser enforces valence (Sandhi), grammar, and ring closure rules
       -- there is no need to reimplement those here.
    """

    def __init__(self):
        self._parser = None  # lazy-initialised to avoid circular import at module load

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def is_valid(self, script_string: str) -> bool:
        """Return True if *script_string* is a valid SCRIPT expression."""
        if not script_string or not script_string.strip():
            return False

        # 1. Bracket / paren / brace balance (O(n), no parser needed)
        if not self._check_balanced(script_string):
            return False

        # 2. 'X' must appear inside square brackets (halogen rule)
        if 'X' in script_string and not self._x_in_brackets(script_string):
            return False

        # 3. Peptide block amino-acid codes
        if not self._check_peptide_codes(script_string):
            return False

        # 4. Full structural parse (authoritative)
        return self._parse_ok(script_string)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _parse_ok(self, script_string: str) -> bool:
        """Attempt a full parse; return True on success."""
        try:
            parser = self._get_parser()
            result = parser.parse(script_string)
            return bool(result.get("success", False))
        except Exception:
            return False

    def _get_parser(self):
        if self._parser is None:
            from .parser import SCRIPTParser
            self._parser = SCRIPTParser()
        return self._parser

    def _check_balanced(self, s: str) -> bool:
        """Return True iff all bracket/paren/brace pairs are balanced."""
        stack = []
        closing = {')': '(', ']': '[', '}': '{'}
        for ch in s:
            if ch in closing.values():
                stack.append(ch)
            elif ch in closing:
                if not stack or stack[-1] != closing[ch]:
                    return False
                stack.pop()
        return not stack

    def _x_in_brackets(self, s: str) -> bool:
        """Return True iff every 'X' occurrence is inside square brackets."""
        depth = 0
        for ch in s:
            if ch == '[':
                depth += 1
            elif ch == ']':
                depth -= 1
            elif ch == 'X' and depth == 0:
                return False
        return True

    def _check_peptide_codes(self, script_string: str) -> bool:
        """Validate single-letter codes inside {..} peptide blocks.

        Spline control-point blocks (~{...}) are skipped because they contain
        numeric/comma content, not amino acid codes.
        """
        valid_amino_acids: Set[str] = {
            'A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I',
            'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'
        }
        for m in re.finditer(r'\{([^}]*)\}', script_string):
            start = m.start()
            # Skip spline control-point blocks (preceded by '~')
            if start > 0 and script_string[start - 1] == '~':
                continue
            inner = m.group(1)
            for part in inner.split('.'):
                # Single-character codes must be valid amino acids;
                # multi-character codes (PTM codes, nucleotide codes) are
                # passed to the full parser for validation.
                if len(part) == 1 and part not in valid_amino_acids:
                    return False
        return True

    @staticmethod
    def _get_valid_elements() -> Set[str]:
        """Full set of IUPAC element symbols (Z=1..118)."""
        return {
            'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
            'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca',
            'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
            'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr',
            'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
            'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
            'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
            'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
            'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
            'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
            'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
            'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og',
        }


def is_valid_SCRIPT(script_string: str) -> bool:
    """Module-level convenience wrapper around SCRIPTValidator.is_valid."""
    return SCRIPTValidator().is_valid(script_string)