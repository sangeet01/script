"""
script/visualizer.py
────────────────────
Terminal-based molecular structure renderer for SCRIPT.

Converts an RDKit mol object into a Unicode line-art grid, with
atom symbols placed at computed 2D positions and bonds drawn
as ─ │ ╱ ╲ lines. Falls back to a text summary when RDKit is absent.
"""

import math
from typing import List, Tuple, Optional

# ── Atom colours (rich markup) ────────────────────────────────────────────────
ATOM_COLORS = {
    "C":  "bright_white",
    "N":  "bright_blue",
    "O":  "bright_red",
    "S":  "bright_yellow",
    "P":  "orange1",
    "F":  "bright_green",
    "Cl": "green",
    "Br": "dark_orange",
    "I":  "magenta",
    "H":  "grey50",
    "Na": "cyan",
    "K":  "cyan",
    "Ca": "cyan",
    "Mg": "cyan",
    "Fe": "bright_red",
    "Cu": "bright_red",
    "Zn": "bright_red",
}
DEFAULT_ATOM_COLOR = "bright_white"

# ── Bond characters ───────────────────────────────────────────────────────────
# (dx, dy) → (single_char, double_char, aromatic_char)
BOND_CHARS = {
    "h":  ("─", "═", "~"),   # horizontal
    "v":  ("│", "║", "¦"),   # vertical
    "du": ("╱", "╱", "╱"),   # diagonal up-right
    "dd": ("╲", "╲", "╲"),   # diagonal down-right
}


def _bond_dir(dx: int, dy: int) -> str:
    """Classify bond direction into h / v / du / dd."""
    adx, ady = abs(dx), abs(dy)
    if adx == 0:
        return "v"
    if ady == 0:
        return "h"
    if (dx > 0 and dy < 0) or (dx < 0 and dy > 0):
        return "du"
    return "dd"


def _bresenham(x1: int, y1: int, x2: int, y2: int) -> List[Tuple[int, int]]:
    """Return all grid cells along a line (exclusive of endpoints)."""
    cells = []
    dx = abs(x2 - x1); dy = abs(y2 - y1)
    sx = 1 if x1 < x2 else -1
    sy = 1 if y1 < y2 else -1
    err = dx - dy
    cx, cy = x1, y1
    while True:
        if (cx, cy) != (x1, y1) and (cx, cy) != (x2, y2):
            cells.append((cx, cy))
        if cx == x2 and cy == y2:
            break
        e2 = 2 * err
        if e2 > -dy:
            err -= dy; cx += sx
        if e2 < dx:
            err += dx; cy += sy
    return cells


def _draw_bond(
    char_grid: List[List[str]],
    x1: int, y1: int, x2: int, y2: int,
    bond_order: float, H: int, W: int,
):
    """Draw one bond onto the character grid."""
    dx = x2 - x1; dy = y2 - y1
    d = _bond_dir(dx, dy)
    single, double, aromatic = BOND_CHARS[d]

    if bond_order == 1.5:
        char = aromatic
    elif bond_order >= 2:
        char = double
    else:
        char = single

    for cx, cy in _bresenham(x1, y1, x2, y2):
        if 0 <= cx < W and 0 <= cy < H:
            if char_grid[cy][cx] == " ":
                char_grid[cy][cx] = char

    # Second line for double bonds (perpendicular offset of 1)
    if bond_order >= 2 and bond_order != 1.5:
        perp_x = 0 if d in ("du", "dd") else (0 if d == "h" else 1)
        perp_y = 1 if d == "h" else (0 if d == "v" else 0)
        for cx, cy in _bresenham(x1, y1, x2, y2):
            nx, ny = cx + perp_x, cy + perp_y
            if 0 <= nx < W and 0 <= ny < H:
                if char_grid[ny][nx] == " ":
                    char_grid[ny][nx] = double


def mol_to_grid(
    mol,
    width: int = 68,
    height: int = 22,
) -> Tuple[List[List[str]], List[List[Optional[str]]]]:
    """
    Render an RDKit mol to (char_grid, color_grid).

    char_grid  – 2D list of single characters
    color_grid – 2D list of rich color strings or None
    """
    try:
        from rdkit.Chem import AllChem
        AllChem.Compute2DCoords(mol)
        conf = mol.GetConformer()
    except Exception:
        return _fallback_grid(mol, width, height)

    n = mol.GetNumAtoms()
    if n == 0:
        return _fallback_grid(mol, width, height)

    # Terminal characters are ~2.1× taller than wide — stretch x
    ASPECT = 2.1
    raw = [(conf.GetAtomPosition(i).x * ASPECT,
            conf.GetAtomPosition(i).y) for i in range(n)]

    xs = [p[0] for p in raw]
    ys = [p[1] for p in raw]
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)

    pad = 3
    rx = max_x - min_x or 1.0
    ry = max_y - min_y or 1.0
    scale = min((width - 2 * pad) / rx, (height - 2 * pad) / ry) * 0.85

    gpos = []
    for x, y in raw:
        gx = int((x - min_x) * scale) + pad
        gy = height - 1 - int((y - min_y) * scale) - pad
        gpos.append((gx, gy))

    char_grid  = [[" "] * width for _ in range(height)]
    color_grid = [[None] * width for _ in range(height)]

    # Bonds first (atoms drawn on top)
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        _draw_bond(
            char_grid,
            gpos[i][0], gpos[i][1],
            gpos[j][0], gpos[j][1],
            bond.GetBondTypeAsDouble(),
            height, width,
        )

    # Atoms on top
    for i, atom in enumerate(mol.GetAtoms()):
        x, y = gpos[i]
        sym = atom.GetSymbol()
        col = ATOM_COLORS.get(sym, DEFAULT_ATOM_COLOR)
        for k, ch in enumerate(sym):
            nx = x + k
            if 0 <= y < height and 0 <= nx < width:
                char_grid[y][nx]  = ch
                color_grid[y][nx] = col

    return char_grid, color_grid


def _fallback_grid(mol, width, height):
    """Plain text grid when 2D coords are unavailable."""
    char_grid  = [[" "] * width for _ in range(height)]
    color_grid = [[None]  * width for _ in range(height)]
    msg = f"[{mol.GetNumAtoms()} atoms, {mol.GetNumBonds()} bonds]"
    cy = height // 2
    cx = (width - len(msg)) // 2
    for k, ch in enumerate(msg):
        if 0 <= cx + k < width:
            char_grid[cy][cx + k] = ch
    return char_grid, color_grid


def grid_to_rich(
    char_grid: List[List[str]],
    color_grid: List[List[Optional[str]]],
) -> str:
    """
    Convert (char_grid, color_grid) to a rich-markup string.
    Atoms are coloured; bond lines are dim grey.
    """
    lines = []
    for row_c, row_col in zip(char_grid, color_grid):
        line = ""
        i = 0
        while i < len(row_c):
            ch = row_c[i]
            col = row_col[i]
            if ch == " ":
                line += " "; i += 1
            elif col is not None:
                # Collect full atom symbol (may be 2 chars)
                sym = ch
                if i + 1 < len(row_c) and row_col[i + 1] == col:
                    sym += row_c[i + 1]; i += 1
                line += f"[bold {col}]{sym}[/bold {col}]"
                i += 1
            else:
                line += f"[dim]{ch}[/dim]"; i += 1
        lines.append(line)
    # Strip trailing blank lines
    while lines and lines[-1].strip() == "":
        lines.pop()
    return "\n".join(lines)


def render_mol(mol, width: int = 68, height: int = 22) -> str:
    """Full pipeline: mol → rich-markup string."""
    g, c = mol_to_grid(mol, width, height)
    return grid_to_rich(g, c)


def mol_properties(mol) -> dict:
    """Return a dict of basic molecular properties."""
    try:
        from rdkit.Chem import Descriptors, rdMolDescriptors
        from rdkit import Chem
        formula = rdMolDescriptors.CalcMolFormula(mol)
        mw      = round(Descriptors.ExactMolWt(mol), 3)
        n_heavy = mol.GetNumHeavyAtoms()
        n_atoms = mol.GetNumAtoms()
        n_bonds = mol.GetNumBonds()
        n_rings = rdMolDescriptors.CalcNumRings(mol)
        n_arom  = rdMolDescriptors.CalcNumAromaticRings(mol)
        n_rot   = rdMolDescriptors.CalcNumRotatableBonds(mol)
        n_stereo= rdMolDescriptors.CalcNumAtomStereoCenters(mol)
        n_hbd   = rdMolDescriptors.CalcNumHBD(mol)
        n_hba   = rdMolDescriptors.CalcNumHBA(mol)
        tpsa    = round(Descriptors.TPSA(mol), 1)
        return {
            "Formula":      formula,
            "Exact MW":     f"{mw} Da",
            "Heavy atoms":  str(n_heavy),
            "Bonds":        str(n_bonds),
            "Rings":        str(n_rings),
            "Aromatic rings": str(n_arom),
            "Rot. bonds":   str(n_rot),
            "Stereocenters":str(n_stereo),
            "HBD / HBA":    f"{n_hbd} / {n_hba}",
            "TPSA":         f"{tpsa} Å²",
        }
    except Exception:
        return {
            "Atoms": str(mol.GetNumAtoms()),
            "Bonds": str(mol.GetNumBonds()),
        }
