"""
script/cli.py
─────────────
SCRIPT command-line interface.

    scr "C:C:C:C:C:C&6:"           draw + info  (default)
    scr draw "CCCCCC&6-"           ASCII art
    scr canon "c1ccccc1" -f smiles canonicalize
    scr validate "CC(=O)O"         Sandhi check
    scr convert "CCO" -f smiles    SMILES → SCRIPT
    scr info "C:C:C:C:C:C&6:"      properties
    echo "CCO" | scr convert -f smiles   pipe support
"""

import sys
import os
import click
from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.text import Text
from rich.columns import Columns
from rich.rule import Rule
from rich.padding import Padding
from rich import box

console = Console()

# ── Brand colours ─────────────────────────────────────────────────────────────
PRIMARY   = "bold cyan"
SUCCESS   = "bold green"
ERROR     = "bold red"
WARNING   = "bold yellow"
DIM       = "dim white"
ACCENT    = "bold bright_white"

LOGO = """[bold cyan]
 ███████╗ ██████╗██████╗ ██╗██████╗ ████████╗
 ██╔════╝██╔════╝██╔══██╗██║██╔══██╗╚══██╔══╝
 ███████╗██║     ██████╔╝██║██████╔╝   ██║   
 ╚════██║██║     ██╔══██╗██║██╔═══╝    ██║   
 ███████║╚██████╗██║  ██║██║██║        ██║   
 ╚══════╝ ╚═════╝╚═╝  ╚═╝╚═╝╚═╝        ╚═╝   
[/bold cyan][dim]  Structural Chemical Representation In Plain Text[/dim]"""


# ── Helpers ───────────────────────────────────────────────────────────────────

def _get_mol_from_input(string: str, from_fmt: str):
    """
    Parse input string to (script_str, rdkit_mol).
    Returns (None, None) on failure — caller handles error display.
    """
    try:
        if from_fmt == "smiles":
            try:
                from rdkit import Chem, RDLogger
                RDLogger.DisableLog("rdApp.*")
                from script.rdkit_bridge import SCRIPTFromMol
                mol = Chem.MolFromSmiles(string)
                if mol is None:
                    return None, None, "RDKit could not parse SMILES"
                script_str = SCRIPTFromMol(mol)
                return script_str, mol, None
            except ImportError:
                return None, None, "RDKit is required for SMILES input (pip install rdkit)"

        elif from_fmt == "script":
            from script.parser import SCRIPTParser
            parser = SCRIPTParser()
            result = parser.parse(string)
            if not result.get("success"):
                return None, None, result.get("error", "Parse failed")
            script_str = string

            # Try to get RDKit mol for visualisation
            mol = None
            try:
                from script.rdkit_bridge import MolFromSCRIPT
                mol = MolFromSCRIPT(string)
            except Exception:
                pass
            return script_str, mol, None

    except Exception as e:
        return None, None, str(e)


def _draw_structure(mol, script_str: str, title: str = ""):
    """Render the molecule and print inside a panel."""
    try:
        from script.visualizer import render_mol
        art = render_mol(mol, width=66, height=20)
        panel_title = f"[{PRIMARY}]{title or script_str[:50]}[/{PRIMARY}]"
        console.print(Panel(
            art,
            title=panel_title,
            border_style="cyan",
            padding=(0, 1),
        ))
    except Exception as e:
        console.print(f"[{WARNING}]Could not render structure: {e}[/{WARNING}]")


def _show_properties(mol, script_str: str):
    """Print a property table for the molecule."""
    from script.visualizer import mol_properties
    props = mol_properties(mol) if mol else {}

    table = Table(box=box.SIMPLE, show_header=False, padding=(0, 2))
    table.add_column("Property", style=DIM, min_width=18)
    table.add_column("Value",    style=ACCENT)

    # Always show the SCRIPT string first
    script_display = script_str if len(script_str) <= 55 else script_str[:52] + "…"
    table.add_row("[bold cyan]SCRIPT[/bold cyan]", f"[cyan]{script_display}[/cyan]")

    for k, v in props.items():
        table.add_row(k, v)

    console.print(table)


def _pipe_input() -> str:
    """Read from stdin if it's being piped."""
    if not sys.stdin.isatty():
        return sys.stdin.read().strip()
    return ""


# ── Main group ────────────────────────────────────────────────────────────────

# ── Custom group: bare molecule strings route to draw ─────────────────────────

KNOWN_COMMANDS = {
    "draw", "canon", "validate", "convert",
    "info", "match", "bench", "about",
}

class SCRGroup(click.Group):
    """Routes bare molecule strings (not subcommand names) to `draw`."""
    def parse_args(self, ctx, args):
        if (args
                and not args[0].startswith("-")
                and args[0] not in KNOWN_COMMANDS):
            args = ["draw"] + args
        return super().parse_args(ctx, args)


@click.group(
    cls=SCRGroup,
    invoke_without_command=True,
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@click.option("-f", "--from", "from_fmt",
              default="script",
              type=click.Choice(["script", "smiles"], case_sensitive=False),
              show_default=True,
              help="Input format.")
@click.option("--no-draw", is_flag=True, help="Skip structure drawing, show info only.")
@click.version_option("3.0.0", "-v", "--version", message="SCRIPT v%(version)s")
@click.pass_context
def main(ctx, from_fmt, no_draw):
    """
    \b
    SCRIPT — canonical molecular notation.

    \b
    Examples:
      scr "C:C:C:C:C:C&6:"               draw benzene
      scr "c1ccccc1" -f smiles            SMILES input
      scr canon "CC(=O)O" -f smiles       canonicalize aspirin
      scr validate "C(C)(C)(C)(C)(C)C"   Sandhi check
      scr convert "CCO" -f smiles         convert to SCRIPT
      scr info "C:C:C:C:C:C&6:"          show properties
      echo "CCO" | scr convert -f smiles  pipe support
    """
    if ctx.invoked_subcommand is None:
        console.print(LOGO)
        console.print(ctx.get_help())


# ── draw ──────────────────────────────────────────────────────────────────────

@main.command()
@click.argument("string", required=False)
@click.option("-f", "--from", "from_fmt", default="script",
              type=click.Choice(["script", "smiles"], case_sensitive=False))
@click.option("--save", "save_path", default=None,
              help="Save PNG to file (requires RDKit).")
@click.option("-W", "--width",  default=68,  show_default=True, help="Terminal width.")
@click.option("-H", "--height", default=22,  show_default=True, help="Terminal height.")
def draw(string, from_fmt, save_path, width, height):
    """Draw a molecule as ASCII art in the terminal."""
    if not string:
        string = _pipe_input()
    if not string:
        console.print(f"[{ERROR}]Provide a molecule string or pipe one.[/{ERROR}]")
        sys.exit(1)

    script_str, mol, err = _get_mol_from_input(string, from_fmt)
    if err:
        console.print(f"[{ERROR}]✗ {err}[/{ERROR}]"); sys.exit(1)

    if mol is not None:
        from script.visualizer import render_mol
        art = render_mol(mol, width=width, height=height)
        console.print(Panel(
            art,
            title=f"[{PRIMARY}]{script_str[:55]}[/{PRIMARY}]",
            border_style="cyan",
            padding=(0, 1),
        ))
    else:
        console.print(f"[{WARNING}]Structure rendering requires RDKit.[/{WARNING}]")

    # Save PNG
    if save_path:
        try:
            from rdkit.Chem import Draw
            img = Draw.MolToImage(mol, size=(600, 400))
            img.save(save_path)
            console.print(f"[{SUCCESS}]✓ Saved to {save_path}[/{SUCCESS}]")
        except Exception as e:
            console.print(f"[{ERROR}]Could not save PNG: {e}[/{ERROR}]")


# ── canon ─────────────────────────────────────────────────────────────────────

@main.command()
@click.argument("string", required=False)
@click.option("-f", "--from", "from_fmt", default="script",
              type=click.Choice(["script", "smiles"], case_sensitive=False))
@click.option("-q", "--quiet", is_flag=True, help="Output canonical string only.")
def canon(string, from_fmt, quiet):
    """Return the canonical SCRIPT string for a molecule."""
    if not string:
        string = _pipe_input()
    if not string:
        console.print(f"[{ERROR}]Provide a molecule string.[/{ERROR}]"); sys.exit(1)

    try:
        if from_fmt == "smiles":
            from rdkit import Chem, RDLogger
            RDLogger.DisableLog("rdApp.*")
            from script.rdkit_bridge import SCRIPTFromMol
            mol = Chem.MolFromSmiles(string)
            if mol is None:
                console.print(f"[{ERROR}]✗ Invalid SMILES[/{ERROR}]"); sys.exit(1)
            canonical = SCRIPTFromMol(mol)
        else:
            from script.canonical import canonicalize_SCRIPT
            canonical = canonicalize_SCRIPT(string)

        if quiet:
            print(canonical)
        else:
            console.print(f"\n[{DIM}]Input   :[/{DIM}] {string}")
            console.print(f"[{PRIMARY}]Canonical:[/{PRIMARY}] [bold]{canonical}[/bold]\n")

    except Exception as e:
        console.print(f"[{ERROR}]✗ {e}[/{ERROR}]"); sys.exit(1)


# ── validate ──────────────────────────────────────────────────────────────────

@main.command()
@click.argument("string", required=False)
@click.option("--batch", "batch_file", default=None,
              help="Validate every line in a file. Exits 1 if any fail.")
@click.option("-q", "--quiet", is_flag=True, help="Suppress output, use exit code only.")
def validate(string, batch_file, quiet):
    """
    Validate a SCRIPT string via Sandhi valence checking.

    Exits 0 if valid, 1 if invalid.
    """
    from script.parser import SCRIPTParser
    parser = SCRIPTParser()

    def _check(s):
        result = parser.parse(s.strip())
        ok = result.get("success", False)
        err = result.get("error", "")
        return ok, err

    if batch_file:
        try:
            lines = open(batch_file).read().splitlines()
        except FileNotFoundError:
            console.print(f"[{ERROR}]File not found: {batch_file}[/{ERROR}]"); sys.exit(1)
        all_ok = True
        for i, line in enumerate(lines, 1):
            if not line.strip() or line.startswith("#"):
                continue
            ok, err = _check(line)
            if not quiet:
                status = f"[{SUCCESS}]✓[/{SUCCESS}]" if ok else f"[{ERROR}]✗[/{ERROR}]"
                note   = f"  [{DIM}]{err[:60]}[/{DIM}]" if not ok else ""
                console.print(f"  {status}  {line[:55]}{note}")
            if not ok:
                all_ok = False
        sys.exit(0 if all_ok else 1)

    if not string:
        string = _pipe_input()
    if not string:
        console.print(f"[{ERROR}]Provide a SCRIPT string.[/{ERROR}]"); sys.exit(1)

    ok, err = _check(string)
    if quiet:
        sys.exit(0 if ok else 1)

    if ok:
        console.print(f"\n  [{SUCCESS}]✓  Valid SCRIPT[/{SUCCESS}]\n")
    else:
        console.print(f"\n  [{ERROR}]✗  Invalid: {err}[/{ERROR}]\n")
        sys.exit(1)


# ── convert ───────────────────────────────────────────────────────────────────

@main.command()
@click.argument("string", required=False)
@click.option("-f", "--from", "from_fmt",
              default="smiles",
              type=click.Choice(["smiles", "script"], case_sensitive=False),
              show_default=True)
@click.option("-t", "--to", "to_fmt",
              default="script",
              type=click.Choice(["script", "smiles", "inchi", "inchikey"], case_sensitive=False),
              show_default=True)
@click.option("--batch", "batch_file", default=None,
              help="Convert every line in a file.")
@click.option("-q", "--quiet", is_flag=True, help="Output converted string only.")
def convert(string, from_fmt, to_fmt, batch_file, quiet):
    """
    Convert between SMILES, SCRIPT, InChI, and InChIKey.

    \b
    Examples:
      scr convert "CCO" -f smiles -t script
      scr convert "CCO" -f smiles -t inchi
      cat mols.smi | scr convert -f smiles
    """
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem import inchi as rdInchi
        RDLogger.DisableLog("rdApp.*")
        from script.rdkit_bridge import SCRIPTFromMol, MolFromSCRIPT
    except ImportError:
        console.print(f"[{ERROR}]RDKit is required for conversion.[/{ERROR}]")
        sys.exit(1)

    def _convert_one(s: str) -> str:
        s = s.strip()
        if not s or s.startswith("#"):
            return ""
        if from_fmt == "smiles":
            mol = Chem.MolFromSmiles(s)
            if mol is None:
                return f"ERROR: invalid SMILES: {s}"
        else:
            mol = MolFromSCRIPT(s)
            if mol is None:
                return f"ERROR: invalid SCRIPT: {s}"

        if to_fmt == "script":
            return SCRIPTFromMol(mol) or "ERROR: conversion failed"
        elif to_fmt == "smiles":
            return Chem.MolToSmiles(mol)
        elif to_fmt == "inchi":
            return rdInchi.MolToInchi(mol) or "ERROR"
        elif to_fmt == "inchikey":
            return rdInchi.InchiToInchiKey(rdInchi.MolToInchi(mol) or "") or "ERROR"
        return "ERROR: unknown target format"

    if batch_file:
        try:
            lines = open(batch_file).read().splitlines()
        except FileNotFoundError:
            console.print(f"[{ERROR}]File not found: {batch_file}[/{ERROR}]"); sys.exit(1)
        for line in lines:
            out = _convert_one(line)
            if out:
                print(out)
        return

    if not string:
        string = _pipe_input()

    if not string:
        console.print(f"[{ERROR}]Provide a molecule string or pipe one.[/{ERROR}]"); sys.exit(1)

    # Handle multi-line piped input
    results = [_convert_one(s) for s in string.splitlines() if s.strip()]

    if quiet or not sys.stdout.isatty():
        for r in results:
            print(r)
        return

    for inp, out in zip(string.splitlines(), results):
        if out.startswith("ERROR"):
            console.print(f"[{DIM}]{inp[:45]}[/{DIM}]  [{ERROR}]{out}[/{ERROR}]")
        else:
            console.print(f"[{DIM}]{inp[:45]}[/{DIM}]  [{SUCCESS}]→[/{SUCCESS}]  [bold]{out}[/bold]")


# ── info ──────────────────────────────────────────────────────────────────────

@main.command()
@click.argument("string", required=False)
@click.option("-f", "--from", "from_fmt", default="script",
              type=click.Choice(["script", "smiles"], case_sensitive=False))
def info(string, from_fmt):
    """Show molecular properties: formula, MW, rings, stereocenters, Lipinski."""
    if not string:
        string = _pipe_input()
    if not string:
        console.print(f"[{ERROR}]Provide a molecule string.[/{ERROR}]"); sys.exit(1)

    script_str, mol, err = _get_mol_from_input(string, from_fmt)
    if err:
        console.print(f"[{ERROR}]✗ {err}[/{ERROR}]"); sys.exit(1)

    _show_properties(mol, script_str)


# ── match ─────────────────────────────────────────────────────────────────────

@main.command()
@click.argument("query")
@click.argument("target")
@click.option("-f", "--from", "from_fmt", default="script",
              type=click.Choice(["script", "smiles"], case_sensitive=False))
def match(query, target, from_fmt):
    """Check if TARGET contains QUERY as a substructure."""
    try:
        from rdkit import Chem, RDLogger
        RDLogger.DisableLog("rdApp.*")
        from script.rdkit_bridge import MolFromSCRIPT

        def to_mol(s):
            if from_fmt == "smiles":
                return Chem.MolFromSmiles(s)
            return MolFromSCRIPT(s)

        q_mol = to_mol(query)
        t_mol = to_mol(target)
        if q_mol is None or t_mol is None:
            console.print(f"[{ERROR}]Could not parse one or both inputs.[/{ERROR}]")
            sys.exit(1)

        hit = t_mol.HasSubstructMatch(q_mol)
        if hit:
            console.print(f"\n  [{SUCCESS}]✓  Substructure match found[/{SUCCESS}]\n")
        else:
            console.print(f"\n  [{ERROR}]✗  No substructure match[/{ERROR}]\n")
        sys.exit(0 if hit else 1)
    except ImportError:
        console.print(f"[{ERROR}]RDKit is required for substructure matching.[/{ERROR}]")
        sys.exit(1)


# ── bench ─────────────────────────────────────────────────────────────────────

@main.command()
@click.option("--quick", is_flag=True, help="Run internal 71-compound suite only.")
def bench(quick):
    """Run the SCRIPT benchmark suite."""
    console.print(Rule("[bold cyan]SCRIPT Benchmark[/bold cyan]"))

    if quick:
        import subprocess, sys
        result = subprocess.run(
            [sys.executable, "benchmark.py"],
            capture_output=True, text=True
        )
        console.print(result.stdout)
        if result.returncode != 0:
            console.print(f"[{ERROR}]{result.stderr}[/{ERROR}]")
    else:
        console.print(f"[{DIM}]Running large benchmark — this takes a few minutes...[/{DIM}]")
        import subprocess, sys
        result = subprocess.run(
            [sys.executable, "tests/large_benchmark.py"],
            capture_output=False, text=True
        )


# ── about ─────────────────────────────────────────────────────────────────────

@main.command()
def about():
    """Show SCRIPT version, design basis, and links."""
    console.print(LOGO)
    console.print()

    table = Table(box=box.SIMPLE, show_header=False, padding=(0, 2))
    table.add_column("", style=DIM,    min_width=16)
    table.add_column("", style=ACCENT)
    table.add_row("Version",     "3.0.0")
    table.add_row("Grammar",     "LALR(1) via Lark")
    table.add_row("Design basis","Paninian linguistics (Ashtadhyayi)")
    table.add_row("Canonicalization","Morgan–Weisfeiler–Lehman")
    table.add_row("Stereo",      "CIP-anchored, parse-order-independent")
    table.add_row("Materials",   "Alloys, phases, surfaces, polymers, spin states")
    table.add_row("Benchmark",   "99.6% on 2,500 ChEMBL drug-like compounds")
    table.add_row("License",     "MIT + Commons Clause")
    table.add_row("Repo",        "https://github.com/sangeet01/script")
    console.print(table)
    console.print()
    console.print(f'[{DIM}]"A linear script to unfold molecular complexity — from the singlet to the surface."[/{DIM}]')
    console.print()


# ── entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    main()
