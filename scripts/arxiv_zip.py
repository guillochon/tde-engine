#!/usr/bin/env python3
"""
Zip the files arXiv needs to compile the manuscript.

The .tex, .bib, figures, and class/bst files sit at the archive root
(not under paper/). The compiled PDF, .bbl, data products, and aux files are
omitted.

    python scripts/arxiv_zip.py              # write tde-engine-arxiv.zip
    python scripts/arxiv_zip.py --out PATH   # choose the zip path
    python scripts/arxiv_zip.py --dry-run    # list files, do not write
"""

import argparse
import sys
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "paper"
DEFAULT_ZIP = ROOT / "tde-engine-arxiv.zip"

# Names as they appear in paper/ and at the zip root.
FILES = [
    "engine.tex",
    "engineNotes.bib",
    "apj.bst",
    "aastex.cls",
    "engine-diagram.pdf",
    "mom-rates.pdf",
    "dominance.pdf",
    "domratio.pdf",
    "window.pdf",
]


def collect():
    missing = [name for name in FILES if not (PAPER / name).is_file()]
    if missing:
        print("Missing from paper/:", file=sys.stderr)
        for name in missing:
            print(f"  {name}", file=sys.stderr)
        sys.exit(1)
    return [(PAPER / name, name) for name in FILES]


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument(
        "--out",
        type=Path,
        default=DEFAULT_ZIP,
        help=f"zip path (default: {DEFAULT_ZIP.name})",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="list files and sizes, do not write",
    )
    args = parser.parse_args()

    items = collect()
    total = 0
    for src, name in items:
        size = src.stat().st_size
        total += size
        print(f"  {name:24s} {size / 1024:8.1f} KB")
    print(f"  {'total':24s} {total / 1024:8.1f} KB")

    if args.dry_run:
        return

    out = args.out.resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(out, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for src, name in items:
            zf.write(src, arcname=name)
    print(f"Wrote {out} ({out.stat().st_size / 1024:.1f} KB)")


if __name__ == "__main__":
    main()
