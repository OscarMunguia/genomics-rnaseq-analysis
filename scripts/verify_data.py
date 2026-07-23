"""Check that rnaseq-01.xlsx exists and has expected FPKM columns."""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA = ROOT / "rnaseq-01.xlsx"

EXPECTED_COLS = [
    "Control-1 (FPKM)",
    "Control-2 (FPKM)",
    "Control-3 (FPKM)",
    "APOEhi-1 (FPKM)",
    "APOEhi-2 (FPKM)",
    "APOEhi-3 (FPKM)",
    "APOElo-1 (FPKM)",
    "APOElo-2 (FPKM)",
    "APOElo-3 (FPKM)",
]


def verify(path: Path = DEFAULT_DATA) -> bool:
    ok = True
    if not path.exists():
        print(f"FAIL: missing {path}")
        print("  See docs/DOWNLOAD_DATA.md")
        return False

    size_mb = path.stat().st_size / (1024 * 1024)
    print(f"OK: found {path.name} ({size_mb:.1f} MB)")

    df = pd.read_excel(path, nrows=2)
    header = df.iloc[0].tolist()
    missing = [c for c in EXPECTED_COLS if c not in header]
    if missing:
        print("FAIL: missing FPKM columns:")
        for col in missing:
            print(f"  - {col}")
        ok = False
    else:
        print(f"OK: all {len(EXPECTED_COLS)} FPKM columns present")

    if "GeneSymbol" not in header:
        print("FAIL: missing GeneSymbol column")
        ok = False
    else:
        print("OK: GeneSymbol column present")

    return ok


def main() -> None:
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_DATA
    sys.exit(0 if verify(path) else 1)


if __name__ == "__main__":
    main()
