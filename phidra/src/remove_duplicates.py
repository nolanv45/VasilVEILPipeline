#!/usr/bin/env python3

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path
from typing import Iterable, Set, Tuple

from Bio import SeqIO


def _canon(seq: str) -> str:
    """Canonicalize sequence for equality."""
    return "".join(seq.split()).upper()


def _md5(seq: str) -> str:
    return hashlib.md5(seq.encode("utf-8")).hexdigest()


def _digest_set(fasta_path: Path) -> Set[str]:
    digs: Set[str] = set()
    for rec in SeqIO.parse(str(fasta_path), "fasta"):
        digs.add(_md5(_canon(str(rec.seq))))
    return digs


def remove_sequences(file1: Path, file2: Path, out_path: Path) -> Tuple[int, int, int]:
    to_remove = _digest_set(file2)

    seen = 0
    kept = 0
    removed = 0
    removed_ids: list[str] = []

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as out_h:
        for rec in SeqIO.parse(str(file1), "fasta"):
            seen += 1
            dig = _md5(_canon(str(rec.seq)))
            if dig in to_remove:
                removed += 1
                removed_ids.append(rec.id)
                continue
            SeqIO.write(rec, out_h, "fasta")
            kept += 1

    if removed_ids:
        with (out_path.with_suffix(out_path.suffix + ".removed_ids.txt")).open(
            "w", encoding="utf-8"
        ) as rm:
            rm.write("\n".join(removed_ids) + "\n")

    return seen, kept, removed


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Remove sequences from file1 whose sequence content occurs in file2."
    )
    ap.add_argument("file1", help="Source FASTA to be filtered")
    ap.add_argument("file2", help="FASTA whose sequences should be removed from file1")
    ap.add_argument("output", help="Filtered output FASTA path")
    args = ap.parse_args()

    f1 = Path(args.file1)
    f2 = Path(args.file2)
    out = Path(args.output)

    if not f1.exists():
        raise SystemExit(f"file1 not found: {f1}")
    if not f2.exists():
        raise SystemExit(f"file2 not found: {f2}")

    seen, kept, removed = remove_sequences(f1, f2, out)
    print(
        f"Input: {seen}; kept: {kept}; removed: {removed}; output: {out}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
