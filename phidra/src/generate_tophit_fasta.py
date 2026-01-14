#!/usr/bin/env python3

from __future__ import annotations

import argparse
import sys
from collections import OrderedDict
from pathlib import Path
from typing import Iterable, List, Tuple

from Bio import SeqIO


def _read_query_ids(tsv_path: Path, header_name: str = "Query_ID") -> List[str]:
    ids: "OrderedDict[str, None]" = OrderedDict()
    with tsv_path.open("r", encoding="utf-8") as f:
        first = f.readline()
        if not first:
            return []
        parts = first.rstrip("\n").split("\t")
        try:
            idx = parts.index(header_name)
        except ValueError:
            # No header: treat the first line as data and assume first column
            idx = 0
            row = parts
            if row and row[0] and not row[0].startswith("#"):
                ids.setdefault(row[idx], None)
        for line in f:
            if not line.strip():
                continue
            row = line.rstrip("\n").split("\t")
            if len(row) <= idx:
                continue
            qid = row[idx]
            if not qid or qid.startswith("#"):
                continue
            ids.setdefault(qid, None)
    return list(ids.keys())


def _write_missing(path: Path, missing: Iterable[str]) -> None:
    miss = list(missing)
    if not miss:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        for mid in miss:
            f.write(mid + "\n")


def extract_fasta(top_tsv: Path, fasta_path: Path, out_fasta: Path) -> Tuple[int, int, int]:
    qids = _read_query_ids(top_tsv)
    if not qids:
        raise SystemExit(f"No query IDs found in: {top_tsv}")

    out_fasta.parent.mkdir(parents=True, exist_ok=True)

    # Avoids loading entire FASTA into memory, scales to very large files
    idx = SeqIO.index(str(fasta_path), "fasta")

    written = 0
    missing: List[str] = []

    with out_fasta.open("w", encoding="utf-8") as out_h:
        for qid in qids:
            rec = idx.get(qid)
            if rec is None:
                alt = qid.split()[0]
                rec = idx.get(alt)
            if rec is None:
                missing.append(qid)
                continue
            SeqIO.write(rec, out_h, "fasta")
            written += 1

    idx.close()

    _write_missing(out_fasta.with_suffix(out_fasta.suffix + ".missing.txt"), missing)
    return len(qids), written, len(missing)


def main() -> None:
    ap = argparse.ArgumentParser(description="Generate FASTA of TopHit queries.")
    ap.add_argument("tophit_file", help="TopHit TSV (expects a 'Query_ID' header, else uses first column)")
    ap.add_argument("fasta_file", help="Original query FASTA")
    ap.add_argument("output_fasta", help="Output FASTA path")
    args = ap.parse_args()

    top = Path(args.tophit_file)
    fa = Path(args.fasta_file)
    out = Path(args.output_fasta)

    if not top.exists():
        raise SystemExit(f"TopHit file not found: {top}")
    if not fa.exists():
        raise SystemExit(f"FASTA not found: {fa}")

    total, written, missing = extract_fasta(top, fa, out)
    print(
        f"TopHit IDs: {total}; written: {written}; missing: {missing}; output: {out}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
