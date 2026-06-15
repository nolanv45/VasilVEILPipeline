#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import sys
import os
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Tuple

from Bio import SeqIO

HEADERS = [
    "Query_ID",
    "Pfam_IDs",
    "Env_Starts",
    "Env_Ends",
    "Pfam_Functions",
    "IDA_ID",
    "IDA_Text",
    "Representative_Domains",
]


@dataclass
class PfamAgg:
    pfam_ids: List[str] = field(default_factory=list)
    pfam_funcs: List[str] = field(default_factory=list)
    env_starts: List[int] = field(default_factory=list)
    env_ends: List[int] = field(default_factory=list)
    ida_id: str = ""
    ida_text: str = ""
    rep_domains: str = ""

    def to_row(self, qid: str) -> Dict[str, str]:
        return {
            "Query_ID": qid,
            "Pfam_IDs": "|".join(self.pfam_ids),
            "Env_Starts": "|".join(map(str, self.env_starts)),
            "Env_Ends": "|".join(map(str, self.env_ends)),
            "Pfam_Functions": "|".join(self.pfam_funcs),
            "IDA_ID": self.ida_id,
            "IDA_Text": self.ida_text,
            "Representative_Domains": self.rep_domains,
        }


def _process_pfam_scan(pfam_scan_file: Path) -> Dict[str, PfamAgg]:
    agg: Dict[str, PfamAgg] = defaultdict(PfamAgg)
    with pfam_scan_file.open("r", encoding="utf-8") as f:
        rdr = csv.DictReader(f, delimiter="\t")
        for row in rdr:
            if not row:
                continue
            if row.get("significance", "0") != "1":
                continue
            sid = row["seq_id"].strip()
            if not sid:
                continue
            agg[sid].pfam_ids.append(row["hmm_acc"].split(".")[0])
            agg[sid].pfam_funcs.append(row["hmm_name"].strip())
            try:
                agg[sid].env_starts.append(int(row["env_start"]))
                agg[sid].env_ends.append(int(row["env_end"]))
            except (KeyError, ValueError):
                # Skip malformed coordinates for this row
                agg[sid].pfam_ids.pop()
                agg[sid].pfam_funcs.pop()
    return agg


def _process_ida_text(txt: str) -> str:
    # Example: "PF00069:Kinase-PF07714:Pkinase_Tyr" -> "PF00069|PF07714"
    parts = []
    for token in txt.split("-"):
        left = token.split(":")[0].strip()
        if left:
            parts.append(left)
    return "|".join(parts)


def _chains_equal(a: str, b: str) -> bool:
    return Counter(a.split("|")) == Counter(b.split("|"))


def _load_ida_valid(ida_file: Path) -> Dict[str, Dict[str, str]]:
    valid: Dict[str, Dict[str, str]] = {}
    with ida_file.open("r", encoding="utf-8") as f:
        rdr = csv.DictReader(f, delimiter="\t")
        for row in rdr:
            chain = _process_ida_text(row.get("IDA Text", ""))
            if not chain:
                continue
            valid[chain] = {
                "IDA_ID": row.get("IDA ID", ""),
                "IDA_Text": row.get("IDA Text", ""),
                "Representative_Domains": row.get("Representative Domains", ""),
            }
    return valid


def _validate(agg: Dict[str, PfamAgg], ida_file: Path) -> Tuple[Dict[str, PfamAgg], Dict[str, PfamAgg]]:
    valid_meta = _load_ida_valid(ida_file)

    validated: Dict[str, PfamAgg] = {}
    unvalidated: Dict[str, PfamAgg] = {}

    for sid, data in agg.items():
        chain = "|".join(data.pfam_ids)
        found = False
        for vchain, meta in valid_meta.items():
            if _chains_equal(chain, vchain):
                data.ida_id = meta["IDA_ID"]
                data.ida_text = meta["IDA_Text"]
                data.rep_domains = meta["Representative_Domains"]
                validated[sid] = data
                found = True
                break
        if not found:
            # Ensure IDA fields exist but are empty placeholders
            data.ida_id = data.ida_id or ""
            data.ida_text = data.ida_text or ""
            data.rep_domains = data.rep_domains or ""
            unvalidated[sid] = data

    return validated, unvalidated


def _write_report(data: Dict[str, PfamAgg], out_tsv: Path) -> None:
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w", encoding="utf-8", newline="") as o:
        w = csv.DictWriter(o, HEADERS, delimiter="\t")
        w.writeheader()
        for sid, d in data.items():
            w.writerow(d.to_row(sid))


def _sanitize(name: str) -> str:
    return "".join(ch if ch.isalnum() or ch in ("_", "-") else "_" for ch in name)


def _extract_seqs(fasta: Path, data: Dict[str, PfamAgg], out_fa: Path, domain_only: bool = True) -> Tuple[int, int]:
    out_fa.parent.mkdir(parents=True, exist_ok=True)
    total = 0
    written = 0
    with out_fa.open("w", encoding="utf-8") as o:
        for rec in SeqIO.parse(str(fasta), "fasta"):
            if rec.id not in data:
                continue
            total += 1
            if domain_only:
                # Sort domains by start for stable output
                for fn, st, en in sorted(
                    zip(data[rec.id].pfam_funcs, data[rec.id].env_starts, data[rec.id].env_ends),
                    key=lambda x: x[1],
                ):
                    # Clamp coordinates to sequence bounds; skip invalid ranges
                    st = max(1, int(st))
                    en = min(len(rec.seq), int(en))
                    if en < st:
                        continue
                    header = f">{rec.id}_{_sanitize(fn)}_{st}_{en}"
                    o.write(header + "\n")
                    o.write(str(rec.seq[st - 1 : en]) + "\n")
                    written += 1
            else:
                o.write(f">{rec.id}\n{rec.seq}\n")
                written += 1
    return total, written


def main() -> None:
    ap = argparse.ArgumentParser(description="Validate Pfam domains with IDA table and emit reports + FASTAs.")
    ap.add_argument("pfam_scan", help="Pfam scan TSV")
    ap.add_argument("fasta", help="Original proteins FASTA")
    ap.add_argument("ida_file", help="IDA TSV")
    ap.add_argument("out_dir", help="Output directory")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    val_dir = out_dir / "validated_ida_report"
    unv_dir = out_dir / "unvalidated_ida_report"
    for d in (val_dir, unv_dir):
        d.mkdir(parents=True, exist_ok=True)

    agg = _process_pfam_scan(Path(args.pfam_scan))
    val, unv = _validate(agg, Path(args.ida_file))

    # reports
    _write_report(val, val_dir / "pfam_validated_report.tsv")
    _write_report(unv, unv_dir / "pfam_unvalidated_report.tsv")

    # sequences
    _extract_seqs(Path(args.fasta), val, val_dir / "domains.fa", True)
    _extract_seqs(Path(args.fasta), val, val_dir / "full_proteins.fa", False)
    _extract_seqs(Path(args.fasta), unv, unv_dir / "domains.fa", True)
    _extract_seqs(Path(args.fasta), unv, unv_dir / "full_proteins.fa", False)

    print(
        f"Validated: {len(val)} | Unvalidated: {len(unv)} | Out: {out_dir}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit("Usage: validate_pfam.py <pfam_scan.tsv> <proteins.fa> <ida.tsv> <out_dir>")
    main()
