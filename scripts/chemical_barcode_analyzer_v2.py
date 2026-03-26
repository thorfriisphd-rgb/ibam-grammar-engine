#!/usr/bin/env python3
"""
chemical_barcode_analyzer_v2.py

Compute machine-readable summaries from projected alignment FASTAs:
- motif width
- per-column occupancy
- per-column dominant residue + freq
- per-column dominant chemical class + freq
- per-column Shannon entropy (AA-level + chemical-level)
- "chemical barcode" (dominant chemical class per column)
- anchor columns (top N by information / dominance)

Designed for projected interface alignments (3D-projected MSAs).

Usage:
  python chemical_barcode_analyzer_v2.py --root /path/to/MyhT_samle_perm --gate core35_chem70 --out summary.tsv
  python chemical_barcode_analyzer_v2.py --fasta /path/to/MG_projected_trimmed_core35_chem70.fa --out single.tsv

Notes:
- Do NOT parse PDFs/SVG. This reads the FASTA used to make them.
- Deterministic file selection (sorted).
"""

from __future__ import annotations
import argparse
import math
import csv
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from collections import Counter

AA_ALPHABET = set(list("ACDEFGHIKLMNPQRSTVWY"))

CHEM_MAP = {
    # Basic
    "K": "B", "R": "B", "H": "B",
    # Acidic
    "D": "A", "E": "A",
    # Polar (incl. Gly; treat as polar/other-friendly)
    "S": "P", "T": "P", "N": "P", "Q": "P", "C": "P", "G": "P",
    # Hydrophobic
    "A": "H", "V": "H", "I": "H", "L": "H", "M": "H", "F": "H", "W": "H", "Y": "H", "P": "H",
}
CHEM_ALPHABET = set(list("BAPH"))  # plus X for unknown

GAP_CHARS = set(["-", ".", "~"])

RUN_DIR_RE = re.compile(r"^(MG_core_.*)$")
GATE_RE = re.compile(r"(core\d+_chem\d+)")
@dataclass
class ColumnStats:
    pos: int
    occupancy: float
    aa_dom: str
    aa_dom_frac: float
    chem_dom: str
    chem_dom_frac: float
    aa_entropy: float
    chem_entropy: float
    chem_fracs: dict

def parse_fasta(path: Path) -> Dict[str, str]:
    seqs: Dict[str, List[str]] = {}
    name: Optional[str] = None
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            else:
                if name is None:
                    raise ValueError(f"FASTA parse error: sequence line before header in {path}")
                seqs[name].append(line)
    out = {k: "".join(v) for k, v in seqs.items()}
    if not out:
        raise ValueError(f"No sequences found in {path}")
    lengths = {len(s) for s in out.values()}
    if len(lengths) != 1:
        raise ValueError(f"Not an alignment: multiple lengths in {path}: {sorted(lengths)}")
    return out

def aa_to_chem(aa: str) -> str:
    aa = aa.upper()
    if aa in GAP_CHARS:
        return "-"
    return CHEM_MAP.get(aa, "X")

def shannon_entropy(counter: Counter, total: int) -> float:
    if total <= 0:
        return 0.0
    H = 0.0
    for c in counter.values():
        p = c / total
        H -= p * math.log2(p)
    return H

def info_content_from_entropy(entropy: float, alphabet_size: int) -> float:
    # Classic sequence-logo-ish: log2(|alphabet|) - H
    return max(0.0, math.log2(alphabet_size) - entropy)

def compute_columns(seqs: Dict[str, str]) -> List[ColumnStats]:
    taxa = list(seqs.keys())
    n = len(taxa)
    L = len(next(iter(seqs.values())))
    cols: List[ColumnStats] = []

    for i in range(L):
        col = [seqs[t][i] for t in taxa]
        nongap = [c for c in col if c not in GAP_CHARS]
        occ = len(nongap) / n

        if not nongap:
            cols.append(ColumnStats(
                pos=i+1, occupancy=0.0,
                aa_dom=".", aa_dom_frac=0.0,
                chem_dom=".", chem_dom_frac=0.0,
                aa_entropy=0.0, chem_entropy=0.0,
                chem_fracs={}
            ))
            continue

        aa_counts = Counter([c.upper() for c in nongap])
        aa_dom, aa_dom_n = aa_counts.most_common(1)[0]
        aa_dom_frac = aa_dom_n / len(nongap)
        aa_H = shannon_entropy(aa_counts, len(nongap))

        chem_col = [aa_to_chem(c) for c in nongap]
        chem_counts = Counter(chem_col)
        chem_dom, chem_dom_n = chem_counts.most_common(1)[0]
        chem_dom_frac = chem_dom_n / len(nongap)
        chem_fracs = {k: v / len(nongap) for k, v in chem_counts.items()}
        chem_H = shannon_entropy(chem_counts, len(nongap))

        cols.append(ColumnStats(
            pos=i+1, occupancy=occ,
            aa_dom=aa_dom, aa_dom_frac=aa_dom_frac,
            chem_dom=chem_dom, chem_dom_frac=chem_dom_frac,
            aa_entropy=aa_H, chem_entropy=chem_H,
            chem_fracs=chem_fracs
        ))
    return cols

def chemical_barcode(cols: List[ColumnStats]) -> str:
    return "".join([c.chem_dom if c.chem_dom not in ["."] else "X" for c in cols])

def anchors(cols: List[ColumnStats], top: int) -> List[Tuple[int, float, float, str]]:
    # Anchor score: chemical info * occupancy * chem_dom_frac
    # (stable for 3D interfaces; less sensitive to AA churn)
    out = []
    for c in cols:
        chem_info = info_content_from_entropy(c.chem_entropy, alphabet_size=4)
        score = chem_info * c.occupancy * c.chem_dom_frac
        out.append((c.pos, score, c.chem_dom_frac, c.chem_dom))
    out.sort(key=lambda x: x[1], reverse=True)
    return out[:top]

def barcode_similarity(b1: str, b2: str) -> float:
    # Hamming similarity over aligned barcodes; if lengths differ, compare to min length.
    m = min(len(b1), len(b2))
    if m == 0:
        return 0.0
    eq = sum(1 for i in range(m) if b1[i] == b2[i])
    return eq / m

def find_fastas(run_dir: Path) -> List[Path]:
    # Broad scan: projection/ and alignment/ for fasta-like.
    pats = [
        run_dir / "projection" / "*.fa*",
        run_dir / "projection" / "*.fasta*",
        run_dir / "alignment" / "*.fa*",
        run_dir / "alignment" / "*.fasta*",
    ]
    files: List[Path] = []
    for pat in pats:
        files.extend(sorted(Path().glob(str(pat))))
    # De-dup + sort deterministically by path
    uniq = sorted({f.resolve() for f in files if f.is_file()})
    return uniq

def pick_fasta(files: List[Path], gate: Optional[str]) -> Optional[Path]:
    if not files:
        return None
    if gate:
        gate_files = [f for f in files if gate in f.name]
        if gate_files:
            return sorted(gate_files)[0]
        return None
    # Otherwise pick the “most projected-looking” file
    def score(p: Path) -> int:
        n = p.name.lower()
        s = 0
        if "project" in n: s += 10
        if "trim" in n: s += 5
        if "mg" in n: s += 3
        if "myh" in n: s += 3
        if n.endswith(".fa") or n.endswith(".fasta"): s += 1
        return s
    return sorted(files, key=score, reverse=True)[0]

def parse_gate_from_name(name: str) -> str:
    m = GATE_RE.search(name)
    return m.group(1) if m else "coreNA_chemNA"

def infer_n_from_name(name: str, seq_count: int) -> int:
    m = re.search(r"_n(\d+)_", name)
    if m:
        return int(m.group(1))
    return seq_count

def analyze_one(fasta: Path, run_label: str, gate_override: Optional[str], top: int):
    seqs = parse_fasta(fasta)
    cols = compute_columns(seqs)
    bc = chemical_barcode(cols)
    anc = anchors(cols, top=top)
    gate = gate_override or parse_gate_from_name(fasta.name)
    return seqs, cols, bc, anc, gate

def write_per_column(path: Path, run_id: str, gate: str, cols: List[ColumnStats]):
    with path.open("w", encoding="utf-8") as out:
        out.write("\t".join([
            "run_id","gate","pos","occupancy",
            "aa_dom","aa_dom_frac","chem_dom","chem_dom_frac",
            "aa_entropy","chem_entropy",
            "aa_info","chem_info"
        ]) + "\n")
        for c in cols:
            aa_info = info_content_from_entropy(c.aa_entropy, 20)
            chem_info = info_content_from_entropy(c.chem_entropy, 4)
            out.write("\t".join(map(str, [
                run_id, gate, c.pos, f"{c.occupancy:.4f}",
                c.aa_dom, f"{c.aa_dom_frac:.4f}",
                c.chem_dom, f"{c.chem_dom_frac:.4f}",
                f"{c.aa_entropy:.4f}", f"{c.chem_entropy:.4f}",
                f"{aa_info:.4f}", f"{chem_info:.4f}"
            ])) + "\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=str, help="Root containing MG_core_* run folders")
    ap.add_argument("--fasta", type=str, help="Analyze a single FASTA")
    ap.add_argument("--gate", type=str, default=None, help="Gate tag to select, e.g. core35_chem70")
    ap.add_argument("--out", type=str, default="barcode_summary.tsv", help="Output TSV (one row per run)")
    ap.add_argument("--per_column_out", type=str, default=None, help="Optional per-column TSV")
    ap.add_argument("--top", type=int, default=5, help="Top N anchors")
    ap.add_argument("--reference_run", type=str, default=None,
                    help="Run_id to use as reference for barcode similarity (default: full n25 if present, else first run).")
    args = ap.parse_args()

    if not args.root and not args.fasta:
        raise SystemExit("Provide --root or --fasta")

    rows = []
    percol_rows: List[Tuple[str,str,List[ColumnStats]]] = []

    # Collect analyses
    if args.fasta:
        fasta = Path(args.fasta).expanduser().resolve()
        run_id = fasta.stem
        seqs, cols, bc, anc, gate = analyze_one(fasta, run_id, args.gate, args.top)
        rows.append((run_id, gate, len(seqs), len(cols), bc, anc))
        percol_rows.append((run_id, gate, cols))
    else:
        root = Path(args.root).expanduser().resolve()
        run_dirs = sorted([p for p in root.iterdir() if p.is_dir() and p.name.startswith("MG_core_")])
        for rd in run_dirs:
            files = find_fastas(rd)
            picked = pick_fasta(files, args.gate)
            if not picked:
                rows.append((rd.name, args.gate or "coreNA_chemNA", "NA", "NA", "NA", []))
                continue
            seqs, cols, bc, anc, gate = analyze_one(picked, rd.name, args.gate, args.top)
            n_taxa = infer_n_from_name(rd.name, len(seqs))
            rows.append((rd.name, gate, n_taxa, len(cols), bc, anc))
            percol_rows.append((rd.name, gate, cols))

    # Choose reference barcode for similarity
    ref_idx = 0
    if args.reference_run:
        for i, r in enumerate(rows):
            if r[0] == args.reference_run:
                ref_idx = i
                break
    else:
        # Prefer n25 if present
        for i, r in enumerate(rows):
            if isinstance(r[2], int) and r[2] == 25:
                ref_idx = i
                break

    ref_bc = rows[ref_idx][4] if rows and rows[ref_idx][4] != "NA" else None

    out_path = Path(args.out).resolve()
    with out_path.open("w", encoding="utf-8") as out:
        out.write("\t".join([
            "run_id","gate","n_taxa","width",
            "mean_occupancy","barcode",
            "anchors",
            "barcode_similarity_to_ref"
        ]) + "\n")

        for run_id, gate, n_taxa, width, bc, anc in rows:
            if bc == "NA" or width == "NA":
                out.write("\t".join(map(str, [run_id, gate, n_taxa, width, "NA", bc, "NA", "NA"])) + "\n")
                continue
            # mean occupancy from barcode is not possible; compute from per-column if available
            cols = next((c for (rid, g, c) in percol_rows if rid == run_id and g == gate), None)
            mean_occ = sum(cc.occupancy for cc in cols)/len(cols) if cols else 0.0
            anc_str = ",".join([f"{pos}{chem}:{domf:.2f}" for (pos, score, domf, chem) in anc])
            sim = barcode_similarity(bc, ref_bc) if ref_bc else 0.0
            out.write("\t".join(map(str, [
                run_id, gate, n_taxa, width,
                f"{mean_occ:.4f}", bc,
                anc_str,
                f"{sim:.4f}"
            ])) + "\n")

    if args.per_column_out:
        per_path = Path(args.per_column_out).resolve()

        # Write concatenated per-column table
        with per_path.open("w", encoding="utf-8") as out:
            out.write("\t".join([
                "run_id","gate","pos","occupancy",
                "aa_dom","aa_dom_frac","chem_dom","chem_dom_frac",
                "aa_entropy","chem_entropy","aa_info","chem_info"
            ]) + "\n")
            for run_id, gate, cols in percol_rows:
                for c in cols:
                    aa_info = info_content_from_entropy(c.aa_entropy, 20)
                    chem_info = info_content_from_entropy(c.chem_entropy, 4)
                    out.write("\t".join(map(str, [
                        run_id, gate, c.pos, f"{c.occupancy:.4f}",
                        c.aa_dom, f"{c.aa_dom_frac:.4f}",
                        c.chem_dom, f"{c.chem_dom_frac:.4f}",
                        f"{c.aa_entropy:.4f}", f"{c.chem_entropy:.4f}",
                        f"{aa_info:.4f}", f"{chem_info:.4f}",
                    ])) + "\n")

        # Write concatenated chemistry-composition table
        chemcomp_out = str(per_path).replace("_percol.tsv", "_chemcomp.tsv")
        with open(chemcomp_out, "w", encoding="utf-8") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(["run_id", "gate", "pos", "chem_class", "frac"])
            for run_id, gate, cols in percol_rows:
                for c in cols:
                    for chem_class, frac in sorted(c.chem_fracs.items()):
                        writer.writerow([run_id, gate, c.pos, chem_class, f"{frac:.4f}"])

    print(f"[OK] Wrote summary: {out_path}")
    if args.per_column_out:
        print(f"[OK] Wrote per-column: {Path(args.per_column_out).resolve()}")
        print(f"[OK] Wrote chemcomp: {chemcomp_out}")

if __name__ == "__main__":
    main()
