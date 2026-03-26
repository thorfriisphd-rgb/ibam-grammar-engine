#!/usr/bin/env python3
"""
complementary_pattern_analyzer_v2.py

Analyze chemical complementarity between two projected interface alignments
(e.g., C12 vs MyhT) across shared taxa.

Modes:
1) contact-map restricted (recommended):
   Provide a TSV with pairs of positions: c12_pos  myht_pos
2) all-pairs exploratory scan:
   Evaluate all position pairs (short motifs only).

Also supports null testing by shuffling columns (preserves per-column chemistry).

Inputs:
- FASTA A (C12 projected)
- FASTA B (MyhT projected)

Outputs:
- TSV of complementarity per pair (and optionally null z-scores)
"""

from __future__ import annotations
import argparse
import random
import math
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from collections import Counter

CHEM_MAP = {
    "K":"B","R":"B","H":"B",
    "D":"A","E":"A",
    "S":"P","T":"P","N":"P","Q":"P","C":"P","G":"P",
    "A":"H","V":"H","I":"H","L":"H","M":"H","F":"H","W":"H","Y":"H","P":"H",
}
GAP = set(["-",".","~"])

def parse_fasta(path: Path) -> Dict[str, str]:
    seqs: Dict[str, List[str]] = {}
    name = None
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line: continue
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            else:
                if name is None:
                    raise ValueError(f"FASTA parse error in {path}")
                seqs[name].append(line)
    out = {k:"".join(v) for k,v in seqs.items()}
    if not out:
        raise ValueError(f"No sequences in {path}")
    Ls = {len(s) for s in out.values()}
    if len(Ls) != 1:
        raise ValueError(f"Not aligned in {path}: lengths={sorted(Ls)}")
    return out

def aa_to_chem(a: str) -> str:
    a = a.upper()
    if a in GAP: return "-"
    return CHEM_MAP.get(a, "X")

def chem_seq(seq: str) -> str:
    return "".join(aa_to_chem(a) for a in seq)

def is_complement(ca: str, cb: str) -> int:
    # Complementarity definition (conservative):
    # Basic↔Acidic is complement, Hydrophobic↔Hydrophobic is complement.
    # You can extend later.
    if ca == "-" or cb == "-": return 0
    if (ca == "B" and cb == "A") or (ca == "A" and cb == "B"): return 1
    if ca == "H" and cb == "H": return 1
    return 0

def read_contact_map(path: Path) -> List[Tuple[int,int]]:
    pairs = []
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"): continue
            parts = line.split()
            if len(parts) < 2:
                raise ValueError(f"Bad contact map line: {line}")
            pairs.append((int(parts[0]), int(parts[1])))
    return pairs

def complement_score_for_pair(chemA: Dict[str,str], chemB: Dict[str,str], posA: int, posB: int, taxa: List[str]) -> float:
    # Positions are 1-based
    hits = 0
    denom = 0
    i = posA - 1
    j = posB - 1
    for t in taxa:
        a = chemA[t][i]
        b = chemB[t][j]
        if a == "-" or b == "-":
            continue
        denom += 1
        hits += is_complement(a, b)
    return hits / denom if denom else 0.0

def shuffle_columns(seqs: Dict[str,str], rng: random.Random) -> Dict[str,str]:
    # Shuffle columns identically across taxa (preserves per-column distributions)
    taxa = list(seqs.keys())
    L = len(next(iter(seqs.values())))
    cols = [[seqs[t][i] for t in taxa] for i in range(L)]
    rng.shuffle(cols)
    out = {t: "".join(cols[i][k] for i in range(L)) for k, t in enumerate(taxa)}
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--c12_fasta", required=True, type=str)
    ap.add_argument("--myht_fasta", required=True, type=str)
    ap.add_argument("--contact_map", default=None, type=str, help="Optional TSV of c12_pos myht_pos (1-based).")
    ap.add_argument("--out", default="complementarity.tsv", type=str)
    ap.add_argument("--null", type=int, default=0, help="Number of null shuffles for z-scores (0 disables).")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--max_pairs", type=int, default=20000, help="Safety cap for all-pairs mode.")
    args = ap.parse_args()

    c12 = parse_fasta(Path(args.c12_fasta).expanduser().resolve())
    myh = parse_fasta(Path(args.myht_fasta).expanduser().resolve())
    taxa = sorted(set(c12.keys()) & set(myh.keys()))
    if len(taxa) < 3:
        raise SystemExit(f"Too few shared taxa ({len(taxa)}).")

    c12c = {t: chem_seq(c12[t]) for t in taxa}
    myhc = {t: chem_seq(myh[t]) for t in taxa}
    L1 = len(next(iter(c12c.values())))
    L2 = len(next(iter(myhc.values())))

    if args.contact_map:
        pairs = read_contact_map(Path(args.contact_map).expanduser().resolve())
    else:
        pairs = [(i+1, j+1) for i in range(L1) for j in range(L2)]
        if len(pairs) > args.max_pairs:
            raise SystemExit(f"All-pairs would be {len(pairs)} pairs; raise --max_pairs or provide --contact_map.")

    rng = random.Random(args.seed)

    # Null model: shuffle columns within each alignment independently.
    null_scores: Dict[Tuple[int,int], List[float]] = {p: [] for p in pairs} if args.null > 0 else {}

    if args.null > 0:
        for _ in range(args.null):
            c12_sh = shuffle_columns(c12c, rng)
            myh_sh = shuffle_columns(myhc, rng)
            for p in pairs:
                s = complement_score_for_pair(c12_sh, myh_sh, p[0], p[1], taxa)
                null_scores[p].append(s)

    outp = Path(args.out).resolve()
    with outp.open("w", encoding="utf-8") as out:
        header = ["c12_pos","myht_pos","score","shared_taxa"]
        if args.null > 0:
            header += ["null_mean","null_sd","z"]
        out.write("\t".join(header) + "\n")

        for p in pairs:
            s = complement_score_for_pair(c12c, myhc, p[0], p[1], taxa)
            row = [p[0], p[1], f"{s:.4f}", len(taxa)]
            if args.null > 0:
                ns = null_scores[p]
                mu = sum(ns)/len(ns) if ns else 0.0
                sd = math.sqrt(sum((x-mu)**2 for x in ns)/len(ns)) if ns else 0.0
                z = (s - mu) / sd if sd > 1e-12 else 0.0
                row += [f"{mu:.4f}", f"{sd:.4f}", f"{z:.3f}"]
            out.write("\t".join(map(str,row)) + "\n")

    print(f"[OK] Wrote: {outp}")
    if args.null > 0:
        print(f"[OK] Included null z-scores with {args.null} shuffles (seed={args.seed}).")

if __name__ == "__main__":
    main()
