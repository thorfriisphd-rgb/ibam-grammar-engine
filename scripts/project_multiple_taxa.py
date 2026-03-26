#!/usr/bin/env python3
import argparse, os, sys, math
from collections import Counter, defaultdict

# ----------------------------
# I/O helpers
# ----------------------------
def read_fasta(path):
    seqs = {}
    h = None
    buf = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if h is not None:
                    seqs[h] = "".join(buf)
                h = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if h is not None:
            seqs[h] = "".join(buf)
    return seqs

def read_watchlist(path):
    hits = set()
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            try:
                hits.add(int(line))
            except ValueError:
                continue
    return hits

# ----------------------------
# Chemistry bins
# ----------------------------
AA_CLASS = {
    # Basic
    "K":"BASIC","R":"BASIC","H":"BASIC",
    # Acidic
    "D":"ACIDIC","E":"ACIDIC",
    # Polar uncharged
    "S":"POLAR","T":"POLAR","N":"POLAR","Q":"POLAR",
    # Hydrophobic
    "A":"HYDRO","V":"HYDRO","I":"HYDRO","L":"HYDRO","M":"HYDRO",
    # Aromatic
    "F":"AROM","Y":"AROM","W":"AROM",
    # Specials
    "G":"SPECIAL","P":"SPECIAL","C":"SPECIAL",
}

def aa_to_class(aa: str) -> str:
    aa = aa.upper()
    if aa in ("-", "X", "?", "*"):
        return "GAP"
    return AA_CLASS.get(aa, "UNK")

def csv_join(items):
    # Stable, readable output
    return ",".join(items)

# ----------------------------
# Main
# ----------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Project watchlist residue positions onto alignment columns, "
                    "trim columns by taxa support and optional chemistry-consensus gate."
    )
    ap.add_argument("--aln", required=True, help="Input alignment FASTA (all taxa).")
    ap.add_argument("--watchdir", required=True, help="Directory containing watchlists.")
    ap.add_argument("--watch-suffix", default="_watch.txt", help="Watchlist filename suffix (default _watch.txt).")

    # Gate 1: coverage threshold
    ap.add_argument("--min-taxa-per-col", type=int, default=None,
                    help="Gate 1: minimum number of taxa supporting a column (watchlist-hit taxa).")
    ap.add_argument("--min-taxa-frac", type=float, default=None,
                    help="Gate 1 alternative: fraction of total taxa, e.g. 0.50 for >=50%%. "
                         "If provided, overrides --min-taxa-per-col.")

    # Gate 2: chemistry consensus among supporting taxa at that column
    ap.add_argument("--chem-thr", type=float, default=0.85,
                    help="Gate 2: chemistry consensus threshold among supporting taxa (default 0.85).")
    ap.add_argument("--disable-chem-gate", action="store_true",
                    help="Disable chemistry gate (keep only Gate 1).")

    # Outputs
    ap.add_argument("--out-fa", required=True, help="Output trimmed FASTA (projected alignment).")
    ap.add_argument("--out-map", required=True, help="Output mapping TSV with per-column stats.")
    ap.add_argument("--out-taxon-report", default=None,
                    help="Optional: write per-taxon minority/consensus statistics TSV (outlier detection aid).")

    args = ap.parse_args()

    aln = read_fasta(args.aln)
    headers = list(aln.keys())
    if not headers:
        sys.exit("ERROR: empty alignment")

    L = len(aln[headers[0]])
    for h in headers:
        if len(aln[h]) != L:
            sys.exit(f"ERROR: alignment length mismatch for {h}")

    n_taxa = len(headers)

    # Gate 1 threshold determination
    #if args.min_talta := False:
    #    pass  # (never used; placeholder to avoid accidental walrus use)
    if args.min_taxa_frac is not None:
        if not (0 < args.min_taxa_frac <= 1.0):
            sys.exit("ERROR: --min-taxa-frac must be in (0,1].")
        min_taxa_per_col = int(math.ceil(args.min_taxa_frac * n_taxa))
    else:
        if args.min_taxa_per_col is None:
            sys.exit("ERROR: must provide either --min-taxa-per-col or --min-taxa-frac")
        min_taxa_per_col = args.min_taxa_per_col

    if not (0.0 < args.chem_thr <= 1.0):
        sys.exit("ERROR: --chem-thr must be in (0,1].")

    # Load watchlists; support binomial->genus fallback
    hit_resids = {}
    for h in headers:
        wl1 = os.path.join(args.watchdir, f"{h}{args.watch_suffix}")
        wl = wl1
        if not os.path.exists(wl):
            genus = h.split("_", 1)[0]
            wl2 = os.path.join(args.watchdir, f"{genus}{args.watch_suffix}")
            if os.path.exists(wl2):
                wl = wl2
            else:
                sys.exit(f"ERROR: missing watchlist for {h}: {wl1} (also tried {wl2})")
        hit_resids[h] = read_watchlist(wl)

    # Map seq positions (1-based, ungapped) to alignment columns for each taxon
    seqpos_to_col = {}
    for h in headers:
        m = {}
        seqpos = 0
        for col, aa in enumerate(aln[h]):
            if aa != "-":
                seqpos += 1
                m[seqpos] = col
        seqpos_to_col[h] = m

    # For each taxon, compute which alignment columns are "hit" by its watchlist
    hit_cols_by_taxon = {}
    for h in headers:
        m = seqpos_to_col[h]
        hit_cols = set()
        for r in hit_resids[h]:
            if r in m:
                hit_cols.add(m[r])
        hit_cols_by_taxon[h] = hit_cols

    # Compute support per alignment column (Gate 1)
    support = [0] * L
    supporters = [[] for _ in range(L)]  # list of taxa supporting each column
    for h in headers:
        for col in hit_cols_by_taxon[h]:
            support[col] += 1
            supporters[col].append(h)

    gate1_cols = [i for i, s in enumerate(support) if s >= min_taxa_per_col]
    if not gate1_cols:
        sys.exit("ERROR: no columns passed Gate 1 (taxa-per-column threshold)")

    kept_cols = []
    col_stats = {}  # col -> dict

    # Per-taxon report counters (based on Gate1 columns only)
    tax_counts = {h: {"support_cols": 0, "majority_class": 0, "minority_class": 0} for h in headers}

    for col in gate1_cols:
        supp_taxa = sorted(supporters[col])
        residues = [aln[h][col] for h in supp_taxa]
        res_counts = Counter(residues)

        classes = [aa_to_class(a) for a in residues]
        class_counts = Counter(classes)

        dom_res, dom_res_n = res_counts.most_common(1)[0]
        dom_cls, dom_cls_n = class_counts.most_common(1)[0]
        n = len(residues)
        dom_res_frac = dom_res_n / n if n else 0.0
        dom_cls_frac = dom_cls_n / n if n else 0.0

        # Identify minority taxa among supporters (chemistry != dominant class)
        minority_taxa = []
        class_by_taxon = {}
        for h in supp_taxa:
            cls = aa_to_class(aln[h][col])
            class_by_taxon[h] = cls
            if cls != dom_cls:
                minority_taxa.append(h)

        # Update per-taxon counters (Gate1 columns only)
        for h in supp_taxa:
            tax_counts[h]["support_cols"] += 1
            if class_by_taxon[h] == dom_cls:
                tax_counts[h]["majority_class"] += 1
            else:
                tax_counts[h]["minority_class"] += 1

        passes_chem = True
        if not args.disable_chem_gate:
            # If dominant class is GAP/UNK, treat as fail (should be rare)
            if dom_cls in ("GAP", "UNK"):
                passes_chem = False
            else:
                passes_chem = (dom_cls_frac >= args.chem_thr)

        if passes_chem:
            kept_cols.append(col)

        col_stats[col] = {
            "support_taxa": support[col],
            "supporters": supp_taxa,
            "n_support": n,
            "dom_res": dom_res,
            "dom_res_frac": dom_res_frac,
            "dom_class": dom_cls,
            "dom_class_frac": dom_cls_frac,
            "chem_gate_pass": passes_chem,
            "minority_taxa": minority_taxa,
        }

    if not kept_cols:
        sys.exit("ERROR: columns passed Gate 1 but none passed chemistry gate")

    # Write trimmed alignment
    with open(args.out_fa, "w") as out:
        for h in headers:
            trimmed = "".join(aln[h][i] for i in kept_cols)
            out.write(f">{h}\n")
            for i in range(0, len(trimmed), 80):
                out.write(trimmed[i:i+80] + "\n")

    # Write mapping table (verbose + supporters/minority taxa lists)
    with open(args.out_map, "w") as out:
        out.write("\t".join([
            "aln_col_1based",
            "support_taxa",
            "n_support_taxa",
            "dom_res",
            "dom_res_frac",
            "dom_class",
            "dom_class_frac",
            "chem_gate_pass",
            "supporters",
            "minority_taxa",
        ]) + "\n")

        for col in kept_cols:
            st = col_stats[col]
            out.write("\t".join([
                str(col + 1),
                str(st["support_taxa"]),
                str(st["n_support"]),
                st["dom_res"],
                f"{st['dom_res_frac']:.3f}",
                st["dom_class"],
                f"{st['dom_class_frac']:.3f}",
                "1" if st["chem_gate_pass"] else "0",
                csv_join(st["supporters"]),
                csv_join(st["minority_taxa"]),
            ]) + "\n")

    # Optional: per-taxon outlier report
    if args.out_taxon_report:
        with open(args.out_taxon_report, "w") as out:
            out.write("\t".join([
                "taxon",
                "support_cols_gate1",
                "majority_class_cols",
                "minority_class_cols",
                "minority_frac",
            ]) + "\n")

            # Sort by most "radical" first: minority fraction descending
            def minority_frac(h):
                sc = tax_counts[h]["support_cols"]
                return (tax_counts[h]["minority_class"] / sc) if sc else 0.0

            for h in sorted(headers, key=minority_frac, reverse=True):
                sc = tax_counts[h]["support_cols"]
                maj = tax_counts[h]["majority_class"]
                mn = tax_counts[h]["minority_class"]
                frac = (mn / sc) if sc else 0.0
                out.write("\t".join([
                    h,
                    str(sc),
                    str(maj),
                    str(mn),
                    f"{frac:.3f}",
                ]) + "\n")

    # Console summary
    print(f"Wrote: {args.out_fa}")
    print(f"Wrote: {args.out_map}")
    if args.out_taxon_report:
        print(f"Wrote: {args.out_taxon_report}")
    print(f"Taxa: {n_taxa}")
    if args.min_taxa_frac is not None:
        print(f"Gate 1: min taxa/col = {min_taxa_per_col} (frac={args.min_taxa_frac})")
    else:
        print(f"Gate 1: min taxa/col = {min_taxa_per_col}")
    if args.disable_chem_gate:
        print("Gate 2: chemistry gate DISABLED")
    else:
        print(f"Gate 2: chemistry consensus >= {args.chem_thr:.2f} among supporting taxa")
    print(f"Kept {len(kept_cols)} / {L} alignment columns "
          f"(passed Gate1={len(gate1_cols)}, passed Gate2={len(kept_cols)}).")

if __name__ == "__main__":
    main()
