#!/usr/bin/env python3
"""
prco_decode_cli.py

Per-residue contact occupancy ("PRCO") decoder for a C12–MyhT complex.

Given a topology + trajectory and two MDAnalysis selection strings (C12 and MyhT),
computes:
  - per-C12 residue occupancy: fraction of frames where residue contacts MyhT
  - top partner residues on MyhT per C12 residue (by occupancy)

Outputs:
  <out>_prco.csv
  <out>_prco_topN.txt

Notes:
- Topology and trajectory MUST have the same atom count/order.
- Prefer protein-only topology + protein-only trajectory for speed.
- Uses heavy atoms by default (excludes H*).
"""

import argparse
import sys
import numpy as np
import pandas as pd

import MDAnalysis as mda
from MDAnalysis.lib.distances import capped_distance

try:
    from tqdm import tqdm
except Exception:
    tqdm = None


def fmt_res(res) -> str:
    return f"{res.resname}{int(res.resid)}"


def parse_args():
    p = argparse.ArgumentParser(
        description="Decode PRCO: per-residue contact occupancy + top Myh partners."
    )
    p.add_argument("--top", required=True, help="Topology file (e.g., md.tpr, protein_only.tpr, protein_only.gro)")
    p.add_argument("--traj", required=True, help="Trajectory file (e.g., md.xtc, md_fit.xtc, md_10ps.xtc)")
    p.add_argument("--c12", required=True, help='MDAnalysis selection for C12 (e.g., "protein and resid 1:309")')
    p.add_argument("--myh", required=True, help='MDAnalysis selection for MyhT (e.g., "protein and resid 678:747")')

    p.add_argument("--cutoff", type=float, default=0.50, help="Contact cutoff in nm (default: 0.50)")
    p.add_argument("--stride", type=int, default=1, help="Frame stride (default: 1 = every frame)")
    p.add_argument("--start", type=int, default=0, help="Start frame index (0-based, default: 0)")
    p.add_argument("--stop", type=int, default=None, help="Stop frame index (exclusive, default: end)")
    p.add_argument("--topN", type=int, default=20, help="Top N residues to write to text report (default: 20)")
    p.add_argument("--topK", type=int, default=2, help="Top K Myh partner residues per C12 residue (default: 2)")
    p.add_argument("--include_h", action="store_true", help="Include hydrogens (default: heavy atoms only)")
    p.add_argument("--no_pbc", action="store_true", help="Ignore PBC/box when computing distances (default: use box)")
    p.add_argument("--no_progress", action="store_true", help="Disable progress bar")
    p.add_argument("--out", required=True, help="Output prefix, e.g., eisenia_cutoff_50ns_dt10ps")

    return p.parse_args()


def main():
    args = parse_args()

    # Build Universe
    try:
        u = mda.Universe(args.top, args.traj)
    except Exception as e:
        print(f"[ERROR] Failed to load Universe with top={args.top} traj={args.traj}\n{e}", file=sys.stderr)
        sys.exit(1)

    # Apply selections; enforce heavy atoms by default
    h_clause = "" if args.include_h else " and not name H*"
    c12_sel = f"({args.c12}){h_clause}"
    myh_sel = f"({args.myh}){h_clause}"

    C12 = u.select_atoms(c12_sel)
    MYH = u.select_atoms(myh_sel)

    if C12.n_atoms == 0:
        print(f"[ERROR] C12 selection returned 0 atoms: {c12_sel}", file=sys.stderr)
        sys.exit(1)
    if MYH.n_atoms == 0:
        print(f"[ERROR] Myh selection returned 0 atoms: {myh_sel}", file=sys.stderr)
        sys.exit(1)

    # Residues
    c12_res = list(C12.residues)
    myh_res = list(MYH.residues)
    n_c12 = len(c12_res)
    n_myh = len(myh_res)

    print(f"[INFO] C12 atoms={C12.n_atoms} residues={n_c12} selection={args.c12}")
    print(f"[INFO] Myh atoms={MYH.n_atoms} residues={n_myh} selection={args.myh}")

    # Map atom -> local residue index
    c12_ix_global = np.array([r.ix for r in c12_res], dtype=int)
    myh_ix_global = np.array([r.ix for r in myh_res], dtype=int)
    c12_ix_to_local = {ix: i for i, ix in enumerate(c12_ix_global)}
    myh_ix_to_local = {ix: i for i, ix in enumerate(myh_ix_global)}

    c12_atom_to_local = np.array([c12_ix_to_local[a.residue.ix] for a in C12.atoms], dtype=int)
    myh_atom_to_local = np.array([myh_ix_to_local[a.residue.ix] for a in MYH.atoms], dtype=int)

    # Frame slicing
    start = max(0, args.start)
    stop = args.stop
    stride = max(1, args.stride)

    # Accumulators
    frames = 0
    c12_contact_counts = np.zeros(n_c12, dtype=int)       # any-contact per C12 residue
    pair_counts = np.zeros((n_c12, n_myh), dtype=int)     # residue-residue contact counts

    traj_iter = u.trajectory[start:stop:stride]
    use_progress = (not args.no_progress) and (tqdm is not None)
    if use_progress:
        # MDAnalysis trajectory slicing can be sized; if not, tqdm still works without total
        try:
            total = len(traj_iter)
        except Exception:
            total = None
        traj_iter = tqdm(traj_iter, total=total)

    cutoff_A = args.cutoff * 10.0  # MDAnalysis uses Å for positions

    for ts in traj_iter:
        frames += 1
        box = None if args.no_pbc else ts.dimensions

        pairs = capped_distance(
            C12.atoms.positions,
            MYH.atoms.positions,
            max_cutoff=cutoff_A,
            box=box,
            return_distances=False
        )

        if len(pairs) == 0:
            continue

        # Atom indices within the two arrays -> residue indices (local)
        c12_res_i = c12_atom_to_local[pairs[:, 0]]
        myh_res_j = myh_atom_to_local[pairs[:, 1]]

        # Unique residue-residue pairs per frame (avoid overcounting multiple atom-atom contacts)
        rr = np.unique(np.stack([c12_res_i, myh_res_j], axis=1), axis=0)

        pair_counts[rr[:, 0], rr[:, 1]] += 1
        c12_contact_counts[np.unique(rr[:, 0])] += 1

    if frames == 0:
        print("[ERROR] No frames processed (check --start/--stop/--stride).", file=sys.stderr)
        sys.exit(1)

    c12_occ = c12_contact_counts / frames
    pair_occ = pair_counts / frames

    # Build output table
    topK = max(1, args.topK)
    rows = []
    for i, r in enumerate(c12_res):
        partners = np.argsort(pair_occ[i, :])[::-1]

        row = {
            "C12_resid": int(r.resid),
            "C12_resname": r.resname,
            "occupancy": float(c12_occ[i]),
        }

        for k in range(min(topK, n_myh)):
            j = partners[k]
            row[f"top_partner_{k+1}"] = fmt_res(myh_res[j])
            row[f"partner{k+1}_occ"] = float(pair_occ[i, j])

        rows.append(row)

    df = pd.DataFrame(rows).sort_values("occupancy", ascending=False)
    csv_path = f"{args.out}_prco.csv"
    txt_path = f"{args.out}_prco_top{args.topN}.txt"

    df.to_csv(csv_path, index=False)

    topN_df = df.head(args.topN)
    with open(txt_path, "w") as fh:
        for _, row in topN_df.iterrows():
            lhs = f"{row['C12_resname']}{int(row['C12_resid'])}\tocc={row['occupancy']:.3f}\t"
            rhs = []
            for k in range(min(topK, n_myh)):
                rhs.append(f"{row[f'top_partner_{k+1}']}({row[f'partner{k+1}_occ']:.3f})")
            fh.write(lhs + "\t".join(rhs) + "\n")

    print(f"[INFO] Frames analyzed: {frames}")
    print(f"[INFO] Wrote: {csv_path}")
    print(f"[INFO] Wrote: {txt_path}")


if __name__ == "__main__":
    main()
