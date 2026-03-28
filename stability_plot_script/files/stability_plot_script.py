#!/usr/bin/env python3
"""
new_analyze_rms.py — One-stop post-MDS analysis + panel exports + Quad composer.

Outputs:
  <prefix>_RMSD.csv/.png
  <prefix>_RMSF.csv/.png
  <prefix>_InterfaceContacts.csv/.png
  <prefix>_ContactOccupancy.csv
  <prefix>_ResidueIndexMap.csv
  <prefix>_Engaged_t0.csv
  <prefix>_Engaged_tend.csv
  <prefix>_ContactChange.csv
  <prefix>_ContactOccupancy.png   (6x4, slimmer bars, Chain B re-indexed from 1)
  <prefix>_StabilityQuad.png      (12x8 composed from four PNGs)
  <prefix>_Summary.txt
"""

import os, argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from pathlib import Path
from datetime import datetime
import MDAnalysis as mda
from MDAnalysis.analysis import rms
from MDAnalysis.lib.distances import distance_array

# ---------------------------
# Helpers
# ---------------------------
def prompt_species():
    while True:
        name = input("> Enter species name: ").strip()
        confirm = input(f"You entered '{name}'. Is this correct? (1 = Yes, 0 = No): ").strip()
        if confirm == "1":
            return name
        elif confirm == "0":
            print("Try again.")
        else:
            print("Please enter 1 or 0.")


def render_occupancy_png_grouped(occ_df: pd.DataFrame, out_png: str,
                                 bar_width: float = 0.9, tick_target: int = 12,
                                 title: str = "Per-Residue Contact Occupancy"):
    """
    Render occupancy panel at 6x4 inches with two indexed blocks:
      - Chain A (C12, seg_0_*): x labels 1..N_A
      - Chain B (MyT): x labels 1..N_B (restart at 1), plotted after Chain A
    Visual divider drawn between blocks.
    """
    # Split by chain (preserve original order)
    is_c12 = occ_df["ChainID"].astype(str).str.startswith("seg_0_")
    dfA = occ_df[is_c12].copy()
    dfB = occ_df[~is_c12].copy()

    # Build plotting positions
    nA, nB = len(dfA), len(dfB)
    xA = np.arange(nA)          # 0..nA-1
    xB = nA + np.arange(nB)     # nA..nA+nB-1

    # Labels: restart from 1 in each block
    labelsA = [str(i) for i in range(1, nA + 1)]
    labelsB = [str(i) for i in range(1, nB + 1)]

    fig, ax = plt.subplots(figsize=(6, 4))

    # (1) Center the panel title
    ax.set_title(title)  # default loc="center"

    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Contact Occupancy")
    ax.set_xlabel("Residues (IBAM block → MyT block)")

    # Bars (width will be auto-guarded below)
    if nA > 0:
        ax.bar(xA, dfA["Occupancy"].values, color="#1f77b4", width=bar_width, label="IBAM")
    if nB > 0:
        ax.bar(xB, dfB["Occupancy"].values, color="#ff9999", width=bar_width, label="MyT")

    # X ticks (sparse) – aim for ~tick_target labels across the whole axis
    N = nA + nB
    if N > 0:
        tick_every = max(1, N // max(1, int(tick_target)))
        ticks = np.arange(0, N, tick_every)

        def label_for_pos(pos):
            if pos < nA:
                return labelsA[pos]
            else:
                idx = pos - nA
                return labelsB[idx] if idx < nB else ""

        xticklabels = [label_for_pos(int(t)) for t in ticks]
        ax.set_xticks(ticks)
        ax.set_xticklabels(xticklabels, rotation=60, ha="right", fontsize=8)
        ax.tick_params(axis="x", pad=2)

    # Vertical divider between blocks (if both present)
    if nA > 0 and nB > 0:
        ax.axvline(x=nA - 0.5, color="0.7", linestyle="--", linewidth=1)

        # Colored labels serving as the legend (no legend box)
        ax.text(0.04, 1.02, "IBAM", ha="left",  va="bottom", transform=ax.transAxes, color="#1f77b4",fontweight="bold")
        ax.text(0.96, 1.02, "MyT",      ha="right", va="bottom", transform=ax.transAxes, color="#ff9999", fontweight="bold")

    # No legend box inside the plot
    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f">>>>>>>>> Wrote occupancy panel PNG: {out_png}")

def compose_stability_quad_from_pngs(
    rmsd_png: str, rmsf_png: str, contacts_png: str, occupancy_png: str,
    out_png: str, suptitle: str, wspace: float = 0.05, hspace: float = 0.05,
    title_y: float = 0.96, top_rect: float = 0.95
):
    for pth in (rmsd_png, rmsf_png, contacts_png, occupancy_png):
        if not Path(pth).is_file():
            raise FileNotFoundError(f"Missing panel PNG: {pth}")

    # Manual layout (don’t use constrained_layout here)
    fig, axs = plt.subplots(
        2, 2, figsize=(12, 8),
        gridspec_kw={"wspace": wspace, "hspace": hspace}
    )

    panels = [
        (axs[0, 0], rmsd_png,     "A  RMSD (Cα)"),
        (axs[0, 1], rmsf_png,     "B  RMSF per Residue (Cα)"),
        (axs[1, 0], contacts_png, "C  Interface Contact Persistence (IBAM ↔ MyT)"),
        (axs[1, 1], occupancy_png,"D  Per-Residue Contact Occupancy"),
    ]
    for ax, png, ttl in panels:
        ax.imshow(mpimg.imread(png))
        ax.set_title(ttl, loc="left")
        ax.axis("off")

    fig.suptitle(suptitle, fontsize=16, y=title_y)
    fig.tight_layout(rect=[0, 0, 1, top_rect])

    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f">>>>>>>>>> Wrote Stability Quad: {out_png}")

# ---------------------------
# CLI
# ---------------------------
parser = argparse.ArgumentParser(description="Post-MDS analysis + Quad composer")
parser.add_argument("topology_file", nargs="?", default="md.tpr")
parser.add_argument("trajectory_file", nargs="?", default="md.xtc")
parser.add_argument("--stride", type=int, default=1, help="analyze every Nth frame")
parser.add_argument("--species", type=str, default=None, help="species name (skip prompt)")
parser.add_argument("--cutoff", type=float, default=4.5, help="contact cutoff Å (default 4.5)")
parser.add_argument("--occ-threshold", type=float, default=0.0, help="min occupancy (0–1) to include in charts")
parser.add_argument("--bar-width", type=float, default=0.9, help="occupancy bar width (default 0.9)")
parser.add_argument("--tick-target", type=int, default=12, help="~number of x tick labels")
parser.add_argument("--no-compose-quad", dest="no_compose_quad", action="store_true",
                    help="skip composing the 2x2 Quad from PNGs")
parser.add_argument("--quad-wspace", type=float, default=0.05, help="horizontal gap between columns")
parser.add_argument("--quad-hspace", type=float, default=0.05, help="vertical gap between rows")
parser.add_argument("--quad-title-y", type=float, default=0.96, help="suptitle y position (0–1)")
parser.add_argument("--quad-top", type=float, default=0.95, help="top margin for tight_layout rect")
args = parser.parse_args()

topology_file   = args.topology_file
trajectory_file = args.trajectory_file
STRIDE          = max(1, int(args.stride))
CUTOFF          = float(args.cutoff)
OCC_THRESHOLD   = float(max(0.0, min(1.0, args.occ_threshold)))
species_name    = args.species or prompt_species()

if not os.path.isfile(topology_file):
    raise FileNotFoundError(f"Topology not found: {topology_file}")
if not os.path.isfile(trajectory_file):
    raise FileNotFoundError(f"Trajectory not found: {trajectory_file}")

output_prefix = os.path.splitext(os.path.basename(trajectory_file))[0]

print(f">> Loading topology:  {topology_file}")
print(f">>> Loading trajectory: {trajectory_file}  (stride={STRIDE})")
u = mda.Universe(topology_file, trajectory_file)

sel_ca = "protein and name CA"
protein_ca = u.select_atoms(sel_ca)
assert protein_ca.n_atoms > 0, "No CA atoms found in selection."


# ---------------------------
# RMSD (Cα) — simple superposition to frame 0
# ---------------------------
print(">>>> Calculating RMSD (Cα)…")
R = rms.RMSD(protein_ca, reference=protein_ca, select=sel_ca, ref_frame=0)
R.run(step=STRIDE)
r = R.results.rmsd
rmsd_df = pd.DataFrame(r, columns=["Frame", "Time(ps)", "RMSD(Å)"])
rmsd_df.to_csv(f"{output_prefix}_RMSD.csv", index=False)

plt.figure(figsize=(6,4))
plt.plot(rmsd_df["Time(ps)"], rmsd_df["RMSD(Å)"], lw=1.5)
plt.xlabel("Time (ps)"); plt.ylabel("RMSD (Å)"); plt.title("RMSD vs Time (Cα)")
plt.grid(True); plt.tight_layout()
plt.savefig(f"{output_prefix}_RMSD.png", dpi=300, bbox_inches="tight"); plt.close()


# ---------------------------
# Interface native-contact persistence & occupancy
# ---------------------------
print(">>>>> Calculating interface contacts…")

# NOTE: adjust segids if your topology names differ
sel_c12 = "protein and segid seg_0_Protein_chain_A and not name H*"
sel_myt = "protein and segid seg_1_Protein_chain_B and not name H*"
ag_c12 = u.select_atoms(sel_c12)
ag_myt = u.select_atoms(sel_myt)

if ag_c12.n_atoms == 0 or ag_myt.n_atoms == 0:
    seginfo = [(str(seg.segid), seg.atoms.n_atoms) for seg in u.segments]
    print("!! Selection returned zero atoms. Available segments (segid, n_atoms):")
    for s, n in seginfo: print(f"   - {s}: {n}")
    raise ValueError("Check segid names in sel_c12 / sel_myt above.")

def residue_info(ag):
    return [(res.resid, res.resname, res.segid) for res in ag.residues]

res_info_c12 = residue_info(ag_c12)
res_info_myt = residue_info(ag_myt)

# Map (resid, segid) -> resname
resname_map = {(int(resid), str(segid)): str(resname)
               for resid, resname, segid in (res_info_c12 + res_info_myt)}

# Atom→residue arrays
a2resid_c12 = ag_c12.atoms.resids; a2segid_c12 = ag_c12.atoms.segids
a2resid_myt = ag_myt.atoms.resids; a2segid_myt = ag_myt.atoms.segids

# Native contacts at frame 0
u.trajectory[0]
D0 = distance_array(ag_c12.positions, ag_myt.positions)
i0, j0 = np.where(D0 <= CUTOFF)
native_pairs = set(zip(i0, j0))

# Time series + occupancy tracking
times, frac_retained = [], []
c12_occ_counts = {(int(resid), str(segid)): 0 for resid, _, segid in res_info_c12}
myt_occ_counts = {(int(resid), str(segid)): 0 for resid, _, segid in res_info_myt}
initial_contact_keys, final_contact_keys = set(), set()

n_total_frames = len(u.trajectory)
n_strided = (n_total_frames + STRIDE - 1) // STRIDE

for frame_idx, ts in enumerate(u.trajectory[::STRIDE]):
    D = distance_array(ag_c12.positions, ag_myt.positions)
    ci, cj = np.where(D <= CUTOFF)
    current_pairs = set(zip(ci, cj))
    retained = len(native_pairs & current_pairs)
    frac = retained / max(1, len(native_pairs))
    times.append(ts.time); frac_retained.append(frac)

    c12_keys_now = set(zip(a2resid_c12[ci], a2segid_c12[ci]))
    myt_keys_now = set(zip(a2resid_myt[cj], a2segid_myt[cj]))
    for k in c12_keys_now: c12_occ_counts[(int(k[0]), str(k[1]))] += 1
    for k in myt_keys_now: myt_occ_counts[(int(k[0]), str(k[1]))] += 1

    if frame_idx == 0: initial_contact_keys = {(int(r), str(s)) for r, s in (c12_keys_now | myt_keys_now)}
    if frame_idx == n_strided - 1: final_contact_keys = {(int(r), str(s)) for r, s in (c12_keys_now | myt_keys_now)}

# Contact-persistence timeseries
contacts_df = pd.DataFrame({"Time(ps)": times, "FracNativeContacts": frac_retained})
contacts_df.to_csv(f"{output_prefix}_InterfaceContacts.csv", index=False)

plt.figure(figsize=(6,4))
plt.plot(contacts_df["Time(ps)"], contacts_df["FracNativeContacts"], lw=1.5)
plt.xlabel("Time (ps)"); plt.ylabel("Fraction of Native Contacts")
plt.ylim(0, 1.05); plt.title("Interface Contact Persistence"); plt.grid(True); plt.tight_layout()
plt.savefig(f"{output_prefix}_InterfaceContacts.png", dpi=300, bbox_inches="tight"); plt.close()

# Per-residue occupancy CSV (C12 then MyT order), thresholded
total_frames_contacts = len(times)
seen, ordered = set(), []
for resid, resname, segid in (res_info_c12 + res_info_myt):
    key = (int(resid), str(segid))
    if key not in seen:
        seen.add(key); ordered.append((int(resid), str(resname), str(segid)))

rows = []
for resid, resname, segid in ordered:
    key = (resid, segid)
    if segid.startswith("seg_0_"):
        occ = c12_occ_counts.get(key, 0) / max(1, total_frames_contacts)
    else:
        occ = myt_occ_counts.get(key, 0) / max(1, total_frames_contacts)
    if occ >= OCC_THRESHOLD:
        rows.append([resid, resname, segid, float(occ)])

occ_df = pd.DataFrame(rows, columns=["Resid", "ResName", "ChainID", "Occupancy"])
occ_df.to_csv(f"{output_prefix}_ContactOccupancy.csv", index=False)

# ResidueIndex map
index_map = occ_df.reset_index().rename(columns={"index": "BarIndex"})
index_map.to_csv(f"{output_prefix}_ResidueIndexMap.csv", index=False)

# t0 / tend engagement
def _keys_to_rows(keys):
    return [(int(r), resname_map.get((int(r), str(s)), "UNK"), str(s))
            for r, s in sorted(keys, key=lambda k: (k[1], k[0]))]

pd.DataFrame(_keys_to_rows(initial_contact_keys), columns=["Resid","ResName","ChainID"]).to_csv(
    f"{output_prefix}_Engaged_t0.csv", index=False)
pd.DataFrame(_keys_to_rows(final_contact_keys), columns=["Resid","ResName","ChainID"]).to_csv(
    f"{output_prefix}_Engaged_tend.csv", index=False)

# Contact change table
change_rows = []
all_keys = initial_contact_keys | final_contact_keys
for resid, segid in sorted(all_keys, key=lambda k: (k[1], k[0])):
    resname = resname_map.get((int(resid), str(segid)), "UNK")
    init = int((resid, segid) in initial_contact_keys)
    fin  = int((resid, segid) in final_contact_keys)
    status = ("kept" if init and fin else
              "lost" if init and not fin else
              "gained" if not init and fin else
              "none")
    change_rows.append([int(resid), str(resname), str(segid), init, fin, status])

change_df = pd.DataFrame(
    change_rows, columns=["Resid", "ResName", "ChainID", "InitialContact", "FinalContact", "Status"])
change_df.to_csv(f"{output_prefix}_ContactChange.csv", index=False)

print(">>>>>> Contacts: time series, occupancy, maps, and change tables written.")


# ---------------------------
# RMSF (Cα) — with MyT interface shading + stats
# ---------------------------
print(">>>>>>> Calculating RMSF (Cα)…")
_ = rms.RMSD(u, select=sel_ca, ref_frame=0).run()  # ensure superposition per frame
coords_list = []
for ts in u.trajectory[::STRIDE]:
    coords_list.append(protein_ca.positions.copy())
coords = np.asarray(coords_list, dtype=np.float32)   # (frames, atoms, 3)
mean_coords = coords.mean(axis=0)
sq_fluct = ((coords - mean_coords) ** 2).sum(axis=2)
rmsf_vals = np.sqrt(sq_fluct.mean(axis=0))

rids = [res.resid for res in protein_ca.residues]
rmsf_df = pd.DataFrame({"Residue": rids, "RMSF(Å)": rmsf_vals})
rmsf_df.to_csv(f"{output_prefix}_RMSF.csv", index=False)

# Robust split by segid (works even if interleaved)
res_segs = [str(res.segid) for res in protein_ca.residues]
idxA = [i for i, s in enumerate(res_segs) if s.startswith("seg_0_")]  # C12
idxB = [i for i, s in enumerate(res_segs) if s.startswith("seg_1_")]  # MyT
xA  = np.arange(len(idxA))
xB  = len(idxA) + np.arange(len(idxB))

fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(xA, rmsf_vals[idxA], lw=1.5, label="IBAM", color="#1f77b4")
ax.plot(xB, rmsf_vals[idxB], lw=1.5, label="MyT",        color="#ff9999")

# Shade MyT interface window using occ_df (occ ≥ 0.5)
iface_idx = np.array([], dtype=int)
try:
    dfB_occ = occ_df[~occ_df["ChainID"].astype(str).str.startswith("seg_0_")].reset_index(drop=True)
    occ_thresh = 0.5
    iface_idx = np.where(dfB_occ["Occupancy"].values >= occ_thresh)[0]
    if iface_idx.size:
        x0 = len(idxA) + iface_idx.min() - 0.5
        x1 = len(idxA) + iface_idx.max() + 0.5
        ax.axvspan(x0, x1, alpha=0.12, linewidth=0)
except Exception as e:
    print(f"[WARN] Could not shade MyT interface in RMSF: {e}")

# Labels, divider (no legend box)
if len(idxA) and len(idxB):
    ax.axvline(x=len(idxA)-0.5, color="0.7", linestyle="--", linewidth=1)
    ax.text(0.04, 1.02, "IBAM", ha="left",  va="bottom", transform=ax.transAxes, color="#1f77b4", fontweight="bold")
    ax.text(0.96, 1.02, "MyT",      ha="right", va="bottom", transform=ax.transAxes, color="#ff9999", fontweight="bold")

ax.set_xlabel("Residues (IBAM block → MyT block)")
ax.set_ylabel("RMSF (Å)")
ax.set_title("RMSF per Residue (Cα)")
ax.grid(True)

# Sparse x-ticks (~12 total)
N = len(idxA) + len(idxB)
if N:
    tick_every = max(1, N // 12)
    ticks = np.arange(0, N, tick_every)
    labA = [str(i) for i in range(1, len(idxA)+1)]
    labB = [str(i) for i in range(1, len(idxB)+1)]
    def lab(pos): return labA[pos] if pos < len(idxA) else labB[pos - len(idxA)]
    ax.set_xticks(ticks)
    ax.set_xticklabels([lab(int(t)) for t in ticks], rotation=60, ha="right", fontsize=8)

# No legend box
fig.tight_layout()
fig.savefig(f"{output_prefix}_RMSF.png", dpi=300, bbox_inches="tight")
plt.close(fig)

# Quick stats for captions (printed and appended to Summary)
c12_mean = float(np.mean(rmsf_vals[idxA])) if len(idxA) else float("nan")
iface_mean = tail_mean = fold = float("nan")
if iface_idx.size:
    iface_mean = float(np.mean(rmsf_vals[idxB][np.array(iface_idx)]))
    tail_mask  = np.arange(len(idxB)) > iface_idx.max()
    if tail_mask.any():
        tail_mean = float(np.mean(rmsf_vals[idxB][tail_mask]))
        fold = tail_mean / iface_mean if iface_mean > 0 else float("nan")
rmsf_msg = (f"[RMSF] IBAM mean = {c12_mean:.2f} Å; "
            f"MyT interface mean = {iface_mean:.2f} Å; "
            f"MyT distal tail mean = {tail_mean:.2f} Å; "
            f"tail/interface = {fold:.1f}×")


# ---------------------------
# Standalone Occupancy PNG (6x4) + Quad
# ---------------------------
print(">>>>>>>> Rendering Occupancy (6x4) and composing Stability Quad…")
occ_png_path = f"{output_prefix}_ContactOccupancy.png"
render_occupancy_png_grouped(
    occ_df, out_png=occ_png_path,
    bar_width=float(args.bar_width),
    tick_target=int(args.tick_target),
    title="Per-Residue Contact Occupancy"
)

quad_out_path = f"{output_prefix}_StabilityQuad.png"
if not args.no_compose_quad:
    compose_stability_quad_from_pngs(
    rmsd_png=f"{output_prefix}_RMSD.png",
    rmsf_png=f"{output_prefix}_RMSF.png",
    contacts_png=f"{output_prefix}_InterfaceContacts.png",
    occupancy_png=occ_png_path,
    out_png=quad_out_path,
    suptitle=f"Stability Quad — {species_name}",
    wspace=float(args.quad_wspace),
    hspace=float(args.quad_hspace),
    title_y=float(args.quad_title_y),
    top_rect=float(args.quad_top),
)
else:
    print("[INFO] Skipped composing Quad (--no-compose-quad).")

# ---------------------------
# Summary file
# ---------------------------
print(">>>>>>>>>>> Writing summary…")
summary_path = f"{output_prefix}_Summary.txt"

frames_analyzed = int(len(rmsd_df))
sim_time_ps     = float(rmsd_df["Time(ps)"].iloc[-1])
final_rmsd      = float(rmsd_df["RMSD(Å)"].iloc[-1])
mean_rmsf       = float(rmsf_df["RMSF(Å)"].mean())
last_q          = max(1, int(0.25 * len(contacts_df)))
contacts_mean_last_q = 100.0 * float(contacts_df["FracNativeContacts"].iloc[-last_q:].mean())

def _res_stats(df: pd.DataFrame):
    if df is None or len(df) == 0: return 0,0,0,0,0.0
    kept   = int(((df["InitialContact"] == 1) & (df["FinalContact"] == 1)).sum())
    gained = int(((df["InitialContact"] == 0) & (df["FinalContact"] == 1)).sum())
    lost   = int(((df["InitialContact"] == 1) & (df["FinalContact"] == 0)).sum())
    total  = int(len(df)); kept_pct = 100.0 * kept / max(1, total)
    return kept, gained, lost, total, kept_pct

overall_kept, overall_gained, overall_lost, overall_total, overall_kept_pct = _res_stats(change_df)
_chain = change_df.get("ChainID", pd.Series(dtype=str)).astype(str)
c12_df = change_df[_chain.str.startswith("seg_0_")]
myt_df = change_df[~_chain.str.startswith("seg_0_")]
c12_kept, c12_gained, c12_lost, c12_total, c12_kept_pct = _res_stats(c12_df)
myt_kept, myt_gained, myt_lost, myt_total, myt_kept_pct = _res_stats(myt_df)

run_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
with open(summary_path, "w") as f:
    f.write(
        f"IBAM–MyT stability analysis — {species_name} (AF3 → explicit-solvent MD)\n"
        f"Run date: {run_date}\n"
        "------------------------------------------------------\n"
        f"Simulation files: {topology_file} + {trajectory_file}\n"
        f"Frames analyzed:  {frames_analyzed}  (stride={STRIDE})\n"
        f"Sim time (ps):    {sim_time_ps:.1f}\n"
        f"RMSD (Cα) final:  {final_rmsd:.2f} Å\n"
        f"RMSF (Cα) mean:   {mean_rmsf:.2f} Å\n"
        f"Interface native-contacts retained (last quartile): {contacts_mean_last_q:.1f}%\n"
        f"Contact cutoff (Å): {CUTOFF:.2f}\n"
        f"Occupancy plot threshold: {OCC_THRESHOLD:.2f}\n"
        f"Residues kept (t0→tend): {overall_kept} / {overall_total} ({overall_kept_pct:.1f}%)\n"
        f"Residues gained/lost: {overall_gained} / {overall_lost}\n"
        f"IBAM kept (t0→tend): {c12_kept} / {c12_total} ({c12_kept_pct:.1f}%)   |   gained/lost: {c12_gained} / {c12_lost}\n"
        f"MyT kept (t0→tend): {myt_kept} / {myt_total} ({myt_kept_pct:.1f}%)   |   gained/lost: {myt_gained} / {myt_lost}\n"
        f"{rmsf_msg}\n"
    )
print(f">>>>>>>>>>>> Summary: {summary_path}")
print(">>>>>>>>>>>>> Done.")