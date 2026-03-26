#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


def load_table(path: Path) -> pd.DataFrame:
    try:
        df = pd.read_csv(path, sep="\t")
        if df.shape[1] == 1:
            df = pd.read_csv(path, sep=r"\s+", engine="python")
    except Exception:
        df = pd.read_csv(path, sep=r"\s+", engine="python")
    return df


def main() -> int:
    ap = argparse.ArgumentParser(description="Plot stacked chemistry barcode from *_chemcomp.tsv")
    ap.add_argument("--infile", required=True, help="Input *_chemcomp.tsv")
    ap.add_argument("--out", required=True, help="Output PNG")
    ap.add_argument(
        "--title",
        default="Panel B — Conserved MG chemical grammar",
        help="Main plot title"
    )
    ap.add_argument(
        "--no-bands",
        action="store_true",
        help="Disable subtle 7-column background shading"
    )
    args = ap.parse_args()

    infile = Path(args.infile)
    if not infile.is_file():
        raise FileNotFoundError(f"Input file not found: {infile}")

    df = load_table(infile)
    required = ["run_id", "gate", "pos", "chem_class", "frac"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns: {missing}; found: {list(df.columns)}")

    df = df.copy()
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce")
    df["frac"] = pd.to_numeric(df["frac"], errors="coerce")
    df = df.dropna(subset=["pos", "frac"])
    df["pos"] = df["pos"].astype(int)
    df["chem_class"] = df["chem_class"].astype(str)

    # Chemistry class order
    class_order = ["H", "P", "B", "A"]
    discovered = [c for c in df["chem_class"].unique() if c not in class_order]
    class_order = [c for c in class_order if c in df["chem_class"].unique()] + sorted(discovered)

    # Color map
    colors = {
        "H": "#555555",  # dark grey
        "P": "#2CA9A1",  # teal
        "B": "#3B6FB6",  # blue
        "A": "#8E44AD",  # purple
    }
    for c in class_order:
        colors.setdefault(c, "#999999")

    mat = (
        df.pivot_table(
            index="pos",
            columns="chem_class",
            values="frac",
            aggfunc="sum",
            fill_value=0.0
        ).sort_index()
    )

    positions = mat.index.tolist()
    if not positions:
        raise ValueError("No valid positions found after parsing input table.")

    fig_w = max(9, 0.72 * len(positions))
    fig_h = 4.8

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    # ------------------------------------------------------------
    # Subtle 7-column background shading to suggest MG periodicity
    # ------------------------------------------------------------
    if not args.no_bands:
        start_pos = min(positions)
        end_pos = max(positions)
        block_idx = 0
        for start in range(start_pos, end_pos + 1, 7):
            if block_idx % 2 == 0:
                ax.axvspan(
                    start - 0.5,
                    min(start + 6.5, end_pos + 0.5),
                    color="0.97",
                    zorder=0
                )
            block_idx += 1

    # ------------------------------------------------------------
    # Determine dominant chemistry class per position
    # ------------------------------------------------------------
    dom = (
        df.sort_values(["pos", "frac"], ascending=[True, False])
          .groupby("pos", as_index=False)
          .first()[["pos", "chem_class"]]
    )
    dom_map = dict(zip(dom["pos"], dom["chem_class"]))

    # ------------------------------------------------------------
    # Plot stacked bars position-by-position so dominant segment
    # can be outlined cleanly
    # ------------------------------------------------------------
    bar_width = 0.68

    for pos in positions:
        bottom = 0.0
        dominant = dom_map.get(pos, None)

        for chem_class in class_order:
            val = float(mat.at[pos, chem_class]) if chem_class in mat.columns else 0.0
            if val <= 0:
                continue

            edgecolor = "black" if chem_class == dominant else "none"
            linewidth = 0.4 if chem_class == dominant else 0.0

            ax.bar(
                pos,
                val,
                bottom=bottom,
                color=colors[chem_class],
                edgecolor=edgecolor,
                linewidth=linewidth,
                width=bar_width,
                zorder=3,
            )
            bottom += val

    # ------------------------------------------------------------
    # Annotate dominant class at the top of each bar
    # ------------------------------------------------------------
    for pos in positions:
        label = dom_map.get(pos, "")
        ax.text(
            pos,
            1.035,
            label,
            ha="center",
            va="bottom",
            fontsize=13,
            fontweight="bold",
            color="black",
            zorder=5,
            clip_on=False
        )

    # ------------------------------------------------------------
    # Axes, labels, title, subtitle
    # ------------------------------------------------------------
    ax.set_ylim(0, 1.08)
    ax.set_xlim(min(positions) - 0.6, max(positions) + 0.6)

    ax.set_xlabel("MG column", fontsize=14)
    ax.set_ylabel("Chemical class fraction", fontsize=14)

    ax.set_xticks(positions)
    ax.set_xticklabels(positions, fontsize=11)

    ax.tick_params(axis="x", width=1.0, length=4, labelsize=11)
    ax.tick_params(axis="y", width=1.0, length=4, labelsize=11)

    ax.set_title(args.title, fontsize=18, pad=14)


    # ------------------------------------------------------------
    # Legend
    # ------------------------------------------------------------
    legend_handles = [
        Patch(facecolor=colors[c], edgecolor="black", label=c)
        for c in class_order
    ]
    ax.legend(
        handles=legend_handles,
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.12),
        ncol=min(len(class_order), 6),
        fontsize=11
    )

    # ------------------------------------------------------------
    # Journal-style axes
    # ------------------------------------------------------------
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.2)
    ax.spines["bottom"].set_linewidth(1.2)

    fig.tight_layout()
    fig.savefig(args.out, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"[OK] Wrote: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
