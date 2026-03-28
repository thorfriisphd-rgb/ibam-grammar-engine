# stability_plot_script.py

One-stop post-MD stability analysis for IBAM–MyT complexes. Takes GROMACS topology (.tpr) and trajectory (.xtc) files and produces a four-panel Stability Quad figure plus per-residue contact data.

## Output

**Stability Quad** (2×2 composite PNG):
- **A** — RMSD (Cα) vs time
- **B** — RMSF per residue (Cα), with IBAM and MyT blocks colour-coded and MyT interface region shaded
- **C** — Interface contact persistence (fraction of native contacts retained over time)
- **D** — Per-residue contact occupancy across the IBAM–MyT interface

**Data files** (CSV): RMSD time series, RMSF per residue, interface contact persistence, per-residue contact occupancy, residue index map, engaged residues at t₀ and t_end, contact change table (kept/gained/lost).

**Summary** (TXT): Key statistics including final RMSD, mean RMSF, interface retention, and per-chain contact turnover.

## Usage

```bash
python stability_plot_script.py md.tpr md.xtc --species "Homo sapiens"
```

The script will prompt for a species name interactively if `--species` is not provided. The species name appears in the Stability Quad title.

## Key options

| Flag | Default | Description |
|------|---------|-------------|
| `--stride` | 1 | Analyse every Nth frame |
| `--cutoff` | 4.5 | Contact distance cutoff (Å) |
| `--occ-threshold` | 0.0 | Minimum occupancy to include in charts |
| `--species` | *(prompt)* | Species name for figure title |
| `--no-compose-quad` | off | Skip composing the 2×2 Quad PNG |
| `--quad-wspace` | 0.05 | Horizontal gap between Quad columns |
| `--quad-hspace` | 0.05 | Vertical gap between Quad rows |

## Requirements

- Python 3.8+
- MDAnalysis
- NumPy, Pandas, Matplotlib

## Performance

Handles 10 ns, 50 ns, and 100 ns trajectories without difficulty on modest hardware. Use `--stride` to downsample large trajectories if needed.

## Note on chain selection

The script expects GROMACS segid naming: `seg_0_Protein_chain_A` (IBAM) and `seg_1_Protein_chain_B` (MyT). If your topology uses different segment identifiers, adjust `sel_c12` and `sel_myt` near line 214.
