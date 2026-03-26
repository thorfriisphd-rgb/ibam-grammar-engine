# IBAM Grammar Engine

**C12orf29/IBAM: PRCO → Projection → WebLogo → Chemical Barcode**

Single-command pipeline wrapping the battle-tested PDWLS and BARWLS
scripts into a reproducible, manifest-tracked workflow.

---

## What it produces

From MD simulation data across multiple taxa:
- PRCO contact decoding per taxon
- SC watchlists (structurally recurrent contact residues)
- MAFFT alignment
- Dual-gate projected FASTA alignments (3D → 1D projection)
- WebLogo sequence logos (per gate threshold)
- Chemical barcode analysis (optional)
- Full manifest with SHA256 checksums and tool versions

---

## Directory layout

```
ibam_pipeline/
├── run_pipeline.sh              ← Single entry point
├── README.md
├── scripts/
│   ├── PDWLS.sh                 ← PRCO → projection → WebLogo
│   ├── prco_decode_cli.py       ← Per-residue contact occupancy
│   ├── project_multiple_taxa.py ← Watchlist → alignment projection
│   ├── BARWLS.sh                ← Chemical barcode orchestrator
│   └── chemical_barcode_analyzer_v2.py
├── results/                     ← Per-run output directories
├── logs/                        ← Per-run logs
└── manifest/                    ← JSON manifests + checksums
```

---

## Quick start

```bash
# Full run: three gate thresholds + barcode
./run_pipeline.sh \
  --data-root ~/Proc-decode_Weblogo_Data_Folder \
  --samples ~/Proc-decode_Weblogo_Data_Folder/Taxon_MDS_data/samples.tsv \
  --gates "50/85,60/90,70/90" \
  --label C12_25taxon \
  --barcode

# Single gate, no barcode
./run_pipeline.sh \
  --data-root ~/Proc-decode_Weblogo_Data_Folder \
  --samples ~/Proc-decode_Weblogo_Data_Folder/Taxon_MDS_data/samples.tsv \
  --gates "60/90" \
  --label C12_60-90_only

# Dry run (validate without executing)
./run_pipeline.sh \
  --data-root ~/Proc-decode_Weblogo_Data_Folder \
  --samples ~/Proc-decode_Weblogo_Data_Folder/Taxon_MDS_data/samples.tsv \
  --dry-run
```

---

## samples.tsv format

Tab-separated, one row per taxon:

```
Branchiostoma   /home/snerx/Proc-decode_Weblogo_Data_Folder/Taxon_MDS_data/Branchiostoma   1   309   310   381
```

Columns: `SAMPLE  WORKDIR  C12_START  C12_END  MYH_START  MYH_END`

Each WORKDIR must contain `md.gro`, `md.tpr`, `md.xtc`.
Lines starting with `#` are skipped.

---

## Environment variables

| Variable | Default | Effect |
|----------|---------|--------|
| `MAKE_PROTEIN_ONLY` | `1` | Generate protein_only.gro via gmx |
| `MAKE_PROT_XTC` | `0` | Also generate md_protein.xtc |
| `OCC_THR` | `0.05` | Minimum occupancy for watchlist inclusion |
| `P1_THR` | (none) | Optional partner 1 occupancy threshold |
| `P2_THR` | (none) | Optional partner 2 occupancy threshold |

---

## What each step does

### Step 1: PDWLS.sh
For each taxon in samples.tsv:
1. Optionally generates protein_only.gro via gmx trjconv
2. Runs prco_decode_cli.py (per-residue contact occupancy from MD trajectory)
3. Generates SC watchlists (residues exceeding occupancy threshold)

Then across all taxa:
4. Concatenates C12 FASTA sequences, runs MAFFT alignment
5. For each gate threshold: projects watchlists onto alignment (project_multiple_taxa.py)
6. Generates WebLogo (PNG, SVG, PDF) from each projected FASTA

### Step 2: BARWLS.sh (optional, --barcode flag)
Scans projection outputs and for each projected FASTA:
1. Classifies residues into chemistry classes (H, P, B, A)
2. Computes per-column occupancy, dominance, entropy
3. Identifies anchor columns
4. Produces summary tables and per-column profiles

### Step 3: Checksums + manifest
- SHA256 checksums of all output files
- JSON manifest with parameters, tool versions, provenance

---

## Output structure

Each run creates:
```
results/<label>_<timestamp>/
├── samples.tsv           ← Copy of input
├── prco_decode/          ← Per-taxon PRCO CSVs
├── alignment/            ← MAFFT alignment (all_taxa.fa, C12_aligned.fa)
├── projection/           ← Projected FASTAs + column maps (per gate)
├── logo/                 ← WebLogo PNGs, SVGs, PDFs (per gate)
├── logs/                 ← MAFFT log
└── barcode/              ← BARWLS output (if --barcode)
```

---

## Dependencies

- Python 3.8+ with: numpy, pandas, MDAnalysis
- MAFFT
- WebLogo 3.7+
- GROMACS 2025.2 (only if protein_only.gro not pre-generated)

---

## Design principle

The scripts in `scripts/` are unmodified working copies of the
analysis tools developed and validated during the IBAM project.
The orchestrator wraps them without altering their internal logic —
it handles directory setup, parameter passing, logging, and
provenance tracking.
