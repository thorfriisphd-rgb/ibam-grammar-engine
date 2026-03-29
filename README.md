# IBAM Grammar Engine

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19241765.svg)](https://doi.org/10.5281/zenodo.19241765)

![GitHub last commit](https://img.shields.io/github/last-commit/thorfriisphd-rgb/ibam-grammar-engine)
![GitHub repo size](https://img.shields.io/github/repo-size/thorfriisphd-rgb/ibam-grammar-engine)



**Structure-based interaction grammar analysis of C12orf29 (IBAM) across 25 taxa**

Reproducibility pipeline for the manuscript: *C12orf29 encodes IBAM (In Between Actin and Myosin), a sarcomeric protein with a conserved actomyosin binding grammar spanning ~1 billion years of evolution.*

---

## Overview

Conventional multiple sequence alignments of C12orf29 orthologs yield 30–40% pairwise identity — sufficient to confirm homology, but too low to make confident claims about conservation of specific functional residues. This pipeline overcomes that limitation by replacing linear sequence comparison with a structure-based projection: instead of asking whether the same residues occupy the same positions in a primary sequence alignment, it asks whether the same residues make the same three-dimensional contacts at the C12–myosin tail (MyhT) binding interface.

The pipeline integrates three analytical stages:

1. **Dynamic contact decoding** from molecular dynamics trajectories (GROMACS)
2. **Evolutionary projection** of interaction residues across taxa into a common positional framework
3. **Chemical barcode analysis** of the resulting interaction grammar

The output is a projected FASTA alignment in which each column represents not a linear sequence position but a structurally recurrent contact site on the binding interface. Conservation signals extracted from this alignment reflect constraint on the three-dimensional interaction surface rather than on primary structure.

---
## Conceptual framework

The IBAM Grammar Engine identifies conserved interaction chemistry at the IBAM–myosin tail interface using a structure-guided evolutionary projection. Instead of relying on conventional multiple sequence alignment alone, the pipeline first decodes persistent IBAM–MyhT contacts from molecular-dynamics trajectories and then projects these interface residues into a shared coordinate system representing positions along the IBAM binding groove. Orthologous sequences from 25 taxa are subsequently mapped onto these interface coordinates, allowing conservation to be evaluated at the level of interaction chemistry rather than primary sequence identity. This approach reveals a conserved pattern of basic, hydrophobic, and aromatic residues that together form an evolutionarily preserved actomyosin interaction grammar spanning ~1 billion years of divergence.


---

## Taxa

The pipeline processes 25 taxa spanning ~1 billion years of eukaryotic divergence, plus a bacterial outgroup:

| Clade | Representative taxa |
|-------|-------------------|
| Bacteria (outgroup) | *Hahella chejuensis* |
| Heterolobosea | *Naegleria gruberi*, *Willaertia magna* |
| Euglenozoa | *Euglena longa*, *Eutreptiella gymnastica* |
| Cnidaria | *Clytia hemispherica*, *Podocoryna carnea* |
| Annelida | *Eisenia fetida*, *Lamellibrachia satsuma* |
| Mollusca | *Haliotis asinina*, *Magallana angulata*, *Octopus vulgaris* |
| Brachiopoda | *Lingula anatina* |
| Phoronida | *Phoronis australis* |
| Myriapoda | *Henia illyrica*, *Lithobius forficatus* |
| Onychophora | *Euperipatoides rowelli* |
| Chordata — Tunicata | *Ciona intestinalis*, *Salpa thompsoni* |
| Chordata — Cephalochordata | *Branchiostoma floridae* |
| Chordata — Vertebrata | *Lampetra planeri*, *Myxine glutinosa*, *Mus musculus*, *Ovis aries* |
| Hemichordata | *Saccoglossus kowalevskii* |

---

## Pipeline Architecture

```
AF3 predicted structures (per taxon × MyhT pair)
        │
        ▼
┌─────────────────────────────────────────────┐
│  STEP 1: PDWLS                              │
│  ├── PRCO decode (contact extraction)       │
│  │   └── 4.5 Å cutoff, ≥30% occupancy       │
│  ├── MAFFT alignment                        │
│  ├── Evolutionary projection                │
│  │   └── Watchlist mapping → projected FASTA│
│  ├── Dual-gate filtering                    │
│  │   └── Occupancy gate × Chemistry gate    │
│  └── WebLogo generation                     │
├─────────────────────────────────────────────┤
│  STEP 2: BARWLS                             │
│  ├── Chemical barcode analysis              │
│  │   └── H / P / B / A class frequencies    │
│  ├── Mini results tables                    │
│  └── Stacked barcode rendering              │
├─────────────────────────────────────────────┤
│  STEP 3: Checksums and manifest             │
│  ├── SHA256 checksums (all inputs/outputs)  │
│  ├── Tool versions                          │
│  └── JSON manifest                          │
└─────────────────────────────────────────────┘
```

---

## Dual-Gate Filtering

The raw projected alignment contains columns with varying levels of phylogenetic support and chemical coherence. To extract the conserved interaction core, a dual-gate filtering scheme is applied:

- **Occupancy gate (core threshold)** — requires a residue to be present in at least *x*% of 25 taxa. Ensures phylogenetic breadth.
- **Chemical dominance gate (chem threshold)** — requires a single physicochemical class (H, P, B, or A) to account for at least *y*% of residues in the column. Ensures functional coherence.

Both gates must be passed simultaneously. The pipeline evaluates multiple gating regimes by default (50/85, 60/90, 70/90) to demonstrate threshold robustness. The core invariant positions are stable across all tested thresholds; progressively stricter gating contracts the cassette by trimming peripheral positions while preserving the functional core.

---

## Quick Start

### Prerequisites

- **GROMACS** ≥ 2025.x (with `gmx` in PATH)
- **Python** ≥ 3.8 with: `numpy`, `pandas`, `MDAnalysis`
- **MAFFT** (sequence alignment)
- **WebLogo** 3.7+ (sequence logo generation)
- pdf2svg (required by WebLogo for SVG output; sudo apt install pdf2svg on Debian/Ubuntu)
- **matplotlib** (barcode rendering)
- **AF3 predicted structures** for each taxon (see Data Availability)

### Installation

```bash
git clone https://github.com/thorfriisphd-rgb/ibam-grammar-engine.git
cd ibam-grammar-engine
```

No additional installation required. All pipeline scripts are self-contained within the `scripts/` directory.

### Download Trajectory Data

Download the MD trajectory dataset from Zenodo:

https://doi.org/10.5281/zenodo.19241765

Extract and point `--data-root` at the extracted directory.

### Running the Pipeline

```bash
# Full run: three gate thresholds + barcode
./run_pipeline.sh \
  --data-root /path/to/data \
  --samples /path/to/data/samples.tsv \
  --gates "50/85,60/90,70/90" \
  --label C12_25taxon \
  --barcode

# Single gate, no barcode
./run_pipeline.sh \
  --data-root /path/to/data \
  --samples /path/to/data/samples.tsv \
  --gates "60/90" \
  --label C12_60-90_only

# Dry run (validate without executing)
./run_pipeline.sh \
  --data-root /path/to/data \
  --samples /path/to/data/samples.tsv \
  --dry-run
```

### Key Flags

| Flag | Description |
|------|-------------|
| `--data-root` | Path to the data directory containing taxon MD trajectories and FASTA files |
| `--samples` | Path to samples.tsv taxon registry |
| `--gates` | Comma-separated dual-gate thresholds (core/chem) |
| `--label` | Run identifier (used in output directory names and manifests) |
| `--barcode` | Enable chemical barcode analysis (BARWLS step) |
| `--dry-run` | Validate configuration without executing |

---

## samples.tsv Format

Tab-separated, one row per taxon:

```
Branchiostoma   /path/to/data/Taxon_MDS_data/Branchiostoma   1   309   310   381
```

Columns: `SAMPLE  WORKDIR  C12_START  C12_END  MYH_START  MYH_END`

Each WORKDIR must contain `md.gro`, `md.tpr`, `md.xtc`.
Lines starting with `#` are skipped.

---

## Output Structure

Each run creates a timestamped results directory:

```
results/<label>_<timestamp>/
├── prco_decode/           ← Per-taxon contact tables
├── alignment/             ← MAFFT alignments
├── projection/            ← Projected FASTA files (one per gating regime)
├── logo/                  ← WebLogo PNGs (one per gating regime)
├── barcode/               ← Chemical barcode analysis
│   └── barcode_*/
│       ├── BAR_*_chemcomp.tsv     ← Long-format chemistry compositions
│       ├── BAR_*_percol.tsv       ← Per-column summary statistics
│       ├── BAR_*_summary.tsv      ← Run-level summaries
│       ├── BAR_*_stacked.png      ← Stacked barcode figures
│       ├── mini_results_table.tsv ← Consolidated results
│       └── MANIFEST.txt           ← Barcode run manifest
└── logs/                  ← Per-taxon processing logs
```

Checksums and manifests:

```
manifest/
├── checksums_*.sha256     ← SHA256 hashes of all input/output files
├── manifest_*.json        ← Run parameters, tool versions, timestamps
└── versions_*.txt         ← Tool version strings
```

---

## Scripts

| Script | Purpose |
|--------|---------|
| `run_pipeline.sh` | Master orchestrator |
| `scripts/PDWLS.sh` | PRCO decode → projection → WebLogo pipeline |
| `scripts/BARWLS.sh` | Chemical barcode analysis and rendering |
| `scripts/prco_decode_cli.py` | Time-integrated contact decoding from MD trajectories |
| `scripts/project_multiple_taxa.py` | Watchlist mapping and projected FASTA generation |
| `scripts/chemical_barcode_analyzer_v2.py` | Physicochemical class frequency analysis |
| `scripts/complementary_pattern_analyzer_v2.py` | Complementary pattern analysis (optional sanity check) |
| `scripts/plot_barcode_stacked.py` | Stacked chemistry barcode figure rendering |

---

## What Each Step Does

### Step 1: PDWLS

For each taxon in samples.tsv:
1. Optionally generates protein_only.gro via `gmx trjconv`
2. Runs `prco_decode_cli.py` — extracts per-residue contact occupancy from the MD trajectory (4.5 Å cutoff, ≥30% occupancy)
3. Generates SC watchlists (structurally recurrent contact residues)

Then across all taxa:
4. Concatenates C12 FASTA sequences, runs MAFFT alignment
5. For each gate threshold: projects watchlists onto alignment (`project_multiple_taxa.py`)
6. Generates WebLogo from each projected FASTA

### Step 2: BARWLS (optional, --barcode flag)

For each projected FASTA:
1. Classifies residues into chemistry classes — H (hydrophobic), P (polar), B (basic), A (acidic)
2. Computes per-column occupancy, dominance, and entropy
3. Identifies anchor columns
4. Produces summary tables and per-column profiles
5. Renders stacked chemistry barcode figures

### Step 3: Checksums and Manifest

- SHA256 checksums of all output files
- JSON manifest with parameters, tool versions, and provenance
- Tool version strings recorded for reproducibility

---

## GROMACS Parameters

All molecular dynamics simulations use standardised parameters:

| Parameter | Value | Reference |
|-----------|-------|-----------|
| Force field | CHARMM36-jul2022 | |
| Water model | TIP3P | |
| Electrostatics | Cutoff (1.0 nm) | Beck et al. 2005; Piana et al. 2012 |
| Production MD | 10 ns (standard); 50/100 ns in reserve | |
| Temperature | 310 K | |
| Contact cutoff | 4.5 Å | |
| Occupancy threshold | ≥30% | |

---

## Environment Variables

| Variable | Default | Effect |
|----------|---------|--------|
| `MAKE_PROTEIN_ONLY` | `1` | Generate protein_only.gro via gmx |
| `MAKE_PROT_XTC` | `0` | Also generate md_protein.xtc |
| `OCC_THR` | `0.05` | Minimum occupancy for watchlist inclusion |
| `P1_THR` | (none) | Optional partner 1 occupancy threshold |
| `P2_THR` | (none) | Optional partner 2 occupancy threshold |

---

---
## Troubleshooting
**WebLogo fails with `OSError: ... requires the program 'pdf2svg'`**
Install the system package: `sudo apt install pdf2svg` (Debian/Ubuntu) or `brew install pdf2svg` (macOS). The Python `weblogo` library depends on this external tool for SVG output but does not install it automatically.
 
**`samples.tsv` contains wrong paths / `WORKDIR not found` errors**
The shipped `samples.tsv` contains absolute paths from the development machine. Edit the `WORKDIR` column so each path points to the corresponding taxon directory inside your local copy of the Zenodo data. The file lives at the root of the Zenodo deposit, not inside `Taxon_MDS_data/`.
 
---

## Reproducibility

This pipeline was designed for full end-to-end reproducibility:

- **Checksummed outputs** — SHA256 hashes for every input and output file
- **Version tracking** — GROMACS, Python, MAFFT, WebLogo versions recorded per run
- **Manifest logging** — JSON manifest captures all parameters, paths, and timestamps
- **Self-contained** — all scripts within the repository; no external dependencies beyond standard bioinformatics tools
- **Deterministic** — fixed random seeds where applicable

To reproduce the complete analysis:

1. Clone this repository
2. Download MD trajectories from Zenodo: https://doi.org/10.5281/zenodo.19241765
3. Extract and set `--data-root` to the extracted directory
4. Run: `./run_pipeline.sh --data-root /path/to/data --samples /path/to/data/Taxon_MDS_data/samples.tsv --gates "50/85,60/90,70/90" --label C12_25taxon --barcode`
5. Compare output checksums against the reference manifest

---

## Data Availability

- **Code**: This repository (MIT license)
- **Data**: AF3 predicted structures and GROMACS MD trajectories for all 25 taxa — https://doi.org/10.5281/zenodo.19241765
- **Control structure**: GCN4 leucine zipper (PDB: 4DMD) used for coiled-coil validation

---

## Citation

If you use this pipeline, please cite:

> Friis, T. (2026). C12orf29 encodes IBAM (In Between Actin and Myosin), a sarcomeric protein with a conserved actomyosin binding grammar. *bioRxiv* [preprint]. DOI: [to be assigned]

---

## Design Principle

The scripts in `scripts/` are unmodified working copies of the analysis tools developed and validated during the IBAM project. The orchestrator wraps them without altering their internal logic — it handles directory setup, parameter passing, logging, and provenance tracking.

---

## License

MIT License. See [LICENSE](LICENSE) for details.

---

## Contact

Thor Friis — [@thorfriisphd-rgb](https://github.com/thorfriisphd-rgb)

Independent researcher, Bodø, Norway.
PhD in Molecular Biology, Queensland University of Technology (QUT).
