# ibam-grammar-engine

**Structure-based interaction grammar analysis of the C12orf29 (IBAM) protein across 25 taxa**

Reproducibility pipeline for the manuscript: *C12orf29 encodes IBAM (Interdigitating Bridge Actin-Myosin), a sarcomeric protein with a conserved actomyosin binding grammar spanning ~1 billion years of evolution.*

---

## Overview

Conventional multiple sequence alignments of C12orf29 orthologs yield 30–40% pairwise identity — sufficient to confirm homology, but too low to make confident claims about conservation of specific functional residues. This pipeline overcomes that limitation by replacing linear sequence comparison with a structure-based projection: instead of asking whether the same residues occupy the same positions in a primary sequence alignment, it asks whether the same residues make the same three-dimensional contacts at the C12–myosin tail (MyhT) binding interface.

The pipeline integrates three analytical stages:

1. **Dynamic contact decoding** from molecular dynamics trajectories (GROMACS)
2. **Evolutionary projection** of interaction residues across taxa into a common positional framework
3. **Chemical barcode analysis** of the resulting interaction grammar

The output is a projected FASTA alignment in which each column represents not a linear sequence position but a structurally recurrent contact site on the binding interface. Conservation signals extracted from this alignment reflect constraint on the three-dimensional interaction surface rather than on primary structure.

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
│  │   └── 4.5 Å cutoff, ≥30% occupancy      │
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
- **Python** ≥ 3.10 with the following packages:
  - `numpy`
  - `pandas`
  - `matplotlib`
  - `mdtraj`
  - `biopython`
- **MAFFT** (sequence alignment)
- **WebLogo** 3.7+ (sequence logo generation)
- **AF3 predicted structures** for each taxon (see Data Availability)

### Installation

```bash
git clone https://github.com/thorfriisphd-rgb/ibam-grammar-engine.git
cd ibam-grammar-engine
```

No additional installation is required. All pipeline scripts are self-contained within the `scripts/` directory.

### Running the Pipeline

The Grammar Engine is invoked via the master orchestrator:

```bash
bash run_pipeline.sh \
  --label C12_25taxon \
  --data_root /path/to/your/data \
  --gates "50/85,60/90,70/90" \
  --barcode
```

#### Key flags

| Flag | Description |
|------|-------------|
| `--label` | Run identifier (used in output directory names and manifests) |
| `--data_root` | Path to the data directory containing taxon MD trajectories and FASTA files |
| `--gates` | Comma-separated dual-gate thresholds (core/chem) |
| `--barcode` | Enable chemical barcode analysis (BARWLS step) |
| `--dry-run` | Validate configuration without executing |

### Output Structure

Each run produces a timestamped results directory:

```
results/C12_25taxon_YYYYMMDD_HHMMSS/
├── alignment/          # MAFFT alignments
├── prco_decode/        # Per-taxon contact tables
├── projection/         # Projected FASTA files (one per gating regime)
├── logo/               # WebLogo PNGs (one per gating regime)
├── barcode/            # Chemical barcode analysis
│   └── barcode_*/
│       ├── BAR_*_chemcomp.tsv     # Long-format chemistry compositions
│       ├── BAR_*_percol.tsv       # Per-column summary statistics
│       ├── BAR_*_summary.tsv      # Run-level summaries
│       ├── BAR_*_stacked.png      # Stacked barcode figures
│       ├── mini_results_table.tsv # Consolidated results
│       └── MANIFEST.txt           # Barcode run manifest
└── logs/               # Per-taxon processing logs
```

Checksums and manifests are written to:

```
manifest/
├── checksums_*.sha256   # SHA256 hashes of all input/output files
└── manifest_*.json      # Run parameters, tool versions, timestamps
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

## GROMACS Parameters

All molecular dynamics simulations use the following standardised parameters:

| Parameter | Value | Reference |
|-----------|-------|-----------|
| Force field | CHARMM36-jul2022 | |
| Water model | TIP3P | |
| Electrostatics | Cutoff (1.2 nm) | Beck et al. 2005; Piana et al. 2012 |
| Production MD | 10 ns (standard); 50/100 ns in reserve | |
| Temperature | 310 K | |
| Contact cutoff | 4.5 Å | |
| Occupancy threshold | ≥30% | |

---

## Key Quantitative Outputs

The pipeline produces the following core metrics for each gating regime:

- **Projected FASTA width** — number of columns surviving dual-gate filtering
- **Mean anchor dominance** — mean fractional dominance of the top chemistry class per column
- **n99 count** — number of positions with >99% single-class dominance
- **Chemistry barcode** — per-column fractional occupancy of H (hydrophobic), P (polar), B (basic), A (acidic) classes

Threshold-robustness is demonstrated by running multiple gating regimes. The core interaction grammar (K/R towers, D punctuation, W anchor, H-rich tail) is stable across all tested thresholds.

---

## Reproducibility

This pipeline was designed for full end-to-end reproducibility:

- **Checksummed outputs** — SHA256 hashes for every input and output file
- **Version tracking** — GROMACS, Python, MAFFT, WebLogo versions recorded per run
- **Manifest logging** — JSON manifest captures all parameters, paths, and timestamps
- **Self-contained** — all scripts are within the repository; no external dependencies beyond standard bioinformatics tools
- **Deterministic** — fixed random seeds where applicable (e.g., permutation tests: seed=42)

To reproduce the complete analysis:

1. Clone this repository
2. Download AF3 structures and MD trajectories from Zenodo (DOI: *[to be assigned]*)
3. Set `--data_root` to the downloaded data directory
4. Run: `bash run_pipeline.sh --label C12_25taxon --gates "50/85,60/90,70/90" --barcode`
5. Compare output checksums against the reference manifest

---

## Data Availability

- **Code**: This repository (MIT license)
- **Data**: AF3 predicted structures and GROMACS MD trajectories for all 25 taxa — Zenodo DOI *[to be assigned]*
- **Control structure**: GCN4 leucine zipper (PDB: 4DMD) used for coiled-coil validation

---

## Citation

If you use this pipeline, please cite:

> Friis, T. (2026). C12orf29 encodes IBAM (In Between Actin and Myosin), a sarcomeric protein with a conserved actomyosin binding grammar. *bioRxiv* [preprint]. DOI: *[to be assigned]*

---

## License

MIT License. See [LICENSE](LICENSE) for details.

---

## Contact

Thor Friis — [@thorfriisphd-rgb](https://github.com/thorfriisphd-rgb)

Independent researcher ("Ronin Researcher"), Bodø, Norway.
PhD in Molecular Biology, Queensland University of Technology (QUT).
