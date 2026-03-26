#!/usr/bin/env bash
set -euo pipefail

############################################
# PDWLS.sh (v3.2 clean)
# samples.tsv -> (optional gmx protein_only) -> PRCO -> SC watchlists
# -> MAFFT -> projection -> WebLogo
############################################

DATA_ROOT="${DATA_ROOT:-}"
TOOLS_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SAMPLES_FILE="samples.tsv"
GATES="40/80"

# Watchlist thresholds
OCC_THR="${OCC_THR:-0.05}"
P1_THR="${P1_THR:-}"   # optional
P2_THR="${P2_THR:-}"   # optional

# Base FASTA directories
FASTA_ROOT="${DATA_ROOT}/FastaFolder"
C12_DIR="${FASTA_ROOT}/C12_Fastafiles"

# Optional gmx preprocessing
MAKE_PROTEIN_ONLY="${MAKE_PROTEIN_ONLY:-1}"  # 1=yes, 0=no
MAKE_PROT_XTC="${MAKE_PROT_XTC:-0}"          # 1=yes, 0=no

while [[ $# -gt 0 ]]; do
  case "$1" in
    --gates)     GATES="$2"; shift 2 ;;
    --data-root) DATA_ROOT="$2"; shift 2 ;;
    *) echo "Unknown argument: $1"; exit 1 ;;
  esac
done

mkdir -p prco_decode alignment projection logo logs

[[ -f "$SAMPLES_FILE" ]] || { echo "[FATAL] samples.tsv not found in WD: $(pwd)"; exit 1; }

echo "=========================================="
echo "PDWLS v3 Run"
echo "WD: $(pwd)"
echo "DATA_ROOT: $DATA_ROOT"
echo "GATES: $GATES"
echo "Watchlist thr: OCC_THR=$OCC_THR P1_THR=${P1_THR:-NA} P2_THR=${P2_THR:-NA}"
echo "GMX: MAKE_PROTEIN_ONLY=$MAKE_PROTEIN_ONLY MAKE_PROT_XTC=$MAKE_PROT_XTC"
echo "=========================================="

# ---------- QC + taxa list ----------
TAXA=()
while IFS=$'\t ' read -r SAMPLE WORKDIR C12S C12E MYHS MYHE; do
  [[ "${SAMPLE:-}" =~ ^# ]] && continue
  [[ -z "${SAMPLE// }" ]] && continue

  [[ -n "${WORKDIR:-}" ]] || { echo "[FATAL] WORKDIR empty for $SAMPLE"; exit 1; }

  for f in md.gro md.tpr md.xtc; do
    [[ -f "$WORKDIR/$f" ]] || { echo "[FATAL] Missing $f in $WORKDIR for $SAMPLE"; exit 1; }
  done

  if ! [[ "$C12S" =~ ^[0-9]+$ && "$C12E" =~ ^[0-9]+$ && "$MYHS" =~ ^[0-9]+$ && "$MYHE" =~ ^[0-9]+$ ]]; then
    echo "[FATAL] Non-integer residue range for $SAMPLE"
    exit 1
  fi
  if (( C12S > C12E || MYHS > MYHE )); then
    echo "[FATAL] Inverted residue range for $SAMPLE"
    exit 1
  fi

  TAXA+=("$SAMPLE")
done < "$SAMPLES_FILE"

echo "Taxa detected (${#TAXA[@]}):"
for t in "${TAXA[@]}"; do echo "  - $t"; done
echo "------------------------------------------"

# ---------- PRCO + SC watchlists ----------
WATCHDIR="${DATA_ROOT}/WatchLists"
mkdir -p "$WATCHDIR"

while IFS=$'\t ' read -r SAMPLE WORKDIR C12S C12E MYHS MYHE; do
  [[ "${SAMPLE:-}" =~ ^# ]] && continue
  [[ -z "${SAMPLE// }" ]] && continue

  # Optional gmx: generate protein_only.gro / md_protein.xtc inside taxon WORKDIR
  if [[ "$MAKE_PROTEIN_ONLY" == "1" ]]; then
    if [[ ! -f "$WORKDIR/protein_only.gro" ]]; then
      echo "[INFO] gmx: generating protein_only.gro for $SAMPLE"
      ( cd "$WORKDIR"
        export GMX_MAXBACKUP=-1
        printf "Protein\n" | gmx trjconv -s md.tpr -f md.gro -o protein_only.gro >/dev/null
      )
    fi
    if [[ "$MAKE_PROT_XTC" == "1" && ! -f "$WORKDIR/md_protein.xtc" ]]; then
      echo "[INFO] gmx: generating md_protein.xtc for $SAMPLE"
      ( cd "$WORKDIR"
        export GMX_MAXBACKUP=-1
        printf "Protein\n" | gmx trjconv -s md.tpr -f md.xtc -o md_protein.xtc >/dev/null
      )
    fi
  fi

  echo "[INFO] PRCO decode: $SAMPLE"

# Choose matching topology/trajectory pair
  TOPO="$WORKDIR/md.tpr"
  TRAJ="$WORKDIR/md.xtc"

  OUT_PREFIX="prco_decode/${SAMPLE}"
  python "$TOOLS_ROOT/prco_decode_cli.py" \
    --top  "$TOPO" \
    --traj "$TRAJ" \
    --c12  "protein and resid ${C12S}:${C12E}" \
    --myh  "protein and resid ${MYHS}:${MYHE}" \
    --out  "$OUT_PREFIX"

  PRCO_CSV="${OUT_PREFIX}_prco.csv"
  [[ -f "$PRCO_CSV" ]] || { echo "[FATAL] Expected PRCO CSV not found: $PRCO_CSV"; exit 1; }

  WATCH_OUT="${WATCHDIR}/${SAMPLE}_watch.txt"
  META_OUT="${WATCHDIR}/${SAMPLE}_watch.meta.tsv"

  python - <<PY
import pandas as pd
import sys
from pathlib import Path

prco = Path(r"$PRCO_CSV")
out  = Path(r"$WATCH_OUT")
thr  = float(r"$OCC_THR")

p1_raw = r"$P1_THR".strip()
p2_raw = r"$P2_THR".strip()
p1_thr = float(p1_raw) if p1_raw else None
p2_thr = float(p2_raw) if p2_raw else None

df = pd.read_csv(prco)

required = ["C12_resid","occupancy","partner1_occ","partner2_occ"]
missing = [c for c in required if c not in df.columns]
if missing:
    sys.stderr.write(f"[FATAL] PRCO CSV missing columns {missing} in {prco}\n")
    sys.stderr.write(f"Found columns: {list(df.columns)}\n")
    sys.exit(3)

keep = df["occupancy"] >= thr
if p1_thr is not None:
    keep = keep & (df["partner1_occ"] >= p1_thr)
if p2_thr is not None:
    keep = keep & (df["partner2_occ"] >= p2_thr)

res = (df.loc[keep, "C12_resid"]
         .dropna()
         .astype(int)
         .drop_duplicates()
         .sort_values()
         .tolist())

out.write_text("\n".join(map(str, res)) + ("\n" if res else ""), encoding="utf-8")

msg = f"[INFO] SC watchlist: {out} (n={len(res)}) occ>={thr}"
if p1_thr is not None: msg += f" p1>={p1_thr}"
if p2_thr is not None: msg += f" p2>={p2_thr}"
print(msg)
PY

  [[ -s "$WATCH_OUT" ]] || { echo "[FATAL] SC watchlist empty for $SAMPLE (OCC_THR=$OCC_THR)"; exit 1; }

  {
    echo -e "timestamp\t$(date -Iseconds)"
    echo -e "source_prco_csv\t$PRCO_CSV"
    echo -e "occ_thr\t$OCC_THR"
    echo -e "p1_thr\t${P1_THR:-}"
    echo -e "p2_thr\t${P2_THR:-}"
    echo -e "top_used\t$TOPO"
    echo -e "traj_used\t$TRAJ"
  } > "$META_OUT"

done < "$SAMPLES_FILE"

echo "[INFO] PRCO + SC watchlists complete."
echo "------------------------------------------"

# ---------- MAFFT alignment ----------
FASTA_LIST="alignment/all_taxa.fa"
> "$FASTA_LIST"

while IFS=$'\t ' read -r SAMPLE WORKDIR C12S C12E MYHS MYHE; do
  [[ "${SAMPLE:-}" =~ ^# ]] && continue
  [[ -z "${SAMPLE// }" ]] && continue

C12_FASTA_DIR="${DATA_ROOT}/FastaFolder/C12_Fastafiles"

FASTA="${C12_FASTA_DIR}/${SAMPLE}_C12.fa"
[[ -f "$FASTA" ]] || { echo "[FATAL] Missing FASTA: $FASTA"; exit 1; }
cat "$FASTA" >> "$FASTA_LIST"
done < "$SAMPLES_FILE"

echo "[INFO] Running MAFFT..."
mafft --auto "$FASTA_LIST" > alignment/C12_aligned.fa 2> logs/mafft.log

# ---------- Projection + WebLogo ----------
IFS=',' read -ra GATE_ARRAY <<< "$GATES"

for G in "${GATE_ARRAY[@]}"; do
  CORE="${G%/*}"
  CHEM="${G#*/}"

  CORE_FRAC="$(python - <<PY
print(float($CORE)/100.0)
PY
)"
  CHEM_FRAC="$(python - <<PY
print(float($CHEM)/100.0)
PY
)"

  echo "[INFO] Projection: core${CORE} chem${CHEM}"

  PROJ_OUT="projection/MG_projected_trimmed_n${#TAXA[@]}_core${CORE}_chem${CHEM}.fa"
  COLMAP_OUT="projection/MG_column_map_n${#TAXA[@]}_core${CORE}_chem${CHEM}.tsv"

  python "$TOOLS_ROOT/project_multiple_taxa.py" \
    --aln alignment/C12_aligned.fa \
    --watchdir "$WATCHDIR" \
    --min-taxa-frac "$CORE_FRAC" \
    --chem-thr "$CHEM_FRAC" \
    --out-fa "$PROJ_OUT" \
    --out-map "$COLMAP_OUT"

  echo "[INFO] WebLogo: core${CORE} chem${CHEM}"
  PYTHONWARNINGS="ignore:pkg_resources is deprecated as an API:UserWarning" \
    weblogo -f "$PROJ_OUT" -o "logo/MG_core_n${#TAXA[@]}_core${CORE}_chem${CHEM}.png" --format png
  PYTHONWARNINGS="ignore:pkg_resources is deprecated as an API:UserWarning" \
    weblogo -f "$PROJ_OUT" -o "logo/MG_core_n${#TAXA[@]}_core${CORE}_chem${CHEM}.svg" --format svg
    PYTHONWARNINGS="ignore:pkg_resources is deprecated as an API:UserWarning" \
  weblogo -f "$PROJ_OUT" -o "logo/MG_core_n${#TAXA[@]}_core${CORE}_chem${CHEM}.pdf" --format pdf
done

echo "=========================================="
echo "PDWLS run complete."
echo "=========================================="
