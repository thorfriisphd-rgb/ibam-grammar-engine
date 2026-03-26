#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# run_pipeline.sh — IBAM Grammar Engine
#
# Single entry point wrapping PDWLS.sh (PRCO → projection →
# WebLogo) and BARWLS.sh (chemical barcode analysis).
#
# Required:
#   --data-root PATH    Data directory root
#   --samples FILE      samples.tsv (taxon, workdir, residue ranges)
#
# Options:
#   --gates STRING      Gate thresholds, comma-separated
#                       (default: "50/85,60/90,70/90")
#   --label STRING      Run label (default: IBAM_grammar)
#   --barcode           Also run BARWLS chemical barcode analysis
#   --barcode-label S   BARWLS label (default: barcode_<label>)
#   --dry-run           Print commands without executing
#   --help              Show this message
#
# Example:
#   ./run_pipeline.sh \
#     --data-root /path/to/data \
#     --samples /path/to/data/Taxon_MDS_data/samples.tsv \
#     --gates "50/85,60/90,70/90" \
#     --label C12_25taxon \
#     --barcode
# ============================================================

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="$PIPELINE_DIR/scripts"

# ── Defaults ───────────────────────────────────────────────
DATA_ROOT=""
SAMPLES_FILE=""
GATES="50/85,60/90,70/90"
LABEL="IBAM_grammar"
RUN_BARCODE=false
BARCODE_LABEL=""
DRY_RUN=false
RUN_ID="$(date +%Y%m%d_%H%M%S)"

# ── Argument parsing ───────────────────────────────────────
while [[ $# -gt 0 ]]; do
  case "$1" in
    --data-root)      DATA_ROOT="$2";      shift 2 ;;
    --samples)        SAMPLES_FILE="$2";   shift 2 ;;
    --gates)          GATES="$2";          shift 2 ;;
    --label)          LABEL="$2";          shift 2 ;;
    --barcode)        RUN_BARCODE=true;    shift   ;;
    --barcode-label)  BARCODE_LABEL="$2";  shift 2 ;;
    --dry-run)        DRY_RUN=true;        shift   ;;
    --help)
      sed -n '2,30p' "$0" | grep '^#' | sed 's/^# \{0,1\}//'
      exit 0 ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

# ── Validate required args ─────────────────────────────────
[[ -n "$DATA_ROOT" ]]   || { echo "[FATAL] --data-root required"; exit 1; }
[[ -n "$SAMPLES_FILE" ]] || { echo "[FATAL] --samples required"; exit 1; }
[[ -d "$DATA_ROOT" ]]   || { echo "[FATAL] data-root not found: $DATA_ROOT"; exit 1; }
[[ -f "$SAMPLES_FILE" ]] || { echo "[FATAL] samples not found: $SAMPLES_FILE"; exit 1; }

# Resolve to absolute paths
SAMPLES_FILE="$(cd "$(dirname "$SAMPLES_FILE")" && pwd)/$(basename "$SAMPLES_FILE")"
DATA_ROOT="$(cd "$DATA_ROOT" && pwd)"
[[ -n "$BARCODE_LABEL" ]] || BARCODE_LABEL="barcode_${LABEL}"

# ── Directories ────────────────────────────────────────────
RUN_DIR="$PIPELINE_DIR/results/${LABEL}_${RUN_ID}"
LOG_DIR="$PIPELINE_DIR/logs"
MANIFEST_DIR="$PIPELINE_DIR/manifest"
mkdir -p "$RUN_DIR" "$LOG_DIR" "$MANIFEST_DIR"
LOG_FILE="$LOG_DIR/${LABEL}_${RUN_ID}.log"

# ── Logging ────────────────────────────────────────────────
log() {
  local level="$1"; shift
  echo "[$(date '+%H:%M:%S')] [$level] $*" | tee -a "$LOG_FILE"
}

# ── Banner ─────────────────────────────────────────────────
log "INFO" "================================================"
log "INFO" "IBAM Grammar Engine"
log "INFO" "Run ID:     $RUN_ID"
log "INFO" "Label:      $LABEL"
log "INFO" "Data root:  $DATA_ROOT"
log "INFO" "Samples:    $SAMPLES_FILE"
log "INFO" "Gates:      $GATES"
log "INFO" "Barcode:    $RUN_BARCODE"
log "INFO" "Work dir:   $RUN_DIR"
[[ "$DRY_RUN" == "true" ]] && log "INFO" "*** DRY RUN ***"
log "INFO" "================================================"

# ── Pre-flight: tools ──────────────────────────────────────
ok=true
for tool in python3 mafft weblogo; do
  if command -v "$tool" &>/dev/null; then
    log "OK" "$tool → $(command -v "$tool")"
  else
    log "FATAL" "Required: $tool"; ok=false
  fi
done

for f in PDWLS.sh prco_decode_cli.py project_multiple_taxa.py; do
  [[ -f "$SCRIPTS/$f" ]] || { log "FATAL" "Missing: $SCRIPTS/$f"; ok=false; }
done

if [[ "$RUN_BARCODE" == "true" ]]; then
  for f in BARWLS.sh chemical_barcode_analyzer_v2.py; do
    [[ -f "$SCRIPTS/$f" ]] || { log "FATAL" "Missing: $SCRIPTS/$f"; ok=false; }
  done
fi

[[ "$ok" == "true" ]] || exit 1

# gmx check (warning only — needed if protein_only.gro is missing)
if command -v gmx &>/dev/null; then
  log "OK" "gmx → $(command -v gmx)"
else
  log "WARN" "gmx not found (needed only if protein_only.gro missing)"
fi

# ── Capture versions ───────────────────────────────────────
VERSIONS="$MANIFEST_DIR/versions_${LABEL}_${RUN_ID}.txt"
{
  echo "=== IBAM Grammar Engine — Versions ==="
  echo "Run: $RUN_ID  Label: $LABEL"
  echo "Date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "Python:  $(python3 --version 2>&1)"
  echo "MAFFT:   $(mafft --version 2>&1 | head -1 || echo unknown)"
  echo "WebLogo: $(weblogo --version 2>&1 | head -1 || echo unknown)"
  echo "GROMACS: $(gmx --version 2>&1 | grep 'GROMACS version' | head -1 || echo 'not found')"
  echo "OS:      $(uname -srm)"
  echo "Host:    $(hostname)"
  echo "Gates:   $GATES"
  echo "Taxa:    $(grep -cvE '^#|^$' "$SAMPLES_FILE" || echo 0)"
} > "$VERSIONS"
log "INFO" "Versions: $VERSIONS"

# ── Taxa count ─────────────────────────────────────────────
N_TAXA=$(grep -cvE '^#|^$' "$SAMPLES_FILE" || echo 0)
log "INFO" "Taxa: $N_TAXA"

# ── STEP 1: PDWLS ──────────────────────────────────────────
log "INFO" "──────────────────────────────────────────"
log "INFO" "STEP 1: PDWLS (PRCO → projection → WebLogo)"
log "INFO" "──────────────────────────────────────────"

# PDWLS reads samples.tsv from CWD, outputs to CWD subdirectories
cp "$SAMPLES_FILE" "$RUN_DIR/samples.tsv"

if [[ "$DRY_RUN" == "false" ]]; then
  (
    cd "$RUN_DIR"
    export DATA_ROOT="$DATA_ROOT"
    export MAKE_PROTEIN_ONLY="${MAKE_PROTEIN_ONLY:-1}"
    export MAKE_PROT_XTC="${MAKE_PROT_XTC:-0}"
    export OCC_THR="${OCC_THR:-0.05}"
    bash "$SCRIPTS/PDWLS.sh" --gates "$GATES" --data-root "$DATA_ROOT"
  ) 2>&1 | tee -a "$LOG_FILE"

  PDWLS_EXIT=${PIPESTATUS[0]:-$?}
  if [[ $PDWLS_EXIT -ne 0 ]]; then
    log "FATAL" "PDWLS.sh exited $PDWLS_EXIT"
    exit $PDWLS_EXIT
  fi
  log "INFO" "PDWLS complete."
else
  log "DRY" "cd $RUN_DIR && bash $SCRIPTS/PDWLS.sh --gates $GATES --data-root $DATA_ROOT"
fi

# ── STEP 2: BARWLS (optional) ──────────────────────────────
if [[ "$RUN_BARCODE" == "true" ]]; then
  log "INFO" "──────────────────────────────────────────"
  log "INFO" "STEP 2: BARWLS (chemical barcode analysis)"
  log "INFO" "──────────────────────────────────────────"

  BARWLS_OUT="$RUN_DIR/barcode"
  mkdir -p "$BARWLS_OUT"

  if [[ "$DRY_RUN" == "false" ]]; then
    bash "$SCRIPTS/BARWLS.sh" \
      --c12_root "$RUN_DIR" \
      --myht_root "$RUN_DIR" \
      --out_root "$BARWLS_OUT" \
      --label "$BARCODE_LABEL" \
      2>&1 | tee -a "$LOG_FILE"
    log "INFO" "BARWLS complete."
  else
    log "DRY" "bash $SCRIPTS/BARWLS.sh --c12_root $RUN_DIR --out_root $BARWLS_OUT --label $BARCODE_LABEL"
  fi
fi

# ── STEP 3: Checksums + manifest ───────────────────────────
log "INFO" "──────────────────────────────────────────"
log "INFO" "STEP 3: Checksums and manifest"
log "INFO" "──────────────────────────────────────────"

CHECKSUMS="$MANIFEST_DIR/checksums_${LABEL}_${RUN_ID}.sha256"
if [[ "$DRY_RUN" == "false" ]]; then
  find "$RUN_DIR" -type f | sort | xargs sha256sum > "$CHECKSUMS" 2>/dev/null || true
  log "INFO" "Checksums: $CHECKSUMS ($(wc -l < "$CHECKSUMS") files)"
fi

# Convert bash boolean to Python boolean for JSON
if [[ "$RUN_BARCODE" == "true" ]]; then
  PY_BARCODE="True"
else
  PY_BARCODE="False"
fi

MANIFEST="$MANIFEST_DIR/manifest_${LABEL}_${RUN_ID}.json"
if [[ "$DRY_RUN" == "false" ]]; then
  python3 - <<PYEOF > "$MANIFEST"
import json, os, platform
from datetime import datetime, timezone

print(json.dumps({
    "run_id": "${RUN_ID}",
    "label": "${LABEL}",
    "pipeline": "IBAM_Grammar_Engine",
    "version": "1.0.0",
    "created_at": datetime.now(timezone.utc).isoformat(),
    "system": {
        "hostname": platform.node(),
        "os": platform.platform(),
        "python": platform.python_version(),
    },
    "parameters": {
        "data_root": "${DATA_ROOT}",
        "samples_file": "${SAMPLES_FILE}",
        "gates": "${GATES}",
        "n_taxa": ${N_TAXA},
        "occ_thr": float(os.environ.get("OCC_THR", "0.05")),
    },
    "outputs": {
        "work_dir": "${RUN_DIR}",
        "barcode": ${PY_BARCODE},
    },
    "description": (
        "Reproducibility manifest for C12orf29/IBAM grammar engine. "
        "PRCO contact decoding, evolutionary projection, dual-gate "
        "filtering, WebLogo generation, and chemical barcode analysis."
    ),
}, indent=2))
PYEOF
  log "INFO" "Manifest: $MANIFEST"
fi

# ── Done ───────────────────────────────────────────────────
log "INFO" "================================================"
log "INFO" "Grammar Engine complete."
log "INFO" "  Results:   $RUN_DIR"
log "INFO" "  Log:       $LOG_FILE"
log "INFO" "  Manifest:  $MANIFEST"
log "INFO" "  Checksums: $CHECKSUMS"
log "INFO" "================================================"

# Quick listing of key outputs
if [[ "$DRY_RUN" == "false" && -d "$RUN_DIR" ]]; then
  echo ""
  echo "── Output directories ──"
  find "$RUN_DIR" -maxdepth 2 -type d | sort | sed "s|$RUN_DIR|.|"
  echo ""
  PROJ_COUNT=$(find "$RUN_DIR" -name "MG_projected_trimmed_*.fa" 2>/dev/null | wc -l)
  LOGO_COUNT=$(find "$RUN_DIR" -name "MG_core_*.png" 2>/dev/null | wc -l)
  echo "Projected FASTAs: $PROJ_COUNT"
  find "$RUN_DIR" -name "MG_projected_trimmed_*.fa" 2>/dev/null | sort | sed "s|$RUN_DIR/||"
  echo ""
  echo "WebLogos: $LOGO_COUNT"
  find "$RUN_DIR" -name "MG_core_*.png" 2>/dev/null | sort | sed "s|$RUN_DIR/||"
fi
