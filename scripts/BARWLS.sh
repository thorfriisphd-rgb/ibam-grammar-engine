#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# BARWLS.sh
# One-command barcode summary builder for C12 + MyhT projected FASTAs.
# One-command barcode summary builder for C12 + MyhT projected FASTAs.
# ------------------------------------------------------------

# -------- defaults (match your current layout) ---------------
HOME_DIR="${HOME}"
TOOLS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

DEFAULT_C12_ROOT=""
DEFAULT_MYHT_ROOT=""
DEFAULT_OUT_ROOT=""

PYTHON_BIN="python3"
TOPN="10"
LABEL="Barcode_run"
RUN_IDENTITY_CHECK="0"
C12_GATE=""
MYHT_GATE=""
NULL_N="200"
SEED="1"

# -------- helper functions -----------------------------------
die() { echo "[FATAL] $*" >&2; exit 1; }
info() { echo "[INFO] $*"; }
ok() { echo "[OK] $*"; }

usage() {
  cat <<EOF
Usage:
  bash BARWLS.sh [options]

Options:
  --c12_root PATH         Root containing C12 Results (default: ${DEFAULT_C12_ROOT})
  --myht_root PATH        Root containing MyhT Results (default: ${DEFAULT_MYHT_ROOT})
  --out_root PATH         Output root (default: ${DEFAULT_OUT_ROOT})
  --label STRING          Label for this run folder (default: ${LABEL})
  --python PATH           Python executable (default: ${PYTHON_BIN})
  --top N                 Top anchors to report (default: ${TOPN})

  --identity_check        Also run identity-map complementarity sanity check (optional)
  --c12_gate STRING       Gate substring to pick C12 FASTA for identity check (e.g. core60_chem90)
  --myht_gate STRING      Gate substring to pick MyhT FASTA for identity check (e.g. core35_chem70)
  --null N                Null shuffles for identity check (default: ${NULL_N})
  --seed N                RNG seed for identity check (default: ${SEED})

Examples:
  # Standard run: barcode everything found under Results_C12 and Results_MyhT
  bash BARWLS.sh --label Barcode_beta-test_04.03.26

  # With optional identity-map complementarity sanity check
  bash BARWLS.sh --label Barcode_beta-test_04.03.26 \\
    --identity_check --c12_gate core60_chem90 --myht_gate core35_chem70

EOF
}

# -------- parse args -----------------------------------------
C12_ROOT="${DEFAULT_C12_ROOT}"
MYHT_ROOT="${DEFAULT_MYHT_ROOT}"
OUT_ROOT="${DEFAULT_OUT_ROOT}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --c12_root) C12_ROOT="$2"; shift 2 ;;
    --myht_root) MYHT_ROOT="$2"; shift 2 ;;
    --out_root) OUT_ROOT="$2"; shift 2 ;;
    --label) LABEL="$2"; shift 2 ;;
    --python) PYTHON_BIN="$2"; shift 2 ;;
    --top) TOPN="$2"; shift 2 ;;
    --identity_check) RUN_IDENTITY_CHECK="1"; shift 1 ;;
    --c12_gate) C12_GATE="$2"; shift 2 ;;
    --myht_gate) MYHT_GATE="$2"; shift 2 ;;
    --null) NULL_N="$2"; shift 2 ;;
    --seed) SEED="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown option: $1 (use --help)" ;;
  esac
done

# -------- validate environment -------------------------------
[[ -d "${TOOLS_DIR}" ]] || die "Tools dir not found: ${TOOLS_DIR}"
[[ -d "${C12_ROOT}" ]] || die "C12 root not found: ${C12_ROOT}"
[[ -d "${MYHT_ROOT}" ]] || die "MyhT root not found: ${MYHT_ROOT}"
command -v "${PYTHON_BIN}" >/dev/null 2>&1 || die "Python not found: ${PYTHON_BIN}"

BARCODE_PY="${TOOLS_DIR}/chemical_barcode_analyzer_v2.py"
COMP_PY="${TOOLS_DIR}/complementary_pattern_analyzer_v2.py"
PLOT_PY="${TOOLS_DIR}/plot_barcode_stacked.py"

[[ -f "${BARCODE_PY}" ]] || die "Missing tool: ${BARCODE_PY}"
[[ -f "${COMP_PY}" ]] || die "Missing tool: ${COMP_PY}"
[[ -f "${PLOT_PY}" ]] || die "Missing tool: ${PLOT_PY}"

# numeric checks
[[ "${TOPN}" =~ ^[0-9]+$ ]] || die "--top must be integer"
[[ "${NULL_N}" =~ ^[0-9]+$ ]] || die "--null must be integer"
[[ "${SEED}" =~ ^[0-9]+$ ]] || die "--seed must be integer"

# -------- create run output folder ---------------------------
TS="$(date +%Y%m%d_%H%M%S)"
RUN_DIR="${OUT_ROOT}/${LABEL}_${TS}"
LOG_DIR="${RUN_DIR}/logs"
mkdir -p "${RUN_DIR}" "${LOG_DIR}"

info "Run dir: ${RUN_DIR}"
info "C12 root: ${C12_ROOT}"
info "MyhT root: ${MYHT_ROOT}"
info "Python: ${PYTHON_BIN}"
info "Top anchors: ${TOPN}"
echo

# -------- write manifest for reproducibility -----------------
MANIFEST="${RUN_DIR}/MANIFEST.txt"
{
  echo "BARWLS run manifest"
  echo "timestamp: ${TS}"
  echo "label: ${LABEL}"
  echo "tools_dir: ${TOOLS_DIR}"
  echo "c12_root: ${C12_ROOT}"
  echo "myht_root: ${MYHT_ROOT}"
  echo "python: ${PYTHON_BIN}"
  echo "topN: ${TOPN}"
  echo
  echo "tool_files_md5:"

  (cd "${TOOLS_DIR}" && md5sum \
  "$(basename "${BARCODE_PY}")" \
  "$(basename "${COMP_PY}")" \
  "$(basename "${PLOT_PY}")" \
  "run_robustness_analysis_v2.sh" 2>/dev/null || true)

  echo "identity_check_enabled: ${RUN_IDENTITY_CHECK}"
  echo "identity_check_c12_gate: ${C12_GATE}"
  echo "identity_check_myht_gate: ${MYHT_GATE}"
  echo "identity_check_null: ${NULL_N}"
  echo "identity_check_seed: ${SEED}"
} > "${MANIFEST}"

ok "Wrote manifest: ${MANIFEST}"
echo

# -------- find projected FASTAs ------------------------------
# We only care about the projection outputs
# pattern: */projection/MG_projected_trimmed_*.fa
C12_LIST="${RUN_DIR}/c12_fastas.txt"
MYHT_LIST="${RUN_DIR}/myht_fastas.txt"

find "${C12_ROOT}" -type f -path "*/projection/MG_projected_trimmed_*.fa" | sort > "${C12_LIST}"
find "${MYHT_ROOT}" -type f -path "*/projection/MG_projected_trimmed_*.fa" | sort > "${MYHT_LIST}"

C12_N=$(wc -l < "${C12_LIST}" | tr -d ' ')
MYHT_N=$(wc -l < "${MYHT_LIST}" | tr -d ' ')

info "Found C12 projected FASTAs: ${C12_N}"
info "Found MyhT projected FASTAs: ${MYHT_N}"

[[ "${C12_N}" -gt 0 ]] || die "No C12 projected FASTAs found under ${C12_ROOT}"
[[ "${MYHT_N}" -gt 0 ]] || die "No MyhT projected FASTAs found under ${MYHT_ROOT}"
echo

# -------- run barcode analyzer over everything ---------------
run_barcode_on_list() {
  local list_file="$1"
  local tag_prefix="$2"  # "C12" or "MYHT"

  while IFS= read -r fa; do
    [[ -f "${fa}" ]] || continue
    base="$(basename "${fa}" .fa)"
    out_summary="${RUN_DIR}/BAR_${base}_summary.tsv"
    out_percol="${RUN_DIR}/BAR_${base}_percol.tsv"

    info "Barcoding ${tag_prefix}: ${base}"
    "${PYTHON_BIN}" "${BARCODE_PY}" \
      --fasta "${fa}" \
      --out "${out_summary}" \
      --per_column_out "${out_percol}" \
      --top "${TOPN}" \
      > "${LOG_DIR}/barcode_${base}.log" 2>&1

    ok "Wrote: $(basename "${out_summary}") and $(basename "${out_percol}")"
  done < "${list_file}"
}

run_barcode_on_list "${C12_LIST}" "C12"
run_barcode_on_list "${MYHT_LIST}" "MYHT"

echo
ok "Barcode analyses complete."
echo

# -------- build mini results table (gate, width, occupancy, mean_anchor_dom, n99) ----
MINI="${RUN_DIR}/mini_results_table.tsv"
MINI_SORT="${RUN_DIR}/mini_results_table_sorted.tsv"

echo -e "gate\twidth\tmean_occupancy\tmean_anchor_dominance\tn99" > "${MINI}"

# We aggregate from the BAR_*_summary.tsv files in RUN_DIR
for f in "${RUN_DIR}"/BAR_*_summary.tsv; do
  [[ -f "$f" ]] || continue

  # These summary TSVs are whitespace-separated in practice (even if named .tsv),
  # and row 2 contains: run_id gate n_taxa width mean_occupancy barcode anchors barcode_similarity_to_ref
  gate=$(awk 'NR==2{print $2}' "$f")
  width=$(awk 'NR==2{print $4}' "$f")
  occ=$(awk 'NR==2{print $5}' "$f")
  anchors=$(awk 'NR==2{print $7}' "$f")

  mean_dom=$(echo "$anchors" | tr ',' '\n' | awk -F: '{s+=$2; c++} END{if(c) printf("%.3f", s/c); else print "NA"}')
  n99=$(echo "$anchors" | tr ',' '\n' | awk -F: '$2+0>=0.99{c++} END{print c+0}')

  echo -e "${gate}\t${width}\t${occ}\t${mean_dom}\t${n99}" >> "${MINI}"
done

{ head -n 1 "${MINI}"; tail -n +2 "${MINI}" | sort -k1,1; } > "${MINI_SORT}"

ok "Wrote: ${MINI}"
ok "Wrote: ${MINI_SORT}"
echo

# -------- split convenience views (optional) -----------------
MYHT_VIEW="${RUN_DIR}/mini_results_MyhT.tsv"
C12_VIEW="${RUN_DIR}/mini_results_C12.tsv"

{
  head -n 1 "${MINI_SORT}"
  grep -E '^core3|^core4' "${MINI_SORT}" || true
} > "${MYHT_VIEW}"

{
  head -n 1 "${MINI_SORT}"
  grep -E '^core5|^core6|^core7' "${MINI_SORT}" || true
} > "${C12_VIEW}"

ok "Wrote: ${MYHT_VIEW}"
ok "Wrote: ${C12_VIEW}"
echo

# -------- render stacked barcode(s) --------------------------
info "Rendering stacked barcode(s)"

for percol in "${RUN_DIR}"/BAR_*_chemcomp.tsv; do
  [[ -f "${percol}" ]] || continue
  base="$(basename "${percol}" _chemcomp.tsv)"
  outpng="${RUN_DIR}/${base}_stacked.png"

  "${PYTHON_BIN}" "${PLOT_PY}" \
    --infile "${percol}" \
    --out "${outpng}" \
    --title "${base}" \
    >> "${LOG_DIR}/plot_barcode_stacked.log" 2>&1

  ok "Rendered: $(basename "${outpng}")"
done

ok "Stacked barcode render complete."
echo

# -------- optional: identity-map complementarity sanity check -
if [[ "${RUN_IDENTITY_CHECK}" == "1" ]]; then
  [[ -n "${C12_GATE}" ]] || die "--identity_check requires --c12_gate"
  [[ -n "${MYHT_GATE}" ]] || die "--identity_check requires --myht_gate"

  # pick first matching FASTA from lists
  C12_FA="$(grep -F "${C12_GATE}" "${C12_LIST}" | head -n 1 || true)"
  MYHT_FA="$(grep -F "${MYHT_GATE}" "${MYHT_LIST}" | head -n 1 || true)"

  [[ -f "${C12_FA}" ]] || die "Could not find C12 FASTA matching gate '${C12_GATE}' in ${C12_LIST}"
  [[ -f "${MYHT_FA}" ]] || die "Could not find MyhT FASTA matching gate '${MYHT_GATE}' in ${MYHT_LIST}"

  info "Identity check using:"
  info "  C12:  ${C12_FA}"
  info "  MyhT: ${MYHT_FA}"

  # compute lengths from first record (counts sequence letters across wrapped lines)
  C12_LEN=$(awk 'BEGIN{L=0} /^>/{if(L>0){print L; exit}} !/^>/{L+=length($0)} END{if(L>0)print L}' "${C12_FA}" | head -n1)
  MYH_LEN=$(awk 'BEGIN{L=0} /^>/{if(L>0){print L; exit}} !/^>/{L+=length($0)} END{if(L>0)print L}' "${MYHT_FA}" | head -n1)
  [[ "${C12_LEN}" =~ ^[0-9]+$ ]] || die "Could not parse C12_LEN"
  [[ "${MYH_LEN}" =~ ^[0-9]+$ ]] || die "Could not parse MYH_LEN"

  MINLEN=$(( C12_LEN<MYH_LEN ? C12_LEN : MYH_LEN ))
  [[ "${MINLEN}" -gt 0 ]] || die "MINLEN computed as 0"

  IDMAP="${RUN_DIR}/identity_map_1to1.tsv"
  seq 1 "${MINLEN}" | awk '{print $1"\t"$1}' > "${IDMAP}"

  COMP_OUT="${RUN_DIR}/complementarity_identitymap_${C12_GATE}_VS_${MYHT_GATE}.tsv"

  "${PYTHON_BIN}" "${COMP_PY}" \
    --c12_fasta "${C12_FA}" \
    --myht_fasta "${MYHT_FA}" \
    --contact_map "${IDMAP}" \
    --null "${NULL_N}" \
    --seed "${SEED}" \
    --out "${COMP_OUT}" \
    > "${LOG_DIR}/complementarity_identitymap.log" 2>&1

  ok "Wrote identity-map complementarity: ${COMP_OUT}"
  echo
fi

ok "BARWLS complete."
info "Run folder ready: ${RUN_DIR}"
