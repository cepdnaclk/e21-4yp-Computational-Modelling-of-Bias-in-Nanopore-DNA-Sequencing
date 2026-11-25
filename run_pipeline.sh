#!/usr/bin/env bash
# Run the Phase 1 pipeline end to end.
#
#   ./run_pipeline.sh                          # config/experiment.yaml, full design
#   ./run_pipeline.sh config/test.yaml         # small test
#   ./run_pipeline.sh config/experiment.yaml --mode ofat
#   ./run_pipeline.sh config/experiment.yaml --select basecall=guppy-6.5.7:hac --select minq=8,9,10
#   (any option after the config goes to bin/design.py: --mode, --select, --strict)
#   Preview without running: python bin/design.py config/experiment.yaml --out ~/tmp/d.csv --select ... --list
#
# Steps: 1) design.py builds design.csv   2) Nextflow runs every combination
# (-resume reuses finished stages)   3) collect_results.py builds the summary tables.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="$(realpath "${1:-$REPO/config/experiment.yaml}")"
shift || true
MODE_ARGS=("$@")

CONDA_BASE="${FYP_CONDA_BASE:-$( (conda info --base) 2>/dev/null || echo "$HOME/anaconda3")}"
export FYP_CONDA_BASE="$CONDA_BASE"   # read by nextflow.config
ENVS="$CONDA_BASE/envs"
export PATH="$CONDA_BASE/condabin:$HOME/bin:$PATH"   # Nextflow activates the envs with conda
PY="$ENVS/env-analysis/bin/python"
NXF="$ENVS/env-workflow/bin/nextflow"
export TMPDIR="$HOME/tmp" NXF_TEMP="$HOME/tmp"
mkdir -p "$TMPDIR"

OUTDIR="$("$PY" -c "import yaml,sys; print(yaml.safe_load(open(sys.argv[1]))['outdir'])" "$CONFIG")"
RUN_ID="$(date -u +%Y%m%dT%H%M%SZ)"
RUN_DIR="$OUTDIR/_runs/$RUN_ID"
WORK_DIR="$(dirname "$OUTDIR")/nf-work"   # same filesystem as outputs, so results are hard-linked, not copied
mkdir -p "$RUN_DIR"

# 1) Design
"$PY" "$REPO/bin/design.py" "$CONFIG" --out "$RUN_DIR/design.csv" "${MODE_ARGS[@]}"
cp "$CONFIG" "$RUN_DIR/experiment.yaml"
cp "$RUN_DIR/design.csv" "$OUTDIR/_runs/latest_design.csv"
{
  echo "run_id: $RUN_ID"
  echo "git_commit: $(git -C "$REPO" rev-parse HEAD 2>/dev/null || echo none)"
  echo "git_dirty: $(git -C "$REPO" status --porcelain 2>/dev/null | wc -l) changed files"
  echo "nextflow: $("$NXF" -version 2>&1 | grep -o 'version [0-9.]*')"
  echo "started_utc: $(date -u +%FT%TZ)"
} > "$RUN_DIR/run_info.yaml"

# 2) Execute (launch dir = run dir's parent keeps .nextflow cache for -resume)
cd "$OUTDIR/_runs"
"$NXF" run "$REPO/main.nf" \
  -params-file "$CONFIG" \
  --design "$RUN_DIR/design.csv" \
  --conda_base "$CONDA_BASE" \
  -w "$WORK_DIR" \
  -resume \
  -with-trace "$RUN_DIR/trace.tsv" \
  -with-report "$RUN_DIR/report.html" \
  -with-timeline "$RUN_DIR/timeline.html" \
  ${FYP_NXF_ARGS:-}        # e.g. FYP_NXF_ARGS="-process.errorStrategy=ignore" to keep going past failures

# 3) Collect every metadata.json + results.tsv into summary tables
"$PY" "$REPO/bin/collect_results.py" --outdir "$OUTDIR"
echo "finished_utc: $(date -u +%FT%TZ)" >> "$RUN_DIR/run_info.yaml"
echo "Done. Summary: $OUTDIR/_summary/"
