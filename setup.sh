#!/usr/bin/env bash
# Rebuild the project toolchain for the current Linux user. No sudo needed.
#
#   ./setup.sh                      everything: conda, environments, Dorado, Guppy, models, ~/.bashrc
#   ./setup.sh --skip-basecallers   environments only (no Dorado / Guppy / models)
#   ./setup.sh --recreate env-qc    delete and rebuild one environment (repeatable)
#   ./setup.sh --no-bashrc          do not touch ~/.bashrc
#   ./setup.sh --check              only print what is installed
#   ./setup.sh --export             after changing an environment: rewrite envs/*.yml and envs/lock/*
#   ./setup.sh --env-suffix -test   build every environment under a new name (env-qc-test, ...)
#                                   to test a rebuild without touching the working ones
#
# Safe to re-run: anything that already exists is kept. Everything goes under
# $HOME (conda, ~/tools, ~/bin, ~/tmp). Run it on a login node with internet.
#
# What it installs
#   1. conda: uses an existing Anaconda/Miniconda/Miniforge, else installs Miniforge to ~/miniforge3
#   2. channels: conda-forge + bioconda, strict priority (~/.condarc)
#   3. every environment, in this order of exactness:
#        a. envs/lock/<env>.linux-64.lock  exact conda builds (same files, md5-checked, no solving)
#           + envs/lock/<env>.pip.txt      exact pip packages (Python modules from PyPI)
#           + SimReadUntil at a fixed git commit
#        b. envs/<env>.yml                 same versions, builds re-solved (other Linux / conda)
#        c. unpinned package list          last resort; prints a warning (compare versions!)
#      then imports every Python module the project's scripts use, in the env that runs them
#   4. Dorado 0.9.6 + R9.4.1 models (~/tools/dorado_models), Guppy 6.5.7 and 6.2.1 (the lab's
#      version) GPU builds, Bonito R9.4.1 models, links in ~/bin
#   5. NCBI taxonomy (taxdump) in the data folder, used by every classifier summary
#   6. ~/.bashrc: conda init, ~/bin on PATH, FYP_CONDA_BASE, DORADO_MODELS, TMPDIR
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DORADO_VERSION=0.9.6          # last release line with R9.4.1 models and FAST5 input
GUPPY_VERSIONS=(6.5.7 6.2.1)     # 6.2.1 = the version the lab used (Method PDF)
BONITO_MODELS=(dna_r9.4.1_e8_fast@v3.4 dna_r9.4.1_e8_hac@v3.3 dna_r9.4.1_e8_sup@v3.3)
TAXONOMY_DIR="${FYP_TAXONOMY:-$HOME/dna_data/FYP/databases/taxonomy}"
DORADO_MODELS=(dna_r9.4.1_e8_fast@v3.4 dna_r9.4.1_e8_hac@v3.3 dna_r9.4.1_e8_sup@v3.3 dna_r9.4.1_e8_hac@v3.3_5mCG_5hmCG@v0)
TOOLS="$HOME/tools"
BIN="$HOME/bin"
export TMPDIR="$HOME/tmp"
LOG="$HOME/.cache/fyp-setup"

# Unpinned fallback for each environment (the packages originally installed).
declare -A FALLBACK=(
  [env-basecall]="python=3.10 samtools hdf5 h5py pip"
  [env-qc]="nanofilt nanoplot nanostat porechop porechop_abi seqkit chopper filtlong cutadapt fastqc python-edlib"
  [env-align]="minimap2 samtools bedtools winnowmap ngmlr lra python=3.10"
  [env-classify]="kraken2 bracken centrifuge krakentools"
  [env-analysis]="python=3.10 pandas numpy scipy matplotlib seaborn biopython blast pysam pyarrow h5py scikit-learn statsmodels lightgbm shap scikit-bio jupyterlab ncbi-genome-download ncbi-datasets-cli seqkit minimap2 pyyaml last mmseqs2 diamond"
  [env-hostile]="hostile"
  [env-taxa]="krakenuniq centrifuger kaiju metabuli ganon taxor kmcp sylph taxonkit"
  [env-protein]="diamond megan"
  [env-bonito]="python=3.10 pip"
  [env-signal]="ont-modkit f5c slow5tools"
  [env-emu]="emu flatten-dict"
  [env-metamaps]="metamaps"
  [env-assembly]="flye"
  [env-sim]="python=3.10 nanosim pip"
  [env-workflow]="nextflow"
)
# pip packages that are not on PyPI under these names; installed separately below.
PIP_SKIP_REGEX='^ *- (simreaduntil)=='
SIMREADUNTIL="git+https://github.com/ratschlab/sim_read_until.git@785beee0b2e3c99d74a89e69241777747b28f79a"
LOCKDIR="$REPO/envs/lock"
PLATFORM="linux-64"

# Python modules each environment must provide (third-party only; the rest is stdlib).
declare -A PY_MODULES=(
  [env-basecall]="h5py numpy pod5 ont_fast5_api"
  [env-qc]="edlib"
  [env-align]=""
  [env-classify]=""
  [env-analysis]="yaml pandas numpy scipy matplotlib seaborn Bio pysam pyarrow h5py sklearn statsmodels lightgbm shap skbio jupyterlab"
  [env-sim]="simreaduntil"
  [env-bonito]="bonito torch"
)

SKIP_BASECALLERS=0; DO_BASHRC=1; CHECK_ONLY=0; EXPORT=0; SUFFIX=""; RECREATE=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --skip-basecallers) SKIP_BASECALLERS=1 ;;
    --no-bashrc) DO_BASHRC=0 ;;
    --check) CHECK_ONLY=1 ;;
    --recreate) RECREATE+=("$2"); shift ;;
    --export) EXPORT=1 ;;
    --env-suffix) SUFFIX="$2"; shift ;;
    -h|--help) sed -n '2,32p' "$0"; exit 0 ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
  shift
done

step() { printf '\n==> %s\n' "$*"; }
warn() { printf '  ! %s\n' "$*" >&2; }
mkdir -p "$TMPDIR" "$BIN" "$TOOLS" "$LOG"

# ---------------------------------------------------------------- 0. platform
step "Checking platform"
[[ "$(uname -s)" == Linux ]] || { echo "Linux only"; exit 1; }
ARCH="$(uname -m)"
[[ "$ARCH" == x86_64 ]] || warn "Architecture $ARCH: envs/*.yml were made on x86_64; expect the unpinned fallback and an arm64 Dorado."
if command -v nvidia-smi >/dev/null; then
  nvidia-smi --query-gpu=name,driver_version --format=csv,noheader | sed 's/^/  GPU: /'
else
  warn "No NVIDIA GPU found: Dorado/Guppy will be installed but basecalling needs a GPU node."
fi
df -h "$HOME" | awk 'NR==2 {print "  Free space in home: " $4}'

# ---------------------------------------------------------------- 1. conda
step "Finding conda"
CONDA_BASE="${FYP_CONDA_BASE:-}"
if [[ -z "$CONDA_BASE" ]]; then
  for c in "$(command -v conda 2>/dev/null || true)" "$HOME/anaconda3/bin/conda" "$HOME/miniconda3/bin/conda" "$HOME/miniforge3/bin/conda"; do
    if [[ -n "$c" && -x "$c" ]]; then CONDA_BASE="$("$c" info --base)"; break; fi
  done
fi
if [[ -z "$CONDA_BASE" ]]; then
  [[ $CHECK_ONLY == 1 ]] && { echo "  conda not installed"; exit 1; }
  echo "  No conda found: installing Miniforge to ~/miniforge3"
  curl -fsSL -o "$TMPDIR/miniforge.sh" "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-${ARCH}.sh"
  bash "$TMPDIR/miniforge.sh" -b -p "$HOME/miniforge3"
  rm -f "$TMPDIR/miniforge.sh"
  CONDA_BASE="$HOME/miniforge3"
fi
CONDA="$CONDA_BASE/bin/conda"
ENVS="$CONDA_BASE/envs"
echo "  conda $("$CONDA" --version | awk '{print $2}') at $CONDA_BASE"

# ---------------------------------------------------------------- export (maintainers)
if [[ $EXPORT == 1 ]]; then
  step "Exporting environments to envs/ and envs/lock/"
  mkdir -p "$LOCKDIR"
  for d in "$ENVS"/env-*; do
    name="$(basename "$d")"; [[ -d "$d/conda-meta" ]] || continue
    "$CONDA" env export -n "$name" --no-builds 2>/dev/null | grep -v '^prefix:' > "$REPO/envs/$name.yml"
    "$CONDA" list -n "$name" --explicit --md5 > "$LOCKDIR/$name.$PLATFORM.lock"
    # pip-installed modules only (conda marks them with channel "pypi"); git installs pinned by commit
    "$CONDA" list -n "$name" | awk '$4 == "pypi" && $1 != "simreaduntil" {print $1 "==" $2}' > "$LOCKDIR/$name.pip.txt"
    if "$CONDA" list -n "$name" | awk '$1 == "simreaduntil"' | grep -q .; then echo "$SIMREADUNTIL" >> "$LOCKDIR/$name.pip.txt"; fi
    echo "  $name: $(grep -c '^https' "$LOCKDIR/$name.$PLATFORM.lock") conda packages, $(wc -l < "$LOCKDIR/$name.pip.txt") pip packages"
  done
  echo "Commit envs/ to Git."
  exit 0
fi

if [[ $CHECK_ONLY == 0 ]]; then
  step "Configuring channels (conda-forge, bioconda, strict)"
  "$CONDA" config --remove channels defaults >/dev/null 2>&1 || true
  "$CONDA" config --add channels bioconda
  "$CONDA" config --add channels conda-forge      # added last = highest priority
  "$CONDA" config --set channel_priority strict
  "$CONDA" config --set solver libmamba 2>/dev/null || true
  "$CONDA" config --show channels | sed 's/^/  /'
fi

# ---------------------------------------------------------------- 2. environments
env_exists() { [[ -d "$ENVS/$1/conda-meta" ]]; }

install_pip() {   # exact pip modules from the lock, then git packages
  local target="$1" src="$2" req="$LOCKDIR/$2.pip.txt"
  [[ -s "$req" ]] || return 0
  grep -v '^git+' "$req" > "$TMPDIR/$target.pip.txt" || true
  if [[ -s "$TMPDIR/$target.pip.txt" ]]; then
    "$ENVS/$target/bin/python" -m pip install -q --no-deps -r "$TMPDIR/$target.pip.txt" >> "$LOG/$target.log" 2>&1 \
      || { warn "$target: pip modules failed (see $LOG/$target.log)"; return 1; }
  fi
  rm -f "$TMPDIR/$target.pip.txt"
  if grep -q '^git+' "$req"; then   # SimReadUntil: its submodule URL uses SSH, rewritten to HTTPS here
    GIT_CONFIG_COUNT=1 GIT_CONFIG_KEY_0=url.https://github.com/.insteadOf GIT_CONFIG_VALUE_0=git@github.com: \
      "$ENVS/$target/bin/python" -m pip install -q --no-deps $(grep '^git+' "$req") >> "$LOG/$target.log" 2>&1 \
      || { warn "$target: git package failed (see $LOG/$target.log)"; return 1; }
  fi
}

create_env() {
  local name="$1" target="$1$SUFFIX" lock="$LOCKDIR/$1.$PLATFORM.lock" yml="$REPO/envs/$1.yml" tmp="$TMPDIR/$1.setup.yml"
  : > "$LOG/$target.log"
  # a) exact builds from the lock file (Linux x86_64 only)
  if [[ "$ARCH" == x86_64 && -f "$lock" ]]; then
    echo "  creating $target from envs/lock (exact builds)..."
    if "$CONDA" create -y -n "$target" --file "$lock" >> "$LOG/$target.log" 2>&1 && install_pip "$target" "$name"; then
      echo "  OK $target (exact)"; return 0
    fi
    warn "$target: exact lock failed (see $LOG/$target.log); trying envs/$name.yml"
    "$CONDA" env remove -y -n "$target" >/dev/null 2>&1 || true
  fi
  # b) same versions, builds re-solved
  grep -v '^prefix:' "$yml" | grep -Ev "$PIP_SKIP_REGEX" > "$tmp"
  if "$CONDA" env create -y -n "$target" -f "$tmp" >> "$LOG/$target.log" 2>&1; then
    rm -f "$tmp"
    grep -q '^git+' "$LOCKDIR/$name.pip.txt" 2>/dev/null && install_pip_git "$target" "$name"
    echo "  OK $target (same versions, builds re-solved)"; return 0
  fi
  rm -f "$tmp"
  warn "$target: pinned versions could not be installed here; trying unpinned packages"
  "$CONDA" env remove -y -n "$target" >/dev/null 2>&1 || true
  # c) unpinned
  # shellcheck disable=SC2086
  "$CONDA" create -y -n "$target" ${FALLBACK[$name]} >> "$LOG/$target.log" 2>&1 \
    || { warn "$target FAILED (see $LOG/$target.log)"; return 1; }
  case "$name" in
    env-basecall) "$ENVS/$target/bin/python" -m pip install -q pod5 ont-fast5-api "setuptools<81" >> "$LOG/$target.log" 2>&1 ;;
    env-bonito)   "$ENVS/$target/bin/python" -m pip install -q ont-bonito==0.7.3 "setuptools<81" >> "$LOG/$target.log" 2>&1 ;;
  esac
  grep -q '^git+' "$LOCKDIR/$name.pip.txt" 2>/dev/null && install_pip_git "$target" "$name"
  echo "  OK $target (UNPINNED: versions may differ from envs/$name.yml)"
}

install_pip_git() {
  GIT_CONFIG_COUNT=1 GIT_CONFIG_KEY_0=url.https://github.com/.insteadOf GIT_CONFIG_VALUE_0=git@github.com: \
    "$ENVS/$1/bin/python" -m pip install -q $(grep '^git+' "$LOCKDIR/$2.pip.txt") >> "$LOG/$1.log" 2>&1 \
    || warn "$1: git package failed (see $LOG/$1.log)"
}

if [[ $CHECK_ONLY == 0 ]]; then
  step "Creating conda environments"
  for yml in "$REPO"/envs/*.yml; do
    name="$(basename "$yml" .yml)"
    if printf '%s\n' "${RECREATE[@]:-}" | grep -qx "$name"; then
      echo "  removing $name$SUFFIX (--recreate)"; "$CONDA" env remove -y -n "$name$SUFFIX" >/dev/null 2>&1 || true
    fi
    if env_exists "$name$SUFFIX"; then echo "  exists: $name$SUFFIX"; else create_env "$name" || true; fi
  done
fi

# ---------------------------------------------------------------- 3. basecallers
if [[ $CHECK_ONLY == 0 && $SKIP_BASECALLERS == 0 ]]; then
  step "Installing Dorado $DORADO_VERSION"
  DARCH=$([[ "$ARCH" == aarch64 ]] && echo linux-arm64 || echo linux-x64)
  DDIR="$TOOLS/dorado-$DORADO_VERSION-$DARCH"
  if [[ ! -x "$DDIR/bin/dorado" ]]; then
    curl -fsSL -o "$TMPDIR/dorado.tar.gz" "https://cdn.oxfordnanoportal.com/software/analysis/dorado-$DORADO_VERSION-$DARCH.tar.gz"
    tar -xzf "$TMPDIR/dorado.tar.gz" -C "$TOOLS" && rm -f "$TMPDIR/dorado.tar.gz"
  else echo "  exists: $DDIR"; fi
  ln -sf "$DDIR/bin/dorado" "$BIN/dorado"

  step "Downloading Dorado R9.4.1 models"
  mkdir -p "$TOOLS/dorado_models"
  for m in "${DORADO_MODELS[@]}"; do
    if [[ -d "$TOOLS/dorado_models/$m" ]]; then echo "  exists: $m"
    else "$BIN/dorado" download --model "$m" --models-directory "$TOOLS/dorado_models" >> "$LOG/dorado.log" 2>&1 && echo "  OK $m" || warn "model $m failed"; fi
  done

  for GUPPY_VERSION in "${GUPPY_VERSIONS[@]}"; do
    step "Installing Guppy $GUPPY_VERSION (GPU build)"
    GDIR="$TOOLS/guppy-$GUPPY_VERSION"
    if [[ "$ARCH" != x86_64 ]]; then warn "Guppy archive is x86_64 only; skipped"
    elif [[ ! -x "$GDIR/bin/guppy_basecaller" ]]; then
      curl -fsSL -o "$TMPDIR/guppy.tar.gz" "https://cdn.oxfordnanoportal.com/software/analysis/ont-guppy_${GUPPY_VERSION}_linux64.tar.gz"
      mkdir -p "$GDIR" && tar -xzf "$TMPDIR/guppy.tar.gz" -C "$GDIR" --strip-components=1 && rm -f "$TMPDIR/guppy.tar.gz"
    else echo "  exists: $GDIR"; fi
  done
  # ~/bin/guppy_basecaller = newest; the pipeline calls each version by its full path
  [[ -x "$TOOLS/guppy-6.5.7/bin/guppy_basecaller" ]] && ln -sf "$TOOLS/guppy-6.5.7/bin/guppy_basecaller" "$BIN/guppy_basecaller"

  step "Downloading Bonito R9.4.1 models"
  BMODELS="$(ls -d "$ENVS"/env-bonito$SUFFIX/lib/python3*/site-packages/bonito/models 2>/dev/null | head -1)"
  if [[ -z "$BMODELS" ]]; then warn "env-bonito$SUFFIX not found; Bonito models skipped"
  else for m in "${BONITO_MODELS[@]}"; do
    if [[ -d "$BMODELS/$m" ]]; then echo "  exists: $m"
    else curl -fsSL -o "$TMPDIR/$m.zip" "https://cdn.oxfordnanoportal.com/software/analysis/bonito/$m.zip" \
         && "$ENVS/env-analysis$SUFFIX/bin/python" -m zipfile -e "$TMPDIR/$m.zip" "$BMODELS" && rm -f "$TMPDIR/$m.zip" && echo "  OK $m" || warn "Bonito model $m failed"; fi
  done; fi
fi

# ---------------------------------------------------------------- taxonomy
if [[ $CHECK_ONLY == 0 ]]; then
  step "NCBI taxonomy (taxdump) in $TAXONOMY_DIR"
  if [[ -s "$TAXONOMY_DIR/nodes.dmp" ]]; then echo "  exists"
  else mkdir -p "$TAXONOMY_DIR"
    curl -fsSL --retry 3 -o "$TMPDIR/taxdump.tar.gz" https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz \
      && tar -xzf "$TMPDIR/taxdump.tar.gz" -C "$TAXONOMY_DIR" nodes.dmp names.dmp merged.dmp delnodes.dmp \
      && echo "NCBI taxdump downloaded $(date -u +%F)" > "$TAXONOMY_DIR/SOURCE.txt" && rm -f "$TMPDIR/taxdump.tar.gz" \
      && echo "  OK" || warn "taxdump download failed"
  fi
fi

# ---------------------------------------------------------------- 4. ~/.bashrc
if [[ $CHECK_ONLY == 0 && $DO_BASHRC == 1 ]]; then
  step "Updating ~/.bashrc"
  touch "$HOME/.bashrc"
  [[ -f "$HOME/.bashrc.bak-fyp" ]] || cp "$HOME/.bashrc" "$HOME/.bashrc.bak-fyp"
  grep -q '>>> conda initialize >>>' "$HOME/.bashrc" || "$CONDA" init bash >/dev/null
  # Replace our block if present, then append the current version.
  sed -i '/^# >>> e21-fyp >>>/,/^# <<< e21-fyp <<</d' "$HOME/.bashrc"
  cat >> "$HOME/.bashrc" <<EOF
# >>> e21-fyp >>>   (managed by $REPO/setup.sh)
case ":\$PATH:" in *":\$HOME/bin:"*) ;; *) export PATH="\$HOME/bin:\$PATH" ;; esac
export FYP_CONDA_BASE="$CONDA_BASE"
export DORADO_MODELS="\$HOME/tools/dorado_models"
export TMPDIR="\$HOME/tmp"
# <<< e21-fyp <<<
EOF
  echo "  backup: ~/.bashrc.bak-fyp; open a new terminal or run: source ~/.bashrc"
fi

# ---------------------------------------------------------------- 5. check
step "Checking Python modules and project scripts"
PY_FAIL=0
for name in "${!PY_MODULES[@]}"; do
  py="$ENVS/$name$SUFFIX/bin/python"
  [[ -x "$py" ]] || { warn "$name$SUFFIX: no python"; PY_FAIL=1; continue; }
  mods="${PY_MODULES[$name]}"
  # every env that runs a bin/ script must import the project helper (stdlib only)
  if out=$(cd "$REPO/bin" && "$py" -c "import fyp_lib; [__import__(m) for m in '$mods'.split()]" 2>&1); then
    printf '  %-14s ok   %s\n' "$name$SUFFIX" "${mods:-(stdlib + fyp_lib)}"
  else
    warn "$name$SUFFIX: $(echo "$out" | tail -1)"; PY_FAIL=1
  fi
done
for script in "$REPO"/bin/*.py; do   # each script's --help runs its imports
  [[ "$(basename "$script")" == fyp_lib.py ]] && continue
  "$ENVS/env-analysis$SUFFIX/bin/python" "$script" --help > /dev/null 2>&1 || { warn "bin/$(basename "$script") --help failed"; PY_FAIL=1; }
done
[[ $PY_FAIL == 0 ]] && echo "  all bin/*.py scripts import cleanly"

step "Installed versions"
v() { printf '  %-14s %-15s ' "$1" "$2"; shift 2; ( "$@" 2>&1 || true ) | grep -v '^\s*$' | head -1; }
run_in() { local e="$1$SUFFIX"; shift; PATH="$ENVS/$e/bin:$PATH" "$@"; }
v - dorado "$BIN/dorado" --version
v - guppy bash -c "'$BIN/guppy_basecaller' --version | grep -o 'Version [^ ,]*'"
v env-basecall pod5 run_in env-basecall pod5 --version
v env-basecall samtools run_in env-basecall samtools --version
v env-qc NanoPlot run_in env-qc NanoPlot --version
v env-qc porechop_abi run_in env-qc porechop_abi --version
v env-qc seqkit run_in env-qc seqkit version
v env-align minimap2 run_in env-align minimap2 --version
v env-classify kraken2 run_in env-classify kraken2 --version
v env-classify bracken run_in env-classify bracken -v
v env-analysis blastn run_in env-analysis blastn -version
v env-analysis python-libs run_in env-analysis python -c "import pandas, pysam, statsmodels, lightgbm, shap, skbio, yaml; print('ok')"
v env-signal modkit run_in env-signal modkit --version
v env-emu emu run_in env-emu emu --version
v env-metamaps metamaps run_in env-metamaps metamaps
v env-assembly flye run_in env-assembly flye --version
v env-sim nanosim run_in env-sim simulator.py --version
v env-sim simreaduntil run_in env-sim python -c "import simreaduntil; print('ok')"
v - guppy-6.2.1 bash -c "'$TOOLS/guppy-6.2.1/bin/guppy_basecaller' --version | grep -o 'Version [^ ,]*'"
v env-bonito bonito run_in env-bonito bash -c "bonito --version 2>&1 | grep -v -i warn | grep -v pkg_resources | tail -1"
v env-qc cutadapt run_in env-qc cutadapt --version
v env-qc fastqc run_in env-qc fastqc --version
v env-align winnowmap run_in env-align winnowmap --version
v env-align ngmlr run_in env-align bash -c "ngmlr --version 2>&1 | grep -i version"
v env-hostile hostile run_in env-hostile hostile --version
v env-taxa kaiju run_in env-taxa bash -c "kaiju -h 2>&1 | grep -i 'kaiju 1'"
v env-taxa metabuli run_in env-taxa metabuli version
v env-taxa centrifuger run_in env-taxa centrifuger -v
v env-taxa krakenuniq run_in env-taxa krakenuniq --version
v env-taxa ganon run_in env-taxa ganon --version
v env-taxa kmcp run_in env-taxa kmcp version
v env-taxa sylph run_in env-taxa sylph --version
v env-protein diamond run_in env-protein diamond version
v env-analysis last run_in env-analysis lastal --version
v env-analysis mmseqs2 run_in env-analysis mmseqs version
v env-workflow nextflow bash -c "PATH='$ENVS/env-workflow$SUFFIX/bin':\$PATH nextflow -version | grep -o 'version [0-9.]*'"
echo
echo "Logs: $LOG"
echo "Next: edit paths in config/experiment.yaml (outdir, fast5, resources) for this account, then ./run_pipeline.sh config/test.yaml"
