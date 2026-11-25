# Computational Modelling of Bias in Nanopore DNA Sequencing

Final-year project (e21, Department of Computer Engineering, University of Peradeniya).
Phase 1 asks where the host–microbe composition bias arises in a MinION run of
green seaweed DNA (R9.4.1, SQK-LSK110): basecalling, trimming, quality filtering,
host mapping or classification. The literature review is in this repository.

## How the pipeline is organised

Every pipeline choice is a **factor** with several **levels**. The pipeline runs
every combination of levels and stores each one in its own folder, so any two
combinations can be compared later.

| Stage | Factor (folder key) | Levels in `config/experiment.yaml` |
|---|---|---|
| Basecall | `basecaller`, `model` | guppy-6.5.7, dorado-0.9.6 × fast, hac, sup (+ optional imported MinKNOW / Guppy 6.2.1 calls) |
| Trim | `trim` | none, porechop_abi (porechop also supported) |
| QC | `minq`, `minlen` | 0, 7, 8, 9, 10 × 0, 500 (0 = no filter) |
| Host mapping | `aligner`, `host_ref` | minimap2 × host references |
| Classification | `classifier`, `db` | kraken2 or centrifuge × databases |
| BLAST of reads | `blast_tool`, `blast_db` | megablast, dc-megablast, blastn, minimap2 × labelled databases |

Fixed settings (max read length 50 kb for artefacts, MAPQ 20, BLAST subsample
2,000 reads with seed 11) are recorded in every `metadata.json`.

Mean read Q is computed the same way for every basecaller (`bin/fyp_lib.py`:
the Phred value of the mean per-base error probability). The basecallers run
with no Q filter (`--disable_qscore_filtering`, `--min-qscore 0`), so the QC
stage alone decides pass and fail, and **failed reads are always kept**.

## Outputs folder

```
<outdir>/                                   (default /home/induwara/dna_data/FYP/outputs)
  basecaller=guppy-6.5.7/model=hac/         reads.fastq.gz, sequencing_summary.tsv.gz, reads_table.tsv.gz, metadata.json, results.tsv
    trim=porechop_abi/                      reads.fastq.gz, reads_table.tsv.gz, metadata.json, results.tsv
      minq=9/minlen=0/                      passed / failed / artefact .fastq.gz, reads_table.tsv.gz, metadata.json, results.tsv
        analysis=host_map/aligner=minimap2/host_ref=<ref>/      passed.bam, failed.bam, host_reads.tsv.gz, metadata.json, results.tsv
        analysis=classify/classifier=kraken2/db=<db>/           *.kreport, top_taxa.tsv, metadata.json, results.tsv
        analysis=blast/blast_tool=megablast/blast_db=<db>/      *.sub.fa, *.hits, best_hits.tsv, metadata.json, results.tsv
  _shared/inventory/                        inventory.csv (per FAST5 file), raw_reads.tsv.gz (channel, mux, time per read)
  _runs/<UTC time>/                         design.csv, experiment.yaml copy, run_info.yaml (git commit), trace, report, timeline
  _summary/                                 nodes.tsv, results_long.tsv / .parquet (every node, every metric)
```

* Each folder name is `factor=level`, so the path alone identifies the combination
  (`pandas`/`pyarrow` can read it as Hive-style partitions).
* Each upstream stage is computed **once** and shared by every downstream combination.
* Every node has `metadata.json` (factors, exact command, tool versions, fixed
  parameters, parent node, md5 of every output) and `results.tsv`
  (`read_set`, `metric`, `value`; read sets are `all`, `passed`, `failed`, `artefact`).
* Files are hard links into Nextflow's work directory (`nf-work/`), so they take no extra space.

## Running

```bash
./run_pipeline.sh config/test.yaml                    # small test: 1 FAST5 file, fast models (~5 min)
./run_pipeline.sh                                     # full design from config/experiment.yaml
./run_pipeline.sh config/experiment.yaml --mode ofat  # baseline + one factor at a time
```

`run_pipeline.sh` (1) builds `design.csv` with `bin/design.py`, (2) runs `main.nf`
with Nextflow using `-resume`, so finished stages are reused, and (3) gathers every
result with `bin/collect_results.py`. Resources that do not exist yet (for
example a database not downloaded) are skipped with a warning, so early stages can
run first. To add a level, add it to `config/experiment.yaml` and rerun: only the
new combinations are computed.

Then estimate effects:

```bash
python bin/effects.py --summary <outdir>/_summary/results_long.tsv --config config/experiment.yaml \
    --stage qc --read-set all --metric pass_rate_reads --logit
```

This prints one-factor-at-a-time contrasts against the baseline (the configuration
of the original run: Guppy HAC, no trimming, Q9, no length filter) and a factorial
model with every factor as a main effect (`--interactions` adds two-way terms).

## Repository layout

| Path | Purpose |
|---|---|
| `main.nf`, `nextflow.config` | Nextflow workflow (execution, caching, publishing) |
| `run_pipeline.sh` | One command: design → run → collect |
| `config/experiment.yaml` | The experiment: factors, levels, baseline, resources |
| `config/test.yaml` | Small test configuration |
| `bin/design.py` | Builds the combination list (full factorial or OFAT) |
| `bin/qc_split.py`, `bin/read_stats.py` | QC split keeping failed reads; per-read Q, length, GC |
| `bin/host_map_summary.py`, `bin/kraken_summary.py`, `bin/blast_summary.py` | Stage summaries → `results.tsv` |
| `bin/write_metadata.py` | `metadata.json` for every node |
| `bin/make_inventory.py` | Dataset inventory and per-read raw metadata from FAST5 |
| `bin/build_blast_db.py` | Labelled BLAST database (`group|id` headers) |
| `bin/collect_results.py`, `bin/effects.py` | Summary tables and effect estimates |
| `envs/*.yml` | Conda environments (`conda env create -f envs/env-qc.yml`) |

## Environments and tools

**New machine or new team member:** clone the repo and run `./setup.sh` (no sudo).
It finds or installs conda (Miniforge), sets the conda-forge + bioconda channels,
builds every environment, installs Dorado, Guppy and the R9.4.1 models, and adds a
marked block to `~/.bashrc` (`~/bin` on PATH, `FYP_CONDA_BASE`, `DORADO_MODELS`,
`TMPDIR`). It is safe to re-run.

Environments are rebuilt from, in order of exactness:

| File | What it pins |
|---|---|
| `envs/lock/<env>.linux-64.lock` | every conda package: exact build and download URL, md5-checked (Linux x86_64) |
| `envs/lock/<env>.pip.txt` | every pip-installed Python module, `name==version`; SimReadUntil at a fixed git commit |
| `envs/<env>.yml` | versions only (fallback for other Linux setups) |

On Linux x86_64 the result is package-for-package identical (tested by rebuilding
all 11 environments and comparing every name, version, build and channel). The
script then imports every Python module the project uses in the environment that
runs it and runs each `bin/*.py --help`.

Options: `--check` (versions only), `--skip-basecallers`, `--recreate <env>`,
`--no-bashrc`, `--env-suffix -test` (rebuild under new names to test).
**After installing or upgrading anything, run `./setup.sh --export` and commit `envs/`.**
Not covered (system level): the NVIDIA driver, git and curl.

Conda environments (`envs/`) plus two binaries in `~/bin`: Dorado 0.9.6 (the last
release line with R9.4.1 models; it reads FAST5 directly) and Guppy 6.5.7. Models
are in `~/tools/dorado_models`. In `env-sim`, SimReadUntil is installed from GitHub
(`pip install git+https://github.com/ratschlab/sim_read_until.git`, with
`git config url.https://github.com/.insteadOf git@github.com:` for its submodule).

Known issue: `pod5 convert fast5` fails on these files (the pip VBZ plugin cannot
decompress them), so POD5 is not used; Dorado reads the FAST5 files directly.
