#!/usr/bin/env bash
# Download / build the reference databases (levels of the refdb factor) relevant to the seaweed sample.
#
#   $FYP_DB/refdb/<refdb>/<tool>/   one folder per classifier that can use that reference database
#   $FYP_DB/refdb/<refdb>/search/   labelled read-search databases (BLAST, LAST, MMseqs2, DIAMOND, FASTA for minimap2)
#
# Levels (config/experiment.yaml, resources.refdbs):
#   refset         host Ulva + other algae + seaweed-associated marine bacteria (bin/refset_taxa.tsv), NCBI
#                  reference genomes, built for EVERY classifier and search tool (fair comparison)
#   refset_nohost  the same without the host group: what classifiers do when the host is missing
#   standard_8gb   Kraken2 Standard (RefSeq archaea, bacteria, viral, plasmid, human) capped at 8 GB (prebuilt)
#   standard_16gb  Kraken2 Standard capped at 16 GB
#   pluspfp_16gb   Kraken2 PlusPFP capped at 16 GB (Standard + protozoa, fungi, plants incl. green algae)
#   core_nt        NCBI core_nt BLAST database (the lab protocol's "nucleotide database"; 263 GB, not in "all")
# plus host_ref ulva_refset: all Ulva genomes of the reference set as one host reference for host_map.
#
# Usage: bin/download_databases.sh <refset|refset_nohost|host_ref|kraken2|core_nt|all> [taxa.tsv]
#   Every step skips what is finished; rerun to resume. Each database gets SOURCE.txt.
#   To run it detached (keeps going after the SSH session closes):
#     setsid nohup bin/download_databases.sh all > ~/dna_data/FYP/databases/download.log 2>&1 < /dev/null &
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DB="${FYP_DB:-$HOME/dna_data/FYP/databases}"
R="$DB/refdb"
BASE="${FYP_CONDA_BASE:-$HOME/anaconda3}"
A="$BASE/envs/env-analysis/bin"
export TMPDIR="${TMPDIR:-$HOME/tmp}"
THREADS="${THREADS:-16}"
mkdir -p "$R" "$TMPDIR"
log() { printf '[%(%F %H:%M:%S)T] %s\n' -1 "$*"; }
K2=https://genome-idx.s3.amazonaws.com/kraken

kraken2_db() {   # kraken2_db LEVEL FILE   (resumable download, then extract)
  local level="$1" file="$2" out="$R/$1/kraken2"
  [[ -s "$out/hash.k2d" ]] && { log "$level/kraken2 exists"; return; }
  log "$level/kraken2: downloading $file"
  curl -fL --retry 10 --retry-delay 30 -C - -o "$TMPDIR/$file" "$K2/$file"
  mkdir -p "$out" && tar -xzf "$TMPDIR/$file" -C "$out" && rm -f "$TMPDIR/$file"
  echo "$K2/$file downloaded $(date -u +%F)" > "$out/SOURCE.txt"
}

kraken2() {
  kraken2_db standard_8gb  k2_standard_08_GB_20261022.tar.gz
  kraken2_db standard_16gb k2_standard_16_GB_20261022.tar.gz
  kraken2_db pluspfp_16gb  k2_pluspfp_16_GB_20261022.tar.gz
}

core_nt() {
  mkdir -p "$R/core_nt/search" && cd "$R/core_nt/search"
  log "core_nt: update_blastdb.pl (263 GB compressed; rerun to resume)"
  "$A/update_blastdb.pl" --decompress --num_threads 4 core_nt taxdb
  echo "NCBI core_nt + taxdb via update_blastdb.pl, $(date -u +%F)" > SOURCE.txt
}

fetch_refset() {   # NCBI reference genomes + proteins for every taxon in the list
  local taxa="${1:-$REPO/bin/refset_taxa.tsv}" S="$DB/refset_genomes"
  [[ -s "$S/refset.tsv" ]] && { log "refset genomes exist ($(wc -l < "$S/refset.tsv"))"; return; }
  mkdir -p "$S/genomes" "$S/proteins" "$S/zips"
  : > "$S/refset.tsv.part"
  while IFS=$'\t' read -r group taxon; do
    [[ -z "$group" || "$group" == \#* ]] && continue
    local z="$S/zips/$(echo "$taxon" | tr ' ' '_').zip" d acc tx
    log "$group: $taxon"
    [[ -s "$z" ]] || "$A/datasets" download genome taxon "$taxon" --reference --include genome,protein --filename "$z" > /dev/null
    rm -rf "$TMPDIR/ds" && unzip -oq "$z" -d "$TMPDIR/ds"
    for d in "$TMPDIR"/ds/ncbi_dataset/data/GC*; do
      [[ -d "$d" ]] || continue
      acc="$(basename "$d")"
      cat "$d"/*.fna > "$S/genomes/$acc.fna"
      [[ -s "$d/protein.faa" ]] && cp "$d/protein.faa" "$S/proteins/$acc.faa"
      tx="$("$A/datasets" summary genome accession "$acc" --as-json-lines | "$A/dataformat" tsv genome --fields organism-tax-id | tail -1)"
      printf '%s\t%s\t%s\t%s\n' "$acc" "$group" "$tx" "$acc" >> "$S/refset.tsv.part"
    done
  done < "$taxa"
  mv "$S/refset.tsv.part" "$S/refset.tsv"
  echo "NCBI datasets reference genomes for $taxa, $(date -u +%F)" > "$S/SOURCE.txt"
}

build_refdb() {    # build_refdb LEVEL REFSET_TSV: every classifier + the labelled search databases
  local level="$1" tsv="$2" S="$DB/refset_genomes" out="$R/$1" args=() g
  log "$level: $(wc -l < "$tsv") genomes -> every classifier"
  "$REPO/bin/build_classifier_dbs.sh" --refset "$tsv" --genomes "$S/genomes" --proteins "$S/proteins" \
    --taxonomy "$DB/taxonomy" --out "$out" --threads "$THREADS"
  log "$level: search databases"
  for g in $(cut -f2 "$tsv" | sort -u); do
    mkdir -p "$out/_by_group/$g/g" "$out/_by_group/$g/p"
    awk -F'\t' -v g="$g" '$2 == g {print $1}' "$tsv" | while read -r acc; do
      ln -sf "$S/genomes/$acc.fna" "$out/_by_group/$g/g/"
      if [[ -s "$S/proteins/$acc.faa" ]]; then ln -sf "$S/proteins/$acc.faa" "$out/_by_group/$g/p/"; fi
    done
    args+=(--group "$g=$out/_by_group/$g/g/*.fna")
    ls "$out/_by_group/$g/p/"*.faa > /dev/null 2>&1 && args+=(--protein-group "$g=$out/_by_group/$g/p/*.faa")
  done
  mkdir -p "$out/search"
  PATH="$A:$PATH" "$A/python" "$REPO/bin/build_blast_db.py" --out "$out/search/refdb" --threads "$THREADS" "${args[@]}"
  cp "$tsv" "$out/refset.tsv"
  log "$level: done"
}

refset() { fetch_refset "${1:-}"; build_refdb refset "$DB/refset_genomes/refset.tsv"; }

refset_nohost() {
  fetch_refset "${1:-}"
  awk -F'\t' '$2 != "host"' "$DB/refset_genomes/refset.tsv" > "$DB/refset_genomes/refset_nohost.tsv"
  build_refdb refset_nohost "$DB/refset_genomes/refset_nohost.tsv"
}

host_ref() {
  fetch_refset "${1:-}"
  local out="$HOME/dna_data/FYP/references/ulva_refset"
  mkdir -p "$out"
  awk -F'\t' '$2 == "host" {print $1}' "$DB/refset_genomes/refset.tsv" | \
    while read -r acc; do cat "$DB/refset_genomes/genomes/$acc.fna"; done > "$out/ulva_refset.fa"
  awk -F'\t' '$2 == "host"' "$DB/refset_genomes/refset.tsv" > "$out/SOURCE.tsv"
  log "host_ref ulva_refset: $(grep -c '>' "$out/ulva_refset.fa") sequences"
}

case "${1:-}" in
  refset|refset_nohost|host_ref) "$1" "${2:-}" ;;
  kraken2) kraken2 ;;
  core_nt) core_nt ;;
  all) kraken2; fetch_refset "${2:-}"; host_ref "${2:-}"; refset "${2:-}"; refset_nohost "${2:-}" ;;
  *) sed -n '2,22p' "$0"; exit 1 ;;
esac
log "finished: $1"
