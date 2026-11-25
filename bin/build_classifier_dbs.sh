#!/usr/bin/env bash
# Build a database for every classifier from ONE reference set, so that a
# difference between classifiers is the classifier, not its reference content.
#
#   build_classifier_dbs.sh --refset refset.tsv --genomes DIR --proteins DIR \
#                           --taxonomy DIR --out DIR [--threads 16] [--tools "kraken2 kaiju ..."]
#
# refset.tsv (tab-separated, no header): name  group  taxid  accession
#   genomes/<name>.fna   nucleotide FASTA      proteins/<name>.faa   protein FASTA (optional)
# --taxonomy: NCBI taxdump directory (nodes.dmp, names.dmp, merged.dmp)
#
# Writes <out>/<tool>/ for each tool, a layout main.nf expects, and prints the
# config/experiment.yaml block to paste (resources.classify_dbs).
# Needs the conda environments env-classify, env-taxa and env-protein.
set -euo pipefail

THREADS=16; TOOLS="kraken2 krakenuniq centrifuge centrifuger kaiju metabuli ganon taxor kmcp sylph diamond megan"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --refset) REFSET="$(realpath "$2")"; shift ;;  --genomes) GEN="$(realpath "$2")"; shift ;;
    --proteins) PROT="$(realpath "$2")"; shift ;;  --taxonomy) TAX="$(realpath "$2")"; shift ;;
    --out) OUT="$2"; shift ;;  --threads) THREADS="$2"; shift ;;  --tools) TOOLS="$2"; shift ;;
    -h|--help) sed -n '2,15p' "$0"; exit 0 ;;  *) echo "unknown option $1"; exit 1 ;;
  esac; shift
done
: "${REFSET:?--refset}" "${GEN:?--genomes}" "${TAX:?--taxonomy}" "${OUT:?--out}"
BASE="${FYP_CONDA_BASE:-$HOME/anaconda3}"
export PATH="$BASE/envs/env-classify/bin:$PATH:$BASE/envs/env-taxa/bin:$BASE/envs/env-protein/bin"
export TMPDIR="${TMPDIR:-$HOME/tmp}"
mkdir -p "$OUT/_prep"; OUT="$(realpath "$OUT")"; P="$OUT/_prep"
log() { printf '[%(%H:%M:%S)T] %s\n' -1 "$*"; }

# ---- shared inputs --------------------------------------------------------------
log "preparing shared inputs"
: > "$P/refs.fa"; : > "$P/kraken.fa"; : > "$P/seqid2taxid.tsv"; : > "$P/genome2taxid.tsv"
: > "$P/prot_kaiju.faa"; : > "$P/prot_plain.faa"
printf 'accession\taccession.version\ttaxid\tgi\n' > "$P/nucl.accession2taxid"
printf 'accession\taccession.version\ttaxid\tgi\n' > "$P/prot.accession2taxid"
mkdir -p "$P/genomes"; pn=0
while IFS=$'\t' read -r name group taxid acc; do
  [[ -z "$name" || "$name" == \#* ]] && continue
  f="$GEN/$name.fna"; [[ -f "$f" ]] || { echo "missing $f"; exit 1; }
  ln -sf "$f" "$P/genomes/$name.fna"
  printf '%s\t%s\n' "$P/genomes/$name.fna" "$taxid" >> "$P/genome2taxid.tsv"
  awk -v t="$taxid" -v k="$P/kraken.fa" -v m="$P/seqid2taxid.tsv" -v a="$P/nucl.accession2taxid" '
    /^>/ { id = substr($1, 2); print ">" id >> k; print ">" id;
           print id "\t" t >> m; acc = id; sub(/\.[0-9]+$/, "", acc); print acc "\t" id "\t" t "\t0" >> a; next }
         { print; print >> k; }' "$f" >> "$P/refs.fa"
  if [[ -n "${PROT:-}" && -f "$PROT/$name.faa" ]]; then
    while IFS= read -r line; do
      if [[ "$line" == ">"* ]]; then pn=$((pn + 1))
        printf '>P%d_%s\n' "$pn" "$taxid" >> "$P/prot_kaiju.faa"; printf '>P%d\n' "$pn" >> "$P/prot_plain.faa"
        printf 'P%d\tP%d\t%s\t0\n' "$pn" "$pn" "$taxid" >> "$P/prot.accession2taxid"
      else printf '%s\n' "$line" >> "$P/prot_kaiju.faa"; printf '%s\n' "$line" >> "$P/prot_plain.faa"; fi
    done < "$PROT/$name.faa"
  fi
done < "$REFSET"
# kraken-style headers: >seqid|kraken:taxid|T
awk -F'\t' 'NR==FNR {t[$1]=$2; next} /^>/ {id=substr($1,2); print ">" id "|kraken:taxid|" t[id]; next} {print}' \
  "$P/seqid2taxid.tsv" "$P/kraken.fa" > "$P/kraken_tagged.fa"
log "$(grep -c '>' "$P/refs.fa") sequences, $(grep -c '>' "$P/prot_plain.faa" || true) proteins"
YAML="$OUT/classify_dbs.yaml"; echo "# paste under resources.classify_dbs in config/experiment.yaml" > "$YAML"
TAG="$(basename "$OUT")"
add_yaml() { printf '  %s_%s:\n    classifier: %s\n    path: %s\n' "$1" "$TAG" "$1" "$OUT/$1" >> "$YAML"; }

build() {
  local tool="$1" d="$OUT/$1"; rm -rf "$d"; mkdir -p "$d"
  log "building $tool"
  case "$tool" in
    kraken2)
      mkdir -p "$d/taxonomy"; cp "$TAX"/{nodes,names}.dmp "$d/taxonomy/"
      kraken2-build --db "$d" --add-to-library "$P/kraken_tagged.fa" --no-masking > "$d/build.log" 2>&1
      kraken2-build --db "$d" --build --threads "$THREADS" >> "$d/build.log" 2>&1 ;;
    krakenuniq)
      mkdir -p "$d/taxonomy" "$d/library"; cp "$TAX"/{nodes,names}.dmp "$d/taxonomy/"
      cp "$P/refs.fa" "$d/library/refs.fna"; cp "$P/seqid2taxid.tsv" "$d/library/refs.map"
      krakenuniq-build --db "$d" --threads "$THREADS" --kmer-len 31 --jellyfish-bin "$BASE/envs/env-taxa/bin/jellyfish" > "$d/build.log" 2>&1 ;;
    centrifuge)
      centrifuge-build -p "$THREADS" --conversion-table "$P/seqid2taxid.tsv" --taxonomy-tree "$TAX/nodes.dmp" \
        --name-table "$TAX/names.dmp" "$P/refs.fa" "$d/index" > "$d/build.log" 2>&1 ;;
    centrifuger)
      centrifuger-build -t "$THREADS" -r "$P/refs.fa" --conversion-table "$P/seqid2taxid.tsv" \
        --taxonomy-tree "$TAX/nodes.dmp" --name-table "$TAX/names.dmp" -o "$d/index" > "$d/build.log" 2>&1 ;;
    kaiju)
      [[ -s "$P/prot_kaiju.faa" ]] || { log "kaiju skipped: no proteins"; return; }
      kaiju-mkbwt -n "$THREADS" -a ACDEFGHIKLMNPQRSTVWY -o "$d/kaiju_db" "$P/prot_kaiju.faa" > "$d/build.log" 2>&1
      kaiju-mkfmi "$d/kaiju_db" >> "$d/build.log" 2>&1 && rm -f "$d/kaiju_db.bwt" "$d/kaiju_db.sa" ;;
    metabuli)
      ls "$P"/genomes/*.fna > "$P/fasta_list.txt"
      metabuli build "$d" "$P/fasta_list.txt" "$P/nucl.accession2taxid" --taxonomy-path "$TAX" \
        --threads "$THREADS" --max-ram 32 > "$d/build.log" 2>&1 ;;
    ganon)
      awk -F'\t' '{n=$1; sub(/.*\//,"",n); sub(/\.fna$/,"",n); print $1 "\t" n "\t" $2}' "$P/genome2taxid.tsv" > "$d/input.tsv"
      cut -f2,3 "$d/input.tsv" > "$d/target2taxid.tsv"
      ganon build-custom --input-file "$d/input.tsv" --db-prefix "$d/db" --taxonomy ncbi \
        --taxonomy-files "$TAX/nodes.dmp" "$TAX/names.dmp" --threads "$THREADS" > "$d/build.log" 2>&1 ;;
    taxor)
      # taxor matches files by accession-style IDs (one underscore), so use REF_<i>
      mkdir -p "$d/seqs"; : > "$d/input.tsv"; : > "$d/ref2name.tsv"; local i=0 gf gt
      while IFS=$'\t' read -r gf gt; do i=$((i + 1)); local n="REF_$i"
        gzip -c "$gf" > "$d/seqs/${n}_genomic.fna.gz"
        printf '%s\t%s\n' "$n" "$(basename "$gf" .fna)" >> "$d/ref2name.tsv"
        lin=$(echo "$gt" | python3 "$(dirname "$0")/taxor_lineage.py" --taxonomy "$TAX" | cut -f2,3)
        printf '%s\t%s\tftp/%s\t%s\t%s\n' "$n" "$gt" "$n" "$n" "$lin" >> "$d/input.tsv"
      done < "$P/genome2taxid.tsv"
      taxor build --input-file "$d/input.tsv" --input-sequence-dir "$d/seqs" --output-filename "$d/taxor.hixf" \
        --threads "$THREADS" --kmer-size 22 --syncmer-size 12 --use-syncmer > "$d/build.log" 2>&1 ;;
    kmcp)
      awk -F'\t' '{n=$1; sub(/.*\//,"",n); sub(/\.fna$/,"",n); print n "\t" $2}' "$P/genome2taxid.tsv" > "$d/taxid.map"
      kmcp compute -I "$P/genomes" -O "$d/tmp" -k 21 -n 10 -l 150 -j "$THREADS" --force > "$d/build.log" 2>&1
      kmcp index -I "$d/tmp" -O "$d/kmcp.db" -j "$THREADS" --force >> "$d/build.log" 2>&1 && rm -rf "$d/tmp" ;;
    sylph)
      sylph sketch -t "$THREADS" -o "$d/db" "$P"/genomes/*.fna > "$d/build.log" 2>&1
      awk -F'\t' '{n=$1; sub(/.*\//,"",n); print n "\t" $2}' "$P/genome2taxid.tsv" > "$d/genome2taxid.tsv" ;;
    diamond|megan)
      [[ -s "$P/prot_plain.faa" ]] || { log "$tool skipped: no proteins"; return; }
      diamond makedb --in "$P/prot_plain.faa" -d "$d/db" -p "$THREADS" --taxonmap "$P/prot.accession2taxid" \
        --taxonnodes "$TAX/nodes.dmp" --taxonnames "$TAX/names.dmp" > "$d/build.log" 2>&1
      cut -f2,3 "$P/prot.accession2taxid" | tail -n +2 > "$d/acc2taxid.tsv" ;;
    *) echo "unknown tool $tool"; exit 1 ;;
  esac
  # verify the database exists before declaring success
  local need
  case "$tool" in
    kraken2) need="$d/hash.k2d" ;;           krakenuniq) need="$d/database.kdb" ;;
    centrifuge) need="$d/index.1.cf" ;;      centrifuger) need="$d/index.1.cfr" ;;
    kaiju) need="$d/kaiju_db.fmi" ;;         metabuli) need="$d/diffIdx" ;;
    ganon) need="$d/db.hibf" ;;              taxor) need="$d/taxor.hixf" ;;
    kmcp) need="$d/kmcp.db" ;;               sylph) need="$d/db.syldb" ;;
    diamond|megan) need="$d/db.dmnd" ;;
  esac
  [[ -e "$need" ]] || { log "FAILED $tool: $need missing; last lines of $d/build.log:"; tail -5 "$d/build.log"; return 1; }
  add_yaml "$tool"
}

FAILED=()
for tool_name in $TOOLS; do build "$tool_name" || FAILED+=("$tool_name"); done
log "done; failed: ${FAILED[*]:-none}"
cat "$YAML"
