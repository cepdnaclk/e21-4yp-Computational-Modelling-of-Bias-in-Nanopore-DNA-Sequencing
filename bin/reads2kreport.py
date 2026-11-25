#!/usr/bin/env python3
"""Convert any classifier's output to a Kraken-style report (kreport).

Every classifier writes its result differently. To compare them fairly, each
one is reduced to the same thing: one NCBI taxid per read (0 = unclassified),
or, for profilers that do not classify single reads (sylph), an abundance per
reference genome. This script rolls those up the NCBI taxonomy and writes a
kreport, which kraken_summary.py then reads exactly as it reads Kraken2's own.

  --assign file.tsv                     read_id <tab> taxid   (one line per read)
  --cami binning.tsv                    CAMI binning format (kmcp, taxor): @@SEQUENCEID ... TAXID
  --sylph profile.tsv --genome-map map  sylph profile + 'genome_file<tab>taxid'
  --total N                             reads in the input (unlisted reads = unclassified)
  --taxonomy DIR                        nodes.dmp, names.dmp (+ merged.dmp) from NCBI taxdump
"""
import argparse
import os
from collections import defaultdict

RANK_CODE = {"superkingdom": "D", "domain": "D", "kingdom": "K", "phylum": "P", "class": "C",
             "order": "O", "family": "F", "genus": "G", "species": "S"}


def load_taxonomy(d):
    parent, rank, name, merged = {}, {}, {}, {}
    with open(os.path.join(d, "nodes.dmp")) as fh:
        for line in fh:
            f = line.split("\t|\t")
            parent[f[0]] = f[1]
            rank[f[0]] = f[2]
    with open(os.path.join(d, "names.dmp")) as fh:
        for line in fh:
            f = line.split("\t|\t")
            if f[3].startswith("scientific name"):
                name[f[0]] = f[1]
    m = os.path.join(d, "merged.dmp")
    if os.path.exists(m):
        with open(m) as fh:
            for line in fh:
                f = line.split("\t|\t")
                merged[f[0]] = f[1].split("\t")[0]
    return parent, rank, name, merged


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--assign")
    ap.add_argument("--cami")
    ap.add_argument("--sylph")
    ap.add_argument("--genome-map")
    ap.add_argument("--total", type=int, required=True)
    ap.add_argument("--taxonomy", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    parent, rank, name, merged = load_taxonomy(a.taxonomy)
    direct = defaultdict(float)
    if a.assign:
        with open(a.assign) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 2 or f[1] in ("", "0", "NA", "-"):
                    continue
                t = merged.get(f[1], f[1])
                direct[t if t in parent else "1"] += 1
    elif a.cami:
        sid, tid, seen_reads = 0, None, set()
        with open(a.cami) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line.startswith("@@"):
                    cols = [c.upper() for c in line[2:].split("\t")]
                    sid, tid = cols.index("SEQUENCEID"), cols.index("TAXID")
                    continue
                if not line or line.startswith(("@", "#")) or tid is None:
                    continue
                f = line.split("\t")
                if f[sid] in seen_reads or len(f) <= tid or f[tid] in ("", "0", "NA"):
                    continue
                seen_reads.add(f[sid])      # one taxon per read (first listed)
                t = merged.get(f[tid], f[tid])
                direct[t if t in parent else "1"] += 1
    elif a.sylph:
        gmap = dict(line.rstrip("\n").split("\t")[:2] for line in open(a.genome_map) if line.strip())
        gmap = {os.path.basename(k): v for k, v in gmap.items()}
        with open(a.sylph) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            gi, ai = header.index("Genome_file"), header.index("Sequence_abundance")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                t = gmap.get(os.path.basename(f[gi]))
                if t:   # sequence abundance (%) of all reads -> pseudo read count
                    direct[merged.get(t, t)] += float(f[ai]) / 100 * a.total
    else:
        raise SystemExit("give --assign, --cami or --sylph")

    clade = defaultdict(float)
    for t, c in direct.items():
        seen = set()
        while t not in seen:
            seen.add(t)
            clade[t] += c
            if t == "1" or t not in parent:
                break
            t = parent[t]
    classified = clade.get("1", 0.0)
    unclassified = max(a.total - classified, 0)
    children = defaultdict(list)
    for t in clade:
        if t != "1" and t in parent:
            children[parent[t]].append(t)

    total = max(a.total, 1)
    with open(a.out, "w") as out:
        out.write(f"{100 * unclassified / total:.2f}\t{round(unclassified)}\t{round(unclassified)}\tU\t0\tunclassified\n")

        def walk(t, depth):
            code = "R" if t == "1" else RANK_CODE.get(rank.get(t, ""), "-")
            out.write(f"{100 * clade[t] / total:.2f}\t{round(clade[t])}\t{round(direct.get(t, 0))}\t{code}\t{t}\t"
                      f"{'  ' * depth}{name.get(t, 'taxid ' + t)}\n")
            for c in sorted(children[t], key=lambda x: -clade[x]):
                walk(c, depth + 1)
        if classified:
            walk("1", 0)


if __name__ == "__main__":
    main()
