#!/usr/bin/env python3
"""Seven-rank lineage strings in the format taxor build expects.

Reads taxids (one per line) on stdin and prints, per taxid:
  taxid <tab> k__..;p__..;c__..;o__..;f__..;g__..;s__.. <tab> id;id;id;id;id;id;id
The top rank is 'superkingdom' or, in NCBI taxonomy from 2024 on, 'domain'.
A missing rank is left empty (k__ with id 0 is never written: the nearest
higher taxon's id is repeated, as taxonkit reformat -F does).

Usage: echo 63186 | taxor_lineage.py --taxonomy DIR
"""
import argparse
import os
import sys

RANKS = [("k", ("superkingdom", "domain", "kingdom")), ("p", ("phylum",)), ("c", ("class",)),
         ("o", ("order",)), ("f", ("family",)), ("g", ("genus",)), ("s", ("species",))]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--taxonomy", required=True)
    a = ap.parse_args()
    wanted = [line.strip() for line in sys.stdin if line.strip()]
    parent, rank, name = {}, {}, {}
    with open(os.path.join(a.taxonomy, "nodes.dmp")) as fh:
        for line in fh:
            f = line.split("\t|\t", 3)
            parent[f[0]], rank[f[0]] = f[1], f[2]
    with open(os.path.join(a.taxonomy, "names.dmp")) as fh:
        for line in fh:
            f = line.split("\t|\t")
            if f[3].startswith("scientific name"):
                name[f[0]] = f[1]
    for t in wanted:
        lineage, x = [], t
        while x in parent and x != "1":
            lineage.append(x)
            x = parent[x]
        by_rank = {rank[x]: x for x in lineage}
        names, ids, last = [], [], "1"
        for prefix, accepted in RANKS:
            hit = next((by_rank[r] for r in accepted if r in by_rank), None)
            names.append(f"{prefix}__{name.get(hit, '') if hit else ''}")
            last = hit or last
            ids.append(last)
        print(f"{t}\t{';'.join(names)}\t{';'.join(ids)}")


if __name__ == "__main__":
    main()
