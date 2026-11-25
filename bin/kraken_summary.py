#!/usr/bin/env python3
"""Summarise Kraken-style reports (Kraken2, or Centrifuge via centrifuge-kreport).

Reports the classified share and the share of reads under key clades, so host
(green algae) and microbial fractions can be compared across combinations.
Shares are of ALL reads, classified or not: if host reads are left
unclassified because the database lacks the host, that shows up here as a
low Viridiplantae share and a high unclassified share.

Usage: kraken_summary.py --report passed=passed.kreport --report failed=failed.kreport
"""
import argparse
import math

from fyp_lib import write_results

CLADES = {
    "bacteria": "2", "archaea": "2157", "viruses": "10239", "eukaryota": "2759",
    "viridiplantae": "33090", "chlorophyta": "3041", "rhodophyta": "2763",
    "stramenopiles": "33634", "fungi": "4751", "human": "9606",
}


def parse_report(path):
    """Return ({taxid: clade_reads}, unclassified, total, [(rank, taxid, name, reads)])."""
    clade, rows, unclassified, root = {}, [], 0, 0
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            reads, rank, taxid, name = int(f[1]), f[-3], f[-2], f[-1].strip()
            clade[taxid] = reads
            rows.append((rank, taxid, name, reads))
            if taxid == "0":
                unclassified = reads
            elif taxid == "1":
                root = reads
    return clade, unclassified, unclassified + root, rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--report", action="append", required=True, help="read_set=path.kreport")
    ap.add_argument("--results", default="results.tsv")
    ap.add_argument("--top", default="top_taxa.tsv")
    ap.add_argument("--top-n", type=int, default=20)
    a = ap.parse_args()

    rows = []
    with open(a.top, "w") as top:
        top.write("read_set\trank\ttaxid\tname\treads\tfrac_of_all_reads\n")
        for spec in a.report:
            read_set, path = spec.split("=", 1)
            clade, uncl, total, taxa = parse_report(path)
            rows.append((read_set, "n_reads", total))
            rows.append((read_set, "frac_classified", round(1 - uncl / total, 5) if total else ""))
            for label, taxid in CLADES.items():
                rows.append((read_set, f"frac_{label}", round(clade.get(taxid, 0) / total, 5) if total else ""))
            host, microbe = clade.get("33090", 0), clade.get("2", 0) + clade.get("2157", 0)
            rows.append((read_set, "log_ratio_viridiplantae_prokaryote",
                         round(math.log((host + 0.5) / (microbe + 0.5)), 5)))
            for rank in ("G", "S"):
                best = sorted((t for t in taxa if t[0] == rank), key=lambda t: -t[3])[: a.top_n]
                for r, taxid, name, reads in best:
                    top.write(f"{read_set}\t{r}\t{taxid}\t{name}\t{reads}\t{reads / total if total else 0:.5f}\n")
    write_results(a.results, rows)


if __name__ == "__main__":
    main()
