#!/usr/bin/env python3
"""Summarise BLAST (outfmt 6) or minimap2 (PAF) hits of subsampled reads.

Database sequences must be labelled 'group|original_id' (build_blast_db.py
does this), so each read's best hit can be assigned to a group such as
'host', 'algae' or 'bacteria'. The key number for RQ1 is the share of
FAILED reads whose best hit is the host, compared with PASSED reads.

Usage:
  blast_summary.py --format blast6 --hits failed=failed.blast6 --queries failed=failed.sub.fa ...
"""
import argparse
from collections import Counter

from fyp_lib import write_results


def n_fasta(path):
    with open(path) as fh:
        return sum(1 for line in fh if line.startswith(">"))


def best_hits(path, fmt):
    """Return {query: (score, subject)} keeping the highest-scoring hit per query."""
    best = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if fmt == "blast6":   # qseqid sseqid pident length evalue bitscore
                q, s, score = f[0], f[1], float(f[5])
            else:                  # PAF: qname ... tname(5) ... matches(9)
                q, s, score = f[0], f[5], float(f[9])
            if q not in best or score > best[q][0]:
                best[q] = (score, s)
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--format", choices=["blast6", "paf"], required=True)
    ap.add_argument("--hits", action="append", required=True, help="read_set=hits file")
    ap.add_argument("--queries", action="append", required=True, help="read_set=query fasta")
    ap.add_argument("--results", default="results.tsv")
    ap.add_argument("--table", default="best_hits.tsv")
    a = ap.parse_args()

    queries = dict(s.split("=", 1) for s in a.queries)
    rows = []
    with open(a.table, "w") as tab:
        tab.write("read_set\tread_id\tgroup\tsubject\tscore\n")
        for spec in a.hits:
            read_set, path = spec.split("=", 1)
            n_q = n_fasta(queries[read_set])
            best = best_hits(path, a.format)
            groups = Counter()
            for q, (score, s) in best.items():
                g = s.split("|", 1)[0] if "|" in s else "unlabelled"
                groups[g] += 1
                tab.write(f"{read_set}\t{q}\t{g}\t{s}\t{score}\n")
            rows.append((read_set, "n_queries", n_q))
            rows.append((read_set, "frac_with_hit", round(len(best) / n_q, 5) if n_q else ""))
            for g, c in sorted(groups.items()):
                rows.append((read_set, f"frac_besthit_{g}", round(c / n_q, 5)))
    write_results(a.results, rows)


if __name__ == "__main__":
    main()
