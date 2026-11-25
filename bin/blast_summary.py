#!/usr/bin/env python3
"""Summarise read-search hits (BLAST-style 12-column tabular, or minimap2 PAF).

Tabular format comes from blastn (megablast, dc-megablast, blastn), DIAMOND
blastx, LAST (BlastTab) and MMseqs2 (m8); PAF from minimap2.

Each read's best hit is put in a group. Labelled databases (build_blast_db.py
writes 'group|original_id') give the group directly. Unlabelled databases such
as NCBI core_nt give a taxid (13th column, staxids); with --taxonomy the
taxid is placed in host (under --host-taxid, default 3118 = Ulva),
other_green_algae (Viridiplantae), other_eukaryota, bacteria, archaea, viruses. The key number for RQ1 is the share of
FAILED reads whose best hit is the host, compared with PASSED reads.

Usage:
  blast_summary.py --format blast6 --hits failed=failed.blast6 --queries failed=failed.sub.fa ...
"""
import argparse
import os
from collections import Counter

from fyp_lib import write_results

GROUP_TAXA = [("bacteria", "2"), ("archaea", "2157"), ("viruses", "10239"),
              ("other_green_algae", "33090"), ("other_eukaryota", "2759")]


def load_parents(tax_dir):
    parent = {}
    with open(os.path.join(tax_dir, "nodes.dmp")) as fh:
        for line in fh:
            f = line.split("\t|\t", 2)
            parent[f[0]] = f[1]
    return parent


def taxid_group(taxid, parent, host):
    lineage, t = set(), taxid
    while t in parent and t not in lineage:
        lineage.add(t)
        if t == "1":
            break
        t = parent[t]
    if host in lineage:
        return "host"
    for name, tid in GROUP_TAXA:
        if tid in lineage:
            return name
    return "other"


def n_fasta(path):
    with open(path) as fh:
        return sum(1 for line in fh if line.startswith(">"))


def best_hits(path, fmt):
    """Return {query: (score, subject)} keeping the highest-scoring hit per query."""
    best = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if fmt == "blast6":   # standard 12-column tabular (blastn, diamond, last BlastTab, mmseqs m8) [+ staxids]
                q, s, score = f[0], f[1], float(f[11] if len(f) >= 12 else f[5])
                if len(f) >= 13 and "|" not in s and f[12] not in ("", "0", "N/A"):
                    s = "taxid:" + f[12].split(";")[0]
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
    ap.add_argument("--taxonomy", help="NCBI taxdump dir, to group taxid hits (unlabelled databases)")
    ap.add_argument("--host-taxid", default="3118", help="host clade taxid (default 3118 = genus Ulva)")
    a = ap.parse_args()
    parent = load_parents(a.taxonomy) if a.taxonomy else {}

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
                if "|" in s:
                    g = s.split("|", 1)[0]
                elif s.startswith("taxid:") and parent:
                    g = taxid_group(s[6:], parent, a.host_taxid)
                else:
                    g = "unlabelled"
                groups[g] += 1
                tab.write(f"{read_set}\t{q}\t{g}\t{s}\t{score}\n")
            rows.append((read_set, "n_queries", n_q))
            rows.append((read_set, "frac_with_hit", round(len(best) / n_q, 5) if n_q else ""))
            for g, c in sorted(groups.items()):
                rows.append((read_set, f"frac_besthit_{g}", round(c / n_q, 5)))
    write_results(a.results, rows)


if __name__ == "__main__":
    main()
