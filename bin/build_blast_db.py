#!/usr/bin/env python3
"""Build a labelled BLAST database from groups of FASTA files.

Each sequence header becomes '>group|original_id', so blast_summary.py can tell
whether a read's best hit is host, other algae or bacteria. The combined FASTA
is kept next to the database because minimap2 (the fast alternative tool)
aligns to the FASTA directly.

Usage:
  build_blast_db.py --out /path/db/algae_bact \
      --group host=ulva/*.fna.gz --group bacteria=bact/*.fna.gz
Then point resources.blast_dbs.<name> in config/experiment.yaml at:
  db: /path/db/algae_bact      fasta: /path/db/algae_bact.fa
"""
import argparse
import glob
import hashlib
import json
import os
import subprocess

from fyp_lib import open_text


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", required=True, help="database prefix (no extension)")
    ap.add_argument("--group", action="append", required=True, help="label=glob (repeatable)")
    a = ap.parse_args()

    fasta = a.out + ".fa"
    os.makedirs(os.path.dirname(os.path.abspath(fasta)), exist_ok=True)
    manifest = {}
    with open(fasta, "w") as out:
        for spec in a.group:
            label, pattern = spec.split("=", 1)
            files = sorted(glob.glob(pattern))
            if not files:
                raise SystemExit(f"No files match {pattern!r} for group {label}")
            manifest[label] = []
            for f in files:
                n = 0
                with open_text(f) as fh:
                    for line in fh:
                        if line.startswith(">"):
                            out.write(f">{label}|{line[1:].split()[0][:40]}\n")
                            n += 1
                        else:
                            out.write(line)
                manifest[label].append({"file": os.path.abspath(f), "n_seqs": n})
    subprocess.run(["makeblastdb", "-in", fasta, "-dbtype", "nucl", "-out", a.out], check=True)
    h = hashlib.md5(open(fasta, "rb").read()).hexdigest()
    with open(a.out + ".manifest.json", "w") as fh:
        json.dump({"fasta": fasta, "fasta_md5": h, "groups": manifest}, fh, indent=2)
    print(f"Built {a.out} ({sum(len(v) for v in manifest.values())} files)")


if __name__ == "__main__":
    main()
