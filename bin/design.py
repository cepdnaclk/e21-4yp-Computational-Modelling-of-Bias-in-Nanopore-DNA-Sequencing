#!/usr/bin/env python3
"""Turn config/experiment.yaml into design.csv: one row per combination (leaf).

The design decides WHICH combinations exist; Nextflow only executes them.
Each row carries its full folder path under outputs/ plus the paths of its
parent stages, so every upstream stage is computed once and shared:

  basecaller=<b>/model=<m>                                    stage: basecall
    trim=<t>                                                  stage: trim
      minq=<q>/minlen=<l>                                     stage: qc
        analysis=host_map/aligner=<a>/host_ref=<r>            leaf
        analysis=classify/classifier=<c>/db=<d>               leaf
        analysis=blast/blast_tool=<t>/blast_db=<d>            leaf

Modes
  full : every combination of every factor level (the default; affordable for one run)
  ofat : the baseline, plus each factor changed one at a time with the others at baseline

Usage: design.py config/experiment.yaml --out design.csv [--mode ofat] [--strict]
"""
import argparse
import csv
import hashlib
import itertools
import os
import sys

import yaml

UPSTREAM = ["basecall", "trim", "minq", "minlen"]
ANALYSIS_FACTORS = {
    "host_map": ["aligner", "host_ref"],
    "classify": ["db"],
    "blast": ["blast_tool", "blast_db"],
}
COLUMNS = [
    "combo_id", "path", "qc_path", "trim_path", "bc_path",
    "basecaller", "model", "bc_type", "model_spec", "import_fastq",
    "trim", "minq", "minlen", "analysis",
    "aligner", "host_ref", "host_ref_fasta",
    "classifier", "db", "db_path",
    "blast_tool", "blast_db", "blast_db_prefix", "blast_db_fasta",
]


def warn(msg):
    print(f"[design] {msg}", file=sys.stderr)


def exists_any(path):
    return bool(path) and (os.path.exists(path) or any(
        os.path.exists(path + ext) for ext in (".nal", ".nin", ".ndb", ".fa")))


def fmt_num(v):
    return str(int(v)) if float(v).is_integer() else str(v)


class Design:
    def __init__(self, cfg, strict):
        self.cfg, self.strict = cfg, strict
        self.res = cfg.get("resources", {})
        self.dropped = set()

    def check(self, ok, what):
        if ok:
            return True
        if self.strict:
            raise SystemExit(f"[design] missing resource: {what}")
        if what not in self.dropped:
            warn(f"skipping (resource not found yet): {what}")
            self.dropped.add(what)
        return False

    # ---- factor levels -------------------------------------------------
    def levels(self):
        f = self.cfg["factors"]
        lv = {
            "basecall": list(f["basecall"]),
            "trim": list(f["trim"]),
            "minq": [fmt_num(q) for q in f["minq"]],
            "minlen": [fmt_num(x) for x in f["minlen"]],
        }
        for analysis, keys in ANALYSIS_FACTORS.items():
            spec = (self.cfg.get("analyses") or {}).get(analysis)
            if spec:
                for k in keys:
                    lv[f"{analysis}.{k}"] = list(spec[k])
        return lv

    def baseline(self, lv):
        b = dict((self.cfg.get("design") or {}).get("baseline") or {})
        out = {}
        for k, levels in lv.items():
            short = k.split(".")[-1]
            v = b.get(k, b.get(short, levels[0]))
            v = fmt_num(v) if k in ("minq", "minlen") else v
            if v not in levels:
                raise SystemExit(f"[design] baseline {k}={v} is not one of the listed levels {levels}")
            out[k] = v
        return out

    # ---- one row -------------------------------------------------------
    def row(self, up, analysis, an):
        basecaller, model = up["basecall"].split(":", 1)
        bc = self.res["basecallers"].get(basecaller)
        if bc is None or model not in bc["models"]:
            raise SystemExit(f"[design] no resources.basecallers entry for {basecaller}:{model}")
        mspec = bc["models"][model]
        r = dict.fromkeys(COLUMNS, "")
        r.update(basecaller=basecaller, model=model, bc_type=bc["type"],
                 trim=up["trim"], minq=up["minq"], minlen=up["minlen"], analysis=analysis)
        if bc["type"] == "import":
            dirs = mspec["fastq"]
            if not self.check(all(os.path.isdir(d) for d in dirs), f"import fastq dirs for {basecaller}:{model}"):
                return None
            r["import_fastq"] = ";".join(dirs)
        else:
            r["model_spec"] = mspec
            if bc["type"] == "dorado":
                mdir = os.path.join(os.path.expanduser(self.cfg["dorado_models"]), mspec)
                if not self.check(os.path.isdir(mdir), f"dorado model {mdir}"):
                    return None
        r["bc_path"] = f"basecaller={basecaller}/model={model}"
        r["trim_path"] = f"{r['bc_path']}/trim={r['trim']}"
        r["qc_path"] = f"{r['trim_path']}/minq={r['minq']}/minlen={r['minlen']}"

        if analysis == "host_map":
            ref = self.res["host_refs"].get(an["host_ref"], {})
            if not self.check(os.path.exists(ref.get("fasta", "")), f"host_ref {an['host_ref']}"):
                return None
            r.update(aligner=an["aligner"], host_ref=an["host_ref"], host_ref_fasta=ref["fasta"])
            leaf = f"analysis=host_map/aligner={an['aligner']}/host_ref={an['host_ref']}"
        elif analysis == "classify":
            db = self.res["classify_dbs"].get(an["db"], {})
            if not self.check(os.path.exists(db.get("path", "")), f"classify db {an['db']}"):
                return None
            r.update(classifier=db["classifier"], db=an["db"], db_path=db["path"])
            leaf = f"analysis=classify/classifier={db['classifier']}/db={an['db']}"
        else:
            db = self.res["blast_dbs"].get(an["blast_db"], {})
            need = db.get("fasta", "") if an["blast_tool"] == "minimap2" else db.get("db", "")
            if not self.check(exists_any(need), f"blast_db {an['blast_db']} for {an['blast_tool']}"):
                return None
            r.update(blast_tool=an["blast_tool"], blast_db=an["blast_db"],
                     blast_db_prefix=db.get("db", ""), blast_db_fasta=db.get("fasta", ""))
            leaf = f"analysis=blast/blast_tool={an['blast_tool']}/blast_db={an['blast_db']}"
        r["path"] = f"{r['qc_path']}/{leaf}"
        r["combo_id"] = hashlib.sha1(r["path"].encode()).hexdigest()[:12]
        return r

    # ---- grids -------------------------------------------------------------
    def analysis_grid(self, lv, analysis, fixed=None):
        keys = [f"{analysis}.{k}" for k in ANALYSIS_FACTORS[analysis] if f"{analysis}.{k}" in lv]
        if not keys:
            return []
        choices = [[fixed[k]] if fixed else lv[k] for k in keys]
        return [{k.split(".")[1]: v for k, v in zip(keys, combo)} for combo in itertools.product(*choices)]

    def build(self, mode):
        lv, rows = self.levels(), []
        analyses = [a for a in ANALYSIS_FACTORS if any(k.startswith(a + ".") for k in lv)]
        if mode == "full":
            for up in itertools.product(*(lv[k] for k in UPSTREAM)):
                upd = dict(zip(UPSTREAM, up))
                for analysis in analyses:
                    for an in self.analysis_grid(lv, analysis):
                        rows.append(self.row(upd, analysis, an))
        else:
            base = self.baseline(lv)
            up_base = {k: base[k] for k in UPSTREAM}
            variants = [dict(up_base)]
            for k in UPSTREAM:
                variants += [dict(up_base, **{k: v}) for v in lv[k] if v != base[k]]
            for upd in variants:   # upstream changes: every analysis at its baseline
                for analysis in analyses:
                    for an in self.analysis_grid(lv, analysis, fixed=base):
                        rows.append(self.row(upd, analysis, an))
            for analysis in analyses:   # analysis changes: upstream at baseline
                for key in [k for k in lv if k.startswith(analysis + ".")]:
                    for v in lv[key]:
                        if v != base[key]:
                            fixed = dict(base, **{key: v})
                            for an in self.analysis_grid(lv, analysis, fixed=fixed):
                                rows.append(self.row(up_base, analysis, an))
        uniq = {}
        for r in rows:
            if r:
                uniq[r["path"]] = r
        return list(uniq.values())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("config")
    ap.add_argument("--out", required=True)
    ap.add_argument("--mode", choices=["full", "ofat"], help="overrides design.mode in the config")
    ap.add_argument("--strict", action="store_true", help="fail instead of skipping missing resources")
    a = ap.parse_args()

    cfg = yaml.safe_load(open(a.config))
    mode = a.mode or (cfg.get("design") or {}).get("mode", "full")
    rows = Design(cfg, a.strict).build(mode)
    if not rows:
        raise SystemExit("[design] no runnable combinations (check resources)")
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS)
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: r["path"]))
    count = lambda col: len({r[col] for r in rows})
    warn(f"mode={mode}: {len(rows)} leaf combinations | basecall nodes {count('bc_path')}, "
         f"trim nodes {count('trim_path')}, qc nodes {count('qc_path')}")


if __name__ == "__main__":
    main()
