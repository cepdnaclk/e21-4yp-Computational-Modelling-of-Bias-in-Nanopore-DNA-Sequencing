#!/usr/bin/env python3
"""Write metadata.json for one node of the outputs/ tree.

The folder path is the single source of truth for the factor levels: every
'key=value' folder name becomes a factor. The metadata adds what the path
cannot hold: the exact command, tool versions, fixed parameters, the parent
node, and checksums of every output file, so any result can be traced back
and reproduced.
"""
import argparse
import datetime
import hashlib
import json
import os
import platform
import socket

from fyp_lib import MEANQ_DEFINITION, parse_combo_path


def md5(path, block=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(block), b""):
            h.update(chunk)
    return h.hexdigest()


def kv_list(items):
    out = {}
    for it in items or []:
        if "=" not in it:
            raise SystemExit(f"Expected key=value, got {it!r}")
        k, v = it.split("=", 1)
        out[k.strip()] = v.strip()
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--stage", required=True, help="basecall | trim | qc | host_map | classify | blast")
    ap.add_argument("--path", required=True, help="combination path relative to outputs/")
    ap.add_argument("--parent", default="", help="parent combination path ('' for the root stage)")
    ap.add_argument("--command-file", help="file holding the exact shell command(s) run")
    ap.add_argument("--tool", action="append", help="name=version (repeatable)")
    ap.add_argument("--param", action="append", help="fixed parameter name=value (repeatable)")
    ap.add_argument("--inputs", nargs="*", default=[], help="input files or descriptions")
    ap.add_argument("--outputs", nargs="*", default=[], help="output files to checksum")
    ap.add_argument("--out", default="metadata.json")
    a = ap.parse_args()

    path = a.path.strip("/")
    meta = {
        "combo_id": hashlib.sha1(path.encode()).hexdigest()[:12],
        "path": path,
        "stage": a.stage,
        "parent_path": a.parent.strip("/"),
        "factors": parse_combo_path(path),
        "params": kv_list(a.param),
        "tools": kv_list(a.tool),
        "command": open(a.command_file).read().strip() if a.command_file else "",
        "inputs": [os.path.basename(p) if os.path.exists(p) else p for p in a.inputs],
        "outputs": {},
        "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "host": socket.gethostname(),
        "python": platform.python_version(),
    }
    if a.stage in ("basecall", "trim", "qc"):
        meta["params"].setdefault("meanq_definition", MEANQ_DEFINITION)
    for p in a.outputs:
        if os.path.exists(p) and os.path.isfile(p):
            meta["outputs"][os.path.basename(p)] = {"bytes": os.path.getsize(p), "md5": md5(os.path.realpath(p))}
    with open(a.out, "w") as fh:
        json.dump(meta, fh, indent=2)
        fh.write("\n")


if __name__ == "__main__":
    main()
