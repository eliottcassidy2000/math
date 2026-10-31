#!/usr/bin/env python3
"""Replay attributed Dauth certificates, or correctly reprice fixed witnesses.

External checkout: https://github.com/Lachy-Dauth/ramsey-neighbourhood-counting
Pinned revision below. This wrapper is not an independent checker rewrite.
Run normal Python, never -O: the external checker validates using assertions.
"""
import argparse
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile

REV = "5770baac746878a3a63b9720fdda97726fe02851"
FILES = {
    "mixed_checker": "ancillary/flag5_r34k/check/check_flag5_mix.py",
    "344": "ancillary/flag5_r34k/certificates/R344_73.json",
    "345": "ancillary/flag5_r34k/certificates/R345_137.json",
    "336_checker": "ancillary/flag_r33k/check_r33k5.py",
    "336": "ancillary/flag_r33k/certificates/five/k6/excl_86.json",
    "444_checker5": "ancillary/flag5_r444/checkers/check5.py",
    "444_checker6": "ancillary/flag5_r444/checkers/check6.py",
    "444": "ancillary/flag5_r444/certificates/cert5x_198.txt",
}


def load(path):
    spec = importlib.util.spec_from_file_location("dauth_check5", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def reprice444(root):
    """Keep all square vectors and multipliers; recompute every valid cap."""
    mod = load(root / FILES["444_checker5"])
    mod.R344 = 73
    mod.DMAX = 72
    d = copy.deepcopy(mod.load_text(root / FILES["444"]))
    d["n"] = 197
    caps = []
    for kind, p4, kap, _, mu in d["caps"]:
        edges = dict(zip(((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)), p4))
        edge = lambda a, b, edges=edges: edges[tuple(sorted((a, b)))]
        known = mod.known_of(kap, edge)
        q = (mod.cap_of(kap, edge) if kind == "up" else 196 - 2 * mod.DMAX) - known
        caps.append((kind, p4, kap, q, mu))
    d["caps"] = caps
    squares = []
    for sq in d["squares"]:
        if sq[0] != "eloc":
            squares.append(sq)
            continue
        kind, b, c, _ = sq[3]
        if kind == "up":
            known = int(b is None and c == 0) + int(c is None and b == 0)
            q = mod.cap_of((b, c, None, None), lambda a, b: 0) - known
        else:
            q = 196 - 2 * mod.DMAX - int(b == 0)
        squares.append(sq[:3] + ((kind, b, c, q),))
    d["squares"] = squares
    return mod.check(d, "fixed 198 witness repriced at197 with R344<=73")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("checkout", type=Path)
    ap.add_argument("mode", choices=["manifest", "344", "345", "336", "444-5", "444-6", "72", "136", "197"])
    a = ap.parse_args()
    root = a.checkout.resolve()
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()
    if revision != REV or not __debug__:
        raise SystemExit("Need pinned revision and Python assertions enabled")
    if a.mode == "manifest":
        print(json.dumps({"revision": revision, "files": {
            key: {"path": path, "sha256": hashlib.sha256((root / path).read_bytes()).hexdigest()}
            for key, path in FILES.items()}}, indent=2))
        return
    if a.mode == "197":
        # Rejection is the observed control result, not a lower-bound proof.
        raise SystemExit(0 if reprice444(root) else 1)
    if a.mode in ("72", "136"):
        key = "344" if a.mode == "72" else "345"
        d = json.loads((root / FILES[key]).read_text())
        d["n"] = int(a.mode)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "repriced.json"
            path.write_text(json.dumps(d))
            cmd = [sys.executable, "-B", str(root / FILES["mixed_checker"]), str(path)]
            if a.mode == "136":
                cmd += ["--input", "3,4,4:73"]
            raise SystemExit(subprocess.call(cmd, cwd=root))
    key = a.mode.split("-")[0]
    checker = "mixed_checker" if key in ("344", "345") else (
        "336_checker" if key == "336" else "444_checker" + a.mode[-1])
    raise SystemExit(subprocess.call([sys.executable, "-B", str(root / FILES[checker]), str(root / FILES[key])], cwd=root))


if __name__ == "__main__":
    main()
