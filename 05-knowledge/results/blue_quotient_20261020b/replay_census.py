"""Replay compact n9 certificates; --rebuild also regenerates their universe.

Requires g++ on PATH, Python 3, and NetworkX for the final Python analysis.
Executables and expanded TSVs stay in an automatically cleaned temp folder.
No Git operation and no installed package changes are performed.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def normalized(raw):
    return raw.replace(b"\r\n", b"\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rebuild", action="store_true")
    args = parser.parse_args()
    compiler = shutil.which("g++")
    need(compiler is not None, "g++ is required")
    manifest = json.loads((HERE / "census_manifest.json").read_text(encoding="utf-8"))
    for kind, base in (("owned_sha256", HERE), ("inherited_sha256", ROOT)):
        for name, expected in manifest[kind].items():
            raw = (base / name).read_bytes()
            if not name.endswith(".gz"):
                raw = normalized(raw)
            need(hashlib.sha256(raw).hexdigest() == expected, "input hash mismatch: " + name)
    with tempfile.TemporaryDirectory(prefix="blue-quotient-n9-replay-") as scratch:
        scratch = Path(scratch)
        prefix = scratch / "n9"
        archived = {}
        for suffix in ("classes.tsv", "blue.tsv", "black_tree.tsv"):
            data = gzip.decompress((HERE / ("n9_" + suffix + ".gz")).read_bytes())
            archived[suffix] = data
            if not args.rebuild:
                (scratch / ("n9_" + suffix)).write_bytes(data)
        if args.rebuild:
            binary = scratch / "producer.exe"
            subprocess.run([compiler, "-O3", "-std=c++17", str(HERE / "census.cpp"), "-o", str(binary)], check=True)
            subprocess.run([str(binary),
                            str(ROOT / "04-computation/experiments/tiling_line_metagraph_census_20261020_classes_n8.tsv"),
                            str(prefix)], check=True)
            for suffix, data in archived.items():
                need(normalized((scratch / ("n9_" + suffix)).read_bytes()) == normalized(data),
                     "rebuilt certificate mismatch: " + suffix)
            need((scratch / "n9_summary.txt").read_text() == (HERE / "n9_summary.txt").read_text(),
                 "rebuilt summary mismatch")
        binary = scratch / "independent.exe"
        subprocess.run([compiler, "-O3", "-std=c++17", str(HERE / "audit.cpp"), "-o", str(binary)], check=True)
        subprocess.run([str(binary), str(prefix)], check=True)
        need((scratch / "n9_independent_audit.txt").read_text() ==
             (HERE / "n9_independent_audit.txt").read_text(), "independent report mismatch")
    outputs = []
    for flags in ([], ["-O"]):
        result = subprocess.run([sys.executable, *flags, str(HERE / "analyse_census.py")],
                                check=True, capture_output=True, text=True)
        outputs.append(result.stdout)
    need(outputs[0] == outputs[1], "Python optimized analysis differs")
    print("PASS: input hashes, independent full-universe/blue/tree replay, and Python normal/-O analysis.")
    if args.rebuild:
        print("PASS: all producer certificates regenerated identically after newline normalization.")


if __name__ == "__main__":
    main()
