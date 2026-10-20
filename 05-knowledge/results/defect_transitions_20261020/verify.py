"""Verify the frozen research packet; optionally replay its exact experiments."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SCRIPTS = ("budget_check.py", "strata_check.py", "blue_check.py", "petersen_check.py")
MANIFEST = HERE / "manifest.json"


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def normalized(data):
    return data.replace(b"\r\n", b"\n")


def digest(path):
    return hashlib.sha256(normalized(path.read_bytes())).hexdigest()


def inputs():
    own = [p for p in HERE.iterdir()
           if p.is_file() and p.suffix in (".py", ".md", ".out")
           and p.name != "verify.out" and ".optimized." not in p.name]
    deps = [
        ROOT / "05-knowledge/results/defect_boundary_20261020/sparse_count.py",
        ROOT / "05-knowledge/results/metagraph_growth_20261020/enumeration.py",
        ROOT / "05-knowledge/results/metagraph_growth_20261020/growth_check.py",
    ]
    for n in range(3, 8):
        for kind in ("classes", "lines"):
            deps.append(ROOT / (
                f"04-computation/experiments/tiling_line_metagraph_census_20261020_{kind}_n{n}.tsv"))
    return sorted(set(own + deps))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", action="store_true", help="record current files after audit")
    parser.add_argument("--replay", action="store_true", help="rerun normal and optimized experiments")
    args = parser.parse_args()
    if args.record:
        record = {"normalization": "CRLF to LF before SHA256",
                  "files": {p.relative_to(ROOT).as_posix(): digest(p) for p in inputs()}}
        MANIFEST.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n",
                            encoding="utf-8")
    record = json.loads(MANIFEST.read_text(encoding="utf-8"))
    for name, expected in record["files"].items():
        path = ROOT / name
        need(path.is_file() and digest(path) == expected, f"hash mismatch: {name}")
    print("Normalized SHA256 inputs verified:", len(record["files"]))
    for name in SCRIPTS:
        tree = ast.parse((HERE / name).read_text(encoding="utf-8"))
        need(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
             f"removable assertion in {name}")
        need((HERE / name.replace(".py", ".out")).is_file(), f"missing transcript: {name}")
    print("All four experiment sources parse; checks contain no removable assertions.")
    if args.replay:
        for name in SCRIPTS:
            expected = normalized((HERE / name.replace(".py", ".out")).read_bytes())
            for flags in ([], ["-O"]):
                result = subprocess.run([sys.executable, *flags, str(HERE / name)],
                                        cwd=ROOT, capture_output=True, check=False)
                need(result.returncode == 0,
                     f"{name} {flags} failed: {result.stderr.decode(errors='replace')}")
                need(normalized(result.stdout) == expected, f"transcript mismatch: {name} {flags}")
            print("Normal and optimized transcript equality:", name, flush=True)
    print("PASS")


if __name__ == "__main__":
    main()
