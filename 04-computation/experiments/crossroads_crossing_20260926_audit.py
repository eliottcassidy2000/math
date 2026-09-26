"""Replay each lane normally and under -O; retain exact LF-normalized outputs."""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
PREFIX = "crossroads_crossing_20260926_"
LANES = ("colour", "potential", "arithmetic", "arithmetic_run", "dyadic", "resource", "phase_audit", "integration")
MANIFEST = ROOT / "05-knowledge/results" / (PREFIX + "manifest.json")


def require(ok, label):
    if not ok:
        raise RuntimeError(label)


def norm(data):
    return data.decode("utf-8").replace("\r\n", "\n").encode("utf-8")


def digest(data):
    return hashlib.sha256(data).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--refresh", action="store_true")
    args = parser.parse_args()
    lanes = LANES
    hashes = {}
    for lane in lanes:
        relative = Path("04-computation/experiments") / (PREFIX + lane + ".py")
        source = ROOT / relative
        require(source.exists(), "missing lane " + lane)
        runs = []
        for optimized in (False, True):
            command = [sys.executable, "-X", "utf8", "-B"]
            if optimized:
                command.append("-O")
            result = subprocess.run(command + [str(source)], cwd=ROOT, capture_output=True, timeout=180)
            require(result.returncode == 0, lane + " failed: " + result.stderr.decode("utf-8", errors="replace"))
            runs.append(norm(result.stdout))
        require(runs[0] == runs[1], lane + " normal/-O mismatch")
        output_relative = Path("05-knowledge/results") / (PREFIX + lane + ".out")
        output = ROOT / output_relative
        if args.refresh:
            if output.exists():
                require(norm(output.read_bytes()) == runs[0], lane + " disagrees with retained output; inspect before refreshing")
            output.write_bytes(runs[0])
        else:
            require(output.exists() and norm(output.read_bytes()) == runs[0], lane + " retained output mismatch")
        hashes[relative.as_posix()] = digest(norm(source.read_bytes()))
        hashes[output_relative.as_posix()] = digest(runs[0])
        print("PASS " + lane + ": normal = optimized = retained")
    deps = [
        "01-canon/theorems/THM-4476-thin-divergent-orbits-reciprocal-sums-finite.md",
        "01-canon/theorems/THM-4503-collatz-poset-bridge-height-selection.md",
        "05-knowledge/results/collatz_blueprint_20260921_energy.md",
        "05-knowledge/results/duck_zeckendorf_20260925.md",
    ]
    for relative in deps:
        hashes[relative] = digest(norm((ROOT / relative).read_bytes()))
    manifest = {"hash_basis": "UTF-8 with LF line endings", "partial": False, "sha256": hashes}
    if args.refresh:
        MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    else:
        require(json.loads(MANIFEST.read_text(encoding="utf-8")) == manifest, "manifest changed")
    print("PASS " + str(len(hashes)) + " script/output/dependency hash pins")


if __name__ == "__main__":
    main()
