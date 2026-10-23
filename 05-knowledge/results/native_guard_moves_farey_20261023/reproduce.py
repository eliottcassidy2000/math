"""Replay exact native guards, moves, Fourier and Farey receipts."""
from hashlib import sha256
from pathlib import Path
import json
import os
import re
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def main():
    provenance = json.loads((HERE / "SOURCE_PROVENANCE.json").read_text(encoding="utf-8"))
    for name, expected in provenance["inherited_sha256"].items():
        if digest(ROOT / name) != expected:
            raise RuntimeError(f"Inherited input changed; review before replay: {name}")
    mode = ["-O"] if "--optimized" in sys.argv else []
    env = dict(os.environ, PYTHONHASHSEED="0", PYTHONIOENCODING="utf-8")
    jobs = ["guards", "moves", "fourier", "farey", "native_moves"]
    for lane in jobs:
        before = digest(HERE / lane / "receipt.json")
        run = subprocess.run([sys.executable, *mode, "-X", "utf8", "-B", str(HERE / lane / "check.py")],
                             cwd=ROOT, env=env, text=True, encoding="utf-8",
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if run.returncode:
            print(run.stdout)
            raise RuntimeError(f"Checker failed: {lane}")
        if before != digest(HERE / lane / "receipt.json"):
            raise RuntimeError(f"Replay mutated frozen receipt: {lane}")
        print(f"PASS {lane}: frozen receipt unchanged", flush=True)
    audits = ["guards/audit_coupled_controls.py", "farey/AUDIT_NATIVE_MOVES_CHECK.py"]
    for program in audits:
        run = subprocess.run([sys.executable, *mode, "-X", "utf8", "-B", str(HERE / program)],
                             cwd=ROOT, env=env, text=True, encoding="utf-8",
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if run.returncode:
            print(run.stdout)
            raise RuntimeError(f"Independent audit failed: {program}")
        print(f"PASS independent {program}", flush=True)
    # This packet contains relative Markdown file links only; URLs and anchors
    # are excluded from the local existence check.
    count = 0
    for page in HERE.rglob("*.md"):
        for target in re.findall(r"\]\(([^)]+)\)", page.read_text(encoding="utf-8")):
            if target.startswith(("http:", "https:", "#")):
                continue
            target = target.split("#")[0]
            if target and not (page.parent / target).exists():
                raise RuntimeError(f"Missing local link in {page.name}: {target}")
            count += 1
    print(f"PASS {len(jobs)} producers + {len(audits)} independent checks, {count} local links, inherited source hashes; mode={mode or ['normal']}")


if __name__ == "__main__":
    main()
