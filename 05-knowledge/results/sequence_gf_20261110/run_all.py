#!/usr/bin/env python3
"""Replay the four new sequence compilers; the heat script refreshes output.json."""
from pathlib import Path
import subprocess
import sys

root = Path(__file__).resolve().parent
scripts = [
    "burnside/cluster_compiler.py",
    "grammar/compiler.py",
    "heat/compile_heat_covariance.py",
    "rays/compile_rays.py",
]
for script in scripts:
    run = subprocess.run([sys.executable, str(root / script)], cwd=root,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if run.returncode:
        print(run.stdout)
        print(run.stderr, file=sys.stderr)
        raise SystemExit(f"FAIL: {script} ({run.returncode})")
    print(f"PASS: {script} ({len(run.stdout.encode())} output bytes)")
print("All four new-formula compilers passed. Prior tournament claims are inherited; Chan-Pak OPEN.")
