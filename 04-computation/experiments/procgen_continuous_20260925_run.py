#!/usr/bin/env python3
"""procgen_continuous_20260925_run.py

Runner for the lane "natural boundaries, Mahler/Cobham and harmonic trees"
(session collatz-procgen-20260922, 2026-09-25).  Runs the four parts one after another
(one process at a time; each stays below 700 MB), and writes
05-knowledge/results/procgen_continuous_20260925.out.  Exits nonzero if any part fails.

  python3 04-computation/experiments/procgen_continuous_20260925_run.py
"""
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
OUT = os.path.join(ROOT, "05-knowledge", "results", "procgen_continuous_20260925.out")
PARTS = [
    ("A", "procgen_continuous_20260925_natural_boundary.py"),
    ("B", "procgen_continuous_20260925_mahler.py"),
    ("C", "procgen_continuous_20260925_tree_harmonic.py"),
    ("D", "procgen_continuous_20260925_hex_game.py"),
]


def main():
    chunks = []
    status = []
    for tag, fn in PARTS:
        t = time.time()
        p = subprocess.run([sys.executable, os.path.join(HERE, fn)], capture_output=True, text=True)
        dt = time.time() - t
        memline = [ln for ln in p.stderr.splitlines() if ln.startswith("[mem] end")]
        status.append((tag, fn, p.returncode, dt, memline[-1] if memline else "[mem] n/a"))
        chunks.append("#" * 78 + f"\n# PART {tag}: {fn}\n" + "#" * 78 + "\n" + p.stdout)
        if p.returncode != 0:
            chunks.append("STDERR:\n" + p.stderr)
    head = ["procgen_continuous_20260925.out -- lane: natural boundaries, Mahler/Cobham, harmonic trees, Hex",
            "session collatz-procgen-20260922, 2026-09-25; produced by procgen_continuous_20260925_run.py",
            ""]
    for tag, fn, rc, dt, mem in status:
        head.append(f"  part {tag}: {fn:<52s} exit {rc}  ({dt:.1f}s)  {mem}")
    head.append("")
    with open(OUT, "w") as f:
        f.write("\n".join(head) + "\n" + "\n".join(chunks))
    print("\n".join(head))
    print(f"wrote {OUT}")
    if any(rc != 0 for _, _, rc, _, _ in status):
        sys.exit(1)


if __name__ == "__main__":
    main()
