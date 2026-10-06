#!/usr/bin/env python3
"""
chessboard_weave_20261006_hp_blocking_theory_efg.py

Runs sections E, F, G of chessboard_weave_20261006_hp_blocking_theory.py (the first full run of
that script was stopped inside section E by a bug: k-path with k >= N has no blocking set, so
the deepening loop never ended; fixed by skipping k >= N).  Sections A-D are in
chessboard_weave_20261006_hp_blocking_theory.out.
"""
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("theory", os.path.join(HERE, "chessboard_weave_20261006_hp_blocking_theory.py"))
theory = importlib.util.module_from_spec(spec)
spec.loader.exec_module(theory)

if __name__ == "__main__":
    print(__doc__.strip())
    print()
    sys.stdout.flush()
    theory.main("EFG")
