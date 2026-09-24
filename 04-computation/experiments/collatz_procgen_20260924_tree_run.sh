#!/bin/sh
# Reproduce 05-knowledge/results/collatz_procgen_20260924_tree.out (lane: inverse tree mod 192).
# Run from the repository root.  Total about 25 s; peak memory about 400 MB; one process at a time.
set -e
OUT=05-knowledge/results/collatz_procgen_20260924_tree.out
{
  echo "collatz_procgen_20260924_tree.out -- inverse tree mod 192 (session collatz-procgen-20260922, 2026-09-24)"
  echo "Produced by 04-computation/experiments/collatz_procgen_20260924_tree_run.sh"
  echo "Sections A1-A8: ..._tree_automaton.py;  C1-C4: ..._tree_controls.py;  K1-K7: ..._tree_counts.py"
  python3 04-computation/experiments/collatz_procgen_20260924_tree_automaton.py
  python3 04-computation/experiments/collatz_procgen_20260924_tree_controls.py
  python3 04-computation/experiments/collatz_procgen_20260924_tree_counts.py
} > "$OUT"
echo "wrote $OUT"
