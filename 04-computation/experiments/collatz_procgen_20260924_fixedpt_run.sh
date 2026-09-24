#!/bin/sh
# collatz_procgen_20260924_fixedpt_run.sh -- fixed points / owner's inequality / arXiv:2502.20642v1 audit lane
# (session collatz-procgen-20260922, 2026-09-24).  Runs the three scripts in sequence (one process at a time)
# and writes 05-knowledge/results/collatz_procgen_20260924_fixedpt.out.  About 35 s; peak memory about 420 MB.
set -e
cd "$(dirname "$0")/../.."
OUT=05-knowledge/results/collatz_procgen_20260924_fixedpt.out
{
  echo "collatz_procgen_20260924_fixedpt.out -- output of 04-computation/experiments/collatz_procgen_20260924_fixedpt_{kawasaki,inequality,chains}.py"
  echo "(session collatz-procgen-20260922, lane: fixed points, the owner's inequality, arXiv:2502.20642v1)"
  echo
  echo "#################### collatz_procgen_20260924_fixedpt_kawasaki.py ####################"
  python3 04-computation/experiments/collatz_procgen_20260924_fixedpt_kawasaki.py
  echo
  echo "#################### collatz_procgen_20260924_fixedpt_inequality.py ####################"
  python3 04-computation/experiments/collatz_procgen_20260924_fixedpt_inequality.py
  echo
  echo "#################### collatz_procgen_20260924_fixedpt_chains.py ####################"
  python3 04-computation/experiments/collatz_procgen_20260924_fixedpt_chains.py
} > "$OUT"
echo "wrote $OUT"
