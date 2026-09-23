#!/usr/bin/env bash
# Hard-class lane (collatz-procgen-20260922): build the LPF engine and regenerate the .out.
# Usage (from the repository root or the worktree root):
#   bash 04-computation/experiments/collatz_procgen_20260922_hard_run.sh          # full run (about 10-15 min, < 200 MB)
#   bash 04-computation/experiments/collatz_procgen_20260922_hard_run.sh --quick  # smoke run (about 1 min)
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
EXP="$ROOT/04-computation/experiments"
SCR="$ROOT/scratch/procgen_hard"
mkdir -p "$SCR"
cc -O2 -o "$SCR/hard_lpf" "$EXP/collatz_procgen_20260922_hard_lpf.c" -lm
OUT="$ROOT/05-knowledge/results/collatz_procgen_20260922_hard.out"
if [[ "${1:-}" == "--quick" ]]; then OUT="$SCR/hard_quick.out"; fi
python3 "$EXP/collatz_procgen_20260922_hard.py" "$@" > "$OUT"
echo "wrote $OUT"
