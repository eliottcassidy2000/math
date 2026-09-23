#!/usr/bin/env bash
# collatz_procgen_20260923_theta_run.sh -- reproduce the theta-beyond-phi lane (HYP-9131, HYP-9130, Euler's
# function / signed pentagonal words, the owner's "defect" principle, and the "1089" test).
#
#   bash 04-computation/experiments/collatz_procgen_20260923_theta_run.sh          # full run
#   bash 04-computation/experiments/collatz_procgen_20260923_theta_run.sh --quick  # smoke run
#
# Full run writes 05-knowledge/results/collatz_procgen_20260923_theta.out
# (quick run writes scratch/procgen_theta/theta_quick.out).
# Requirements: python3 with gmpy2, python-flint, sympy, mpmath. One process at a time; small memory.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
EXP="$ROOT/04-computation/experiments"
mkdir -p "$ROOT/scratch/procgen_theta"
if [[ "${1:-}" == "--quick" ]]; then
  OUT="$ROOT/scratch/procgen_theta/theta_quick.out"; FLAG="--quick"
else
  OUT="$ROOT/05-knowledge/results/collatz_procgen_20260923_theta.out"; FLAG=""
fi
{
  echo "# collatz_procgen_20260923_theta -- HYP-9131 (theta beyond phi), Euler's function, HYP-9130 (restricted), defects, '1089'"
  echo "# python: $(python3 -c 'import sys,gmpy2,flint,sympy,mpmath; print(sys.version.split()[0], "gmpy2", gmpy2.version(), "python-flint", flint.__version__, "sympy", sympy.__version__, "mpmath", mpmath.__version__)')"
  python3 "$EXP/collatz_procgen_20260923_theta_hankel.py" $FLAG
  python3 "$EXP/collatz_procgen_20260923_theta_euler.py" $FLAG
  python3 "$EXP/collatz_procgen_20260923_theta_defect.py" $FLAG
  python3 "$EXP/collatz_procgen_20260923_theta_zero.py" $FLAG
  python3 "$EXP/collatz_procgen_20260923_theta_1089.py" $FLAG
} > "$OUT" 2>&1
echo "wrote $OUT"
