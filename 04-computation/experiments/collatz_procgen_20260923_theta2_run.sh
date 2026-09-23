#!/usr/bin/env bash
# collatz_procgen_20260923_theta2_run.sh -- reproduce the theta round-2 lane: HYP-9132 (2-adic
# non-quadraticity of theta at rho = 2^10/3^9, Theorem NQ) and HYP-9127 (cube determinant families, no-go).
#
#   bash 04-computation/experiments/collatz_procgen_20260923_theta2_run.sh          # full run (~2 min)
#   bash 04-computation/experiments/collatz_procgen_20260923_theta2_run.sh --quick  # smoke run (~15 s)
#
# Full run writes 05-knowledge/results/collatz_procgen_20260923_theta2.out
# (quick run writes scratch/procgen_theta2/theta2_quick.out).
# Requirements: python3 with gmpy2, python-flint, sympy, mpmath, numpy. One process at a time; < 200 MB.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
EXP="$ROOT/04-computation/experiments"
mkdir -p "$ROOT/scratch/procgen_theta2"
if [[ "${1:-}" == "--quick" ]]; then
  OUT="$ROOT/scratch/procgen_theta2/theta2_quick.out"; FLAG="--quick"
else
  OUT="$ROOT/05-knowledge/results/collatz_procgen_20260923_theta2.out"; FLAG=""
fi
{
  echo "# collatz_procgen_20260923_theta2 -- HYP-9132 (non-quadraticity, Theorem NQ) and HYP-9127 (cube determinant no-go)"
  echo "# python: $(python3 -c 'import sys,gmpy2,flint,sympy,mpmath,numpy; print(sys.version.split()[0], "gmpy2", gmpy2.version(), "python-flint", flint.__version__, "sympy", sympy.__version__, "mpmath", mpmath.__version__, "numpy", numpy.__version__)')"
  python3 "$EXP/collatz_procgen_20260923_theta2_nq.py" $FLAG
  python3 "$EXP/collatz_procgen_20260923_theta2_cube.py" $FLAG
} > "$OUT" 2>&1
echo "wrote $OUT"
