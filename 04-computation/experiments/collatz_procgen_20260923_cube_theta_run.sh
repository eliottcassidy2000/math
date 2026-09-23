#!/usr/bin/env bash
# collatz_procgen_20260923_cube_theta_run.sh -- reproduce the cube-theta lane (HYP-9127).
#
#   bash 04-computation/experiments/collatz_procgen_20260923_cube_theta_run.sh          # full run
#   bash 04-computation/experiments/collatz_procgen_20260923_cube_theta_run.sh --quick  # smoke run
#
# Full run writes 05-knowledge/results/collatz_procgen_20260923_cube_theta.out
# (quick run writes scratch/procgen_cube/cube_theta_quick.out).
# Requirements: python3 with gmpy2 and python-flint; PARI/GP (gp) optional for the 10^6/10^7-bit
# reconstruction certificates. One process at a time; peak RSS stays well below 900 MB.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
EXP="$ROOT/04-computation/experiments"
mkdir -p "$ROOT/scratch/procgen_cube"
if [[ "${1:-}" == "--quick" ]]; then
  OUT="$ROOT/scratch/procgen_cube/cube_theta_quick.out"; FLAG="--quick"
else
  OUT="$ROOT/05-knowledge/results/collatz_procgen_20260923_cube_theta.out"; FLAG=""
fi
{
  echo "# collatz_procgen_20260923_cube_theta -- HYP-9127 lane output"
  echo "# python: $(python3 -c 'import sys,gmpy2,flint; print(sys.version.split()[0], "gmpy2", gmpy2.version(), "python-flint", flint.__version__)')"
  echo "# gp: $(command -v gp >/dev/null && echo 'print(version())' | gp -q 2>/dev/null || echo 'absent')"
  python3 "$EXP/collatz_procgen_20260923_cube_theta_ledger.py" $FLAG
  python3 "$EXP/collatz_procgen_20260923_cube_theta_lattice.py" $FLAG
  python3 "$EXP/collatz_procgen_20260923_cube_theta_pade.py" $FLAG
} > "$OUT" 2>&1
echo "wrote $OUT"
