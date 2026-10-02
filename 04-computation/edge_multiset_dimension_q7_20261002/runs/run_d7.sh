#!/bin/bash
# d=7 exhaustive runs, engine 0 (bulk), 2 workers, resumable (re-running skips completed chunks)
cd "$(dirname "$0")/.."
KMIN=${1:-1}; KMAX=${2:-12}
for k in $(seq $KMIN $KMAX); do
  python3 py/driver.py 7 $k --workers 2 --target 120 --engine 0 --run runs/d7_e0 >> logs/d7_e0.out 2>&1 || { echo "FAILED k=$k" >> logs/d7_e0.out; exit 1; }
done
echo "ALLDONE $KMIN..$KMAX" >> logs/d7_e0.out
