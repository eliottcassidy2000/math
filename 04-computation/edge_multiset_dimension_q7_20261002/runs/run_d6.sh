#!/bin/bash
# d=6 validation: k=1..15 with engine 0 (bulk), then k=14,15 again with engine 1 (per-leaf) as a cross-check
cd "$(dirname "$0")/.."
for k in $(seq 1 15); do python3 py/driver.py 6 $k --workers 2 --target 60 --engine 0 --run runs/d6_e0 >> logs/d6_e0.out 2>&1 || echo "FAILED k=$k" >> logs/d6_e0.out; done
for k in 14 15; do python3 py/driver.py 6 $k --workers 2 --target 60 --engine 1 --run runs/d6_e1 >> logs/d6_e1.out 2>&1 || echo "FAILED k=$k" >> logs/d6_e1.out; done
echo ALLDONE >> logs/d6_e0.out
