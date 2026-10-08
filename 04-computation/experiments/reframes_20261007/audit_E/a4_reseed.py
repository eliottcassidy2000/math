#!/usr/bin/env python3
"""Re-run of check A4 (flip-coin fairness on actual orbits) with fresh seeds, to see whether a -2.7 sigma cell is noise."""
import random, sys
sys.argv = sys.argv[:1]
import a_table_check as A
for seed in (1, 2, 3):
    A.rnd = random.Random(1000 + seed)
    print(f"--- seed {1000+seed}")
    A.check_A4(npaths=1500, nsteps=1000)
