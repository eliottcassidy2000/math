#!/usr/bin/env python3
"""Brute-force long runs of the eleven rule systems (partial loop acceleration), per-machine time budget.
Reports: halted?, macro steps, final parameter bit lengths.  Session opus-2026-10-25-S56."""
import sys, time, json, os
from deepmind_cryptids_20261025 import step, START, HALT, MACHINES
HERE = os.path.dirname(os.path.abspath(__file__))
budget_s = float(sys.argv[1]) if len(sys.argv) > 1 else 60
cap_bits = int(sys.argv[2]) if len(sys.argv) > 2 else 400000
res = {}
for nm in MACHINES:
    p = START[nm]; n = 0; t0 = time.time(); halted = False
    while time.time() - t0 < budget_s:
        for _ in range(1000):
            if p == HALT:
                halted = True; break
            p = step(nm, p); n += 1
        if halted or p == HALT:
            halted = True; break
        if max(x.bit_length() for x in p) > cap_bits:
            break
    res[nm] = {'halted': halted, 'macro_steps': n, 'final_bits': None if halted else [x.bit_length() for x in p],
               'seconds': round(time.time() - t0, 1)}
    print(nm, res[nm], flush=True)
json.dump(res, open(os.path.join(HERE, f'long_run_{int(budget_s)}s.json'), 'w'), indent=1)
