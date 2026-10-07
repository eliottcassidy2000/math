#!/usr/bin/env python3
"""Extend the +1 barrier (THM-4601 (ii) / HYP-9240) to D <= DMAX with exact rational limit chains."""
import sys
from twoanchor_core import limit_chain
DMAX = int(sys.argv[1]); D0 = int(sys.argv[2]) if len(sys.argv) > 2 else 1
inph = outph = to0 = 0; maxN = 0; mink = None; absorbed = []
for D in range(D0, DMAX + 1):
    r = limit_chain(D)
    maxN = max(maxN, r['N'])
    if r.get('kind') == 'to0': to0 += 1; continue
    if r['inphase']:
        inph += 1
        if r['k'] == 0: absorbed.append(D)
        if mink is None or r['k'] < mink[0]: mink = (r['k'], D)
    else: outph += 1
print(f"D in [{D0},{DMAX}]: absorbed {absorbed}; to0 {to0}; in-phase {inph}; out-of-phase {outph}; min in-phase debt {mink}; max landing N {maxN}")
