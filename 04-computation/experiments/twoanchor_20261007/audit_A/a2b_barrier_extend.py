#!/usr/bin/env python3
"""Audit A: optional extension of the +1 barrier with the audit's own scaled limit chain (a2_barrier.limit_chain_scaled),
D in [D0, D1]; cross-checks the session's barrier_extend.out (D in [3001, 8000])."""
import sys
from a2_barrier import limit_chain_scaled
D0, D1 = int(sys.argv[1]), int(sys.argv[2])
absorbed = []; to0 = []; inph = outph = 0; maxN = (0, None); mink = None; rat = []
for D in range(D0, D1 + 1):
    r = limit_chain_scaled(D)
    assert r['kind'] in ('cycle', 'to0'), r
    if r['N'] > maxN[0]: maxN = (r['N'], D)
    if r['kind'] == 'to0': to0.append(D); continue
    if r['absorbed']: absorbed.append(D)
    if r['inphase']:
        inph += 1; rat.append(r['k']/D)
        if mink is None or r['k'] < mink[0]: mink = (r['k'], D)
    else: outph += 1
print(f"D in [{D0},{D1}]: absorbed {absorbed}; to0 {to0}; in-phase {inph}; out-of-phase {outph}; "
      f"min in-phase debt (k, D) {mink}; max landing (N, D) {maxN}; in-phase k/D in [{min(rat):.3f}, {max(rat):.3f}]")
