#!/usr/bin/env python3
"""Generalized +1 barrier census. Residual-type sources whose run end x shadows a periodic point c' (x -> c' 2-adically,
c' = 1 mod 4 so that x is a legal run end). The deletion child h_D has run end y_D = (x+1)/3^D - 1; the D-chain starts at
(D, 3^D - 1) [x = 3^D y + 3^D - 1]. In the limit x* = c' the chain is the exact rational pair (c', (c'+1)/3^D - 1).
Outcome per (c', D): ABSORB (run-uniform certificate!), ANCHOR (universal state analogue), SHIFT, DRIFT, UNRESOLVED."""
import sys
sys.path.insert(0, '../runcompress_20261007')
from transitions import limit_transition, cycle_reps, density, par, T
from fractions import Fraction as Fr
from collections import Counter
cycles = cycle_reps(9, 4)
anchors = []
for w, c, cy in cycles:
    for p in cy:                      # every point of the cycle that is a legal run end (odd, = 1 mod 4 2-adically)
        num, den = p.numerator, p.denominator
        if num % 2 == 0: continue
        if (num * pow(den, -1, 4)) % 4 != 1: continue
        anchors.append((w, p, density(cy)))
print(f"{len(anchors)} legal anchor points (cycle points = 1 mod 4) from {len(cycles)} cycles (U-word sum <= 9, length <= 4)")
summary = Counter(); absorbs = []; best = {}
for w, p, dens in anchors:
    for D in range(1, 41):
        r = limit_transition(D, Fr(3**D - 1), 'source', p, maxit=20000)
        kind = r['kind']
        summary[kind] += 1
        if kind == 'ABSORB': absorbs.append((w, p, D, r['time']))
        if kind == 'ANCHOR':
            k = abs(r['k_entry'])
            if (w, p) not in best or k < best[(w, p)][0]: best[(w, p)] = (k, D, r['entry_time'])
print("outcome counts over (anchor, D <= 40):", dict(summary))
print("ABSORB (run-uniform certificates at a non-trivial anchor):", absorbs[:20], len(absorbs))
print("minimal |debt| of re-anchoring D per anchor (first 20):")
for (w, p), (k, D, t) in sorted(best.items(), key=lambda z: z[1][0])[:20]:
    print(f"   anchor {str(p):>10s} (word {w}): debt {k} via D={D}, entry t={t}")
# --- per-anchor summary: barrier anchors (no ABSORB for any D <= 40) ---
from collections import defaultdict
per = defaultdict(set)
for w, p, D, t in absorbs: per[(w, p)].add(D)
barriers = [(w, p, dens) for (w, p, dens) in anchors if (w, p) not in per]
print(f"barrier anchors (no absorption for D <= 40): {len(barriers)} of {len(anchors)}")
for w, p, dens in barriers[:40]:
    print(f"   anchor {str(p):>12s}  word {w}  odd density {float(dens):.3f}")
print("absorbing anchors with their absorbing depths (first 25):")
for (w, p), Ds in list(per.items())[:25]:
    print(f"   anchor {str(p):>12s}  word {w}: D in {sorted(Ds)}")
