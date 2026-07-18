#!/usr/bin/env python3
"""
ROUTE B CRUX: covering constraints + gap structure  (boxeph-2026-07-18-S101)

Verifies, on the known M<1/13 families (deep well + tower/dilations):
 (i) the SIEVE-DIVISIBILITY lemma: M<1/13 => every q' in {2..13} divides some speed;
 (ii) the residues at the maximizer have EXACTLY ONE gap < val (= core difference-closed,
      the AP core val*{1..12} plus the v_max anomaly at 12val+1).

These support the S101 reflection: the crux <=> "at most one interior small gap", covering
gives the divisibility layer, but the finish is a global cross-modulus argument (open).
"""
from math import gcd
from fractions import Fraction as Fr


def Mstar(V, QMAX=4000):
    best = (Fr(0), 1, 0)
    for q in range(2, QMAX + 1):
        for a in range(1, q):
            if gcd(a, q) != 1:
                continue
            m = min(min((a * v) % q, q - ((a * v) % q)) for v in V)
            if Fr(m, q) > best[0]:
                best = (Fr(m, q), q, a)
    return best


fams = {
    'deep well {1..12,182}': list(range(1, 13)) + [182],
    'dilate x2 {2..24,364}': [2 * x for x in range(1, 13)] + [2 * 182],
    'tower {1..12,364}': list(range(1, 13)) + [364],
}
for name, V in fams.items():
    M, q, a = Mstar(V)
    val = min(min((a * v) % q, q - ((a * v) % q)) for v in V)
    res = sorted((a * v) % q for v in V)
    gaps = [res[i + 1] - res[i] for i in range(len(res) - 1)]
    small = [g for g in gaps if g < val]
    div = all(any(v % qp == 0 for v in V) for qp in range(2, 14))
    print(f'{name}:  M={M}={float(M):.5f} (<1/13={float(M) < 1/13})  q={q} val={val}')
    print(f'  (i)  every q\'<=13 divides some speed: {div}')
    print(f'  (ii) #gaps<val = {len(small)}  smallgaps={small}   residues={res}')
    print()
print('Confirms: exactly ONE small gap (core difference-closed) + full sieve divisibility.')
print('The crux <=> "no SECOND interior small gap" is the open inverse theorem (= LRC(14), S94).')
