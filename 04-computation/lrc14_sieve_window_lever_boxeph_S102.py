#!/usr/bin/env python3
"""
REFUTING the S101 "free sieve window" lever  (boxeph-2026-07-18-S102)

S101 named the finishing lever: two interior small gaps (at the maximizer modulus q)
force a free SIEVE window at some q'<=13, hence a witness M>=1/13, contradiction.

This script REFUTES it: a 2-small-gap configuration that is SIEVE-COMPLETE (no witness
at any q'<=13) but is beaten at q=24 (M=1/12 > 1/13). So the crux's forcing lives at a
MEDIUM modulus 13 < q' < q, not at the sieve moduli. The lever, as named, is false.
"""
from math import gcd
from fractions import Fraction as Fr


def maxband(V, qs):
    best = (Fr(0), 0, 0)
    for q in qs:
        for a in range(1, q):
            if gcd(a, q) != 1:
                continue
            m = min(min((a * v) % q, q - ((a * v) % q)) for v in V)
            if Fr(m, q) > best[0]:
                best = (Fr(m, q), q, a)
    return best


# deep-well residue AP {14*1..14*12, 169} with a twin inserted (+15, -168) => TWO small gaps at q=183
twin2 = [14, 15, 28, 42, 56, 70, 84, 98, 112, 126, 140, 154, 169]
print('twin@2 family (two small gaps at q=183):', twin2)
lo = maxband(twin2, range(2, 14))        # sieve moduli q'<=13
mid = maxband(twin2, range(14, 27))
allm = maxband(twin2, range(2, 600))
print(f'  best over q\'<=13 (SIEVE):  M={lo[0]}={float(lo[0]):.5f} at q={lo[1]}  ->  >=1/13? {float(lo[0]) >= 1/13}')
print(f'  best over 14<=q\'<=26     :  M={mid[0]}={float(mid[0]):.5f} at q={mid[1]}')
print(f'  true max over all q\'     :  M={allm[0]}={float(allm[0]):.5f} at q={allm[1]}')
print(f'  sieve killed by: 13|169={169 % 13 == 0}, 11|154={154 % 11 == 0}, 7|14={14 % 7 == 0}, ...')
print()
print('CONCLUSION: sieve-complete (no q\'<=13 witness) yet M=1/12>1/13 via q=24.')
print('=> two small gaps do NOT force a sieve window; the crux forces at a MEDIUM modulus q>13.')
print('   The S101 "free sieve window" lever is REFUTED. Crux = medium-modulus inverse theorem (open).')
