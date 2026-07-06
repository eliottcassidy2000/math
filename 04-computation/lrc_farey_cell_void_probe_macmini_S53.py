#!/usr/bin/env python3
"""
mac-mini-2026-07-05-S53 -- HYP-4100: the FAREY-CELL reframe of the second-value
gap (G) + the falsification probe.

LEMMA (one line, verified here): 1/13 and 2/25 are Farey neighbors
(det = -1), so every rational strictly inside (1/13, 2/25) has reduced
numerator >= 3 and denominator >= 38 (minimum = the mediant 3/38).
Hence (G) <=> no 12-set's M lands in the cell <=> no 12-set sustains a
clearance->=3 witness at q in (12.5c, 13c).

PROBE (enumerate-don't-sample): structured 12-set families searched for
M strictly inside (1/13, 2/25); plus targeted search for M = 3/38 exactly.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import itertools, sys
sys.path.insert(0, '04-computation')
from lonely_profile import profile

LO, HI, MED = F(1, 13), F(2, 25), F(3, 38)
print(f"Farey adjacency: 1*25 - 13*2 = {1*25 - 13*2} (neighbors ✓); mediant = {F(1+2, 13+25)}")
def M_of(S):
    for cap in (10, 6, 4, 3):
        m = profile(sorted(S), F(1, cap)).M()
        if m is not None: return m
    return None

in_gap = []
tested = 0
def probe(S, fam):
    global tested, in_gap
    S = sorted(set(S))
    if len(S) != 12 or reduce(gcd, S) != 1: return
    tested += 1
    M = M_of(S)
    if M is not None and LO < M < HI:
        in_gap.append((fam, S, M))
        print(f"  !! IN-GAP: [{fam}] {S}: M = {M}")

# families: drops of {1..13}; single lifts of {1..12}; the 14-fold ladder bases;
# dilation-mixes; deep-well descendants; two-lift combos (small grid)
for d in range(1, 14):
    probe([v for v in range(1, 14) if v != d], 'drop13')
for r in range(1, 13):
    for k in range(1, 14):
        probe([v for v in range(1, 13) if v != r] + [r + 13 * k], 'lift')
for r in range(1, 13):
    probe([v for v in range(1, 13) if v != r] + [14 * r], '14fold')
for c in (2, 3):
    for x in range(1, 30):
        probe([c * v for v in range(1, 12)] + [x], 'dilmix')
for r1, r2 in itertools.combinations(range(1, 13), 2):
    for k1 in (1, 2):
        for k2 in (1, 2):
            probe([v for v in range(1, 13) if v not in (r1, r2)] + [r1 + 13*k1, r2 + 13*k2], '2lift')
# targeted: can ANY 12-set hit 3/38 exactly? residues mod 38 must clear +-2;
# try sets built from the 38-grid structure
for a in range(1, 19):
    if gcd(a, 38) != 1: continue
    # runners at residues 3..35 (|res| >= 3): pick 12 with small representatives
    S = []
    for res in [3,4,5,6,7,8,9,10,11,12,13,14]:
        # smallest v with a*v = res mod 38
        ainv = pow(a, -1, 38)
        v = (res * ainv) % 38
        S.append(v if v > 0 else 38)
    probe(S, f'grid38-a{a}')
print(f"\ntested {tested} primitive 12-sets; IN-GAP hits: {len(in_gap)}")
print("(G) Farey-cell reframe:", "SUPPORTED (void empty in all families)" if not in_gap else "FALSIFIED CANDIDATES ABOVE")
