#!/usr/bin/env python3
"""
klein-2026-07-05-S134 (HYP-4099) - THE TOWER LEAF SCAN: which compressed covering
primitive families survive EVERY (base, tops)-split of the multi-killer window?

Split criterion (lonely_of_window_multi + cite_margin): choose tops T (1 <= |T| <= 6),
base = rest (13-|T| runners, cited at beta = 1/(14-|T|)), B = max base |v|:
  CLOSES the family iff  sum_{j in T} 3/(7 v_j)  <  2*((beta-1/14)/B)*((7-|T|)/7).
A family is a TOWER LEAF iff NO split closes it (and no repeat: repeats -> citation).

Scan structured compressed covering primitive families; report the leaf population and
its max spread (predicted <= ~C(13,6) = 1716 by the threshold telescope).
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd

def covering(V):
    return all(any(x % q == 0 for x in V) for q in range(2, 15))

def compressed(V):
    W = sorted(V)
    return W[-1] <= 13 * W[-2]

def closes(V):
    """Does SOME split close V? V sorted ascending, distinct (no repeats here)."""
    n = len(V)
    for tsize in range(1, 7):
        beta = F(1, 14 - tsize)
        # best split of this size: tops = the tsize LARGEST (maximizes v_j, minimizes B)
        tops = V[n - tsize:]
        base = V[:n - tsize]
        B = base[-1]
        lhs = sum(F(3, 7 * vj) for vj in tops)
        rhs = 2 * ((beta - F(1, 14)) / B) * F(7 - tsize, 7)
        if lhs < rhs:
            return True, tsize
    return False, None

# structured seeds: dilated APs + killers, blocks, mixed
leaves = []
checked = 0
seeds = set()
for c in range(1, 30):
    for extra in range(1, 400):
        V = tuple(sorted(set([c * j for j in range(1, 13)] + [extra])))
        if len(V) == 13: seeds.add(V)
for N in range(2, 200, 7):
    seeds.add(tuple(range(N, N + 13)))          # blocks
    seeds.add(tuple(sorted(set(list(range(1, 13)) + [N])))[:13])
for c in (1, 2, 3, 5, 7):
    for k in (2, 3, 14, 26, 84, 156, 182):
        V = tuple(sorted(set([c * j for j in range(1, 14) if j != 12] + [k * 12])))
        if len(V) == 13: seeds.add(V)
for V in sorted(seeds):
    if len(set(V)) < 13: continue      # repeats -> citation leg
    if not covering(V) or not compressed(V) or gcd(*V) != 1: continue
    checked += 1
    ok, tsize = closes(V)
    if not ok:
        leaves.append(V)

print(f"checked {checked} structured compressed covering primitive families")
print(f"TOWER LEAVES (no split closes): {len(leaves)}")
sp = 0
for V in leaves[:20]:
    s = V[-1] / V[0]
    sp = max(sp, s)
    print(f"  {V}  spread {s:.1f}")
if leaves:
    allsp = max(v[-1]/v[0] for v in leaves)
    print(f"max leaf spread: {allsp:.1f} (predicted bound ~C(13,6) = 1716)")
