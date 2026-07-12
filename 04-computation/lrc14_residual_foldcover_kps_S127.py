# -*- coding: utf-8 -*-
# kps-2026-07-11-S127cont46: attacking the RESIDUAL 8.5% of Route B as a bounded fold-class covering on the
# small COPRIME-to-q sub-family (owner's named next step).
#
# From cont.45: the coprime-reduction => clearing at q via a UNIT multiplier p  <=>  [q nmid v_i for all i]
# AND the COPRIME-to-q sub-family misses a unit-fold-class. The elementary guarantee (#coprime <= phi(q)/2-1)
# clears 91.5%; the residual 8.5% has #coprime >= phi(q)/2 at every easy q, so it clears (when it does) only
# by a genuine FOLD-CLASS MISS (coprime runners COLLIDE into fewer than phi(q)/2 classes). THIS characterizes
# the residual precisely and looks for the bounded mechanism:
#  (1) exact per-q "unit-clearing" = (# missed fold-classes by coprime sub-family) -- the composite analog of
#      klein THM-718; verify it matches true clearing (up to non-unit-p clearing).
#  (2) the OBSTRUCTION (never unit-clears) = coprime-to-q sub-family HITS ALL phi(q)/2 fold-classes at EVERY
#      window q. Characterize: does DC make this impossible? find the always-missed q.
#  (3) collision structure: why do the coprime runners collide (miss) even when >= phi(q)/2 of them?
#  (4) backstop: non-unit-p clearing for any family the unit-analysis misses.
import random
from math import gcd
from functools import reduce
from collections import Counter

def is_DC(v):  return all(any(x % d == 0 for x in v) for d in range(2, 15))
def prim(v):   return reduce(gcd, v) == 1
def longest_run(v):
    v = sorted(set(v)); b = m = 1
    for i in range(1, len(v)):
        if v[i] == v[i-1] + 1: m += 1; b = max(b, m)
        else: m = 1
    return b
def phi(q): return sum(1 for a in range(1, q) if gcd(a, q) == 1)
def fold(a, q): return min(a % q, (-a) % q)   # antipodal fold-class rep (0..q//2)

def coprime_foldclasses(v, q):
    return {fold(vi % q, q) for vi in v if gcd(vi, q) == 1}   # unit-fold-classes hit by coprime runners
def unit_clears(v, q):
    # clears via a UNIT multiplier: q nmid every v_i AND coprime sub-family misses a unit-fold-class
    if any(vi % q == 0 for vi in v): return False
    hit = coprime_foldclasses(v, q)
    total = phi(q) // 2
    return len(hit) < total
def clears_any_p(v, q):           # TRUE clearing (any multiplier, unit or not)
    lo = -(-q // 14)
    return any(all(lo <= (vi * p) % q <= q - lo for vi in v) for p in range(1, q))

def rand_spread_DC(rng, vmax):
    for _ in range(20000):
        v = sorted(rng.sample(range(1, vmax + 1), 13))
        if prim(v) and is_DC(v) and longest_run(v) <= 7:
            return v
    return None

def main():
    rng = random.Random(20260712)
    window = [q for q in range(15, 30) if q % 14]
    fams = []
    for vmax in (40, 60, 90, 130, 200):
        for _ in range(900):
            v = rand_spread_DC(rng, vmax)
            if v: fams.append(v)
    fams = [list(t) for t in {tuple(v) for v in fams}]
    print(f"spread primitive DC families: {len(fams)}")

    # (1) does UNIT-clearing (fold-class miss) match TRUE clearing? where does non-unit-p add clearing?
    only_nonunit = 0
    for v in fams:
        u = any(unit_clears(v, q) for q in window)
        t = any(clears_any_p(v, q) for q in window)
        if t and not u: only_nonunit += 1
    print(f"(1) families that clear but NOT via any unit-multiplier fold-class miss: {only_nonunit}/{len(fams)}")
    print(f"    => unit-clearing (coprime fold-class miss) {'CAPTURES ALL' if only_nonunit==0 else 'misses '+str(only_nonunit)} true clearing on the window")

    # (2) the residual = families where the elementary guarantee (#coprime<=phi/2-1) fails at every q
    def guar(v, q): return (all(vi % q != 0 for vi in v)) and (sum(1 for vi in v if gcd(vi, q) == 1) <= phi(q)//2 - 1)
    resid = [v for v in fams if not any(guar(v, q) for q in window)]
    print(f"(2) residual (elementary guarantee fails everywhere): {len(resid)} ({100*len(resid)/len(fams):.1f}%)")
    # for the residual, the MIN over window of (#fold-classes MISSED by coprime sub-family) -- must be >=1 to clear
    miss_profile = []
    for v in resid:
        best = max((phi(q)//2 - len(coprime_foldclasses(v, q))) if all(vi%q for vi in v) else -99 for q in window)
        miss_profile.append(best)
    print(f"(3) residual: max-over-window (#missed fold-classes by coprime sub-family) distribution: {dict(sorted(Counter(miss_profile).items()))}")
    print(f"    (>=1 means SOME window q has a fold-class miss => unit-clears. -99 excluded via q|v_i)")
    # which window q carries the residual's clearing? (the miss happens where?)
    carry = Counter()
    for v in resid:
        for q in window:
            if unit_clears(v, q): carry[q] += 1; break
    print(f"(4) residual first-unit-clearing-q distribution: {dict(sorted(carry.items()))}")

    # (5) the COLLISION mechanism: at the residual's clearing q, #coprime runners vs #distinct fold-classes
    colls = []
    for v in resid:
        for q in window:
            if unit_clears(v, q):
                nc = sum(1 for vi in v if gcd(vi, q) == 1)
                nf = len(coprime_foldclasses(v, q))
                colls.append(nc - nf)   # collisions = coprime runners sharing a fold-class
                break
    print(f"(5) at the clearing q: collisions (#coprime - #distinct-foldclasses) distribution: {dict(sorted(Counter(colls).items()))}")
    print(f"    => residual clears BECAUSE coprime runners collide (>=1 collision frees a fold-class); median collisions = {sorted(colls)[len(colls)//2] if colls else 'NA'}")

    # (6) worst case: any family clearing NOWHERE in the window? (would be an LRC(14) counterexample)
    fail = [v for v in fams if not any(clears_any_p(v, q) for q in window)]
    print(f"(6) families clearing NOWHERE in [15,29]: {len(fail)} (must be 0 -- else LRC(14) counterexample)")
    if fail: print("    !!!", fail[0])

if __name__ == "__main__":
    main()
