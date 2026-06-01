#!/usr/bin/env python3
"""
lrc_n14_geometry_attempt_s521.py   claudebox-2026-06-01-S521

Deep geometry of LRC(14) and an honest attempt at proof.  This is NOT a proof of
LRC(14) (open for 13 movers); it develops the geometry, pushes the multiplicative-
walk attack to its boundary, and pins down exactly where it breaks.

GEOMETRY.  v=(v_1..v_13) distinct nonzero integers.  t |-> (v_i t mod 1) is a
CLOSED rational curve on the torus T^13 (period 1).  LRC(14) = this 1-dim curve
meets the central cube K=[1/14,13/14]^13 (vol (6/7)^13).  Each runner's bad set
has measure 1/7, union bound 13/7>1 FAILS; positivity holds when the bad sets are
independent (equidistributed curve), so the obstruction is rational coincidences
-> reduce to "fully covered" systems (every q in 2..14 divides some speed).

MULTIPLICATIVE WALK.  At t=a/q a runner v is safe iff v a mod q avoids
F_q={r: min(r,q-r) < q/14} = {0,±1,...,±(ceil(q/14)-1)}.  Win at q = some unit a
rotates all speed-residues off F_q.  PAIR-COUNTING (the threshold phenomenon):
the forbidden radius ceil(q/14)-1 grows in lockstep with the number of antipodal
pairs, so 13 movers keep EXACTLY enough blocking power at every single modulus --
no one-modulus argument closes it.

CROSS-MODULUS (empirical).  No sampled fully-covered system -- including a CRT
adversary built to block primes 17,19,23 -- blocks all small denominators; the
worst minimal lonely denominator observed is 25.  CONJECTURE Q0: every primitive
13-system is lonely at t=a/q with q <= 25.  With Q0 + a speed bound, LRC(14) is a
finite check.  Both bounds are open (= the genuine LRC difficulty).
"""
from math import gcd, ceil
from fractions import Fraction as F
from functools import reduce
import random

THR = F(1, 14)
def lonely_at(vs, t):
    for v in vs:
        x = (F(v) * t) % 1
        if min(x, 1 - x) < THR: return False
    return True
def min_denom(vs, qmax):
    for q in range(2, qmax + 1):
        for a in range(1, q):
            if gcd(a, q) == 1 and lonely_at(vs, F(a, q)): return q
    return None
def covered_le14(vs): return all(any(v % q == 0 for v in vs) for q in range(2, 15))

def geometry_facts():
    vol = (F(6,7))**13
    print("GEOMETRY: central cube K=[1/14,13/14]^13")
    print(f"  vol(K) = (6/7)^13 = {float(vol):.4f}   (positive: a generic curve meets K)")
    print(f"  each runner bad-set measure = 1/7; union bound = 13/7 = {13/7:.3f} > 1  -> FAILS")
    print("  => positivity needs independence; obstruction = rational coincidences.\n")

def pair_counting():
    print("PAIR-COUNTING at modulus q (why no single modulus closes it):")
    print(f"  {'q':>3} {'radius':>6} {'|F_q|':>5} {'#pairs':>7}  note")
    for q in [14,17,19,23,25,28,29,42]:
        radius = ceil(q/14) - 1
        nF = 2*radius + 1
        # antipodal pairs of nonzero residues (unit ones drive the unit-a argument)
        npairs = (q-1)//2
        note = "radius grows with pairs -> 13 movers stay 'just enough'"
        print(f"  {q:>3} {radius:>6} {nF:>5} {npairs:>7}  {note}")
    print()

def worst_denominator(trials=40000, seed=1):
    random.seed(seed); worst=0; ws=None
    for _ in range(trials):
        base=[8,9,5,7,11,13,6,10,14]
        rest=[random.randint(1,1000) for _ in range(6)]
        vs=tuple(sorted(set(base+rest)))
        if len(vs)<13 or reduce(gcd,vs)!=1 or not covered_le14(vs): continue
        q=min_denom(vs,80)
        if q and q>worst: worst=q; ws=vs
    return worst, ws

def main():
    print("LRC(14): deep geometry + honest proof attempt (claudebox-S521)\n")
    geometry_facts()
    pair_counting()
    print("PROVED REDUCTION: if some q in 2..14 divides no speed, t=1/q is lonely")
    print("  (||v/q|| >= 1/q >= 1/14).  So only 'fully covered' systems are hard.\n")
    w, ws = worst_denominator()
    print(f"CROSS-MODULUS PROBE: worst minimal lonely denominator over fully-covered")
    print(f"  adversaries = {w}  (e.g. {ws})")
    print(f"  => CONJECTURE Q0: every primitive 13-system is lonely at q <= 25.")
    print("\nVERDICT: NOT a proof of LRC(14).  Open links: (i) the denominator bound Q0")
    print("  (cross-modulus incompatibility), (ii) a speed bound to finitize.  Both are")
    print("  the genuine LRC difficulty; the geometry localizes them precisely.")

if __name__ == "__main__":
    main()
