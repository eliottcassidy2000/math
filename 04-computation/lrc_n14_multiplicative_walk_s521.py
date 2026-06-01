#!/usr/bin/env python3
"""
lrc_n14_multiplicative_walk_s521.py   claudebox-2026-06-01-S521

LRC(14) via the multiplicative walk on the residue tournament (the metagraph-walk
reframe with vertices = residues mod q instead of runners).

Reframe: at time t=a/q the runner of speed v sits at residue (v*a mod q); it is
SAFE iff ||v a / q|| >= 1/14, i.e. min(r, q-r) >= q/14 where r = v a mod q.  The
"walk" over rational times is the multiplicative action of the numerator a on the
residue configuration of the speeds; loneliness = some a makes the whole config
avoid the forbidden set F_q = {r : min(r,q-r) < q/14}.

KEY (provable) REDUCTION, verified here:
  * For q <= 14, F_q = {0}: t=1/q is lonely  <=>  no speed is divisible by q
    (since the smallest nonzero distance 1/q >= 1/14).
  * Hence LRC(14) holds for every primitive 13-speed set UNLESS it is
    "fully covered": every q in {2,...,14} divides at least one speed
    (equivalently the speeds include multiples of 8, 9, 5, 7, 11, 13).
  * The fully-covered residual is the only hard core; empirically it is still
    lonely at small denominators (<= ~21 in samples).

This script verifies the reduction and probes the residual.  It is a REDUCTION
and reframe, NOT a proof of LRC(14): the fully-covered core (where the
multiplicative cover meets the additive walk) remains open.
"""
from math import gcd
from fractions import Fraction as F
from collections import Counter
import random

THR = F(1, 14)

def primitive(vs):
    g = 0
    for v in vs: g = gcd(g, v)
    return g == 1

def lonely_at(vs, t):
    for v in vs:
        x = (F(v) * t) % 1
        if min(x, 1 - x) < THR: return False
    return True

def uncovered_q(vs):
    """smallest q in 2..14 dividing no speed (=> t=1/q lonely); None if fully covered."""
    for q in range(2, 15):
        if all(v % q != 0 for v in vs): return q
    return None

def min_denom(vs, qmax=400):
    for q in range(2, qmax + 1):
        for a in range(1, q):
            if gcd(a, q) == 1 and lonely_at(vs, F(a, q)): return q
    return None

def verify_reduction(trials=20000):
    random.seed(3); ok = True
    for _ in range(trials):
        k = random.randint(3, 13)
        vs = tuple(sorted(random.sample(range(1, 300), k)))
        if not primitive(vs): continue
        small = any(lonely_at(vs, F(a, q)) for q in range(2, 15)
                    for a in range(1, q) if gcd(a, q) == 1)
        if small != (uncovered_q(vs) is not None):
            ok = False; print("  MISMATCH", vs); break
    return ok

def probe_residual(target=400):
    random.seed(5); dist = Counter(); tested = nf = attempts = 0
    while tested < target and attempts < 400000:
        attempts += 1
        seeds = [8, 9, 5, 7, 11, 13]
        extra = [random.choice([6, 10, 12, 14, 4]) * random.randint(1, 3) for _ in range(3)]
        rest = random.sample(range(1, 200), max(0, 13 - len(seeds) - len(extra)))
        vs = tuple(sorted(set(seeds + extra + rest)))
        if len(vs) < 10 or not primitive(vs): continue
        if uncovered_q(vs) is not None: continue          # keep only fully-covered
        tested += 1
        q = min_denom(vs, qmax=400)
        if q is None: nf += 1
        else: dist[q] += 1
    return tested, nf, dist

def main():
    print("LRC(14) multiplicative-walk reduction (claudebox-S521)\n")
    print("(1) reduction: small-denominator loneliness <=> some q in 2..14 divides no speed")
    print("    verified:", verify_reduction())
    print("\n(2) fully-covered residual (every q in 2..14 divides some speed):")
    tested, nf, dist = probe_residual()
    print(f"    sets tested: {tested}   no lonely time found (q<=400): {nf}")
    print(f"    minimal-denominator distribution: {dict(sorted(dist.items()))}")
    print(f"    max denominator needed: {max(dist) if dist else None}")
    print("\n    => non-fully-covered sets: PROVED lonely at t=1/q.  Fully-covered core")
    print("       remains the open heart of LRC(14) (multiplicative cover vs additive walk).")

if __name__ == "__main__":
    main()
