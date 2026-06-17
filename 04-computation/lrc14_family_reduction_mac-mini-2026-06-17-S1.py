#!/usr/bin/env python3
"""
LRC(14) PROVE route: reduce the open family.
Identify large sub-families of primitive 13-sets where M(S) >= 1/14 is provable,
shrinking the open problem to a small "hard core".

Exact gap tool M(S) is validated against a dense grid (see prompt).
THRESH = 1/14. A counterexample is a primitive 13-set with M(S) < 1/14.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random

def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r

def g(S, t):
    return min(nrm(v * t) for v in S)

def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2 * k + 1, 2 * v) <= F(1, 2):
            C.add(F(2 * k + 1, 2 * v)); k += 1
    for i in range(len(S)):
        for j in range(i + 1, len(S)):
            for d in (S[i] + S[j], S[j] - S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2):
                        C.add(F(k, d)); k += 1
    C.add(F(1, 2)); return C

def M(S):
    b = F(0); at = None
    for t in cand(S):
        v = g(S, t)
        if v > b: b = v; at = t
    return b, at

def gcd_all(S):
    return reduce(gcd, S, 0)

# Fast float screen for broad searches
def gf(S, t):
    best = 2.0
    for v in S:
        r = (v * t) % 1.0
        d = r if r <= 0.5 else 1.0 - r
        if d < best: best = d
    return best

THRESH = F(1, 14)
T14 = 1.0 / 14.0

def report_family(name, gen, n_trials, seed):
    """gen(rng) -> sorted primitive 13-set (or None to skip)."""
    rng = random.Random(seed)
    worst = F(10); worstS = None; cnt = 0; below = 0; near = 0
    NEAR = F(1, 13)  # within ~10% above 1/14
    for _ in range(n_trials):
        S = gen(rng)
        if S is None: continue
        if len(set(S)) != 13: continue
        if gcd_all(S) != 1: continue
        cnt += 1
        # float screen first
        # exact M
        m, at = M(S)
        if m < THRESH: below += 1; print("  !!! BELOW 1/14:", S, "M=", m, "at", at)
        if m <= NEAR: near += 1
        if m < worst: worst = m; worstS = S
    print(f"[{name}] tested {cnt} primitive 13-sets; below 1/14: {below}; "
          f"near(<=1/13): {near}")
    if worstS:
        mm, att = M(worstS)
        print(f"  worst M = {mm} = {float(mm):.6f}  at tau={att}")
        print(f"  worstS = {worstS}")
    return worst, worstS

if __name__ == "__main__":
    # sanity
    m, at = M(list(range(1, 14)))
    print("SANITY tight {1..13}: M=", m, "at", at, " (expect 1/14 at 5/14)")
    print()

    # FAMILY (a): contains a multiple of 14
    def gen_a(rng):
        mult = 14 * rng.randint(1, 6)
        rest = set()
        while len(rest) < 12:
            x = rng.randint(1, 70)
            if x != mult: rest.add(x)
        return sorted(rest | {mult})
    report_family("a: contains multiple of 14", gen_a, 5000, 1)
    print()

    # FAMILY (c): spread out (each speed >= 2x previous, geometric-ish)
    def gen_c(rng):
        S = [1]
        for _ in range(12):
            nxt = S[-1] * 2 + rng.randint(0, 3)
            S.append(nxt)
        return sorted(set(S))
    report_family("c: spread out (>=2x growth)", gen_c, 3000, 2)
    print()

    # FAMILY: generic random small speeds (the dangerous regime)
    def gen_rand(rng):
        return sorted(rng.sample(range(1, 30), 13))
    report_family("rand: small random 1..29", gen_rand, 5000, 3)
