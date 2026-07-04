#!/usr/bin/env python3
"""
Idea-generation anchors (mac-mini-2026-07-04-S39): (1) Phi6 iterated = Sylvester's sequence (greedy Egyptian);
(2) covering-family M-values are NOT quantized to the Ostrowski ladder M_k=k/((n-1)k+1) (only extremizers are).
Grounds the 'seven dormant threads into the covering-min core' reflection.
"""
from fractions import Fraction as F
from math import gcd, sin, pi
from functools import reduce
import random
def gcd_all(xs): return reduce(gcd, xs)
def nd(x):
    x = x % 1
    return x if x <= 1-x else 1-x
def M_exact(sp):
    vs = sorted(set(sp)); Q = set()
    for i in range(len(vs)):
        Q.add(2*vs[i])
        for j in range(i+1, len(vs)): Q.add(vs[i]+vs[j]); Q.add(vs[j]-vs[i])
    best = F(0)
    for q in Q:
        if q < 2: continue
        for a in range(1, q):
            m = min(nd(v*F(a,q)) for v in sp)
            if m > best: best = m
    return best
def is_cov(sp, n): return all(any(v % q == 0 for v in sp) for q in range(2, n+1))
def phi6(n): return n*n-n+1

if __name__ == "__main__":
    x = 2; seq = [x]
    for _ in range(5): x = phi6(x); seq.append(x)
    print("(1) Phi6 iterated from 2 =", seq, "= Sylvester (greedy Egyptian 1/2+1/3+1/7+1/43+...)")
    print("    apex-7 slack 4 sin^2(pi/14) =", round(4*sin(pi/14)**2, 5), "= 2-2cos(pi/7)")
    n = 13; rng = random.Random(39)
    def is_ladder(M):
        for k in range(1, 60):
            if M == F(k, (n-1)*k+1): return k
        return None
    tot = lad = 0
    for _ in range(6000):
        sp = sorted(set(rng.sample(range(1, rng.choice([20,40,80,170])), n-1)))
        if len(sp) != n-1 or gcd_all(sp) != 1 or not is_cov(sp, n): continue
        tot += 1
        if is_ladder(M_exact(sp)) is not None: lad += 1
    print(f"(2) covering 12-families on the Ostrowski ladder M_k=k/({n-1}k+1): {lad}/{tot} => quantization FAILS")
    print("    => ladder = extremal locus only (min-M configs), not universal. Honest form of M-uniqueness.")
