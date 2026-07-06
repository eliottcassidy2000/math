#!/usr/bin/env python3
"""lrc_cluster_gcd_kps_S17.py -- HYP-4217 part 2: the cluster-gcd ladder.
Planted big-gcd clusters (d > 50*Sum(C)/(25-4|C|)) must be loose, witness at a
periodic copy t0 + j/d of the T-citation point. See the draft
03-artifacts/drafts/cluster-gcd-ladder-kps-S17.md for the theorem + proof."""
from fractions import Fraction
from math import gcd
from functools import reduce
import random
rng = random.Random(17)
def margin_at(ws, p, q):
    return Fraction(min(min((w*p) % q, q - (w*p) % q) for w in ws), q)
def witness_at_copies(W, T, d, Q0=250):
    for q in range(2, Q0+1):
        for a in range(1, q):
            if gcd(a,q) > 1: continue
            if margin_at(T, a, q) >= Fraction(1,12):
                for j in range(d):
                    p2, q2 = a*d + j*q, q*d
                    g = gcd(p2, q2)
                    if margin_at(W, p2//g, q2//g) >= Fraction(2,25):
                        return (a, q, j)
    return None
viol = tested = 0
for trial in range(25):
    nc = rng.randint(1, 5)
    C = sorted(rng.sample(range(1, 15), nc))
    bound = 50*sum(C)//(25 - 4*nc)
    d = bound + rng.randint(1, 50)
    T = [d*x for x in sorted(rng.sample(range(1, 40), 12-nc))]
    W = T + C
    if reduce(gcd, W) != 1: continue
    tested += 1
    if witness_at_copies(W, T, d) is None:
        viol += 1; print("FAIL", C, d)
print(f"{tested} planted families, {viol} failures (0 = ladder mechanism verified)")
