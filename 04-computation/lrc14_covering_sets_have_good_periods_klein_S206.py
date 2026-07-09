#!/usr/bin/env python3
"""
klein-2026-07-09-S206: DO COVERING SETS ALWAYS HAVE A GOOD PERIOD?  (and does the tight AP never?)

THM-523: LRC(14) reduces to COVERING 13-sets (a multiple of every q in {2..14}); every NON-covering set
has the explicit witness tau=1/q.  S206 verified: the TIGHT locus (M=1/n exactly) is entirely
NON-covering (n=4..7 exhaustive; {1..13} at n=14).

CONSEQUENCE TO TEST.  The covering-case good-period machinery analyses the co-offset cluster
E = {Vmax - v : v in S} (so 0 in E, the observer anchor) and asks for a STRICT good period
    exists j in 1..Vmax-1 :  7 * maxCircGap{e*j mod Vmax}  >  Vmax.
klein-S201 found the cluster E={0..12} at Vmax=13 has NO good period -- and every averaging route died on
it.  But that cluster comes ONLY from the velocity set {1..13}, which is NOT covering.

HYPOTHESIS: restricted to clusters arising from COVERING velocity sets, a strict good period ALWAYS
exists (with margin).  Then the S201 pathology never belonged to the covering case at all.
"""
import random
from math import gcd

def is_cov(S): return all(any(s % q == 0 for s in S) for q in range(2, 15))
def prim(S):
    g = 0
    for s in S: g = gcd(g, s)
    return g == 1

def best_margin(S):
    """max over j of 7*maxCircGap/Vmax  for the co-offset cluster E={Vmax-v}. >1 iff strict good period."""
    Vmax = max(S)
    E = sorted({Vmax - v for v in S})
    best = 0
    for j in range(1, Vmax):
        r = sorted({(e * j) % Vmax for e in E})
        g = max([r[i+1] - r[i] for i in range(len(r)-1)] + [Vmax - r[-1] + r[0]])
        best = max(best, 7 * g / Vmax)
    return best

random.seed(206)
print("margin(S) := max_j 7*maxCircGap/Vmax  on the co-offset cluster.  >1 <=> STRICT good period.\n")

# (1) exhaustive small covering sets
from itertools import combinations
print("(1) EXHAUSTIVE primitive covering 13-sets, speeds <= 18")
n_cov = 0; fails = 0; minmarg = 9e9; worst = None
for S in combinations(range(1, 19), 13):
    if not prim(S) or not is_cov(S): continue
    n_cov += 1
    m = best_margin(S)
    if m <= 1: fails += 1
    if m < minmarg: minmarg, worst = m, S
print(f"    #covering = {n_cov};  #without strict good period = {fails};  min margin = {minmarg:.4f}")
print(f"    worst-margin set: {worst}\n")

# (2) random covering sets at larger speeds
print("(2) RANDOM primitive covering 13-sets, speeds <= 60")
tries = 0; got = 0; fails = 0; minmarg = 9e9; worst = None
while got < 400 and tries < 400000:
    tries += 1
    S = sorted(random.sample(range(1, 61), 13))
    if not prim(S) or not is_cov(S): continue
    got += 1
    m = best_margin(S)
    if m <= 1: fails += 1
    if m < minmarg: minmarg, worst = m, S
print(f"    sampled {got} covering sets;  #without strict good period = {fails};  min margin = {minmarg:.4f}")
print(f"    worst-margin set: {worst}\n")

# (3) the NON-covering tight locus: does it fail?
print("(3) The NON-covering TIGHT sets (the ones every averaging route died on)")
for name, S in [("{1..13} AP (tight, M=1/14)", list(range(1, 14))),
                ("{1,3,4,5,9,...} GW-like",   [1,3,4,5,9,11,12,13,17,19,20,21,25])]:
    S = sorted(S)
    print(f"    {name:>28}: covering={is_cov(S)!s:>5}  margin={best_margin(S):.4f}  "
          f"{'NO strict good period' if best_margin(S)<=1 else 'has good period'}")

print("\nREADING: if #without-strict-good-period == 0 over covering sets while the NON-covering tight AP")
print("fails, then the S201 pathology (and the tight-AP counterexamples that broke arc-count / c<D3 /")
print("smooth-mean / kissing) lives OUTSIDE the covering case. The covering branch -- the ONLY branch")
print("LRC(14) actually needs (THM-523) -- has comfortable good periods.")
