#!/usr/bin/env python3
"""
TARGETED search: does a PRIMITIVE tight family with q*=28 exist? (mac-mini-2026-07-03-S31, numpy-fast)
Even block {2,4,..,26} is tight at q*=28 but imprimitive. Try to primitivize => primitive tight q*=28?
If YES => tight locus > {AP,GW}, 'primitive tight => q*=14' FALSE (override THM-523). If NO (thorough) => support.
numpy float-M pre-filter (|M-1/14|<2e-4) then EXACT rational confirm; only near-tight get exact-M.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations
import numpy as np, random

def gcd_all(xs): return reduce(gcd, xs)
_G = np.arange(1, 8400)/8400.0   # grid resolves 1/8400 << 1/28
def approxM(sp):
    v = np.array(sp, dtype=np.float64)
    ph = np.outer(v, _G) % 1.0
    return np.minimum(ph, 1.0-ph).min(axis=0).max()
def nd(x):
    x = x % 1
    return x if x <= 1-x else 1-x
def M_exact(sp):
    vs = sorted(set(sp)); Q = set()
    for i in range(len(vs)):
        Q.add(2*vs[i])
        for j in range(i+1, len(vs)):
            Q.add(vs[i]+vs[j]); Q.add(vs[j]-vs[i])
    best = F(0)
    for q in Q:
        if q < 2: continue
        for a in range(1, q):
            m = min(nd(v*F(a,q)) for v in sp)
            if m > best: best = m
    return best
def qstar_set(sp, M):
    vs = sorted(set(sp)); Q = set()
    for i in range(len(vs)):
        Q.add(2*vs[i])
        for j in range(i+1, len(vs)):
            Q.add(vs[i]+vs[j]); Q.add(vs[j]-vs[i])
    ds = set()
    for q in sorted(Q):
        for a in range(1, q):
            if gcd(a,q)!=1: continue
            if min(nd(v*F(a,q)) for v in sp) == M: ds.add(F(a,q).denominator)
    return sorted(ds)

target = F(1,14); tf = 1.0/14
def confirm(S, hits):
    if abs(approxM(S) - tf) > 2e-4: return
    if M_exact(S) == target:
        hits.append((sorted(S), qstar_set(S, target)))

if __name__ == "__main__":
    import sys
    def out(*a): print(*a); sys.stdout.flush()
    hits = []
    EB = [2*i for i in range(1,14)]
    out("(a) replace 1 even block runner with an odd one (minimal primitivization):")
    n1 = 0
    for drop in range(13):
        for w in range(1, 400, 2):
            S = EB[:drop] + EB[drop+1:] + [w]
            if len(set(S)) != 13 or gcd_all(S) != 1: continue
            n1 += 1; confirm(S, hits)
    out(f"   tested {n1}; primitive tight found: {len(hits)}   {hits[:6]}")

    out("(b) replace 2 even runners with 2 odds (float-prefiltered):")
    h0 = len(hits); n2 = 0
    for d1,d2 in combinations(range(13),2):
        remain = [EB[i] for i in range(13) if i not in (d1,d2)]
        for w1,w2 in combinations(range(1,120,2),2):
            S = remain + [w1,w2]
            if len(set(S)) != 13 or gcd_all(S) != 1: continue
            n2 += 1; confirm(S, hits)
    out(f"   tested {n2}; new primitive tight found: {len(hits)-h0}   {[h for h in hits[h0:]][:6]}")

    out("(c) 2*U + odd tighteners, and broad random mixed-parity (float-prefiltered):")
    rng = random.Random(31); h1 = len(hits); n3 = 0
    for _ in range(60000):
        ne = rng.choice([9,10,11,12])   # #even runners
        evn = rng.sample([2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32], ne)
        odd = rng.sample([1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35], 13-len(set(evn)))
        S = sorted(set(evn) | set(odd))
        if len(S) != 13 or gcd_all(S) != 1: continue
        n3 += 1; confirm(S, hits)
    out(f"   tested {n3}; new primitive tight found: {len(hits)-h1}   {[h for h in hits[h1:]][:6]}")

    out("\n" + "="*82)
    q_gt14 = [(S,ds) for S,ds in hits if any(d>14 for d in ds)]
    out(f"TOTAL primitive tight families found: {len(hits)}; with q*>14: {len(q_gt14)}")
    if q_gt14:
        out("  !!! tight locus BIGGER than {AP,GW} — primitive tight WITH q*>14 exists:")
        for S,ds in q_gt14[:12]: out(f"     {S}  q*={ds}  gcd={gcd_all(S)}")
    else:
        out("  => NONE with q*>14. Every primitive tight family found sits at q*=14.")
        out("     Supports the CONFINEMENT: primitive tight => q*=14 (rigidity = finite mod-14 problem).")
    if hits:
        out(f"  (all primitive tight found, q* denominators): {sorted(set(tuple(ds) for _,ds in hits))}")
