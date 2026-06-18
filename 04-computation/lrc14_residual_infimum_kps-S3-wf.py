#!/usr/bin/env python3
"""
LRC(14) S3 residual -- (a) exact M infimum over the residual (looseness),
(b) bound on r (#sub-bands) and the structural multi-band theorem,
(c) whether a global witness always exists (the real claim we can defend).

We use the EXACT M tool (binding-pair candidate set) on the residual sets to confirm
M >= 1/14 with margin, and report the minimum exact M. This is the decisive looseness
check: if min exact M over thousands of residual S3 sets is comfortably > 1/14, the
residual is non-dangerous (as THM-526 claims, min ~ 4/31).
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations
import random
from collections import Counter

def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r
def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2*k+1, 2*v) <= F(1, 2): C.add(F(2*k+1, 2*v)); k += 1
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for d in (S[i]+S[j], S[j]-S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2): C.add(F(k, d)); k += 1
    C.add(F(1, 2)); return C
def Mval(S):
    b = F(0); at = None
    for t in cand(S):
        v = min(nrm(x * t) for x in S)
        if v > b: b = v; at = t
    return b, at

def merge(ivs):
    ivs = sorted(ivs); m = []
    for a, b in ivs:
        if m and a <= m[-1][1]: m[-1] = (m[-1][0], max(m[-1][1], b))
        else: m.append((a, b))
    return m
def small_safe_arcs(P, h=F(1, 14)):
    iv = []
    for u in P:
        for j in range(0, u):
            c = F(j, u); a = (c - h/u) % 1; b = (c + h/u) % 1
            if a < b: iv.append((a, b))
            else: iv.append((a, F(1))); iv.append((F(0), b))
    m = merge(iv); safe = []; prev = F(0)
    for a, b in m:
        if a > prev: safe.append((prev, a))
        prev = max(prev, b)
    if prev < 1: safe.append((prev, F(1)))
    return safe
def teeth_in(u, lo, hi, h=F(1, 14)):
    out = []
    jmin = int(lo * u) - 1; jmax = int(hi * u) + 1
    for j in range(jmin, jmax + 1):
        c = F(j, u); a = c - h/u; b = c + h/u
        aa = max(a, lo); bb = min(b, hi)
        if aa < bb: out.append((aa, bb))
    return out
def find_global_witness(S):
    P = [u for u in S if u <= 13]; L = [u for u in S if u > 13]
    safe = small_safe_arcs(P); best = None
    for (a, b) in safe:
        dgr = merge([t for u in L for t in teeth_in(u, a, b)]); prev = a
        for (x, y) in dgr:
            if x > prev:
                w = x - prev
                if best is None or w > best[0]: best = (w, prev, x)
            prev = max(prev, y)
        if b > prev:
            w = b - prev
            if best is None or w > best[0]: best = (w, prev, b)
    if best is None: return None
    w, lo, hi = best
    return ((lo + hi) / 2, w)

def gcd_all(S): return reduce(gcd, S, 0)
def is_covering(S): return all(any(v % q == 0 for v in S) for q in range(2, 15))
def classify(S):
    S = sorted(set(S)); Vmin, Vmax = S[0], S[-1]
    k = sum(1 for v in S if v > 13)
    if k <= 1: return 'S1'
    if Vmax < 13 * Vmin: return 'S2'
    return 'S3'
def missing_q(P): return [q for q in range(2, 15) if not any(v % q == 0 for v in P)]
def single_gap_closable(S):
    P = [u for u in S if u <= 13]; L = [u for u in S if u > 13]
    Vmin, Vmax = min(L), max(L); s = Vmax - Vmin
    safe = small_safe_arcs(P); K = 0
    while (13 * Vmin - Vmax > 14 * K * s) if s > 0 else (K == 0):
        lo = F(14*K+1, 14)/Vmin; hi = F(14*K+13, 14)/Vmax
        if lo < hi:
            for (a, b) in safe:
                if max(a, lo) < min(b, hi): return True
        if hi >= 1 or (s == 0 and K > 0): break
        K += 1
        if K > 14*Vmax: break
    return False

def gen_constructive(seed=0, target=2000, Vrange=(50, 2000)):
    rng = random.Random(seed); out = []; tries = 0
    smalls = []
    base = list(range(1, 14))
    for sz in (11, 10, 9, 8, 7):
        for P in combinations(base, sz): smalls.append(list(P))
    while len(out) < target and tries < target * 300:
        tries += 1
        P = rng.choice(smalls); c = 13 - len(P)
        if c < 2: continue
        miss = missing_q(P); V = rng.randint(*Vrange)
        spread = rng.choice([14, 20, 28, 40, 56, 70, 90])
        window = list(range(V, V + spread + 1))
        if len(window) < c: continue
        cluster = rng.sample(window, c)
        if not all(any(v % q == 0 for v in cluster) for q in miss): continue
        S = sorted(set(P) | set(cluster))
        if len(S) != 13 or gcd_all(S) != 1: continue
        if not is_covering(S) or classify(S) != 'S3': continue
        out.append(S)
    return out

if __name__ == "__main__":
    # Build residual (not single-gap closable) and check exact M + global witness.
    allS = []
    for (lo, hi) in [(50, 200), (200, 800)]:   # keep Vmax modest so exact M is fast
        allS += gen_constructive(seed=lo+13, target=1500, Vrange=(lo, hi))
    residual = [S for S in allS if not single_gap_closable(S)]
    print(f"n S3 = {len(allS)}, residual (not single-gap closable) = {len(residual)}")
    minM = F(10); minS = None; gw_ok = 0; rdist = Counter()
    nM = 0
    for S in residual:
        gw = find_global_witness(S)
        if gw and gw[1] > 0: gw_ok += 1
        m, at = Mval(S)
        nM += 1
        if m < minM: minM = m; minS = S
    print(f"global witness (meas>0) on residual: {gw_ok}/{len(residual)}")
    print(f"exact M computed on {nM} residual sets")
    print(f"min exact M over residual = {minM} = {float(minM):.5f}")
    print(f"   threshold 1/14 = {float(F(1,14)):.5f}; ratio min M/(1/14) = {float(minM*14):.3f}")
    print(f"   at S = {minS}")
