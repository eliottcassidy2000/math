#!/usr/bin/env python3
"""
LRC(14) S3 residual -- the PROVABLE THEOREM and the honest bounded residual.

We establish a clean, rigorous SUFFICIENT lemma (the "two-gap / cluster-collapse
window" lemma) that closes a sub-regime of S3 with a genuine proof, and we delimit
exactly what is left.

================================ THE LEMMA ==================================
LEMMA (cluster-collapse global witness).  Let S=P∪L, P⊆{1..13}, L a cluster with
Vmin=min L, Vmax=max L, s=Vmax-Vmin (spread), c=|L|. Suppose there is an integer k and
a closed interval J=[t0-r, t0+r] ⊆ (0,1) such that:
  (C1) every cluster runner is safe on J:  for all u∈L and all t∈J, ||u t|| >= 1/14.
       A clean sufficient condition (single shared gap): there is an integer K with
            K + 1/14 <= Vmin*t  and  Vmax*t <= K + 13/14   for all t∈J,
       i.e. J ⊆ ( (K+1/14)/Vmin , (K+13/14)/Vmax ).  Call this the cluster window
       W_K := ( (K+1/14)/Vmin , (K+13/14)/Vmax ), nonempty iff
            (K+13/14)*Vmin > (K+1/14)*Vmax   <=>   13*Vmin - Vmax > 14*K*s.        (*)
  (C2) the small part P has a safe point in J.
Then M(S) >= 1/14.   [This is LEMMA 2 / band-fit, with the band = the WHOLE cluster
collapsed into the single gap index K; J's midpoint is a global witness.]

The cluster window W_K is a single-gap band fit for the ENTIRE cluster at gap index K.
Its width |W_K| = (K+13/14)/Vmax - (K+1/14)/Vmin.  The SMALL part is safe on a
fraction ~ mu_P of the circle; we need W_K to contain a small-safe point.

================================ THIS SCRIPT ==================================
(1) For every S3 set, enumerate ALL cluster windows W_K (K=0,1,2,... while nonempty by
    (*)), and check whether ANY of them contains a small-part-safe sub-interval. Report
    the success rate of this PURELY single-gap (provable, no widest-vs-average) lemma.
(2) For the sets it CLOSES, that is a rigorous proof. For the rest, report what's left.
(3) Also test the union: cluster windows W_K exist for K=0..Kmax with
    Kmax = floor((13 Vmin - Vmax)/(14 s)). Count them; the small part safe set has
    measure mu_P, and the windows are spaced ~1/Vmin apart -- so a pigeonhole over
    windows vs small-danger gives a provable count.
"""
from fractions import Fraction as F
from math import gcd, floor
from functools import reduce
from itertools import combinations
import random

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
    m = merge(iv)
    safe = []; prev = F(0)
    for a, b in m:
        if a > prev: safe.append((prev, a))
        prev = max(prev, b)
    if prev < 1: safe.append((prev, F(1)))
    return safe

def interval_has_small_safe(J, safe_arcs):
    """does open interval J=(lo,hi) intersect any small safe arc? return witness sub-int."""
    lo, hi = J
    for (a, b) in safe_arcs:
        aa = max(a, lo); bb = min(b, hi)
        if aa < bb: return (aa, bb)
    return None

def cluster_windows(Vmin, Vmax, s):
    """yield (K, lo, hi) for each nonempty cluster window W_K (single shared gap index K).
    W_K = ((K+1/14)/Vmin, (K+13/14)/Vmax). Nonempty iff 13 Vmin - Vmax > 14 K s.
    Note if s==0 (single value) all K up to ... but cluster has >=2 distinct, s>=1 here
    unless duplicate -- our S3 clusters have distinct values so s>=1; but if c could be
    such that Vmin==Vmax impossible. For the bound we also allow the per-runner version."""
    K = 0
    while True:
        if not (13 * Vmin - Vmax > 14 * K * s if s > 0 else True):
            # for s=0 window nonempty iff 13Vmin>Vmax i.e. (*) with s=0 always-ish; cap K
            if s == 0:
                if K > Vmin: break
            else:
                break
        lo = F(14 * K + 1, 14) / Vmin
        hi = F(14 * K + 13, 14) / Vmax
        if lo < hi and lo < 1:
            yield (K, lo, hi)
        if hi >= 1: break
        K += 1
        if K > 14 * Vmax: break

def gcd_all(S): return reduce(gcd, S, 0)
def is_covering(S): return all(any(v % q == 0 for v in S) for q in range(2, 15))
def classify(S):
    S = sorted(set(S)); Vmin, Vmax = S[0], S[-1]
    k = sum(1 for v in S if v > 13)
    if k <= 1: return 'S1'
    if Vmax < 13 * Vmin: return 'S2'
    return 'S3'
def missing_q(P): return [q for q in range(2, 15) if not any(v % q == 0 for v in P)]

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
        spread = rng.choice([2, 3, 5, 7, 10, 14, 20, 28, 40, 56])
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
    allS = []
    for (lo, hi) in [(50, 200), (200, 800), (800, 3000)]:
        allS += gen_constructive(seed=lo+1, target=1500, Vrange=(lo, hi))
    print(f"n S3 = {len(allS)}")
    closed = 0; not_closed = []
    nwindows = []
    cluster_window_nonempty = 0
    for S in allS:
        P = [u for u in S if u <= 13]; L = [u for u in S if u > 13]
        Vmin, Vmax = min(L), max(L); s = Vmax - Vmin
        safe = small_safe_arcs(P)
        wins = list(cluster_windows(Vmin, Vmax, s))
        nwindows.append(len(wins))
        if wins: cluster_window_nonempty += 1
        ok = False
        for (K, a, b) in wins:
            w = interval_has_small_safe((a, b), safe)
            if w is not None: ok = True; break
        if ok: closed += 1
        else: not_closed.append((S, len(wins), Vmin, Vmax, s, c if False else len(L)))
    n = len(allS)
    print(f"\nSINGLE-GAP cluster-collapse lemma closes {closed}/{n} S3 sets "
          f"({100*closed/n:.1f}%)  -- these are RIGOROUSLY proved M>=1/14.")
    print(f"cluster windows nonempty (at least one K): {cluster_window_nonempty}/{n}")
    print(f"# cluster windows W_K: min={min(nwindows)} mean={sum(nwindows)/len(nwindows):.2f} max={max(nwindows)}")
    print(f"\nNOT closed by single-gap lemma: {len(not_closed)} sets")
    for (S, nw, vmin, vmax, s, c) in not_closed[:12]:
        print(f"   S={S}  #wins={nw} Vmin={vmin} Vmax={vmax} s={s} c={c}")
