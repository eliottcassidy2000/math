#!/usr/bin/env python3
"""
LRC(14) S3 -- CONSTRUCTIVE adversarial generator + exact infimum of the
global-witness margin.

Constructive S3 generator: choose small part P (subset of {1..13}), compute the q in
2..14 that P fails to cover, then build a tight cluster L (near V, spread s) that covers
those q. This reaches the genuinely tight S3 sets fast (the rejection sampler stalls
because covering 14 with a 2-element tight cluster is rare).

We compute the global-witness margin G*·7·Vmax exactly and report the infimum + the
worst sets, sweeping V upward to test whether the margin stays > 1 as Vmax -> infinity.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations
import random

def merge(ivs):
    ivs = sorted(ivs); m = []
    for a, b in ivs:
        if m and a <= m[-1][1]: m[-1] = (m[-1][0], max(m[-1][1], b))
        else: m.append((a, b))
    return m

def teeth_in(u, lo, hi, h=F(1, 14)):
    out = []
    jmin = int(lo * u) - 1; jmax = int(hi * u) + 1
    for j in range(jmin, jmax + 1):
        c = F(j, u); a = c - h/u; b = c + h/u
        aa = max(a, lo); bb = min(b, hi)
        if aa < bb: out.append((aa, bb))
    return out

def safe_components(A, h=F(1, 14)):
    iv = []
    for u in A:
        for j in range(0, u):
            c = F(j, u); a = (c - h/u) % 1; b = (c + h/u) % 1
            if a < b: iv.append((a, b))
            else: iv.append((a, F(1))); iv.append((F(0), b))
    iv.sort(); merged = []
    for a, b in iv:
        if merged and a <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else: merged.append((a, b))
    safe = []; prev = F(0)
    for a, b in merged:
        if a > prev: safe.append((prev, a))
        prev = max(prev, b)
    if prev < 1: safe.append((prev, F(1)))
    return safe

def largest_free_gap(lo, hi, cluster):
    dgr = []
    for u in cluster: dgr.extend(teeth_in(u, lo, hi))
    dgr = merge(dgr)
    best = F(0); prev = lo
    for a, b in dgr:
        if a > prev and a - prev > best: best = a - prev
        prev = max(prev, b)
    if hi > prev and hi - prev > best: best = hi - prev
    return best

def global_witness_margin(S):
    small = [u for u in S if u <= 13]; cluster = [u for u in S if u > 13]
    Vmax = max(cluster)
    sc = safe_components(small)
    best = F(0)
    for (a, b) in sc:
        g = largest_free_gap(a, b, cluster)
        if g > best: best = g
    return best * 7 * Vmax, best

def gcd_all(S): return reduce(gcd, S, 0)
def is_covering(S): return all(any(v % q == 0 for v in S) for q in range(2, 15))
def classify(S):
    S = sorted(set(S)); Vmin, Vmax = S[0], S[-1]
    k = sum(1 for v in S if v > 13)
    if k <= 1: return 'S1'
    if Vmax < 13 * Vmin: return 'S2'
    return 'S3'

def missing_q(P):
    return [q for q in range(2, 15) if not any(v % q == 0 for v in P)]

def gen_constructive(seed=0, target=3000, Vrange=(50, 3000)):
    rng = random.Random(seed); out = []; tries = 0
    smalls = []
    # all small parts of size 7..11 from {1..13}
    base = list(range(1, 14))
    for sz in (11, 10, 9, 8, 7):
        for P in combinations(base, sz):
            smalls.append(list(P))
    while len(out) < target and tries < target * 300:
        tries += 1
        P = rng.choice(smalls)
        c = 13 - len(P)
        if c < 2: continue
        miss = missing_q(P)   # cluster must cover these
        V = rng.randint(*Vrange)
        spread = rng.choice([2, 3, 5, 7, 10, 14, 20, 28, 40, 56])
        # build cluster of size c near V covering miss
        # pick c values in [V, V+spread]; check covers miss
        window = list(range(V, V + spread + 1))
        if len(window) < c: continue
        cluster = rng.sample(window, c)
        if not all(any(v % q == 0 for v in cluster) for q in miss): continue
        S = sorted(set(P) | set(cluster))
        if len(S) != 13: continue
        if gcd_all(S) != 1: continue
        if not is_covering(S): continue
        if classify(S) != 'S3': continue
        out.append(S)
    return out

if __name__ == "__main__":
    allmargins = []
    for (lo, hi) in [(50, 200), (200, 800), (800, 3000), (3000, 12000)]:
        S3 = gen_constructive(seed=lo, target=1500, Vrange=(lo, hi))
        ms = [(global_witness_margin(S)[0], S) for S in S3]
        ms.sort(key=lambda x: x[0])
        allmargins += ms
        mn = float(ms[0][0]) if ms else None
        mean = float(sum(m for m, _ in ms)/len(ms)) if ms else None
        print(f"V in [{lo},{hi}): n={len(S3)}  min margin={mn:.4f}  mean={mean:.3f}  "
              f"all>1: {ms[0][0] > 1 if ms else '?'}")
    allmargins.sort(key=lambda x: x[0])
    print(f"\nOVERALL min global-witness margin = {float(allmargins[0][0]):.4f}")
    print(f"all > 1: {'YES' if allmargins[0][0] > 1 else 'NO -- FAILURE'}")
    print("\n12 worst S3 sets (smallest global-witness margin):")
    for m, S in allmargins[:12]:
        cl = [u for u in S if u > 13]; sm = [u for u in S if u <= 13]
        print(f"   margin={float(m):.4f}  small={sm} cluster={cl} "
              f"(c={len(cl)},s={max(cl)-min(cl)},Vmax={max(cl)})")
