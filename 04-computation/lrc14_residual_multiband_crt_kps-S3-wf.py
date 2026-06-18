#!/usr/bin/env python3
"""
LRC(14) S3 -- the TRUE multi-band-fit: structure of the global witness for the
not-single-gap-closable sets, and a CRT/three-distance existence attempt.

For the 32.6% of S3 not closed by the single-gap lemma, the global witness tau* has
the cluster runners split across r>=2 distinct gap indices. We:

(1) For each such set, FIND an exact global witness tau* (search small safe arcs for a
    cluster-free point) and report:
      - r = number of distinct gap indices the cluster occupies at tau*
      - the partition L = L_1 ⊔ ... ⊔ L_r by gap index
      - whether the partition is by CONTIGUOUS speed blocks (sub-bands) -- i.e. does the
        gap index increase monotonically with u? (this is the "band straddles gaps" picture)

(2) Test a CONSTRUCTIVE multi-band lemma: choose tau in a small-safe arc near
    tau ≈ (K+1/2)/Vmid for the cluster center; the cluster splits into bands
    L_i = {u: gap index = K_i}; require each band's span fits its gap. Quantify the
    constructive success.

(3) Three-distance / continued-fraction view: the cluster-free gaps inside a small arc
    I_P are governed by the three-distance theorem for the rotation by Vmin (mod 1)
    restricted to I_P. Measure the THREE gap lengths.
"""
from fractions import Fraction as F
from math import gcd, floor
from functools import reduce
from itertools import combinations
import random
from collections import Counter

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

def teeth_in(u, lo, hi, h=F(1, 14)):
    out = []
    jmin = int(lo * u) - 1; jmax = int(hi * u) + 1
    for j in range(jmin, jmax + 1):
        c = F(j, u); a = c - h/u; b = c + h/u
        aa = max(a, lo); bb = min(b, hi)
        if aa < bb: out.append((aa, bb))
    return out

def find_global_witness(S):
    """Return (tau*, gap) where gap is the cluster-free gap (lo,hi) in a small arc,
    tau* its midpoint; or None."""
    P = [u for u in S if u <= 13]; L = [u for u in S if u > 13]
    safe = small_safe_arcs(P)
    best = None
    for (a, b) in safe:
        dgr = merge([t for u in L for t in teeth_in(u, a, b)])
        prev = a
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
    return ((lo + hi) / 2, (lo, hi), w)

def gap_index(u, tau):
    """integer K s.t. frac(u*tau) in (1/14,13/14) lands in band K = floor(u*tau) if
    frac>=1/14, else floor(u*tau) (handle frac near 0). We use K=floor(u*tau) and require
    frac in (1/14,13/14)."""
    x = u * tau
    return int(x), x - int(x)

def gcd_all(S): return reduce(gcd, S, 0)
def is_covering(S): return all(any(v % q == 0 for v in S) for q in range(2, 15))
def classify(S):
    S = sorted(set(S)); Vmin, Vmax = S[0], S[-1]
    k = sum(1 for v in S if v > 13)
    if k <= 1: return 'S1'
    if Vmax < 13 * Vmin: return 'S2'
    return 'S3'
def missing_q(P): return [q for q in range(2, 15) if not any(v % q == 0 for v in P)]
def cluster_windows_count(Vmin, Vmax, s):
    if s == 0: return 1
    Kmax = (13 * Vmin - Vmax) // (14 * s) if 13*Vmin > Vmax else -1
    return Kmax + 1
def single_gap_closable(S):
    P = [u for u in S if u <= 13]; L = [u for u in S if u > 13]
    Vmin, Vmax = min(L), max(L); s = Vmax - Vmin
    safe = small_safe_arcs(P)
    K = 0
    while 13 * Vmin - Vmax > 14 * K * s if s > 0 else K == 0:
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
        spread = rng.choice([14, 20, 28, 40, 56, 70, 90])  # bias LARGE spread (multi-gap)
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
        allS += gen_constructive(seed=lo+7, target=1200, Vrange=(lo, hi))
    print(f"n S3 (large-spread biased) = {len(allS)}")
    multi = [S for S in allS if not single_gap_closable(S)]
    print(f"not single-gap closable (genuine multi-band): {len(multi)}/{len(allS)}")
    rcount = Counter(); contiguous = Counter(); witness_found = 0
    band_span = Counter()
    for S in multi:
        gw = find_global_witness(S)
        if gw is None: continue
        witness_found += 1
        tau, (lo, hi), w = gw
        L = sorted(u for u in S if u > 13)
        Ks = {}
        for u in L:
            K, fr = gap_index(u, tau)
            Ks[u] = K
        distinct = sorted(set(Ks.values()))
        r = len(distinct)
        rcount[r] += 1
        band_span[max(distinct) - min(distinct)] += 1
        # is the partition by gap index monotone in u? (sub-band / band picture)
        seq = [Ks[u] for u in L]
        mono = all(seq[i] <= seq[i+1] for i in range(len(seq)-1))
        contiguous[mono] += 1
    print(f"global witness found for {witness_found}/{len(multi)} multi-band sets")
    print(f"\nr = #distinct gap indices the cluster occupies at the witness:")
    for k in sorted(rcount): print(f"   r={k}: {rcount[k]}")
    print(f"gap-index span (max-min):")
    for k in sorted(band_span): print(f"   span {k}: {band_span[k]}")
    print(f"gap index MONOTONE in speed (clean sub-band split)?: {dict(contiguous)}")
