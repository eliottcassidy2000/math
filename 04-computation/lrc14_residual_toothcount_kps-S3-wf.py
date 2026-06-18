#!/usr/bin/env python3
"""
LRC(14) S3 -- RIGOROUS tooth-counting lower bound on the global-witness gap.

Setup: P = small part, safe arc I_P=(alpha,beta), width W_P. Cluster L={u_1<...<u_c},
all in [Vmin,Vmax], Vmax=max S.

A single runner u plants teeth centered at j/u, half-width 1/(14u). The number of
tooth CENTERS j/u inside an interval of length W_P is at most floor(W_P*u)+1.
So the number of teeth of u meeting I_P is N_u <= floor(W_P*u)+1 <= W_P*u + 1.
Total cluster teeth meeting I_P: N <= sum_u (W_P*u+1) = W_P*(sum_L u) + c.
After merging, #merged-danger-components <= N, so #free-gaps <= N+1.
Total cluster danger measure in I_P <= sum_u N_u*(1/(7u)) <= sum_u (W_P*u+1)/(7u)
    = sum_u (W_P/7 + 1/(7u)) = c*W_P/7 + (1/7)*sum_u 1/u.
Free measure in I_P >= W_P - [c*W_P/7 + (1/7)sum_L 1/u]
                     = W_P*(1 - c/7) - (1/7) sum_L 1/u.
Largest gap >= free_measure / (#gaps) >= free_measure/(N+1).

PROBLEM: 1-c/7 can be <= 0 for c>=7. And N+1 can be large. So this CRUDE bound is too
weak. We need a SHARPER count exploiting cluster TIGHTNESS: the c runners' teeth nearly
COINCIDE, so the merged danger has far fewer components than N.

SHARPER MODEL (cluster collapse): all u in L are within s of Vmin. On I_P, consider the
'cluster phase' theta(tau)=Vmin*tau. As tau ranges over I_P (width W_P), theta ranges over
width Vmin*W_P. Each cluster runner u=Vmin+d (0<=d<=s) has u*tau = theta + d*tau, so within
I_P its phase is theta + d*tau, and d*tau in [d*alpha, d*beta], a shift of width d*W_P<=s*W_P.
So all c runners' phase-combs are SHIFTS of the Vmin comb by <= s*W_P (in tau units: the
teeth of u and Vmin differ in position by <= (d/Vmin)* (local), bounded). If s*W_P < tooth
spacing... we measure the EFFECTIVE number of merged comb-passes.

THIS SCRIPT measures, in the worst regime, the merged-danger COMPONENT COUNT vs the crude
N bound, to see how much tightness helps, and tests a refined provable bound:
   Largest gap >= W_P/(K)  -  (tooth widths)
where K = number of merged cluster passes ~ ceil(Vmax*W_P)+something.
"""
from fractions import Fraction as F
from math import gcd, floor, ceil
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

def gap_stats(lo, hi, cluster):
    dgr = []
    for u in cluster: dgr.extend(teeth_in(u, lo, hi))
    dgr = merge(dgr)
    gaps = []; prev = lo
    for a, b in dgr:
        if a > prev: gaps.append(a - prev)
        prev = max(prev, b)
    if hi > prev: gaps.append(hi - prev)
    return (max(gaps) if gaps else F(0)), len(dgr), len(gaps), sum(gaps)

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
    S3 = gen_constructive(seed=99, target=2000, Vrange=(50, 3000))
    print(f"n S3 = {len(S3)}")
    # For each set, for its BEST small arc (the one achieving widest cluster-free gap),
    # report: W_P, Vmax, #merged danger comps n_merged, crude N bound, ratio.
    rows = []
    worst_margin = (F(10), None)
    n_le_widthV = 0
    for S in S3:
        small = [u for u in S if u <= 13]; cluster = [u for u in S if u > 13]
        Vmax = max(cluster); c = len(cluster); s = Vmax - min(cluster)
        sc = safe_components(small)
        best = None
        for (a, b) in sc:
            G, nmerged, ngaps, free = gap_stats(a, b, cluster)
            cand = (G, (b - a), nmerged, ngaps, free, a, b)
            if best is None or G > best[0]: best = cand
        G, WP, nmerged, ngaps, free, a, b = best
        margin = G * 7 * Vmax
        # crude tooth count: N <= W_P*sum(cluster)+c ; actual nmerged
        Ncrude = float(WP) * sum(cluster) + c
        # tightness-aware: number of distinct cluster-passes ~ ceil(Vmax*WP)+1 (merged)
        passes = ceil(float(Vmax * WP)) + 1
        rows.append((float(margin), c, s, float(Vmax*WP), nmerged, Ncrude, passes))
        if margin < worst_margin[0]: worst_margin = (margin, S, c, s, float(Vmax*WP), nmerged)
        # is merged comp count <= c*(ceil(Vmax*WP)+1)?  (each runner <= that many passes)
    # Test the sharp claim: n_merged <= c*(ceil(Vmax*WP)+1) and also n_merged <= ceil((Vmax)*WP)+c?
    viol1 = sum(1 for (m,c,s,vw,nm,Nc,p) in rows if nm > c*(ceil(vw)+1))
    viol2 = sum(1 for (m,c,s,vw,nm,Nc,p) in rows if nm > ceil(vw)+c)  # tightness: ~one comb + c offsets
    print(f"\nClaim A: n_merged <= c*(ceil(Vmax*W_P)+1) : violations {viol1}/{len(rows)}")
    print(f"Claim B: n_merged <= ceil(Vmax*W_P)+c     : violations {viol2}/{len(rows)}")
    margins = [r[0] for r in rows]
    print(f"\nglobal-witness margin: min={min(margins):.4f} mean={sum(margins)/len(margins):.3f}")
    m, S, c, s, vw, nm = worst_margin
    print(f"worst: margin={float(m):.4f}  S={S}")
    print(f"       c={c} s={s} Vmax*W_P={vw:.3f} n_merged={nm}")
    # Among worst 30, show Vmax*W_P and n_merged to see the tightness lever
    rows.sort(key=lambda x: x[0])
    print("\n worst 15: margin  c  s  Vmax*W_P  n_merged  (ceil(Vmax*WP)+c)")
    for (m,c,s,vw,nm,Nc,p) in rows[:15]:
        print(f"   {m:.3f}  c={c} s={s}  VW={vw:.2f}  nm={nm}  bnd={ceil(vw)+c}")
