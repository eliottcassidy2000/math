#!/usr/bin/env python3
r"""
LRC(14) residual case S3 -- ANGLE "two-gap-lemma2plus"  (CONSOLIDATED).

Generalize LEMMA 2 (single-gap band fit) to the genuine S3 mixed sets where the cluster
band straddles several gap indices. This file consolidates the full investigation and
records the PROVED results, the STRUCTURAL theorem, and the honest residual.

================================ SUMMARY OF RESULTS =============================

S = P (small part, subset of {1..13}) ∪ L (cluster, all speeds > 13), Vmin=min L,
    Vmax=max L=max S, s=Vmax-Vmin spread, c=|L|.

[STRUCTURAL THEOREM, PROVED]  (multi-band monotone split)
  At any t>0, the gap-index map u -> floor(u t) is monotone non-decreasing in u. So at a
  global witness tau*, the cluster L partitions into CONTIGUOUS speed sub-bands
  L = L_1 ⊔ ... ⊔ L_r, each landing in a distinct gap index K_1<...<K_r. Verified 2445/2445.
  Empirically r in {2,...,6}; r <= c always.

[LEMMA A, PROVED -- cluster-collapse single-gap (r=1)]  (generalizes LEMMA 2 to whole cluster)
  For integer K, the cluster window W_K = ((K+1/14)/Vmin, (K+13/14)/Vmax) is level-1/14
  SAFE for EVERY speed u in [Vmin,Vmax] (hence all of L). It is nonempty iff
  13 Vmin - Vmax > 14 K s. If some W_K contains a small-part-safe point, M(S) >= 1/14.
  Both claims verified EXACTLY (200000 / 1256 cases, 0 violations). CLOSES ~67% of S3.

[LEMMA B, PROVED -- short-arc tooth count]
  In a small-part safe arc I_P=(a,b) of width W_P, runner u plants <= floor(W_P u)+1 teeth;
  total cluster danger <= sum_{u in L}(floor(W_P u)+1)/(7u); largest free gap
  >= (W_P - danger_ub)/(N+1), N=sum(floor(W_P u)+1). If this exceeds 1/(7 Vmax), M(S)>=1/14.
  RIGOROUS. Adds a few % (the average-gap bound is weak for big clusters).

[GLOBAL-WITNESS EXISTENCE, VERIFIED (not yet PROVED in general)]
  Every S3 set has a single tau* safe for ALL of S (M(S)>=1/14 outright). Verified on
  thousands incl. adversarial; min margin G*·7·Vmax ~ 1.28 (>1), robust as Vmax->inf.

[LOOSENESS, VERIFIED]  min exact M over residual S3 (not Lemma-A-closable) = 2/19 ~ 0.105
  = 1.47x threshold. Residual is non-dangerous (consistent with THM-526's 4/31).

[HONEST NEGATIVE]  The crude measure bound mu_P > meas(cluster danger) holds 0/800: the
  cluster danger (<= c/7) exceeds small-part safe measure mu_P; positivity of the global
  safe set is due to ALIGNMENT (cluster teeth overlap small teeth), invisible to a union
  bound. The general multi-band (r>=2) existence needs the widest-vs-average lever, which
  has NO closed form here. The residual is INFINITE (Vmax unbounded), so Lemmas A+B do
  NOT finitize S3. A complete proof of S3 via this elementary angle is NOT achievable.

This script re-runs the decisive checks. Other scripts in the family hold the detailed
probes (anatomy, mechanism, multibandfit, toothcount, infimum, lemma_verify).
"""
from fractions import Fraction as F
from math import gcd, floor
from functools import reduce
from itertools import combinations
import random
from collections import Counter

# ---- tools ----
def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r
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
def short_arc_closable(S):
    P = [u for u in S if u <= 13]; L = [u for u in S if u > 13]
    Vmax = max(L)
    for (a, b) in small_safe_arcs(P):
        WP = b - a
        N = sum(floor(WP * u) + 1 for u in L)
        danger_ub = sum(F(floor(WP * u) + 1, 7 * u) for u in L)
        free_lb = WP - danger_ub
        if free_lb <= 0: continue
        if free_lb / (N + 1) > F(1, 7 * Vmax): return True
    return False
def gen_constructive(seed=0, target=2000, Vrange=(50, 2000), spreads=None):
    rng = random.Random(seed); out = []; tries = 0
    smalls = []
    base = list(range(1, 14))
    for sz in (11, 10, 9, 8, 7):
        for P in combinations(base, sz): smalls.append(list(P))
    if spreads is None: spreads = [2,3,5,7,10,14,20,28,40,56,70,90]
    while len(out) < target and tries < target * 300:
        tries += 1
        P = rng.choice(smalls); c = 13 - len(P)
        if c < 2: continue
        miss = missing_q(P); V = rng.randint(*Vrange)
        spread = rng.choice(spreads)
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
    print("="*78)
    print("LRC(14) S3 residual -- two-gap-lemma2plus consolidated checks")
    print("="*78)

    # build a wide S3 sample
    allS = []
    for (lo, hi) in [(50, 200), (200, 800), (800, 3000), (3000, 12000)]:
        allS += gen_constructive(seed=lo+1, target=1300, Vrange=(lo, hi))
    n = len(allS)
    print(f"\nS3 sample size: {n}  (Vmax up to ~12000)")

    # rigorous lemma coverage
    sgA = sum(1 for S in allS if single_gap_closable(S))
    sbB = sum(1 for S in allS if short_arc_closable(S))
    either = sum(1 for S in allS if single_gap_closable(S) or short_arc_closable(S))
    print(f"\nRIGOROUS coverage (genuine proofs of M(S)>=1/14):")
    print(f"  LEMMA A (cluster-collapse single-gap): {sgA}/{n} = {100*sgA/n:.1f}%")
    print(f"  LEMMA B (short-arc tooth count):       {sbB}/{n} = {100*sbB/n:.1f}%")
    print(f"  A or B:                                {either}/{n} = {100*either/n:.1f}%")

    # global witness existence + margin
    margins = []
    for S in allS:
        gw = find_global_witness(S)
        Vmax = max(u for u in S if u > 13)
        margins.append(float(gw[1] * 7 * Vmax) if gw else 0.0)
    print(f"\nGLOBAL WITNESS (VERIFIED, not proved in general):")
    print(f"  exists (margin>0) on {sum(1 for m in margins if m>0)}/{n}")
    print(f"  margin G*·7·Vmax: min={min(margins):.4f} mean={sum(margins)/len(margins):.3f}"
          f" (all>1: {min(margins)>1})")

    # multi-band structure on the genuine multi-band sets
    multi = [S for S in allS if not single_gap_closable(S)]
    rcount = Counter(); mono_ok = 0
    for S in multi:
        gw = find_global_witness(S)
        if not gw: continue
        tau = gw[0]; L = sorted(u for u in S if u > 13)
        Ks = [int(u * tau) for u in L]
        rcount[len(set(Ks))] += 1
        if all(Ks[i] <= Ks[i+1] for i in range(len(Ks)-1)): mono_ok += 1
    print(f"\nMULTI-BAND STRUCTURE (PROVED: monotone sub-band split):")
    print(f"  genuine multi-band sets (not single-gap closable): {len(multi)}")
    print(f"  gap-index monotone in speed: {mono_ok}/{len(multi)}")
    print(f"  r=#distinct gap indices: {dict(sorted(rcount.items()))}")
    print(f"  max r observed: {max(rcount) if rcount else 0}  (always <= c)")

    print("\nSee companion scripts for: anatomy, mechanism, multibandfit, toothcount,")
    print("infimum (min exact M=2/19~0.105 over residual), and lemma_verify (exact proofs).")
