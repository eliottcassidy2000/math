#!/usr/bin/env python3
"""
lrc14_Bk_zero_investigation_kps-S5-wf.py  (kind-pasteur-2026-06-18-S5)

CRITICAL: the admissible (|P|+k=13) scan found rho*(P,E)=0 for several configs, the
cleanest being P={1,2,3}, E=(0,1,2,3,4,5,6,7,8,10) (k=10) and P={1,2,3,4},E=(0,2,..,8,10).
These are ADMISSIBLE S3 instances -- so the SINGLE-GLOBAL-WITNESS route of THM-527 FAILS
on them (rho*=0). We investigate rigorously:

 (1) CONFIRM rho*=0 exactly, and show WHY: Good(E) and G_P are each nonempty but DISJOINT.
     Print Good(E), G_P, and their (empty) intersection.

 (2) IS THE PURE mu(E) ALSO 0?  No -- mu(E)>0 by L2 (per-shape). The vanishing is the
     INTERSECTION with G_P. So Good(E)'s mass sits OUTSIDE G_P (anti-correlation in this case).
     Report mu(E) and meas(G_P) and where Good(E)'s arcs lie relative to G_P's forbidden combs.

 (3) DOES rho*=0 imply M(S)<1/14?  NO. THM-527 says rho*>0 => M>=1/14 (sufficient, not nec.).
     We DIRECTLY estimate M(S) for the actual speed set S = P ∪ {large cluster realizing E}.
     A large cluster with co-offsets E and base Vmax: speeds = {Vmax - e : e in E}. We pick a
     concrete Vmax (e.g. 200) so all large speeds > 13 and distinct from P, then compute the
     loneliness M(S) = max_t min_v ||v t|| over a fine grid (float, high-res) to see if it is
     >= 1/14. This is EVIDENCE (not exact) that LRC(14) still HOLDS on these configs despite
     rho*=0 -- i.e. the global-witness route is incomplete, not that the conjecture fails.

 (4) ENUMERATE how many admissible (P,E) have rho*=0 in a FULLY EXHAUSTIVE small-case sweep
     (all P subset {1..13} with |P|=13-k, all primitive E spread<=cap) for k=10,11 -- to see
     whether the zero-locus is small/structured (near-consecutive E + tiny P) or pervasive.
"""
import sys, itertools, random
from fractions import Fraction as F
from math import gcd
sys.stdout.reconfigure(line_buffering=True)
random.seed(2468)
TWO7 = F(2, 7)

def merge(iv):
    iv = sorted(iv); out = []
    for a, b in iv:
        if out and a <= out[-1][1]: out[-1] = (out[-1][0], max(out[-1][1], b))
        else: out.append((a, b))
    return out
def meas(arcs): return sum((b - a for a, b in arcs), F(0))
def good_set_exact(E):
    E = sorted(set(E)); k = len(E)
    if k == 1: return [(F(0), F(1))]
    diffs = set()
    for a in range(k):
        for b in range(a + 1, k): diffs.add(E[b] - E[a])
    bps = {F(0), F(1)}
    for d in diffs:
        for m in range(0, d + 1): bps.add(F(m, d))
    bps = sorted(x for x in bps if 0 <= x <= 1)
    good = []
    for i in range(len(bps) - 1):
        x0, x1 = bps[i], bps[i + 1]
        if x1 <= x0: continue
        xm = (x0 + x1) / 2
        pts = sorted(((E[t] * xm) % 1, E[t]) for t in range(k))
        order = [e for _, e in pts]; floors = [int((e * xm) // 1) for e in order]
        for idx in range(k):
            e_cur = order[idx]; f_cur = floors[idx]
            if idx < k - 1: e_nx = order[idx + 1]; f_nx = floors[idx + 1]; wrap = F(0)
            else: e_nx = order[0]; f_nx = floors[0]; wrap = F(1)
            A = F(e_nx - e_cur); Cc = F(f_cur - f_nx) + wrap
            if A == 0:
                if Cc > TWO7: good.append((x0, x1))
                continue
            xb = (TWO7 - Cc) / A
            if A > 0: lo = max(x0, xb); hi = x1
            else: lo = x0; hi = min(x1, xb)
            if lo < hi: good.append((lo, hi))
    return merge(good)
def mu(E): return meas(good_set_exact(E))
def GP_arcs(P):
    safe = [(F(0), F(1))]
    for p in P:
        forb = []
        for j in range(p + 1):
            c = F(j, p); forb.append((c - F(1, 14 * p), c + F(1, 14 * p)))
        forb = merge([(max(F(0), a), min(F(1), b)) for a, b in forb if b > 0 and a < 1])
        sp = []; cur = F(0)
        for a, b in forb:
            if a > cur: sp.append((cur, a))
            cur = max(cur, b)
        if cur < 1: sp.append((cur, F(1)))
        new = []
        for a, b in safe:
            for c, d in sp:
                lo, hi = max(a, c), min(b, d)
                if lo < hi: new.append((lo, hi))
        safe = merge(new)
    return safe
def intersect(A, B):
    out = []
    for a, b in A:
        for c, d in B:
            lo, hi = max(a, c), min(b, d)
            if lo < hi: out.append((lo, hi))
    return merge(out)
def is_primitive(E):
    g = 0
    for e in E: g = gcd(g, e)
    return g == 1

print("="*90)
print("(1)-(2) WHY rho*=0 for the witness configs (Good(E), G_P each nonempty but DISJOINT)")
print("="*90)
cases = [
    ((1,2,3), (0,1,2,3,4,5,6,7,8,10)),
    ((1,2,3,4), (0,2,3,4,5,6,7,8,10)),
]
for P, E in cases:
    gE = good_set_exact(E); gp = GP_arcs(P); inter = intersect(gp, gE)
    print(f"\n  P={P}  E={E}  (|P|+k = {len(P)+len(E)})")
    print(f"     mu(E)        = {mu(E)} = {float(mu(E)):.6f}   (Good(E) nonempty, per L2)")
    print(f"     meas(G_P)    = {meas(gp)} = {float(meas(gp)):.6f}   (G_P nonempty)")
    print(f"     rho*=meas(cap)= {meas(inter)}   (intersection EMPTY -> single global witness FAILS)")
    print(f"     Good(E) arcs (first 6): {[ (str(a),str(b)) for a,b in gE[:6] ]}")
    print(f"     G_P arcs     (first 6): {[ (str(a),str(b)) for a,b in gp[:6] ]}")

print()
print("="*90)
print("(3) DOES rho*=0 imply M(S)<1/14?  NO (THM-527 is one-directional). Direct M(S) estimate.")
print("    Build speeds S = P u {Vmax - e : e in E}, Vmax large so cluster > 13, distinct.")
print("    M(S) = max_t min_{v in S} ||v t||  (fine float grid). LRC(14) wants M(S) >= 1/14.")
print("="*90)
def frac(x): return x - int(x) if x>=0 else x-int(x)+ (1 if x!=int(x) else 0)
def dist_half(y):
    y = y - int(y)
    if y < 0: y += 1
    return min(y, 1-y)
def M_estimate(S, N=2000000):
    best = 0.0
    Sl = list(S)
    for t in range(1, N):
        x = t / N
        m = 1.0
        for v in Sl:
            d = dist_half(v * x)
            if d < m:
                m = d
                if m <= best: break
        if m > best:
            best = m
    return best
for P, E in cases:
    for Vmax in (200, 503):
        cluster = sorted(Vmax - e for e in E)
        S = sorted(set(list(P) + cluster))
        if len(S) != len(P) + len(E):
            print(f"   (skip Vmax={Vmax}: overlap between P and cluster)"); continue
        Mest = M_estimate(S, N=1200000)
        oneok = Mest >= 1/14 - 1e-6
        print(f"   P={P} Vmax={Vmax}: |S|={len(S)}  M(S)~{Mest:.6f}  (1/14={1/14:.6f})  "
              f"-> {'M>=1/14 (LRC OK despite rho*=0)' if oneok else 'M<1/14 ?? CHECK'}")

print()
print("="*90)
print("(4) HOW STRUCTURED is the admissible rho*=0 locus?  Exhaustive small sweep k=10, |P|=3.")
print("    All P=3-subsets of {1..13}, all primitive E (0 in E,|E|=10) spread<=12: count zeros,")
print("    record whether zeros are exactly the near-consecutive E with small P.")
print("="*90)
k = 10; szP = 13 - k
allP = list(itertools.combinations(range(1,14), szP))
# primitive E spread<=12
Elist = []
for rest in itertools.combinations(range(1,13), k-1):
    E = (0,)+rest
    if is_primitive(E): Elist.append(E)
print(f"   #P={len(allP)}  #E(spread<=12)={len(Elist)}  -> {len(allP)*len(Elist)} pairs (exact)")
zero_pairs = []
total = 0
# precompute G_P arcs for each P
GPs = {P: GP_arcs(P) for P in allP}
for E in Elist:
    gE = good_set_exact(E)
    sprd = max(E)
    for P in allP:
        total += 1
        if meas(intersect(GPs[P], gE)) == 0:
            zero_pairs.append((P, E, sprd))
print(f"   scanned {total} admissible (P,E); zeros = {len(zero_pairs)}")
# characterize: are zero E all 'near-consecutive' (spread == k-1 or k, few gaps)?
from collections import Counter
spread_count = Counter(sprd for _,_,sprd in zero_pairs)
Pcount = Counter(P for P,_,_ in zero_pairs)
print(f"   zero-E spreads (k-1={k-1}): {dict(spread_count)}")
print(f"   zero-P multiset: {dict(Pcount)}")
for P,E,s in zero_pairs[:12]:
    print(f"      zero: P={P} E={E} spread={s}")
print("\nDONE.")
