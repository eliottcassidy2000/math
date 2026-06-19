#!/usr/bin/env python3
"""
LRC(14) -- Adversarial verification of the "gp-intersection-uniform" angle (kps-S5-wf).

Goal: rigorously re-derive every claimed EXACT rational, and HUNT for counterexamples
to the claimed floors. EXACT rationals only (fractions.Fraction).

Claims under test:
 (0) m_P := min_{|P|<=10} meas(G_P) = 14249/252252, with per-size minima listed.
 (1) via-max refutation: rho*_{2/7} = 0 EXACT on E={0,2,3,..,k}, k>=7, e.g.
     k=7 P={1,2,3,6,12,13} E={0,2,3,4,5,6,8}: meas(G_P)=515/1092, mu_{2/7}=13/35, rho*=0.
 (2) k<=7 pigeonhole: mu_{1/7}(E)=1 for ALL E (a.e.).
 (3) k>=8 union bound: rho*_{1/7} >= meas(G_P)+mu_{1/7}(E)-1; consecutive minimizes mu_{1/7};
     consecutive mu_{1/7}(k): k=8 691/735,... ; thr_k; union floor 1891/5880 (k=8 P={1,5,7,8,9}).
 (4) quasi-independence R_{1/7} >= 67053/84241.

G_P = {x in [0,1): ||p x|| >= 1/14 for all p in P}.  Danger of p = {x: ||p x|| < 1/14}.
maxgap is computed on the circle of the points {frac(e_i x)}.
"""

from fractions import Fraction as F
from itertools import combinations
import sys

# ---------------------------------------------------------------------------
# EXACT measure of a union of arcs given as list of (a,b) with 0<=a<b<=1.
# We represent sets on [0,1) as sorted disjoint (a,b) intervals (rationals).
# ---------------------------------------------------------------------------

def union_intervals(ivs):
    ivs = sorted([(a, b) for (a, b) in ivs if a < b])
    out = []
    for a, b in ivs:
        if out and a <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return out

def meas(ivs):
    return sum((b - a for a, b in union_intervals(ivs)), F(0))

def intersect_intervals(A, B):
    A = union_intervals(A); B = union_intervals(B)
    out = []
    i = j = 0
    while i < len(A) and j < len(B):
        a0, a1 = A[i]; b0, b1 = B[j]
        lo = max(a0, b0); hi = min(a1, b1)
        if lo < hi:
            out.append((lo, hi))
        if a1 < b1:
            i += 1
        else:
            j += 1
    return out

def complement(ivs):
    ivs = union_intervals(ivs)
    out = []
    prev = F(0)
    for a, b in ivs:
        if a > prev:
            out.append((prev, a))
        prev = max(prev, b)
    if prev < F(1):
        out.append((prev, F(1)))
    return out

# ---------------------------------------------------------------------------
# G_P : safe set of small speeds.  ||p x|| < 1/14 danger.
#   For integer p>=1: ||p x|| < 1/14 means dist(p x, Z) < 1/14.
#   On [0,1): p x in [k - 1/14, k + 1/14] for integer k.  x in [(k-1/14)/p,(k+1/14)/p].
# Build danger as union of arcs over k=0..p, then clip to [0,1) (with wrap handled by
# taking k from 0..p and intersecting with [0,1)).
# ---------------------------------------------------------------------------

def danger_p(p):
    th = F(1, 14)
    arcs = []
    # p x mod 1 < 1/14 or > 1 - 1/14.  Solve over x in [0,1):
    # p x in (k - th, k + th) for k = 0,1,...,p  (k=0 gives negative half from wrap)
    for k in range(0, p + 1):
        lo = (F(k) - th) / p
        hi = (F(k) + th) / p
        lo = max(lo, F(0)); hi = min(hi, F(1))
        if lo < hi:
            arcs.append((lo, hi))
    return union_intervals(arcs)

def safe_GP(P):
    if not P:
        return [(F(0), F(1))]
    dang = []
    for p in P:
        dang.extend(danger_p(p))
    dang = union_intervals(dang)
    return complement(dang)

# ---------------------------------------------------------------------------
# Good_E^theta = {x in [0,1): maxgap of points {frac(e_i x)} > theta}.
# Exact via breakpoints:
#   collisions: (e_i - e_j) x in Z  -> x = m/(e_i-e_j)
#   gap=theta crossings: for an ordered adjacent pair on the circle, the gap is a
#     piecewise-linear function of x; we add ALL candidate x where some difference
#     ((e_i - e_j) x mod 1) equals theta or 1-theta (these are the only places a gap
#     between two specific points can equal theta).  We over-generate breakpoints
#     (safe): any x where (e_i-e_j) x = integer +/- theta.
# Then on each open cell, maxgap is smooth & monotone-affine in the *active* adjacency;
# we sample the midpoint and ALSO require that no gap=theta crossing happens strictly
# inside (guaranteed because we added all such crossings).  Midpoint decision is exact.
# ---------------------------------------------------------------------------

def good_E(E, theta):
    E = sorted(set(E))
    k = len(E)
    th = theta
    bps = set([F(0), F(1)])
    diffs = set()
    for i in range(k):
        for j in range(i + 1, k):
            d = E[j] - E[i]
            if d != 0:
                diffs.add(d)
    for d in diffs:
        # collisions: x = m/d
        for m in range(0, d + 1):
            x = F(m, d)
            if 0 <= x <= 1:
                bps.add(x)
        # gap=theta crossings: (d x mod 1) = th or = 1-th
        # d x = m + th  -> x=(m+th)/d ; d x = m + (1-th) -> x=(m+1-th)/d
        m = 0
        while True:
            cands = [(F(m) + th) / d, (F(m) + (1 - th)) / d]
            done = True
            for x in cands:
                if 0 <= x <= 1:
                    bps.add(x)
                    done = False
                if x <= 1:
                    done = False
            if (F(m)) / d > 1:
                break
            m += 1
            if m > d + 2:
                break
    bps = sorted(b for b in bps if F(0) <= b <= F(1))
    good = []
    for a, b in zip(bps, bps[1:]):
        if a >= b:
            continue
        mid = (a + b) / 2
        pts = sorted(set((e * mid) % 1 for e in E))
        if len(pts) == 1:
            mg = F(1)
        else:
            gaps = [pts[t + 1] - pts[t] for t in range(len(pts) - 1)]
            gaps.append(pts[0] + 1 - pts[-1])
            mg = max(gaps)
        if mg > th:
            good.append((a, b))
    return union_intervals(good)

def mu_exact(E, theta):
    return meas(good_E(E, theta))

def rho_star(P, E, theta):
    GP = safe_GP(P)
    GE = good_E(E, theta)
    return meas(intersect_intervals(GP, GE))

# ---------------------------------------------------------------------------
# Independent cross-check of good_E via fine rational sampling refinement:
# (sanity, not used for final exact numbers) -- skip heavy.
# ---------------------------------------------------------------------------

def report(label, val):
    print(f"{label}: {val} = {float(val):.6f}")

# ===========================================================================
if __name__ == "__main__":
    out = []
    def P(*a):
        s = " ".join(str(x) for x in a)
        print(s); out.append(s)

    P("="*70)
    P("CLAIM (0): m_P = min_{|P|<=10} meas(G_P), per-size minima")
    P("="*70)
    claimed = {0:F(1), 1:F(6,7), 2:F(66,91), 3:F(55,91), 4:F(1979,4004),
               5:F(2243,5880), 6:F(3029,10780), 7:F(45107,229320),
               8:F(2479,17640), 9:F(10601,114660), 10:F(14249,252252)}
    speeds = list(range(1,14))
    per_size_min = {}
    per_size_arg = {}
    for sz in range(0,11):
        best = None; bestP=None
        for Pset in combinations(speeds, sz):
            m = meas(safe_GP(list(Pset)))
            if best is None or m < best:
                best = m; bestP = Pset
        per_size_min[sz]=best; per_size_arg[sz]=bestP
        match = "OK" if best==claimed[sz] else f"MISMATCH claimed {claimed[sz]}"
        P(f"  |P|={sz}: min meas(G_P)={best} ({float(best):.6f}) at {bestP}  [{match}]")
    mP = per_size_min[10]
    P(f"  => m_P = {mP} = {float(mP):.6f}  (claimed 14249/252252={F(14249,252252)})")
    P(f"  match m_P: {mP==F(14249,252252)}")
    P(f"  argmin |P|=10: {per_size_arg[10]} (claimed (1,2,3,5,7,8,9,11,12,13))")

    P("")
    P("="*70)
    P("CLAIM (1): via-max refutation rho*_{2/7}=0 EXACT")
    P("="*70)
    E7=[0,2,3,4,5,6,8]; P7=[1,2,3,6,12,13]
    gE7 = meas(safe_GP(P7)); muE7 = mu_exact(E7, F(2,7)); rE7 = rho_star(P7,E7,F(2,7))
    P(f"  k=7 P={P7} E={E7}")
    P(f"    meas(G_P)={gE7} (claimed 515/1092={F(515,1092)}) match={gE7==F(515,1092)}")
    P(f"    mu_2/7(E)={muE7} (claimed 13/35={F(13,35)}) match={muE7==F(13,35)}")
    P(f"    rho*_2/7={rE7} (claimed 0) match={rE7==F(0)}")
    # other witnesses
    for (kk,Pw,Ew) in [
        (8,[1,2,3,12,13],[0,2,3,4,5,6,7,9]),
        (9,[1,2,3,13],[0,2,3,4,5,6,7,8,10]),
        (10,[1,2,3],[0,2,3,4,5,6,7,8,9,11]),
    ]:
        r = rho_star(Pw,Ew,F(2,7))
        P(f"  k={kk} P={Pw} E={Ew}: rho*_2/7={r} match0={r==F(0)}")

    P("")
    P("="*70)
    P("CLAIM (2): k<=7 pigeonhole mu_{1/7}(E)=1 (a.e.)  -- EXHAUSTIVE bounded spread")
    P("="*70)
    th17 = F(1,7)
    for k in [3,4,5,6,7]:
        # exhaustive over shapes 0 in E, spread<=k+ (bound). spread cap chosen modest.
        cap = k+6
        worst = None; worstE=None
        cnt=0
        for rest in combinations(range(1,cap+1), k-1):
            E=[0]+list(rest)
            m = mu_exact(E, th17)
            cnt+=1
            if worst is None or m<worst:
                worst=m; worstE=E
        P(f"  k={k}: min mu_1/7 over {cnt} shapes (spread<={cap}) = {worst} at {worstE}")

    P("")
    P("="*70)
    P("CLAIM (3): k>=8 consecutive mu_{1/7} & union floor")
    P("="*70)
    cons_claim = {8:F(691,735),9:F(247,294),10:F(38,49),11:F(1381,2205),
                  12:F(13823,24255),13:F(477,1078)}
    thr_claim = {8:F(3637,5880),9:F(2025,4004),10:F(36,91),11:F(25,91),
                 12:F(1,7),13:F(0)}
    for k in range(8,14):
        E=list(range(0,k))  # consecutive 0..k-1
        m = mu_exact(E, th17)
        cm = cons_claim[k]
        P(f"  k={k} consecutive E={E}: mu_1/7={m} (claimed {cm}) match={m==cm}")
    # thr_k from per-size meas(G_P): for k, |P|=13-k, thr_k = 1 - min_{|P|=13-k} meas(G_P)
    P("  thr_k = 1 - min_{|P|=13-k} meas(G_P):")
    for k in range(8,14):
        sz = 13-k
        thr = 1 - per_size_min[sz]
        ck = thr_claim[k]
        P(f"    k={k} |P|={sz}: thr={thr} (claimed {ck}) match={thr==ck}; consec mu>=thr? {mu_exact(list(range(k)),th17) >= thr}")
    # union floor: min over k>=8, P of meas(G_P)+mu_1/7(consec_k)-1
    P("  union floor min[meas(G_P)+mu_1/7(consec_k)-1]:")
    best=None; bestinfo=None
    for k in range(8,14):
        sz=13-k
        muc = mu_exact(list(range(k)), th17)
        for Pset in combinations(speeds, sz):
            v = meas(safe_GP(list(Pset))) + muc - 1
            if best is None or v<best:
                best=v; bestinfo=(k,Pset)
    P(f"    min = {best} = {float(best):.6f} at k={bestinfo[0]} P={bestinfo[1]} (claimed 1891/5880={F(1891,5880)}) match={best==F(1891,5880)}")

    P("")
    P("="*70)
    P("ADVERSARIAL HUNT: mu_{1/7}(E) < thr_k over NON-consecutive shapes (k>=8)")
    P("="*70)
    # exhaustive bounded spread for k=8,9; broad for k>=10
    viol=0
    for k in [8,9]:
        sz=13-k; thr = 1 - per_size_min[sz]
        cap=k+4
        worst=None; worstE=None; cnt=0
        for rest in combinations(range(1,cap+1), k-1):
            E=[0]+list(rest)
            m=mu_exact(E,th17); cnt+=1
            if worst is None or m<worst:
                worst=m; worstE=E
            if m < thr:
                viol+=1
                if viol<=5:
                    P(f"    VIOLATION k={k} E={E} mu={m} < thr={thr}")
        P(f"  k={k}: exhaustive {cnt} shapes spread<={cap}; min mu_1/7={worst} at {worstE}; thr={thr}; min>=thr? {worst>=thr}")
    # k=10..13 broad structured + larger spread sweeps
    import random
    random.seed(7)
    for k in range(10,14):
        sz=13-k; thr = 1 - per_size_min[sz]
        worst=None; worstE=None; cnt=0
        # consecutive + perforated-AP + common-factor + random
        shapes=[]
        shapes.append(list(range(k)))
        for d in [2,3,5,7]:
            shapes.append([0]+[d*i for i in range(1,k)])  # AP gcd d (will scale)
        for drop in range(1,k):
            s=[i for i in range(k+1) if i!=drop]; shapes.append(s)  # perforated
        for _ in range(2000):
            cap=random.choice([k+2,k+5,k+10,2*k,4*k])
            rest=random.sample(range(1,cap+1), k-1)
            shapes.append([0]+sorted(rest))
        for E in shapes:
            m=mu_exact(E,th17); cnt+=1
            if worst is None or m<worst:
                worst=m; worstE=E
            if m<thr:
                viol+=1
                if viol<=10:
                    P(f"    VIOLATION k={k} E={E} mu={m} < thr={thr}")
        P(f"  k={k}: {cnt} shapes (structured+random); min mu_1/7={worst} at {worstE}; thr={thr}; min>=thr? {worst>=thr}")
    P(f"  TOTAL thr violations found: {viol}")

    P("")
    P("="*70)
    P("ADVERSARIAL HUNT (pure mu floor at 1/7 and 2/7): any mu=0?")
    P("="*70)
    # at 1/7, k<=7 must be 1; hunt mu_1/7=0 for any k (should be impossible since >0 per shape)
    found0=0
    for k in range(3,14):
        for _ in range(300):
            cap=random.choice([k+3,2*k,5*k,40])
            rest=random.sample(range(1,cap+1), min(k-1,cap))
            if len(rest)<k-1: continue
            E=[0]+sorted(rest)
            if mu_exact(E,th17)==0:
                found0+=1; P(f"   mu_1/7=0 at k={k} E={E}")
    P(f"  mu_1/7=0 occurrences: {found0}")

    P("")
    P("="*70)
    P("CLAIM (4): quasi-independence R_1/7 over consecutive E (k=8..13)")
    P("="*70)
    Rmin=None; Rarg=None
    for k in range(8,14):
        sz=13-k
        E=list(range(k)); muc=mu_exact(E,th17)
        for Pset in combinations(speeds,sz):
            gp=meas(safe_GP(list(Pset)))
            rho=rho_star(list(Pset),E,th17)
            denom=gp*muc
            if denom>0:
                R=rho/denom
                if Rmin is None or R<Rmin:
                    Rmin=R; Rarg=(k,Pset)
    P(f"  min R_1/7 = {Rmin} = {float(Rmin):.6f} at {Rarg} (claimed >=67053/84241={float(F(67053,84241)):.6f})")
    P(f"  Rmin >= 67053/84241 ? {Rmin >= F(67053,84241)}")

    P("")
    P("DONE.")
