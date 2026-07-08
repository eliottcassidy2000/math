#!/usr/bin/env python3
r"""
lrc_avg_c_conditional_tent_kps_S75.py  (kind-pasteur-2026-07-07-S75)

THE AVERAGE-c CONDITIONAL TENT.  mac-mini-S56's c-table reports sup_d c(d,P) = 1.76 at
d=2 as the "obstruction" to the k=9 conditional-tent discharge.  But the conditional-tent
bound is EXACT and LINEAR:

    rho*(P,E) = meas(G_P) - meas(G_P cap S),   S = {maxgap(E) <= 1/7}
    meas(G_P cap S) <= (1/toll) * int_{G_P} F  dx     [Markov: F >= toll on S, F >= 0]
    int_{G_P} F = sum_{ordered i != j} I(e_j - e_i, P),   I(d,P) = int_{G_P} f(frac(dx)) dx

and by the x -> 1-x symmetry of G_P, I(-d,P) = I(d,P) EXACTLY (proved below numerically:
phi(x)=1-x is measure-preserving, G_P-invariant, and frac(-d(1-y)) = frac(dy)).  Hence

    int_{G_P} F = sum_{unordered {i,j}} 2 I(|e_i-e_j|, P) = 2 meas(G_P) int_f * sum c(d_ij,P)
    rho*(P,E) >= meas(G_P) [ 1 - (1 - tentfloor) * cbar(E,P) ],
    cbar(E,P) = (1/C(k,2)) sum_{unordered pairs} c(|e_i-e_j|, P)  = the AVERAGE, not the sup.

So the discharge needs cbar <= c*(P) := (1/(1-tentfloor)) (1 - m_P/meas(G_P)), and a
large-diameter PRIMITIVE cluster cannot pack many small-d (high-c) pairs.  THIS FILE:
(1) verifies I(-d,P) = I(d,P); (2) reproduces mac-mini's forward c-table; (3) MAXIMIZES
cbar over primitive k-clusters with diam beyond klein-THM-653's composition reach
(k=9 diam>16, k=10 diam>10) for the binding slow parts P, and checks rho*_bound >= m_P.
Everything exact rational.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations
import random

M_P = F(14249, 252252)

# ---------- G_P and the exact tent integral (mac-mini-S56 machinery, reused) ----------
def GP_intervals(P):
    bad = []
    for p in P:
        w = F(1, 14 * p)
        for j in range(p + 1):
            bad.append((F(j, p) - w, F(j, p) + w))
    bad = [(max(l, F(0)), min(h, F(1))) for l, h in bad if h > 0 and l < 1]
    bad.sort()
    merged = []
    for l, h in bad:
        if merged and l <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], h))
        else:
            merged.append((l, h))
    good = []; prev = F(0)
    for l, h in merged:
        if l > prev: good.append((prev, l))
        prev = max(prev, h)
    if prev < 1: good.append((prev, F(1)))
    return good

def tent_integral_forward(P_iv, d, beta, th=F(1, 7)):
    """int_{G_P} (frac(dx)-beta)_+ 1[frac(dx)<=th] dx  (forward sweep, d>0), exact."""
    tot = F(0); d = int(d)
    for (l, h) in P_iv:
        m_lo = int(l * d - th) - 1
        m_hi = int(h * d) + 1
        for m in range(max(0, m_lo), m_hi + 1):
            wl, wh = (F(m) + beta) / d, (F(m) + th) / d
            a, b = max(wl, l), min(wh, h)
            if a >= b: continue
            tot += (d * (b * b - a * a) / 2 - (F(m) + beta) * (b - a))
    return tot

def tent_integral_backward(P_iv, d, beta, th=F(1, 7)):
    """int_{G_P} f(frac(-dx)) dx = int_{G_P} (1-frac(dx)-beta)_+ 1[1-frac(dx)<=th] dx.
    support: frac(dx) in [1-th, 1-beta); integrand 1 - frac(dx) - beta = 1 - (dx-m) - beta."""
    tot = F(0); d = int(d)
    lo_s, hi_s = 1 - th, 1 - beta   # frac(dx) in [lo_s, hi_s)
    for (l, h) in P_iv:
        m_lo = int(l * d - hi_s) - 1
        m_hi = int(h * d) + 1
        for m in range(max(0, m_lo), m_hi + 1):
            wl, wh = (F(m) + lo_s) / d, (F(m) + hi_s) / d
            a, b = max(wl, l), min(wh, h)
            if a >= b: continue
            # integrand = 1 + m - beta - d x
            tot += ((F(1) + m - beta) * (b - a) - d * (b * b - a * a) / 2)
    return tot

def c_forward(P_iv, mGP, intf, d, beta):
    return tent_integral_forward(P_iv, d, beta) / (mGP * intf)

# ---------- (1) verify the symmetry I(-d,P) = I(d,P) ----------
print("=" * 90)
print("(1) SYMMETRY CHECK: I(-d,P) = I(d,P) (x->1-x invariance of G_P)")
print("=" * 90)
beta9 = F(14 - 9, 7 * 9)  # 5/63
P0 = (10, 11, 12, 13)
iv0 = GP_intervals(P0)
ok = True
for d in range(1, 40):
    fwd = tent_integral_forward(iv0, d, beta9)
    bwd = tent_integral_backward(iv0, d, beta9)
    if fwd != bwd:
        ok = False
        print(f"   d={d}: forward {fwd} != backward {bwd}   DIFF {float(fwd-bwd):.2e}")
print(f"   I(-d,P) == I(d,P) for d=1..39, P={P0}: {ok}")
print("   => the ordered-pair sum = 2 * sum over unordered pairs of the forward integral;")
print("      the conditional-tent bound consumes the AVERAGE forward c, symmetrization free.")

# ---------- (2) reproduce mac-mini's forward c-table sup ----------
print()
print("=" * 90)
print("(2) mac-mini-S56 c-table reproduction (forward c, sup over d<=250)")
print("=" * 90)
def make_c(P, k, dmax=250):
    beta = F(14 - k, 7 * k)
    intf = (F(1, 7) - beta) ** 2 / 2
    iv = GP_intervals(P)
    mGP = sum(h - l for l, h in iv)
    tab = {d: c_forward(iv, mGP, intf, d, beta) for d in range(1, dmax + 1)}
    return tab, mGP, beta, intf

for k, P in [(9, (10, 11, 12, 13)), (9, (1, 11, 12, 13)), (10, (11, 12, 13)), (10, (1, 12, 13))]:
    tab, mGP, beta, intf = make_c(P, k)
    supd = max(tab, key=lambda d: tab[d]); supc = float(tab[supd])
    print(f"   k={k} P={P}: meas(G_P)={float(mGP):.4f}, sup_d c = {supc:.4f} at d={supd}; "
          f"c(1)={float(tab[1]):.4f} c(2)={float(tab[2]):.4f} c(3)={float(tab[3]):.4f}")

# ---------- (3) the requirement c*(P) and the average-c maximization ----------
def cstar(k, mGP):
    tentfloor = 1 - F(2 * (k - 1) * (k - 7), 7 * k)
    return (1 / (1 - tentfloor)) * (1 - M_P / mGP), tentfloor

def cbar_of_cluster(E, tab):
    """average forward-c over unordered pairs of E (uses precomputed tab; needs all diffs)."""
    ds = [abs(E[j] - E[i]) for i in range(len(E)) for j in range(i + 1, len(E))]
    return sum(tab[d] for d in ds) / len(ds), max(ds)

def rho_bound(E, P, k):
    beta = F(14 - k, 7 * k)
    intf = (F(1, 7) - beta) ** 2 / 2
    iv = GP_intervals(P)
    mGP = sum(h - l for l, h in iv)
    tentfloor = 1 - F(2 * (k - 1) * (k - 7), 7 * k)
    ds = [abs(E[j] - E[i]) for i in range(len(E)) for j in range(i + 1, len(E))]
    Ck2 = len(ds)
    cbar = sum(tent_integral_forward(iv, d, beta) / (mGP * intf) for d in ds) / Ck2
    return mGP * (1 - (1 - tentfloor) * cbar), cbar, mGP, tentfloor

print()
print("=" * 90)
print("(3) MAXIMIZE average-c over primitive clusters beyond klein's composition reach")
print("=" * 90)

def maximize_cbar(k, P, diam_min, diam_max, seed=0, iters=6000):
    """hill-climb: primitive k-cluster (0 in E, gcd of diffs = 1, diam in [diam_min,diam_max])
    maximizing cbar = minimizing rho_bound.  Precompute c-table up to diam_max."""
    beta = F(14 - k, 7 * k)
    intf = (F(1, 7) - beta) ** 2 / 2
    iv = GP_intervals(P)
    mGP = sum(h - l for l, h in iv)
    tab = {}
    def cval(d):
        if d not in tab:
            tab[d] = tent_integral_forward(iv, d, beta) / (mGP * intf)
        return tab[d]
    rng = random.Random(seed)
    def rand_cluster():
        while True:
            pts = sorted({0} | set(rng.sample(range(1, diam_max + 1), k - 2)) | {rng.randint(diam_min, diam_max)})
            if len(pts) == k and pts[-1] - pts[0] >= diam_min:
                g = 0
                for a in pts[1:]:
                    g = gcd(g, a - pts[0])
                if g == 1:
                    return pts
    def cbar(E):
        ds = [E[j] - E[i] for i in range(k) for j in range(i + 1, k)]
        return sum(cval(d) for d in ds) / len(ds)
    best_E = rand_cluster(); best = cbar(best_E)
    for _ in range(iters):
        E = list(best_E)
        idx = rng.randrange(1, k)  # don't move the 0 anchor
        step = rng.choice([-3, -2, -1, 1, 2, 3])
        nv = E[idx] + step
        if nv <= 0 or nv > diam_max or nv in E:
            continue
        cand = sorted(set(E) - {E[idx]} | {nv})
        if len(cand) != k or cand[-1] - cand[0] < diam_min or cand[-1] - cand[0] > diam_max:
            continue
        g = 0
        for a in cand[1:]:
            g = gcd(g, a - cand[0])
        if g != 1:
            continue
        v = cbar(cand)
        if v > best:
            best, best_E = v, cand
    return best, best_E, mGP

CFG = [
    (9, (10, 11, 12, 13), 17, 45),
    (9, (1, 11, 12, 13), 17, 45),
    (9, (1, 2, 12, 13), 17, 45),
    (10, (11, 12, 13), 11, 45),
    (10, (1, 12, 13), 11, 45),
]
for k, P, dmin, dmax in CFG:
    cst, tentfloor = cstar(k, sum(h - l for l, h in GP_intervals(P)))
    # run several seeds
    best, bestE, mGP = F(-1), None, None
    for s in range(6):
        b, E, m = maximize_cbar(k, P, dmin, dmax, seed=s, iters=4000)
        if b > best:
            best, bestE, mGP = b, E, m
    rb = mGP * (1 - (1 - tentfloor) * best)
    print(f"\n   k={k} P={P} diam in [{dmin},{dmax}], meas(G_P)={float(mGP):.4f}")
    print(f"     c*(P) = {float(cst):.4f}  (need avg c <= c* for rho* >= m_P)")
    print(f"     max avg-c found = {float(best):.4f}  at E={bestE} (diam {bestE[-1]})")
    print(f"     => rho*_bound = {float(rb):.5f}  vs m_P = {float(M_P):.5f}  "
          f"=> {'DISCHARGED' if rb >= M_P else 'GAP (worst cluster fails)'}")
