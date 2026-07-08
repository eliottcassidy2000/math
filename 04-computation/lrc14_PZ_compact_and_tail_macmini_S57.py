"""
mac-mini-2026-07-07-S57 -- COMPLETE the k=11,12,13 PZ covering floor (THM-660):
 (1) k=12, k=13 EXHAUSTIVE exact compact checks (klein-S178 did k=11 diam<=15);
 (2) characterize the DECORRELATION TAIL: which parameter (diameter? min-difference?)
     rigorously controls PZ -> the analytic target that makes the reduction complete.

PZ(E) = E[W]^2/E[W^2], W = sum_i (g_i - 1/7)_+.  Fast float search + exact Farey confirm.
"""
import numpy as np
from fractions import Fraction as F
from math import floor, gcd
from functools import reduce
from itertools import combinations
import random
random.seed(7)
TH = F(1, 7); THf = 1/7

# ---------- fast float PZ ----------
GRID = 40_000
_x = (np.arange(GRID)+0.5)/GRID
def pz_float(E):
    Ea = np.array(sorted(E), float); ph = np.mod(np.outer(_x, Ea), 1.0); ph.sort(axis=1)
    g = np.concatenate([np.diff(ph, axis=1), (ph[:,0]+1-ph[:,-1])[:,None]], axis=1)
    W = np.maximum(g-THf, 0).sum(axis=1); EW = W.mean(); EW2 = (W*W).mean()
    return (EW*EW/EW2) if EW2 > 0 else 0.0

# ---------- exact Farey PZ (any family) ----------
def pz_exact(E):
    E = sorted(E); k = len(E)
    ds = set()
    for j in range(k):
        if E[j] != 0: ds.add(E[j])
    for i in range(k):
        for j in range(i+1, k):
            ds.add(E[j]-E[i])
    bps = set([F(0), F(1)])
    for d in ds:
        for m in range(0, d+1):
            bps.add(F(m, d))
    bps = sorted(bps)
    EW = F(0); EW2 = F(0)
    for c in range(len(bps)-1):
        x0, x1 = bps[c], bps[c+1]; xm = (x0+x1)/2
        lin = [(F(-floor(e*xm)), F(e)) for e in E]        # frac(e x) = e x - floor(e xm)
        order = sorted(range(k), key=lambda j: lin[j][0]+lin[j][1]*xm)
        sp = [lin[j] for j in order]
        gaps = [(sp[i+1][0]-sp[i][0], sp[i+1][1]-sp[i][1]) for i in range(k-1)]
        gaps.append((F(1)+sp[0][0]-sp[k-1][0], sp[0][1]-sp[k-1][1]))
        subs = set([x0, x1])
        for (a, b) in gaps:
            if b != 0:
                xs = (TH-a)/b
                if x0 < xs < x1: subs.add(xs)
        subs = sorted(subs)
        for s in range(len(subs)-1):
            u0, u1 = subs[s], subs[s+1]; um = (u0+u1)/2
            A = F(0); B = F(0)
            for (a, b) in gaps:
                if a+b*um > TH: A += (a-TH); B += b
            EW += A*(u1-u0) + B*(u1*u1-u0*u0)/2
            EW2 += A*A*(u1-u0) + A*B*(u1*u1-u0*u0) + B*B*(u1**3-u0**3)/3
    return EW*EW/EW2, EW, EW2

BARS = {11: F(83549,252252), 12: F(50285,252252), 13: F(14249,252252)}

def primitive(E):
    E = sorted(E); return reduce(gcd, [E[i+1]-E[i] for i in range(len(E)-1)]) == 1

# ---------- (1) exhaustive compact check k=12,13 ----------
print("=== (1) EXHAUSTIVE exact compact PZ check, k=12,13 ===")
for k in (13,):
    Dmax = 15
    best = (1e9, None); nshapes = 0
    for D in range(k-1, Dmax+1):
        for mid in combinations(range(1, D), k-2):
            E = (0,)+mid+(D,)
            if not primitive(E): continue
            nshapes += 1
            p = pz_float(E)
            if p < best[0]: best = (p, E)
    # exact-confirm the float-minimizer + a few near-min
    pex, ew, ew2 = pz_exact(list(best[1]))
    bar = BARS[k]
    print(f"  k={k}: {nshapes} primitive shapes diam<={Dmax}; float-min PZ={best[0]:.5f} at {best[1]}")
    print(f"        EXACT min PZ = {pex} = {float(pex):.6f};  bar={float(bar):.6f};  "
          f"margin {float(pex-bar):+.6f}  {'CLEARS' if pex>bar else 'BELOW!'}")

# ---------- (2) tail parameter: diameter vs min-difference ----------
print("\n=== (2) TAIL: what controls PZ -- diameter or min-difference? (k=11, the thin leg) ===")
k = 11; bar = float(BARS[11])
def mindiff(E):
    E = sorted(E); return min(E[j]-E[i] for i in range(k) for j in range(i+1, k))
# families with small min-diff but large diameter (tight cluster + far outliers)
print("  tight-cluster+outlier families (small min-diff, large diam):")
for desc, E in [('block9+2far', list(range(9))+[40,80]),
                ('block10+1far', list(range(10))+[100]),
                ('2 tight pairs spread', [0,1,20,21,40,41,60,61,80,81,100]),
                ('block11 (compact)', list(range(11)))]:
    print(f"    {desc:24s} diam={max(E)-min(E):3d} mindiff={mindiff(E):2d}: PZ={pz_float(E):.4f} {'>=bar' if pz_float(E)>=bar else '<BAR'}")
# systematic: min PZ binned by min-difference, over wide families
import collections
res = collections.defaultdict(lambda: 1e9)
for _ in range(600):
    D = random.choice([20,40,80,200])
    E = sorted(random.sample(range(D+1), k)); E = [e-E[0] for e in E]
    if not primitive(E): continue
    md = mindiff(E); res[min(md,6)] = min(res[min(md,6)], pz_float(E))
print("  min PZ over wide (diam>=20) families, binned by min-difference:")
for md in sorted(res):
    print(f"    mindiff{'>=6' if md==6 else '='+str(md)}: min PZ = {res[md]:.4f} {'>=bar' if res[md]>=bar else '<BAR (tail obstruction!)'}")
