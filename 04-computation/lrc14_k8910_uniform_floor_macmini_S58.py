"""
mac-mini-2026-07-08-S58 -- close the k=8,9,10 (A') legs at the UNIFORM (min over
families) level, mirroring the k=11,12,13 D3 closure. THM-661's table gives B_d
for the BLOCK; the union bound (THM-530) needs min_E >= bar. Check whether the
clean degree-3 D3 clears k=8,9,10 uniformly; if not, use degree-4 B_4.
Float pre-filter (fast) + exact-confirm (Farey) of candidates near the min.
"""
import sys
from fractions import Fraction as F
from math import floor, gcd
from functools import reduce
from itertools import combinations
import numpy as np
from scipy.optimize import linprog

TH = F(1, 7); M = F(6, 7); Mf = 6.0/7.0
BARS = {8: 0.6750, 9: 0.5622, 10: 0.4521}
BARS_EXACT = {8: None, 9: None, 10: None}   # decimals suffice (margins are large)

# ---- fast float moments/D3 ----
GRIDf = 8000
_xf = (np.arange(GRIDf) + 0.5) / GRIDf
def moments_float(E):
    Ea = np.array(sorted(E), float)
    ph = np.mod(np.outer(_xf, Ea), 1.0); ph.sort(axis=1)
    g = np.concatenate([np.diff(ph, axis=1), (ph[:, 0]+1-ph[:, -1])[:, None]], axis=1)
    W = np.maximum(g - 1.0/7.0, 0).sum(axis=1)
    return [None, W.mean(), (W*W).mean(), (W**3).mean(), (W**4).mean()]
def d3_float(m):
    den = m[2] - m[3]/Mf
    return (m[1]/Mf + (m[1]-m[2]/Mf)**2/den) if den > 1e-12 else m[1]/Mf

# ---- exact moments/D3 (Farey) ----
def int_pow(A, B, p, lo, hi):
    if B == 0: return A**p * (hi - lo)
    return ((A + B*hi)**(p+1) - (A + B*lo)**(p+1)) / (B * (p+1))
def W_moments(E, maxp=4):
    E = sorted(E); k = len(E)
    ds = set(e for e in E if e != 0)
    for i in range(k):
        for j in range(i+1, k): ds.add(E[j]-E[i])
    bps = set([F(0), F(1)])
    for d in ds:
        for m in range(0, d+1): bps.add(F(m, d))
    bps = sorted(bps)
    mom = [F(0)]*(maxp+1)
    for c in range(len(bps)-1):
        x0, x1 = bps[c], bps[c+1]; xm = (x0+x1)/2
        lin = [(F(-floor(e*xm)), F(e)) for e in E]
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
            A = B = F(0)
            for (a, b) in gaps:
                if a + b*um > TH: A += (a-TH); B += b
            for p in range(1, maxp+1):
                mom[p] += int_pow(A, B, p, u0, u1)
    return mom
def D3_exact(mom):
    m1, m2, m3 = mom[1], mom[2], mom[3]
    den = m2 - m3/M
    return (m1/M + (m1 - m2/M)**2/den) if den > 0 else m1/M

_wgrid = np.linspace(1e-4, Mf, 800)
_Aub = np.column_stack([_wgrid, _wgrid**2, _wgrid**3, _wgrid**4])
def B4(mom):   # degree-4 moment LP lower bound on mu
    m = [float(mom[i]) for i in range(5)]
    res = linprog([-m[1], -m[2], -m[3], -m[4]], A_ub=_Aub, b_ub=np.ones(len(_wgrid)),
                  bounds=[(None, None)]*4, method='highs')
    return -res.fun if res.success else None

def primitive(E):
    E = sorted(E); return reduce(gcd, [E[i+1]-E[i] for i in range(len(E)-1)]) == 1

print("k=8,9,10 UNIFORM density floor: exhaustive compact min (float pre-filter + exact)\n", flush=True)
for k in (8, 9, 10):
    bar = BARS[k]
    Dmax = k + 4
    blk = tuple(range(k)); mblk = W_moments(blk)
    d3blk = float(D3_exact(mblk)); b4blk = B4(mblk)
    print(f"k={k} (bar={bar:.4f}):  D3(block)={d3blk:.4f}  B4(block)={b4blk:.4f}", flush=True)
    # float scan -> keep candidates with low D3 OR low B4
    cands = []; nshapes = 0; fmin = 9.9
    for D in range(k-1, Dmax+1):
        for mid in combinations(range(1, D), k-2):
            E = (0,) + mid + (D,)
            if not primitive(E): continue
            nshapes += 1
            mf = moments_float(E); d3f = d3_float(mf)
            fmin = min(fmin, d3f)
            cands.append((d3f, E))
    cands.sort()
    # exact-confirm the lowest ~40 by float D3, and B4 on each
    bestD3 = (F(10), None); bestB4 = (9.9, None)
    for d3f, E in cands[:40]:
        mom = W_moments(list(E))
        d3 = D3_exact(mom); b4 = B4(mom)
        if d3 < bestD3[0]: bestD3 = (d3, E)
        if b4 < bestB4[0]: bestB4 = (b4, E)
    print(f"   {nshapes} shapes prim-diam<= {Dmax}; float D3 min = {fmin:.4f}", flush=True)
    print(f"   EXACT min_E D3 = {float(bestD3[0]):.6f} at {bestD3[1]}  "
          f"margin {float(bestD3[0])-bar:+.4f}  {'CLEARS' if float(bestD3[0])>bar else 'BELOW bar (need deg4)'}", flush=True)
    print(f"   min_E B4 (deg-4) = {bestB4[0]:.6f} at {bestB4[1]}  "
          f"margin {bestB4[0]-bar:+.4f}  {'CLEARS' if bestB4[0]>bar else 'BELOW!!'}\n", flush=True)
