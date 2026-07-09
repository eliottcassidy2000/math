"""
mac-mini-2026-07-08-S59 -- THM-527-A large-spread half via WEYL (first-moment route).

W(y) = sum_i (g_i(y)-1/7)_+ = uncovered measure >= 0.  A good period (j in 1..Vmax-1 with
maxgap>1/7) exists  <=>  SUM_{j=1}^{Vmax-1} W(j/Vmax) > 0  <=>  Vmax*E_grid[W] - W(0) > 0,
W(0)=6/7.  So it suffices:  E_grid[W] := (1/Vmax) sum_{j=0}^{Vmax-1} W(j/Vmax)  >  6/(7 Vmax).

Inclusion-exclusion:  E_grid[W] = sum_{S subseteq [k]} (-1)^|S| E_grid[O_S],  O_S = overlap of the
|S| arcs = (1/7 - span_S)_+.  Main term (iid, O_S->(1/7)^|S|):  E_grid[W] -> (6/7)^k > 0.
The Weyl corrections are the grid-vs-iid deviations E_grid[O_S]-(1/7)^|S|, each an exponential
sum over resonances n.e ≡ 0 (mod Vmax) -> 0 as the cluster decorrelates (large spread).

This computes E_grid[W] EXACTLY and checks (a) E_grid[W] > 6/(7Vmax) (=> good period), and
(b) how close E_grid[W] is to the iid (6/7)^k, on the WORST (structured) large-spread clusters.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce

def W_at(E, j, Vmax):
    """W(j/Vmax) EXACT: uncovered measure = sum (gap-1/7)_+, gaps in units of 1, as Fraction."""
    ph = sorted({(e*j) % Vmax for e in E})     # phases * Vmax (ints), distinct
    m = len(ph)
    TH = F(1, 7)
    if m == 1:
        return F(6, 7)                          # single phase: one gap =1 => W = 1-1/7
    W = F(0)
    for i in range(m):
        if i < m-1:
            g = F(ph[i+1]-ph[i], Vmax)
        else:
            g = F(ph[0]+Vmax-ph[-1], Vmax)
        if g > TH:
            W += g - TH
    return W

def egrid_W(E, Vmax):
    """exact E_grid[W] and SUM_{j>=1} W(j/Vmax)."""
    tot = F(0)
    for j in range(Vmax):
        tot += W_at(E, j, Vmax)
    eg = tot / Vmax
    sum_pos = tot - F(6, 7)                      # subtract j=0 term W(0)=6/7
    return eg, sum_pos

def primitive(E):
    E = sorted(E); return reduce(gcd, [E[i+1]-E[i] for i in range(len(E)-1)]) == 1

def twoblock(k, s):
    a = k//2; b = k-a
    return sorted(set(list(range(a)) + [s-b+1+i for i in range(b)]))
def apshape(k, s):
    d = max(1, s//(k-1)); return sorted(set([d*i for i in range(k-1)] + [s]))

print("WEYL first-moment: E_grid[W] vs iid (6/7)^k vs threshold 6/(7Vmax), structured large-spread\n")
for k in (11, 13):
    iid = (6/7)**k
    print(f"k={k}: iid (6/7)^k = {iid:.5f}")
    for name, mk in [('2-block', twoblock), ('AP', apshape)]:
        for s in [30, 60, 120]:
            E = mk(k, s)
            if len(E) != k or not primitive(E): continue
            for Vmax in [s+1, s+7, 2*s+1]:        # large-spread: spread ~ Vmax
                if Vmax <= max(E): continue
                eg, spos = egrid_W(E, Vmax)
                thr = F(6, 7*Vmax)
                ok = eg > thr
                print(f"   {name:8s} s={s:3d} Vmax={Vmax:4d}: E_grid[W]={float(eg):.5f} "
                      f"(iid {iid:.4f}, dev {float(eg)-iid:+.4f}); sum_{{j>=1}}W={float(spos):.3f}>0? {spos>0}; "
                      f">6/(7Vmax)? {ok}")
    print()
