#!/usr/bin/env python3
"""
klein-2026-07-08-S192: the WORST-CASE test for the large-spread half of THM-527-A.

The finite-Vmax glue needs: for the actual covering set, the ruler grid
x_j=(j+1/2)/Vmax hits G*={x: maxgap{frac(e_i x)}>1/7} (a good period exists).
The pigeonhole is #good_j >= rho* * Vmax - #arcs, so it needs #arcs < rho*·Vmax.

WORST CASE for the pigeonhole: spread ~ Vmax (Vmin as small as possible).  A real
cluster has u_i = Vmax - e_i > 13, so with e in [0,s] (min e=0 => the Vmax elt),
the tightest is Vmax = s + 14 (Vmin=14).  We DIRECTLY count good j in {0..Vmax-1}
for primitive co-offset sets e at Vmax = s+14, over many random e and spreads s,
and report the MIN #good (want >> 1) and min maxarc*Vmax (want > 1).

Also imposes a realistic G_P: P = {1,...,13-k} (the small part), x in G_P iff
||p x||>=1/14 for all p in P.  Good period must be in G_P too.
"""
import numpy as np
from math import gcd
from functools import reduce

INV7 = 1.0/7.0

def gcd_list(xs):
    return reduce(gcd, xs)

def good_count_on_ruler(e, Vmax, P):
    """count j in {0..Vmax-1} with x_j=(j+1/2)/Vmax in G_P and maxgap{frac(e_i x_j)}>1/7."""
    j = np.arange(Vmax)
    x = (j + 0.5)/Vmax
    # G_P membership
    inGP = np.ones(Vmax, dtype=bool)
    for p in P:
        d = np.abs(((p*x + 0.5) % 1.0) - 0.5)   # ||p x||
        inGP &= (d >= 1.0/14 - 1e-12)
    # maxgap of cluster phases
    e = np.asarray(e)
    Ph = np.mod(np.outer(x, e), 1.0)
    Ph.sort(axis=1)
    g = np.diff(Ph, axis=1)
    g = np.concatenate([g, (1.0 - Ph[:,-1] + Ph[:,0])[:,None]], axis=1)
    mg = g.max(axis=1)
    good = inGP & (mg > INV7 + 1e-12)
    return int(good.sum()), float(inGP.mean())

def arcs_and_maxarc(e, N=200000):
    e=np.asarray(e); x=(np.arange(N)+0.5)/N
    Ph=np.mod(np.outer(x,e),1.0); Ph.sort(axis=1)
    g=np.diff(Ph,axis=1); g=np.concatenate([g,(1.0-Ph[:,-1]+Ph[:,0])[:,None]],axis=1)
    good=(g.max(axis=1)>INV7+1e-12).astype(int)
    edges=np.diff(np.concatenate([good,good[:1]]))
    nc=int((edges==1).sum());
    if good.all(): nc=1
    gg=np.concatenate([good,good]); runs=[];c=0
    for v in gg:
        if v: c+=1
        else:
            if c>0: runs.append(c);c=0
    if c>0: runs.append(c)
    maxarc=(max(runs) if runs else 0)/N
    return nc, maxarc, good.mean()

rng = np.random.default_rng(31415926)
print("Worst-case (Vmax = spread+14, Vmin=14) large-spread test.")
print("Report MIN over random primitive clusters of: #good_j on the ruler, and maxarc*Vmax.\n")

for k in (11, 12, 13):
    P = list(range(1, 13-k+1))   # |P| = 13-k ; k=13 => P empty
    print(f"===== k={k}, P={P} =====")
    for s in (30, 60, 120, 250, 500, 1000):
        Vmax = s + 14
        min_good = 10**9; min_ma_v = 1e9; min_rho=1.0; n_fail=0; trials=25
        for _ in range(trials):
            # primitive co-offsets: 0 plus k-1 distinct in [1,s], gcd=1, and s present
            e = rng.choice(np.arange(1, s), size=k-2, replace=False)
            e = np.concatenate([[0], e, [s]])   # ensure spread exactly s
            e = e // gcd_list([int(v) for v in e]) if gcd_list([int(v) for v in e])>1 else e
            # ensure cluster valid: u_i = Vmax - e_i > 13  => e_i < Vmax-13 = s+1 ok
            ng, rhoGP = good_count_on_ruler(e, Vmax, P)
            if ng < min_good: min_good = ng
            if ng == 0: n_fail += 1
            nc, ma, rho = arcs_and_maxarc(e, N=120000)
            min_ma_v = min(min_ma_v, ma*Vmax)
            min_rho = min(min_rho, rho)
        print(f"  s={s:5d} Vmax={Vmax:5d}: MIN #good_j={min_good:6d}  (fails:{n_fail}/{trials})  "
              f"min maxarc*Vmax={min_ma_v:5.2f}  min rho*={min_rho:.3f}")
    print()

print("If MIN #good_j >= 1 (no fails) with big margin and maxarc*Vmax>1 across all")
print("spreads at the WORST Vmax=s+14, the large-spread half holds: a good ruler")
print("period always exists.  The pigeonhole #good >= rho*Vmax-#arcs then just needs")
print("the rigorous linear arc-count bound #arcs<=c*spread with c<rho* (~0.2 observed).")
