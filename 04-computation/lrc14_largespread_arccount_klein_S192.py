#!/usr/bin/env python3
"""
klein-2026-07-08-S192: does #arcs(G*) stay BOUNDED (Vmax-free) for LARGE-SPREAD
clusters?  This is the sole open half of THM-527-A (the finite-Vmax glue).

G*(E) = { x in [0,1) : maxgap{frac(e_i x)} > 1/7 }  (the good set; the cluster
E fits in a <5/7 arc complement -> a covering witness).  mac-mini-S58 proved
#arcs depends only on cluster-INTERNAL differences and MEASURED #arcs ~ k+1 for
small-spread blocks/near-2-APs, but bounded #arcs = O(k^2 spread^2) grows with
spread.  If #arcs is ACTUALLY bounded (~O(k)) even for spread ~ Vmax, the
large-spread half closes by the same finite check: a good period exists once
Vmax > #arcs/m_P, and Vmax <= that is finite.

Test: k=11..13 clusters with spread ranging from O(k) up to several thousand
(random + structured), count connected components of G* on a fine grid, and
watch whether #arcs grows with spread.  Also report max-arc-length * Vmax_proxy.
"""
import numpy as np
from math import gcd

INV7 = 1.0/7.0
m_P  = 14249.0/252252.0     # the uniform rho* floor (THM-530)

def maxgap(phases):
    p = np.sort(np.mod(phases, 1.0))
    g = np.diff(p); g = np.append(g, 1.0 - p[-1] + p[0])
    return g.max()

def arc_stats(E, N=400000):
    """#connected components of {x: maxgap{frac(e_i x)}>1/7} and max arc length."""
    E = np.asarray(E, dtype=np.float64)
    x = (np.arange(N)+0.5)/N
    # phases: N x k ; maxgap per row
    Ph = np.mod(np.outer(x, E), 1.0)
    Ph.sort(axis=1)
    g = np.diff(Ph, axis=1)
    g = np.concatenate([g, (1.0 - Ph[:,-1] + Ph[:,0])[:,None]], axis=1)
    mg = g.max(axis=1)
    good = mg > INV7 + 1e-12
    # count components on the circle (wrap)
    d = np.diff(good.astype(int))
    ups = np.sum(d==1); downs=np.sum(d==-1)
    # wrap correction: treat as circular
    ncomp = ups
    if good[0] and good[-1]:
        # first and last are same component if both good
        ncomp = max(ncomp, downs)  # components = number of up-edges on circle
        ncomp = ups if ups>0 else (1 if good.all() else 0)
    else:
        ncomp = ups + (1 if (good[0] and not good[-1]) else 0)
    # robust circular component count:
    gi = good.astype(int)
    edges = np.diff(np.concatenate([gi, gi[:1]]))  # circular
    ncomp = int(np.sum(edges==1))
    if gi.all(): ncomp = 1
    frac = good.mean()
    # max arc length (circular run of good)
    if ncomp==0:
        maxarc=0.0
    else:
        # find run lengths circularly
        gg = np.concatenate([gi,gi])
        runs=[]; c=0
        for v in gg:
            if v: c+=1
            else:
                if c>0: runs.append(c); c=0
        if c>0: runs.append(c)
        maxarc = (max(runs) if runs else 0)/N
    return ncomp, frac, maxarc

print(f"m_P (rho* floor) = {m_P:.5f}   1/m_P = {1/m_P:.1f}")
print("Testing #arcs(G*) vs spread.  If #arcs stays ~O(k) for large spread,")
print("the large-spread half of THM-527-A closes (Vmax > #arcs/m_P finite check).\n")

rng = np.random.default_rng(20260708)

for k in (11, 12, 13):
    print(f"===== k={k} =====")
    # (a) blocks / near-blocks at increasing scale (spread grows, structured)
    for scale in (1, 10, 100, 1000):
        E = list(scale*np.arange(k))
        nc,fr,ma = arc_stats(E)
        print(f"  block*{scale:<5d} spread={scale*(k-1):6d}: #arcs={nc:3d} frac={fr:.4f} maxarc={ma:.5f}")
    # (b) RANDOM large-spread clusters (co-offsets e_i random in [0,S])
    for S in (30, 100, 300, 1000, 3000):
        ncs=[]; frs=[]; mas=[]
        for _ in range(12):
            E = rng.choice(np.arange(1,S+1), size=k-1, replace=False)
            E = np.concatenate([[0], E])
            nc,fr,ma = arc_stats(E, N=300000)
            ncs.append(nc); frs.append(fr); mas.append(ma)
        ncs=np.array(ncs)
        print(f"  random spread<= {S:5d}: #arcs mean={ncs.mean():5.1f} max={ncs.max():3d} "
              f"| frac mean={np.mean(frs):.4f} min={np.min(frs):.4f} | min maxarc={np.min(mas):.5f}")
    print()

print("READ: if #arcs 'max' stays small (bounded, not growing ~S^2) and frac stays")
print(">= ~m_P, then Vmax > #arcs/frac ~ O(k)/m_P is a UNIFORM threshold => good period")
print("exists for Vmax above it, finite check below. That closes THM-527-A large-spread.")
