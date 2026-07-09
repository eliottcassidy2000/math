"""
mac-mini-2026-07-08-S58 -- backstop for the k=8,9,10 UNIFORM density floor:
confirm the TAIL (large prim-diam) never dips below the BLOCK value of B4
(deg-4 moment bound), so min_E B4 = B4(block) >= bar. Block is the compact
minimizer (exhaustive); decorrelation raises the moments, so the block is global.
"""
import numpy as np
from math import gcd
from functools import reduce
import random
random.seed(89)
from scipy.optimize import linprog

BARS = {8: 0.6750, 9: 0.5622, 10: 0.4521}
B4_BLOCK = {8: 0.761132, 9: 0.644603, 10: 0.553111}   # exact-confirmed compact min (block)
Mf = 6.0/7.0
GRIDf = 8000
_xf = (np.arange(GRIDf) + 0.5) / GRIDf
def moments_float(E):
    Ea = np.array(sorted(E), float)
    ph = np.mod(np.outer(_xf, Ea), 1.0); ph.sort(axis=1)
    g = np.concatenate([np.diff(ph, axis=1), (ph[:, 0]+1-ph[:, -1])[:, None]], axis=1)
    W = np.maximum(g - 1.0/7.0, 0).sum(axis=1)
    return [W.mean(), (W*W).mean(), (W**3).mean(), (W**4).mean()]
_wgrid = np.linspace(1e-4, Mf, 800)
_Aub = np.column_stack([_wgrid, _wgrid**2, _wgrid**3, _wgrid**4])
def B4_float(E):
    m = moments_float(E)
    res = linprog([-m[0], -m[1], -m[2], -m[3]], A_ub=_Aub, b_ub=np.ones(len(_wgrid)),
                  bounds=[(None, None)]*4, method='highs')
    return -res.fun if res.success else None
def primitive(E):
    E = sorted(E); return reduce(gcd, [E[i+1]-E[i] for i in range(len(E)-1)]) == 1

print("k=8,9,10 backstop: does any large-diameter tail shape undercut B4(block)?\n", flush=True)
for k in (8, 9, 10):
    bar = BARS[k]; blkv = B4_BLOCK[k]
    gmin = (9.9, None); n = 0
    # structured near-AP + random large-diameter primitive shapes
    for d in range(1, 9):
        ap = [d*i for i in range(k-1)]; hi = ap[-1]
        for p in list(range(1, hi)) + [hi+1, hi+2, 2*hi, 3*hi]:
            E = sorted(set(ap + [p]))
            if len(E) != k or not primitive(E): continue
            n += 1; v = B4_float(E)
            if v is not None and v < gmin[0]: gmin = (v, tuple(E))
    for _ in range(4000):
        D = random.choice([15, 20, 30, 50, 100, 200])
        E = sorted(random.sample(range(1, D+1), k-1)); E = [0]+E
        if len(set(E)) != k or not primitive(E): continue
        n += 1; v = B4_float(E)
        if v is not None and v < gmin[0]: gmin = (v, tuple(E))
    ok = gmin[0] >= blkv - 0.003
    print(f"k={k} (bar={bar:.4f}, B4 block={blkv:.4f}): scanned {n} tail shapes", flush=True)
    print(f"   min tail B4 = {gmin[0]:.5f} at {gmin[1]}", flush=True)
    print(f"   tail {'>=' if ok else '<'} block (tol 3e-3): "
          f"{'YES -- block is global min' if ok else 'NO -- INVESTIGATE'};  "
          f"tail margin vs bar {gmin[0]-bar:+.4f}\n", flush=True)
