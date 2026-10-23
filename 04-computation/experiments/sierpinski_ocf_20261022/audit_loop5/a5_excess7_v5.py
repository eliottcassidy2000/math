# Audit (HYP-9371 section 8, "Next step"): is excess 7 needed at z^12?  At L = 12 an excess-7 cluster has 12 - 7 = 5
# blocks, so its homomorphism sums see the 5-vertex census, which varies across doubly regular tournaments of one
# order (affine in TT5).  Excess-7 cores with V <= 4 are universal (<= 4-census); with V = 5 they have E = 12 edges,
# so at z^12 every path is a single arc, and a loop then gives S_vv = 0.  So the non-universal part of [z^12] X_7 is
#     Delta X_7^(12) = sum over LOOPLESS excess-7 cores on 5 vertices of K(core, 1..1) Delta hom(core; S) / (|Aut| 23^12).
# Own code: all compositions of 12 edges into the 10 vertex pairs (exhaustive), minimum degree >= 3, connected,
# brute-force canonical forms over the 120 permutations; K from a5_kernel.c (first principles); exact integer hom
# numbers on the 14 DRT(23)s.  Usage: python3 a5_excess7_v5.py KERNEL_BINARY
import os, sys, json, itertools, time
from fractions import Fraction
from math import factorial
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); D = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, D)
import a5_z10_lib as Z
from drt_switch import code_to_S

binary = sys.argv[1]
V, E = 5, 12
pairs = [(a, b) for a in range(V) for b in range(a + 1, V)]             # loopless
t0 = time.time()
P = np.array(list(itertools.permutations(range(V))), dtype=np.int64)
iu = np.triu_indices(V, 1)
canon = {}; nlab = 0
for bars in itertools.combinations(range(E + len(pairs) - 1), len(pairs) - 1):   # stars and bars
    mv = []; prev = -1
    for b in bars: mv.append(b - prev - 1); prev = b
    mv.append(E + len(pairs) - 1 - prev - 1)
    M = np.zeros((V, V), dtype=np.int64)
    for (a, b), m in zip(pairs, mv): M[a, b] = M[b, a] = m
    if M.sum(1).min() < 3: continue
    # connectivity
    seen = {0}; st = [0]
    while st:
        x = st.pop()
        for y in range(V):
            if M[x, y] and y not in seen: seen.add(y); st.append(y)
    if len(seen) < V: continue
    nlab += 1
    vecs = M[P[:, :, None], P[:, None, :]][:, iu[0], iu[1]]
    c = tuple(int(x) for x in vecs[np.lexsort(vecs.T[::-1])[0]])
    if c in canon: continue
    vaut = int(np.sum(np.all(vecs == M[iu], axis=1)))
    f = 1
    for m in c: f *= factorial(m)
    canon[c] = vaut * f
# orbit-counting check: sum over classes of V!/|vertex Aut| = number of labelled graphs
orbit_sum = sum(factorial(V) // (aut // np.prod([factorial(m) for m in c])) for c, aut in canon.items())
print(f"loopless connected multigraphs on 5 vertices, 12 edges, min degree >= 3: {nlab} labelled, {len(canon)} classes "
      f"(orbit-stabiliser check {orbit_sum} = {nlab}: {orbit_sum == nlab}) ({time.time() - t0:.0f}s)", flush=True)
cores = []
for c, aut in canon.items():
    edges = []
    for (a, b), m in zip(pairs, c): edges += [(a, b)] * m
    cores.append((edges, aut))
Ks = [Z.kernels(binary, V, edges, [tuple([1] * E)])[0] for edges, aut in cores]
print(f"kernels K(core, 1..1): {sum(1 for k in Ks if k)} nonzero of {len(Ks)} ({time.time() - t0:.0f}s)", flush=True)
cls = json.load(open(os.path.join(D, 'drt23_classes.json'))); dbl = json.load(open(os.path.join(D, 'drt23_doubled.json')))
tours = [(c['label'], c['code'], c['TT5']) for c in cls] + [('D23_dbl', dbl['code'], dbl['TT5'])]
vals = {}
for lab, code, tt in tours:
    H = Z.Homs(code_to_S(code), 1)
    s = Fraction(0)
    for (edges, aut), K in zip(cores, Ks):
        if K: s += Fraction(K * H.hom(V, edges, tuple([1] * E)), aut)
    vals[lab] = (s / 23 ** 12, tt)
P23 = vals['P_23']
print(f"P_23: [z^12] X_7 restricted to the V = 5 cores = {float(P23[0]):+.6e}")
print("Delta X_7^(12) (exact, V = 5 cores) relative to P_23, and per unit Delta TT5:")
for lab, (v, tt) in vals.items():
    if lab == 'P_23': continue
    dv = v - P23[0]
    print(f"   {lab:8s} Delta TT5 = {tt - P23[1]:4d}: Delta X_7^(12) = {float(dv):+.6e}  (per TT5: {float(dv / (tt - P23[1])):+.6e})")
print(f"({time.time() - t0:.0f}s)")
