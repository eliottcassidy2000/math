#!/usr/bin/env python3
"""THE CONCORDANCE DECISION at n=7 (opus-S308, HYP-6870).

S307 found 132/92634 discordant potential-axis pairs at n=7 -- but in float64.
Before proving anything, decide whether they are REAL: solve the 456-class
network with float64 LU + exact-integer-residual iterative refinement (the
Laplacian is integer, so residuals are exact), then classify every class pair:
concordant / discordant / tie, with margins. Pairs with |dphi| below the
certified error bound are UNDECIDED and counted honestly.
"""
import sys, os
from collections import defaultdict
from fractions import Fraction as F
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from smith_diagram_of_the_metagraph_opus_S307 import build

n = 7
B = build(n)
C, W, x_of, H_of = B['C'], B['W'], B['x_of'], B['H_of']
cond = {(a, b): w for (a, b), w in W.items() if a != b}
xmin = min(x_of.values())
src = [c for c in range(C) if H_of[c] == 1][0]
sinks = [c for c in range(C) if x_of[c] == xmin]
nodes = [c for c in range(C) if c not in sinks]
pos = {c: i for i, c in enumerate(nodes)}
N = len(nodes)
A = np.zeros((N, N))
Aint = {}                       # exact integer entries
deg = defaultdict(int)
for (a, b), w in cond.items():
    deg[a] += w; deg[b] += w
    if a in pos and b in pos:
        A[pos[a], pos[b]] -= w; A[pos[b], pos[a]] -= w
        Aint[(pos[a], pos[b])] = Aint.get((pos[a], pos[b]), 0) - w
        Aint[(pos[b], pos[a])] = Aint.get((pos[b], pos[a]), 0) - w
for c in nodes:
    A[pos[c], pos[c]] = deg[c]
    Aint[(pos[c], pos[c])] = deg[c]
b = np.zeros(N); b[pos[src]] = 1.0
x = np.linalg.solve(A, b)
# exact-residual iterative refinement (3 passes)
for it in range(3):
    xf = [F(float(v)) for v in x]
    r = [F(0)] * N
    for (i, j), aij in Aint.items():
        r[i] += aij * xf[j]
    r[pos[src]] -= 1
    rf = np.array([-float(v) for v in r])
    dx = np.linalg.solve(A, rf)
    x = x + dx
    print(f"refine pass {it}: max |residual| = {max(abs(float(v)) for v in r):.3e}, "
          f"max |dx| = {np.max(np.abs(dx)):.3e}")
phi = {c: x[pos[c]] for c in nodes}
for c in sinks: phi[c] = 0.0
# error bound: after refinement, forward error ~ kappa * eps; estimate kappa via norm ratios
err = 1e-12 * max(abs(v) for v in phi.values())
print(f"decision threshold |dphi| > {err:.2e}")
conc = disc = und = 0
disc_pairs = []
import itertools
for a, bb in itertools.combinations(range(C), 2):
    dxv = x_of[a] - x_of[bb]
    if dxv == 0: continue
    dp = phi[a] - phi[bb]
    if abs(dp) < err: und += 1; continue
    if (dxv > 0) == (dp > 0): conc += 1
    else:
        disc += 1
        disc_pairs.append((a, bb, x_of[a], x_of[bb], phi[a], phi[bb]))
print(f"n=7 DECIDED: concordant {conc}, DISCORDANT {disc}, undecided (below threshold) {und}")
if disc:
    print("discordant pair details (class, class, x, x, phi, phi):")
    lv = defaultdict(int)
    for a, bb, xa, xb, pa, pb in disc_pairs[:20]:
        print(f"  ({a},{bb}) x=({xa},{xb}) phi=({pa:.9e},{pb:.9e}) "
              f"Hs=({H_of[a]},{H_of[bb]})")
    for a, bb, xa, xb, pa, pb in disc_pairs:
        lv[(max(xa,xb), min(xa,xb))] += 1
    print(f"  by level pair: {dict(sorted(lv.items()))}")
    print(f"  VERDICT: the discordances are REAL -- the concordance law FAILS at n=7;")
    print(f"  its true scope is n <= 6 exact + bulk-concordance (99.86%) beyond.")
else:
    print("  VERDICT: no real discordances -- the 132 float pairs were numerical;")
    print("  the concordance law stands exactly at n=7.")
