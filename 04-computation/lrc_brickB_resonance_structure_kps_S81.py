#!/usr/bin/env python3
r"""
lrc_brickB_resonance_structure_kps_S81.py  (kind-pasteur-2026-07-08-S81, HYP-5357)

STRUCTURAL clarification of brick (B) = the resonance lemma Var(W) <= c*R2 (the open analytic
mile for the k=11 covering tail).  NOT a proof -- it pinpoints the residual.

THE EXACT DECOMPOSITION.  W = uncovered measure = int_0^1 prod_i (1 - u_i(y)) dy, u_i(y) =
1[e_i x-phase covers y].  Expand: W = sum_{S subset [k]} (-1)^|S| L_S, L_S = meas(cap_{i in S}
A_i) = overlap of the S-arcs.
  * |S|=0: L = 1 (const).  |S|=1: L_i = 1/7 (CONSTANT) => both contribute 0 to Var(W).
  * So Var(W) = sum_{|S|,|T| >= 2} (-1)^{|S|+|T|} Cov(L_S, L_T) -- STARTS AT PAIRS.
The leading (|S|=|T|=2) term is Var(sum_{i<j} ov_ij), additive-energy-linear (THM-641 gives it
EXACTLY: the pair overlap masses).  The higher terms (|S| or |T| >= 3) are triple/quad overlap
correlations -- NO proved law yet.

THE RESONANCE = the (large, ~96%) cancellation of the pair-sum variance by the triple+ terms:
for the k=11 block, Var(sum ov_ij) ~ 1.26 but Var(W) = 0.047 (3.7% of it).  So brick (B) reframes
from 'prove Var(W) <= c*R2 (the resonance sign)' to 'derive the 3-arc and 4-arc overlap mass
laws (the THM-641 analog for triples/quads)' -- a concrete, derivable target.
"""
import numpy as np
from collections import Counter

def W_of(E, x):
    ph = sorted((e*x) % 1.0 for e in E); n = len(ph)
    g = [ph[(i+1) % n] - ph[i] if i < n-1 else ph[0] + 1 - ph[-1] for i in range(n)]
    return sum(max(gi - 1/7, 0) for gi in g)
def ov(d, x):
    t = (d*x) % 1.0; dd = min(t, 1-t); return max(1/7 - dd, 0)
def R2(E):
    r = Counter()
    for i in range(len(E)):
        for j in range(len(E)):
            if i != j: r[E[i]-E[j]] += 1
    return sum(v*v for v in r.values())

for E in [list(range(11)), [0,1,2,3,4,5,6,7,8,9,16]]:
    res = 200000; xs = (np.arange(res)+0.5)/res
    W = np.array([W_of(E, x) for x in xs]); VarW = W.var()
    pairs = [(i, j) for i in range(11) for j in range(i+1, 11)]
    Ssum = np.zeros(res)
    for (i, j) in pairs:
        d = E[j] - E[i]; Ssum += np.array([ov(d, x) for x in xs])
    Varpair = Ssum.var(); r2 = R2(E)
    print(f"E={E}: R2={r2}")
    print(f"  Var(W) = {VarW:.6f} = R2*{VarW/r2:.2e}")
    print(f"  Var(pair-sum) = {Varpair:.6f}; Var(W)/Var(pair-sum) = {VarW/Varpair:.4f}"
          f" => {100*(1-VarW/Varpair):.1f}% higher-order cancellation (= the RESONANCE)")
print()
print("|S|=1 terms (L_i=1/7 const) => 0 to Var; THM-641 gives the pair term exactly;")
print("residual = the triple/quad arc-overlap mass laws (THM-641's method should extend).")
