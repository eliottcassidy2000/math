#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_threadB_mobius_constant_degeneracy_opus_0621.py   (opus, 2026-06-21, THREAD B honest close)

THE HONEST CLOSE.  The crux run revealed the consec-max "certificate" the slack-LP found is the
TRIVIAL CONSTANT functional  F(E) = measS7(consec)  (y_r = [measS7(consec), 0, 0, ...]).  This is
CIRCULAR: F(E)=const is valid (V: F>=measS7) iff measS7(E)<=measS7(consec) for all E, which IS the
consec-max claim.  So the "certificate exists / non-vacuity" results were detecting consec-as-max
TAUTOLOGICALLY, not proving it from the subset lattice.

This script makes the degeneracy explicit and asks the only NON-circular question the Mobius
route can pose:

  Q.  Is there a degree-R functional that certifies consec-max WITHOUT absorbing validity into a
      constant, i.e. a functional whose validity (V) F(E)>=measS7(E) is provable PER-E from
      structure (the actual CJJ aim), AND that is consec-maximal?

  The clean operationalization: the ONLY per-E-provable validity we have on file is the Bonferroni
  even-truncation  B_R = sum_{r<=R, even} ... actually  B_R(E) = sum_{r=0}^R (-1)^r S_r(E) >= measS7(E)
  for EVEN R (Bonferroni, PROVED per-E, THM-534/HYP-2726).  So the NON-circular degree-R Mobius/
  Delsarte certificate must be a NONNEGATIVE combination of {B_2, B_4, ...} (and the Krawtchouk
  c_j>=0 duals), NOT the free constant.  We therefore RE-RUN the slack-LP but restrict the
  functional to the CONE of per-E-valid Bonferroni/Krawtchouk functionals (the genuine level-1
  Delsarte dual cone), and ask: is consec the MAX of such a functional over the window?

  Concretely: F = sum_j c_j M_j(E) with c_j>=0 (M_j = Krawtchouk moments, the Delsarte-valid cone;
  F>=measS7 per-E is then automatic).  min over c>=0 (normalized) of [max_E F(E) - F(consec)].
  If 0 -> consec is the argmax of a Delsarte-cone functional => clean consec-max.  If >0 -> even
  the full Delsarte cone cannot make consec the argmax (the bound's consec-max is genuinely
  aggregate, not a single Delsarte inequality).
"""
import sys, itertools
from math import comb, gcd
from fractions import Fraction as F
from functools import reduce, lru_cache
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)
H = F(1, 14)
INNER = list(range(1, 7))
ALL_SUBSETS = [frozenset(b for b in INNER if (mask >> (b - 1)) & 1) for mask in range(64)]
from scipy.optimize import linprog
import numpy as np

def Kraw(j, t, n=6):
    return sum((-1) ** i * comb(t, i) * comb(n - t, j - i) for i in range(j + 1))
KTAB = [[Kraw(j, t) for t in range(7)] for j in range(7)]

def occupancy_p(E):
    """p_t = meas{#missed inner sectors = t}, t=0..6."""
    E = sorted(set(E)); Enz = [e for e in E if e != 0]
    bps = {F(0), F(1)}
    for e in Enz:
        for m in range(0, 7 * e + 1):
            bps.add(F(m, 7 * e))
    bps = sorted(bps)
    p = [F(0)] * 7
    for x0, x1 in zip(bps, bps[1:]):
        if x1 <= x0: continue
        xm = (x0 + x1) / 2
        hit = set(int(7 * e * xm) % 7 for e in E)
        t = 6 - len([s for s in hit if s != 0])
        p[t] += x1 - x0
    return p

def Mvec(p): return [sum(KTAB[j][t] * p[t] for t in range(7)) for j in range(7)]
def Svec(p): return [sum(comb(t, r) * p[t] for t in range(7)) for r in range(7)]
def measS7(p): return p[0]
def primitive(E): return reduce(gcd, [e for e in E if e != 0], 0) == 1
def consec(k): return list(range(k))
WINDOWS = {8: 16, 9: 15, 10: 14}

def build(k):
    recs = []
    for rest in itertools.combinations(range(1, WINDOWS[k] + 1), k - 1):
        E = [0] + list(rest)
        if not primitive(E): continue
        p = occupancy_p(E)
        recs.append((tuple(E), p))
    return recs

print("=" * 100)
print("HONEST CLOSE: the slack-LP 'certificate' is the TRIVIAL CONSTANT F=measS7(consec) (circular).")
print("Real question: is consec the ARGMAX over the window of a PER-E-VALID Delsarte-cone functional")
print("  F = sum_j c_j M_j, c_j>=0 (g(t)=sum c_j K_j(t)>=1[t=0] => F>=measS7 per-E, no constant trick)?")
print("=" * 100)

for k in [8, 9, 10]:
    recs = build(k)
    pc = next(p for (E, p) in recs if E == tuple(consec(k)))
    Mc = Mvec(pc); s7c = measS7(pc)
    # Delsarte cone: g(t) = sum_j c_j K_j(t),  need g(t) >= 1[t=0] (t=0..6), c arbitrary sign but g>=...
    # Standard level-1 Delsarte dual: variables g(t) free with g(t)>=1[t=0]; F(E)=sum_t g(t) p_t(E).
    # We work in g-space directly (7 vars g_0..g_6), constraint g_t >= 1[t=0].
    # min s s.t.  F(E)-F(consec) <= s for all E ;  F(consec)=measS7(consec) ;  g_t>=1[t=0].
    # NOTE: validity F(E)>=measS7(E) is AUTOMATIC from g(t)>=1[t=0] (Bonferroni), so we don't add it;
    #       that's the whole point -- validity is per-E structural, NOT a free constraint.
    n = 7
    pcv = np.array([float(x) for x in pc])
    c = [0.0] * n + [1.0]  # vars g_0..g_6, s
    A_ub = []; b_ub = []
    for (E, p) in recs:
        pv = np.array([float(x) for x in p])
        A_ub.append(list(pv - pcv) + [-1.0]); b_ub.append(0.0)   # F(E)-F(consec)<=s
    # g_t >= 1[t=0]  -> -g_t <= -1[t=0]
    for t in range(7):
        row = [0.0] * n + [0.0]; row[t] = -1.0
        A_ub.append(row); b_ub.append(-1.0 if t == 0 else 0.0)
    A_eq = [list(pcv) + [0.0]]; b_eq = [float(s7c)]              # F(consec)=measS7(consec) [tightness]
    res = linprog(c, A_ub=np.array(A_ub), b_ub=np.array(b_ub), A_eq=np.array(A_eq), b_eq=np.array(b_eq),
                  bounds=[(None, None)] * n + [(0.0, None)], method="highs")
    if not res.success:
        print(f"\n  k={k}: Delsarte-cone slack-LP infeasible ({res.message})")
        print(f"        => NO per-E-valid Delsarte functional is tight at consec => consec-max is NOT a")
        print(f"           single Delsarte inequality (genuinely aggregate).")
        continue
    s = res.x[-1]; g = res.x[:n]
    print(f"\n  k={k} (#shapes={len(recs)}): min slack over Delsarte cone (g(t)>=1[t=0], F(consec)=measS7) = {s:.6f}")
    if s < 1e-7:
        print(f"        => consec IS the argmax of a per-E-valid Delsarte functional. g(t)~{[round(x,4) for x in g]}")
        print(f"           CLEAN consec-max via ONE Delsarte inequality (not circular).")
    else:
        # report the beater
        Fc = float(np.dot(g, pcv))
        beaters = []
        for (E, p) in recs:
            pv = np.array([float(x) for x in p])
            fe = float(np.dot(g, pv))
            if fe > Fc + 1e-9: beaters.append((fe, list(E)))
        beaters.sort(reverse=True)
        print(f"        => consec is NOT the argmax of the best tight Delsarte functional (slack {s:.5f}).")
        if beaters[:3]:
            print(f"           top beaters of F over window: {beaters[:3]}")
        print(f"           CONFIRMS: consec-max is genuinely aggregate -- NOT a single Delsarte/Mobius")
        print(f"           inequality.  The per-E-valid (non-circular) cone CANNOT pin consec at degree 1.")

print("\n" + "=" * 100)
print("FINAL THREAD B VERDICT")
print("  - The per-subset (Mobius) BOUND strictly improves the level-1 moment bound (Part B), but")
print("  - the per-subset Mobius EXTREMALITY 'certificate' the slack-LP returns is the TRIVIAL")
print("    constant F=measS7(consec) -- CIRCULAR (validity = the claim).  The genuinely per-subset")
print("    weights collapse to symmetric/constant (Prop 1.2 COLLAPSE for the extremizer).")
print("  - The only NON-circular degree-1 object, the per-E-valid Delsarte cone, is tested here:")
print("    its consec-tight slack tells whether consec-max is ONE Delsarte inequality or aggregate.")
print("=" * 100)
print("\nDONE.")
