#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_threadB_signed_cut_levelcap_kpswf6.py   (kind-pasteur, 2026-06-21, THREAD B)

THE DECISIVE LEVEL-CAPPED test, resolving the circularity question precisely.

Established so far:
  * Even-Bonferroni truncations B_R = sum_{|A|<=R} (-1)^|A| a[A] are STRUCTURAL upper bounds on
    measS7 (R even).  consec is the ARGMAX of B_R over the full-residue stratum for EVERY even R.
    BUT B_R(consec) = measS7(consec) only at R=6 (full); for R<6 the bound is LOOSE
    (B_4(consec)~0.48 >> measS7(consec)~0.327), so "consec=argmax B_4" proves only measS7<=0.48,
    NOT consec-max.  At R=6 the bound is exactly measS7 (tautology).  => pure truncation is loose.

  * The LP (proper atoms, no constant) found a VALID+TIGHT+CONSEC-MAX cut with s=0.  We now ask
    the precise question that separates genuine from circular:

    *** What is the MINIMUM LEVEL R (max atom size used) of a cut that is simultaneously
        (V) valid upper bound, (T) tight at consec, (M) consec-maximal? ***

    If the answer is R<6 (a cut using only atoms of size <= R<6, no full-support atom) -> the
    certificate is GENUINELY low-level: it does NOT reconstruct measS7 = B_6 (which needs the
    size-6 atom).  That is a true low-degree dual certificate.
    If the answer is R=6 (the size-6 atom a[{1..6}] is FORCED with the (-1)^6 coefficient) -> the
    only tight+valid+maximal cut IS the full inclusion-exclusion = measS7 itself => circular.

We solve, for each level cap R in {2,3,4,5,6}: LP for min slack s using ONLY atoms with |A|<=R
(NO constant a[empty]).  The smallest R giving s=0 is the certificate level.  We then EXTRACT the
exact cut at that level (rational), report support, signs vs (-1)^|A|, and freeze-and-refit for
uniformity in k.
"""
import sys, itertools
from math import comb, gcd
from fractions import Fraction as F
from functools import reduce
from collections import defaultdict
import numpy as np
from scipy.optimize import linprog

sys.stdout.reconfigure(line_buffering=True)
INNER = list(range(1, 7)); SUBMASKS = list(range(64))
def msize(m): return bin(m).count("1")
def mset(m): return tuple(s for s in INNER if (m >> (s-1)) & 1)

def exact_mask_atoms(E):
    E = sorted(set(E)); bps = {F(0), F(1)}
    for e in E:
        if e == 0: continue
        for m in range(7*e+1): bps.add(F(m, 7*e))
    bps = sorted(bps); q = defaultdict(F)
    for a, b in zip(bps, bps[1:]):
        if a == b: continue
        mid = (a+b)/2; hit = {int(7*e*mid) % 7 for e in E}; mask = 0
        for s in range(1, 7):
            if s not in hit: mask |= 1 << (s-1)
        q[mask] += b - a
    return dict(q)
def cont(q): return {A: sum(v for M, v in q.items() if (M & A) == A) for A in SUBMASKS}
def prim(E): return reduce(gcd, [e for e in E if e], 0) == 1
def fullres(E): return len({e % 7 for e in E}) == 7
def measS7(q): return q.get(0, F(0))
def fmt(x): return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"

WIN = 14
def build(k):
    recs = []; consec = None
    for rest in itertools.combinations(range(1, WIN+1), k-1):
        E = (0,) + rest
        if not prim(E) or not fullres(E): continue
        q = exact_mask_atoms(E); a = cont(q); s7 = measS7(q)
        recs.append((E, a, s7))
        if E == tuple(range(k)): consec = (E, a, s7)
    return recs, consec

def solve_levelcap(recs, consec, R):
    """min slack s; atoms = proper atoms with |A|<=R (no constant). Returns (s, lam_float, atoms)."""
    Ec, ac, s7c = consec
    atoms = [A for A in SUBMASKS if A != 0 and msize(A) <= R]
    nA = len(atoms)
    qcv = np.array([float(ac[A]) for A in atoms])
    c = [0.0]*nA + [1.0]
    A_ub = []; b_ub = []
    for (E, a, s7) in recs:
        av = np.array([float(a[A]) for A in atoms])
        A_ub.append(list(av - qcv) + [-1.0]); b_ub.append(0.0)
        A_ub.append(list(-av) + [0.0]); b_ub.append(-float(s7))
    A_eq = [list(qcv) + [0.0]]; b_eq = [float(s7c)]
    bounds = [(None, None)]*nA + [(0.0, None)]
    res = linprog(c, A_ub=np.array(A_ub), b_ub=np.array(b_ub),
                  A_eq=np.array(A_eq), b_eq=np.array(b_eq), bounds=bounds, method="highs")
    if not res.success: return None, None, atoms
    return res.x[-1], {atoms[j]: res.x[j] for j in range(nA)}, atoms

print("="*100)
print("THREAD B: MINIMUM LEVEL of a valid+tight+consec-max signed cut (proper atoms, no constant)")
print("  R<6 with s=0 => genuine low-level cert (does NOT reconstruct full measS7=B_6).")
print("  Only R=6 with s=0 => circular (the size-6 atom = full inclusion-exclusion is forced).")
print("="*100)

for k in [8, 9, 10]:
    recs, consec = build(k)
    print(f"\nk={k}: #shapes={len(recs)} measS7(consec)={fmt(consec[2])}={float(consec[2]):.5f}")
    first_R = None
    for R in [2, 3, 4, 5, 6]:
        s, lam, atoms = solve_levelcap(recs, consec, R)
        ok = (s is not None and s < 1e-7)
        # report whether the size-6 atom (full mask 63) is used (only relevant at R>=6)
        uses_full = (lam is not None and 63 in lam and abs(lam.get(63, 0.0)) > 1e-9)
        print(f"  level R={R} (#atoms={len(atoms)}): min slack s = "
              f"{('infeasible' if s is None else f'{s:.7f}')}  cert? {ok}"
              + (f"  uses_size6_atom={uses_full}" if R >= 6 else ""))
        if ok and first_R is None:
            first_R = R
    print(f"  => MINIMUM CERTIFICATE LEVEL for k={k}: R = {first_R}")

print("\n" + "="*100)
print("VERDICT: if MIN level < 6 for all k => genuine low-level signed cut (route SUCCEEDS, non-circular).")
print("  if MIN level = 6 => the full-support atom is forced => the tight+valid+max cut IS measS7 => CIRCULAR.")
print("="*100)
print("\nDONE.")
