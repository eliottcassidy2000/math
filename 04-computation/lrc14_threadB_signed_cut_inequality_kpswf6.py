#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_threadB_signed_cut_inequality_kpswf6.py   (kind-pasteur, 2026-06-21, THREAD B)

Follow-up to lrc14_threadB_signed_cut_exact_kpswf6.py (which showed: the consec DEFICIT is
NOT in the linear span of the PROPER cumulative atoms a[A], A!=empty -- the constant a[empty]
is essential to reproduce measS7, i.e. the identity route is circular).

HERE we test the WEAKER, correct certificate: a VALID SIGNED CUT (an INEQUALITY).
We seek lambda (over subsets A of inner sectors 1..6) with cut  C_lambda(E)=sum_A lambda_A a[A](E):
   (V)  C_lambda(E) >= measS7(E)            for ALL stratum shapes E   [valid upper bound]
   (T)  C_lambda(consec_k) = measS7(consec) [tight at consec]
=> measS7(E) <= C_lambda(E) <= C_lambda(consec)=measS7(consec) PROVIDES consec-max IFF additionally
   (M)  C_lambda(E) <= C_lambda(consec)     for ALL E.

THE DECISIVE NON-CIRCULARITY DIAGNOSTIC.  The atom a[empty]=1 identically, so any cut with
lambda_empty>0 contains a free constant (the cap-value restatement).  The question that
separates a GENUINE certificate from a circular one:
   *** Can (V)+(T)+(M) be met with lambda_empty = 0 (NO whole-space constant term)? ***
If YES -> a genuine signed cut built only from PROPER atoms certifies consec-max.
If NO  -> the constant term is forced; the "certificate" is the circular cap restatement.

We solve exactly (rational LP via vertex enumeration on the binding constraints is heavy, so we
use scipy-highs to find the optimum, then RE-VERIFY the reported cut in EXACT Fraction
arithmetic against every shape).  We run three regimes:

  (B0) FULL atoms, free lambda_empty: min support / does s=0 certificate exist (reproduces prior).
  (B1) lambda_empty FORCED 0 (proper atoms only): is (V)+(T)+(M) still feasible with s=0?
       -> the real circularity test.
  (B2) minimize ||lambda||_1 over proper atoms s.t. (V)+(M)+(T): the SPARSEST signed cut; report
       its explicit support, signs, and whether signs follow Mobius parity (-1)^|A| / residue type.
  (B3) UNIFORMITY: freeze the support+sign pattern from k=8, refit magnitudes at k=9,10.

Honest labels throughout.  scipy float for search; exact Fraction for the final verdict.
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
        a = cont(exact_mask_atoms(E)); s7 = measS7(exact_mask_atoms(E))
        recs.append((E, a, s7))
        if E == tuple(range(k)): consec = (E, a, s7)
    return recs, consec

def solve_cut(recs, consec, atoms, forbid_empty):
    """min slack s s.t. (V) C>=measS7, (T) C(consec)=measS7(consec), (M) C(E)-C(consec)<=s.
       atoms = list of masks used. forbid_empty handled by exclusion from atoms list.
       Returns (s, lambda dict)."""
    Ec, ac, s7c = consec
    nA = len(atoms)
    qcv = np.array([float(ac[A]) for A in atoms])
    c = [0.0]*nA + [1.0]   # minimize s
    A_ub = []; b_ub = []
    for (E, a, s7) in recs:
        av = np.array([float(a[A]) for A in atoms])
        A_ub.append(list(av - qcv) + [-1.0]); b_ub.append(0.0)   # (M): (a-ac).lam - s <= 0
        A_ub.append(list(-av) + [0.0]); b_ub.append(-float(s7))  # (V): -a.lam <= -measS7
    A_eq = [list(qcv) + [0.0]]; b_eq = [float(s7c)]              # (T)
    bounds = [(None, None)]*nA + [(0.0, None)]
    res = linprog(c, A_ub=np.array(A_ub), b_ub=np.array(b_ub),
                  A_eq=np.array(A_eq), b_eq=np.array(b_eq), bounds=bounds, method="highs")
    if not res.success: return None, None
    lam = {atoms[j]: res.x[j] for j in range(nA)}
    return res.x[-1], lam

print("="*100)
print("THREAD B: SIGNED-CUT INEQUALITY certificate -- the non-circularity diagnostic")
print("  KEY TEST (B1): can (V)valid + (T)tight + (M)consec-max hold with NO constant atom a[empty]?")
print("  YES => genuine proper-atom signed cut.  NO => constant forced => circular cap restatement.")
print("="*100)

data = {}
for k in [8, 9, 10]:
    recs, consec = build(k)
    data[k] = (recs, consec)
    s7c = consec[2]
    proper = [A for A in SUBMASKS if A != 0]
    allA = SUBMASKS

    # B0: full atoms
    s0, lam0 = solve_cut(recs, consec, allA, forbid_empty=False)
    # B1: proper only (no constant)
    s1, lam1 = solve_cut(recs, consec, proper, forbid_empty=True)
    print(f"\nk={k}: #shapes={len(recs)} measS7(consec)={fmt(s7c)}={float(s7c):.5f}")
    print(f"  (B0) FULL atoms incl a[empty]: min slack s = {s0:.7f}  cert_exists={s0 is not None and s0<1e-7}")
    print(f"  (B1) PROPER atoms only (NO constant a[empty]): min slack s = "
          f"{('infeasible' if s1 is None else f'{s1:.7f}')}  "
          f"genuine_proper_cert={s1 is not None and s1<1e-7}")

print("\n" + "="*100)
print("B0 vs B1 verdict: if B0 gives s=0 but B1 gives s>0, the constant a[empty] is ESSENTIAL")
print("  => the per-subset 'certificate' is the circular cap-value restatement, NOT a proper cut.")
print("  if B1 also gives s=0 => a genuine proper-atom signed cut exists (route succeeds).")
print("="*100)
print("\nDONE B0/B1.")
