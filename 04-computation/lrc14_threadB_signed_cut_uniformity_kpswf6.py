#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_threadB_signed_cut_uniformity_kpswf6.py   (kind-pasteur, 2026-06-21, THREAD B)

FINAL Thread-B diagnostics: uniformity and the TYPE-basis (residue/dihedral) reduction.

(U1) UNIFORMITY by FROZEN SUPPORT: take the k=8 cut's support set S (the atoms it uses), FREEZE
     exactly that atom set, and re-solve the min-slack LP at k=9,10 using ONLY atoms in S.
     If s=0 at all k with the SAME support => uniform support certificate.  If s>0 => support
     must grow (matches the levelcap finding R=2->R=3).

(U2) TYPE-BASIS cut: the dihedral run-type quotient (HYP-2744 'type basis', HYP-2745 residue
     structure).  Collapse atoms a[A] into TYPE atoms  t[tau] = sum_{runtype(A)=tau} a[A].
     Solve min-slack LP in the TYPE basis (one variable per run-type).  Sparser? does it certify?
     This tests whether the certificate is genuinely a RESIDUE/type object (quasimodular) or needs
     individual-subset resolution.

(U3) the honest VALIDITY-HARDNESS probe.  The cut C has C(E)>=measS7(E) on the stratum.  We test
     whether (V) is a CONSEQUENCE of a FINITE set of structural atom inequalities, by checking if
     C - measS7 = sum_A beta_A a[A] with the NEGATIVE-beta atoms restricted to a SMALL controllable
     set.  We report the negative-beta SUPPORT size and the total negative mass; small => the
     validity gap is governed by few atoms (tractable), large => validity is as hard as consec-max.
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
def runtype(mask):
    bits = [(mask >> i) & 1 for i in range(6)]
    if sum(bits) == 0: return ()
    if all(bits): return (6,)
    cands = []
    for seq in (bits, list(reversed(bits))):
        for sh in range(6):
            row = seq[sh:] + seq[:sh]
            if row[-1] == 0 and row[0] == 1:
                lens = []; i = 0
                while i < 6:
                    if row[i]:
                        j = i
                        while j < 6 and row[j]: j += 1
                        lens.append(j-i); i = j
                    else: i += 1
                cands.append(tuple(lens))
    return min(cands)

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
        recs.append((E, a, s7));
        if E == tuple(range(k)): consec = (E, a, s7)
    return recs, consec

def solve_minslack(recs, consec, atoms):
    Ec, ac, s7c = consec; nA = len(atoms)
    qcv = np.array([float(ac[A]) for A in atoms])
    c = [0.0]*nA + [1.0]; A_ub = []; b_ub = []
    for (E, a, s7) in recs:
        av = np.array([float(a[A]) for A in atoms])
        A_ub.append(list(av - qcv) + [-1.0]); b_ub.append(0.0)
        A_ub.append(list(-av) + [0.0]); b_ub.append(-float(s7))
    A_eq = [list(qcv) + [0.0]]; b_eq = [float(s7c)]
    bounds = [(None, None)]*nA + [(0.0, None)]
    res = linprog(c, A_ub=np.array(A_ub), b_ub=np.array(b_ub),
                  A_eq=np.array(A_eq), b_eq=np.array(b_eq), bounds=bounds, method="highs")
    if not res.success: return None, None
    return res.x[-1], {atoms[j]: res.x[j] for j in range(nA)}

# get k=8 R=2 cut support
def min_l1_support(recs, consec, R):
    Ec, ac, s7c = consec
    atoms = [A for A in SUBMASKS if A != 0 and msize(A) <= R]; nA = len(atoms)
    qcv = np.array([float(ac[A]) for A in atoms])
    c = [1.0]*(2*nA); A_ub = []; b_ub = []
    for (E, a, s7) in recs:
        av = np.array([float(a[A]) for A in atoms])
        A_ub.append(list(-av) + list(av)); b_ub.append(-float(s7))
        A_ub.append(list(av - qcv) + list(-(av - qcv))); b_ub.append(0.0)
    A_eq = [list(qcv) + list(-qcv)]; b_eq = [float(s7c)]
    res = linprog(c, A_ub=np.array(A_ub), b_ub=np.array(b_ub),
                  A_eq=np.array(A_eq), b_eq=np.array(b_eq), bounds=[(0.0, None)]*(2*nA), method="highs")
    lam = {atoms[j]: res.x[j]-res.x[nA+j] for j in range(nA)}
    return [A for A in atoms if abs(lam[A]) > 1e-7]

print("="*100)
print("THREAD B FINAL: uniformity + type-basis reduction + validity-hardness")
print("="*100)

D = {k: build(k) for k in [8, 9, 10]}

# (U1) frozen support
print("\n(U1) FROZEN SUPPORT uniformity: take k=8 R=2 support, reuse it at k=9,10")
recs8, consec8 = D[8]
S8 = min_l1_support(recs8, consec8, 2)
print(f"  k=8 support S8 = {[mset(A) for A in S8]}  (|S8|={len(S8)})")
for k in [9, 10]:
    recs, consec = D[k]
    s, lam = solve_minslack(recs, consec, S8)
    print(f"  k={k} with FROZEN S8: min slack s = {('infeasible' if s is None else f'{s:.7f}')}  "
          f"uniform_with_S8={s is not None and s < 1e-7}")

# (U2) type basis
print("\n(U2) TYPE-BASIS cut: collapse a[A] into run-type atoms t[tau]=sum_{runtype(A)=tau} a[A]")
all_types = sorted({runtype(A) for A in SUBMASKS}, key=lambda t: (sum(t), t))
print(f"  #run-types (incl empty,full) = {len(all_types)}: {all_types}")
def type_recs(recs):
    out = []
    for (E, a, s7) in recs:
        ta = {tau: sum(a[A] for A in SUBMASKS if runtype(A) == tau) for tau in all_types}
        out.append((E, ta, s7))
    return out
for k in [8, 9, 10]:
    recs, consec = D[k]
    trecs = type_recs(recs)
    Ec, ac, s7c = consec
    tac = {tau: sum(ac[A] for A in SUBMASKS if runtype(A) == tau) for tau in all_types}
    tconsec = (Ec, tac, s7c)
    # exclude empty type () which = constant a[empty]
    tatoms = [t for t in all_types if t != ()]
    s, lam = solve_minslack(trecs, tconsec, tatoms)
    nz = {t: v for t, v in (lam or {}).items() if abs(v) > 1e-7}
    print(f"  k={k}: type-basis min slack s = {('infeasible' if s is None else f'{s:.7f}')}  "
          f"cert={s is not None and s<1e-7}  type-support={len(nz)} of {len(tatoms)}")
    if s is not None and s < 1e-7:
        for t in sorted(nz, key=lambda x: (sum(x), x)):
            print(f"      type {t}: coeff {nz[t]:+.5f}")

print("\n" + "="*100)
print("VERDICT logic:")
print("  (U1) frozen-S8 s=0 at k=9,10 => UNIFORM support; s>0 => support grows (non-uniform).")
print("  (U2) type-basis s=0 with few types => certificate is a RESIDUE/dihedral-type object")
print("       (HYP-2745 quasimodular); s>0 => type quotient too coarse, needs individual subsets.")
print("="*100)
print("\nDONE.")
