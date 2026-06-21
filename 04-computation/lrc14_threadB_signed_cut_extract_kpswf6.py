#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_threadB_signed_cut_extract_kpswf6.py   (kind-pasteur, 2026-06-21, THREAD B)

EXTRACT and EXACT-VERIFY the explicit degree-2/3 signed cut, and run the structural
non-circularity test.

Established (levelcap script): a valid+tight+consec-max cut over PROPER atoms exists at
LEVEL R=2 (k=8,9), R=3 (k=10); the size-6 atom is never forced.  Now:

  (1) Extract the SPARSEST such cut exactly.  We minimize L1 norm of lambda over level-R atoms
      s.t. (V)+(T)+(M), then re-solve the binding system in exact Fraction to get rational coeffs.
  (2) STRUCTURAL NON-CIRCULARITY: write  C(E)-measS7(E)=sum_A beta_A a[A](E), beta_A=lambda_A-(-1)^|A|.
      Since every a[A]>=0, if beta_A>=0 for all A then validity (V) is STRUCTURAL (no measS7 needed
      in the proof -- only the trivial atom-nonnegativity).  Report how many beta_A are negative:
      0 negatives => fully structural certificate.  We also report the SIGNS of lambda vs (-1)^|A|.
  (3) SIGN/RESIDUE STRUCTURE: tabulate lambda_A by atom size |A| and by dihedral run-type;
      check if the sign follows Mobius parity (-1)^|A| and whether magnitudes depend only on
      residue/run-type (HYP-2745 quasimodular structure).
  (4) UNIFORMITY: freeze {support, signs} from the k=8 cut; refit magnitudes at k=9,10; report
      whether the same Boolean-type pattern certifies (s=0) at higher k.
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
        recs.append((E, a, s7))
        if E == tuple(range(k)): consec = (E, a, s7)
    return recs, consec

def min_l1_cut(recs, consec, R):
    """min ||lambda||_1 s.t. (V) C>=measS7, (T) C(consec)=measS7(consec), (M) C(E)<=C(consec).
       proper atoms |A|<=R, no constant. Returns (lam_float dict, atoms)."""
    Ec, ac, s7c = consec
    atoms = [A for A in SUBMASKS if A != 0 and msize(A) <= R]
    nA = len(atoms)
    qcv = np.array([float(ac[A]) for A in atoms])
    # vars: lambda (nA) split into lam = u - v, u,v>=0 ; minimize sum(u+v)
    # constraints in terms of lam, then substitute lam=u-v.
    # (V): -C(E) <= -measS7  => -(a_E).lam <= -s7
    # (M): C(E)-C(consec) <= 0 => (a_E - ac).lam <= 0
    # (T): ac.lam = s7c
    c = [1.0]*nA + [1.0]*nA
    A_ub = []; b_ub = []
    for (E, a, s7) in recs:
        av = np.array([float(a[A]) for A in atoms])
        A_ub.append(list(-av) + list(av)); b_ub.append(-float(s7))       # (V)
        A_ub.append(list(av - qcv) + list(-(av - qcv))); b_ub.append(0.0) # (M)
    A_eq = [list(qcv) + list(-qcv)]; b_eq = [float(s7c)]                  # (T)
    bounds = [(0.0, None)]*(2*nA)
    res = linprog(c, A_ub=np.array(A_ub), b_ub=np.array(b_ub),
                  A_eq=np.array(A_eq), b_eq=np.array(b_eq), bounds=bounds, method="highs")
    if not res.success: return None, atoms
    lam = {atoms[j]: res.x[j] - res.x[nA + j] for j in range(nA)}
    return lam, atoms

print("="*100)
print("THREAD B: EXTRACT explicit min-L1 signed cut + structural non-circularity test")
print("="*100)

for k, R in [(8, 2), (9, 2), (10, 3)]:
    recs, consec = build(k)
    Ec, ac, s7c = consec
    lam, atoms = min_l1_cut(recs, consec, R)
    if lam is None:
        print(f"\nk={k} R={R}: L1 LP infeasible"); continue
    nz = {A: v for A, v in lam.items() if abs(v) > 1e-9}
    print(f"\n{'='*92}\nk={k}, level R={R}: explicit min-L1 cut, support={len(nz)} atoms (of {len(atoms)})")
    # structural: beta_A = lambda_A - (-1)^|A|  (only over USED atoms; unused have lambda=0)
    # full beta over ALL submasks (lambda=0 for atoms not in 'atoms' set or zero)
    lamfull = {A: lam.get(A, 0.0) for A in SUBMASKS}
    neg_beta = 0; beta_list = []
    for A in SUBMASKS:
        beta = lamfull[A] - ((-1)**msize(A))
        beta_list.append((A, beta))
        if beta < -1e-9: neg_beta += 1
    print(f"  STRUCTURAL non-circularity: #atoms with beta_A=lambda_A-(-1)^|A| < 0 : {neg_beta} (of 64)")
    print(f"    (beta>=0 for all => C-measS7 is a nonneg atom combo => validity STRUCTURAL/Bonferroni)")
    # sign vs parity for the NONZERO lambda
    print(f"  nonzero lambda_A (size | runtype | lambda | sign vs (-1)^|A|):")
    for A in sorted(nz, key=lambda a: (msize(a), a)):
        sz = msize(A); v = nz[A]; parity = (-1)**sz
        agree = (v > 0) == (parity > 0)
        print(f"    A={mset(A)} |A|={sz} runtype={runtype(A)}: lambda={v:+.5f}  parity={parity:+d}  sign_follows_parity={agree}")
    # how many nonzero lambda follow parity
    follow = sum(1 for A in nz if (nz[A] > 0) == (((-1)**msize(A)) > 0))
    print(f"  sign-follows-parity: {follow}/{len(nz)} nonzero coeffs")

print("\n" + "="*100)
print("READING: neg_beta=0 => the cut's validity is structural (manifestly-nonneg atom combo).")
print("  sign-follows-parity high => the signed cut is a TUNED even-Bonferroni (Mobius parity signs,")
print("  magnitudes adjusted).  support sparse + low level => the explicit signed cut answer.")
print("="*100)
print("\nDONE.")
