#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_threadB_k8_3type_verify_kpswf6.py   (kind-pasteur, 2026-06-21, THREAD B)

EXACT VERIFICATION of the cleanest Thread-B deliverable: a 3-TYPE signed cut that certifies
consec-max at k=8 on the full-residue stratum.

The uniformity script found (float): a degree-(<=3) cut in the dihedral RUN-TYPE basis using only
THREE type-atoms certifies consec-max at k=8 with slack 0.  Here we:
  (1) recompute the type atoms t[tau](E) = sum_{runtype(A)=tau} a[A](E) in EXACT Fraction,
  (2) RE-SOLVE the certificate LP restricted to type-support {(1,),(2,),(1,1,1)} exactly enough
      to get rational coeffs, then
  (3) VERIFY, in exact Fraction over ALL 319 stratum shapes, the sandwich
        measS7(E) <= C(E) <= C(consec) = measS7(consec),     C(E)=sum_tau lam_tau t[tau](E)
      i.e. (V) C(E)>=measS7(E) for all E, (T) C(consec)=measS7(consec), (M) C(E)<=C(consec) all E.
  If all three hold exactly => a verified-exact 3-type consec-max certificate for k=8.

NOTE on the basis: a type-atom t[tau] = sum over subsets A with the given cyclic run-type of the
cumulative atom a[A].  This is a permutation/dihedral-symmetrized Mobius atom (CJJ 'view (b)/(d)'
symmetrization).  The cut is therefore an explicit SIGNED (here all-positive) aggregate in the
type basis -- exactly HYP-2744's requested object, at k=8.
"""
import sys, itertools
from math import gcd
from fractions import Fraction as F
from functools import reduce
from collections import defaultdict

sys.stdout.reconfigure(line_buffering=True)
INNER = list(range(1, 7)); SUBMASKS = list(range(64))
def msize(m): return bin(m).count("1")
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

TYPES = [(1,), (2,), (1, 1, 1)]
def typevec(a):  # exact type atoms for the 3 chosen types
    return {tau: sum(a[A] for A in SUBMASKS if runtype(A) == tau) for tau in TYPES}

WIN = 14; k = 8
recs = []; consec = None
for rest in itertools.combinations(range(1, WIN+1), k-1):
    E = (0,) + rest
    if not prim(E) or not fullres(E): continue
    q = exact_mask_atoms(E); a = cont(q); s7 = measS7(q)
    tv = typevec(a)
    recs.append((E, tv, s7))
    if E == tuple(range(k)): consec = (E, tv, s7)
Ec, tc, s7c = consec
print("="*100)
print(f"THREAD B: EXACT 3-type certificate at k=8 (full-residue stratum, #shapes={len(recs)})")
print(f"  types = {TYPES};  measS7(consec)={fmt(s7c)}")
print("="*100)
print(f"  consec type values: " + ", ".join(f"t[{t}]={fmt(tc[t])}" for t in TYPES))

# Solve exact: find lam (3 rationals) with (T) sum lam_t tc[t] = s7c, and minimize the slack
# s = max_E (C(E)-C(consec)) subject to (V) C(E)>=measS7(E).  We do an exact rational LP via
# enumerating candidate bases is overkill; instead we take the float optimum's basis and solve the
# exact 3x3 binding system: (T) + two tight (V) constraints (the binding shapes).  We detect them
# by scipy then verify exactly.
import numpy as np
from scipy.optimize import linprog
tatoms = TYPES; nA = 3
qcv = np.array([float(tc[t]) for t in tatoms])
c = [0.0]*nA + [1.0]; A_ub = []; b_ub = []; rows = []
for (E, tv, s7) in recs:
    av = np.array([float(tv[t]) for t in tatoms])
    A_ub.append(list(av - qcv) + [-1.0]); b_ub.append(0.0); rows.append(('M', E))
    A_ub.append(list(-av) + [0.0]); b_ub.append(-float(s7)); rows.append(('V', E))
A_eq = [list(qcv) + [0.0]]; b_eq = [float(s7c)]
res = linprog(c, A_ub=np.array(A_ub), b_ub=np.array(b_ub), A_eq=np.array(A_eq), b_eq=np.array(b_eq),
              bounds=[(None, None)]*nA + [(0.0, None)], method="highs")
print(f"\n  float LP: slack s = {res.x[-1]:.8f}  lam = {[f'{v:.5f}' for v in res.x[:3]]}")

# identify binding (V) shapes (tight) to set up exact 3x3 system with (T)
lam_f = res.x[:3]
def Cf(tv): return sum(lam_f[i]*float(tv[t]) for i, t in enumerate(tatoms))
binding = []
for (E, tv, s7) in recs:
    if abs(Cf(tv) - float(s7)) < 1e-6:
        binding.append((E, tv, s7))
# pick consec is in (T); choose two binding V shapes that with (T) give a nonsingular 3x3
def solve3(rowsM, rhs):
    # exact 3x3 solve with Fraction
    M = [r[:] for r in rowsM]; b = rhs[:]
    n = 3
    for i in range(n):
        # pivot
        p = next(j for j in range(i, n) if M[j][i] != 0)
        M[i], M[p] = M[p], M[i]; b[i], b[p] = b[p], b[i]
        inv = M[i][i]; M[i] = [x/inv for x in M[i]]; b[i] = b[i]/inv
        for j in range(n):
            if j != i and M[j][i] != 0:
                f = M[j][i]; M[j] = [a-f*c for a, c in zip(M[j], M[i])]; b[j] = b[j]-f*b[i]
    return b

found = None
for (E1, tv1, s71), (E2, tv2, s72) in itertools.combinations(binding, 2):
    rowsM = [[tc[t] for t in tatoms], [tv1[t] for t in tatoms], [tv2[t] for t in tatoms]]
    rhs = [s7c, s71, s72]
    # check nonsingular
    det = (rowsM[0][0]*(rowsM[1][1]*rowsM[2][2]-rowsM[1][2]*rowsM[2][1])
           - rowsM[0][1]*(rowsM[1][0]*rowsM[2][2]-rowsM[1][2]*rowsM[2][0])
           + rowsM[0][2]*(rowsM[1][0]*rowsM[2][1]-rowsM[1][1]*rowsM[2][0]))
    if det == 0: continue
    lam = solve3(rowsM, rhs)
    # verify exact over ALL shapes
    def C(tv): return sum(lam[i]*tv[t] for i, t in enumerate(tatoms))
    Cc = C(tc)
    okV = all(C(tv) >= s7 for (_, tv, s7) in recs)
    okT = (Cc == s7c)
    okM = all(C(tv) <= Cc for (_, tv, _) in recs)
    if okV and okT and okM:
        found = (lam, E1, E2); break

if found:
    lam, E1, E2 = found
    print("\n  EXACT 3-type certificate FOUND and VERIFIED over all 319 stratum shapes:")
    for i, t in enumerate(tatoms):
        print(f"    lambda[{t}] = {fmt(lam[i])}")
    print(f"    binding shapes (besides consec): {E1}, {E2}")
    Cc = sum(lam[i]*tc[t] for i, t in enumerate(tatoms))
    print(f"    C(consec) = {fmt(Cc)} = measS7(consec) = {fmt(s7c)}  [TIGHT: {Cc==s7c}]")
    print("    (V) C(E)>=measS7(E) all E: VERIFIED   (M) C(E)<=C(consec) all E: VERIFIED")
    print("    => measS7(E) <= C(E) <= C(consec) = measS7(consec) for ALL stratum E  => CONSEC-MAX (k=8).")
else:
    print("\n  No exact 3-type certificate from binding pairs (float slack may be >0).")
print("\nDONE.")
