#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_threadB_k8_3type_exactLP_kpswf6.py   (kind-pasteur, 2026-06-21, THREAD B)

Exact rational feasibility of the 3-type consec-max certificate at k=8.
With slack s=0 the certificate conditions are (3 unknowns lam over types [(1,),(2,),(1,1,1)]):
   (T)  sum_t lam_t * tc[t]   = measS7(consec)
   (V)  sum_t lam_t * tv[t](E) >= measS7(E)        for all E
   (M)  sum_t lam_t * tv[t](E) <= sum_t lam_t*tc[t] for all E      (= C(E) <= C(consec))
We solve this exact rational feasibility by 3-variable vertex enumeration:
  use (T) to eliminate one variable, leaving a 2D polytope in (lam1,lam2); each remaining
  constraint is a half-plane; enumerate all vertices (intersections of pairs of boundary lines)
  with exact Fraction, test feasibility of each vertex, return the first feasible one.
This is a rigorous EXACT certificate search (no float).
"""
import sys, itertools
from math import gcd
from fractions import Fraction as F
from functools import reduce
from collections import defaultdict

sys.stdout.reconfigure(line_buffering=True)
INNER = list(range(1, 7)); SUBMASKS = list(range(64))
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
WIN = 14; k = 8
recs = []; consec = None
for rest in itertools.combinations(range(1, WIN+1), k-1):
    E = (0,) + rest
    if not prim(E) or not fullres(E): continue
    a = cont(exact_mask_atoms(E)); s7 = measS7(exact_mask_atoms(E))
    tv = [sum(a[A] for A in SUBMASKS if runtype(A) == t) for t in TYPES]
    recs.append((E, tv, s7))
    if E == tuple(range(k)): consec = (E, tv, s7)
Ec, tc, s7c = consec

print("="*100)
print(f"THREAD B: EXACT rational 3-type certificate at k=8  (#shapes={len(recs)}, types={TYPES})")
print("="*100)

# Eliminate lam2 (coeff of TYPES[2]) using (T): lam0*tc0 + lam1*tc1 + lam2*tc2 = s7c
# => lam2 = (s7c - lam0*tc0 - lam1*tc1)/tc2
tc0, tc1, tc2 = tc
assert tc2 != 0
# Each constraint a row over (lam0,lam1,lam2): coefficient vector w, threshold.
# (V): w.lam >= m   (M): w.lam <= C(consec)=s7c (since C(consec)=s7c by (T))
# Substitute lam2 -> express each as alpha*lam0 + beta*lam1 >= / <= rhs  in 2 vars.
def sub(w):  # w=(w0,w1,w2) -> linear in (lam0,lam1): A*lam0 + B*lam1 + const
    A = w[0] - w[2]*tc0/tc2
    B = w[1] - w[2]*tc1/tc2
    const = w[2]*s7c/tc2
    return A, B, const   # value = A*lam0 + B*lam1 + const

halfplanes = []  # each (A,B,c): A*lam0+B*lam1 >= c   (we store as >= form)
for (E, tv, s7) in recs:
    A, B, const = sub(tv)
    # (V): tv.lam >= s7  => A*lam0+B*lam1+const >= s7 => A*lam0+B*lam1 >= s7-const
    halfplanes.append((A, B, s7 - const))
    # (M): tv.lam <= s7c => A*lam0+B*lam1+const <= s7c => -A*lam0-B*lam1 >= -(s7c-const)
    halfplanes.append((-A, -B, -(s7c - const)))

# vertex enumeration: intersect each pair of boundary lines (A*x+B*y=c), test all halfplanes
def line_int(h1, h2):
    A1, B1, C1 = h1; A2, B2, C2 = h2
    det = A1*B2 - A2*B1
    if det == 0: return None
    x = (C1*B2 - C2*B1)/det
    y = (A1*C2 - A2*C1)/det
    return (x, y)

def feasible(pt):
    x, y = pt
    return all(A*x + B*y >= c for (A, B, c) in halfplanes)

found = None
n = len(halfplanes)
for i in range(n):
    for j in range(i+1, n):
        pt = line_int(halfplanes[i], halfplanes[j])
        if pt is None: continue
        if feasible(pt):
            found = pt; break
    if found: break

if found:
    lam0, lam1 = found
    lam2 = (s7c - lam0*tc0 - lam1*tc1)/tc2
    lam = [lam0, lam1, lam2]
    # final exact verification
    def C(tv): return sum(lam[i]*tv[i] for i in range(3))
    Cc = C(tc)
    okT = (Cc == s7c)
    okV = all(C(tv) >= s7 for (_, tv, s7) in recs)
    okM = all(C(tv) <= Cc for (_, tv, _) in recs)
    print("\n  EXACT 3-type certificate FOUND (rational vertex), verified over ALL stratum shapes:")
    for i, t in enumerate(TYPES):
        print(f"    lambda[{t}] = {fmt(lam[i])}  = {float(lam[i]):.6f}")
    print(f"    C(consec) = {fmt(Cc)} ; measS7(consec) = {fmt(s7c)} ; TIGHT: {okT}")
    print(f"    (V) C(E)>=measS7(E) all E: {okV}")
    print(f"    (M) C(E)<=C(consec) all E: {okM}")
    if okT and okV and okM:
        print("\n    *** VERIFIED EXACT: measS7(E) <= C(E) <= C(consec)=measS7(consec) for all 319 stratum E.")
        print("        C(E) = lam0*t[(1,)] + lam1*t[(2,)] + lam2*t[(1,1,1)],  t[tau]=sum_{runtype(A)=tau} a[A].")
        print("        This is an explicit SIGNED AGGREGATE CUT in the dihedral run-type (Mobius) basis,")
        print("        level 3 (largest type has 3 sectors), support 3, certifying consec-max at k=8.")
else:
    print("\n  No feasible rational vertex among pairwise line intersections (degenerate/unbounded).")
print("\nDONE.")
