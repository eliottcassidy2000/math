#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_threadB_signed_cut_exact_kpswf6.py   (kind-pasteur, 2026-06-21, THREAD B)

GOAL (HYP-2744 explicit open target): find the EXPLICIT SPARSEST SIGNED CUT in the
Boolean/type (Mobius) basis that certifies LRC consec-max -- or prove it is still circular.

SETUP (CJJ "view d" Mobius inversion, exact rational arithmetic):
  Inner sectors U={1..6}.  For a speed set E:
    q[M] = meas{x : exact missed inner-sector mask = M}        (Boolean atoms)
    a[A] = meas{A subset missed set} = sum_{M superset A} q[M] (cumulative / containment atoms)
    measS7(E) = q[empty] = sum_A (-1)^|A| a[A]                 (inclusion-exclusion = Bonferroni)
  consec-max claim:  measS7(consec_k) >= measS7(E) for all primitive E in the window.

WHY PRIOR CERTIFICATES WERE CIRCULAR.  The earlier opus LP (certificate_existence) found a
functional F_lambda = sum lambda_A q_A with F>=measS7 (validity) and F(consec)=measS7(consec).
But its optimum was the CONSTANT  F = (481/1470) q_empty  (q_empty=1 here is the WHOLE-space
atom a[empty]=1, so F = measS7(consec) identically) -- it just restates the cap value.  And its
validity constraint  F(E) >= measS7(E)  literally feeds measS7 into the LP, so "tight at consec
+ valid" is information-free about WHY consec wins.

THE GENUINE DUAL CERTIFICATE (non-circular).  We want a signed cut
        C(E) = sum_A lambda_A a[A](E)
that expresses the consec DEFICIT as a MANIFESTLY-NONNEGATIVE combination of atoms:
   (D)  C(E) = measS7(consec_k) - measS7(E)     for ALL E in the window          [deficit-identity]
   (P)  C(E) >= 0                               for ALL E in the window          [nonneg => consec-max]
The cut is NON-CIRCULAR exactly when (P) holds STRUCTURALLY -- i.e. each a[A] is itself a measure
(a[A] >= 0 always, by definition), and the lambda_A pattern with the RIGHT SIGNS makes the whole
combination nonneg.  A POSITIVE-coefficient cut (all lambda_A >= 0 on a set of atoms that are
themselves >=0) would be a trivial monotone certificate; the open question is whether a SIGNED
(mixed +/-) low-support cut works and whether the signs follow the residue/Mobius structure.

We pose three exact LPs (rational, via fractions; we use scipy only to get a vertex then
RE-SOLVE the binding system exactly):

  LP-A (deficit-identity feasibility): does there exist lambda with
        sum_A lambda_A a[A](E) = measS7(consec) - measS7(E)   for ALL E (exact linear system)?
     The atoms a[A](E) over the window form a matrix; we ask if the consec-deficit vector is in
     its column space, and find the SPARSEST such lambda (min-support via L1 / greedy).
     If yes -> the deficit is an EXACT atom combination; we then check sign(lambda) structure.

  LP-B (the real certificate): minimize ||lambda||_0 (support) subject to
        (D-tight)  C(consec)=0  (auto), and  C(E) = deficit(E) is NOT imposed pointwise; instead
        (V)  sum_A lambda_A a[A](E) >= measS7(consec) - measS7(E)  for all E   [valid LOWER cut on deficit]
        (P-built-in) we FORBID lambda_empty (no constant/whole-space term) and FORBID using the
                     single atom that equals measS7 -- the certificate must be built from
                     proper higher atoms.  Among feasible lambda, find min support.

  LP-C (uniformity in k): take the support/sign PATTERN found at k=8, FREEZE it (same A's, signs),
     re-fit only magnitudes at k=9,10 -- does the SAME Boolean-type pattern certify all k?

OUTPUT: the explicit signed cut (subsets + signed rational coeffs), its support size (level),
whether signs follow residue/Mobius parity, and whether it is uniform in k.  Honest 'still
circular' if the only feasible cut is the constant / measS7-atom.
"""
import sys, itertools
from math import comb, gcd
from fractions import Fraction as F
from functools import reduce, lru_cache
from collections import defaultdict

sys.stdout.reconfigure(line_buffering=True)
H = F(1, 14)
INNER = list(range(1, 7))
# subsets of the 6 inner sectors, indexed by mask
SUBMASKS = list(range(64))
def mask_size(m): return bin(m).count("1")
def mask_set(m): return frozenset(s for s in INNER if (m >> (s-1)) & 1)

def exact_mask_atoms(E):
    """q[M] for M a mask over inner sectors 1..6 (bit s-1)."""
    E = sorted(set(E))
    bps = {F(0), F(1)}
    for e in E:
        if e == 0: continue
        for m in range(7*e+1): bps.add(F(m, 7*e))
    bps = sorted(bps)
    q = defaultdict(F)
    for a, b in zip(bps, bps[1:]):
        if a == b: continue
        mid = (a + b) / 2
        hit = {int(7*e*mid) % 7 for e in E}
        mask = 0
        for s in range(1, 7):
            if s not in hit: mask |= 1 << (s-1)
        q[mask] += b - a
    return dict(q)

def containment_atoms(q):
    """a[A] = sum_{M superset A} q[M], indexed by mask A."""
    a = {}
    for A in SUBMASKS:
        a[A] = sum(v for M, v in q.items() if (M & A) == A)
    return a

def measS7(q): return q.get(0, F(0))   # q[empty mask]

def primitive(E): return reduce(gcd, [e for e in E if e != 0], 0) == 1
def full_residue(E): return len({e % 7 for e in E}) == 7   # HYP-2749 full stratum

WINDOWS = {8: 14, 9: 14, 10: 14}   # E={0}+ (k-1)-subset of [1..WINDOW]; span<=14 captures binding cases

def build_stratum(k, window):
    """All primitive E={0}+rest hitting all 7 residues (full stratum). Return list of (E, a-atoms-vec, measS7)."""
    recs = []
    consec_rec = None
    for rest in itertools.combinations(range(1, window+1), k-1):
        E = (0,) + rest
        if not primitive(E): continue
        if not full_residue(E): continue   # restrict to full-residue stratum (HYP-2749)
        q = exact_mask_atoms(E)
        a = containment_atoms(q)
        s7 = measS7(q)
        recs.append((E, a, s7))
        if E == tuple(range(k)): consec_rec = (E, a, s7)
    return recs, consec_rec

def fmt(x):
    if x.denominator == 1: return str(x.numerator)
    return f"{x.numerator}/{x.denominator}"

# ---------- exact rational linear algebra ----------
def rref(M, ncols):
    """Reduced row echelon, M list of lists of Fraction. Returns (M, pivot_cols)."""
    M = [row[:] for row in M]
    nrows = len(M)
    piv = []
    r = 0
    for c in range(ncols):
        # find pivot
        pr = None
        for i in range(r, nrows):
            if M[i][c] != 0:
                pr = i; break
        if pr is None: continue
        M[r], M[pr] = M[pr], M[r]
        inv = M[r][c]
        M[r] = [x / inv for x in M[r]]
        for i in range(nrows):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [a - f*b for a, b in zip(M[i], M[r])]
        piv.append(c); r += 1
        if r == nrows: break
    return M, piv

def in_column_space(A_cols, target):
    """Is target (list of Fraction, len = #rows) in span of columns A_cols (list of column-vectors)?
       Returns (bool, coeffs or None) using augmented rref on rows."""
    nrows = len(target)
    ncolsA = len(A_cols)
    # Build augmented matrix [A | target] row-wise: row i = [A_cols[0][i],...,A_cols[ncolsA-1][i], target[i]]
    M = [[A_cols[j][i] for j in range(ncolsA)] + [target[i]] for i in range(nrows)]
    R, piv = rref(M, ncolsA + 1)
    # consistent iff no pivot in last column
    if (ncolsA) in piv:
        return False, None
    # back out a particular solution (free vars = 0)
    sol = [F(0)] * ncolsA
    # map: pivot col c -> the row that has leading 1 there
    for ridx, c in enumerate(piv):
        if c < ncolsA:
            sol[c] = R[ridx][ncolsA]
    return True, sol


print("="*100)
print("THREAD B: EXPLICIT SIGNED CUT in the Boolean-Mobius (a[A]) basis -- consec-max dual certificate")
print("  Atoms a[A](E)=meas{A subset missed}. measS7=q[empty]=sum_A (-1)^|A| a[A].")
print("  Target: sparsest SIGNED lambda with  C(E)=sum_A lambda_A a[A](E) = measS7(consec)-measS7(E),")
print("          C(E)>=0 structurally => consec-max.  Honest check vs circular (constant/measS7-atom).")
print("="*100)

results_by_k = {}
for k in [8, 9, 10]:
    win = WINDOWS[k]
    recs, consec_rec = build_stratum(k, win)
    assert consec_rec is not None, f"consec not in stratum for k={k}"
    Ec, ac, s7c = consec_rec
    nrec = len(recs)
    print(f"\n{'='*92}\nk={k}: full-residue stratum window [1..{win}], #shapes={nrec}, "
          f"measS7(consec)={fmt(s7c)}={float(s7c):.5f}")

    # ---- LP-A: deficit-identity in atom column space ----
    # target vector: deficit(E) = measS7(consec) - measS7(E) over all recs (rows = recs)
    target = [s7c - s7 for (_, _, s7) in recs]
    # columns = a[A] over recs, for A in SUBMASKS. We DROP A=empty (a[empty]=1 constant => the
    # circular whole-space term) to force a non-constant cut. Keep all proper A (size>=1).
    proper = [A for A in SUBMASKS if A != 0]
    A_cols = [[recs[i][1][A] for i in range(nrec)] for A in proper]
    feasible, sol = in_column_space(A_cols, target)
    print(f"  LP-A (deficit in span of PROPER atoms, A!=empty): feasible={feasible}")
    if feasible:
        # report this particular solution's support; then minimize support by trying small-support subsets
        nz = [(proper[j], sol[j]) for j in range(len(proper)) if sol[j] != 0]
        print(f"    particular solution support={len(nz)}")
    results_by_k[k] = dict(recs=recs, consec=consec_rec, proper=proper, A_cols=A_cols,
                            target=target, feasible_A=feasible, sol_A=sol if feasible else None)

print("\nDONE LP-A pass.")
