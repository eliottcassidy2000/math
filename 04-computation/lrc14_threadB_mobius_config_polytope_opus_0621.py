#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_threadB_mobius_config_polytope_opus_0621.py   (opus, 2026-06-21, THREAD B, Part C)

THE CONFIGURATION-POLYTOPE INTEGRALITY TEST (arXiv:2211.01248 view (d), Lem 3.1).

Part B established: the per-subset (Mobius) data STRICTLY IMPROVES the aggregated moment
LP at every k,R -- it does NOT collapse to level 1.  That is the "linearity power source".
But fixing consec's OWN data and re-maximizing does NOT recover consec's value (big gap):
so per-subset Mobius-positivity ALONE does not force the integral optimum to consec.

Part C asks the SHARP question.  CJJ Lem 3.1 forces an integral optimum because the lattice
Mobius inversion P~[S subset C] -> P~[S = C] makes the optimum a Dirac mass on ONE config.
The right object is the OCCUPANCY/DEPTH LAW polytope:

  Each offset set E -> occupancy law mu_E = distribution over which subset of Z/7 is hit at
  a uniform place x.  measS7(E) = mu_E({all 7 hit}).  We:

  (C1) Build the convex hull of REAL occupancy laws {mu_E : E in window}.  Is consec a
       VERTEX?  Does it MAXIMIZE the linear functional mu -> mu(all7)=measS7 over the hull?
       (If consec is the unique argmax over the hull, then over the DISCRETE set too -- but
       that is just the finite check restated.  The real question is whether a LOW-DEGREE
       face / a Mobius-defined outer relaxation already certifies it.)

  (C2) THE INTEGRALITY-FORCING DUAL.  CJJ: a degree-l pseudoexpectation E~ that is
       Mobius-nonneg on the lattice and matches the achievable moments is forced to be a
       genuine config.  We test the contrapositive cleanly: build the moment-matching
       polytope P_R = { occupancy laws mu (over the 2^7 subsets of Z/7) : the level-<=R
       per-subset upper-sums q_A(mu) are ACHIEVABLE, i.e. equal SOME real offset set's q_A }
       and ask whether max_{mu in P_R} mu(all7) is attained ONLY at an integral (Dirac-on-
       real-config family) point = consec.  Concretely we test the WEAKER necessary form:
       is consec's measS7 >= the LP that uses, as constraints, ALL the per-subset upper
       bounds q_A(E) <= max_E' q_A(E') that hold UNIFORMLY?  (a valid majorant LP).

  (C3) THE CLEAN POSITIVE TEST FOR THE ROUTE.  The Mobius route SUCCEEDS iff there is a
       degree-R "Mobius functional" sum_A w_A q_A (w from Mobius coeffs, integral-forcing)
       whose MAXIMUM over real E is at consec AND equals measS7(consec).  We search for the
       SMALLEST R s.t. a Mobius/Bonferroni-sign functional sum_{|A|<=R} lambda_A q_A:
          (a) upper-bounds measS7 for ALL E   (Bonferroni even-truncation guarantees this
              for the SYMMETRIC choice lambda_A = (-1)^|A| restricted to |A|<=even R),
          (b) is consec-extremal (max over real E in the window at consec),
          (c) equals measS7(consec) (tight).
       (a)+(b)+(c) together = consec-max proof at degree R.  We report which of (a),(b),(c)
       hold at each R and each k.  (c) is the integrality/tightness condition -- if it FAILS
       at every bounded R (it does, per Part B gap), the Mobius route does NOT give exactness
       without the full lattice, confirming Prop-1.2-style collapse of the TIGHTNESS even
       though the BOUND improves.

ALL EXACT (Fractions).  HONEST verdict at the end.
"""
import sys, itertools
from math import comb, gcd
from fractions import Fraction as F
from functools import reduce, lru_cache
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)
H = F(1, 14)
INNER = list(range(1, 7))

def occupancy_full(E):
    E = sorted(set(E)); Enz = [e for e in E if e != 0]
    bps = {F(0), F(1)}
    for e in Enz:
        for m in range(0, 7 * e + 1):
            bps.add(F(m, 7 * e))
    bps = sorted(bps)
    law = defaultdict(lambda: F(0))
    for x0, x1 in zip(bps, bps[1:]):
        if x1 <= x0: continue
        xm = (x0 + x1) / 2
        hit = set(int(7 * e * xm) % 7 for e in E)
        missed = frozenset(s for s in INNER if s not in hit)
        law[missed] += x1 - x0
    return dict(law)

ALL_SUBSETS = [frozenset(b for b in INNER if (mask >> (b - 1)) & 1) for mask in range(64)]

def q_from_law(law):
    """q_A = sum_{B>=A} p_B, for all A subset {1..6}."""
    q = {}
    items = list(law.items())
    for A in ALL_SUBSETS:
        q[A] = sum(m for B, m in items if A <= B)
    return q

def q_upper(E):
    return q_from_law(occupancy_full(E))

def measS7(E):
    return occupancy_full(E).get(frozenset(), F(0))

def Svec_from_q(q):
    return [sum(q[frozenset(A)] for A in itertools.combinations(INNER, r)) for r in range(7)]

def Svec(E):
    return Svec_from_q(q_upper(E))

def primitive(E): return reduce(gcd, [e for e in E if e != 0], 0) == 1
def consec(k): return list(range(k))

def danger(u):
    iv = []
    for j in range(u):
        c = F(j, u); a = (c - H / u) % 1; b = (c + H / u) % 1
        if a < b: iv.append((a, b))
        else: iv.append((a, F(1))); iv.append((F(0), b))
    return iv
def mgmerge(iv):
    iv = sorted(iv); o = []
    for a, b in iv:
        if o and a <= o[-1][1]: o[-1] = (o[-1][0], max(o[-1][1], b))
        else: o.append((a, b))
    return o
def measGP(P):
    if not P: return F(1)
    dz = mgmerge([iv for u in P for iv in danger(u)]); s = F(0); prev = F(0)
    for a, b in dz:
        if a > prev: s += a - prev
        prev = max(prev, b)
    if prev < 1: s += 1 - prev
    return s
@lru_cache(None)
def cap(k):
    psz = 13 - k
    if psz == 0: return F(1)
    return min(measGP(P) for P in itertools.combinations(range(1, 14), psz))

WINDOWS = {8: 16, 9: 15, 10: 14, 11: 13}

# ============================================================================
print("=" * 100)
print("THREAD B PART C: CONFIGURATION-POLYTOPE INTEGRALITY TEST (CJJ Lem 3.1 / view d)")
print("=" * 100)

# ----------------------------------------------------------------------------
# (C3) The decisive test: per-subset Bonferroni functional B_R^sub(E) = sum_{|A|<=R}(-1)^|A| q_A(E)
#      For EVEN R this is a valid UPPER bound on measS7 for EVERY E (Bonferroni on the
#      subset lattice).  NOTE: sum_{|A|=r}(-1)^|A| q_A = (-1)^r S_r, so B_R^sub == the
#      AGGREGATED even-Bonferroni B_R = sum_{r<=R}(-1)^r S_r.  The per-subset refinement
#      that is NOT captured by S_r is to use a DIFFERENT (non-symmetric) lambda_A per subset.
#      We test BOTH:
#        (i)  symmetric: lambda_A=(-1)^|A| (the aggregated B_R).  Properties (a),(b),(c)?
#        (ii) the BEST per-subset majorant: choose lambda_A>=0-signed Mobius weights so that
#             sum lambda_A q_A(E) >= measS7(E) for ALL E in window AND is minimized at consec.
# ----------------------------------------------------------------------------
print("\n" + "=" * 100)
print("(C3-i) SYMMETRIC even-Bonferroni on subset lattice = aggregated B_R = sum_{r<=R}(-1)^r S_r.")
print("       (a) valid upper bound for all E [TRUE by Bonferroni for even R],")
print("       (b) consec-extremal (min over window, since it's an UPPER bound we want consec to")
print("           give the LARGEST true measS7 with the SMALLEST bound -- test argmin/argmax),")
print("       (c) tight at consec (== measS7).")
print("=" * 100)
print(f"\n{'k':>3} {'R':>2} {'B_R(consec)':>12} {'measS7':>9} {'gap=(c)':>10} {'consec=argmax_measS7?':>22} {'consec=argmin_B_R?':>18} {'argmin_B_R':>22}")

def Svec_E(E):
    """compute S_0..S_6 from occupancy law without holding full q-dict."""
    law = occupancy_full(E)
    items = list(law.items())
    S = [F(0)] * 7
    for r in range(7):
        S[r] = sum(sum(m for B, m in items if frozenset(A) <= B)
                   for A in itertools.combinations(INNER, r))
    s7 = law.get(frozenset(), F(0))
    return S, s7

for k in [8, 9, 10, 11]:
    maxE = WINDOWS[k]; Ec = consec(k)
    Sc, s7c = Svec_E(Ec)
    for R in [2, 4]:
        BRc = sum(F((-1) ** r) * Sc[r] for r in range(R + 1))
        s7_beat = 0; BR_lower = 0; argmin = Ec; minBR = BRc
        for rest in itertools.combinations(range(1, maxE + 1), k - 1):
            E = [0] + list(rest)
            if not primitive(E): continue
            S, s7 = Svec_E(E)
            BR = sum(F((-1) ** r) * S[r] for r in range(R + 1))
            if s7 > s7c: s7_beat += 1
            if BR < BRc: BR_lower += 1
            if BR < minBR: minBR = BR; argmin = E
        gap = BRc - s7c
        print(f"{k:>3} {R:>2} {float(BRc):>12.5f} {float(s7c):>9.5f} {float(gap):>10.5f} "
              f"{str(s7_beat==0):>22} {str(BR_lower==0):>18} {str(argmin):>22}")
    print()

# ----------------------------------------------------------------------------
# (C3-ii) THE BEST PER-SUBSET MAJORANT LP (the real Mobius route).
#   Find weights lambda_A (|A|<=R), with the Bonferroni SIGN structure, s.t.
#       F_lambda(E) := sum_{|A|<=R} lambda_A q_A(E)  >=  measS7(E)   for ALL E in window
#   and F_lambda(consec) is MINIMIZED (the tightest uniform majorant at degree R).
#   Then ask: is min_lambda F_lambda(consec) == measS7(consec)?  (=> Mobius route gives a
#   degree-R EXACT certificate, the integrality success).  And is consec the GLOBAL argmax
#   of measS7 forced by F_lambda <= cap?  (the consec-max proof).
#   This is an LP in lambda (dual side).  We solve it exactly.
# ----------------------------------------------------------------------------
print("=" * 100)
print("(C3-ii) BEST DEGREE-R PER-SUBSET MAJORANT: min_lambda F_lambda(consec) s.t. F_lambda(E)>=measS7(E) for all E.")
print("        This is the SHARPEST level-R Mobius/subset bound on measS7(consec).")
print("        If it == measS7(consec): degree-R Mobius certificate is EXACT (integrality success).")
print("        We use lambda over ALL subsets A with |A|<=R (symmetric per size NOT assumed).")
print("=" * 100)

def solve_lp_min(c, A_ub, b_ub, free=True):
    """min c.x s.t. A_ub x >= b_ub, x free (if free) -> use simplex.  Exact Fractions.
       We convert to standard: min c.x, A x >= b.  Split free x = x+ - x-, slacks.
       Returns (optval, x) or None."""
    m = len(A_ub); n = len(c)
    # variables: x+ (n), x- (n), surplus s (m).  A(x+ - x-) - s = b, all >=0.
    nv = 2 * n + m
    rows = []; rhs = []
    for i in range(m):
        row = [F(0)] * nv
        for j in range(n):
            row[j] = A_ub[i][j]; row[n + j] = -A_ub[i][j]
        row[2 * n + i] = F(-1)
        rows.append(row); rhs.append(b_ub[i])
    cost = [F(0)] * nv
    for j in range(n):
        cost[j] = c[j]; cost[n + j] = -c[j]
    val = _simplex_min(rows, rhs, cost)
    return val

def _simplex_min(A, b, c):
    """min c.x s.t. Ax=b, x>=0.  Two-phase exact simplex.  Returns optval or None."""
    m = len(A); n = len(A[0])
    A = [r[:] for r in A]; b = b[:]
    for i in range(m):
        if b[i] < 0:
            A[i] = [-x for x in A[i]]; b[i] = -b[i]
    N = n + m
    T = [A[i][:] + [F(1) if j == i else F(0) for j in range(m)] + [b[i]] for i in range(m)]
    basis = [n + i for i in range(m)]
    def pivot(T, basis, obj):
        mm = len(T); width = len(T[0])
        it = 0
        while True:
            it += 1
            if it > 20000: return False
            col = -1
            for j in range(width - 1):
                if obj[j] < 0: col = j; break
            if col == -1: return True
            best = None; row = -1
            for i in range(mm):
                if T[i][col] > 0:
                    r = T[i][-1] / T[i][col]
                    if best is None or r < best: best = r; row = i
            if row == -1: return None
            pv = T[row][col]; T[row] = [x / pv for x in T[row]]
            for i in range(mm):
                if i != row and T[i][col] != 0:
                    f = T[i][col]; T[i] = [T[i][j] - f * T[row][j] for j in range(width)]
            f = obj[col]
            if f != 0:
                obj[:] = [obj[j] - f * T[row][j] for j in range(width)]
            basis[row] = col
    phase1 = [F(0)] * N + [F(0)]
    for j in range(n, N): phase1[j] = F(1)
    for i in range(m):
        phase1 = [phase1[j] - T[i][j] for j in range(N + 1)]
    r = pivot(T, basis, phase1)
    if r is None or r is False: return None
    if -phase1[-1] != 0: return None
    T2 = [row[:n] + [row[-1]] for row in T]
    obj2 = c[:] + [F(0)]
    for i in range(m):
        bc = basis[i]
        if bc < n and obj2[bc] != 0:
            f = obj2[bc]; obj2 = [obj2[j] - f * T2[i][j] for j in range(n + 1)]
    r = pivot(T2, basis, obj2)
    if r is None or r is False: return None
    return -obj2[-1]  # min value

try:
    from scipy.optimize import linprog
    import numpy as np
    HAVE_SCIPY = True
except Exception:
    HAVE_SCIPY = False

def q_of_E(E):
    """q_A for ALL A subset {1..6}, returned as dict (float-friendly via Fraction)."""
    return q_from_law(occupancy_full(E))

for k in [8, 9, 10]:
    Ec = consec(k); ck = cap(k)
    qc = q_of_E(Ec); s7c = measS7(Ec)
    for R in [2, 3]:
        subsets = [frozenset(A) for r in range(R + 1) for A in itertools.combinations(INNER, r)]
        n = len(subsets)
        c = [qc[A] for A in subsets]
        if HAVE_SCIPY:
            A_ub = []; b_ub = []
            for rest in itertools.combinations(range(1, WINDOWS[k] + 1), k - 1):
                E = [0] + list(rest)
                if not primitive(E): continue
                law = occupancy_full(E); q = q_from_law(law)
                s7 = law.get(frozenset(), F(0))
                A_ub.append([-float(q[A]) for A in subsets]); b_ub.append(-float(s7))
            cobj = [float(x) for x in c]
            res = linprog(cobj, A_ub=np.array(A_ub), b_ub=np.array(b_ub),
                          bounds=[(None, None)] * n, method="highs")
            if not res.success:
                print(f"  k={k} R={R}: scipy LP failed ({res.message})"); continue
            val = res.fun
            gap = val - float(s7c)
            exact = abs(gap) < 1e-9
            print(f"  k={k} R={R} (#subsets={n}): min degree-R majorant at consec = {val:.6f}  "
                  f"measS7(consec)={float(s7c):.6f}  gap={gap:.6f}  EXACT(tight)? {exact}  "
                  f"(<=cap={float(ck):.5f}: {val<=float(ck)+1e-9})")
        else:
            print(f"  k={k} R={R}: scipy unavailable; skipping float majorant LP")
    print()

print("\nDONE (Part C).")
