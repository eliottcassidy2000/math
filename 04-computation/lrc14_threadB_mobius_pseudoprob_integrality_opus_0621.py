#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_threadB_mobius_pseudoprob_integrality_opus_0621.py   (opus, 2026-06-21, THREAD B)

THE MOBIUS-INVERSION / PSEUDOPROBABILITY INTEGRALITY ROUTE (arXiv:2211.01248 Sec 3, Lem 3.1).

CJJ's level-l hierarchy improves on Delsarte (level 1) ONLY via LINEARITY of the optimizer
(closed under linear combinations).  For NON-linear optimizers it COLLAPSES to level 1 at
every level (Prop 1.2).  View (d): Mobius inversion on the subspace lattice turns
P~[S subset C] -> P~[S = C], forcing an integral optimum WITHOUT whole-polytope integrality.

OUR TRANSLATION.  The "sectors hit" object lives on the subset lattice of the 6 inner
sectors {1..6} (sector 0 always hit by e=0).  For an offset set E:
   N(x) = #missed inner sectors,  measS7(E) = P(N=0) = p_0.
The factorial moments S_r = E[C(N,r)] = sum_{|A|=r, A subset {1..6}} q_A(E),
   q_A(E) := meas{x : ALL sectors in A are MISSED} = P~[missed >= A]  (UPPER/zeta).
Mobius inversion on the subset lattice (= Bonferroni = our inclusion-exclusion):
   p_B(E) = meas{x : missed set EXACTLY = B} = sum_{A >= B} (-1)^{|A\B|} q_A(E)  (>= 0 always).
   measS7 = p_emptyset = sum_{all A} (-1)^{|A|} q_A.

THE QUESTION (deliverable): reinterpret the LP variables as the UPPER pseudoprobabilities
q_A (one per subset A of {1..6}, 2^6=64 vars), with the Mobius-positivity p_B>=0 (64
constraints).  Does this FORCE the optimal pseudo-distribution to be a GENUINE single
offset configuration (integral) => consec exact, WITHOUT the finite atlas?  Or does it
COLLAPSE to the level-1 moment LP (Prop 1.2, a valuable negative)?

We BUILD and SOLVE the small LPs exactly (Fractions), three nested relaxations:

  LP-MOMENT  (level-1 Delsarte, what's already studied): vars p_t (t=0..6),
             constraints sum C(t,r) p_t = S_r(consec) for r=0..R, max p_0.  Subset
             structure AGGREGATED into level sums.

  LP-SUBSET  (the Mobius route): vars q_A for ALL A subset {1..6} (64), Mobius-positivity
             p_B(q) = sum_{A>=B}(-1)^{|A\B|} q_A >= 0 for all B, q_emptyset=1, and the
             ACHIEVABLE single-subset upper bounds 0<=q_A<=cap-type, max p_emptyset.
             We test: with the EXACT consec q_A as the only data fixed at level r<=R,
             does the per-subset Mobius LP have a SMALLER feasible max than the aggregated
             moment LP?  (If equal -> COLLAPSE per Prop 1.2.  If strictly smaller and the
             max is achieved only by a genuine offset law -> INTEGRALITY route works.)

  LP-CONFIG  (the integral target): vars are pseudoprobabilities over the actual depth-law
             / occupancy configurations of REAL offset sets.  We check whether the convex
             hull of real offset depth-laws, intersected with the moment constraints,
             pins measS7 to consec (i.e. whether consec is a vertex/the maximizer of the
             configuration polytope).

HONEST verdict at the end: COLLAPSES vs IMPROVES, and whether a low level proves consec-max.
"""
import sys, itertools
from math import comb, gcd
from fractions import Fraction as F
from functools import reduce, lru_cache
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)
H = F(1, 14)
INNER = list(range(1, 7))  # the 6 inner sectors; sector 0 always hit

# ----------------------------------------------------------------------------
# Core exact occupancy: per-x which sectors hit, breakpoint pass.
# q_A(E) = meas{ all sectors in A missed }.  p_B(E)=meas{ missed set EXACTLY B }.
# ----------------------------------------------------------------------------
def occupancy_full(E):
    """Return dict B(frozenset of MISSED inner sectors) -> measure, exact."""
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
        hit = set(int(7 * e * xm) % 7 for e in E)  # includes 0
        missed = frozenset(s for s in INNER if s not in hit)
        law[missed] += x1 - x0
    return dict(law)

def q_upper(E):
    """q_A = meas{A subset missed} for ALL A subset {1..6}.  q_A = sum_{B>=A} p_B."""
    law = occupancy_full(E)
    q = {}
    for r in range(7):
        for A in itertools.combinations(INNER, r):
            As = frozenset(A)
            q[As] = sum(m for B, m in law.items() if As <= B)
    return q, law

def Svec_from_q(q):
    S = [F(0)] * 7
    for r in range(7):
        S[r] = sum(q[frozenset(A)] for A in itertools.combinations(INNER, r))
    return S

def measS7(E):
    law = occupancy_full(E)
    return law.get(frozenset(), F(0))

def primitive(E): return reduce(gcd, [e for e in E if e != 0], 0) == 1
def consec(k): return list(range(k))

# ----------------------------------------------------------------------------
# caps (gp danger-zone)
# ----------------------------------------------------------------------------
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

# ----------------------------------------------------------------------------
# exact LP solver via vertex enumeration on a small system.  Returns max of objective.
# We use the moment-LP and a generic "max c.x s.t. Ax=b, x>=0" vertex enumerator.
# ----------------------------------------------------------------------------
def solve_square(rows, b):
    n = len(rows)
    M = [rows[i][:] + [b[i]] for i in range(n)]
    for c in range(n):
        piv = None
        for r in range(c, n):
            if M[r][c] != 0: piv = r; break
        if piv is None: return None
        M[c], M[piv] = M[piv], M[c]
        pv = M[c][c]; M[c] = [x / pv for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]; M[r] = [M[r][j] - f * M[c][j] for j in range(n + 1)]
    return [M[i][n] for i in range(n)]

def moment_lp_max_p0(S, R):
    """max p_0 over prob measures on t=0..6 with C(t,r) moments = S_r, r=0..R."""
    ts = list(range(7)); best = None
    for sz in range(1, R + 2):
        for supp in itertools.combinations(ts, sz):
            rows = [[F(comb(t, r)) for t in supp] for r in range(min(sz, R + 1))]
            b = [S[r] for r in range(min(sz, R + 1))]
            if len(rows) != sz:  # underdetermined; pad with more moment rows up to sz
                extra = sz - len(rows)
                for r in range(min(sz, R + 1), min(sz, R + 1) + extra):
                    if r > 6: break
                    rows.append([F(comb(t, r)) for t in supp]); b.append(S[r])
            if len(rows) != sz: continue
            sol = solve_square(rows, b)
            if sol is None or any(x < 0 for x in sol): continue
            # verify ALL moments r=0..R
            ok = all(sum(sol[i] * F(comb(supp[i], r)) for i in range(sz)) == S[r] for r in range(R + 1))
            if not ok: continue
            p0 = sum(sol[i] for i in range(sz) if supp[i] == 0)
            if best is None or p0 > best: best = p0
    return best

# ============================================================================
print("=" * 100)
print("THREAD B: MOBIUS / PSEUDOPROBABILITY INTEGRALITY ROUTE  (arXiv:2211.01248 Lem 3.1)")
print("=" * 100)

# ----------------------------------------------------------------------------
# (A) The Mobius dictionary made explicit, and exact for consec.
# ----------------------------------------------------------------------------
print("\n" + "=" * 100)
print("(A) MOBIUS DICTIONARY: q_A (UPPER, =zeta) <-> p_B (EXACT, =mobius) on subset lattice of {1..6}")
print("    measS7 = p_emptyset = sum_A (-1)^|A| q_A  (full inclusion-exclusion = Mobius at bottom).")
print("=" * 100)
for k in [8]:
    E = consec(k); q, law = q_upper(E); S = Svec_from_q(q)
    print(f"  consec_{k}:  measS7 = {measS7(E)} = {float(measS7(E)):.5f}")
    # verify Mobius inversion bottom element
    mob0 = sum(F((-1) ** len(A)) * q[frozenset(A)] for r in range(7) for A in itertools.combinations(INNER, r))
    print(f"  sum_A (-1)^|A| q_A = {mob0}   matches p_emptyset: {mob0 == measS7(E)}")
    # the per-subset q_A by size
    print("  q_A by size r (these are the per-subset pseudoprobabilities, NOT just S_r):")
    for r in range(4):
        vals = sorted(set(q[frozenset(A)] for A in itertools.combinations(INNER, r)))
        print(f"    r={r}: distinct q_A values = {[str(v) for v in vals]}  (count of distinct = {len(vals)})  S_{r}={S[r]}")

# ----------------------------------------------------------------------------
# (B) THE COLLAPSE TEST.  Compare:
#     (i)  moment-LP_R(consec)  = max p_0 given AGGREGATED level sums S_0..S_R.
#     (ii) subset-LP_R(consec)  = max p_emptyset given the EXACT per-subset q_A for |A|<=R
#          plus Mobius-positivity p_B>=0 for all 64 subsets B.
#     If subset-LP == moment-LP for all R -> COLLAPSE (Prop 1.2).
#     If subset-LP < moment-LP -> the per-subset (Mobius) data IMPROVES the bound.
# ----------------------------------------------------------------------------
print("\n" + "=" * 100)
print("(B) COLLAPSE TEST: subset-LP (per-subset q_A, |A|<=R, + Mobius p_B>=0) vs moment-LP (level sums).")
print("    Fix consec's exact data up to level R; maximize measS7 over feasible pseudo-laws.")
print("=" * 100)

def subset_lp_max(E, R):
    """max p_emptyset over pseudo-laws p_B>=0 (sum=1) whose UPPER sums q_A = sum_{B>=A} p_B
       match consec's EXACT q_A for ALL A with |A|<=R.  Vertex/LP over the 64-dim p_B simplex
       with the constraints q_A(p)=q_A(E) (|A|<=R).  Solve via LP duality on the small system
       by enumerating: since #constraints = sum_{r<=R} C(6,r), the optimum has support of that
       size.  For tractability we instead solve the LP directly with a simple exact simplex."""
    q, _ = q_upper(E)
    # variables: p_B for all B subset {1..6} (64).  We index by the 6-bit mask.
    Bs = [frozenset(b for b in INNER if (mask >> (b - 1)) & 1) for mask in range(64)]
    # constraint rows: for each A with |A|<=R, sum_{B>=A} p_B = q_A.  (A=emptyset gives sum p=1=q_empty.)
    As = [frozenset(A) for r in range(R + 1) for A in itertools.combinations(INNER, r)]
    rowsA = []
    bA = []
    for A in As:
        row = [F(1) if A <= B else F(0) for B in Bs]
        rowsA.append(row); bA.append(q[A])
    # objective: max p_{emptyset} = coefficient on B=emptyset
    obj = [F(1) if len(B) == 0 else F(0) for B in Bs]
    return exact_lp_max(rowsA, bA, obj)

def exact_lp_max(A, b, c):
    """Exact rational LP: max c.x s.t. A x = b, x>=0.  Two-phase simplex with Fractions.
       A: list of rows (each length n), b: rhs (length m), c: objective (length n)."""
    m = len(A); n = len(A[0])
    # Make b>=0
    A = [row[:] for row in A]; b = b[:]
    for i in range(m):
        if b[i] < 0:
            A[i] = [-x for x in A[i]]; b[i] = -b[i]
    # Phase 1: add artificials
    # Tableau columns: n structural + m artificial.  Basis = artificials.
    N = n + m
    T = [A[i][:] + [F(1) if j == i else F(0) for j in range(m)] + [b[i]] for i in range(m)]
    basis = [n + i for i in range(m)]
    def pivot(T, basis, obj_row):
        m = len(T); width = len(T[0])
        while True:
            # find entering col with negative reduced cost (we minimize obj_row dot)
            col = -1
            for j in range(width - 1):
                if obj_row[j] < 0: col = j; break
            if col == -1: break
            # ratio test
            best = None; row = -1
            for i in range(m):
                if T[i][col] > 0:
                    r = T[i][-1] / T[i][col]
                    if best is None or r < best: best = r; row = i
            if row == -1: return None  # unbounded
            pv = T[row][col]; T[row] = [x / pv for x in T[row]]
            for i in range(m):
                if i != row and T[i][col] != 0:
                    f = T[i][col]; T[i] = [T[i][j] - f * T[row][j] for j in range(width)]
            f = obj_row[col]
            if f != 0:
                obj_row[:] = [obj_row[j] - f * T[row][j] for j in range(width)]
            basis[row] = col
        return True
    # phase-1 objective: minimize sum of artificials
    phase1 = [F(0)] * N + [F(0)]
    for j in range(n, N): phase1[j] = F(1)
    # reduce phase1 over basis
    for i in range(m):
        phase1 = [phase1[j] - T[i][j] for j in range(N + 1)]
    r = pivot(T, basis, phase1)
    if r is None: return None
    # feasibility: artificial sum 0
    if -phase1[-1] != 0:
        return None  # infeasible
    # drop artificial columns
    T2 = [row[:n] + [row[-1]] for row in T]
    # Phase 2: maximize c.x  -> minimize -c
    obj2 = [-ci for ci in c] + [F(0)]
    for i in range(m):
        bcol = basis[i]
        if bcol < n and obj2[bcol] != 0:
            f = obj2[bcol]
            obj2 = [obj2[j] - f * T2[i][j] for j in range(n + 1)]
    r = pivot(T2, basis, obj2)
    if r is None: return None
    return obj2[-1]  # = max c.x (since we tracked -(-min) ... obj2[-1] holds value)

# Build the table for consec, k=8..13, R=1,2,3.
print(f"\n{'k':>3} {'R':>2} {'moment-LP':>12} {'subset-LP':>12} {'measS7':>10} {'cap':>10}  verdict")
for k in [8, 9, 10, 11]:
    E = consec(k); q, _ = q_upper(E); S = Svec_from_q(q); s7 = measS7(E); ck = cap(k)
    for R in [1, 2, 3]:
        mlp = moment_lp_max_p0(S, R)
        slp = subset_lp_max(E, R)
        verdict = "COLLAPSE" if (mlp is not None and slp is not None and mlp == slp) else \
                  ("subset IMPROVES" if (mlp is not None and slp is not None and slp < mlp) else "?")
        mstr = f"{float(mlp):.5f}" if mlp is not None else "infeas"
        sstr = f"{float(slp):.5f}" if slp is not None else "infeas"
        print(f"{k:>3} {R:>2} {mstr:>12} {sstr:>12} {float(s7):>10.5f} {float(ck):>10.5f}  {verdict}")
    print()

print("\nDONE (Part B). See (C) for the configuration-polytope integrality test.")
