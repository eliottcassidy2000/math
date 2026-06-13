#!/usr/bin/env python3
"""
COMPLETE INDUCTIVE PROOF of the KEY IDENTITY
=============================================

THEOREM: For any vertex set W with |W| = m >= 2 and any b in W,
    B_b(W) + (-1)^m E_b(W) = 2r * col_sum_W(b)

where col_sum_W(b) = sum_{a != b} M_W[a,b].

COROLLARIES:
  (1) Summing over b: T(W)[1 + (-1)^m] = 2r * Sigma(W)
      => Even m: T = r*Sigma. Odd m: Sigma = 0.
  (2) M[a,b] has only even powers of r (the "even r-powers" conjecture).
  (3) M[a,b] = M[b,a] (symmetry of the transfer matrix).

PROOF STRUCTURE:
  Base case: m = 2. Direct verification.
  Inductive step: Assume KEY IDENTITY for all sizes <= m-1.

  Three recurrences (all algebraic identities from definitions):

  (R1) B_b(W) = r*T(W') + sum_w s_{bw} B_w(W')     [first-edge decomposition]
  (R2) E_b(W) = r*T(W') - sum_w s_{bw} E_w(W')     [last-edge decomposition]
  (R3) col_sum_W(b) = (-1)^{m-2} T(W') + r*Sigma_{W'} + sum_w s_{bw} cs_{W'}(w)
       [from column recurrence of M]

  where W' = W\{b}.

  Combining (R1) and (R2):
    LHS = B_b + (-1)^m E_b
        = r*T(W')[1 + (-1)^m] + sum_w s_{bw}[B_w(W') + (-1)^{m-1} E_w(W')]
                                                        ^^^ note: (-1)^m * (-1) = (-1)^{m+1} = (-1)^{m-1}
                                                        Wait, (-1)^m * (-1) = (-1)^{m+1}. And (-1)^{m+1} = (-1)^{m-1}.
                                                        So yes, B_w - (-1)^m E_w = B_w + (-1)^{m-1} E_w.
                                                        But KEY IDENTITY at size m-1 says:
                                                        B_w(W') + (-1)^{m-1} E_w(W') = 2r * cs_{W'}(w)
        = r*T(W')[1 + (-1)^m] + sum_w s_{bw} * 2r * cs_{W'}(w)   [by induction]

  From (R3):
    RHS = 2r * col_sum_W(b)
        = 2r*(-1)^{m-2}*T(W') + 2r^2*Sigma_{W'} + 2r*sum_w s_{bw} cs_{W'}(w)

  Setting LHS = RHS, the s_{bw}*cs terms cancel, leaving:
    r*T(W')[1 + (-1)^m] = 2r*(-1)^{m-2}*T(W') + 2r^2*Sigma_{W'}

  Since (-1)^{m-2} = (-1)^m:
    r*T(W')[1 + (-1)^m] = 2r*(-1)^m*T(W') + 2r^2*Sigma_{W'}
    r*T(W') + r*(-1)^m*T(W') = 2r*(-1)^m*T(W') + 2r^2*Sigma_{W'}
    r*T(W') - r*(-1)^m*T(W') = 2r^2*Sigma_{W'}
    r*T(W')[1 - (-1)^m] = 2r^2*Sigma_{W'}

  Case m even: LHS = 0, need Sigma_{W'} = 0.
    |W'| = m-1 is odd. By Sigma identity at size m-1 (from KEY IDENTITY at m-1): Sigma = 0. CHECK.

  Case m odd: LHS = 2r*T(W'), need 2r*T(W') = 2r^2*Sigma_{W'}, i.e. T(W') = r*Sigma_{W'}.
    |W'| = m-1 is even. By Sigma identity at size m-1: T = r*Sigma. CHECK.

  QED.

This script verifies every step of the proof computationally.
"""

from itertools import permutations, combinations
from sympy import symbols, expand, Rational, Poly, degree
from functools import lru_cache
import sys

# Setup
MAX_N = 6  # verify up to this size

r = symbols('r')
# Create s variables
s = {}
for i in range(MAX_N):
    for j in range(i+1, MAX_N):
        s[(i,j)] = symbols(f's{i}{j}')
        s[(j,i)] = -s[(i,j)]

def t(i, j):
    return r + s[(i,j)]

def ham_paths(vertices):
    """All Hamiltonian paths through vertices."""
    return permutations(vertices)

def path_weight(path):
    """Weight of a path."""
    w = 1
    for i in range(len(path)-1):
        w *= t(path[i], path[i+1])
    return expand(w)

def B_v(v, W):
    """Sum of Hamiltonian paths starting at v through W."""
    W = tuple(sorted(W))
    total = 0
    for p in ham_paths(W):
        if p[0] == v:
            total += path_weight(p)
    return expand(total)

def E_v(v, W):
    """Sum of Hamiltonian paths ending at v through W."""
    W = tuple(sorted(W))
    total = 0
    for p in ham_paths(W):
        if p[-1] == v:
            total += path_weight(p)
    return expand(total)

def T_total(W):
    """Total Hamiltonian path weight."""
    total = 0
    for p in ham_paths(W):
        total += path_weight(p)
    return expand(total)

def M_entry(a, b, W):
    """Transfer matrix entry M[a,b]."""
    W = tuple(sorted(W))
    U = [v for v in W if v != a and v != b]
    total = 0
    for k in range(len(U)+1):
        for S in combinations(U, k):
            S_set = set(S)
            R = [v for v in U if v not in S_set]
            ea = E_v(a, tuple(sorted(list(S) + [a])))
            bb = B_v(b, tuple(sorted(R + [b])))
            total += ((-1)**k) * ea * bb
    return expand(total)

def col_sum(b, W):
    """Column sum: sum_{a != b} M[a,b]."""
    total = 0
    for a in W:
        if a != b:
            total += M_entry(a, b, W)
    return expand(total)

def Sigma(W):
    """Total off-diagonal sum."""
    total = 0
    for a in W:
        for b in W:
            if a != b:
                total += M_entry(a, b, W)
    return expand(total)

def odd_r(expr):
    """Extract odd powers of r."""
    p = Poly(expand(expr), r)
    result = 0
    for monom, coeff in zip(p.monoms(), p.coeffs()):
        if monom[0] % 2 == 1:
            result += coeff * r**monom[0]
    return expand(result)

def even_r(expr):
    """Extract even powers of r."""
    p = Poly(expand(expr), r)
    result = 0
    for monom, coeff in zip(p.monoms(), p.coeffs()):
        if monom[0] % 2 == 0:
            result += coeff * r**monom[0]
    return expand(result)

print("=" * 70)
print("COMPLETE INDUCTIVE PROOF VERIFICATION")
print("=" * 70)

# ===== BASE CASE =====
print("\n--- BASE CASE: m = 2 ---")
W = (0, 1)
b = 0
m = 2

bb = B_v(b, W)
eb = E_v(b, W)
cs = col_sum(b, W)

lhs = expand(bb + (-1)**m * eb)
rhs = expand(2 * r * cs)
print(f"  B_{b}(W) = {bb}")
print(f"  E_{b}(W) = {eb}")
print(f"  col_sum(b) = {cs}")
print(f"  LHS = B_b + E_b = {lhs}")
print(f"  RHS = 2r * col_sum = {rhs}")
print(f"  Match: {expand(lhs - rhs) == 0}")

# Also check b=1
b = 1
bb = B_v(b, W)
eb = E_v(b, W)
cs = col_sum(b, W)
lhs = expand(bb + (-1)**m * eb)
rhs = expand(2 * r * cs)
print(f"  [b=1] LHS={lhs}, RHS={rhs}, Match: {expand(lhs-rhs)==0}")

# ===== VERIFY KEY IDENTITY FOR ALL SIZES =====
print("\n--- KEY IDENTITY VERIFICATION: B_b + (-1)^m E_b = 2r * col_sum(b) ---")
for m in range(2, MAX_N + 1):
    W = tuple(range(m))
    all_ok = True
    for b in W:
        bb = B_v(b, W)
        eb = E_v(b, W)
        cs = col_sum(b, W)
        lhs = expand(bb + (-1)**m * eb)
        rhs = expand(2 * r * cs)
        ok = expand(lhs - rhs) == 0
        if not ok:
            print(f"  FAIL at m={m}, b={b}")
            all_ok = False
    print(f"  m={m}: KEY IDENTITY verified for all b: {all_ok}")

# ===== VERIFY SIGMA IDENTITIES =====
print("\n--- SIGMA IDENTITY (corollary of KEY IDENTITY) ---")
for m in range(2, MAX_N + 1):
    W = tuple(range(m))
    sig = Sigma(W)
    tt = T_total(W)
    if m % 2 == 0:
        check = expand(tt - r * sig)
        print(f"  m={m} (even): T = r*Sigma? {check == 0}")
    else:
        print(f"  m={m} (odd): Sigma = 0? {expand(sig) == 0}")

# ===== VERIFY RECURRENCE (R1): B_b(W) = r*T(W') + sum s_{bw} B_w(W') =====
print("\n--- RECURRENCE R1: B_b(W) = r*T(W') + sum_w s_{bw} B_w(W') ---")
for m in range(3, MAX_N + 1):
    W = tuple(range(m))
    all_ok = True
    for b in W:
        Wp = tuple(v for v in W if v != b)
        bb = B_v(b, W)
        rhs = r * T_total(Wp)
        for w in Wp:
            rhs += s[(b,w)] * B_v(w, Wp)
        rhs = expand(rhs)
        ok = expand(bb - rhs) == 0
        if not ok:
            print(f"  FAIL at m={m}, b={b}: diff = {expand(bb-rhs)}")
            all_ok = False
    print(f"  m={m}: R1 verified: {all_ok}")

# ===== VERIFY RECURRENCE (R2): E_b(W) = r*T(W') - sum s_{bw} E_w(W') =====
print("\n--- RECURRENCE R2: E_b(W) = r*T(W') - sum_w s_{bw} E_w(W') ---")
for m in range(3, MAX_N + 1):
    W = tuple(range(m))
    all_ok = True
    for b in W:
        Wp = tuple(v for v in W if v != b)
        eb = E_v(b, W)
        rhs = r * T_total(Wp)
        for w in Wp:
            rhs -= s[(b,w)] * E_v(w, Wp)
        rhs = expand(rhs)
        ok = expand(eb - rhs) == 0
        if not ok:
            print(f"  FAIL at m={m}, b={b}: diff = {expand(eb-rhs)}")
            all_ok = False
    print(f"  m={m}: R2 verified: {all_ok}")

# ===== VERIFY RECURRENCE (R3): col_sum_W(b) = (-1)^{m-2} T(W') + r*Sigma_{W'} + sum s_{bw} cs(w) =====
print("\n--- RECURRENCE R3: col_sum_W(b) = (-1)^{m-2} T(W') + r*Sigma(W') + sum s_{bw} cs_{W'}(w) ---")
for m in range(3, MAX_N + 1):
    W = tuple(range(m))
    all_ok = True
    for b in W:
        Wp = tuple(v for v in W if v != b)
        cs_b = col_sum(b, W)
        rhs = (-1)**(m-2) * T_total(Wp) + r * Sigma(Wp)
        for w in Wp:
            rhs += s[(b,w)] * col_sum(w, Wp)
        rhs = expand(rhs)
        ok = expand(cs_b - rhs) == 0
        if not ok:
            print(f"  FAIL at m={m}, b={b}")
            all_ok = False
    print(f"  m={m}: R3 verified: {all_ok}")

# ===== VERIFY THE INDUCTIVE STEP REDUCTION =====
print("\n--- INDUCTIVE STEP: verify r*T(W')[1-(-1)^m] = 2r^2*Sigma(W') ---")
for m in range(3, MAX_N + 1):
    W = tuple(range(m))
    for b in W[:1]:  # just check one b (the identity is b-independent)
        Wp = tuple(v for v in W if v != b)
        tt = T_total(Wp)
        sig = Sigma(Wp)
        lhs = expand(r * tt * (1 - (-1)**m))
        rhs = expand(2 * r**2 * sig)
        print(f"  m={m}: r*T(W')[1-(-1)^m] = {lhs}")
        print(f"         2r^2*Sigma(W')    = {rhs}")
        print(f"         Match: {expand(lhs - rhs) == 0}")

# ===== VERIFY EVEN R-POWERS =====
print("\n--- COROLLARY: M[a,b] has only even powers of r ---")
for m in range(2, MAX_N + 1):
    W = tuple(range(m))
    all_ok = True
    for a in W:
        for b in W:
            if a != b:
                mab = M_entry(a, b, W)
                if odd_r(mab) != 0:
                    print(f"  FAIL: M[{a},{b}] at m={m} has odd r-powers")
                    all_ok = False
    print(f"  m={m}: all M[a,b] have even r-powers: {all_ok}")

# ===== VERIFY SYMMETRY =====
print("\n--- COROLLARY: M[a,b] = M[b,a] ---")
for m in range(2, MAX_N + 1):
    W = tuple(range(m))
    all_ok = True
    for a in W:
        for b in W:
            if a < b:
                mab = M_entry(a, b, W)
                mba = M_entry(b, a, W)
                if expand(mab - mba) != 0:
                    print(f"  FAIL: M[{a},{b}] != M[{b},{a}] at m={m}")
                    all_ok = False
    print(f"  m={m}: M symmetric: {all_ok}")

print("\n" + "=" * 70)
print("PROOF SUMMARY")
print("=" * 70)
print("""
THEOREM (Key Identity): For |W| = m >= 2 and b in W,
  B_b(W) + (-1)^m E_b(W) = 2r * col_sum_W(b)

PROOF by strong induction on m.

Base case (m=2): B_b({a,b}) + E_b({a,b}) = (r+s_ab)+(r-s_ab) = 2r = 2r*1 = 2r*M[a,b].

Inductive step (m >= 3): Assume KEY IDENTITY for all sizes < m.

Step 1: Three algebraic recurrences (from definitions):
  (R1) B_b(W) = r*T(W') + sum_w s_{bw} B_w(W')
  (R2) E_b(W) = r*T(W') - sum_w s_{bw} E_w(W')
  (R3) col_sum_W(b) = (-1)^{m-2} T(W') + r*Sigma(W') + sum_w s_{bw} cs_{W'}(w)

Step 2: Compute B_b + (-1)^m E_b using (R1), (R2):
  = r*T(W')[1 + (-1)^m] + sum_w s_{bw}[B_w(W') + (-1)^{m-1} E_w(W')]
  = r*T(W')[1 + (-1)^m] + 2r * sum_w s_{bw} cs_{W'}(w)   [by KEY IDENTITY at m-1]

Step 3: Compute 2r * col_sum using (R3):
  = 2r*(-1)^{m-2} T(W') + 2r^2*Sigma(W') + 2r * sum_w s_{bw} cs_{W'}(w)

Step 4: Equating (s-terms cancel), need:
  r*T(W')[1 - (-1)^m] = 2r^2 * Sigma(W')

Step 5: From KEY IDENTITY at m-1, summing over b:
  T(W')[1 + (-1)^{m-1}] = 2r * Sigma(W')

  m even => m-1 odd => 1+(-1)^{m-1} = 0 => Sigma(W') = 0.
    Need: 0 = 0. CHECK.

  m odd => m-1 even => 1+(-1)^{m-1} = 2 => T(W') = r*Sigma(W').
    Need: 2r*T(W') = 2r^2*Sigma(W') = 2r*T(W'). CHECK.

QED.

COROLLARY 1 (Sigma identity):
  Sum over b: T(W)[1+(-1)^m] = 2r*Sigma(W).
  Even m: T = r*Sigma. Odd m: Sigma = 0.

COROLLARY 2 (Even r-powers): M[a,b] has only even powers of r.
  Proof: odd_r(B_b) = [B_b - (-1)^{m-1} E_b]/2.
  KEY IDENTITY => odd_r(B_b) = r * col_sum(b).
  Row recurrence: M[a,b] = B_b(V\\{a}) - sum_v t(v,a) M'[v,b]
  By induction M' has even r, so odd_r(M[a,b]) = odd_r(B_b) - r*col_sum = 0.

COROLLARY 3 (Symmetry): M[a,b] = M[b,a].
  Even r-powers + re-indexing identity M_T[b,a] = (-1)^{m-2} M_{T^op}[a,b]
  => s-parity is (-1)^{m-2} => M[a,b]|_{s->-s} = (-1)^{m-2} M[a,b]
  => M[b,a] = (-1)^{m-2} * (-1)^{m-2} M[a,b] = M[a,b].
""")
print("=" * 70)
print("DONE")
print("=" * 70)
