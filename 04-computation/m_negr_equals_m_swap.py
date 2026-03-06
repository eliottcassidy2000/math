#!/usr/bin/env python3
"""
ALGEBRAIC PROOF: M[a,b](-r) = M[b,a](r)

This is the bridge between the two forms of the symmetry conjecture:
  (i) M[a,b] = M[b,a]   (transfer matrix symmetry)
  (ii) M[a,b] has only even powers of r   (even-r-powers)

PROOF:
  1. T(-r)[i,j] = -r + s_ij = -(r - s_ij) = -T(r)[j,i] = -T(r)^T[i,j]
     So T(-r) = -T(r)^T.

  2. For a path (v_1, ..., v_k, a) with k arcs, weight under -T^T is:
     prod (-T^T)[v_i, v_{i+1}] = prod (-T[v_{i+1}, v_i])
     = (-1)^k * prod T[v_{i+1}, v_i]
     = (-1)^k * [weight of REVERSED path (a, v_k, ..., v_1) under T]

  3. So: E_a(S+a; -T^T) = (-1)^|S| * B_a(S+a; T)
         B_b(R+b; -T^T) = (-1)^|R| * E_b(R+b; T)

  4. M[a,b](-r) = sum_S (-1)^|S| * E_a(S+a; T(-r)) * B_b(R+b; T(-r))
                = sum_S (-1)^|S| * [(-1)^|S| B_a(S+a; T)] * [(-1)^|R| E_b(R+b; T)]
                = sum_S (-1)^{|S| + |S| + |R|} * B_a(S+a; T) * E_b(R+b; T)
                = sum_S (-1)^{|R|} * B_a(S+a; T) * E_b(R+b; T)   [since 2|S| is even]

  5. Substituting S' = R (relabel the subset):
     = sum_{S'} (-1)^{|S'|} * B_a(R'+a; T) * E_b(S'+b; T)   [where R' = U\S']

  6. But M[b,a](r) = sum_S (-1)^|S| * E_b(S+b; T) * B_a(R+a; T)

  7. Comparing: M[a,b](-r) = M[b,a](r). QED.

CONSEQUENCE: M[a,b] = M[b,a] iff M[a,b](-r) = M[a,b](r) iff only even r-powers.

This script verifies the identity computationally.

kind-pasteur-2026-03-06-S23
"""

from itertools import permutations
from sympy import symbols, expand, Poly

def setup(n):
    r = symbols('r')
    sv = {}
    for i in range(n):
        for j in range(i+1, n):
            sv[(i,j)] = symbols(f's{i}{j}')
    def s(i, j):
        if i == j: return 0
        if i < j: return sv[(i,j)]
        return -sv[(j,i)]
    def t(i, j):
        if i == j: return 0
        return r + s(i, j)
    def t_neg(i, j):
        """T(-r)[i,j] = -r + s_ij"""
        if i == j: return 0
        return -r + s(i, j)
    return r, sv, s, t, t_neg

def hp(t_fn, vset, start=None, end=None):
    vl = sorted(vset)
    if len(vl) == 0: return 0
    if len(vl) == 1:
        if start is not None and vl[0] != start: return 0
        if end is not None and vl[0] != end: return 0
        return 1
    total = 0
    for p in permutations(vl):
        if start is not None and p[0] != start: continue
        if end is not None and p[-1] != end: continue
        prod = 1
        for i in range(len(p)-1):
            prod *= t_fn(p[i], p[i+1])
        total += prod
    return expand(total)

def transfer_M(t_fn, n, a, b):
    U = [v for v in range(n) if v != a and v != b]
    result = 0
    for mask in range(1 << len(U)):
        S = [U[i] for i in range(len(U)) if mask & (1 << i)]
        R = [U[i] for i in range(len(U)) if not (mask & (1 << i))]
        sign = (-1)**len(S)
        ea = hp(t_fn, set(S)|{a}, end=a)
        bb = hp(t_fn, set(R)|{b}, start=b)
        result += sign * ea * bb
    return expand(result)

print("=" * 70)
print("VERIFICATION: M[a,b](-r) = M[b,a](r)")
print("=" * 70)

for n in [4, 5, 6]:
    r, sv, s, t, t_neg = setup(n)
    a, b = 0, 1

    # Compute M[a,b](r)
    M_ab = transfer_M(t, n, a, b)

    # Compute M[b,a](r)
    M_ba = transfer_M(t, n, b, a)

    # Compute M[a,b](-r) by substitution
    M_ab_neg = expand(M_ab.subs(r, -r))

    # Check M[a,b](-r) = M[b,a](r)
    diff1 = expand(M_ab_neg - M_ba)
    print(f"\n  n={n}: M[a,b](-r) - M[b,a](r) = {diff1}")
    assert diff1 == 0, f"FAILED at n={n}!"

    # Check even-r-powers iff symmetry
    diff2 = expand(M_ab - M_ba)
    p_ab = Poly(M_ab, r)
    has_odd = False
    for power in range(p_ab.degree() + 1):
        coeff = expand(p_ab.nth(power))
        if power % 2 == 1 and coeff != 0:
            has_odd = True
    if diff2 == 0:
        print(f"         M[a,b] = M[b,a] (symmetric)")
        print(f"         Even r-powers only: {not has_odd}")
    else:
        print(f"         M[a,b] != M[b,a] (difference has {len(diff2.as_ordered_terms())} terms)")
        print(f"         Even r-powers only: {not has_odd}")

print("\n" + "=" * 70)
print("ALGEBRAIC PROOF SUMMARY")
print("=" * 70)
print("""
PROVED: M[a,b](-r, s) = M[b,a](r, s) for all c-tournaments.

Proof uses only:
  (1) T(-r) = -T(r)^T  (from t_ij + t_ji = c = 2r)
  (2) Path reversal: reversing a path of k arcs under -T^T gives
      (-1)^k times the reversed path weight under T.
  (3) The inclusion-exclusion sum structure of M.

COROLLARY: The following are equivalent:
  (i)   M[a,b] = M[b,a]  (transfer matrix symmetry)
  (ii)  M[a,b](r,s) has only even powers of r
  (iii) M[a,b](r,s) = M[a,b](r,-s) * (-1)^{n-2}  (s-parity)

The equivalence (i) <=> (ii) follows from M[a,b](-r) = M[b,a](r):
  M[a,b] = M[b,a]  iff  M[a,b](-r) = M[a,b](r)  iff  only even r-powers.

WHAT REMAINS: Proving any ONE of (i), (ii), or (iii).
""")
