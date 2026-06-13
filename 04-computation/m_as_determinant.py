#!/usr/bin/env python3
"""
Can M[a,b] be expressed as a determinant of some matrix?

If M[a,b] = det(some (n-2)x(n-2) matrix depending on a,b,T),
and this matrix is symmetric under a<->b swap, then M[a,b] = M[b,a] follows.

We test: for small n, find a matrix Q such that M[a,b] = det(Q).
At n=4: M[a,b] is degree 2 in arc weights (2 arcs in each cover).
  A 2x2 determinant ad-bc would have degree 2. Can we find Q?

At n=5: M[a,b] is degree 3. A 3x3 determinant has degree 3. Can we find Q?

kind-pasteur-2026-03-06-S23
"""
from itertools import permutations
from sympy import symbols, expand, Poly, Matrix, det, solve
from collections import defaultdict

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
    return r, sv, s, t

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
print("CAN M[a,b] BE A DETERMINANT?")
print("=" * 70)

# n=4, a=0, b=1
n = 4
r, sv, s, t = setup(n)
a, b = 0, 1

M_ab = transfer_M(t, n, a, b)
print(f"\nn=4: M[0,1] = {M_ab}")

# M[a,b] at n=4 has degree 2 in arc weights. Let's check.
p_r = Poly(M_ab, r)
print(f"  r-degree: {p_r.degree()}")
for k in range(p_r.degree() + 1):
    coeff = expand(p_r.nth(k))
    print(f"  r^{k}: {coeff}")

# Check: is M[0,1] = det of a 2x2 matrix?
# M = ad - bc where a,b,c,d are linear in arc weights
# r^2 coefficient: 1 (constant). r^0 coefficient: s_02*s_13 + s_03*s_12 + s_02*s_31 + ...
# wait, let me just check if it's det of t-matrix entries

# Actually, let me check some natural candidates
U = [2, 3]
# Candidate 1: M[a,b] = det [[t(a,c), t(a,d)], [t(b,c), t(b,d)]] ?
# = t(0,2)*t(1,3) - t(0,3)*t(1,2)
cand1 = expand(t(0,2)*t(1,3) - t(0,3)*t(1,2))
print(f"\n  Candidate 1: t(0,2)*t(1,3) - t(0,3)*t(1,2) = {cand1}")
print(f"  Match: {expand(M_ab - cand1) == 0}")

# Candidate 2: det [[t(2,0), t(3,0)], [t(2,1), t(3,1)]]
# = t(2,0)*t(3,1) - t(3,0)*t(2,1)
cand2 = expand(t(2,0)*t(3,1) - t(3,0)*t(2,1))
print(f"\n  Candidate 2: t(2,0)*t(3,1) - t(3,0)*t(2,1) = {cand2}")
print(f"  Match: {expand(M_ab - cand2) == 0}")

# Let me compute M[a,b] more explicitly
M_ab_explicit = transfer_M(t, n, a, b)
print(f"\n  M[0,1] expanded: {expand(M_ab_explicit)}")
M_ba_explicit = transfer_M(t, n, b, a)
print(f"  M[1,0] expanded: {expand(M_ba_explicit)}")
print(f"  M[0,1] = M[1,0]: {expand(M_ab_explicit - M_ba_explicit) == 0}")

# Let me try the FULL transfer matrix for n=4
print("\n" + "=" * 70)
print("FULL TRANSFER MATRIX at n=4")
print("=" * 70)

M_full = {}
for aa in range(n):
    for bb in range(n):
        if aa == bb: continue
        M_full[(aa,bb)] = transfer_M(t, n, aa, bb)

# Check symmetry for all pairs
for aa in range(n):
    for bb in range(aa+1, n):
        diff = expand(M_full[(aa,bb)] - M_full[(bb,aa)])
        print(f"  M[{aa},{bb}] - M[{bb},{aa}] = {'0 OK' if diff == 0 else 'NONZERO!'}")

# Check: is the full matrix M a function of the SKEW-SYMMETRIC part?
# At r=0: M[a,b](0) depends only on s_ij
print("\n  r^0 of M[0,1]:", expand(Poly(M_full[(0,1)], r).nth(0)))
print("  r^2 of M[0,1]:", expand(Poly(M_full[(0,1)], r).nth(2)))

# What is the r^2 coefficient?
r2_coeff = expand(Poly(M_full[(0,1)], r).nth(2))
print(f"\n  r^2 coefficient of M[0,1]: {r2_coeff}")
print(f"  This should be a constant (since total degree is 2, r^2 * s^0): {r2_coeff}")

# Now let's check: is M[a,b] = per(something) - det(something)?
# Or: is M[a,b] related to a PFAFFIAN?
print("\n" + "=" * 70)
print("STRUCTURE OF M[a,b] at n=5")
print("=" * 70)

n = 5
r, sv, s, t = setup(n)
a, b = 0, 1

M_ab = transfer_M(t, n, a, b)
p_r = Poly(M_ab, r)
print(f"  r-degree of M[0,1]: {p_r.degree()}")

for k in range(p_r.degree() + 1):
    coeff = expand(p_r.nth(k))
    if coeff == 0:
        print(f"  r^{k}: 0")
    else:
        terms = coeff.as_ordered_terms()
        print(f"  r^{k}: {len(terms)} terms")

# The r^0 coefficient has s-degree 3
# Can it be a 3x3 determinant?
r0 = expand(p_r.nth(0))
print(f"\n  r^0 = {r0}")

# Check: is r0 = det of some 3x3 matrix of s-values?
# A 3x3 determinant of linear s-entries would have s-degree 3
# Let's see if r0 matches a specific pattern
# Internal vertices at n=5, a=0, b=1: U = {2,3,4}

# Try: det of [[s(2,3), s(2,4), s(2,0)], [s(3,2), s(3,4), s(3,0)], ...]
# This is too speculative. Let me instead check the STRUCTURE of the r^0 term.

# Group terms by which s-variables appear
s_vars = list(sv.values())
ps = Poly(r0, *s_vars)
monomial_data = ps.as_dict()
print(f"\n  r^0 has {len(monomial_data)} monomials:")
for monom, c in sorted(monomial_data.items(), key=lambda x: x[0]):
    s_names = []
    for i, (pair, var) in enumerate(sorted(sv.items())):
        if monom[i] > 0:
            s_names.append(f"s{pair[0]}{pair[1]}^{monom[i]}" if monom[i] > 1 else f"s{pair[0]}{pair[1]}")
    print(f"    {c:+d} * {'*'.join(s_names) if s_names else '1'}")
