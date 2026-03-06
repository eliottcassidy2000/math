#!/usr/bin/env python3
"""
Symbolic proof of OCF at n=8 using SymPy.

27 arc variables, 40320 permutations. Expected to take hours.

delta_I = -2*sum(s_x*H(B_x)) + 2*(D5-C5) + 2*(D7-C7)

Note: No 9-cycle term needed since n=8 < 9.

Instance: kind-pasteur-2026-03-05-S7
"""

from sympy import symbols, expand
from itertools import permutations, combinations
import time

n = 8
I, J = 0, 1
others = [2, 3, 4, 5, 6, 7]

# Variables: 6 p's + 6 q's + C(6,2)=15 internals = 27 total
p_vars = {x: symbols(f'p{x}') for x in others}
q_vars = {x: symbols(f'q{x}') for x in others}
t_vars = {}
for a in others:
    for b in others:
        if a < b:
            t_vars[(a,b)] = symbols(f't{a}{b}')

n_vars = len(p_vars) + len(q_vars) + len(t_vars)
print(f"Variables: {len(p_vars)}+{len(q_vars)}+{len(t_vars)} = {n_vars}")

def T(a, b):
    if a == b: return 0
    if a == I and b == J: return 1
    if a == J and b == I: return 0
    if a in others and b == I: return p_vars[a]
    if a == I and b in others: return 1 - p_vars[b]
    if a == J and b in others: return q_vars[b]
    if a in others and b == J: return 1 - q_vars[a]
    if a in others and b in others:
        if a < b: return t_vars[(a,b)]
        else: return 1 - t_vars[(b,a)]
    raise ValueError(f"Unknown ({a},{b})")

def Tp(a, b):
    if a == I and b == J: return 0
    if a == J and b == I: return 1
    return T(a, b)


print(f"=== n={n} Symbolic Proof ===")
start = time.time()

# H(T): sum over 40320 permutations
print(f"Computing H(T) ({n}! = 40320 perms)...")
HT = 0
nonzero = 0
for i_perm, perm in enumerate(permutations(range(n))):
    w = 1
    for k in range(n - 1):
        arc = T(perm[k], perm[k+1])
        if arc == 0:
            w = 0
            break
        w *= arc
    if w != 0:
        HT += w
        nonzero += 1
    if (i_perm + 1) % 5000 == 0:
        print(f"  {i_perm+1}/40320 perms, {nonzero} nonzero, {time.time()-start:.1f}s")

print(f"Expanding H(T)...")
t1 = time.time()
HT = expand(HT)
print(f"  H(T): {len(HT.as_coefficients_dict())} monomials ({time.time()-t1:.1f}s)")

print(f"Computing H(T') (40320 perms)...")
t1 = time.time()
HTp = 0
nonzero_p = 0
for i_perm, perm in enumerate(permutations(range(n))):
    w = 1
    for k in range(n - 1):
        arc = Tp(perm[k], perm[k+1])
        if arc == 0:
            w = 0
            break
        w *= arc
    if w != 0:
        HTp += w
        nonzero_p += 1
    if (i_perm + 1) % 5000 == 0:
        print(f"  {i_perm+1}/40320 perms, {nonzero_p} nonzero, {time.time()-t1:.1f}s")

print(f"Expanding H(T')...")
HTp = expand(HTp)
print(f"  H(T'): {len(HTp.as_coefficients_dict())} monomials ({time.time()-t1:.1f}s)")

print("Computing delta_H = H(T) - H(T')...")
delta_H = expand(HT - HTp)
n_terms_H = len(delta_H.as_coefficients_dict())
print(f"  delta_H: {n_terms_H} monomials")

# delta_I = -2*sum(s_x*H(B_x)) + 2*(D5-C5) + 2*(D7-C7)
print("Computing delta_I components...")

s = {x: 1 - p_vars[x] - q_vars[x] for x in others}

def h_sub(verts):
    total = 0
    for perm in permutations(verts):
        w = 1
        for k in range(len(perm) - 1):
            w *= T(perm[k], perm[k+1])
        total += w
    return expand(total)

# H(B_x) for 5-vertex sub-tournaments
formula_sum = 0
for x in others:
    Bx = [v for v in others if v != x]
    print(f"  Computing H(B_{x}) ({len(Bx)}! = {len(list(permutations(Bx)))} perms)...", end="", flush=True)
    t1 = time.time()
    hbx = h_sub(Bx)
    print(f" {len(hbx.as_coefficients_dict())} terms ({time.time()-t1:.1f}s)")
    formula_sum += s[x] * hbx
formula_sum = expand(formula_sum)

# 5-cycles: choose 3 from others
print("Computing 5-cycle terms...")
t1 = time.time()
D5 = C5 = 0
for subset in combinations(others, 3):
    for perm in permutations(subset):
        v1, v2, v3 = perm
        D5 += T(J, v1) * T(v1, v2) * T(v2, v3) * T(v3, I)
        C5 += T(I, v1) * T(v1, v2) * T(v2, v3) * T(v3, J)
D5 = expand(D5)
C5 = expand(C5)
print(f"  D5: {len(D5.as_coefficients_dict())} terms, C5: {len(C5.as_coefficients_dict())} terms ({time.time()-t1:.1f}s)")

# 7-cycles: choose 5 from others
print("Computing 7-cycle terms...")
t1 = time.time()
D7 = C7 = 0
for subset in combinations(others, 5):
    for perm in permutations(subset):
        v1, v2, v3, v4, v5 = perm
        D7 += T(J, v1) * T(v1, v2) * T(v2, v3) * T(v3, v4) * T(v4, v5) * T(v5, I)
        C7 += T(I, v1) * T(v1, v2) * T(v2, v3) * T(v3, v4) * T(v4, v5) * T(v5, J)
D7 = expand(D7)
C7 = expand(C7)
print(f"  D7: {len(D7.as_coefficients_dict())} terms, C7: {len(C7.as_coefficients_dict())} terms ({time.time()-t1:.1f}s)")

delta_I = expand(-2 * formula_sum + 2 * (D5 - C5) + 2 * (D7 - C7))
n_terms_I = len(delta_I.as_coefficients_dict())
print(f"  delta_I: {n_terms_I} monomials")

print(f"\nChecking delta_H == delta_I...")
diff = expand(delta_H - delta_I)

elapsed = time.time() - start
if diff == 0:
    print(f"\n*** PROVED: H(T)-H(T') = delta_I as POLYNOMIAL IDENTITY at n={n} ***")
    print(f"    (Computed in {elapsed:.1f}s)")
    print(f"    {n_vars} arc variables, {n_terms_H} monomial terms")
    print(f"    Formula: -2*sum(s_x*H(B_x)) + 2*(D5-C5) + 2*(D7-C7)")
else:
    nd = len(diff.as_coefficients_dict())
    print(f"\nDifference has {nd} terms (computed in {elapsed:.1f}s)")
    if nd <= 10:
        print(f"diff = {diff}")
