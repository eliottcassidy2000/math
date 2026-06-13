#!/usr/bin/env python3
"""
Symbolic proof of OCF at n=6 using SymPy.

H(T)-H(T') = delta_I where
delta_I = -2*sum(s_x*H(B_x)) + 2*(D5-C5)

14 arc variables, 720 permutations. Should be feasible.

Instance: kind-pasteur-2026-03-05-S7
"""

from sympy import symbols, expand
from itertools import permutations, combinations
import time

# Variables: 4 p's, 4 q's, C(4,2)=6 internals = 14 total
p2, p3, p4, p5 = symbols('p2 p3 p4 p5')
q2, q3, q4, q5 = symbols('q2 q3 q4 q5')
t23, t24, t25, t34, t35, t45 = symbols('t23 t24 t25 t34 t35 t45')

n = 6
I, J = 0, 1
others = [2, 3, 4, 5]

# Build arc function
p = {2: p2, 3: p3, 4: p4, 5: p5}
q = {2: q2, 3: q3, 4: q4, 5: q5}
t = {(2,3): t23, (2,4): t24, (2,5): t25, (3,4): t34, (3,5): t35, (4,5): t45}

def T(a, b):
    if a == b: return 0
    if a == I and b == J: return 1
    if a == J and b == I: return 0
    if a in others and b == I: return p[a]
    if a == I and b in others: return 1 - p[b]
    if a == J and b in others: return q[b]
    if a in others and b == J: return 1 - q[a]
    if a in others and b in others:
        if a < b: return t[(a,b)]
        else: return 1 - t[(b,a)]
    raise ValueError(f"Unknown arc ({a},{b})")

def Tp(a, b):
    if a == I and b == J: return 0
    if a == J and b == I: return 1
    return T(a, b)


print("=== n=6 Symbolic Proof ===")
start = time.time()

# Compute H(T) and H(T') by summing over all 720 permutations
print("Computing H(T) (720 perms)...")
HT = 0
count = 0
for perm in permutations(range(n)):
    w = 1
    for k in range(n - 1):
        arc = T(perm[k], perm[k+1])
        if arc == 0:
            w = 0
            break
        w *= arc
    if w != 0:
        HT += w
        count += 1

HT = expand(HT)
print(f"  {count} nonzero terms, H(T) has {len(HT.as_coefficients_dict())} monomials ({time.time()-start:.1f}s)")

print("Computing H(T') (720 perms)...")
t1 = time.time()
HTp = 0
for perm in permutations(range(n)):
    w = 1
    for k in range(n - 1):
        arc = Tp(perm[k], perm[k+1])
        if arc == 0:
            w = 0
            break
        w *= arc
    if w != 0:
        HTp += w

HTp = expand(HTp)
print(f"  H(T') has {len(HTp.as_coefficients_dict())} monomials ({time.time()-t1:.1f}s)")

print("Computing delta_H...")
delta_H = expand(HT - HTp)
print(f"  delta_H has {len(delta_H.as_coefficients_dict())} monomials")

# Compute delta_I = -2*sum(s_x*H(B_x)) + 2*(D5-C5)
print("Computing delta_I...")

s = {x: 1 - p[x] - q[x] for x in others}

# H(B_x) for 3-vertex sub-tournaments B_x = others \ {x}
def h_sub(verts):
    total = 0
    for perm in permutations(verts):
        w = 1
        for k in range(len(perm) - 1):
            w *= T(perm[k], perm[k+1])
        total += w
    return expand(total)

formula_sum = 0
for x in others:
    Bx = [v for v in others if v != x]
    hbx = h_sub(Bx)
    formula_sum += s[x] * hbx
formula_sum = expand(formula_sum)

# 5-cycle counts
D5 = 0  # lost (destroyed)
C5 = 0  # gained (created)
for subset in combinations(others, 3):
    for perm in permutations(subset):
        v1, v2, v3 = perm
        D5 += T(J, v1) * T(v1, v2) * T(v2, v3) * T(v3, I)
        C5 += T(I, v1) * T(v1, v2) * T(v2, v3) * T(v3, J)

D5 = expand(D5)
C5 = expand(C5)

delta_I = expand(-2 * formula_sum + 2 * (D5 - C5))
print(f"  delta_I has {len(delta_I.as_coefficients_dict())} monomials")

# Check
print("\nChecking delta_H == delta_I...")
diff = expand(delta_H - delta_I)
print(f"delta_H - delta_I = {diff}")

elapsed = time.time() - start
if diff == 0:
    print(f"\n*** PROVED: H(T)-H(T') = delta_I as POLYNOMIAL IDENTITY at n=6 ***")
    print(f"    (Computed in {elapsed:.1f}s)")
    print(f"    14 arc variables, {len(delta_H.as_coefficients_dict())} monomial terms")
    print(f"    Formula: -2*sum(s_x*H(B_x)) + 2*(D5-C5)")
else:
    print(f"\nFAILED in {elapsed:.1f}s")
    # Check number of differing terms
    d = diff.as_coefficients_dict()
    print(f"Difference has {len(d)} terms")
    for mono, coeff in list(d.items())[:5]:
        print(f"  {coeff} * {mono}")
