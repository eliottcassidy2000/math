#!/usr/bin/env python3
"""
Symbolic proof of OCF at n=5 using SymPy.

Compute H(T)-H(T') as a polynomial in the 9 arc variables,
and verify it equals delta_I = -2*sum(s_x) + 2*(D5-C5).

This proves the identity as a POLYNOMIAL identity, not just on {0,1}.

Instance: kind-pasteur-2026-03-05-S7
"""

from sympy import symbols, expand, simplify, collect, Poly
from itertools import permutations

# Variables
# p_a = T[a][0] (a beats 0=i), q_a = T[1][a] (1=j beats a)
# t_ab = T[a][b] for a<b in {2,3,4}
p2, p3, p4 = symbols('p2 p3 p4')
q2, q3, q4 = symbols('q2 q3 q4')
t23, t24, t34 = symbols('t23 t24 t34')

# Tournament matrix T as function
# Vertices: 0=i, 1=j, 2=a, 3=b, 4=c
def T(a, b):
    if a == b: return 0
    if a == 0 and b == 1: return 1  # i->j
    if a == 1 and b == 0: return 0
    # Interface
    if a == 2 and b == 0: return p2
    if a == 0 and b == 2: return 1 - p2
    if a == 3 and b == 0: return p3
    if a == 0 and b == 3: return 1 - p3
    if a == 4 and b == 0: return p4
    if a == 0 and b == 4: return 1 - p4
    if a == 1 and b == 2: return q2
    if a == 2 and b == 1: return 1 - q2
    if a == 1 and b == 3: return q3
    if a == 3 and b == 1: return 1 - q3
    if a == 1 and b == 4: return q4
    if a == 4 and b == 1: return 1 - q4
    # Internal
    if a == 2 and b == 3: return t23
    if a == 3 and b == 2: return 1 - t23
    if a == 2 and b == 4: return t24
    if a == 4 and b == 2: return 1 - t24
    if a == 3 and b == 4: return t34
    if a == 4 and b == 3: return 1 - t34
    raise ValueError(f"Unknown arc ({a},{b})")


def Tp(a, b):
    """T' = T with i->j flipped to j->i."""
    if a == 0 and b == 1: return 0
    if a == 1 and b == 0: return 1
    return T(a, b)


n = 5
all_verts = list(range(n))


def ham_path_count(Tf):
    """Compute H(T) as symbolic polynomial."""
    total = 0
    for perm in permutations(all_verts):
        weight = 1
        for k in range(n - 1):
            weight *= Tf(perm[k], perm[k + 1])
        total += weight
    return expand(total)


print("Computing H(T)...")
HT = ham_path_count(T)
print(f"H(T) has {len(HT.as_coefficients_dict())} terms")

print("Computing H(T')...")
HTp = ham_path_count(Tp)
print(f"H(T') has {len(HTp.as_coefficients_dict())} terms")

print("Computing delta_H = H(T) - H(T')...")
delta_H = expand(HT - HTp)
print(f"delta_H has {len(delta_H.as_coefficients_dict())} terms")

# Now compute delta_I = -2*sum(s_x) + 2*(D5-C5)
# s_x = 1 - p_x - q_x
s2 = 1 - p2 - q2
s3 = 1 - p3 - q3
s4 = 1 - p4 - q4

# 3-cycle term: D3-C3 = -sum(s_x)
term_3 = -(s2 + s3 + s4)

# 5-cycle term: D5-C5
# Lost 5-cycle (I,J,sigma(1),sigma(2),sigma(3),I):
# T[J][sig1] * T[sig1][sig2] * T[sig2][sig3] * T[sig3][I]
others = [2, 3, 4]
D5 = 0  # lost (destroyed)
C5 = 0  # gained (created)
for perm in permutations(others):
    v1, v2, v3 = perm
    # Lost: (0,1,v1,v2,v3) cycle back to 0
    D5 += T(1, v1) * T(v1, v2) * T(v2, v3) * T(v3, 0)
    # Gained: (1,0,v1,v2,v3) cycle back to 1
    C5 += Tp(0, v1) * T(v1, v2) * T(v2, v3) * T(v3, 1)

D5 = expand(D5)
C5 = expand(C5)
# Note: Tp(0,v1) = T(0,v1) since 0->v1 doesn't involve 0->1 arc
# Actually Tp(0,v1) for v1 in {2,3,4}: T'[0][v1] = T[0][v1] = 1-p_{v1}. So C5 uses the same arcs as T.

term_5 = D5 - C5

delta_I = expand(2 * (term_3 + term_5))
print(f"\ndelta_I has {len(delta_I.as_coefficients_dict())} terms")

# Check
diff = expand(delta_H - delta_I)
print(f"\ndelta_H - delta_I = {diff}")

if diff == 0:
    print("\n*** PROVED: H(T)-H(T') = delta_I as POLYNOMIAL IDENTITY at n=5 ***")
    print("This holds for ALL values of the 9 arc variables, not just {0,1}!")
    print("\nThe polynomial identity is:")
    print(f"  H(T) - H(T') = 2*(-sum(s_x)) + 2*(D5-C5)")
    print(f"  = -2*(s2+s3+s4) + 2*(D5-C5)")

    # Express delta_H in a clean form
    print("\n--- Simplified delta_H ---")
    delta_H_collected = collect(delta_H, [p2, p3, p4, q2, q3, q4])
    # Try to factor nicely
    print(f"delta_H = {delta_H}")

    # Show the 5-cycle correction explicitly
    print(f"\n5-cycle correction 2*(D5-C5) = {expand(2*term_5)}")
    print(f"3-cycle term -2*sum(s) = {expand(-2*(s2+s3+s4))}")
else:
    print(f"\nFAILED: difference = {diff}")
    # Debug
    print(f"delta_H = {delta_H}")
    print(f"delta_I = {delta_I}")
