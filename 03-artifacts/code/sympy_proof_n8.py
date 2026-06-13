#!/usr/bin/env python3
"""
Symbolic proof of OCF at n=8 using SymPy.

27 arc variables, 40320 permutations. Expected to take hours.

Uses the FULL A-clique formula (not the simplified n<=7 version):
  delta_I = 2 * sum_C [gained(C) - lost(C)] * H(comp(C))

The simplified formula -2*sum(s_x*H(B_x)) + 2*(D5-C5) + 2*(D7-C7)
FAILS at n=8 because 5-cycle complements have 3 vertices and H(3-vert)
can be 1 or 3 (not always 1 as at n<=7).

Instance: kind-pasteur-2026-03-05-S7 (fixed kind-pasteur-2026-03-05-S8)
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


print(f"=== n={n} Symbolic Proof (Full A-clique Formula) ===")
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

# === FULL A-CLIQUE FORMULA ===
# delta_I = 2 * sum_C [gained(C) - lost(C)] * H(comp(C))
# where C ranges over all odd cycle lengths L = 3, 5, 7
# and comp(C) = others \ (vertices of C besides i,j)

print("Computing delta_I via full A-clique formula...")

def h_sub(verts):
    """Compute H(T[verts]) as symbolic polynomial."""
    if len(verts) <= 1:
        return 1
    total = 0
    for perm in permutations(verts):
        w = 1
        for k in range(len(perm) - 1):
            w *= T(perm[k], perm[k+1])
        total += w
    return expand(total)

delta_I = 0

for L in range(3, n + 1, 2):  # odd cycle lengths: 3, 5, 7
    num_inter = L - 2  # intermediate vertices (from others)
    print(f"  Processing {L}-cycles ({num_inter} intermediate vertices)...")
    t1 = time.time()
    cycle_count = 0

    for subset in combinations(others, num_inter):
        complement = [v for v in others if v not in subset]
        H_comp = h_sub(complement)

        for perm in permutations(subset):
            # Lost cycle: i=0 -> j=1 -> perm[0] -> ... -> perm[-1] -> i=0
            # Uses arc i->j (which exists in T, not in T')
            # Arcs to check: T[J, perm[0]], T[perm[k], perm[k+1]], T[perm[-1], I]
            lost_w = T(J, perm[0])
            for k in range(len(perm) - 1):
                lost_w *= T(perm[k], perm[k+1])
            lost_w *= T(perm[-1], I)

            # Gained cycle: j=1 -> i=0 -> perm[0] -> ... -> perm[-1] -> j=1
            # Uses arc j->i (which exists in T', not in T)
            # Arcs to check: T[I, perm[0]], T[perm[k], perm[k+1]], T[perm[-1], J]
            gained_w = T(I, perm[0])
            for k in range(len(perm) - 1):
                gained_w *= T(perm[k], perm[k+1])
            gained_w *= T(perm[-1], J)

            delta_I += 2 * (gained_w - lost_w) * H_comp
            cycle_count += 1

    print(f"    {cycle_count} cycle patterns ({time.time()-t1:.1f}s)")

print(f"Expanding delta_I...")
t1 = time.time()
delta_I = expand(delta_I)
n_terms_I = len(delta_I.as_coefficients_dict())
print(f"  delta_I: {n_terms_I} monomials ({time.time()-t1:.1f}s)")

print(f"\nChecking delta_H == delta_I...")
diff = expand(delta_H - delta_I)

elapsed = time.time() - start
if diff == 0:
    print(f"\n*** PROVED: H(T)-H(T') = delta_I as POLYNOMIAL IDENTITY at n={n} ***")
    print(f"    (Computed in {elapsed:.1f}s)")
    print(f"    {n_vars} arc variables, {n_terms_H} monomial terms")
    print(f"    Formula: full A-clique: 2*sum_C [gained-lost]*H(comp(C))")
else:
    nd = len(diff.as_coefficients_dict())
    print(f"\nDifference has {nd} terms (computed in {elapsed:.1f}s)")
    if nd <= 10:
        print(f"diff = {diff}")
