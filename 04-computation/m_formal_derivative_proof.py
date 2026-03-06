#!/usr/bin/env python3
"""
FORMAL DERIVATIVE / MATRIX IDENTITY APPROACH to M[a,b] = M[b,a]

Key idea: Express M[a,b] as a matrix operation on A = (t_{ij}),
then show the operation is symmetric.

APPROACH 1: Permanent derivative
  d/d(t_{ab}) per(A) = sum_{sigma: sigma(a)=b} prod_{i != a} t(i, sigma(i))
  This is the (a,b) "permanent cofactor" = per(A with row a, col b deleted).
  NOT equal to M[a,b] but related.

APPROACH 2: (I - A)^{-1} connection
  For a nilpotent A (self-loops = 0), (I-A)^{-1} = I + A + A^2 + ... + A^{n-1}.
  The (a,b) entry of (I-A)^{-1} counts walks from a to b.
  But M[a,b] involves Hamiltonian paths, not general walks.
  Connection: M[a,b] = (-1)^{n-2} * [coeff of prod_{i in U} lambda_i] in
  det(I - A + diag(lambda))? (Hafnian-like formula?)

APPROACH 3: Characteristic polynomial
  det(lambda I - A) = sum_k (-1)^k e_k(eigenvalues) lambda^{n-k}
  The coefficient of lambda^2 in det(lambda I - A) involves...
  M[a,b] appears as a SPECIFIC minor / cofactor structure?

APPROACH 4: Direct proof via multilinear expansion
  Each path weight = product of (r + s_e). Expand:
  M[a,b] = sum_{k=0}^{n-2} r^k * C_k
  where C_k = [r^k] sum_S (-1)^|S| sum_{paths} e_{n-2-k}(s_edges).

  Prove C_k = 0 for odd k by showing the sum has a sign-reversing involution.

kind-pasteur-2026-03-06-S24
"""

from itertools import permutations, combinations
from sympy import symbols, expand, Symbol, Matrix, det, Poly, factor
from collections import defaultdict

def make_tournament(n):
    """Set up symbolic c-tournament with n vertices."""
    r = Symbol('r')
    s = {}
    for i in range(n):
        for j in range(n):
            if i < j:
                s[(i,j)] = Symbol(f's{i}{j}')
                s[(j,i)] = -s[(i,j)]

    def t(i, j):
        if i == j: return 0
        return r + s[(i,j)]

    return r, s, t

def compute_M(a, b, n, r, s, t):
    """Compute M[a,b] via inclusion-exclusion."""
    U = [v for v in range(n) if v != a and v != b]
    total = 0
    for mask in range(2**len(U)):
        S = [U[i] for i in range(len(U)) if (mask >> i) & 1]
        R = [u for u in U if u not in S]
        sign = (-1)**len(S)
        Sa = set(S + [a])
        Rb = set(R + [b])

        ea = 0
        for p in permutations(sorted(Sa)):
            if p[-1] != a: continue
            w = 1
            for k in range(len(p)-1): w *= t(p[k], p[k+1])
            ea += w
        if len(Sa) == 1: ea = 1

        bb = 0
        for p in permutations(sorted(Rb)):
            if p[0] != b: continue
            w = 1
            for k in range(len(p)-1): w *= t(p[k], p[k+1])
            bb += w
        if len(Rb) == 1: bb = 1

        total += sign * ea * bb
    return expand(total)

# ==============================================================
# TEST 1: Characteristic polynomial connection
# ==============================================================
print("=" * 70)
print("TEST 1: Characteristic polynomial det(lambda*I - A)")
print("=" * 70)

for n in [4, 5]:
    r, s, t = make_tournament(n)
    lam = Symbol('lam')

    # Build A matrix
    A = Matrix(n, n, lambda i, j: t(i, j) if i != j else 0)

    # Characteristic polynomial
    char_mat = lam * Matrix.eye(n) - A
    char_poly = expand(det(char_mat))

    # Extract coefficient of lam^2
    p = Poly(char_poly, lam)

    # M[a,b] for a=0, b=1
    M_01 = compute_M(0, 1, n, r, s, t)
    M_10 = compute_M(1, 0, n, r, s, t)

    print(f"\nn={n}:")
    print(f"  M[0,1] = {M_01}")
    print(f"  M[1,0] = {M_10}")
    print(f"  M[0,1] - M[1,0] = {expand(M_01 - M_10)}")

    # Check various coefficients of char poly
    for k in range(n+1):
        coeff_k = expand(p.nth(k))
        # Check if M[0,1] divides or relates to any coefficient
        if k <= 3:
            print(f"  [lam^{k}] char_poly = {coeff_k}")

    # The (0,1) cofactor of (lam*I - A) is a poly in lam
    cofactor_01 = expand((-1)**(0+1) * det(char_mat.minor_submatrix(0, 1)))
    p_cof = Poly(cofactor_01, lam)

    print(f"  Cofactor(0,1) of (lam*I - A):")
    for k in range(p_cof.degree() + 1):
        ck = expand(p_cof.nth(k))
        print(f"    [lam^{k}] = {ck}")

    # Check: is M[0,1] the coefficient of lam^0 in cofactor(0,1)?
    cof_lam0 = expand(p_cof.nth(0))
    print(f"  [lam^0] cofactor = {cof_lam0}")
    print(f"  Matches M[0,1]? {expand(cof_lam0 - M_01) == 0}")
    print(f"  Matches -M[0,1]? {expand(cof_lam0 + M_01) == 0}")


# ==============================================================
# TEST 2: M as residue of resolvent (I - zA)^{-1}
# ==============================================================
print("\n" + "=" * 70)
print("TEST 2: Resolvent / generating function for walks")
print("=" * 70)

n = 4
r, s, t = make_tournament(n)
z = Symbol('z')
A = Matrix(n, n, lambda i, j: t(i, j) if i != j else 0)

# (I - zA)^{-1} = adj(I - zA) / det(I - zA)
# The (a,b) entry of adj(I - zA) = (-1)^{a+b} det((I-zA) with row b, col a deleted)

# For walks: (I - zA)^{-1}_{ab} = sum_k z^k A^k_{ab}
# = sum of walks from a to b of length k, weighted by z^k

# M[a,b] is the Hamiltonian path part. It's the coefficient of z^{n-2}
# in some modified version of the resolvent.

# Actually, the inclusion-exclusion in M[a,b] is an EXPONENTIAL formula
# that extracts the "connected" (Hamiltonian) part from the walk sum.

# Key: M[a,b] = [z^{n-2}] sum_S (-1)^|S| z^|S| E_a(S+a) * z^|R| B_b(R+b)?
# No, that doesn't work because z is not natural here.

# Let me instead test: does adj(I - A)_{ba} = M[a,b]?
# adj(I-A)_{ba} = (-1)^{a+b} det((I-A) with row a, col b deleted)

I_minus_A = Matrix.eye(n) - A
adj_ba = expand((-1)**(0+1) * det(I_minus_A.minor_submatrix(0, 1)))
M_01 = compute_M(0, 1, n, r, s, t)

print(f"\nn=4:")
print(f"  adj(I-A)[1,0] = {adj_ba}")
print(f"  M[0,1] = {M_01}")
print(f"  Match? {expand(adj_ba - M_01) == 0}")

# What about det(I - A)?
det_IA = expand(det(I_minus_A))
print(f"  det(I-A) = {det_IA}")

# Try: adj(I - A) itself (the full adjugate)
# adj(I-A)_{ij} = (-1)^{i+j} det(minor(I-A, j, i))
# Is adj(I-A) symmetric?
adj_10 = expand((-1)**(1+0) * det(I_minus_A.minor_submatrix(1, 0)))
print(f"  adj(I-A)[0,1] = {adj_10}")
print(f"  adj symmetric? {expand(adj_ba - adj_10) == 0}")


# ==============================================================
# TEST 3: The REAL connection - M as Hamilton path IE formula
# Each E_a(S+a) is a weighted sum of Ham paths ending at a.
# These paths can be expressed as permanents of submatrices.
# ==============================================================
print("\n" + "=" * 70)
print("TEST 3: E_a as permanent of path matrix")
print("=" * 70)

n = 5
r, s, t = make_tournament(n)

# E_a(W) = sum over Ham paths through W ending at a.
# For W = {u1, u2, ..., uk, a}, a Ham path is u_{sigma(1)} -> ... -> u_{sigma(k)} -> a
# with weight prod_{i=1}^{k-1} t(u_{sigma(i)}, u_{sigma(i+1)}) * t(u_{sigma(k)}, a)

# This is the permanent of a specific matrix:
# Consider the (k+1) x (k+1) matrix P where:
# P_{ij} = t(w_i, w_j) for the "path" from position i to position j.
# But a Hamiltonian path isn't a permutation matrix, it's a linear ordering.

# Actually, E_a(W) = per of the (k x k) matrix M' where:
# Rows = positions 1, ..., k (position in path, from start to penultimate)
# Cols = vertices W \ {a}
# M'_{i, u} = ... this doesn't quite work because each position depends on neighbors.

# E_a(W) is actually NOT a simple permanent. It's a permanent-like sum
# but each term is a product of CONSECUTIVE vertex pair weights.

# Let me instead verify: is the full M[a,b] for n=5 equal to
# something involving the adj(I-A)?

a, b = 0, 1
A = Matrix(n, n, lambda i, j: t(i, j) if i != j else 0)
I_minus_A = Matrix.eye(n) - A
M_01 = compute_M(0, 1, n, r, s, t)

# Test: cofactor of (I - A) at various positions
cof_01 = expand((-1)**(0+1) * det(I_minus_A.minor_submatrix(0, 1)))
cof_10 = expand((-1)**(1+0) * det(I_minus_A.minor_submatrix(1, 0)))

print(f"\nn=5:")
print(f"  M[0,1] = {M_01}")
print(f"  cof(I-A)[0,1] = {cof_01}")
print(f"  cof(I-A)[1,0] = {cof_10}")

# The cofactor has degree n-1 = 4, M has degree n-2 = 3.
# Can we extract M from the cofactor?

# Express cofactor as polynomial in r
p_cof = Poly(cof_01, r)
p_M = Poly(M_01, r)

print(f"\n  cofactor r-expansion:")
for k in range(p_cof.degree() + 1):
    print(f"    [r^{k}] = {expand(p_cof.nth(k))}")

print(f"\n  M[0,1] r-expansion:")
for k in range(p_M.degree() + 1):
    print(f"    [r^{k}] = {expand(p_M.nth(k))}")


# ==============================================================
# TEST 4: Direct [r^1] analysis at n=5 — term-by-term cancellation
# ==============================================================
print("\n" + "=" * 70)
print("TEST 4: [r^1] coefficient detailed structure at n=5")
print("=" * 70)

n = 5
r, s, t = make_tournament(n)
a, b = 0, 1
U = [v for v in range(n) if v != a and v != b]  # [2, 3, 4]

# For each (S, path_a, path_b), compute the r^1 contribution
# Group by the "s-monomial" (the product of s-values from all edges except the promoted one)
r1_by_monomial = defaultdict(list)
r1_total = 0

for mask in range(2**len(U)):
    S = [U[i] for i in range(len(U)) if (mask >> i) & 1]
    R = [u for u in U if u not in S]
    sign = (-1)**len(S)
    Sa = set(S + [a])
    Rb = set(R + [b])

    # All path pairs
    paths_a = []
    for p in permutations(sorted(Sa)):
        if p[-1] == a:
            edges = [(p[k], p[k+1]) for k in range(len(p)-1)]
            paths_a.append(edges)
    if len(Sa) == 1:
        paths_a = [[]]  # trivial path

    paths_b = []
    for p in permutations(sorted(Rb)):
        if p[0] == b:
            edges = [(p[k], p[k+1]) for k in range(len(p)-1)]
            paths_b.append(edges)
    if len(Rb) == 1:
        paths_b = [[]]

    for pa_edges in paths_a:
        for pb_edges in paths_b:
            all_edges = pa_edges + pb_edges
            m = len(all_edges)  # = n-2 = 3

            if m == 0: continue

            # r^1 coefficient: for each edge j, the product of s-values of all other edges
            for j in range(m):
                s_product = sign  # include the (-1)^|S| sign
                edge_labels = []
                for i in range(m):
                    if i == j:
                        continue  # this edge contributes r, not s
                    u, v = all_edges[i]
                    s_product *= s[(u,v)]
                    edge_labels.append(f"s{min(u,v)}{max(u,v)}" + ("" if u < v else "(-)"))

                promoted_edge = all_edges[j]
                key = (frozenset(all_edges) - {promoted_edge}, promoted_edge)

                s_product = expand(s_product)
                r1_total = expand(r1_total + s_product)

                # Record which S produced this
                r1_by_monomial[str(expand(s_product))].append({
                    'S': tuple(S), 'sign': sign,
                    'promoted': promoted_edge,
                    'other_edges': [e for i, e in enumerate(all_edges) if i != j],
                    'value': s_product
                })

print(f"\n  Total [r^1] = {r1_total}")
print(f"  (Should be 0: {r1_total == 0})")

# Group by s-monomial value and show cancellations
monomial_groups = defaultdict(list)
for key, entries in r1_by_monomial.items():
    for entry in entries:
        val = entry['value']
        monomial_groups[str(val)].append(entry)

print(f"\n  Number of distinct s-monomial groups: {len(monomial_groups)}")

# Show a few representative cancellations
shown = 0
for mono_str, entries in sorted(monomial_groups.items(), key=lambda x: len(x[1]), reverse=True):
    if shown >= 5: break
    total = sum(e['value'] for e in entries)
    if total != 0:
        print(f"\n  NON-CANCELING group: {mono_str}")
        for e in entries:
            print(f"    S={e['S']}, promoted={e['promoted']}, others={e['other_edges']}")
    shown += 1

# More informative: group by the s-variable monomial (up to sign)
print("\n  Detailed cancellation analysis:")
# Each r^1 term is sign * s_e1 * s_e2 (product of 2 s-values, since m=3, pick 1 for r)
# Let's collect by the PAIR of s-indices
pair_contributions = defaultdict(list)
for mask in range(2**len(U)):
    S = [U[i] for i in range(len(U)) if (mask >> i) & 1]
    R = [u for u in U if u not in S]
    sign = (-1)**len(S)
    Sa = set(S + [a])
    Rb = set(R + [b])

    paths_a = []
    for p in permutations(sorted(Sa)):
        if p[-1] == a:
            paths_a.append([(p[k], p[k+1]) for k in range(len(p)-1)])
    if len(Sa) == 1:
        paths_a = [[]]

    paths_b = []
    for p in permutations(sorted(Rb)):
        if p[0] == b:
            paths_b.append([(p[k], p[k+1]) for k in range(len(p)-1)])
    if len(Rb) == 1:
        paths_b = [[]]

    for pa_edges in paths_a:
        for pb_edges in paths_b:
            all_edges = pa_edges + pb_edges
            m = len(all_edges)
            if m == 0: continue

            for j in range(m):
                # The two s-edges (not promoted)
                s_edges = [all_edges[i] for i in range(m) if i != j]
                promoted = all_edges[j]

                s_val = sign
                for u, v in s_edges:
                    s_val *= s[(u,v)]

                pair_key = tuple(sorted([
                    (min(e), max(e)) for e in s_edges
                ]))
                pair_contributions[pair_key].append({
                    'S': tuple(S), 'sign': sign,
                    'promoted': promoted,
                    's_edges': s_edges,
                    'value': expand(s_val)
                })

print(f"\n  Edge-pair groups: {len(pair_contributions)}")
for pair_key in sorted(pair_contributions.keys()):
    entries = pair_contributions[pair_key]
    total = expand(sum(e['value'] for e in entries))
    if total != 0:
        print(f"\n  NONZERO pair {pair_key}: total = {total}")
        for e in entries:
            print(f"    S={e['S']}, sign={e['sign']:+d}, promoted={e['promoted']}, "
                  f"s_edges={e['s_edges']}, value={e['value']}")
    else:
        # Count how many terms cancel
        n_terms = len(entries)
        if n_terms <= 4:
            print(f"\n  Zero pair {pair_key} ({n_terms} terms):")
            for e in entries:
                print(f"    S={e['S']}, sign={e['sign']:+d}, promoted={e['promoted']}, "
                      f"s_edges={e['s_edges']}, value={e['value']}")

print("\n" + "=" * 70)
print("TEST 5: The KEY identity - does M relate to a SYMMETRIC matrix function?")
print("=" * 70)

# M[a,b] for ALL pairs at n=4, check if M matrix is symmetric
n = 4
r, s, t = make_tournament(n)
M_matrix = {}
for a in range(n):
    for b in range(n):
        if a != b:
            M_matrix[(a,b)] = compute_M(a, b, n, r, s, t)

print(f"\nn=4: Full transfer matrix:")
for a in range(n):
    for b in range(n):
        if a != b:
            val = M_matrix[(a,b)]
            print(f"  M[{a},{b}] = {val}")

# Check symmetry
sym_check = True
for a in range(n):
    for b in range(a+1, n):
        diff = expand(M_matrix[(a,b)] - M_matrix[(b,a)])
        if diff != 0:
            print(f"  M[{a},{b}] - M[{b},{a}] = {diff} != 0")
            sym_check = False
print(f"\n  Full matrix symmetric? {sym_check}")

# Check: is M = adj(something)?
# If M = adj(L) for some matrix L, then M is symmetric iff adj(L) is symmetric
# iff L is symmetric (for invertible L) or L has special structure.

# What if M = adj(S) where S is the skew-symmetric matrix (s_{ij})?
S_mat = Matrix(n, n, lambda i, j: s[(i,j)] if i != j else 0)
adj_S = {}
for a in range(n):
    for b in range(n):
        if a != b:
            adj_S[(a,b)] = expand((-1)**(a+b) * det(S_mat.minor_submatrix(b, a)))

print(f"\nn=4: adj(S) vs M at r=0:")
for a in range(n):
    for b in range(a+1, n):
        M_at_0 = expand(M_matrix[(a,b)].subs(r, 0))
        adj_at = adj_S[(a,b)]
        match = expand(M_at_0 - adj_at) == 0
        if not match:
            ratio_check = expand(M_at_0 + adj_at) == 0
            print(f"  M[{a},{b}](0) = {M_at_0}, adj(S)[{a},{b}] = {adj_at}, "
                  f"match={match}, neg_match={ratio_check}")
        else:
            print(f"  M[{a},{b}](0) = adj(S)[{a},{b}] = {M_at_0}")

# What is adj(S) exactly? For skew-symmetric S, adj(S) = Pf(S_{ab}) * something?
print(f"\n  adj(S) matrix:")
for a in range(n):
    for b in range(n):
        if a != b:
            print(f"  adj(S)[{a},{b}] = {adj_S[(a,b)]}")

# Check if adj(S) is symmetric
adj_sym = True
for a in range(n):
    for b in range(a+1, n):
        if expand(adj_S[(a,b)] - adj_S[(b,a)]) != 0:
            adj_sym = False
            print(f"  adj(S)[{a},{b}] - adj(S)[{b},{a}] = {expand(adj_S[(a,b)] - adj_S[(b,a)])}")
print(f"  adj(S) symmetric? {adj_sym}")

# For a 4x4 skew-symmetric matrix, adj(S) = Pf(S) * S^{-1} (if S invertible)
# and Pf(S) * S^{-1} = Pf(S_{ij}) (signed Pfaffian minors)
# Since S is skew-symmetric, adj(S)^T = adj(S^T) = adj(-S) = (-1)^{n-1} adj(S)
# For n=4: adj(-S) = (-1)^3 adj(S) = -adj(S). So adj(S)^T = -adj(S)!
# This means adj(S) is ANTISYMMETRIC for n=4!

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
