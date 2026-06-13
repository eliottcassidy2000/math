#!/usr/bin/env python3
"""
cross_domain_connections.py — Exploring F(T,x) through lenses from different
areas of mathematics, looking for unexpected isomorphisms.

IDEA: F(T,x) is palindromic of degree n-1. What other objects in mathematics
produce palindromic polynomials with similar structural properties?

CONNECTIONS TO EXPLORE:
1. Ehrhart theory: F(T,x) as h*-vector
2. Hilbert series: F(T,x) as h-vector of a Cohen-Macaulay ring
3. Poincare duality: palindrome ↔ Poincare duality on a "tournament manifold"
4. Zeta polynomial of a poset: does the tournament poset have a zeta polynomial related to F?
5. Species / exponential structures
6. Hard-core lattice gas: F(T,x) as partition function
7. Kazhdan-Lusztig: palindromic KL polynomials
8. Jones polynomial: palindromic knot invariants
9. Fibonacci cubes: palindromic rank sequences

Author: opus-2026-03-07-S46
"""
from itertools import permutations, combinations
from math import comb, factorial
import numpy as np

def tournament_from_bits(n, bits):
    adj = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            idx += 1
    return adj

def compute_F(adj, n):
    F = [0]*n
    for P in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if adj[P[i]][P[i+1]])
        F[fwd] += 1
    return F

# ============================================================
# CONNECTION 1: F(T,x) AS CHROMATIC-LIKE POLYNOMIAL
# ============================================================
print("=" * 60)
print("CONNECTION 1: F(T,x) AS CHROMATIC-LIKE INVARIANT")
print("=" * 60)

# The chromatic polynomial chi_G(k) counts proper k-colorings.
# For the complete graph K_n: chi(k) = k(k-1)...(k-n+1) = falling factorial
# Worpitzky of chi gives Stirling numbers.
#
# For tournaments, what is the "coloring" interpretation of a_m?
# a_m = sum_k F_k C(m+n-1-k, n-1)
#
# C(m+n-1-k, n-1) = C(m+n-1-k, m-k) = # of multisets of size m-k from [n-1]
# (or lattice paths, or balls-in-bins)
#
# So a_m = sum_k F_k * (# weak compositions of m-k into n parts)
# This is the number of ways to:
# 1. Choose a permutation P of [n] contributing x^k to F
# 2. Choose a weak composition (c_0,...,c_{n-1}) with sum m-k
# 3. "Assign" the composition to the permutation

# What does this COUNT combinatorially?
# A_m might count (P, c) pairs where P is a permutation and c is a
# weak composition satisfying some condition...

# Actually, sum_k F_k C(m+n-1-k, n-1) counts the number of words
# w_1 w_2 ... w_n where w_i ∈ {1,...,m} and w_{P[i]} ≤ w_{P[i+1]}
# whenever P[i]→P[i+1] in T... no, that's not quite right.

# The STANDARD interpretation (from Worpitzky): a_m counts
# the number of weakly increasing sequences 1 ≤ w_1 ≤ ... ≤ w_n ≤ m+1
# where w_i < w_{i+1} at certain positions determined by the permutation.
# But here we have tournament structure, not just descent structure.

# Alternative: P-partition theory! A P-partition is a map f: [n] → Z_+
# such that f(i) ≥ f(j) if i >_P j (and strict if i covers j from left).
# The order polynomial Ω_P(m) counts P-partitions with values in {1,...,m}.
# For a NATURALLY labeled poset, Ω_P(m) = sum_{σ compatible} C(m+n-1-des(σ), n)
# where des counts descents relative to the linear extension.

# For a tournament T, define the "tournament order" i >_T j iff i→j.
# This is NOT a poset (not transitive), but we can still define a
# "tournament partition": f: [n] → {1,...,m} with f(i) ≥ f(j) if i→j.
# Since T is a tournament, this is a TOTAL ordering constraint on each pair.
# f(i) ≥ f(j) if i→j AND f(j) ≥ f(i) if j→i.
# If i→j and j→i... wait, in a tournament exactly one holds.
# So: for each i,j with i→j, we need f(i) ≥ f(j).
# This defines f as a "weak homomorphism" from T to the chain [m].

# The number of such f is exactly the order polynomial of the
# transitive closure of T! But wait — T may have cycles, so the
# transitive closure is the complete graph, making f constant.
# Unless T is acyclic (which tournaments are iff they're transitive).

# So this interpretation only works for the transitive tournament.
# For cyclic tournaments, we need to allow "violations" — exactly what
# F(T,x) counts!

# ACTUAL INTERPRETATION: a_m counts surjections f: [n] → {0,...,m} such that
# the number of "forward steps" (consecutive elements i, f(i) in the
# permutation ordering) equals... hmm, this is getting complicated.

# Let me just verify the P-partition connection for transitive tournament.
n = 5
adj_trans = [[1 if i < j else 0 for j in range(n)] for i in range(n)]
F_trans = compute_F(adj_trans, n)
print(f"\nTransitive T5: F = {F_trans}")

# P-partition of the chain 1 < 2 < 3 < 4 < 5:
# f: [5] → {1,...,m} with f(1) ≤ f(2) ≤ f(3) ≤ f(4) ≤ f(5)
# This is C(m+4, 4) = C(m+n-1, n-1) ... no, it's C(m+n-1, n-1)
# for weakly increasing sequences, which is DIFFERENT from a_m.
# a_m for transitive = (m+1)^5 - m^5.
# C(m+4, 4) for m=0: 1, m=1: 5, m=2: 15, m=3: 35
# a(0)=1, a(1)=31, a(2)=211, a(3)=781
# These are very different!

# So the P-partition interpretation doesn't directly apply.
# But there IS a connection through Stanley's theory of P-partitions
# and the DESCENT polynomial of a poset, which generalizes to tournaments
# via the "tournament descent polynomial" = F(T,x).

print("\nP-partition counts (weakly increasing f:[n]->{1..m}) for chain:")
for m in range(1, 6):
    pp = comb(m + n - 1, n - 1)
    print(f"  m={m}: P-partition count = {pp}, a_m = {(m+1)**n - m**n}")

# ============================================================
# CONNECTION 2: BERNOULLI NUMBERS AND POWER SUMS
# ============================================================
print("\n" + "=" * 60)
print("CONNECTION 2: BERNOULLI NUMBERS")
print("=" * 60)

# For transitive tournament: a_m = (m+1)^n - m^n
# sum_{m=0}^{M-1} a_m = M^n (telescoping!)
# So sum_{m=0}^{M-1} a_m = M^n = 1^n + ... well, no.
# (1^n - 0^n) + (2^n - 1^n) + ... + (M^n - (M-1)^n) = M^n. Trivially.
# But if we use the Bernoulli number formula for power sums:
# sum_{k=1}^{M} k^{n-1} = B_n(M+1)/n - B_n(1)/n (Bernoulli polynomial)
#
# For general tournaments: sum_{m=0}^{M-1} a_m(T) = ?
# This should ALSO telescope if we think of a_m as a forward difference.

n = 5
for bits in [0]:  # transitive
    adj = tournament_from_bits(n, bits)
    F = compute_F(adj, n)

    # Check telescoping
    from math import comb as C
    def worpitzky_a(F, n, m):
        return sum(F[k] * C(m + n - 1 - k, n - 1) for k in range(n))

    print(f"\nn={n}: F={F}")
    running = 0
    for m in range(6):
        a = worpitzky_a(F, n, m)
        running += a
        print(f"  sum a_0..a_{m} = {running}, (m+1)^n = {(m+1)**n}")

# For non-transitive tournament:
for bits in [7]:  # some non-transitive
    adj = tournament_from_bits(n, bits)
    F = compute_F(adj, n)
    print(f"\nn={n}: F={F}")
    running = 0
    for m in range(6):
        a = worpitzky_a(F, n, m)
        running += a
        print(f"  sum a_0..a_{m} = {running}")

# ============================================================
# CONNECTION 3: MAHONIAN STATISTICS
# ============================================================
print("\n" + "=" * 60)
print("CONNECTION 3: MAHONIAN STATISTICS AND INVERSIONS")
print("=" * 60)

# The inversion polynomial I(T,q) = sum_P q^{inv_T(P)}
# where inv_T(P) = #{(i,j): i<j, P[i]→P[j]}
# = # forward edges in the arrangement ordering, not permutation ordering
#
# For the identity permutation, inv_T(id) = # edges i→j with i<j
# = # upper-triangular 1s in adjacency matrix
#
# QUESTION: How does I(T,q) relate to F(T,x)?
# They are DIFFERENT statistics — F counts consecutive forward edges,
# I counts ALL forward edges (not just consecutive).
#
# For the transitive tournament (i→j iff i<j):
# inv_T(P) = # pairs (i,j) with i<j and P[i]<P[j] = coinversions = C(n,2) - inv(P)
# So I(transitive, q) = q^{C(n,2)} * sum_P q^{-inv(P)} = sum_P q^{coinv(P)}
# = [n]_q! (q-factorial)
#
# F(transitive, x) = Eulerian polynomial A_n(x)
#
# [n]_q! and A_n(x) are related by: [n]_q! = A_n(q) evaluated differently?
# No, they're different polynomials entirely.
# [5]! = 1+q+q^2+... = (1+q)(1+q+q^2)(1+q+q^2+q^3)(1+q+q^2+q^3+q^4)

n = 5
adj_trans = [[1 if i < j else 0 for j in range(n)] for i in range(n)]

# Compute F and I for transitive tournament
F_trans = compute_F(adj_trans, n)
I_trans = [0] * (n*(n-1)//2 + 1)
for P in permutations(range(n)):
    inv = sum(1 for i in range(n) for j in range(i+1, n) if adj_trans[P[i]][P[j]])
    I_trans[inv] += 1

print(f"\nTransitive T5:")
print(f"  F(T,x) = {F_trans} (degree {n-1})")
print(f"  I(T,q) = {I_trans} (degree {n*(n-1)//2})")

# Now for a non-transitive tournament
for bits in [7, 15, 31]:
    adj = tournament_from_bits(n, bits)
    F = compute_F(adj, n)
    I = [0] * (n*(n-1)//2 + 1)
    for P in permutations(range(n)):
        inv = sum(1 for i in range(n) for j in range(i+1, n) if adj[P[i]][P[j]])
        I[inv] += 1

    if tuple(F) != tuple(F_trans):
        print(f"\n  bits={bits}: F={F}")
        print(f"  I = {I}")
        # Is I also palindromic?
        is_pal = all(I[k] == I[len(I)-1-k] for k in range(len(I)//2 + 1))
        print(f"  I palindromic? {is_pal}")
        break

# ============================================================
# CONNECTION 4: TRANSFER MATRIX AND MARKOV CHAINS
# ============================================================
print("\n" + "=" * 60)
print("CONNECTION 4: TRANSFER MATRIX / MARKOV CHAIN")
print("=" * 60)

# Think of a tournament T as defining a Markov chain:
# Transition matrix M[i][j] = 1/(n-1) if i→j, 0 otherwise
# (each vertex transitions uniformly to one of its out-neighbors)
# Wait, that's not right — tournaments have exactly one direction per edge,
# but out-degree varies.
#
# Better: Use the adjacency matrix A. The matrix x*A + (1-x)*(J-I-A)
# = x*A + (1-x)*A^T (for tournaments where A^T = J-I-A)
# = (2x-1)*A + (1-x)*(J-I)
#
# At x=1: this is A (follow forward edges)
# At x=0: this is A^T (follow backward edges)
# At x=1/2: this is (J-I)/2 (uniform on all neighbors)
#
# The "transfer matrix" for F(T,x) counts Hamiltonian paths weighted by x.
# F(T,x) = sum_P prod_{i=0}^{n-2} [x*A + (1-x)*A^T]_{P[i],P[i+1]}
#         = sum_P prod_i (x if A[P[i]][P[i+1]]=1 else (1-x))
# Wait, that's not right either. F(T,x) = sum_P x^{fwd(P)}
# = sum_P x^{fwd} * 1^{back} where back = (n-1) - fwd.
# So each forward step contributes x and each backward step contributes 1.
# But then F(T,1) = n! (every step contributes 1), not just H(T).
# Actually F(T,1) = sum_P 1^{fwd} = n!, yes.

# The "weighted adjacency matrix" is W(x)[i][j] = x if i→j, 1 if j→i (i≠j)
# = x*A[i][j] + A^T[i][j] = x*A[i][j] + (1-A[i][j]) for i≠j
# = 1 + (x-1)*A[i][j] for i≠j
# So W(x) = (J-I) + (x-1)*A

# F(T,x) is the PERMANENT of the Hamiltonian path matrix?
# No, F(T,x) = sum over all n! orderings of prod of consecutive W entries.
# This is like a "path permanent" — not the usual matrix permanent.

# The transfer matrix formulation:
# Let T_x be the n×n matrix with T_x[i][j] = x if i→j, 1 if j→i, 0 if i=j.
# Then F(T,x) = sum_{paths i_0→i_1→...→i_{n-1}} prod T_x[i_k][i_{k+1}]
# where the sum is over all Hamiltonian paths (visiting each vertex exactly once).
#
# This is NOT the matrix permanent (which sums over all matchings, not paths).
# But it IS related to the permanent of a modified matrix...

# Actually, the Hamiltonian path sum IS related to the permanent via
# the result of Irving-Omar (THM-082): H(T) = sum_S det(A^T[S]) per(A[S^c])

# Let's compute det(W(x)) — the determinant of the weighted adjacency matrix
n = 5
print(f"\nn={n}: det(W(x)) for various tournaments")
seen = set()
for bits in range(1 << (n*(n-1)//2)):
    adj = tournament_from_bits(n, bits)
    F = compute_F(adj, n)
    key = tuple(F)
    if key in seen:
        continue
    seen.add(key)

    # det(W(x)) = det(J-I + (x-1)A) at several x values
    A = np.array(adj, dtype=float)
    J = np.ones((n,n))
    I = np.eye(n)

    det_vals = {}
    for x_val in [0, 0.5, 1, 2, -1]:
        W = (J - I) + (x_val - 1) * A
        det_vals[x_val] = np.linalg.det(W)

    print(f"  F={F}")
    print(f"  det(W(x)) at x=0,0.5,1,2,-1: {[f'{det_vals[x]:.1f}' for x in [0, 0.5, 1, 2, -1]]}")

# ============================================================
# CONNECTION 5: COXETER GROUPS AND DESCENT ALGEBRAS
# ============================================================
print("\n" + "=" * 60)
print("CONNECTION 5: DESCENT ALGEBRA")
print("=" * 60)

# In the symmetric group S_n, the descent set of a permutation σ is
# Des(σ) = {i : σ(i) > σ(i+1)}.
# The descent algebra (Solomon's descent algebra) has basis elements
# D_S = sum_{Des(σ)=S} σ for each subset S ⊆ [n-1].
#
# For a tournament T, we can define the "T-descent set":
# Des_T(σ) = {i : σ(i) →_T σ(i+1)} = positions of forward edges
# This is NOT a subset of [n-1] in the usual sense, but it IS a function
# from permutations to subsets of {0,...,n-2}.
#
# F_k(T) = #{σ : |Des_T(σ)| = k}
#
# QUESTION: Does {D_T := sum_{Des_T(σ)=S} σ : S ⊆ {0,...,n-2}} form a subalgebra?
# Almost certainly not for general T, but maybe for certain tournaments?

# Let's check: for the transitive tournament, Des_T(σ) = complement of descent set
# (since σ(i) <_T σ(i+1) iff σ(i) < σ(i+1)).
# So Des_T(σ) = Asc(σ) = [n-1] \ Des(σ).
# The ascent algebra is isomorphic to the descent algebra.

# For the regular tournament (n=5, all scores 2):
# Des_T is a more exotic function.

# Instead of the full algebra, let's look at the element
# R_T(x) = sum_σ x^{|Des_T(σ)|} * σ
# This is an element of the group algebra Q[S_n][x].
# Its "character" is F(T,x) (evaluating σ at 1).
#
# But we can also evaluate at other representations!
# If we take the sign character: sum_σ sgn(σ) x^{|Des_T(σ)|}
# This gives a signed version of F(T,x).

n = 5
print(f"\nn={n}: Signed F(T,x) = sum sgn(σ) x^{{fwd_T(σ)}}")
seen = set()
for bits in range(1 << (n*(n-1)//2)):
    adj = tournament_from_bits(n, bits)
    F = compute_F(adj, n)
    key = tuple(F)
    if key in seen:
        continue
    seen.add(key)

    # Signed version
    SF = [0]*n
    for P in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if adj[P[i]][P[i+1]])
        inv = sum(1 for i in range(n) for j in range(i+1, n) if P[i] > P[j])
        sgn = (-1)**inv
        SF[fwd] += sgn

    print(f"  F={F}, SF={SF}")

# ============================================================
# CONNECTION 6: QUANTUM GROUPS / q-ANALOGUE
# ============================================================
print("\n" + "=" * 60)
print("CONNECTION 6: q-ANALOGUE OF F(T,x)")
print("=" * 60)

# Define F(T,x,q) = sum_P x^{fwd(P)} * q^{inv(P)}
# where inv(P) is the usual inversion count of the permutation.
# At q=1: F(T,x,1) = F(T,x)
# At x=1: F(T,1,q) = [n]_q! (q-factorial, independent of T? No!)
# Actually at x=1: F(T,1,q) = sum_P q^{inv(P)} = [n]_q! regardless of T.
# Because we're summing over ALL permutations.
# Wait, is that true? sum_P q^{inv(P)} = product_{k=1}^n [k]_q.
# This is indeed independent of T (inv counts inversions in P, not in T).

# So F(T,x,q) is a bivariate polynomial that interpolates between
# F(T,x) (at q=1) and [n]_q! (at x=1).

n = 4
print(f"\nn={n}: F(T,x,q) = sum x^fwd * q^inv")

# Compute for a few tournaments and display as grid
seen = set()
for bits in range(1 << (n*(n-1)//2)):
    adj = tournament_from_bits(n, bits)
    F = compute_F(adj, n)
    key = tuple(F)
    if key in seen:
        continue
    seen.add(key)

    # F(T,x,q) as 2D array indexed by (fwd, inv)
    Fxq = [[0] * (n*(n-1)//2 + 1) for _ in range(n)]
    for P in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if adj[P[i]][P[i+1]])
        inv = sum(1 for i in range(n) for j in range(i+1, n) if P[i] > P[j])
        Fxq[fwd][inv] += 1

    print(f"\n  F={F}")
    print(f"  F(T,x,q) table (rows=fwd, cols=inv):")
    for k in range(n):
        row = Fxq[k]
        if sum(row) > 0:
            nonzero = [(j, row[j]) for j in range(len(row)) if row[j] > 0]
            print(f"    x^{k}: {nonzero}")

    # Check: is each row palindromic in q?
    for k in range(n):
        row = Fxq[k]
        max_inv = max((j for j in range(len(row)) if row[j] > 0), default=0)
        is_pal = all(row[j] == row[max_inv - j] for j in range(max_inv//2 + 1))
        if sum(row) > 0:
            pass  # We'll check below

    # Marginals
    q_marginal = [sum(Fxq[k][j] for k in range(n)) for j in range(n*(n-1)//2 + 1)]
    x_marginal = [sum(row) for row in Fxq]
    print(f"    q-marginal (inv dist): {[q_marginal[j] for j in range(len(q_marginal)) if q_marginal[j]>0]}")
    print(f"    x-marginal (fwd dist): {x_marginal} = F(T,x)")
