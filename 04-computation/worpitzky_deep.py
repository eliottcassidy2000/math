#!/usr/bin/env python3
"""
DEEP DIVE: Worpitzky expansion of F(T,x)

Key observations from initial run:
1. Worpitzky coefficients are ALWAYS integers
2. sum(w_k) = 0 for all non-transitive tournaments
3. For transitive: w = [A(n-1,0), A(n-1,1), ..., A(n-1,n-2), 0] (Eulerian numbers!)
4. w_{n-1} = 0 always

The sum(w_k) = F(T,1) - n! property needs investigation.
Actually: the Worpitzky basis at x=1 gives C(1+k, n-1).
sum_k w_k C(1+k, n-1) = F(T,1) = H(T).

Wait, let me think about what's actually going on.

Worpitzky identity: x^{n-1} = sum_{k=0}^{n-2} A(n-1,k) C(x+k, n-1)

So for the "all permutations" case (n! permutations, each weighted by x^{fwd}):
sum_P x^{fwd(P)} = sum_k #{P with k forward edges} x^k

For the COMPLETE symmetric group S_n with the standard word statistic:
sum_{sigma} x^{des(sigma)} = A_{n}(x) = Eulerian polynomial

The Worpitzky identity says:
x^n = sum_k A(n,k) C(x+k, n)

This means x^n expanded in Worpitzky basis has Eulerian numbers as coefficients.

For our F(T,x), we're restricting to Hamiltonian paths of T (a subset of S_n),
and "forward edges" plays the role of "descents" (or ascents, depending on convention).

ACTUALLY: Let me reconsider. The standard Worpitzky identity relates
the monomial basis to the binomial basis. What I should investigate is:

HYPOTHESIS A: F(T,x) = sum_k w_k(T) C(x+k, n-1)
where w_k(T) are integers (CONFIRMED).

HYPOTHESIS B: The w_k(T) have a combinatorial interpretation.

HYPOTHESIS C: F(T,x) - F(transitive, x) has special structure in Worpitzky basis.

HYPOTHESIS D: The W-polynomial W(T,r) has a clean Worpitzky expansion.

Also: the connection to Ehrhart theory. If F(T,x) is the h*-vector of some
polytope, then Worpitzky non-negativity would follow from Ehrhart theory.
The ORDER POLYTOPE of a poset has h*-vector related to descent counts.
Tournaments induce partial orders (via transitive closure), so this might connect.
"""
from itertools import permutations
from math import comb, factorial
import random
import numpy as np

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

def transitive_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            A[i][j] = 1
    return A

def all_tournaments(n):
    """Generate all labeled tournaments on n vertices"""
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        yield A

def forward_edge_poly(A, n):
    F = [0] * n
    for P in permutations(range(n)):
        if not all(A[P[i]][P[i+1]] for i in range(n-1)):
            continue
        fwd = sum(1 for i in range(n-1) if P[i] < P[i+1])
        F[fwd] += 1
    return F

def worpitzky_basis_matrix(d):
    n = d + 1
    M = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            M[i][j] = comb(i + j, d)
    return M

def expand_in_worpitzky(F_coeffs):
    d = len(F_coeffs) - 1
    n = d + 1
    F_vals = [sum(F_coeffs[k] * x**k for k in range(n)) for x in range(n)]
    M = worpitzky_basis_matrix(d)
    w = np.linalg.solve(M, F_vals)
    return [int(round(x)) for x in w]

def eulerian_number(n, k):
    return sum((-1)**j * comb(n+1, j) * (k+1-j)**n for j in range(k+2))

random.seed(42)

# ===== EXHAUSTIVE AT n=4,5 =====
for n in [4, 5]:
    print(f"\n{'='*60}")
    print(f"EXHAUSTIVE n={n}")
    print(f"{'='*60}")

    w_distribution = {}
    for A in all_tournaments(n):
        F = forward_edge_poly(A, n)
        w = tuple(expand_in_worpitzky(F))
        H = sum(F)
        key = (H, w)
        w_distribution[key] = w_distribution.get(key, 0) + 1

    print(f"\nDistinct (H, w) pairs: {len(w_distribution)}")
    for (H, w), count in sorted(w_distribution.items()):
        print(f"  H={H:3d}, w={list(w)}, count={count}")

    # Check: is w always non-negative when shifted?
    # i.e., is there a universal shift c such that w_k + c >= 0?
    print(f"\n  Check alternating signs:")
    for (H, w), count in sorted(w_distribution.items()):
        signs = ['+' if x > 0 else ('-' if x < 0 else '0') for x in w]
        print(f"    H={H:3d}: signs={''.join(signs)}")

# ===== DEEPER: W-POLYNOMIAL IN WORPITZKY BASIS =====
print(f"\n\n{'='*60}")
print(f"W-POLYNOMIAL IN WORPITZKY BASIS")
print(f"{'='*60}")

# The W-polynomial W(T,r) = sum_P prod(r + s_e) where s_e = A[P_i][P_{i+1}] - 1/2
# This equals sum_k F_k (r+1/2)^k (r-1/2)^{n-1-k}
# Let's expand W(T,r) in the Worpitzky basis for polynomials in r.

def W_poly_coeffs(A, n):
    """Compute W(T,r) as polynomial in r. Return list of coefficients."""
    from numpy.polynomial import polynomial as P
    W = np.zeros(n)  # degree n-1, so n coefficients
    for perm in permutations(range(n)):
        if not all(A[perm[i]][perm[i+1]] for i in range(n-1)):
            continue
        # Product of (r + s_e) for each edge in the path
        poly = np.array([1.0])
        for i in range(n-1):
            s_e = A[perm[i]][perm[i+1]] - 0.5  # +0.5 or -0.5
            factor = np.array([s_e, 1.0])  # s_e + r
            poly = np.convolve(poly, factor)
        # Add to W
        for k in range(len(poly)):
            if k < n:
                W[k] += poly[k]
    return W

for n in [3, 4, 5]:
    print(f"\nn={n}:")
    A = random_tournament(n)
    W = W_poly_coeffs(A, n)
    print(f"  W(T,r) coefficients (constant to r^{n-1}): {[round(x,4) for x in W]}")

    # Expand in Worpitzky basis
    # Since W is degree n-1, use same approach
    d = n - 1
    # Evaluate W at r=0,1,...,d
    W_vals = [sum(W[k] * r**k for k in range(n)) for r in range(n)]
    M = worpitzky_basis_matrix(d)
    w = np.linalg.solve(M, W_vals)
    print(f"  Worpitzky coeffs: {[round(x,4) for x in w]}")

# ===== KEY INSIGHT: DIFFERENCE FROM TRANSITIVE =====
print(f"\n\n{'='*60}")
print(f"F(T,x) - F(transitive,x) IN WORPITZKY BASIS")
print(f"{'='*60}")

for n in [4, 5]:
    print(f"\nn={n}:")
    E = [eulerian_number(n-1, k) for k in range(n-1)] + [0]

    for A in all_tournaments(n):
        F = forward_edge_poly(A, n)
        w = expand_in_worpitzky(F)
        # Difference from transitive (Eulerian numbers)
        diff = [w[k] - E[k] for k in range(n)]
        H = sum(F)
        if H > 1 and any(d != 0 for d in diff):
            print(f"  H={H}, w={w}, diff_from_Eulerian={diff}")
            # Check: is diff a multiple of something?
            g = 0
            for d in diff:
                if d != 0:
                    from math import gcd
                    g = gcd(g, abs(d))
            if g > 0:
                print(f"    GCD of nonzero diffs: {g}, normalized: {[d//g for d in diff]}")

# ===== THE BIG QUESTION: What does w_k count? =====
print(f"\n\n{'='*60}")
print(f"COMBINATORIAL INTERPRETATION OF w_k")
print(f"{'='*60}")

# For the identity S_n: sum_sigma x^{des(sigma)} = A_n(x) = sum_k A(n,k) x^k
# Worpitzky: x^n = sum_k A(n,k) C(x+k,n)
# For tournaments: F(T,x) = sum_k F_k x^k = sum_k w_k C(x+k, n-1)
#
# From Worpitzky, w_k = sum_j F_j * [coefficient of C(x+k,n-1) in x^j]
# The inverse Worpitzky transform uses the "second kind Eulerian numbers" or
# Stirling-like numbers.
#
# Actually, the change of basis from {x^k} to {C(x+k,n-1)} involves
# evaluating both bases and inverting. Let me compute the change-of-basis matrix.

n = 5
d = n - 1  # degree 4
print(f"\nChange of basis matrix for n={n} (degree {d}):")
print(f"From monomial basis x^k to Worpitzky basis C(x+k,{d}):")

# M[i][j] = C(i+j, d) evaluated at x=i, for basis element j
M = worpitzky_basis_matrix(d)
print(f"  Basis matrix (eval points 0..{d}):\n{M.astype(int)}")
Minv = np.linalg.inv(M)
print(f"  Inverse (gives Worpitzky coeffs from values):\n{np.round(Minv, 4)}")

# The monomial x^j evaluated at x=0,...,d gives V = Vandermonde
V = np.zeros((d+1, d+1))
for i in range(d+1):
    for j in range(d+1):
        V[i][j] = i**j
print(f"\n  Vandermonde:\n{V.astype(int)}")

# Change of basis: w = Minv @ V @ f  where f = [F_0, ..., F_d]
# So the change-of-basis matrix is C = Minv @ V
C = Minv @ V
print(f"\n  Change of basis (monomial -> Worpitzky):")
print(f"  C = Minv @ V:\n{np.round(C, 4)}")

# Check against Eulerian numbers
print(f"\n  Eulerian numbers: {[eulerian_number(d,k) for k in range(d+1)]}")
# For x^d: the monomial x^d should map to Eulerian numbers
f_xd = [0]*d + [1]  # x^d = [0,...,0,1]
w_xd = [int(round(x)) for x in C @ f_xd]
print(f"  x^{d} in Worpitzky basis: {w_xd}")
print(f"  Should be Eulerian: {[eulerian_number(d,k) for k in range(d)] + [0]}")

# For x^0 = 1: should be [1,0,0,...,0] since C(0+0,d) = 0 for d>0... hmm
# C(0,4) = 0, C(1,4) = 0, ... only C(d,d) = 1.
# So 1 = w_0 C(0+0,d) + ... but C(k,d)=0 for k<d. Only C(d,d)=1.
# Wait, C(x+k,d) at x=0: C(k,d). This is 0 for k<d, 1 for k=d.
# So F(T,0) = w_d * C(d,d) = w_d.
# And F(T,0) = F_0 = #{HPs with 0 forward edges}.
print(f"\n  w_{{n-1}} should equal F_0:")
for trial in range(3):
    A = random_tournament(n)
    F = forward_edge_poly(A, n)
    w = expand_in_worpitzky(F)
    print(f"    F_0={F[0]}, w_{n-1}={w[n-1]}")

# ===== SHIFTED WORPITZKY: F(T,x+1) =====
print(f"\n\n{'='*60}")
print(f"SHIFTED F(T,x+1) AND WORPITZKY")
print(f"{'='*60}")
# Since F is palindromic, F(T,x) = x^{n-1} F(T, 1/x).
# Let's try expanding F(T, x+1) in Worpitzky basis.
# For the Eulerian polynomial: A_n(x+1) = sum_k A(n,k)(x+1)^k
# But Worpitzky says x^n = sum A(n,k) C(x+k,n)

for n in [4, 5]:
    print(f"\nn={n}:")
    for trial in range(5):
        A = random_tournament(n)
        F = forward_edge_poly(A, n)
        H = sum(F)

        # F(T, x+1): shift x -> x+1
        # F(T, x+1) = sum_k F_k (x+1)^k
        from numpy.polynomial import polynomial as P
        # Build polynomial coefficients of F(T, x+1)
        F_shifted = [0.0] * n
        for k in range(n):
            # F_k * (x+1)^k: expand using binomial theorem
            for j in range(k+1):
                F_shifted[j] += F[k] * comb(k, j)

        F_shifted_int = [int(round(x)) for x in F_shifted]
        w = expand_in_worpitzky(F_shifted_int)
        print(f"  Trial {trial}: H={H}, F(x+1)={F_shifted_int}, w={w}")
        print(f"    Non-negative? {all(x >= 0 for x in w)}")
