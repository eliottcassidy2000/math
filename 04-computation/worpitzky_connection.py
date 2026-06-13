#!/usr/bin/env python3
"""
WORPITZKY POLYNOMIALS AND TOURNAMENT FORWARD-EDGE POLYNOMIALS

Worpitzky's identity: x^n = sum_{k=0}^{n-1} A(n,k) * C(x+k, n)

where A(n,k) = Eulerian number = #{perm of [n] with k descents}.

Equivalently: (x+1)^n = sum_{k=0}^{n-1} A(n,k) * C(x+n-k, n)

The Eulerian polynomial: A_n(t) = sum_k A(n,k) * t^k

Key properties of A_n(t):
1. Palindromic: A(n,k) = A(n,n-1-k)
2. A_n(1) = n!
3. A_n(t) = (1-t)^{n+1} * sum_{j>=0} (j+1)^n * t^j

CONNECTION TO TOURNAMENTS:
Our F(T,x) = sum_{k=0}^{n-1} F_k * x^k where F_k = #{HPs with k forward edges}.
F(T,x) is palindromic: F_k = F_{n-1-k}.
F(T,1) = H(T) = total HP count.

For the TRANSITIVE tournament: every HP is the unique topological sort,
so F_k = 0 for k != n-1 and F_{n-1} = 1. (Only 1 HP, all edges forward.)
Wait, that's not right — transitive has exactly 1 HP.

For a RANDOM permutation sigma of [n], the HP is sigma, and the number of
"forward edges" (descents in the permutation viewed as a sequence) equals
the descent count. So for the COMPLETE TOURNAMENT (all edges both ways = not
a tournament), the forward-edge polynomial would be the Eulerian polynomial.

But tournaments have exactly ONE direction per edge. So F(T,x) for a
specific tournament T is a "restricted descent polynomial" — it counts
only those permutations that are Hamiltonian paths of T, weighted by descents.

HYPOTHESIS: F(T,x) expanded in the Worpitzky basis C(x+k, n-1) gives
coefficients that are always non-negative (or always alternate in sign,
or have some other structural property).

The Worpitzky basis for degree n-1 polynomials is:
{C(x, n-1), C(x+1, n-1), ..., C(x+n-2, n-1)}

So F(T,x) = sum_{k=0}^{n-2} w_k(T) * C(x+k, n-1)

Let's compute w_k(T) for various tournaments and see what patterns emerge.
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

def forward_edge_poly(A, n):
    """Compute F_k = #{HPs with exactly k forward edges}"""
    F = [0] * n  # F[k] for k=0,...,n-1
    for P in permutations(range(n)):
        if not all(A[P[i]][P[i+1]] for i in range(n-1)):
            continue
        fwd = sum(1 for i in range(n-1) if P[i] < P[i+1])
        F[fwd] += 1
    return F

def worpitzky_basis_value(k, x, d):
    """C(x+k, d) — the k-th Worpitzky basis polynomial of degree d evaluated at x"""
    return comb(x + k, d)

def worpitzky_basis_matrix(d):
    """Matrix M where M[i][j] = C(j + i, d) for i=0..d, j=0..d
    (evaluating basis poly i at x=j)
    Actually we need to evaluate at enough points to solve for coefficients.
    F(T,x) is degree n-1, and basis has n elements, so we evaluate at x=0,1,...,n-1.
    """
    n = d + 1  # number of basis elements = d+1
    M = np.zeros((n, n))
    for i in range(n):  # evaluation point x=i
        for j in range(n):  # basis element j: C(x+j, d)
            M[i][j] = comb(i + j, d)
    return M

def expand_in_worpitzky(F_coeffs):
    """Given F(T,x) = sum F_k x^k, expand in Worpitzky basis C(x+j, d) where d=len(F)-1.
    Returns w_0, w_1, ..., w_d such that F(T,x) = sum w_j C(x+j, d).
    """
    d = len(F_coeffs) - 1  # degree
    n = d + 1  # number of coefficients

    # Evaluate F at x=0,1,...,d
    F_vals = []
    for x in range(n):
        val = sum(F_coeffs[k] * x**k for k in range(n))
        F_vals.append(val)

    # Build Worpitzky basis matrix
    M = worpitzky_basis_matrix(d)

    # Solve M @ w = F_vals
    w = np.linalg.solve(M, F_vals)
    return w

def eulerian_number(n, k):
    """A(n,k) = number of permutations of [n] with exactly k descents"""
    return sum((-1)**j * comb(n+1, j) * (k+1-j)**n for j in range(k+2))

random.seed(42)

print("=== WORPITZKY EXPANSION OF F(T,x) ===\n")

for n in [3, 4, 5, 6, 7]:
    print(f"\n{'='*60}")
    print(f"n = {n}")
    print(f"{'='*60}")

    # Eulerian numbers for reference
    E = [eulerian_number(n-1, k) for k in range(n-1)]
    print(f"Eulerian numbers A({n-1},k): {E}")
    print(f"Sum = {sum(E)} = {factorial(n-1)}")

    # Transitive tournament
    A = transitive_tournament(n)
    F = forward_edge_poly(A, n)
    w = expand_in_worpitzky(F)
    print(f"\nTransitive: F = {F}")
    print(f"  Worpitzky coeffs w = {[round(x) for x in w]}")

    # Random tournaments
    num_trials = 10 if n <= 5 else 5
    all_w = []
    for trial in range(num_trials):
        A = random_tournament(n)
        F = forward_edge_poly(A, n)
        w = expand_in_worpitzky(F)
        all_w.append(w)

        if trial < 3:
            H = sum(F)
            print(f"\nTrial {trial}: F = {F}, H = {H}")
            print(f"  Worpitzky coeffs w = {[round(x,2) for x in w]}")
            print(f"  All non-negative? {all(x >= -0.01 for x in w)}")
            print(f"  Sum(w) = {sum(w):.2f}")
            print(f"  Palindromic? {all(abs(w[i] - w[len(w)-1-i]) < 0.01 for i in range(len(w)))}")

    # Check properties across all trials
    all_nonneg = sum(1 for w in all_w if all(x >= -0.01 for x in w))
    all_palindromic = sum(1 for w in all_w
                          if all(abs(w[i] - w[len(w)-1-i]) < 0.01 for i in range(len(w))))
    all_integer = sum(1 for w in all_w if all(abs(x - round(x)) < 0.01 for x in w))

    print(f"\n  Summary over {num_trials} trials:")
    print(f"    All non-negative: {all_nonneg}/{num_trials}")
    print(f"    All palindromic: {all_palindromic}/{num_trials}")
    print(f"    All integer: {all_integer}/{num_trials}")

print("\n\n=== ALTERNATIVE: h*-VECTOR EXPANSION ===")
print("The h*-vector uses the basis {x^k * (1-x)^{n-1-k}} or similar.")
print("For Ehrhart theory, the h*-vector in the Worpitzky basis gives non-negative coeffs")
print("for lattice polytopes. Our F(T,x) is palindromic like the h-vector of a Gorenstein polytope.")

# Test: expand F(T,x) in the basis t^k (1+t)^{n-1-2k} (for palindromic polys)
print("\n\n=== GAMMA EXPANSION ===")
print("For palindromic polynomials, the gamma expansion is:")
print("F(T,x) = sum_{k} gamma_k * x^k * (1+x)^{n-1-2k}")
print("gamma-positivity means all gamma_k >= 0.")
print("This is a STRONGER property than palindromicity + non-negativity.")

for n in [3, 4, 5, 6, 7]:
    print(f"\nn = {n}:")

    def gamma_expand(F):
        """Expand palindromic F in gamma basis x^k(1+x)^{d-2k} where d=len(F)-1"""
        d = len(F) - 1  # degree
        m = d // 2  # max gamma index

        # Build basis matrix: evaluate x^k(1+x)^{d-2k} at x=0,...,d
        M = np.zeros((d+1, m+1))
        for j in range(d+1):  # evaluation point
            for k in range(m+1):  # gamma index
                M[j][k] = j**k * (1+j)**(d-2*k)

        # Least squares since system may be overdetermined
        gamma, _, _, _ = np.linalg.lstsq(M, [sum(F[i]*j**i for i in range(d+1)) for j in range(d+1)], rcond=None)

        # Verify
        F_reconstructed = [0] * (d+1)
        for k in range(m+1):
            # x^k(1+x)^{d-2k} expanded
            for j in range(d-2*k+1):
                if k+j <= d:
                    F_reconstructed[k+j] += gamma[k] * comb(d-2*k, j)

        err = max(abs(F[i] - F_reconstructed[i]) for i in range(d+1))
        return gamma, err

    # A few examples
    num_tests = 20 if n <= 5 else 5
    gamma_pos_count = 0
    for trial in range(num_tests):
        A = random_tournament(n)
        F = forward_edge_poly(A, n)
        gamma, err = gamma_expand(F)

        if all(g >= -0.01 for g in gamma) and err < 0.1:
            gamma_pos_count += 1

        if trial < 2:
            print(f"  Trial {trial}: F={F}, gamma={[round(g,2) for g in gamma]}, err={err:.2e}")
            print(f"    gamma-positive? {all(g >= -0.01 for g in gamma)}")

    print(f"  gamma-positive: {gamma_pos_count}/{num_tests}")
