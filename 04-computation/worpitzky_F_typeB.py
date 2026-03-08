#!/usr/bin/env python3
"""
F(T,x) IN TYPE B WORPITZKY BASIS

The standard Worpitzky identity: x^n = sum_k A(n,k) C(x+k,n)
The type B Worpitzky identity: (2x+1)^{n} = sum_k T_B(n+1,k+1) C(x+k,n)

F(T,x) is palindromic: F_k = F_{n-1-k}. It lives in the space of
polynomials of degree n-1 that satisfy p(x) = x^{n-1} p(1/x).

What if we expand F(T,x) in the TYPE B Worpitzky basis instead of
the standard Worpitzky basis?

The type B basis is {C(x+k, n-1)} for k=0,...,n-1 (same basis!).
The difference is in what the coefficients represent.

For the standard Worpitzky expansion of F(T,x):
  F(T,x) = sum_k w_k C(x+k, n-1)
We found w_k are integers but NOT palindromic, NOT non-negative.

KEY INSIGHT: Perhaps the right variable substitution is x = (t-1)/2
or x = (2t+1)/(2t-1) (Mobius), giving a different basis expansion.

Actually, the Mobius transform THM-K says:
  W(T,r) = (r-1/2)^{n-1} * F(T, (2r+1)/(2r-1))

At r = (x-1)/2: 2r+1 = x, 2r-1 = x-2, (2r+1)/(2r-1) = x/(x-2)
  W(T, (x-1)/2) = ((x-1)/2 - 1/2)^{n-1} * F(T, x/(x-2))
                 = ((x-2)/2)^{n-1} * F(T, x/(x-2))

Hmm, let me try a different approach. The transitive tournament has
F(T,x) = x^{n-1} and this maps to A(n-1,k) in Worpitzky. For a general
tournament, F(T,x) = sum F_k x^k is palindromic.

Since F is palindromic, it can be written in terms of the symmetric
function of x and 1/x... or more usefully, in the gamma basis.

But here's another angle: the TYPE B EULERIAN POLYNOMIAL is defined as
B_n(x) = sum_k T_B(n,k+1) x^k. It satisfies:
  B_n(x) * (1-x)^{-n} = sum_{m>=0} (2m+1)^{n-1} x^m

This is the type B analog of
  A_n(x) * (1-x)^{-n-1} = sum_{m>=0} (m+1)^n x^m

KEY QUESTION: Does the F-polynomial of a tournament relate to the
type B Eulerian polynomial through some substitution?

F(T,x) at x=1: F(T,1) = H(T)
B_n(1) = (2n)!! / 2^n = n! * 2^n / 2^n = n!
Wait, row sums of T_B(n,k) are (2n)!!/n! actually... let me check.
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

def all_tournaments(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

def forward_edge_poly(A, n):
    F = [0] * n
    for P in permutations(range(n)):
        if not all(A[P[i]][P[i+1]] for i in range(n-1)):
            continue
        fwd = sum(1 for i in range(n-1) if P[i] < P[i+1])
        F[fwd] += 1
    return F

def typeB_eulerian(m, k):
    """T_B(m, k) for k=1,...,m"""
    return sum((-1)**(k-i) * comb(m, k-i) * (2*i-1)**(m-1) for i in range(1, k+1))

random.seed(42)

# ===== Understand F(T,x) structure =====
print("=" * 70)
print("F(T,x) PALINDROMIC STRUCTURE AND WORPITZKY")
print("=" * 70)

# For palindromic polynomials, a natural basis is:
# q_k(x) = x^k + x^{d-k} for k = 0, ..., floor(d/2)
# and x^{d/2} if d is even.
# But this doesn't connect to Worpitzky/Eulerian.

# ALTERNATIVE: The "h-vector" interpretation.
# If F(T,x) = h_0 + h_1 x + ... + h_d x^d is palindromic with h_k = h_{d-k},
# then f(x) = (1-x)^{d+1} * sum_{m>=0} F(T,m+1) x^m
# where F(T,m+1) = sum_k h_k * C(m+k, d)
# This is the Ehrhart-Macdonald-Stanley formula.
# The h-vector is nonneg iff it's the h-vector of a simplicial polytope.

# For the STANDARD Worpitzky expansion x^d = sum A(d,k) C(x+k,d):
# This is the h-vector when counting lattice points.

# For OUR F(T,x):
# F(T, x) = sum_k F_k x^k = sum_j w_j C(x+j, d)
# Then evaluating at x=m (nonneg integer):
# F(T, m) = sum_j w_j C(m+j, d)
# This counts... lattice points of some polytope?

# Let me check what F(T, m) counts for small m.
print("\n=== F(T, m) for integer m ===")

for n in [4, 5]:
    d = n - 1
    print(f"\nn={n}:")
    for trial in range(3):
        A = random_tournament(n)
        F = forward_edge_poly(A, n)
        H = sum(F)
        print(f"  Trial {trial}: F={F}, H={H}")
        for m in range(6):
            Fm = sum(F[k] * m**k for k in range(n))
            print(f"    F(T,{m}) = {Fm}")
        # F(T,0) = F_0, F(T,1) = H
        # F(T,-1) = sum F_k (-1)^k = ... for palindromic, F(T,-1) = (-1)^{d} F(T,-1)
        # so for even d: F(T,-1) = F(T,-1) (trivial)
        # for odd d: F(T,-1) = -F(T,-1), so F(T,-1) = 0.
        Fm1 = sum(F[k] * (-1)**k for k in range(n))
        print(f"    F(T,-1) = {Fm1} (should be 0 for odd d={d})")

# ===== The TYPE B connection =====
print("\n\n" + "=" * 70)
print("TYPE B EULERIAN POLYNOMIAL CONNECTION")
print("=" * 70)

# The type B Eulerian polynomial:
# B_n(x) = sum_{k=0}^{n-1} T_B(n, k+1) x^k
# These are palindromic! T_B(n,k) = T_B(n, n+1-k).
# B_n(1) = row sum = 2^n * n! / (wait, let me compute)

print("\nType B Eulerian polynomials:")
for m in range(2, 7):
    coeffs = [typeB_eulerian(m, k+1) for k in range(m)]
    total = sum(coeffs)
    print(f"  B_{m}(x) coeffs = {coeffs}, sum = {total}")
    print(f"    2^{m-1} * (m-1)! = {2**(m-1) * factorial(m-1)}")
    print(f"    (2m-1)!! = {int(np.prod(range(1, 2*m, 2)))}")

# So B_n(1) = (2n-1)!! (double factorial of odd numbers)

# Now: F(T,x) / H gives a "normalized descent polynomial".
# For the transitive tournament (H=1), F = [0,...,0,1], i.e., x^{n-1}.
# For the "anti-transitive" (reverse labeling) with H=1, F = [1,0,...,0], i.e., 1.
# For a random tournament, F is somewhere in between.

# KEY QUESTION: Is F(T,x) related to a "partial type B Eulerian polynomial"?
# I.e., can we decompose B_n(x) into tournament-indexed pieces?

# The standard Eulerian polynomial A_n(x) sums over ALL permutations.
# F(T,x) sums over Hamiltonian paths of T (a subset of permutations).
# Sum over ALL tournaments T of F(T,x) = sum over all perms of x^{asc(perm)}
# since each perm is a HP of exactly one tournament (determined by the edge directions).

# Actually NO: each perm is a valid HP when A[P_i][P_{i+1}]=1 for all i.
# For a fixed permutation P, it's a HP of T iff all edges P_i->P_{i+1} are in T.
# The number of tournaments for which P is a HP: we need A[P_i][P_{i+1}]=1 for
# the n-1 edges on the path, and the remaining C(n,2)-(n-1) edges can be anything.
# So it's 2^{C(n,2)-(n-1)}.

print("\n\n" + "=" * 70)
print("SUM OVER ALL TOURNAMENTS OF F(T,x)")
print("=" * 70)

for n in [3, 4, 5]:
    d = n - 1
    total_F = [0] * n
    count = 0
    for A in all_tournaments(n):
        F = forward_edge_poly(A, n)
        for k in range(n):
            total_F[k] += F[k]
        count += 1

    print(f"\nn={n}: sum F(T,x) over all {count} tournaments:")
    print(f"  Total F_k = {total_F}")
    print(f"  Sum = {sum(total_F)} = {count} * n! / ? ")
    # Each permutation P is a HP of 2^{C(n,2)-(n-1)} tournaments
    # Total HPs = n! * 2^{C(n,2)-(n-1)}
    extra = comb(n,2) - (n-1)
    expected_total = factorial(n) * 2**extra
    print(f"  Expected total = {n}! * 2^{extra} = {expected_total}")

    # The distribution of total_F_k: since each perm contributes to F_{asc(perm)},
    # total_F_k = #{perms with k ascents} * 2^{C(n,2)-(n-1)}
    # = A(n, k) * 2^{extra}
    E = [sum((-1)**j * comb(n+1, j) * (k+1-j)**n for j in range(k+2)) for k in range(n)]
    expected_F = [E[k] * 2**extra for k in range(n)]
    print(f"  Eulerian * 2^{extra} = {expected_F}")
    print(f"  Match: {total_F == expected_F}")

# ===== F(T,x) as restricted descent polynomial =====
print("\n\n" + "=" * 70)
print("F(T,x) = RESTRICTED EULERIAN POLYNOMIAL")
print("=" * 70)
print("F(T,x) counts permutations (HPs of T) by ascent number.")
print("It's the 'descent polynomial' of the sub-poset of permutations")
print("compatible with tournament T.")
print()
print("KEY INSIGHT: Sum_T F(T,x) = A_n(x) * 2^{C(n,2)-(n-1)}")
print("This means the AVERAGE F(T,x) over all tournaments is")
print("proportional to the Eulerian polynomial A_n(x).")
print()

for n in [4, 5]:
    d = n - 1
    extra = comb(n,2) - (n-1)
    print(f"\nn={n}:")
    # Average F(T,x) = A_n(x) * n! / total_tournaments... actually:
    # Average F(T,x) = A_n(x) * 2^{extra} / 2^{C(n,2)}
    #                = A_n(x) / 2^{n-1}
    # since extra = C(n,2)-(n-1) = (n(n-1)/2)-(n-1) = (n-1)(n-2)/2
    # and C(n,2) = n(n-1)/2
    # So 2^extra / 2^{C(n,2)} = 1 / 2^{n-1}
    E = [sum((-1)**j * comb(n+1, j) * (k+1-j)**n for j in range(k+2)) for k in range(n)]
    avg_F = [e / 2**(n-1) for e in E]
    print(f"  Average F(T,x) coeffs = {avg_F}")
    print(f"  = A_{n}(x) / 2^{n-1}")
    print(f"  Average H = {sum(avg_F)} = {factorial(n)} / 2^{n-1}")

# ===== Does F(T,x) / H have a universal Worpitzky expansion? =====
print("\n\n" + "=" * 70)
print("F(T,x) / H IN WORPITZKY BASIS")
print("=" * 70)
print("We know from the wrong-W calculation that this is NOT universal.")
print("But maybe the SECOND Worpitzky coefficient (w_1) has a universal part?")

from collections import defaultdict

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

for n in [5]:
    d = n - 1
    print(f"\nn={n}:")
    by_H = defaultdict(list)

    for A in all_tournaments(n):
        F = forward_edge_poly(A, n)
        w = expand_in_worpitzky(F)
        H = sum(F)
        by_H[H].append(w)

    for H in sorted(by_H.keys()):
        ws = by_H[H]
        # Check which components are universal (same for all T with same H)
        for k in range(n):
            vals = set(w[k] for w in ws)
            if len(vals) == 1:
                print(f"  H={H}, w_{k} = {vals.pop()} (UNIVERSAL)")
        # Check w_{n-2}
        vals_n2 = set(w[n-2] for w in ws)
        if len(vals_n2) == 1:
            print(f"  H={H}, w_{{n-2}} = {vals_n2.pop()} (UNIVERSAL)")

    # The last coefficient w_{n-1} = F_0 and w_0 = ? Let me check.
    print(f"\n  Checking if w_{{n-2}} = (n-1)*H:")
    for H in sorted(by_H.keys()):
        ws = by_H[H]
        vals = set(w[n-2] for w in ws)
        if len(vals) == 1:
            v = vals.pop()
            print(f"    H={H}: w_{{n-2}} = {v}, (n-1)*H = {(n-1)*H}, H = {H}")
