#!/usr/bin/env python3
"""
f_poly_catalan_connection.py — Look for hidden combinatorial structures
in F(T,x) coefficients.

IDEA 1: The RECIPROCAL POLYNOMIAL. Since F is palindromic of degree n-1,
define R(T,x) = F(T,x) / (1+x)^? or F(T,x) / H.
The "gamma expansion" of palindromic polynomials:
  F(x) = sum gamma_i * x^i * (1+x)^{n-1-2i}
where the gamma_i are the "gamma coefficients".
If all gamma_i >= 0, the polynomial is "gamma-positive" (implies unimodality).

IDEA 2: Compare F with the EULERIAN POLYNOMIAL of the symmetric group:
  A_n(x) = sum_{sigma in S_n} x^{des(sigma)} = sum A(n,k) x^k
where des = number of descents. For the "trivial tournament" (total order),
F(T,x) should equal A_n(x) up to normalization.

IDEA 3: Is there a q-analog? Replace x by q and look at F as a q-analog of n!.

IDEA 4: RECIPROCITY: F(T, 1/x) * x^{n-1} = F(T, x) (palindrome).
What about F(T, -x)?

Author: opus-2026-03-07-S44
"""
from itertools import permutations, combinations
import math
import random

def tournament_from_bits(bits, n):
    A = [[0]*n for _ in range(n)]
    pos = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> pos) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            pos += 1
    return A

def compute_F(A, n):
    F = [0] * n
    for P in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if A[P[i]][P[i+1]])
        F[fwd] += 1
    return F

def gamma_expansion(F):
    """Compute gamma coefficients of palindromic polynomial F.
    F(x) = sum gamma_i * x^i * (1+x)^{d-2i} where d = deg(F).
    """
    d = len(F) - 1  # degree
    gamma = [0] * (d // 2 + 1)

    # Work from lowest gamma upward
    # gamma_0 * (1+x)^d + gamma_1 * x * (1+x)^{d-2} + ...
    remaining = list(F)

    for i in range(d // 2 + 1):
        # gamma_i contributes to coefficients i through d-i
        # via x^i * (1+x)^{d-2i}
        # The coefficient of x^i in the expansion gives gamma_i = remaining[i]
        # after subtracting previous gamma contributions.
        gamma[i] = remaining[i]

        # Subtract gamma_i * x^i * (1+x)^{d-2i} from remaining
        deg_binom = d - 2*i
        for j in range(deg_binom + 1):
            remaining[i + j] -= gamma[i] * math.comb(deg_binom, j)

    return gamma

# ============================================================
# GAMMA-POSITIVITY
# ============================================================
print("=" * 60)
print("GAMMA EXPANSION OF F(T,x)")
print("=" * 60)

for n in [3, 5]:
    m = n*(n-1)//2
    seen_F = set()

    for bits in range(1 << m):
        A = tournament_from_bits(bits, n)
        F = compute_F(A, n)
        key = tuple(F)
        if key in seen_F:
            continue
        seen_F.add(key)

        H = F[n-1]
        gamma = gamma_expansion(F)
        all_nonneg = all(g >= 0 for g in gamma)
        print(f"  n={n} H={H:3d} F={F} gamma={gamma} {'POSITIVE' if all_nonneg else 'NEGATIVE'}")

# n=7 sampling
print(f"\n  n=7 sampling:")
n = 7
m = n*(n-1)//2
random.seed(42)
pos_count = 0
neg_count = 0
neg_examples = []
seen_F = set()

for trial in range(300):
    bits = random.getrandbits(m)
    A = tournament_from_bits(bits, n)
    F = compute_F(A, n)
    key = tuple(F)
    if key in seen_F:
        continue
    seen_F.add(key)

    gamma = gamma_expansion(F)
    if all(g >= 0 for g in gamma):
        pos_count += 1
    else:
        neg_count += 1
        neg_examples.append((F[n-1], gamma))

total = pos_count + neg_count
print(f"  {pos_count}/{total} gamma-positive, {neg_count}/{total} have negative gamma")
if neg_examples:
    for H, g in neg_examples[:3]:
        print(f"    H={H}: gamma={g}")

# ============================================================
# COMPARISON WITH EULERIAN NUMBERS
# ============================================================
print("\n" + "=" * 60)
print("COMPARISON WITH EULERIAN NUMBERS A(n,k)")
print("=" * 60)

# The Eulerian polynomial A_n(x) = sum A(n,k) x^k
# where A(n,k) = # permutations of [n] with exactly k descents.
# A descent at position i means sigma(i) > sigma(i+1).

# For the TRANSITIVE tournament (total order 0 < 1 < ... < n-1):
# A "descent" w.r.t. the tournament is: A[P_i][P_{i+1}] = 0, i.e., P_i > P_{i+1}
# So fwd(P) = (n-1) - des(P) where des counts "tournament descents" = standard descents.
# F_k = #{P : fwd = k} = #{P : des = n-1-k} = A(n, n-1-k)

# Let's verify:
def eulerian_numbers(n):
    """Compute Eulerian numbers A(n,k) for k=0..n-1."""
    A = [0] * n
    A[0] = 1
    for i in range(1, n):
        new_A = [0] * n
        for k in range(n):
            new_A[k] = (k+1) * A[k] + (i - k) * (A[k-1] if k > 0 else 0)
        A = new_A
    return A

for n in [3, 5, 7]:
    eul = eulerian_numbers(n)
    eul_rev = list(reversed(eul))  # A(n, n-1-k) for k=0..n-1

    # Transitive tournament
    A_trans = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            A_trans[i][j] = 1
    F_trans = compute_F(A_trans, n)

    match = F_trans == eul_rev
    print(f"  n={n}: Eulerian A(n,k) reversed = {eul_rev}")
    print(f"  n={n}: F(transitive) = {F_trans}")
    print(f"  n={n}: Match = {match}")
    print()

# ============================================================
# F(T,x) / EULERIAN RATIO
# ============================================================
print("=" * 60)
print("F(T,x) / A_n(x) RATIO — is it structured?")
print("=" * 60)

n = 5
eul = eulerian_numbers(n)
eul_rev = list(reversed(eul))

m = n*(n-1)//2
seen_F = set()

for bits in range(1 << m):
    A = tournament_from_bits(bits, n)
    F = compute_F(A, n)
    key = tuple(F)
    if key in seen_F:
        continue
    seen_F.add(key)

    ratio = [F[k] / eul_rev[k] if eul_rev[k] != 0 else None for k in range(n)]
    H = F[n-1]
    print(f"  H={H:3d} F/A_rev = {[f'{r:.3f}' if r else 'N/A' for r in ratio]}")

# ============================================================
# INTERLACING WITH CATALAN / NARAYANA
# ============================================================
print("\n" + "=" * 60)
print("NARAYANA NUMBER CONNECTION")
print("=" * 60)
# Narayana N(n,k) = C(n,k) * C(n,k-1) / n
# The Narayana polynomial N_n(x) = sum N(n,k) x^k is also palindromic!
# Does F(T,x) "interlace" with N_n(x) for some n?

for n in [5, 7]:
    # Narayana polynomial of degree n-1: N(n-1, k) for k=0..n-2
    # Actually Narayana N(m,k) for m vertices: sum_{k=1}^{m} N(m,k) x^{k-1}
    nar = [math.comb(n-1,k) * math.comb(n-1,k+1) // (n-1) for k in range(n-1)]
    # Pad to n terms
    nar_full = nar + [0] * (n - len(nar))

    eul = eulerian_numbers(n)
    eul_rev = list(reversed(eul))

    print(f"  n={n}:")
    print(f"    Eulerian (reversed): {eul_rev}")
    print(f"    Narayana(n-1):      {nar}")
    print(f"    Sum Narayana = Catalan(n-1) = {sum(nar)}")
    print(f"    Sum Eulerian = n! = {sum(eul_rev)}")
