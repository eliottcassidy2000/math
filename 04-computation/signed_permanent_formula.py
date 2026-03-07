#!/usr/bin/env python3
"""
signed_permanent_formula.py — Explore the formula for S(T)/2^{n-1} at n=5,7.

Known: S/16 = H - 3*t3 at n=5 (verified).
Question: What is the formula at n=7?

S(T) = sum_P prod (2*A[P_i][P_{i+1}] - 1)
     = sum_{k=0}^{n-1} (-1)^{n-1-k} * 2^k * D_k

At n=7: S = -D_0 + 2*D_1 - 4*D_2 + 8*D_3 - 16*D_4 + 32*D_5 - 64*D_6
       S/64 = D_6 - D_5/2 + D_4/4 - D_3/8 + D_2/16 - D_1/32 + D_0/64

D_6 = H (number of Hamiltonian paths = permutations with all edges forward)
D_5 = H*(n-1) - ... (permutations with exactly one backward edge)

Actually D_k = sum_P C(fwd(P), k). And D_{n-1} = H.

So S = sum_{k} (-1)^{n-1-k} 2^k D_k.

For n=7 (odd, n-1=6):
S = D_0 - 2*D_1 + 4*D_2 - 8*D_3 + 16*D_4 - 32*D_5 + 64*D_6
  = D_0 - 2*D_1 + 4*D_2 - 8*D_3 + 16*D_4 - 32*D_5 + 64*H

S/64 = H + (D_0 - 2*D_1 + 4*D_2 - 8*D_3 + 16*D_4 - 32*D_5)/64

For n=5:
S = D_0 - 2*D_1 + 4*D_2 - 8*D_3 + 16*D_4
  = D_0 - 2*D_1 + 4*D_2 - 8*D_3 + 16*H

S/16 = H + (D_0 - 2*D_1 + 4*D_2 - 8*D_3)/16

We know S/16 = H - 3*t3, so:
(D_0 - 2*D_1 + 4*D_2 - 8*D_3)/16 = -3*t3

D_0 = n! = 120
D_1 = n!*(n-1)/2 = 120*4/2 = 240 (always, since C(n,2) edges total, each perm
      has n-1 positions, sum of forward edges over all perms = (n-1)*n!/2 by symmetry)
Wait, D_1 = sum_P fwd(P). By symmetry each of the C(n,2) edges is in the
forward direction in exactly n!*(n-1)/(2*C(n,2))... no.

Actually D_1 = sum_P fwd(P). Each edge (i,j) with A[i][j]=1 contributes to
position p of permutation P whenever P_p = i and P_{p+1} = j. The number of
such perms is (n-2)!. So D_1 = C(n,2) * (n-2)! = C(n,2) * (n-2)!.

Wait that's wrong too. Each ORDERED pair (i,j) with A[i][j]=1 appears in
position p iff P_p=i, P_{p+1}=j. Number of perms with P_p=i, P_{p+1}=j
is (n-2)!. Number of positions is n-1. But each perm P only has one occurrence
of the consecutive pair (i,j). So:

D_1 = sum_P sum_{p=0}^{n-2} A[P_p][P_{p+1}] = sum_{(i,j): A[i][j]=1}
      #{perms with (i,j) consecutive} = sum_{(i,j): A[i][j]=1} (n-1)*(n-2)!

Wait no. The number of perms where i immediately precedes j is (n-1)!...
Hmm. Fix i at position p, j at position p+1. There are (n-2)! arrangements
of remaining n-2 elements. There are n-1 choices of p. So total = (n-1)*(n-2)!
= (n-1)!. But the SAME perm can have (i,j) at different positions? No, i
appears exactly once, j appears exactly once.

So #{perms with i immediately before j} = (n-1)!. (Fix i before j as
consecutive pair, permute the rest and the pair's position.)

Actually: think of it as inserting the pair (i,j) as a block into n-1 possible
positions among n-2 other elements: (n-1) * (n-2)! = (n-1)!. Yes.

D_1 = sum_{(i,j): A[i][j]=1} (n-1)! = C(n,2) * (n-1)!

For n=5: D_1 = 10 * 24 = 240. Let me verify.

Author: opus-2026-03-07-S43b
"""
from itertools import permutations, combinations
from collections import Counter
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

def compute_all(A, n):
    """Compute H, S, D_k, t3, t5."""
    D = [0] * n
    H = 0
    S = 0

    for P in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if A[P[i]][P[i+1]])
        for k in range(n):
            D[k] += math.comb(fwd, k)

        prod_b = 1
        for i in range(n-1):
            prod_b *= (2*A[P[i]][P[i+1]] - 1)
        S += prod_b

    H = D[n-1]  # All edges forward = Hamiltonian path

    # Count 3-cycles
    t3 = 0
    for triple in combinations(range(n), 3):
        i, j, k = triple
        if (A[i][j] and A[j][k] and A[k][i]) or (A[i][k] and A[k][j] and A[j][i]):
            t3 += 1

    # Count 5-cycles
    t5 = 0
    for quint in combinations(range(n), 5):
        for perm in permutations(quint):
            is_cycle = all(A[perm[i]][perm[(i+1)%5]] for i in range(5))
            if is_cycle:
                t5 += 1
        # Each 5-cycle counted 5 times (rotations), 2 directions = 10 times
    t5 //= 10

    return H, S, D, t3, t5

# n=5 verification
print("=== n=5: Dk values and S formula ===")
n = 5
m = n*(n-1)//2

# Collect data by isomorphism class (use (H,t3,t5) as proxy)
data5 = {}
for bits in range(1 << m):
    A = tournament_from_bits(bits, n)
    H, S, D, t3, t5 = compute_all(A, n)
    key = (H, t3, t5)
    if key not in data5:
        data5[key] = (H, S, D[:], t3, t5)

for key in sorted(data5.keys()):
    H, S, D, t3, t5 = data5[key]
    c0 = S / 2**(n-1)
    # Check alternating sum
    alt_sum = sum((-1)**(n-1-k) * (2**k) * D[k] for k in range(n))
    print(f"H={H:3d} t3={t3} t5={t5} S={S:5d} c0={c0:5.1f} D={D}")
    print(f"  alt_sum={alt_sum} (should = S={S}), D0={D[0]} D1={D[1]} D2={D[2]} D3={D[3]} D4={D[4]}")
    # Verify D_0 = n! = 120, D_1 = C(n,2)*(n-1)! = 10*24 = 240
    assert D[0] == math.factorial(n), f"D0 wrong: {D[0]}"

print(f"\nD_0 = {math.factorial(n)} (universal)")
print(f"D_1 = C(n,2)*(n-1)! = {math.comb(n,2)*math.factorial(n-1)} (universal)")

# Now n=7 sampling
print("\n=== n=7: S/64 formula exploration (sampling) ===")
n = 7
m = n*(n-1)//2
random.seed(42)

print(f"{'H':>5} {'S':>8} {'S/64':>7} {'c0':>6} {'t3':>3} {'t5':>3} {'H-3t3':>7} {'diff':>7}")

results7 = []
for trial in range(30):
    bits = random.getrandbits(m)
    A = tournament_from_bits(bits, n)
    H, S, D, t3, t5 = compute_all(A, n)
    c0 = S / 64
    diff = c0 - (H - 3*t3)
    results7.append((H, S, c0, t3, t5, D[:]))
    print(f"{H:5d} {S:8d} {c0:7.2f} {c0:6.2f} {t3:3d} {t5:3d} {H-3*t3:7d} {diff:7.2f}")

# Check what S/64 depends on
print("\n=== Regression: S/64 vs H, t3, t5 ===")
print("Trying: S/64 = a*H + b*t3 + c*t5 + d")
# Use first 20 for fitting, last 10 for testing
import numpy as np

try:
    X = []
    Y = []
    for H, S, c0, t3, t5, D in results7:
        X.append([H, t3, t5, 1])
        Y.append(c0)
    X = np.array(X, dtype=float)
    Y = np.array(Y, dtype=float)

    # Least squares
    coeffs, residuals, rank, sv = np.linalg.lstsq(X, Y, rcond=None)
    print(f"  Coefficients: H*{coeffs[0]:.6f} + t3*{coeffs[1]:.6f} + t5*{coeffs[2]:.6f} + {coeffs[3]:.6f}")
    pred = X @ coeffs
    max_err = max(abs(Y - pred))
    print(f"  Max error: {max_err:.6f}")

    # Try H, t3, t5, t7, bc33
    # Count t7
    print("\n=== With D_k values ===")
    # Check: S = sum (-1)^{n-1-k} * 2^k * D_k
    for H, S, c0, t3, t5, D in results7[:5]:
        alt = sum((-1)**(n-1-k) * (2**k) * D[k] for k in range(n))
        print(f"  H={H}, S={S}, alt_sum={alt}, match={S==alt}")
        # D values
        print(f"    D = {D}")
        # c0 = S/64 = H + (lower)/64 where lower = S - 64*H
        lower = S - 64*H
        print(f"    S - 64H = {lower}, (S-64H)/64 = {lower/64:.4f}")

except ImportError:
    print("  numpy not available, skipping regression")
