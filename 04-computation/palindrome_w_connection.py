#!/usr/bin/env python3
"""
HYPOTHESIS G: The palindromic forward-edge polynomial F(T,x) and the
W-polynomial W(T,r) are related by a CHANGE OF BASIS.

F(T,x) = sum_{k=0}^{n-1} D_k x^k,  palindromic: D_k = D_{n-1-k}
W(T,r) = sum_P prod_{i=0}^{n-2} (r + s_{P_i,P_{i+1}}) where s = ±1/2

Both are polynomials of degree n-1 computed from the same data.
What is their exact relationship?

The key: F(T,x) counts by NUMBER of forward edges.
W(T,r) evaluates a PRODUCT with weights.

For a permutation P with fwd(P) = k and bwd(P) = n-1-k:
  prod(r + s_e) where s_e = +1/2 if forward, -1/2 if backward
  = (r+1/2)^k * (r-1/2)^{n-1-k}

So: W(T,r) = sum_P (r+1/2)^{fwd(P)} * (r-1/2)^{n-1-fwd(P)}
           = sum_{k=0}^{n-1} D_k * (r+1/2)^k * (r-1/2)^{n-1-k}

This IS the change of basis! W is F evaluated at... well, it's more like:
  W(T,r) = sum_k D_k * (r+1/2)^k * (r-1/2)^{n-1-k}

Let y = (r+1/2)/(r-1/2). Then:
  W(T,r) = (r-1/2)^{n-1} * sum_k D_k * y^k
          = (r-1/2)^{n-1} * F(T, y)

where y = (r+1/2)/(r-1/2) = (2r+1)/(2r-1).

So: W(T,r) = (r-1/2)^{n-1} * F(T, (2r+1)/(2r-1))

Let's verify this numerically!

And at r=0: W(T,0) = (-1/2)^{n-1} * F(T, -1) = (-1/2)^{n-1} * S(T)
where S(T) = sum_P (-1)^{bwd(P)} = sum_P prod(-1)^{(1-A)/1}...

Wait, F(T,-1) = sum_P (-1)^{fwd(P)}, not the signed permanent.
S(T) = sum_P prod(2A[P_i][P_{i+1}]-1) = sum_P prod(1 if fwd, -1 if bwd)
     = sum_P (-1)^{bwd(P)} = sum_P (-1)^{n-1-fwd(P)}
     = (-1)^{n-1} sum_P (-1)^{fwd(P)} * (-1)^{n-1}... no.

Actually: S(T) = sum_P prod_{i=0}^{n-2} (2A[P_i][P_{i+1}]-1)
If A=1 (forward): factor = +1
If A=0 (backward): factor = -1

So prod = (+1)^{fwd} * (-1)^{bwd} = (-1)^{bwd} = (-1)^{n-1-fwd}

S(T) = sum_P (-1)^{n-1-fwd(P)} = (-1)^{n-1} * sum_P (-1)^{fwd(P)}
     = (-1)^{n-1} * F(T, -1)

And W(T,0) = (-1/2)^{n-1} * F(T, (0+1/2)/(0-1/2))
           = (-1/2)^{n-1} * F(T, -1)
           = (1/2)^{n-1} * (-1)^{n-1} * F(T, -1)
           = (1/2)^{n-1} * S(T)

So c_0 = W(T,0) = S(T)/2^{n-1}.

And H(T) = D_{n-1} = coefficient of x^{n-1} in F(T,x).
Also H(T) = F(T,1) - ... no, F(T,1) = n!.

H(T) = D_{n-1}. From the W relation:
W(T,r) = sum_k D_k * (r+1/2)^k * (r-1/2)^{n-1-k}

At r=1/2: W(T,1/2) = sum_k D_k * 1^k * 0^{n-1-k} = D_{n-1} = H(T)
At r=-1/2: W(T,-1/2) = sum_k D_k * 0^k * (-1)^{n-1-k} = D_0 = H(T^op)

By palindromicity, D_0 = D_{n-1}, so W(1/2) = W(-1/2) = H(T).

Let's verify all this.
"""
from itertools import permutations
import numpy as np
from math import factorial

def tiling_to_tournament(bits, n):
    A = [[0]*n for _ in range(n)]
    for i in range(1, n):
        A[i][i-1] = 1
    tiles = [(a,b) for a in range(n) for b in range(a) if a-b >= 2]
    tiles.sort()
    for idx, (a, b) in enumerate(tiles):
        if (bits >> idx) & 1: A[b][a] = 1
        else: A[a][b] = 1
    return A

def random_tournament(n):
    import random
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

def compute_F(A, n):
    """Forward-edge polynomial F(T,x) = sum_P x^fwd(P)"""
    D = [0] * n
    for P in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if A[P[i]][P[i+1]])
        D[fwd] += 1
    return D

def compute_W_eval(A, n, r):
    """W(T,r) = sum_P prod(r + s_e)"""
    W = 0
    for P in permutations(range(n)):
        prod_val = 1
        for i in range(n-1):
            s = 0.5 if A[P[i]][P[i+1]] else -0.5
            prod_val *= (r + s)
        W += prod_val
    return W

def eval_F(D, x):
    return sum(D[k] * x**k for k in range(len(D)))

# ========== Verification ==========
print("=== Verify W(T,r) = (r-1/2)^{n-1} * F(T, (2r+1)/(2r-1)) ===\n")

import random
random.seed(42)

for n in [4, 5, 6]:
    print(f"n={n}:")
    for trial in range(3):
        A = random_tournament(n)
        D = compute_F(A, n)

        # Check palindromicity
        palindrome = all(D[k] == D[n-1-k] for k in range(n))

        # Check W(r) = (r-1/2)^{n-1} * F(y) for several r values
        max_err = 0
        for r in [0.0, 0.3, 0.7, 1.0, 1.5, 2.0, -0.3, -1.0]:
            W_actual = compute_W_eval(A, n, r)
            if abs(r - 0.5) < 1e-10:
                continue  # y -> infinity
            y = (2*r + 1) / (2*r - 1)
            W_predicted = (r - 0.5)**(n-1) * eval_F(D, y)
            err = abs(W_actual - W_predicted)
            max_err = max(max_err, err)

        print(f"  trial {trial}: palindrome={palindrome}, max_err={max_err:.2e}")

    # Special evaluations
    A = random_tournament(n)
    D = compute_F(A, n)
    H = D[n-1]
    S = compute_W_eval(A, n, 0) * 2**(n-1)  # c_0 * 2^{n-1}

    print(f"  H = D_{n-1} = {H}")
    print(f"  D_0 = {D[0]} (should equal H by palindromicity)")
    print(f"  W(1/2) = {compute_W_eval(A, n, 0.5):.1f} (should equal H = {H})")
    print(f"  W(-1/2) = {compute_W_eval(A, n, -0.5):.1f} (should equal H = {H})")
    print(f"  W(0) = {compute_W_eval(A, n, 0):.4f}, S/2^(n-1) = {S/2**(n-1):.4f}")
    print(f"  F(-1) = {eval_F(D, -1)}, (-1)^(n-1)*S = {(-1)**(n-1) * S:.0f}")
    print()

# ========== DEEP: What does perpendicularity mean in F-polynomial terms? ==========
print("\n=== Perpendicularity in F-polynomial terms ===")
n = 5
tiles = [(a,b) for a in range(n) for b in range(a) if a-b >= 2]
m = len(tiles)

# For each tiling, compute F-polynomial coefficients D_0,...,D_4
all_D = []
for bits in range(1 << m):
    A = tiling_to_tournament(bits, n)
    D = compute_F(A, n)
    all_D.append(D)

all_D = np.array(all_D, dtype=float)

# By palindromicity: D_0=D_4, D_1=D_3, and D_0+D_1+D_2+D_3+D_4 = 120
# So we have 2 free variables: D_0 (=H) and D_1 (or equivalently D_2 = 120-2D_0-2D_1)
# The "W-polynomial layers" are:
#   c_0 = W(0) = (-1/2)^4 * F(-1) = (1/16) * (D_0 - D_1 + D_2 - D_3 + D_4)
#       = (1/16) * (2D_0 - 2D_1 + D_2)
#       = (1/16) * (2D_0 - 2D_1 + 120 - 2D_0 - 2D_1)
#       = (1/16) * (120 - 4D_1) = (30 - D_1)/4
#
#   H = D_0 = D_4

print("Correlation matrix of D_k:")
corr = np.corrcoef(all_D.T)
for i in range(n):
    print(f"  D_{i}: ", end="")
    for j in range(n):
        print(f"{corr[i,j]:>7.3f}", end=" ")
    print()

print(f"\n  c_0 = (30 - D_1)/4")
print(f"  H = D_0 = D_4")
print(f"  Corr(D_0, D_1) = {corr[0,1]:.4f}")
print(f"  Corr(D_0, D_2) = {corr[0,2]:.4f}")
print(f"  Corr(D_1, D_2) = {corr[1,2]:.4f}")

# c_0 and H are perpendicular means c_0 and D_0 are nearly uncorrelated
# c_0 = (30 - D_1)/4, so Cov(c_0, H) = Cov((30-D_1)/4, D_0) = -Cov(D_1, D_0)/4
# Perpendicularity means Corr(D_0, D_1) ≈ 0!

c0_arr = (30 - all_D[:,1]) / 4
H_arr = all_D[:,0]
print(f"\n  Corr(c_0, H) = {np.corrcoef(c0_arr, H_arr)[0,1]:.4f}")
print(f"  This equals Corr(D_0, D_1) with sign flip: "
      f"{-np.corrcoef(all_D[:,0], all_D[:,1])[0,1]:.4f}")

# So PERPENDICULARITY AT n=5 means: D_0 and D_1 are nearly uncorrelated!
# i.e., the number of Hamiltonian paths (D_0=D_4) is nearly independent of
# the number of paths with exactly 1 forward edge (D_1=D_3).

print("\n\n=== D_0 vs D_1 scatter ===")
# Show the joint distribution
from collections import Counter
pairs = Counter(zip(all_D[:,0].astype(int), all_D[:,1].astype(int)))
for (d0, d1), count in sorted(pairs.items()):
    print(f"  D_0={d0:>2}, D_1={d1:>2}: count={count:>2}")

# ========== At n=7 (sampled) ==========
print("\n\n=== Perpendicularity at n=7 ===")
n = 7
random.seed(42)

D0_list = []
D1_list = []
D2_list = []
c0_list = []

for trial in range(200):
    A = random_tournament(n)
    D = compute_F(A, n)
    D0_list.append(D[0])  # = H
    D1_list.append(D[1])  # = D_5
    D2_list.append(D[2])  # = D_4
    # c_0 = W(0) = (-1/2)^6 * F(-1) = (1/64) * F(-1)
    F_neg1 = eval_F(D, -1)
    c0_list.append(F_neg1 / 64)

D0_arr = np.array(D0_list)
D1_arr = np.array(D1_list)
D2_arr = np.array(D2_list)
c0_arr = np.array(c0_list)

print(f"  Corr(D_0, D_1) = {np.corrcoef(D0_arr, D1_arr)[0,1]:.4f}")
print(f"  Corr(D_0, D_2) = {np.corrcoef(D0_arr, D2_arr)[0,1]:.4f}")
print(f"  Corr(D_1, D_2) = {np.corrcoef(D1_arr, D2_arr)[0,1]:.4f}")
print(f"  Corr(c_0, H=D_0) = {np.corrcoef(c0_arr, D0_arr)[0,1]:.4f}")
print(f"  Corr(c_0, D_1) = {np.corrcoef(c0_arr, D1_arr)[0,1]:.4f}")
print(f"  Corr(c_0, D_2) = {np.corrcoef(c0_arr, D2_arr)[0,1]:.4f}")

# ========== HYPOTHESIS H: H and c_0 are weakly correlated because they live in
# different "spectral bands" of the tournament ===
print(f"\n  Key insight: perpendicularity = near-independence of D_0 and D_1")
print(f"  D_0 = HP count (global connectivity)")
print(f"  D_1 = paths with exactly 1 forward edge (local orientation)")
print(f"  These measure DIFFERENT aspects of tournament structure!")
