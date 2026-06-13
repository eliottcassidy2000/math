"""
Combined multi-degree Walsh energy as β₁ detector at n=6.
opus-2026-03-15-S72b
"""
import numpy as np
from itertools import combinations
import sys
sys.path.insert(0, '04-computation')
from path_homology_v2 import path_betti_numbers

def tournament_from_bits(n, bits):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def count_hamiltonian_paths(A):
    n = A.shape[0]
    dp = np.zeros((1 << n, n), dtype=np.int64)
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return int(sum(dp[full]))

n = 6
num_edges = n*(n-1)//2
N = 1 << num_edges

print(f"Computing H and β₁ for all {N} tournaments at n={n}...")

all_H = []
all_beta1 = []

for bits in range(N):
    A = tournament_from_bits(n, bits)
    H = count_hamiltonian_paths(A)
    all_H.append(H)
    betti = path_betti_numbers(A, n, max_dim=2)
    all_beta1.append(betti[1] if len(betti) > 1 else 0)
    if (bits+1) % 10000 == 0:
        print(f"  ...{bits+1}/{N}")

all_bits = list(range(N))
all_H_arr = np.array(all_H, dtype=float)
beta1_arr = np.array(all_beta1)

print(f"\nComputing Walsh degree-2 and degree-4 energies...")

# Compute H_hat[S] for subsets of size 2 and 4
edge_indices = list(range(num_edges))
E2 = np.zeros(N)
E4 = np.zeros(N)

# Precompute all bit arrays for sign computation
all_bits_arr = np.array(all_bits)

print("  Degree 2...")
for S in combinations(edge_indices, 2):
    mask = (1 << S[0]) | (1 << S[1])
    # Vectorized popcount mod 2
    masked = all_bits_arr & mask
    # popcount of 2-bit mask: count bits
    bit0 = (masked >> S[0]) & 1
    bit1 = (masked >> S[1]) & 1
    parity = bit0 ^ bit1
    signs = 1 - 2*parity  # (-1)^parity
    H_hat_S = np.dot(all_H_arr, signs) / N
    if abs(H_hat_S) > 1e-10:
        E2 += H_hat_S * signs

print("  Degree 4...")
cnt = 0
total = len(list(combinations(edge_indices, 4)))
for S in combinations(edge_indices, 4):
    mask = 0
    for s in S:
        mask |= (1 << s)
    masked = all_bits_arr & mask
    # popcount mod 2 for 4 bits
    parity = np.zeros(N, dtype=int)
    for s in S:
        parity ^= (masked >> s) & 1
    signs = 1 - 2*parity
    H_hat_S = np.dot(all_H_arr, signs) / N
    if abs(H_hat_S) > 1e-10:
        E4 += H_hat_S * signs
    cnt += 1
    if cnt % 500 == 0:
        print(f"    {cnt}/{total}")

print("Done.")

idx0 = beta1_arr == 0
idx1 = beta1_arr > 0

print(f"\nn=6: β₁=0 count={idx0.sum()}, β₁=1 count={idx1.sum()}")
print(f"\nE₂ by β₁:")
print(f"  β₁=0: mean={E2[idx0].mean():.4f}, range=[{E2[idx0].min():.2f}, {E2[idx0].max():.2f}]")
print(f"  β₁=1: mean={E2[idx1].mean():.4f}, range=[{E2[idx1].min():.2f}, {E2[idx1].max():.2f}]")

print(f"\nE₄ by β₁:")
print(f"  β₁=0: mean={E4[idx0].mean():.4f}, range=[{E4[idx0].min():.2f}, {E4[idx0].max():.2f}]")
print(f"  β₁=1: mean={E4[idx1].mean():.4f}, range=[{E4[idx1].min():.2f}, {E4[idx1].max():.2f}]")

# Try logistic regression
try:
    from sklearn.linear_model import LogisticRegression
    X = np.column_stack([E2, E4])
    y = beta1_arr
    clf = LogisticRegression(max_iter=1000)
    clf.fit(X, y)
    pred = clf.predict(X)
    acc = (pred == y).mean()
    print(f"\nLogistic regression (E₂, E₄) → β₁:")
    print(f"  Accuracy: {acc:.4f}")
    print(f"  Coefficients: E₂={clf.coef_[0][0]:.4f}, E₄={clf.coef_[0][1]:.4f}")

    # With H too
    X2 = np.column_stack([E2, E4, all_H_arr])
    clf2 = LogisticRegression(max_iter=1000)
    clf2.fit(X2, y)
    pred2 = clf2.predict(X2)
    acc2 = (pred2 == y).mean()
    print(f"\nLogistic regression (E₂, E₄, H) → β₁:")
    print(f"  Accuracy: {acc2:.4f}")
except ImportError:
    print("\nsklearn not available, skipping logistic regression")
    # Manual: check if E2 > threshold works
    for thresh in [-6, -3, 0, 3, 6]:
        tp = (idx1 & (E2 > thresh)).sum()
        fp = (idx0 & (E2 > thresh)).sum()
        fn = (idx1 & (E2 <= thresh)).sum()
        print(f"  E₂ > {thresh}: TP={tp}, FP={fp}, FN={fn}, precision={tp/(tp+fp+1e-10):.4f}, recall={tp/(tp+fn+1e-10):.4f}")

# Score sequence → β₁ rate
from collections import Counter
score_beta = Counter()
score_count = Counter()
for bits in range(N):
    A = tournament_from_bits(n, bits)
    score = tuple(sorted(A.sum(axis=1).astype(int)))
    score_beta[score] += all_beta1[bits]
    score_count[score] += 1

print(f"\nScore sequence → β₁ rate:")
for score in sorted(score_count.keys()):
    rate = score_beta[score] / score_count[score]
    print(f"  {score}: {score_beta[score]}/{score_count[score]} = {rate:.4f}")

# β₁=1 with lowest E₂
low_E2_beta1 = np.where(idx1)[0]
sorted_by_E2 = low_E2_beta1[np.argsort(E2[low_E2_beta1])]
print(f"\nβ₁=1 tournaments with LOWEST E₂:")
for bits in sorted_by_E2[:10]:
    A = tournament_from_bits(n, bits)
    score = tuple(sorted(A.sum(axis=1).astype(int)))
    t3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                s = A[i][j]*A[j][k]*A[k][i] + A[i][k]*A[k][j]*A[j][i]
                t3 += s
    print(f"  bits={bits}: E₂={E2[bits]:.1f}, E₄={E4[bits]:.1f}, H={all_H[bits]}, t3={t3}, score={score}")
