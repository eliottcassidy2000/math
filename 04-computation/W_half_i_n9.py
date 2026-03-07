"""
W_half_i_n9.py -- Test whether W(i/2) = 0 characterizes H-maximizers at n=9.

W(r) = sum over permutations P of prod_{j=0}^{n-2} (r + A[P(j),P(j+1)] - 1/2)

We evaluate W at r=1/2 (gives H(T)) and r=i/2 (complex).
"""

import random
import time

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
                A[j][i] = 0
            else:
                A[i][j] = 0
                A[j][i] = 1
    return A

def circulant_tournament(n, gen_set):
    """Circulant tournament: vertex i beats i+g mod n for g in gen_set."""
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for g in gen_set:
            j = (i + g) % n
            A[i][j] = 1
    return A

def compute_W(A, n, r):
    """DP approach for W(r) = sum over Hamiltonian paths of product of edge weights."""
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = complex(1, 0) if isinstance(r, complex) else 1.0
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            val = dp.get((mask, v), 0)
            if val == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                new_mask = mask | (1 << u)
                dp[(new_mask, u)] = (
                    dp.get((new_mask, u), 0) + val * (r + A[v][u] - 0.5)
                )
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def tournament_key(A, n):
    """Canonical key for deduplication (just the upper triangle bits)."""
    bits = []
    for i in range(n):
        for j in range(i+1, n):
            bits.append(str(A[i][j]))
    return ''.join(bits)

def score_sequence(A, n):
    scores = sorted([sum(A[i]) for i in range(n)])
    return tuple(scores)

def is_regular(A, n):
    target = (n - 1) / 2
    return all(abs(sum(A[i]) - target) < 0.01 for i in range(n))

n = 9
half_i = complex(0, 0.5)  # i/2
half = 0.5

print("=" * 70)
print("W(i/2) = 0 characterization of H-maximizers at n=9")
print("=" * 70)

# --- Part 1: Specific circulant tournaments ---
print("\n--- Circulant tournaments on 9 vertices ---")
# All possible generating sets of size 4 from {1,2,3,4,5,6,7,8}
# For a tournament, if g is in gen_set, then n-g must NOT be (they'd conflict)
# So gen_set picks one from each pair {1,8},{2,7},{3,6},{4,5}
from itertools import product as iprod

circulant_results = []
for choices in iprod([0, 1], repeat=4):
    pairs = [(1,8),(2,7),(3,6),(4,5)]
    gen_set = [pairs[i][choices[i]] for i in range(4)]
    A = circulant_tournament(n, gen_set)
    H = int(round(compute_W(A, n, half).real))
    W_hi = compute_W(A, n, half_i)
    reg = is_regular(A, n)
    circulant_results.append((gen_set, H, W_hi, reg))
    print(f"  gen={gen_set}, H={H}, W(i/2)={W_hi.real:.4f}+{W_hi.imag:.4f}i, regular={reg}")

best_circ_H = max(r[1] for r in circulant_results)
print(f"\nBest circulant H = {best_circ_H}")

# --- Part 2: Random sampling ---
print("\n--- Random sampling: 200 tournaments at n=9 ---")
N_SAMPLES = 200

results = []
t0 = time.time()

for trial in range(N_SAMPLES):
    A = random_tournament(n)
    H_val = compute_W(A, n, half)
    H = int(round(H_val.real))
    W_hi = compute_W(A, n, half_i)
    results.append((H, W_hi, A))
    if (trial + 1) % 50 == 0:
        elapsed = time.time() - t0
        print(f"  {trial+1}/{N_SAMPLES} done, elapsed={elapsed:.1f}s")

# Also add circulant tournaments to the pool
for gen_set, H, W_hi, reg in circulant_results:
    A = circulant_tournament(n, gen_set)
    results.append((H, W_hi, A))

elapsed = time.time() - t0
print(f"Total time: {elapsed:.1f}s")

# --- Analysis ---
max_H = max(r[0] for r in results)
print(f"\nMax H found: {max_H}")

# Maximizers
maximizers = [(H, W, A) for H, W, A in results if H == max_H]
print(f"Number of max-H tournaments: {len(maximizers)}")

for i, (H, W, A) in enumerate(maximizers[:5]):
    sc = score_sequence(A, n)
    print(f"  Maximizer {i}: H={H}, score={sc}, W(i/2)={W.real:.6f}+{W.imag:.6f}i, |W(i/2)|={abs(W):.6f}")

# W(i/2) near zero
tol = 0.1
zero_W = [(H, W, A) for H, W, A in results if abs(W) < tol]
print(f"\nTournaments with |W(i/2)| < {tol}: {len(zero_W)} out of {len(results)}")
if zero_W:
    zero_H_vals = sorted(set(r[0] for r in zero_W))
    print(f"  Their H values: {zero_H_vals}")
    print(f"  Max H among W(i/2)~0: {max(r[0] for r in zero_W)}")
    # Check: does W(i/2)=0 imply H=max?
    all_max = all(r[0] == max_H for r in zero_W)
    print(f"  All W(i/2)~0 have max H? {all_max}")

# Check: do all maximizers have W(i/2)=0?
max_W_vals = [abs(W) for H, W, _ in results if H == max_H]
print(f"\n|W(i/2)| for maximizers: min={min(max_W_vals):.6f}, max={max(max_W_vals):.6f}")
all_max_zero = all(v < tol for v in max_W_vals)
print(f"All maximizers have |W(i/2)| < {tol}? {all_max_zero}")

# Correlation between H and |W(i/2)|
H_vals = [r[0] for r in results]
absW_vals = [abs(r[1]) for r in results]
realW_vals = [r[1].real for r in results]

mean_H = sum(H_vals) / len(H_vals)
mean_absW = sum(absW_vals) / len(absW_vals)
mean_realW = sum(realW_vals) / len(realW_vals)

cov_H_absW = sum((H_vals[i] - mean_H) * (absW_vals[i] - mean_absW) for i in range(len(results))) / len(results)
var_H = sum((h - mean_H)**2 for h in H_vals) / len(results)
var_absW = sum((w - mean_absW)**2 for w in absW_vals) / len(results)
if var_H > 0 and var_absW > 0:
    corr_H_absW = cov_H_absW / (var_H**0.5 * var_absW**0.5)
else:
    corr_H_absW = 0
print(f"\nCorrelation(H, |W(i/2)|) = {corr_H_absW:.4f}")

cov_H_realW = sum((H_vals[i] - mean_H) * (realW_vals[i] - mean_realW) for i in range(len(results))) / len(results)
var_realW = sum((w - mean_realW)**2 for w in realW_vals) / len(results)
if var_H > 0 and var_realW > 0:
    corr_H_realW = cov_H_realW / (var_H**0.5 * var_realW**0.5)
else:
    corr_H_realW = 0
print(f"Correlation(H, Re(W(i/2))) = {corr_H_realW:.4f}")

# Top-10 by H
print("\n--- Top 10 by H ---")
sorted_results = sorted(results, key=lambda r: -r[0])
seen = set()
count = 0
for H, W, A in sorted_results:
    k = tournament_key(A, n)
    if k in seen:
        continue
    seen.add(k)
    sc = score_sequence(A, n)
    reg = is_regular(A, n)
    print(f"  H={H}, |W(i/2)|={abs(W):.4f}, W(i/2)=({W.real:.4f}, {W.imag:.4f}), score={sc}, reg={reg}")
    count += 1
    if count >= 10:
        break

# Bottom-10 by H
print("\n--- Bottom 10 by H ---")
seen = set()
count = 0
for H, W, A in reversed(sorted_results):
    k = tournament_key(A, n)
    if k in seen:
        continue
    seen.add(k)
    sc = score_sequence(A, n)
    print(f"  H={H}, |W(i/2)|={abs(W):.4f}, W(i/2)=({W.real:.4f}, {W.imag:.4f}), score={sc}")
    count += 1
    if count >= 10:
        break

# --- Part 3: Try QR tournament (doubly regular) ---
print("\n--- Quadratic residue tournament on Z/9Z (if applicable) ---")
print("9 is not prime, so no standard Paley/QR tournament exists.")
print("The circulant {1,2,3,4} is the natural candidate for 'most regular'.")

# Summarize
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"n = {n}")
print(f"Samples: {len(results)} (200 random + {len(circulant_results)} circulants)")
print(f"Max H found: {max_H}")
print(f"Best circulant H: {best_circ_H}")
print(f"|W(i/2)| < 0.1 count: {len(zero_W)}")
if zero_W:
    print(f"W(i/2)~0 => max H? {all(r[0] == max_H for r in zero_W)}")
print(f"All maximizers have W(i/2)~0? {all_max_zero}")
print(f"Corr(H, |W(i/2)|) = {corr_H_absW:.4f}")
print(f"Corr(H, Re(W(i/2))) = {corr_H_realW:.4f}")
