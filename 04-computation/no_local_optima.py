#!/usr/bin/env python3
"""
no_local_optima.py — opus-2026-03-13-S67j

VERIFYING: The H landscape on tournament flip graph has NO local optima.
Every local maximum of H(T) is a global maximum.

This is conjectured to hold for all n. We verify up to n=6 (m=15, 2^15 tournaments).

THEORETICAL FRAMEWORK (DPO Rewriting Connection):
- Each arc reversal is a double-pushout rewrite rule (Bajaj 2024)
- Two arc reversals COMMUTE iff their arcs don't share a vertex
- Non-commuting reversals form "causal conflicts"
- The interaction graph of causal conflicts = the LINE GRAPH L(K_n)
- This is exactly the support structure of degree-2 Fourier coefficients!

CONJECTURE (HYP-732): For all n, every local max of H on the tournament
flip graph is a global max. Equivalently, H has no "spurious" local optima.

CONJECTURE (HYP-733): H is "pseudo-convex" on the tournament polytope:
the set {T : H(T) ≥ c} is connected in the flip graph for every c.

This has implications for:
- Optimization: H-maximization can be solved by ANY local search
- Spin glass physics: the energy landscape is "unfrustrated"
- Coding theory: the code of H-maximizers has good covering properties
"""

import numpy as np
from itertools import permutations
import math
import time

def all_tournaments(n):
    edges = [(i,j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)
    tournaments = []
    for bits in range(2**m):
        A = np.zeros((n,n), dtype=int)
        for k, (i,j) in enumerate(edges):
            if (bits >> k) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        tournaments.append(A)
    return tournaments, edges

def count_ham_paths(A):
    n = A.shape[0]
    count = 0
    for perm in permutations(range(n)):
        valid = True
        for i in range(n-1):
            if A[perm[i]][perm[i+1]] != 1:
                valid = False
                break
        if valid:
            count += 1
    return count

def score_seq(A):
    return tuple(sorted(A.sum(axis=1).astype(int)))

# =====================================================================
# PART 1: VERIFY NO LOCAL OPTIMA FOR n=3,4,5
# =====================================================================
print("=" * 70)
print("NO LOCAL OPTIMA IN H LANDSCAPE")
print("=" * 70)

for n in [3, 4, 5]:
    t0 = time.time()
    tournaments, edges = all_tournaments(n)
    m = len(edges)
    N = 2**m

    H = np.array([count_ham_paths(A) for A in tournaments])
    global_max = np.max(H)
    global_min = np.min(H)

    # Check each tournament for local optimality
    local_max_count = 0
    local_max_not_global = 0
    local_min_count = 0
    local_min_not_global = 0

    # Also check: is every superlevel set {T : H(T) >= c} connected?
    for bits in range(N):
        is_local_max = True
        is_local_min = True
        for k in range(m):
            neighbor = bits ^ (1 << k)
            if H[neighbor] > H[bits]:
                is_local_max = False
            if H[neighbor] < H[bits]:
                is_local_min = False

        if is_local_max:
            local_max_count += 1
            if H[bits] < global_max:
                local_max_not_global += 1

        if is_local_min:
            local_min_count += 1
            if H[bits] > global_min:
                local_min_not_global += 1

    elapsed = time.time() - t0
    print(f"\n  n={n} (m={m}, N={N}, {elapsed:.1f}s):")
    print(f"    Global max: H={global_max}")
    print(f"    Global min: H={global_min}")
    print(f"    Local maxima: {local_max_count} (all global? {local_max_not_global == 0})")
    print(f"    Local minima: {local_min_count} (all global? {local_min_not_global == 0})")
    print(f"    Spurious local max: {local_max_not_global}")
    print(f"    Spurious local min: {local_min_not_global}")

    # Also analyze the H-value distribution
    unique_H = sorted(set(H))
    print(f"    Distinct H values: {unique_H}")

    # Superlevel set connectivity via BFS
    for threshold in unique_H[1:]:  # skip minimum
        superlevel = set(bits for bits in range(N) if H[bits] >= threshold)
        if len(superlevel) == 0:
            continue
        # BFS from first element
        start = next(iter(superlevel))
        visited = {start}
        queue = [start]
        while queue:
            current = queue.pop(0)
            for k in range(m):
                neighbor = current ^ (1 << k)
                if neighbor in superlevel and neighbor not in visited:
                    visited.add(neighbor)
                    queue.append(neighbor)
        connected = len(visited) == len(superlevel)
        if not connected:
            print(f"    *** DISCONNECTED superlevel at H>={threshold}: "
                  f"{len(visited)}/{len(superlevel)} reachable")
        # Only report disconnected ones
    print(f"    All superlevel sets connected? (checked all {len(unique_H)-1} thresholds)")

# =====================================================================
# PART 2: n=6 CHECK (use det(I+A) proxy first, then verify)
# =====================================================================
print("\n" + "=" * 70)
print("n=6 CHECK (m=15, N=32768)")
print("=" * 70)

n = 6
edges = [(i,j) for i in range(n) for j in range(i+1, n)]
m = len(edges)
N = 2**m

print(f"  Building all {N} tournaments and computing H via det(I+A)...")
t0 = time.time()

# det(I+A) is a PROXY for H (they're correlated but not identical at n>=6)
# For exact check at n=6, we'd need actual H computation
# det(I+A) = sum over all subsets S of det(A[S,S]) = number of "spanning subgraphs"
# but at n=6 we can actually compute H exactly since 6! = 720 permutations

H_6 = np.zeros(N, dtype=int)
detIpA_6 = np.zeros(N)

for bits in range(N):
    A = np.zeros((n,n), dtype=int)
    for k, (i,j) in enumerate(edges):
        if (bits >> k) & 1:
            A[i][j] = 1
        else:
            A[j][i] = 1

    # Exact H computation
    count = 0
    for perm in permutations(range(n)):
        valid = True
        for i in range(n-1):
            if A[perm[i]][perm[i+1]] != 1:
                valid = False
                break
        if valid:
            count += 1
    H_6[bits] = count
    detIpA_6[bits] = abs(np.linalg.det(np.eye(n) + A.astype(float)))

elapsed = time.time() - t0
print(f"  Computed in {elapsed:.1f}s")

global_max = np.max(H_6)
global_min = np.min(H_6)

# Check local optima
local_max_count = 0
local_max_not_global = 0
local_min_count = 0

for bits in range(N):
    is_local_max = True
    is_local_min = True
    for k in range(m):
        neighbor = bits ^ (1 << k)
        if H_6[neighbor] > H_6[bits]:
            is_local_max = False
        if H_6[neighbor] < H_6[bits]:
            is_local_min = False

    if is_local_max:
        local_max_count += 1
        if H_6[bits] < global_max:
            local_max_not_global += 1
    if is_local_min:
        local_min_count += 1

print(f"  Global max: H={global_max}")
print(f"  Global min: H={global_min}")
print(f"  Local maxima: {local_max_count} (all global? {local_max_not_global == 0})")
print(f"  Local minima: {local_min_count}")
print(f"  SPURIOUS LOCAL MAXIMA: {local_max_not_global}")

# Distinct H values
unique_H = sorted(set(H_6))
print(f"  Distinct H values: {unique_H}")

# Gradient flow structure: steepest ascent from transitive
# Find transitive tournament
for bits in range(N):
    A = np.zeros((n,n), dtype=int)
    for k, (i,j) in enumerate(edges):
        if (bits >> k) & 1:
            A[i][j] = 1
        else:
            A[j][i] = 1
    ss = score_seq(A)
    if ss == tuple(range(n)):
        trans_bits = bits
        break

# Steepest ascent
current = trans_bits
path = [current]
while True:
    best_k = -1
    best_delta = 0
    for k in range(m):
        neighbor = current ^ (1 << k)
        delta = H_6[neighbor] - H_6[current]
        if delta > best_delta:
            best_delta = delta
            best_k = k
    if best_k == -1:
        break
    current = current ^ (1 << best_k)
    path.append(current)

H_path = [int(H_6[b]) for b in path]
print(f"  Steepest ascent from transitive: {len(path)-1} steps, H path: {H_path}")

# =====================================================================
# PART 3: CAUSAL STRUCTURE OF ARC REVERSALS
# =====================================================================
print("\n" + "=" * 70)
print("CAUSAL STRUCTURE: ARC REVERSAL COMMUTATIVITY")
print("=" * 70)
print("  Two arc reversals COMMUTE iff their arcs share NO vertex.")
print("  The non-commutativity graph = LINE GRAPH L(K_n)")

for n in [3, 4, 5, 6, 7]:
    edges = [(i,j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)

    # Build line graph: edges in L(K_n) connect pairs sharing a vertex
    line_edges = 0
    independent_pairs = 0
    for a in range(m):
        for b in range(a+1, m):
            if edges[a][0] in edges[b] or edges[a][1] in edges[b]:
                line_edges += 1
            else:
                independent_pairs += 1

    print(f"\n  n={n}: m={m} arcs, L(K_n) has {line_edges} edges, {independent_pairs} independent pairs")
    print(f"    Independent fraction: {independent_pairs}/{m*(m-1)//2} = {independent_pairs/(m*(m-1)//2):.4f}")

    # Maximum independent set in L(K_n) = maximum matching in K_n = floor(n/2)
    max_matching = n // 2
    print(f"    Max commuting set: {max_matching} (= maximum matching in K_n = floor(n/2))")
    print(f"    => Can simultaneously reverse {max_matching} arcs without causal conflict")

# =====================================================================
# PART 4: LANDSCAPE DIAMETER AND GRADIENT STATISTICS
# =====================================================================
print("\n" + "=" * 70)
print("LANDSCAPE GEOMETRY: GRADIENT STATISTICS AT n=5")
print("=" * 70)

n = 5
tournaments, edges = all_tournaments(n)
m = len(edges)
N = 2**m

H = np.array([count_ham_paths(A) for A in tournaments])

# For each H-value level, compute average gradient magnitude
for h_val in sorted(set(H)):
    indices = [bits for bits in range(N) if H[bits] == h_val]
    grad_norms = []
    for bits in indices:
        grad = np.array([H[bits ^ (1 << k)] - H[bits] for k in range(m)])
        grad_norms.append(np.linalg.norm(grad))
    upward = []
    for bits in indices:
        up = sum(1 for k in range(m) if H[bits ^ (1 << k)] > H[bits])
        upward.append(up)

    print(f"  H={h_val:2d}: {len(indices):4d} tournaments, "
          f"|∇H|={np.mean(grad_norms):.2f}±{np.std(grad_norms):.2f}, "
          f"upward arcs: {np.mean(upward):.1f}±{np.std(upward):.1f}")

# At global max, gradient points INWARD (all neighbors have lower H)
# => regular tournaments are at a "peak" with no way up
# At transitive tournament, MOST arcs improve H

print("\n" + "=" * 70)
print("PATTERN: LOCAL MAX COUNT = GLOBAL MAX COUNT")
print("=" * 70)
for n in [3, 4, 5]:
    tournaments, edges = all_tournaments(n)
    m = len(edges)
    H = np.array([count_ham_paths(A) for A in tournaments])
    global_max = np.max(H)
    count_max = np.sum(H == global_max)
    # Compare with: n! * (number of regular tournaments) / |Aut|
    print(f"  n={n}: {count_max} global max tournaments (H={global_max})")

print(f"  n=6: {np.sum(H_6 == np.max(H_6))} global max tournaments (H={np.max(H_6)})")


print("\n\nDONE — no_local_optima.py complete")
