#!/usr/bin/env python3
"""
full_cg_iminus1.py — opus-2026-03-14-S76

CRITICAL TEST: Does I(-1) ≤ 1 hold for the FULL conflict graph
(all odd cycles), or only for the 3-cycle subgraph?

Background:
- The OCF says H(T) = I(Ω(T), 2) where Ω(T) has ALL odd cycles as vertices
- α₁ = total odd cycles, α₂ = disjoint pairs, etc.
- Previous work showed α₁(3-cyc) ≥ α₂(3,3 pairs) for n ≤ 9
- But the 5-cycle extension script showed α₂(full) >> α₁(full) at n=8!
- This script checks: does the FULL I(-1) = 1 - α₁ + α₂ - α₃ + ... ≤ 1?

Method: exhaustive at n=5,6, sampling at n=7,8,9
"""

from itertools import combinations, permutations
import random

def find_all_odd_cycles(adj, n, max_length=None):
    """Find all chordless directed odd cycles up to given length.
    Returns list of frozensets of vertices."""
    if max_length is None:
        max_length = n

    all_cycles = []

    for length in range(3, max_length + 1, 2):  # 3, 5, 7, ...
        if length > n:
            break
        for combo in combinations(range(n), length):
            verts = list(combo)
            fs = frozenset(combo)
            # Try to find a directed cycle on these vertices
            found = False
            for perm in permutations(verts):
                if found:
                    break
                # Check if perm is a directed cycle
                is_cycle = True
                for idx in range(length):
                    if not (adj[perm[idx]] & (1 << perm[(idx+1) % length])):
                        is_cycle = False
                        break
                if is_cycle:
                    # Check chordless: no "short" arcs
                    chordless = True
                    for idx in range(length):
                        for step in range(2, length - 1):
                            v1 = perm[idx]
                            v2 = perm[(idx + step) % length]
                            if adj[v1] & (1 << v2):
                                chordless = False
                                break
                        if not chordless:
                            break
                    if chordless:
                        all_cycles.append(fs)
                        found = True

    return all_cycles

def find_3cycles(adj, n):
    """Fast 3-cycle finder."""
    cycles = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (adj[i] & (1<<j)) and (adj[j] & (1<<k)) and (adj[k] & (1<<i)):
                    cycles.append(frozenset([i,j,k]))
                elif (adj[i] & (1<<k)) and (adj[k] & (1<<j)) and (adj[j] & (1<<i)):
                    cycles.append(frozenset([i,j,k]))
    return cycles

def compute_independence_polynomial(cycles_list):
    """Compute the independence polynomial coefficients α₀, α₁, α₂, ...
    by finding all independent sets in the conflict graph.

    The conflict graph has cycles_list as vertices,
    edges between cycles sharing a vertex.
    """
    m = len(cycles_list)
    if m == 0:
        return [1]  # just α₀ = 1

    # Build adjacency for conflict graph
    # conflict[i][j] = True if cycles i and j share a vertex
    conflict = [[False]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1, m):
            if len(cycles_list[i] & cycles_list[j]) > 0:
                conflict[i][j] = True
                conflict[j][i] = True

    # Find all independent sets by enumeration
    # For small m, use bitmask
    if m <= 20:
        max_k = 0
        alpha = {}  # alpha[k] = count of independent sets of size k
        alpha[0] = 1

        for mask in range(1, 1 << m):
            # Check if mask represents an independent set
            bits = []
            temp = mask
            idx = 0
            while temp:
                if temp & 1:
                    bits.append(idx)
                temp >>= 1
                idx += 1

            is_indep = True
            for i in range(len(bits)):
                for j in range(i+1, len(bits)):
                    if conflict[bits[i]][bits[j]]:
                        is_indep = False
                        break
                if not is_indep:
                    break

            if is_indep:
                k = len(bits)
                alpha[k] = alpha.get(k, 0) + 1
                max_k = max(max_k, k)

        result = [alpha.get(k, 0) for k in range(max_k + 1)]
        return result
    else:
        # For larger m, just compute α₁, α₂, α₃ directly
        alpha = [1, m]  # α₀=1, α₁=m

        # α₂ = disjoint pairs
        a2 = 0
        for i in range(m):
            for j in range(i+1, m):
                if not conflict[i][j]:
                    a2 += 1
        alpha.append(a2)

        # α₃ = disjoint triples
        a3 = 0
        for i in range(m):
            for j in range(i+1, m):
                if conflict[i][j]:
                    continue
                for k in range(j+1, m):
                    if not conflict[i][k] and not conflict[j][k]:
                        a3 += 1
        alpha.append(a3)

        return alpha

def ham_paths(adj, n):
    """Count Hamiltonian paths."""
    count = 0
    for perm in permutations(range(n)):
        is_path = True
        for idx in range(n-1):
            if not (adj[perm[idx]] & (1 << perm[idx+1])):
                is_path = False
                break
        if is_path:
            count += 1
    return count

def random_tournament(n):
    adj = [0] * n
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                adj[i] |= (1 << j)
            else:
                adj[j] |= (1 << i)
    return adj

# ====================================================================
print("=" * 70)
print("PART 1: EXHAUSTIVE CHECK AT n=5")
print("=" * 70)
print()

n = 5
edges = [(i,j) for i in range(n) for j in range(i+1,n)]
num_edges = len(edges)
total = 2 ** num_edges

max_iminus1 = -float('inf')
min_iminus1 = float('inf')
violations = 0
h_mismatch = 0

for bits in range(total):
    adj = [0] * n
    for idx, (i,j) in enumerate(edges):
        if bits & (1 << idx):
            adj[i] |= (1 << j)
        else:
            adj[j] |= (1 << i)

    # Find ALL odd cycles (at n=5: only 3 and 5)
    cycles = find_all_odd_cycles(adj, n)

    # Compute independence polynomial
    alpha = compute_independence_polynomial(cycles)

    # I(-1) = Σ (-1)^k α_k
    iminus1 = sum((-1)**k * alpha[k] for k in range(len(alpha)))

    # I(2) = Σ 2^k α_k
    i2 = sum(2**k * alpha[k] for k in range(len(alpha)))

    # H
    H = ham_paths(adj, n)

    if H != i2:
        h_mismatch += 1
        if h_mismatch <= 3:
            print(f"  MISMATCH: H={H}, I(2)={i2}, alpha={alpha}, cycles={len(cycles)}")

    if iminus1 > 1:
        violations += 1

    max_iminus1 = max(max_iminus1, iminus1)
    min_iminus1 = min(min_iminus1, iminus1)

print(f"n=5: {total} tournaments")
print(f"  H = I(2) mismatches: {h_mismatch}")
print(f"  I(-1) range: [{min_iminus1}, {max_iminus1}]")
print(f"  I(-1) > 1 violations: {violations}")
print(f"  I(-1) ≤ 1 for all? {violations == 0}")

# ====================================================================
print()
print("=" * 70)
print("PART 2: EXHAUSTIVE CHECK AT n=6")
print("=" * 70)
print()

n = 6
edges = [(i,j) for i in range(n) for j in range(i+1,n)]
num_edges = len(edges)
total = 2 ** num_edges

max_iminus1 = -float('inf')
min_iminus1 = float('inf')
violations = 0
h_mismatch = 0

chunk = total // 10
for bits in range(total):
    if bits % chunk == 0:
        print(f"  Progress: {bits}/{total}...", flush=True)

    adj = [0] * n
    for idx, (i,j) in enumerate(edges):
        if bits & (1 << idx):
            adj[i] |= (1 << j)
        else:
            adj[j] |= (1 << i)

    cycles = find_all_odd_cycles(adj, n)
    alpha = compute_independence_polynomial(cycles)
    iminus1 = sum((-1)**k * alpha[k] for k in range(len(alpha)))
    i2 = sum(2**k * alpha[k] for k in range(len(alpha)))
    H = ham_paths(adj, n)

    if H != i2:
        h_mismatch += 1
        if h_mismatch <= 3:
            c3 = [c for c in cycles if len(c) == 3]
            c5 = [c for c in cycles if len(c) == 5]
            print(f"  MISMATCH: H={H}, I(2)={i2}, #3cyc={len(c3)}, #5cyc={len(c5)}, alpha={alpha}")

    if iminus1 > 1:
        violations += 1

    max_iminus1 = max(max_iminus1, iminus1)
    min_iminus1 = min(min_iminus1, iminus1)

print(f"n=6: {total} tournaments")
print(f"  H = I(2) mismatches: {h_mismatch}")
print(f"  I(-1) range: [{min_iminus1}, {max_iminus1}]")
print(f"  I(-1) > 1 violations: {violations}")

# ====================================================================
print()
print("=" * 70)
print("PART 3: RANDOM SAMPLING AT n=7,8,9")
print("=" * 70)
print()

random.seed(42)
for n in [7, 8, 9]:
    nsamples = 200 if n <= 8 else 100
    max_iminus1 = -float('inf')
    min_iminus1 = float('inf')
    violations = 0
    h_mismatch = 0
    a1_gt_a2_count = 0

    for trial in range(nsamples):
        adj = random_tournament(n)
        cycles = find_all_odd_cycles(adj, n, max_length=n)

        # For efficiency at n=9, limit polynomial computation
        if len(cycles) <= 20:
            alpha = compute_independence_polynomial(cycles)
        else:
            alpha = compute_independence_polynomial(cycles)  # will use truncated version

        iminus1 = sum((-1)**k * alpha[k] for k in range(len(alpha)))
        i2 = sum(2**k * alpha[k] for k in range(len(alpha)))

        if n <= 8:
            H = ham_paths(adj, n)
            if H != i2:
                h_mismatch += 1
                if h_mismatch <= 2:
                    c_by_len = {}
                    for c in cycles:
                        l = len(c)
                        c_by_len[l] = c_by_len.get(l, 0) + 1
                    print(f"  n={n} MISMATCH at trial {trial}: H={H}, I(2)={i2}")
                    print(f"    Cycles by length: {c_by_len}")
                    print(f"    Alpha: {alpha}")

        if iminus1 > 1:
            violations += 1

        if len(alpha) >= 3 and alpha[1] >= alpha[2]:
            a1_gt_a2_count += 1

        max_iminus1 = max(max_iminus1, iminus1)
        min_iminus1 = min(min_iminus1, iminus1)

    print(f"n={n}: {nsamples} random tournaments")
    print(f"  H = I(2) mismatches: {h_mismatch}" + (" (not checked)" if n > 8 else ""))
    print(f"  I(-1) range: [{min_iminus1}, {max_iminus1}]")
    print(f"  I(-1) > 1 violations: {violations}/{nsamples}")
    print(f"  α₁ ≥ α₂ (full CG): {a1_gt_a2_count}/{nsamples}")
    print()

# ====================================================================
print()
print("=" * 70)
print("PART 4: WHAT IS GOING WRONG? CHORDLESS VS ALL CYCLES")
print("=" * 70)
print()
print("KEY QUESTION: Does the OCF use CHORDLESS odd cycles or ALL odd cycles?")
print()
print("The definition says 'directed odd cycles' — which could mean:")
print("  (a) All directed odd cycles (including those with chords)")
print("  (b) Only CHORDLESS directed odd cycles")
print()
print("Testing both interpretations at n=5:")

n = 5

def find_all_odd_cycles_with_chords(adj, n):
    """Find ALL directed odd cycles (including non-chordless)."""
    all_cycles = []
    for length in range(3, n + 1, 2):
        if length > n:
            break
        for combo in combinations(range(n), length):
            verts = list(combo)
            fs = frozenset(combo)
            found = False
            for perm in permutations(verts):
                if found:
                    break
                is_cycle = True
                for idx in range(length):
                    if not (adj[perm[idx]] & (1 << perm[(idx+1) % length])):
                        is_cycle = False
                        break
                if is_cycle:
                    all_cycles.append(fs)
                    found = True
    return all_cycles

# Test on a few tournaments
random.seed(99)
for trial in range(5):
    adj = random_tournament(n)
    H = ham_paths(adj, n)

    # Chordless
    cycles_cl = find_all_odd_cycles(adj, n)
    alpha_cl = compute_independence_polynomial(cycles_cl)
    i2_cl = sum(2**k * alpha_cl[k] for k in range(len(alpha_cl)))
    im1_cl = sum((-1)**k * alpha_cl[k] for k in range(len(alpha_cl)))

    # All (with chords)
    cycles_all = find_all_odd_cycles_with_chords(adj, n)
    alpha_all = compute_independence_polynomial(cycles_all)
    i2_all = sum(2**k * alpha_all[k] for k in range(len(alpha_all)))
    im1_all = sum((-1)**k * alpha_all[k] for k in range(len(alpha_all)))

    print(f"  T{trial}: H={H}")
    print(f"    Chordless: {len(cycles_cl)} cycles, I(2)={i2_cl}, I(-1)={im1_cl}, alpha={alpha_cl}")
    print(f"    All:       {len(cycles_all)} cycles, I(2)={i2_all}, I(-1)={im1_all}, alpha={alpha_all}")
    print(f"    Match chordless? {H == i2_cl}  Match all? {H == i2_all}")
    print()

print()
print("=" * 70)
print("SUMMARY")
print("=" * 70)
print()
print("This script determines:")
print("1. Whether H = I(Ω(T), 2) uses chordless or all odd cycles")
print("2. Whether I(-1) ≤ 1 holds for the FULL conflict graph")
print("3. Whether α₁ ≥ α₂ holds for the full CG")
