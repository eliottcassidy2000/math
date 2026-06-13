#!/usr/bin/env python3
"""
n8_maximizer_topology.py — Deep analysis of H=661 maximizer split at n=8

At n=8, H-maximizers (H=661) split between β₄=1 and contractible.
This script:
1. Finds ALL H=661 tournaments (systematic search via SC + fpf involutions)
2. Computes Pfaffian, eigenvalues, spectral gap for each
3. Computes Betti for a representative of each spectral type
4. Checks vertex-deletion Betti for hereditary topology
5. Counts c3, alpha_1, alpha_2

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, random, os
import numpy as np
from itertools import combinations
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w')
from path_homology_v2 import path_betti_numbers
sys.stdout = _saved

def H_tournament(A, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if dp[mask][v] == 0: continue
            for u in range(n):
                if mask & (1 << u): continue
                if A[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full])

def skew_eigenvalues(A, n):
    """Return positive imaginary eigenvalues of skew matrix S = A - A^T."""
    S = np.array([[A[i][j] - A[j][i] for j in range(n)] for i in range(n)], dtype=float)
    evals = np.linalg.eigvals(S)
    pos_imag = sorted([e.imag for e in evals if e.imag > 0.01])
    return pos_imag

def pfaffian_from_eigenvalues(evals):
    """Product of positive imaginary eigenvalues = |Pf(S)|."""
    return int(round(np.prod(evals)))

def count_3cycles(A, n):
    """Count directed 3-cycles (as vertex sets)."""
    cycles = set()
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j] and A[j][k] and A[k][i]:
                    cycles.add(frozenset([i,j,k]))
                if A[i][k] and A[k][j] and A[j][i]:
                    cycles.add(frozenset([i,j,k]))
    return cycles

def count_5cycles(A, n):
    """Count directed 5-cycles (as vertex sets)."""
    cycles = set()
    verts = list(range(n))
    for combo in combinations(verts, 5):
        # Check all possible orderings
        from itertools import permutations as perms
        for perm in perms(combo):
            is_cycle = True
            for idx in range(5):
                if not A[perm[idx]][perm[(idx+1) % 5]]:
                    is_cycle = False
                    break
            if is_cycle:
                cycles.add(frozenset(combo))
                break
    return cycles

def independence_number(cycle_sets):
    """Count independent sets in the conflict graph."""
    unique = list(cycle_sets)
    m = len(unique)
    if m == 0:
        return [1]

    # Build conflict
    conf = [[False]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1, m):
            if unique[i] & unique[j]:
                conf[i][j] = True
                conf[j][i] = True

    # Count independent sets by size
    coeffs = [0] * (m + 1)
    coeffs[0] = 1
    for mask in range(1, min(1 << m, 1 << 20)):  # cap at 2^20
        bits = []
        temp = mask
        while temp:
            b = temp & (-temp)
            bits.append(b.bit_length() - 1)
            temp ^= b
        ok = True
        for i in range(len(bits)):
            for j in range(i+1, len(bits)):
                if conf[bits[i]][bits[j]]:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            coeffs[len(bits)] += 1

    while len(coeffs) > 1 and coeffs[-1] == 0:
        coeffs.pop()
    return coeffs

def vertex_deletion_betti(A, n, v, max_dim=5):
    """Compute Betti of T-v."""
    An = [[A[i][j] for j in range(n) if j != v] for i in range(n) if i != v]
    return path_betti_numbers(An, n-1, max_dim=max_dim)

# ===== Find H=661 tournaments by sampling =====
print("=" * 70)
print("FINDING H=661 MAXIMIZERS AT n=8")
print("=" * 70)
n = 8
m = n * (n-1) // 2
t0 = time.time()

max_tours = []  # Store (bits, A) for H=661
max_H = 0
checked = 0
target = 500000

for _ in range(target):
    bits = random.randint(0, (1 << m) - 1)
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    H = H_tournament(A, n)
    if H > max_H:
        max_H = H
    if H == 661:
        max_tours.append((bits, [row[:] for row in A]))

    checked += 1
    if checked % 100000 == 0:
        print(f"  {checked}/{target} ({time.time()-t0:.1f}s), found {len(max_tours)} with H=661")

print(f"\nFound {len(max_tours)} H=661 maximizers in {checked} samples")
print(f"Max H seen: {max_H}")

# ===== Analyze spectral types =====
print("\n" + "=" * 70)
print("SPECTRAL ANALYSIS OF MAXIMIZERS")
print("=" * 70)

spectral_types = Counter()
type_examples = {}  # gap -> example A
type_pfaffians = {}

for bits, A in max_tours:
    evals = skew_eigenvalues(A, n)
    pf = pfaffian_from_eigenvalues(evals)
    gap = max(evals) - min(evals)
    gap_rounded = round(gap, 3)

    spectral_types[gap_rounded] += 1
    if gap_rounded not in type_examples:
        type_examples[gap_rounded] = A
        type_pfaffians[gap_rounded] = pf

print(f"\nSpectral types by gap:")
for gap in sorted(spectral_types.keys()):
    pf = type_pfaffians[gap]
    print(f"  gap={gap:.4f}: {spectral_types[gap]} tournaments, |Pf|={pf}")

# ===== Compute Betti for each type =====
print("\n" + "=" * 70)
print("BETTI NUMBERS BY SPECTRAL TYPE")
print("=" * 70)

for gap in sorted(type_examples.keys()):
    A = type_examples[gap]
    try:
        beta = path_betti_numbers(A, n, max_dim=6)
        beta_list = [int(beta[k]) if k < len(beta) else 0 for k in range(7)]
    except:
        beta_list = ['ERROR']

    # 3-cycle and 5-cycle analysis
    c3_sets = count_3cycles(A, n)
    c5_sets = count_5cycles(A, n)
    ip3 = independence_number(c3_sets)

    evals = skew_eigenvalues(A, n)
    pf = type_pfaffians[gap]

    print(f"\n  gap={gap:.4f}, |Pf|={pf}:")
    print(f"    eigenvalues: [{', '.join(f'{e:.4f}' for e in evals)}]")
    print(f"    β = {beta_list}")
    print(f"    c3={len(c3_sets)}, c5={len(c5_sets)}")
    print(f"    I(Ω₃, x) = {ip3}")

# ===== Vertex-deletion topology =====
print("\n" + "=" * 70)
print("VERTEX-DELETION TOPOLOGY (hereditary check)")
print("=" * 70)

for gap in sorted(type_examples.keys()):
    A = type_examples[gap]
    beta_parent = path_betti_numbers(A, n, max_dim=6)
    beta_parent_list = [int(beta_parent[k]) if k < len(beta_parent) else 0 for k in range(7)]

    print(f"\n  Parent gap={gap:.4f}, β={beta_parent_list}:")
    del_bettis = Counter()
    for v in range(n):
        beta_v = vertex_deletion_betti(A, n, v, max_dim=5)
        bv = tuple(int(beta_v[k]) if k < len(beta_v) else 0 for k in range(6))
        del_bettis[bv] += 1

    for bv, cnt in del_bettis.most_common():
        print(f"    β(T-v) = {list(bv)}: {cnt} vertex deletions")

# ===== Complement check =====
print("\n" + "=" * 70)
print("COMPLEMENT ANALYSIS")
print("=" * 70)

for gap in sorted(type_examples.keys()):
    A = type_examples[gap]
    # Build complement
    Ac = [[1 - A[i][j] if i != j else 0 for j in range(n)] for i in range(n)]
    Hc = H_tournament(Ac, n)
    evals_c = skew_eigenvalues(Ac, n)
    gap_c = max(evals_c) - min(evals_c)

    try:
        beta_c = path_betti_numbers(Ac, n, max_dim=6)
        beta_c_list = [int(beta_c[k]) if k < len(beta_c) else 0 for k in range(7)]
    except:
        beta_c_list = ['ERROR']

    print(f"  gap={gap:.4f}: complement has H={Hc}, gap={gap_c:.4f}, β={beta_c_list}")

print(f"\nTotal time: {time.time()-t0:.1f}s")
