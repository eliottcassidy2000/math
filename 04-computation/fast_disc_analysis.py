#!/usr/bin/env python3
"""
Fast discriminant analysis using DFS-based cycle finding.

Key questions:
1. Find more negative-discriminant n=9 tournaments
2. Characterize the inter-group structure
3. n=10 Omega_3 analysis (3-cycles only, fast)

Author: opus-2026-03-06-S19
"""
import numpy as np
from collections import Counter
import random

random.seed(42)
np.random.seed(42)

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def find_odd_cycles_dfs(A, n, max_len=None):
    """DFS-based directed odd cycle finder. Much faster than permutation enumeration."""
    if max_len is None:
        max_len = n
    cycles = set()

    def dfs(path, start):
        cur = path[-1]
        L = len(path)
        if L > 1 and A[cur][start] and L % 2 == 1 and L <= max_len:
            min_idx = path.index(min(path))
            normalized = tuple(path[min_idx:] + path[:min_idx])
            cycles.add(normalized)
        if L < max_len:
            for nxt in range(n):
                if nxt != start and nxt not in path and A[cur][nxt]:
                    dfs(path + [nxt], start)

    for v in range(n):
        dfs([v], v)
    return list(cycles)

def find_3cycles(A, n):
    """Find all directed 3-cycles."""
    cycles = set()
    for i in range(n):
        for j in range(n):
            if j == i or not A[i][j]: continue
            for k in range(n):
                if k == i or k == j or not A[j][k]: continue
                if A[k][i]:
                    c = (i, j, k)
                    mi = c.index(min(c))
                    cycles.add(c[mi:] + c[:mi])
    return list(cycles)

def independence_poly_dc(adj_dict):
    memo = {}
    def solve(verts):
        if verts in memo:
            return memo[verts]
        if not verts:
            return [1]
        v = max(verts, key=lambda u: len(adj_dict[u] & verts))
        p1 = solve(verts - {v})
        p2 = solve(verts - (adj_dict[v] & verts) - {v})
        maxlen = max(len(p1), len(p2) + 1)
        result = [0] * maxlen
        for i in range(len(p1)):
            result[i] += p1[i]
        for i in range(len(p2)):
            result[i + 1] += p2[i]
        memo[verts] = result
        return result
    return solve(frozenset(adj_dict.keys()))

def build_omega_poly(cycles):
    m = len(cycles)
    if m == 0:
        return [1]
    vsets = [frozenset(c) for c in cycles]
    adj = {}
    for i in range(m):
        adj[i] = frozenset(j for j in range(m) if j != i and vsets[i] & vsets[j])
    return independence_poly_dc(adj)

def cubic_disc(coeffs):
    if len(coeffs) != 4:
        return None
    a0, a1, a2, a3 = coeffs
    return 18*a3*a2*a1*a0 - 4*a2**3*a0 + a2**2*a1**2 - 4*a3*a1**3 - 27*a3**2*a0**2

def newton_check(coeffs):
    deg = len(coeffs) - 1
    while deg > 0 and coeffs[deg] == 0:
        deg -= 1
    if deg < 2:
        return True
    for k in range(1, deg):
        if coeffs[k]**2 < coeffs[k-1] * coeffs[k+1] * (k+1) / k:
            return False
    return True

# ============================================================
print("=" * 70)
print("FAST DISCRIMINANT ANALYSIS")
print("=" * 70)

# ============================================================
# Part 1: Hunt negative discriminants at n=9 using full Omega (DFS-based)
# ============================================================
print("\n--- Part 1: Full Omega at n=9, targeted search ---")
n = 9
found = 0
tested = 0

for trial in range(3000):
    T = random_tournament(n)
    scores = sorted(sum(row) for row in T)
    score_var = sum((s - 4)**2 for s in scores) / n

    # Only test high-variance tournaments
    if score_var < 2.5:
        continue

    cycles = find_odd_cycles_dfs(T, n)
    m = len(cycles)
    if m == 0:
        continue
    tested += 1

    coeffs = build_omega_poly(cycles)

    newton_ok = newton_check(coeffs)
    d = cubic_disc(coeffs)

    if not newton_ok or (d is not None and d < 0):
        found += 1
        print(f"\n  FOUND #{found} at trial {trial}:")
        print(f"    Scores: {tuple(scores)}, var={score_var:.3f}")
        print(f"    I(Omega,x) = {coeffs}")
        print(f"    Newton: {'ok' if newton_ok else 'FAIL'}")
        print(f"    disc = {d}")
        print(f"    # cycles: {m}, by length: {dict(Counter(len(c) for c in cycles))}")
        if found >= 10:
            break

print(f"\nTested {tested} high-variance n=9 tournaments, found {found}")

# ============================================================
# Part 2: Between-group structure characterization
# ============================================================
print("\n--- Part 2: Between-group structure of counterexample ---")

T_CE = [
    [0, 1, 0, 1, 0, 0, 1, 1, 0],
    [0, 0, 0, 1, 0, 0, 0, 0, 0],
    [1, 1, 0, 0, 1, 1, 1, 1, 0],
    [0, 0, 1, 0, 0, 1, 0, 1, 0],
    [1, 1, 0, 1, 0, 0, 0, 1, 0],
    [1, 1, 0, 0, 1, 0, 1, 1, 1],
    [0, 1, 0, 1, 1, 0, 0, 1, 0],
    [0, 1, 0, 0, 0, 0, 0, 0, 0],
    [1, 1, 1, 1, 1, 0, 1, 1, 0],
]

# The three groups: {0,4,6}, {2,5,8}, {1,3,7}
# Between-group arcs form a 3x3 "group tournament"
groups = [[0, 4, 6], [2, 5, 8], [1, 3, 7]]

print("  Inter-group arc matrix (3x3):")
group_scores = [[0, 0, 0] for _ in range(3)]
for gi, G1 in enumerate(groups):
    for gj, G2 in enumerate(groups):
        if gi == gj:
            continue
        fwd = sum(T_CE[u][v] for u in G1 for v in G2)
        group_scores[gi][gj] = fwd

for gi in range(3):
    row = [group_scores[gi][gj] for gj in range(3)]
    print(f"    Group {gi} ({groups[gi]}): {row} (total out = {sum(row)})")

# The "group tournament" is: G0->G2 (9-0), G1->G0 (9-0), G1->G2 (7-2)
# So G1 > G0 > G2, with G1 also beating G2 directly
print("\n  Group tournament is ALMOST transitive: G1 >> G0 >> G2, with G1 >> G2")
print("  The extreme dominance (9-0) means crossing arcs only go one way,")
print("  forcing cycles to route through specific vertices.")

# ============================================================
# Part 3: Quantify "how extreme" the counterexample's group structure is
# ============================================================
print("\n--- Part 3: Group imbalance in random 3d3c tournaments ---")

# For tournaments with 3 vertex-disjoint 3-cycles, measure the
# maximum inter-group domination (max of |fwd - bwd| over group pairs)
imbalance_data = []

for trial in range(2000):
    T = random_tournament(n)
    c3 = find_3cycles(T, n)

    # Check for 3 disjoint 3-cycles
    found_triple = None
    vsets3 = [frozenset(c) for c in c3]
    for i in range(len(c3)):
        for j in range(i+1, len(c3)):
            if vsets3[i] & vsets3[j]:
                continue
            for k in range(j+1, len(c3)):
                if not (vsets3[i] & vsets3[k]) and not (vsets3[j] & vsets3[k]):
                    found_triple = (c3[i], c3[j], c3[k])
                    break
            if found_triple:
                break
        if found_triple:
            break

    if not found_triple:
        continue

    g = [list(found_triple[0]), list(found_triple[1]), list(found_triple[2])]

    # Compute max imbalance
    max_imb = 0
    for gi in range(3):
        for gj in range(gi+1, 3):
            fwd = sum(T[u][v] for u in g[gi] for v in g[gj])
            bwd = sum(T[v][u] for u in g[gi] for v in g[gj])
            max_imb = max(max_imb, abs(fwd - bwd))

    # Full Omega analysis
    cycles = find_odd_cycles_dfs(T, n)
    coeffs = build_omega_poly(cycles)
    newton_ok = newton_check(coeffs)
    d = cubic_disc(coeffs)

    fail = not newton_ok or (d is not None and d < 0)
    imbalance_data.append((max_imb, fail, coeffs, len(cycles)))

print(f"  {len(imbalance_data)} 3d3c tournaments analyzed")

# Group by max imbalance
imb_groups = {}
for imb, fail, coeffs, nc in imbalance_data:
    if imb not in imb_groups:
        imb_groups[imb] = {'total': 0, 'fail': 0}
    imb_groups[imb]['total'] += 1
    if fail:
        imb_groups[imb]['fail'] += 1

print(f"\n  Max imbalance | Total | Failures | Rate")
for imb in sorted(imb_groups.keys()):
    g = imb_groups[imb]
    rate = 100 * g['fail'] / g['total'] if g['total'] > 0 else 0
    marker = " ***" if g['fail'] > 0 else ""
    print(f"  {imb:>13d} | {g['total']:>5d} | {g['fail']:>8d} | {rate:>5.2f}%{marker}")

# Show the failure cases
failures = [(imb, fail, coeffs, nc) for imb, fail, coeffs, nc in imbalance_data if fail]
print(f"\n  Total failures: {len(failures)}")
for imb, _, coeffs, nc in failures:
    d = cubic_disc(coeffs)
    print(f"    max_imb={imb}, coeffs={coeffs}, disc={d}, #cyc={nc}")

# ============================================================
# Part 4: n=10 Omega_3 analysis
# ============================================================
print("\n" + "=" * 70)
print("PART 4: n=10 Omega_3 ANALYSIS")
print("=" * 70)

n10 = 10
newton_fail_10 = 0
disc_neg_10 = 0
total_10 = 0
deg_dist = Counter()
all_coeffs_10 = []

for trial in range(500):
    T = random_tournament(n10)
    c3 = find_3cycles(T, n10)
    m3 = len(c3)
    if m3 == 0:
        continue
    total_10 += 1

    coeffs = build_omega_poly(c3)
    deg = len(coeffs) - 1
    deg_dist[deg] += 1

    newton_ok = newton_check(coeffs)
    if not newton_ok:
        newton_fail_10 += 1

    d = cubic_disc(coeffs)
    if d is not None and d < 0:
        disc_neg_10 += 1

    all_coeffs_10.append((coeffs, newton_ok, d))

print(f"  n=10, {total_10} tournaments with 3-cycles out of 500")
print(f"  Degree distribution: {dict(sorted(deg_dist.items()))}")
print(f"  Newton failures: {newton_fail_10}")
print(f"  Disc < 0 (degree 3): {disc_neg_10}")

# At n=10, max independent set of 3-cycles = 3 (need 9 vertices for 3 disjoint 3-cycles)
# Can we get 3 disjoint 3-cycles? Yes, at n=10 just need 3 vertex-disjoint ones out of 10 vertices
# alpha(Omega_3) can be at most 3 since 4 * 3 = 12 > 10
print(f"\n  At n=10: alpha(Omega_3) <= 3 (12 > 10)")
print(f"  Degree 3 polynomials: {deg_dist.get(3, 0)}")

# ============================================================
# Part 5: When is alpha_3(Omega_3) = 1 vs 2 vs 3?
# ============================================================
print("\n" + "=" * 70)
print("PART 5: INDEPENDENCE NUMBER DISTRIBUTION AT n=10")
print("=" * 70)

alpha_dist = Counter()
for coeffs, _, _ in all_coeffs_10:
    deg = len(coeffs) - 1
    alpha_dist[deg] += 1

print(f"  alpha(Omega_3) distribution:")
for deg in sorted(alpha_dist.keys()):
    pct = 100 * alpha_dist[deg] / total_10
    print(f"    alpha = {deg}: {alpha_dist[deg]} ({pct:.1f}%)")

# For degree-3 cases at n=10, show coefficient ranges
deg3_10 = [(c, n, d) for c, n, d in all_coeffs_10 if len(c) == 4]
if deg3_10:
    a1s = [c[1] for c, _, _ in deg3_10]
    a2s = [c[2] for c, _, _ in deg3_10]
    a3s = [c[3] for c, _, _ in deg3_10]
    print(f"\n  Degree-3 coefficient ranges at n=10:")
    print(f"    a1: [{min(a1s)}, {max(a1s)}], avg={np.mean(a1s):.1f}")
    print(f"    a2: [{min(a2s)}, {max(a2s)}], avg={np.mean(a2s):.1f}")
    print(f"    a3: [{min(a3s)}, {max(a3s)}], avg={np.mean(a3s):.1f}")

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
