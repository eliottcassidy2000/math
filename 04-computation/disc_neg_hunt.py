#!/usr/bin/env python3
"""
Hunt for the 2 negative-discriminant tournaments found in 3d3c analysis.
Also test n=10 failure rate.

Author: opus-2026-03-06-S19
"""
import random
from itertools import combinations, permutations
from collections import Counter

random.seed(42)

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def find_all_directed_odd_cycles(A, n):
    cycles = set()
    for length in range(3, n + 1, 2):
        for subset in combinations(range(n), length):
            for perm in permutations(subset[1:]):
                cycle = (subset[0],) + perm
                if all(A[cycle[i]][cycle[(i+1) % length]] for i in range(length)):
                    min_idx = cycle.index(min(cycle))
                    rotated = cycle[min_idx:] + cycle[:min_idx]
                    cycles.add(rotated)
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

def cubic_disc(a0, a1, a2, a3):
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
print("HUNTING NEGATIVE-DISCRIMINANT TOURNAMENTS AT n=9")
print("=" * 70)

n = 9
found = 0

for trial in range(2000):
    T = random_tournament(n)
    cycles = find_all_directed_odd_cycles(T, n)
    m = len(cycles)
    if m == 0:
        continue

    vsets = [frozenset(c) for c in cycles]
    adj = {}
    for i in range(m):
        adj[i] = frozenset(j for j in range(m) if j != i and vsets[i] & vsets[j])
    coeffs = independence_poly_dc(adj)

    if len(coeffs) == 4:
        a0, a1, a2, a3 = coeffs
        d = cubic_disc(a0, a1, a2, a3)
        newton_ok = newton_check(coeffs)
        if d < 0:
            found += 1
            scores = tuple(sorted(sum(row) for row in T))
            score_var = sum((s - 4)**2 for s in scores) / n
            print(f"\n  FOUND #{found} at trial {trial}:")
            print(f"    Scores: {scores}, var={score_var:.3f}")
            print(f"    I(Omega,x) = {coeffs}")
            print(f"    disc = {d}")
            print(f"    Newton: {'ok' if newton_ok else 'FAIL'}")
            print(f"    # cycles: {m}")
            lens = Counter(len(c) for c in cycles)
            print(f"    Cycle lengths: {dict(lens)}")

            # Between-group structure
            # Find the independent triple
            ind_pairs = []
            for i in range(m):
                for j in range(i+1, m):
                    if j not in adj[i]:
                        ind_pairs.append((i, j))
            print(f"    Independent pairs: {len(ind_pairs)}")
            # Independent triples
            ind_triples = []
            for i, j in ind_pairs:
                for k in range(j+1, m):
                    if k not in adj[i] and k not in adj[j]:
                        ind_triples.append((i, j, k))
            print(f"    Independent triples: {len(ind_triples)}")

            if found >= 5:
                break

print(f"\nTotal found: {found} in 2000 random n=9 tournaments")

# ============================================================
# Focused search: high-variance scores
# ============================================================
print("\n" + "=" * 70)
print("TARGETED SEARCH: High-variance scores at n=9")
print("=" * 70)

found_targeted = 0
for trial in range(5000):
    # Generate high-variance tournament
    T = random_tournament(n)
    scores = sorted(sum(row) for row in T)
    score_var = sum((s - 4)**2 for s in scores) / n
    if score_var < 3.0:
        continue

    cycles = find_all_directed_odd_cycles(T, n)
    m = len(cycles)
    if m == 0:
        continue

    vsets = [frozenset(c) for c in cycles]
    adj = {}
    for i in range(m):
        adj[i] = frozenset(j for j in range(m) if j != i and vsets[i] & vsets[j])
    coeffs = independence_poly_dc(adj)

    if len(coeffs) >= 4:
        newton_ok = newton_check(coeffs)
        if not newton_ok:
            found_targeted += 1
            d = cubic_disc(*coeffs[:4]) if len(coeffs) == 4 else None
            print(f"  FOUND: scores={tuple(scores)}, var={score_var:.3f}, "
                  f"coeffs={coeffs}, Newton=FAIL, disc={d}, #cyc={m}")
            if found_targeted >= 5:
                break

print(f"Total Newton failures in targeted search: {found_targeted}")

# ============================================================
# n=10 exploration (very quick)
# ============================================================
print("\n" + "=" * 70)
print("n=10 EXPLORATION (Newton check on Omega_3 only)")
print("=" * 70)

n10 = 10
newton_fail_10 = 0
disc_neg_10 = 0
total_10 = 0

for trial in range(500):
    T = random_tournament(n10)
    # Only check Omega_3 (3-cycles) for speed
    c3 = set()
    for i in range(n10):
        for j in range(n10):
            if j == i or not T[i][j]:
                continue
            for k in range(n10):
                if k == i or k == j or not T[j][k]:
                    continue
                if T[k][i]:
                    c = (i, j, k)
                    mi = c.index(min(c))
                    c3.add(c[mi:] + c[:mi])
    c3 = list(c3)
    m3 = len(c3)
    if m3 == 0:
        continue

    vsets3 = [frozenset(c) for c in c3]
    adj3 = {}
    for i in range(m3):
        adj3[i] = frozenset(j for j in range(m3) if j != i and vsets3[i] & vsets3[j])
    coeffs3 = independence_poly_dc(adj3)
    total_10 += 1

    newton_ok = newton_check(coeffs3)
    if not newton_ok:
        newton_fail_10 += 1

    if len(coeffs3) == 4:
        d = cubic_disc(*coeffs3)
        if d < 0:
            disc_neg_10 += 1

print(f"  n=10, 500 random tournaments, {total_10} with 3-cycles:")
print(f"  Omega_3 Newton fails: {newton_fail_10}")
print(f"  Omega_3 disc < 0 (degree 3): {disc_neg_10}")
print(f"  Max degree of I(Omega_3): {max(len(coeffs3)-1 for _ in [1])}")

# Degree distribution at n=10
deg_dist_10 = Counter()
for trial in range(200):
    T = random_tournament(n10)
    c3 = set()
    for i in range(n10):
        for j in range(n10):
            if j == i or not T[i][j]:
                continue
            for k in range(n10):
                if k == i or k == j or not T[j][k]:
                    continue
                if T[k][i]:
                    c = (i, j, k)
                    mi = c.index(min(c))
                    c3.add(c[mi:] + c[:mi])
    c3 = list(c3)
    m3 = len(c3)
    if m3 == 0:
        continue

    vsets3 = [frozenset(c) for c in c3]
    adj3 = {}
    for i in range(m3):
        adj3[i] = frozenset(j for j in range(m3) if j != i and vsets3[i] & vsets3[j])
    coeffs3 = independence_poly_dc(adj3)
    deg_dist_10[len(coeffs3)-1] += 1

print(f"\n  Degree distribution of I(Omega_3, x) at n=10: {dict(sorted(deg_dist_10.items()))}")
print(f"  At n=10, max independent set of 3-cycles can be 3 (covering 9 of 10 vertices)")
print(f"  or 4? No: 4 disjoint 3-cycles need 12 > 10 vertices. So alpha <= 3.")
