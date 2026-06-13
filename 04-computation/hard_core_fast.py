#!/usr/bin/env python3
"""
INV-028: Hard-core lattice gas — fast version (n<=6 only).
Tests properties of the conflict graph Omega(T) relevant to
the statistical mechanics interpretation H(T) = Z(Omega(T), lambda=2).

Instance: opus-2026-03-05-S4b
"""
import random
from itertools import combinations

def all_tournaments(n):
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    for bits in range(1 << len(pairs)):
        T = [[0]*n for _ in range(n)]
        for k, (i,j) in enumerate(pairs):
            if bits & (1 << k):
                T[i][j] = 1
            else:
                T[j][i] = 1
        yield T

def find_odd_cycles_fast(T):
    """Find all directed odd cycles using bitmask DP."""
    n = len(T)
    cycles = []
    for length in range(3, n+1, 2):
        for combo in combinations(range(n), length):
            start = combo[0]
            others = list(combo[1:])
            mo = len(others)
            full = (1 << mo) - 1
            dp = [0] * ((1 << mo) * mo)
            for k, o in enumerate(others):
                if T[start][o]:
                    dp[(1<<k)*mo + k] = 1
            for mask in range(1, full+1):
                for li in range(mo):
                    c = dp[mask*mo + li]
                    if c == 0: continue
                    for ni in range(mo):
                        if mask & (1<<ni): continue
                        if T[others[li]][others[ni]]:
                            dp[(mask|(1<<ni))*mo + ni] += c
            cnt = sum(dp[full*mo + li] for li in range(mo) if T[others[li]][start])
            if cnt > 0:
                cycles.append(frozenset(combo))
    return cycles

def ham_count(T):
    n = len(T)
    FULL = (1 << n) - 1
    dp = [[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v] = 1
    for mask in range(1, 1<<n):
        for last in range(n):
            c = dp[mask][last]
            if c == 0: continue
            for nxt in range(n):
                if mask & (1<<nxt): continue
                if T[last][nxt]:
                    dp[mask|(1<<nxt)][nxt] += c
    return sum(dp[FULL])

print("=" * 60)
print("HARD-CORE LATTICE GAS: Omega(T) PROPERTIES")
print("=" * 60)

for n in [4, 5, 6]:
    print(f"\n--- n = {n} (exhaustive) ---")

    max_deg_dist = {}
    density_sum = 0
    density_count = 0
    ncycles_dist = {}
    alpha_dist = {}  # independence number distribution
    total_T = 0
    lambda2_always_above_shearer = True

    for T in all_tournaments(n):
        total_T += 1
        cycles = find_odd_cycles_fast(T)
        nc = len(cycles)
        ncycles_dist[nc] = ncycles_dist.get(nc, 0) + 1

        if nc == 0:
            continue

        # Build adjacency
        adj = [set() for _ in range(nc)]
        for i in range(nc):
            for j in range(i+1, nc):
                if cycles[i] & cycles[j]:
                    adj[i].add(j)
                    adj[j].add(i)

        degrees = [len(adj[i]) for i in range(nc)]
        max_deg = max(degrees) if degrees else 0
        max_deg_dist[max_deg] = max_deg_dist.get(max_deg, 0) + 1

        edges = sum(len(a) for a in adj) // 2
        max_edges = nc * (nc-1) // 2
        if max_edges > 0:
            density_sum += edges / max_edges
            density_count += 1

        # Independence number (exact for small nc)
        alpha = 0
        for size in range(nc+1):
            found = False
            for subset in combinations(range(nc), size):
                indep = True
                for a in range(len(subset)):
                    for b in range(a+1, len(subset)):
                        if subset[b] in adj[subset[a]]:
                            indep = False
                            break
                    if not indep: break
                if indep:
                    alpha = size
                    found = True
            if not found:
                break
        alpha_dist[alpha] = alpha_dist.get(alpha, 0) + 1

        # Check if lambda=2 is above Shearer bound
        if max_deg >= 2:
            shearer = 1.0 / (max_deg - 1)
            if 2 <= shearer:
                lambda2_always_above_shearer = False

    print(f"  Total tournaments: {total_T}")
    print(f"  #cycles distribution: {dict(sorted(ncycles_dist.items()))}")
    print(f"  Max degree distribution: {dict(sorted(max_deg_dist.items()))}")
    if density_count > 0:
        print(f"  Average density: {density_sum/density_count:.3f}")
    print(f"  Independence number distribution: {dict(sorted(alpha_dist.items()))}")
    print(f"  lambda=2 always above Shearer: {lambda2_always_above_shearer}")

    # Is Omega always chordal?
    print(f"\n  Checking chordality and perfectness...")
    non_chordal = 0
    has_odd_hole = 0
    for T in all_tournaments(n):
        cycles = find_odd_cycles_fast(T)
        nc = len(cycles)
        if nc < 4: continue

        adj = [set() for _ in range(nc)]
        for i in range(nc):
            for j in range(i+1, nc):
                if cycles[i] & cycles[j]:
                    adj[i].add(j)
                    adj[j].add(i)

        # Check for induced 4-cycle (chordless cycle of length 4)
        for a, b, c, d in combinations(range(nc), 4):
            # Check if a-b-c-d-a is an induced 4-cycle
            edges_present = set()
            for x, y in [(a,b),(b,c),(c,d),(d,a),(a,c),(b,d)]:
                if y in adj[x]:
                    edges_present.add((x,y))
            if len(edges_present) == 4:
                # Check it's a 4-cycle
                if ({(a,b),(b,c),(c,d),(d,a)} <= edges_present and
                    (a,c) not in edges_present and (b,d) not in edges_present):
                    non_chordal += 1
                    break

        # Check for induced 5-cycle (odd hole)
        if nc >= 5:
            for combo in combinations(range(nc), 5):
                # Check all 12 possible 5-cycle orderings
                from itertools import permutations
                for perm in permutations(combo):
                    is_5cycle = True
                    for i in range(5):
                        if perm[(i+1)%5] not in adj[perm[i]]:
                            is_5cycle = False
                            break
                    if not is_5cycle: continue
                    # Check no chords
                    has_chord = False
                    for i in range(5):
                        for j in range(i+2, 5):
                            if j == (i+4)%5: continue
                            if perm[j] in adj[perm[i]]:
                                has_chord = True
                                break
                        if has_chord: break
                    if is_5cycle and not has_chord:
                        has_odd_hole += 1
                        break
                if has_odd_hole: break

    print(f"  Non-chordal Omega(T): {non_chordal}/{total_T}")
    print(f"  Omega(T) with odd hole (5-cycle): {has_odd_hole}/{total_T}")

print("\n" + "=" * 60)
print("CLUSTER EXPANSION CONVERGENCE BOUNDS")
print("=" * 60)
print("lambda=2 vs critical fugacity bounds:")
for Delta in range(1, 12):
    if Delta == 1:
        tree_bound = float('inf')
    else:
        tree_bound = ((Delta-1)**(Delta-1)) / (Delta**Delta)
    kp = 1.0 / (2.718 * (Delta + 1))
    print(f"  Delta={Delta:2d}: tree_bound={tree_bound:.6f}, KP={kp:.6f}, "
          f"lambda=2 {'BELOW' if 2 < tree_bound else 'ABOVE'} tree threshold")

print("\nKey finding: lambda=2 is above ALL convergence thresholds for Delta >= 2.")
print("This means the hard-core model on Omega(T) at fugacity 2 is in the")
print("non-perturbative / strong-coupling regime. The partition function")
print("Z(Omega, 2) = H(T) cannot be computed via cluster expansion.")
print("This is consistent with the difficulty of proving OCF — it's a")
print("non-perturbative identity requiring exact cancellations.")
