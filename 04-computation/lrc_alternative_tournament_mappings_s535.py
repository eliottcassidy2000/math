#!/usr/bin/env python3
"""Alternative tournament mappings for LRC — searching for the most restricted class set.

opus-2026-06-01-S535

The current mapping (half-turn on runner positions) gives A000016 iso classes.
Can we find a mapping that:
(a) still captures the LRC lonely condition (observer = source)
(b) gives FEWER realizable iso classes
(c) makes the reachability question easier

MAPPINGS EXPLORED:

1. DISTANCE-RANK TOURNAMENT: orient by ||v_i t|| vs ||v_j t||
   (who is farther from observer). Lonely = observer ranks last.

2. GAP-ORDER TOURNAMENT: orient by the circular gap i→j vs j→i.
   Captures the gap structure directly.

3. CRT CHANNEL TOURNAMENT: collapse runners into CRT classes.
   Tournament on O(sqrt(n)) vertices instead of n.

4. SAFE/BLOCK PARTITION TOURNAMENT: tournament on {observer, safe-group,
   block-group}. The simplest possible reduction.

5. RAPIDITY TOURNAMENT: orient by hyperbolic distance (formal group).
   The lonely condition = observer has min rapidity.

6. APEX-REDUCED TOURNAMENT: condition on c=1 (apex aligned), then
   build tournament on the remaining runners.

7. DIFFERENCE-PAIR TOURNAMENT: vertices = speed differences.
   Tournament captures the pairwise relative dynamics.

8. LAYERED TOURNAMENT: partition runners by ||v_i t|| into layers
   [0,1/n), [1/n,2/n), ..., [(n-1)/2n, 1/2]. Tournament on layers.

For each: compute realizable iso classes, compare with A000016, check reachability.
"""

from __future__ import annotations

from fractions import Fraction
from math import gcd, comb
from functools import reduce
from itertools import combinations, permutations
from collections import Counter, defaultdict


ONE = Fraction(1)
ZERO = Fraction(0)


def frac(x: Fraction) -> Fraction:
    return x - Fraction(x.numerator // x.denominator)


def dist0(x: Fraction) -> Fraction:
    f = frac(x)
    return min(f, ONE - f)


def canonicalize(adj, n):
    """Canonicalize tournament under all vertex permutations."""
    best = adj
    for perm in permutations(range(n)):
        new = tuple(tuple(adj[perm[i]][perm[j]] for j in range(n)) for i in range(n))
        if new < best:
            best = new
    return best


def canonicalize_pointed(adj, n):
    """Canonicalize with vertex 0 fixed (observer)."""
    best = adj
    for perm in permutations(range(1, n)):
        full = (0,) + perm
        new = tuple(tuple(adj[full[i]][full[j]] for j in range(n)) for i in range(n))
        if new < best:
            best = new
    return best


def compute_walls(speeds, n):
    """All wall times for LRC."""
    walls = set()
    walls.add(ZERO)
    for v in speeds:
        for a in [1, n - 1]:
            for k in range(v):
                walls.add(Fraction(k * n + a, v * n))
        for k in range(v):
            walls.add(Fraction(k, v))
    for i, vi in enumerate(speeds):
        for j, vj in enumerate(speeds):
            if i >= j:
                continue
            diff = abs(vi - vj)
            for k in range(diff):
                walls.add(Fraction(k, diff))
                t = Fraction(2 * k + 1, 2 * diff)
                if ZERO <= t < ONE:
                    walls.add(t)
    return sorted(walls)


# ═══════════════════════════════════════════════════════════════
# MAPPING 1: DISTANCE-RANK TOURNAMENT
# ═══════════════════════════════════════════════════════════════

def distance_rank_tournament(speeds, n, t):
    """Tournament where i→j iff ||v_i t|| > ||v_j t|| (farther = winner).

    The observer is at distance 0 from itself.
    Lonely iff the observer is the MINIMUM (closest to itself = 0 distance).
    But in this tournament: observer → everyone iff all distances > 0 = always!
    Wait — the observer's distance is 0, so everyone beats the observer.
    Lonely iff all distances ≥ 1/n... this mapping doesn't capture loneliness
    as a SOURCE condition.

    REVISED: orient by SAFETY MARGIN: s_i = ||v_i t|| - 1/n.
    Positive = safe, negative = blocking.
    i→j iff s_i > s_j (more safe = winner).
    Observer has margin -1/n (always blocking itself? No, observer IS the reference.)

    BETTER REVISED: vertices = runners only (exclude observer).
    Orient i→j iff ||v_i t|| > ||v_j t||.
    This is a TOTAL ORDER on runners by distance.
    The tournament is ALWAYS TRANSITIVE (total orders give transitive tournaments).
    Iso class count = 1. Too trivial.

    EVEN BETTER: Include the observer as a vertex with "distance" = 1/n (the threshold).
    i→observer iff ||v_i t|| > 1/n (runner is safe → runner beats observer).
    observer→i iff ||v_i t|| < 1/n (runner is blocking → observer beats runner).
    Runner-runner: i→j iff ||v_i t|| > ||v_j t||.

    Now: observer lonely iff all runners beat observer iff all are safe.
    And the runner sub-tournament is TRANSITIVE (total order by distance).
    The observer can be anywhere in this total order.
    Observer at position k iff exactly k runners are closer to observer.
    Iso class determined by k = observer position in the ranking.
    n iso classes: k=0 (observer at top=lonely) to k=n-1 (all beat observer).

    Wait — observer lonely means all runners BEAT observer (i→obs for all i),
    so observer is at the BOTTOM. In this tournament, lonely = observer is SINK.

    The number of iso classes = n (one per observer position).
    This is MORE restricted than A000016! And the reachability question is:
    can the walk on {0,...,n-1} (observer position) reach 0 (or n-1)?
    """
    m = len(speeds)
    N = m + 1  # n vertices: observer + runners
    adj = [[0]*N for _ in range(N)]

    # Distances
    dists = [Fraction(1, n)]  # observer "distance" = threshold
    for v in speeds:
        dists.append(dist0(Fraction(v) * t))

    # Orient by distance (larger beats smaller)
    for i in range(N):
        for j in range(N):
            if i == j:
                continue
            if dists[i] > dists[j]:
                adj[i][j] = 1
            elif dists[i] == dists[j]:
                adj[i][j] = 1 if i < j else 0  # tie-break

    return tuple(tuple(r) for r in adj)


# ═══════════════════════════════════════════════════════════════
# MAPPING 2: GAP-ORDER TOURNAMENT (on gaps, not runners)
# ═══════════════════════════════════════════════════════════════

def gap_order_tournament(speeds, n, t):
    """Tournament on the n GAPS (between consecutive points on circle).

    The n points (observer + n-1 runners) create n circular gaps.
    Orient gap_i → gap_j iff gap_i > gap_j (larger gap wins).

    Lonely iff the observer's two adjacent gaps are both ≥ 1/n.
    In this tournament: lonely iff gaps adjacent to observer are "winners"
    (they beat many other gaps).

    Number of iso classes: determined by the sorted gap multiset.
    For the regular polygon: all gaps = 1/n → all ties → one class.
    For other configurations: gaps differ → various transitive-like tournaments.
    """
    # Compute positions and gaps
    positions = [ZERO] + [frac(Fraction(v) * t) for v in speeds]
    sorted_pos = sorted(set(positions))  # remove duplicates
    if len(sorted_pos) < n:
        return None  # degenerate (ties)

    # Compute gaps
    gaps = []
    for i in range(len(sorted_pos)):
        gap = sorted_pos[(i+1) % len(sorted_pos)] - sorted_pos[i]
        if gap < 0:
            gap += 1
        gaps.append(gap)

    N = len(gaps)
    adj = [[0]*N for _ in range(N)]
    for i in range(N):
        for j in range(N):
            if i == j:
                continue
            if gaps[i] > gaps[j]:
                adj[i][j] = 1
            elif gaps[i] == gaps[j]:
                adj[i][j] = 1 if i < j else 0

    return tuple(tuple(r) for r in adj)


# ═══════════════════════════════════════════════════════════════
# MAPPING 3: LAYERED TOURNAMENT (partition by distance shell)
# ═══════════════════════════════════════════════════════════════

def layered_tournament(speeds, n, t, num_layers=3):
    """Tournament on distance layers.

    Partition runners by distance from observer:
    Layer 0: ||v_i t|| < 1/n (blocking)
    Layer 1: 1/n ≤ ||v_i t|| < 2/n (borderline safe)
    Layer 2: ||v_i t|| ≥ 2/n (deeply safe)

    Build tournament on {observer, layer0, layer1, layer2}.
    Observer → layer_k iff k ≥ 1 (safe layers).
    Layer_i → layer_j by group size comparison or some other rule.

    Lonely iff layer0 is EMPTY (all runners in layer1 or layer2).
    """
    thr = Fraction(1, n)
    layers = [[] for _ in range(num_layers)]

    for v in speeds:
        d = dist0(Fraction(v) * t)
        layer = 0
        for L in range(num_layers):
            if d >= Fraction(L, n) and (L == num_layers - 1 or d < Fraction(L + 1, n)):
                layer = L
                break
        layers[layer].append(v)

    # Build tournament on num_layers+1 vertices (observer + layers)
    N = num_layers + 1
    adj = [[0]*N for _ in range(N)]

    # Observer = vertex 0
    # Observer → layer_k iff k ≥ 1
    for k in range(1, num_layers):
        adj[0][k + 1] = 1 if len(layers[k]) > 0 else 0

    # Layer_i → layer_j: i beats j if i has more runners (size dominance)
    for i in range(num_layers):
        for j in range(num_layers):
            if i == j:
                continue
            if len(layers[i]) > len(layers[j]):
                adj[i + 1][j + 1] = 1
            elif len(layers[i]) == len(layers[j]):
                adj[i + 1][j + 1] = 1 if i > j else 0

    return tuple(tuple(r) for r in adj), tuple(len(l) for l in layers)


# ═══════════════════════════════════════════════════════════════
# MAPPING 4: APEX-REDUCED TOURNAMENT
# ═══════════════════════════════════════════════════════════════

def apex_reduced_tournament(speeds, n, t):
    """Condition on apex aligned (c=1), then build tournament on remaining runners.

    The apex connects observer to runner_1 (slowest).
    If apex is aligned (runner_1 safe): build the half-turn tournament
    on runners 2,...,n-1 only (n-2 vertices).
    If apex is anti-aligned: return None (not in the conditioned space).

    The number of iso classes of this reduced tournament is smaller
    because we've removed one runner and conditioned on the apex.
    """
    thr = Fraction(1, n)

    # Check apex: is runner_1 (speed = min(speeds)) safe?
    if dist0(Fraction(speeds[0]) * t) < thr:
        return None  # apex not aligned

    # Build tournament on remaining runners (excluding runner_1)
    remaining = speeds[1:]
    m = len(remaining)
    adj = [[0]*m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            if i == j:
                continue
            diff = frac(Fraction(remaining[i] - remaining[j]) * t)
            if ZERO < diff < Fraction(1, 2):
                adj[i][j] = 1

    return tuple(tuple(r) for r in adj)


# ═══════════════════════════════════════════════════════════════
# MAPPING 5: MODULAR RESIDUE TOURNAMENT
# ═══════════════════════════════════════════════════════════════

def modular_tournament(speeds, n, t, modulus=None):
    """Tournament on residue classes mod some modulus.

    For n=6, modulus=3: classes {1,4}, {2,5}, {3}.
    For n=14, modulus=7: classes {1,8},...,{6,13},{7}.

    Arc between classes: class A → class B if MAJORITY of runner arcs
    from A-runners to B-runners go A→B.
    """
    if modulus is None:
        # Find smallest non-trivial factor of n
        for p in range(2, n):
            if n % p == 0:
                modulus = n // p
                break
        else:
            modulus = n

    thr = Fraction(1, n)

    # Group runners by residue mod modulus
    classes = defaultdict(list)
    for v in speeds:
        classes[v % modulus].append(v)

    class_list = sorted(classes.keys())
    m = len(class_list) + 1  # +1 for observer

    adj = [[0]*m for _ in range(m)]

    # Observer = vertex 0
    for ci, c in enumerate(class_list):
        runners = classes[c]
        # Observer → class_c iff ALL runners in c are safe
        all_safe = all(dist0(Fraction(v) * t) >= thr for v in runners)
        if all_safe:
            adj[0][ci + 1] = 1
        else:
            adj[ci + 1][0] = 1

    # Class vs class: majority vote
    for ci, c1 in enumerate(class_list):
        for cj, c2 in enumerate(class_list):
            if ci == cj:
                continue
            wins = 0
            total = 0
            for v1 in classes[c1]:
                for v2 in classes[c2]:
                    total += 1
                    diff = frac(Fraction(v1 - v2) * t)
                    if ZERO < diff < Fraction(1, 2):
                        wins += 1
            if wins > total / 2:
                adj[ci + 1][cj + 1] = 1
            elif wins < total / 2:
                adj[cj + 1][ci + 1] = 1
            else:
                adj[ci + 1][cj + 1] = 1 if ci < cj else 0

    return tuple(tuple(r) for r in adj), m


# ═══════════════════════════════════════════════════════════════
# MAIN: Compare all mappings
# ═══════════════════════════════════════════════════════════════

def compare_mappings(n_values=[4, 5, 6]):
    """For each n and mapping: count realizable iso classes."""
    print("=" * 70)
    print("COMPARING TOURNAMENT MAPPINGS FOR LRC")
    print("=" * 70)
    print()

    for n in n_values:
        speeds_initial = tuple(range(1, n))
        max_speed = {4: 15, 5: 10, 6: 9}[n]

        print(f"{'='*60}")
        print(f"n={n}")
        print(f"{'='*60}")
        print()

        # Collect all speed sets
        speed_sets = []
        for combo in combinations(range(1, max_speed + 1), n - 1):
            if reduce(gcd, combo) == 1:
                speed_sets.append(combo)

        print(f"Speed sets tested: {len(speed_sets)}")
        print()

        # For each mapping, count distinct classes
        mappings = {}

        # Standard half-turn (baseline)
        classes_standard = set()
        classes_standard_pointed = set()
        for speeds in speed_sets:
            walls = compute_walls(speeds, n)
            walls_ext = walls + [ONE]
            for idx in range(len(walls)):
                if walls_ext[idx + 1] <= walls_ext[idx]:
                    continue
                t_mid = (walls_ext[idx] + walls_ext[idx + 1]) / 2
                thr = Fraction(1, n)
                m = n
                adj = [[0]*m for _ in range(m)]
                for i in range(m):
                    for j in range(m):
                        if i == j: continue
                        if i == 0:
                            adj[0][j] = 1 if dist0(Fraction(speeds[j-1]) * t_mid) >= thr else 0
                        elif j == 0:
                            adj[i][0] = 1 if dist0(Fraction(speeds[i-1]) * t_mid) < thr else 0
                        else:
                            diff = frac(Fraction(speeds[i-1] - speeds[j-1]) * t_mid)
                            adj[i][j] = 1 if ZERO < diff < Fraction(1, 2) else 0
                adj_t = tuple(tuple(r) for r in adj)
                if n <= 6:
                    classes_standard.add(canonicalize(adj_t, m))
                    classes_standard_pointed.add(canonicalize_pointed(adj_t, m))

        print(f"  STANDARD (half-turn, n={n} vertices):")
        if n <= 6:
            print(f"    unpointed iso classes: {len(classes_standard)}")
            print(f"    pointed iso classes: {len(classes_standard_pointed)}")
        print()

        # Distance-rank tournament
        classes_dist = set()
        for speeds in speed_sets[:50]:  # sample
            walls = compute_walls(speeds, n)
            walls_ext = walls + [ONE]
            for idx in range(len(walls)):
                if walls_ext[idx + 1] <= walls_ext[idx]:
                    continue
                t_mid = (walls_ext[idx] + walls_ext[idx + 1]) / 2
                adj = distance_rank_tournament(speeds, n, t_mid)
                if n <= 6:
                    classes_dist.add(canonicalize_pointed(adj, n))

        print(f"  DISTANCE-RANK (n={n} vertices, 50 speed sets):")
        if n <= 6:
            print(f"    pointed iso classes: {len(classes_dist)}")
        print()

        # Apex-reduced tournament (conditioned on c=1)
        classes_apex = set()
        apex_count = 0
        for speeds in speed_sets[:50]:
            walls = compute_walls(speeds, n)
            walls_ext = walls + [ONE]
            for idx in range(len(walls)):
                if walls_ext[idx + 1] <= walls_ext[idx]:
                    continue
                t_mid = (walls_ext[idx] + walls_ext[idx + 1]) / 2
                adj = apex_reduced_tournament(speeds, n, t_mid)
                if adj is not None:
                    apex_count += 1
                    if n - 1 <= 6:
                        classes_apex.add(canonicalize(adj, n - 2))

        print(f"  APEX-REDUCED ({n-2} runners, conditioned on c=1, 50 sets):")
        if n - 1 <= 6:
            print(f"    unpointed iso classes: {len(classes_apex)}")
            print(f"    apex-aligned cells: {apex_count}")
        print()

        # Modular tournament
        for modulus_name, modulus in [("mod-2", 2), ("mod-3", 3)]:
            if n % modulus != 0 and modulus != 2:
                continue
            classes_mod = set()
            mod_n_vertices = 0
            for speeds in speed_sets[:50]:
                walls = compute_walls(speeds, n)
                walls_ext = walls + [ONE]
                for idx in range(len(walls)):
                    if walls_ext[idx + 1] <= walls_ext[idx]:
                        continue
                    t_mid = (walls_ext[idx] + walls_ext[idx + 1]) / 2
                    adj, mod_nv = modular_tournament(speeds, n, t_mid, modulus)
                    mod_n_vertices = mod_nv
                    if mod_nv <= 6:
                        classes_mod.add(canonicalize_pointed(adj, mod_nv))

            print(f"  MODULAR ({modulus_name}, {mod_n_vertices} vertices, 50 sets):")
            if mod_n_vertices <= 6:
                print(f"    pointed iso classes: {len(classes_mod)}")
            print()

        # Layered tournament
        classes_layer = set()
        for speeds in speed_sets[:50]:
            walls = compute_walls(speeds, n)
            walls_ext = walls + [ONE]
            for idx in range(len(walls)):
                if walls_ext[idx + 1] <= walls_ext[idx]:
                    continue
                t_mid = (walls_ext[idx] + walls_ext[idx + 1]) / 2
                adj, layer_sizes = layered_tournament(speeds, n, t_mid, 3)
                classes_layer.add((canonicalize_pointed(adj, 4), layer_sizes))

        print(f"  LAYERED (3 layers + observer = 4 vertices, 50 sets):")
        print(f"    (class, layer_sizes) pairs: {len(classes_layer)}")
        print()

    # Summary
    print("=" * 60)
    print("SUMMARY: Iso class counts by mapping")
    print("=" * 60)
    print()
    print("The MOST RESTRICTED mapping = fewest iso classes.")
    print("The MOST USEFUL mapping = captures LRC condition + few classes.")
    print()
    print("RANKINGS (fewer = better for proof):")
    print("  1. Distance-rank: n pointed classes (trivially transitive)")
    print("     BUT: runner structure lost. Walk = walk on {0,...,n-1}.")
    print("  2. Layered: O(n) classes (partition by distance shell)")
    print("     Captures coarse structure. Lonely = layer0 empty.")
    print("  3. Modular: O(A000568(sqrt(n))) classes")
    print("     CRT structure preserved. Lonely = observer beats all classes.")
    print("  4. Apex-reduced: A000568(n-2) classes (one fewer runner)")
    print("     Conditioned on apex. Captures the INTERIOR structure.")
    print("  5. Standard: A000016(n-1) classes (current best)")
    print()
    print("THE KEY INSIGHT: the distance-rank mapping collapses to n classes")
    print("because the tournament is ALWAYS TRANSITIVE (total order by distance).")
    print("This is the SIMPLEST possible mapping. The walk on {0,...,n-1}")
    print("(observer position in the distance ranking) IS the observer outdegree walk.")
    print("LRC = the walk reaches position n-1 (observer at the bottom = all safe).")
    print()
    print("But the distance-rank walk is NOT just a random walk on {0,...,n-1}.")
    print("It's constrained by the runner dynamics. The CONSTRAINTS come from")
    print("the runner sub-tournament — which the distance-rank mapping LOSES.")
    print()
    print("The IDEAL mapping preserves enough runner structure to constrain the walk")
    print("while having few enough classes for a tractable reachability analysis.")
    print("The MODULAR mapping at n=14 with modulus 7 gives O(A000568(8))=O(6880)")
    print("classes — but the CRT restriction makes most unrealizable.")
    print("The ACTUAL count of realizable CRT-modular classes at n=14 may be << 6880.")
    print()


def main():
    print("Alternative Tournament Mappings for LRC — opus-S535")
    print()

    compare_mappings()


if __name__ == "__main__":
    main()
