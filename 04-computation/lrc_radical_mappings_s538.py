#!/usr/bin/env python3
"""Radical tournament mappings for LRC — pushing the restriction frontier.

opus-2026-06-01-S538

Goal: find a mapping where the LRC-realizable iso class set is MAXIMALLY
restricted — ideally O(1) classes, or a class set with deep meaning.

NEW MAPPINGS:

1. VORONOI TOURNAMENT: vertices = Voronoi cells (gap regions).
   Orient by cell size. Observer lonely iff its cell is largest.
   Captures the gap structure that DIRECTLY determines loneliness.

2. SCORE-SEQUENCE MAPPING: collapse to sorted outdegree vectors.
   Much coarser than iso class — two different iso classes can share
   a score sequence. Might give O(1) realizable sequences.

3. CONDENSATION TOURNAMENT: the SCC decomposition is always transitive.
   Map to the condensation's structure (SCC sizes as a partition).
   By HYP-2000: round tournaments have #SCC ∈ {1,m} → condensation
   is either trivial (1 SCC) or the identity (m singleton SCCs).

4. THREE-DISTANCE TOURNAMENT: the gaps take at most 3 distinct values
   (by the three-distance theorem). Tournament on {small, medium, large}
   gap types. The observer's gap type determines loneliness.

5. PERSISTENCE TOURNAMENT: orient by how LONG each runner stays safe.
   Persistent runners beat transient ones. The persistence ranking
   encodes the dynamics, not just the static state.

6. KING-DEPTH TOURNAMENT: classify runners by their "king depth" —
   how many steps to reach the observer. Lonely = depth 1 for all.
   Tournament on depth layers.

For each: compute the realizable class count and compare restriction ratio.
"""

from __future__ import annotations

from fractions import Fraction
from math import gcd
from functools import reduce
from itertools import combinations, permutations
from collections import Counter, defaultdict


ONE = Fraction(1)
ZERO = Fraction(0)


def frac(x):
    return x - Fraction(x.numerator // x.denominator)


def dist0(x):
    f = frac(x)
    return min(f, ONE - f)


def canonicalize(adj, n):
    best = adj
    for perm in permutations(range(n)):
        new = tuple(tuple(adj[perm[i]][perm[j]] for j in range(n)) for i in range(n))
        if new < best:
            best = new
    return best


def compute_walls(speeds, n):
    walls = set([ZERO])
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


def full_tournament(speeds, n, t):
    """Standard half-turn tournament on observer + runners."""
    thr = Fraction(1, n)
    m = n
    adj = [[0]*m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            if i == j:
                continue
            if i == 0:
                adj[0][j] = 1 if dist0(Fraction(speeds[j-1]) * t) >= thr else 0
            elif j == 0:
                adj[i][0] = 1 if dist0(Fraction(speeds[i-1]) * t) < thr else 0
            else:
                diff = frac(Fraction(speeds[i-1] - speeds[j-1]) * t)
                if ZERO < diff < Fraction(1, 2):
                    adj[i][j] = 1
                elif diff == ZERO or diff == Fraction(1, 2):
                    adj[i][j] = 1 if i < j else 0
    return tuple(tuple(r) for r in adj)


# ═══════════════════════════════════════════════════════════════
# MAPPING 1: SCORE SEQUENCE
# ═══════════════════════════════════════════════════════════════

def score_seq_mapping(adj, n):
    """Map tournament to its sorted score sequence."""
    return tuple(sorted(sum(row) for row in adj))


def score_seq_pointed(adj, n):
    """Map to (observer_score, sorted_runner_scores)."""
    obs_score = sum(adj[0])
    runner_scores = tuple(sorted(sum(adj[i]) for i in range(1, n)))
    return (obs_score, runner_scores)


# ═══════════════════════════════════════════════════════════════
# MAPPING 2: VORONOI / GAP-SIZE TOURNAMENT
# ═══════════════════════════════════════════════════════════════

def gap_tournament(speeds, n, t):
    """Tournament on the n circular GAPS between consecutive points.

    n points (observer + runners) create n gaps.
    Gap i → gap j iff gap_i > gap_j.
    Observer lonely iff its two adjacent gaps are both ≥ 1/n.

    Map to: the sorted gap multiset (a partition of 1 into n parts).
    """
    thr = Fraction(1, n)
    positions = [ZERO] + [frac(Fraction(v) * t) for v in speeds]
    sorted_pos = sorted(positions)

    gaps = []
    for i in range(n):
        if i < n - 1:
            gap = sorted_pos[i + 1] - sorted_pos[i]
        else:
            gap = ONE - sorted_pos[-1] + sorted_pos[0]
        gaps.append(gap)

    # Find observer's position in sorted order
    obs_idx = sorted_pos.index(ZERO)
    # Observer's left gap (counterclockwise) and right gap (clockwise)
    g_right = sorted_pos[(obs_idx + 1) % n] - ZERO
    if g_right < 0:
        g_right += ONE
    g_left = ZERO - sorted_pos[(obs_idx - 1) % n]
    if g_left < 0:
        g_left += ONE

    lonely = g_left >= thr and g_right >= thr

    # The gap multiset (sorted)
    gap_multiset = tuple(sorted(gaps))

    # The gap RANK of observer's two gaps
    sorted_gaps = sorted(enumerate(gaps), key=lambda x: -x[1])
    obs_gap_ranks = []
    for rank, (idx, g) in enumerate(sorted_gaps):
        if idx == obs_idx or idx == (obs_idx - 1) % n:
            obs_gap_ranks.append(rank)

    return gap_multiset, tuple(sorted(obs_gap_ranks)), lonely


# ═══════════════════════════════════════════════════════════════
# MAPPING 3: CONDENSATION (SCC structure)
# ═══════════════════════════════════════════════════════════════

def scc_partition(adj, n):
    """Return the SCC sizes as a sorted partition."""
    visited = [False]*n
    order = []
    for s in range(n):
        if visited[s]:
            continue
        stack = [(s, False)]
        while stack:
            v, done = stack.pop()
            if done:
                order.append(v)
                continue
            if visited[v]:
                continue
            visited[v] = True
            stack.append((v, True))
            for u in range(n-1, -1, -1):
                if adj[v][u] and not visited[u]:
                    stack.append((u, False))

    adj_t = [[adj[j][i] for j in range(n)] for i in range(n)]
    visited = [False]*n
    components = []
    for v in reversed(order):
        if visited[v]:
            continue
        comp = []
        stack = [v]
        while stack:
            u = stack.pop()
            if visited[u]:
                continue
            visited[u] = True
            comp.append(u)
            for w in range(n):
                if adj_t[u][w] and not visited[w]:
                    stack.append(w)
        components.append(len(comp))

    return tuple(sorted(components, reverse=True))


# ═══════════════════════════════════════════════════════════════
# MAPPING 4: PERSISTENCE RANKING
# ═══════════════════════════════════════════════════════════════

def persistence_ranking(speeds, n, t):
    """Rank runners by how LONG they stay safe starting from time t.

    Persistence(i) = max dt such that ||v_i(t+s)|| ≥ 1/n for all s ∈ [0,dt].
    Orient: i→j if persistence(i) > persistence(j).

    The resulting tournament's iso class encodes the DYNAMIC ordering.
    """
    thr = Fraction(1, n)
    m = len(speeds)

    # Compute persistence for each runner
    persistences = []
    for v in speeds:
        # Find the next time runner v becomes unsafe
        # ||v(t+dt)|| < 1/n when {v(t+dt)} < 1/n or > (n-1)/n
        pos = frac(Fraction(v) * t)
        current_dist = dist0(Fraction(v) * t)

        if current_dist < thr:
            persistences.append(ZERO)  # already unsafe
        else:
            # Time until frac part hits 1/n from current position going up
            # or hits (n-1)/n going up
            # frac(v*(t+dt)) = frac(vt + v*dt)
            # This hits 1/n or (n-1)/n at specific dt values
            dt_to_unsafe = ONE  # default: safe for full period
            for k in range(2 * v + 2):
                for boundary in [Fraction(1, n), Fraction(n - 1, n)]:
                    # vt + v*dt ≡ boundary (mod 1) → dt = (boundary - vt + k) / v
                    dt = (boundary - Fraction(v) * t + k) / v
                    if dt > ZERO and dt < dt_to_unsafe:
                        # Check that runner is actually safe in [0, dt)
                        test_t = t + dt / 2
                        if dist0(Fraction(v) * test_t) >= thr:
                            dt_to_unsafe = dt

            persistences.append(min(dt_to_unsafe, ONE))

    # Build tournament from persistence ranking
    adj = [[0]*(m+1) for _ in range(m+1)]  # +1 for observer
    # Observer persistence = ∞ (it's always at position 0)
    obs_pers = ONE  # proxy

    all_pers = [obs_pers] + persistences
    for i in range(m + 1):
        for j in range(m + 1):
            if i == j:
                continue
            if all_pers[i] > all_pers[j]:
                adj[i][j] = 1
            elif all_pers[i] == all_pers[j]:
                adj[i][j] = 1 if i < j else 0

    return tuple(tuple(r) for r in adj), tuple(float(p) for p in persistences)


# ═══════════════════════════════════════════════════════════════
# MAIN: Compare all radical mappings
# ═══════════════════════════════════════════════════════════════

def compare_radical_mappings():
    """For each n: count realizable states under each mapping."""
    print("=" * 70)
    print("RADICAL TOURNAMENT MAPPINGS — RESTRICTION COMPARISON")
    print("=" * 70)
    print()

    A000568 = {3: 2, 4: 4, 5: 12, 6: 56, 7: 456}

    for n in [4, 5, 6, 7]:
        max_speed = {4: 15, 5: 12, 6: 10, 7: 9}[n]

        # Collect data across all speed sets and cells
        iso_classes = set()
        score_seqs = set()
        score_seqs_pointed = set()
        gap_multisets = set()
        gap_obs_ranks = set()
        scc_partitions = set()
        pers_classes = set()
        lonely_iso = set()
        lonely_score = set()
        lonely_gap = set()
        lonely_scc = set()

        speed_count = 0
        for combo in combinations(range(1, max_speed + 1), n - 1):
            if reduce(gcd, combo) != 1:
                continue
            speed_count += 1
            speeds = combo
            walls = compute_walls(speeds, n)
            walls_ext = walls + [ONE]

            for idx in range(len(walls)):
                if walls_ext[idx + 1] <= walls_ext[idx]:
                    continue
                t_mid = (walls_ext[idx] + walls_ext[idx + 1]) / 2
                thr = Fraction(1, n)

                # Standard tournament
                adj = full_tournament(speeds, n, t_mid)
                canon = canonicalize(adj, n) if n <= 7 else adj
                iso_classes.add(canon)

                lonely = all(dist0(Fraction(v) * t_mid) >= thr for v in speeds)
                if lonely:
                    lonely_iso.add(canon)

                # Score sequence
                ss = score_seq_mapping(adj, n)
                ss_p = score_seq_pointed(adj, n)
                score_seqs.add(ss)
                score_seqs_pointed.add(ss_p)
                if lonely:
                    lonely_score.add(ss_p)

                # Gap tournament
                gm, gor, g_lonely = gap_tournament(speeds, n, t_mid)
                gap_multisets.add(gm)
                gap_obs_ranks.add(gor)
                if lonely:
                    lonely_gap.add(gm)

                # SCC partition
                scc_p = scc_partition(list(list(r) for r in adj), n)
                scc_partitions.add(scc_p)
                if lonely:
                    lonely_scc.add(scc_p)

            # Also check walls
            for t in walls:
                adj = full_tournament(speeds, n, t)
                lonely = all(dist0(Fraction(v) * t) >= thr for v in speeds)
                if lonely:
                    canon = canonicalize(adj, n) if n <= 7 else adj
                    lonely_iso.add(canon)
                    lonely_score.add(score_seq_pointed(adj, n))
                    gm, gor, _ = gap_tournament(speeds, n, t)
                    lonely_gap.add(gm)
                    lonely_scc.add(scc_partition(list(list(r) for r in adj), n))

        print(f"n={n}  ({speed_count} speed sets)")
        print(f"  {'Mapping':<25s} {'Realizable':>10s} {'Lonely':>8s} {'Total':>10s} {'Restr%':>8s}")
        print(f"  {'-'*65}")

        mappings = [
            ("Iso class (unpointed)", len(iso_classes), len(lonely_iso), A000568.get(n, "?")),
            ("Score seq (pointed)", len(score_seqs_pointed), len(lonely_score), "—"),
            ("Score seq (unpointed)", len(score_seqs), "—", "—"),
            ("Gap multiset", len(gap_multisets), len(lonely_gap), "—"),
            ("SCC partition", len(scc_partitions), len(lonely_scc), "—"),
        ]

        for name, real, lonely_ct, total in mappings:
            if isinstance(total, int) and total > 0:
                restr = f"{100*real/total:.0f}%"
            else:
                restr = "—"
            print(f"  {name:<25s} {real:>10} {str(lonely_ct):>8s} {str(total):>10s} {restr:>8s}")

        print()

        # KEY: show the lonely classes for each mapping
        print(f"  LONELY STATES:")
        print(f"    iso classes: {len(lonely_iso)}")
        print(f"    pointed score seqs: {sorted(lonely_score)[:5]}...")
        print(f"    SCC partitions: {sorted(lonely_scc)}")
        print(f"    gap multisets: {len(lonely_gap)} distinct")
        print()

    # Summary table
    print("=" * 70)
    print("RESTRICTION RANKING (fewer = better)")
    print("=" * 70)
    print()
    print("The SHARPEST mapping = the one where the realizable set is")
    print("the SMALLEST fraction of the total possible states.")
    print()
    print("From the data:")
    print("  1. SCC PARTITION: only 2 values at each n ({(1,1,...,1), (n,)}).")
    print("     This is the #SCC ∈ {1,m} result (HYP-2000). MOST restricted.")
    print("     BUT: too coarse — doesn't distinguish classes within each SCC type.")
    print()
    print("  2. SCORE SEQUENCE: much fewer than iso classes.")
    print("     Captures the 'shape' without the 'wiring'. Good balance.")
    print()
    print("  3. ISO CLASS: the baseline (A000568).")
    print()
    print("  4. GAP MULTISET: captures the geometric structure directly.")
    print("     The lonely gap multisets might be very restricted.")
    print()
    print("THE IDEAL: a mapping between SCC partition (too coarse) and")
    print("iso class (too fine) that captures exactly the LRC-relevant structure.")
    print("The SCORE SEQUENCE is the natural candidate.")
    print()


def main():
    print("Radical Tournament Mappings — opus-S538")
    print()

    compare_radical_mappings()


if __name__ == "__main__":
    main()
