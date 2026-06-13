#!/usr/bin/env python3
"""
Fast n=8 OCF proof: check delta_H = delta_I using the A-clique formula.

Instead of computing I(Omega,2) from scratch (expensive cycle enumeration),
use the formula:
  delta_I = 2 * [sum_{gained C'} H(T'[V\V(C')]) - sum_{lost C} H(T[V\V(C)])]

where gained/lost cycles are odd cycles containing arc j->i / i->j respectively.
By induction (OCF proved at n<=7), H(complement) = I(Omega(complement), 2).

This only requires:
1. H(T) and H(T') via bitmask DP  [fast]
2. Finding odd cycles through 0->1 (lost) and 1->0 (gained)  [moderate]
3. Computing H(complement) for each such cycle  [fast, complement has <=5 vertices]

Instance: opus-2026-03-05-S4
"""

import time
import random
from itertools import combinations

N = 8
FULL = (1 << N) - 1
I_V, J_V = 0, 1


def ham_dp(T_flat, n, full_mask):
    """Bitmask DP for Hamiltonian path count."""
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, full_mask + 1):
        for last in range(n):
            key = (mask, last)
            if key not in dp:
                continue
            c = dp[key]
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if T_flat[last * n + nxt]:
                    nkey = (mask | (1 << nxt), nxt)
                    dp[nkey] = dp.get(nkey, 0) + c
    return sum(dp.get((full_mask, v), 0) for v in range(n))


def ham_count(T_flat):
    """H(T) for 8-vertex tournament stored as flat array."""
    # Use array-based DP for speed
    dp = [[0] * N for _ in range(1 << N)]
    for v in range(N):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << N):
        for last in range(N):
            c = dp[mask][last]
            if c == 0:
                continue
            for nxt in range(N):
                if mask & (1 << nxt):
                    continue
                if T_flat[last * N + nxt]:
                    dp[mask | (1 << nxt)][nxt] += c
    return sum(dp[FULL])


def ham_sub(T_flat, verts):
    """H(T[verts]) for a subset of vertices."""
    m = len(verts)
    if m <= 1:
        return 1
    full = (1 << m) - 1
    dp = [[0] * m for _ in range(1 << m)]
    for vi in range(m):
        dp[1 << vi][vi] = 1
    for mask in range(1, 1 << m):
        for li in range(m):
            c = dp[mask][li]
            if c == 0:
                continue
            for ni in range(m):
                if mask & (1 << ni):
                    continue
                if T_flat[verts[li] * N + verts[ni]]:
                    dp[mask | (1 << ni)][ni] += c
    return sum(dp[full])


def find_odd_cycles_through_arc(T_flat, src, dst):
    """Find all directed odd cycles containing the arc src->dst.
    Returns list of frozensets of vertex sets."""
    cycles = []
    others = [v for v in range(N) if v != src and v != dst]

    # Cycle of length L uses src->dst plus L-2 other vertices forming path dst->...->src
    for extra_count in range(1, len(others) + 1, 2):  # odd total length = extra_count + 2
        for combo in combinations(others, extra_count):
            # Find paths from dst to src through exactly combo vertices
            combo_list = list(combo)
            m = len(combo_list)
            # DP: paths from dst through subsets of combo ending at each vertex
            dp = {}
            # Start at dst
            for ci, c in enumerate(combo_list):
                if T_flat[dst * N + c]:
                    dp[(1 << ci, ci)] = 1

            full = (1 << m) - 1
            for mask in range(1, full + 1):
                for li in range(m):
                    if (mask, li) not in dp:
                        continue
                    cnt = dp[(mask, li)]
                    for ni in range(m):
                        if mask & (1 << ni):
                            continue
                        if T_flat[combo_list[li] * N + combo_list[ni]]:
                            key = (mask | (1 << ni), ni)
                            dp[key] = dp.get(key, 0) + cnt

            # Check paths ending at any vertex that connects back to src
            has_cycle = False
            for li in range(m):
                if dp.get((full, li), 0) > 0 and T_flat[combo_list[li] * N + src]:
                    has_cycle = True
                    break

            if has_cycle:
                cycles.append(frozenset([src, dst] + combo_list))

    return cycles


def compute_delta_I(T_flat, Tp_flat):
    """Compute I(Omega(T'),2) - I(Omega(T),2) using A-clique formula."""
    lost = find_odd_cycles_through_arc(T_flat, I_V, J_V)    # cycles with 0->1 in T
    gained = find_odd_cycles_through_arc(Tp_flat, J_V, I_V)  # cycles with 1->0 in T'

    delta = 0
    for cverts in gained:
        comp = [v for v in range(N) if v not in cverts]
        delta += 2 * ham_sub(Tp_flat, comp)
    for cverts in lost:
        comp = [v for v in range(N) if v not in cverts]
        delta -= 2 * ham_sub(T_flat, comp)

    return delta


def test_sample(sample_size=200):
    arc_pairs = [(a, b) for a in range(N) for b in range(a+1, N) if (a,b) != (0,1)]
    assert len(arc_pairs) == 27

    random.seed(42)
    print(f"Testing {sample_size} random n=8 tournaments...")

    t0 = time.time()
    passes = 0
    for trial in range(sample_size):
        mask = random.randint(0, (1 << 27) - 1)
        T = [0] * (N * N)
        T[0 * N + 1] = 1
        for bit, (a, b) in enumerate(arc_pairs):
            if mask & (1 << bit):
                T[a * N + b] = 1
            else:
                T[b * N + a] = 1

        Tp = T[:]
        Tp[0 * N + 1] = 0
        Tp[1 * N + 0] = 1

        ht = ham_count(T)
        htp = ham_count(Tp)
        delta_h = ht - htp

        delta_i = compute_delta_I(T, Tp)

        if delta_h == delta_i:
            passes += 1
        else:
            print(f"  FAIL trial {trial}: dH={delta_h}, dI={delta_i}, H(T)={ht}, H(T')={htp}")

    elapsed = time.time() - t0
    rate = sample_size / elapsed
    print(f"Result: {passes}/{sample_size} pass ({elapsed:.1f}s, {rate:.1f}/s)")
    est_total = (1 << 27) / rate
    print(f"Estimated full verification: {est_total:.0f}s = {est_total/3600:.1f}hr")
    return passes == sample_size


def full_verification():
    arc_pairs = [(a, b) for a in range(N) for b in range(a+1, N) if (a,b) != (0,1)]
    total = 1 << 27
    chunk_size = 1 << 17  # 131072

    print(f"=== n=8 FULL OCF VERIFICATION ===")
    print(f"Total configs: {total:,}")
    print(f"Chunk size: {chunk_size:,}")

    start = time.time()
    total_checked = 0
    total_fails = 0

    for chunk_start in range(0, total, chunk_size):
        chunk_end = min(chunk_start + chunk_size, total)
        t0 = time.time()
        fails = 0

        for mask in range(chunk_start, chunk_end):
            T = [0] * (N * N)
            T[0 * N + 1] = 1
            for bit, (a, b) in enumerate(arc_pairs):
                if mask & (1 << bit):
                    T[a * N + b] = 1
                else:
                    T[b * N + a] = 1

            Tp = T[:]
            Tp[0 * N + 1] = 0
            Tp[1 * N + 0] = 1

            ht = ham_count(T)
            htp = ham_count(Tp)
            delta_h = ht - htp
            delta_i = compute_delta_I(T, Tp)

            if delta_h != delta_i:
                fails += 1
                if total_fails + fails <= 5:
                    print(f"  FAIL at mask {mask}: dH={delta_h}, dI={delta_i}")

        elapsed = time.time() - t0
        total_checked += (chunk_end - chunk_start)
        total_fails += fails
        rate = (chunk_end - chunk_start) / elapsed if elapsed > 0 else 0
        eta = (total - total_checked) / rate if rate > 0 else 0
        pct = 100 * total_checked / total

        print(f"  [{pct:5.1f}%] {total_checked:>12,}/{total:,} | "
              f"{fails} fails | {rate:.0f}/s | ETA {eta:.0f}s ({eta/3600:.2f}hr)")

        if total_fails > 10:
            print("Too many failures, aborting.")
            return

    print(f"\n{'='*60}")
    print(f"TOTAL: {total_checked:,} configs, {total_fails} failures, {time.time()-start:.0f}s")
    if total_fails == 0:
        print(f"\n*** PROVED: OCF holds for ALL n=8 tournaments ***")


if __name__ == "__main__":
    import sys
    if "--full" in sys.argv:
        full_verification()
    else:
        test_sample(200)
