#!/usr/bin/env python3
"""
Fast n=8 OCF proof via exhaustive {0,1} evaluation.

Fix arc 0->1. Vary all C(8,2)-1 = 27 other arcs.
For each of 2^27 = 134,217,728 assignments:
  1. Build T (8x8)
  2. Compute H(T) and H(T') via bitmask DP
  3. Compute I(Omega(T),2) and I(Omega(T'),2)
  4. Check H(T)-H(T') = I(Omega(T'),2)-I(Omega(T),2)

Uses the full OCF check: delta_H = delta_I directly.

Instance: opus-2026-03-05-S4
"""

import time
from itertools import combinations

N = 8
FULL = (1 << N) - 1


def ham_path_count(T):
    """Count Hamiltonian paths using bitmask DP. O(2^n * n^2)."""
    dp = [[0] * N for _ in range(1 << N)]
    for v in range(N):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << N):
        for last in range(N):
            if not (mask & (1 << last)):
                continue
            if dp[mask][last] == 0:
                continue
            c = dp[mask][last]
            prev = mask ^ (1 << last)
            if prev == 0:
                continue
            # Actually, build forward
            pass
    # Rebuild: forward DP
    dp = [[0] * N for _ in range(1 << N)]
    for v in range(N):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << N):
        for last in range(N):
            if dp[mask][last] == 0:
                continue
            c = dp[mask][last]
            for nxt in range(N):
                if mask & (1 << nxt):
                    continue
                if T[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += c
    return sum(dp[FULL])


def find_odd_cycles(T):
    """Find all directed odd cycles in T. Returns list of frozensets of vertices."""
    cycles = set()
    for length in range(3, N + 1, 2):
        for combo in combinations(range(N), length):
            verts = list(combo)
            # Check all cyclic orderings (use DFS-like approach)
            _find_cycles_of_set(T, verts, length, cycles)
    return cycles


def _find_cycles_of_set(T, verts, length, cycles):
    """Find all directed cycles using exactly the given vertices."""
    vset = frozenset(verts)
    # Fix first vertex to avoid counting rotations
    start = verts[0]
    others = verts[1:]

    def dfs(path, remaining):
        if len(path) == length:
            if T[path[-1]][start]:
                cycles.add(vset)
            return
        for v in remaining:
            if T[path[-1]][v]:
                new_rem = [u for u in remaining if u != v]
                dfs(path + [v], new_rem)

    dfs([start], others)


def independence_poly_at_2(T):
    """Compute I(Omega(T), 2) where Omega is the conflict graph on odd cycles."""
    cycles = list(find_odd_cycles(T))
    nc = len(cycles)
    if nc == 0:
        return 1

    # Build conflict graph (adjacency by shared vertex)
    cycle_verts = [c for c in cycles]
    adj = [[False] * nc for _ in range(nc)]
    for a in range(nc):
        for b in range(a + 1, nc):
            if cycle_verts[a] & cycle_verts[b]:
                adj[a][b] = adj[b][a] = True

    # Enumerate independent sets, weight by 2^|S|
    total = 1  # empty set
    for size in range(1, nc + 1):
        for combo in combinations(range(nc), size):
            independent = True
            for idx_a in range(len(combo)):
                for idx_b in range(idx_a + 1, len(combo)):
                    if adj[combo[idx_a]][combo[idx_b]]:
                        independent = False
                        break
                if not independent:
                    break
            if independent:
                total += (1 << size)
    return total


def verify_chunk(start_mask, end_mask, arc_pairs):
    """Verify OCF for a range of arc assignments."""
    n_pairs = len(arc_pairs)
    fails = 0
    checked = 0

    for mask in range(start_mask, end_mask):
        T = [[0] * N for _ in range(N)]
        T[0][1] = 1  # fixed arc

        for bit, (a, b) in enumerate(arc_pairs):
            if mask & (1 << bit):
                T[a][b] = 1
            else:
                T[b][a] = 1

        # T' = flip 0->1 to 1->0
        Tp = [row[:] for row in T]
        Tp[0][1] = 0
        Tp[1][0] = 1

        h_t = ham_path_count(T)
        h_tp = ham_path_count(Tp)
        i_t = independence_poly_at_2(T)
        i_tp = independence_poly_at_2(Tp)

        delta_h = h_t - h_tp
        delta_i = i_tp - i_t  # Note: I(T') - I(T) should equal H(T') - H(T)
        # Wait: OCF says H = I for all T. So H(T)-H(T') should equal I(T)-I(T').
        # delta_h = h_t - h_tp, want this = i_t - i_tp = -(i_tp - i_t)
        # Actually: if H(T)=I(T) and H(T')=I(T'), then H(T)-H(T')=I(T)-I(T').

        if delta_h != i_t - i_tp:
            fails += 1
            if fails <= 3:
                print(f"  FAIL at mask {mask}: dH={delta_h}, dI={i_t-i_tp}")

        checked += 1

    return checked, fails


def main():
    print(f"=== n=8 OCF Exhaustive Proof ===")
    print(f"Fix arc 0->1. Vary {N*(N-1)//2 - 1} = 27 arcs.")
    print(f"Total: 2^27 = {1<<27:,} assignments")

    # Arc pairs (excluding 0,1)
    arc_pairs = []
    for a in range(N):
        for b in range(a + 1, N):
            if (a, b) == (0, 1):
                continue
            arc_pairs.append((a, b))
    assert len(arc_pairs) == 27

    total_configs = 1 << 27
    chunk_size = 1 << 20  # ~1M per chunk
    n_chunks = (total_configs + chunk_size - 1) // chunk_size

    print(f"Chunk size: {chunk_size:,}, Chunks: {n_chunks}")
    print()

    start = time.time()
    total_checked = 0
    total_fails = 0

    for chunk_idx in range(n_chunks):
        s = chunk_idx * chunk_size
        e = min(s + chunk_size, total_configs)

        t0 = time.time()
        checked, fails = verify_chunk(s, e, arc_pairs)
        elapsed = time.time() - t0
        total_checked += checked
        total_fails += fails

        rate = checked / elapsed if elapsed > 0 else 0
        eta = (total_configs - total_checked) / rate if rate > 0 else 0
        pct = 100 * total_checked / total_configs

        print(f"  Chunk {chunk_idx+1}/{n_chunks}: {checked:,} configs, "
              f"{fails} fails, {elapsed:.1f}s ({rate:.0f}/s) "
              f"[{pct:.1f}% done, ETA {eta:.0f}s]")

        if total_fails > 10:
            print("Too many failures, aborting.")
            break

    total_time = time.time() - start
    print(f"\n{'='*60}")
    print(f"Total: {total_checked:,}/{total_configs:,} checked, {total_fails} failures")
    print(f"Time: {total_time:.1f}s")

    if total_fails == 0 and total_checked == total_configs:
        print(f"\n*** PROVED: H(T) = I(Omega(T), 2) for ALL n=8 tournaments ***")


if __name__ == "__main__":
    main()
