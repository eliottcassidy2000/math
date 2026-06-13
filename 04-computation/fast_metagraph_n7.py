"""
fast_metagraph_n7.py — opus-2026-04-04-S8

Build the n=7 metagraph (456 classes) using a FAST invariant.

The (scores, c3, H) invariant gives 242 groups (not enough — 456 classes).
We need to add MORE invariants to split the remaining groups.

Strategy:
  1. Compute (scores, c3, c5, H) — adding 5-cycle count
  2. If still under-distinguishing, add sorted list of vertex-deleted H values
  3. Verify against known class count A000568(7) = 456
"""
import numpy as np
from collections import defaultdict
import time

def make_tiles(n):
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            if x - y >= 2:
                tiles.append((x, y))
    return tiles

def tiling_to_adj(n, tiles, bits):
    T = np.zeros((n, n), dtype=np.int8)
    for k in range(n-1):
        T[k+1][k] = 1
    for i, (x, y) in enumerate(tiles):
        a, b = x-1, y-1
        if bits & (1 << i):
            T[b][a] = 1
        else:
            T[a][b] = 1
    return T

def count_hp(T, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    full = (1 << n) - 1
    total = 0
    for mask in range(1, 1 << n):
        for v in range(n):
            cnt = dp.get((mask, v), 0)
            if cnt == 0: continue
            if mask == full: total += cnt; continue
            for u in range(n):
                if mask & (1 << u): continue
                if T[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + cnt
    return total

def count_cycles(T, n, length):
    """Count directed cycles of given length via trace of A^k / k."""
    A = T.astype(np.float64)
    Ak = np.eye(n)
    for _ in range(length):
        Ak = Ak @ A
    return int(round(np.trace(Ak) / length))

def rich_invariant(T, n):
    """Compute a rich invariant for tournament classification."""
    scores = tuple(sorted(int(T[i].sum()) for i in range(n)))
    c3 = count_cycles(T, n, 3)
    h = count_hp(T, n)

    # Add vertex-deleted H values for finer discrimination
    del_h = []
    for v in range(n):
        # Delete vertex v
        remaining = [u for u in range(n) if u != v]
        T_del = np.zeros((n-1, n-1), dtype=np.int8)
        for i, u in enumerate(remaining):
            for j, w in enumerate(remaining):
                T_del[i][j] = T[u][w]
        del_h.append(count_hp(T_del, n-1))
    del_h_sorted = tuple(sorted(del_h))

    return (scores, c3, h, del_h_sorted)

def build_fast_n7():
    n = 7
    tiles = make_tiles(n)
    m = len(tiles)

    print(f"Building n={n} metagraph with rich invariant")
    print(f"m={m} tiles, 2^m={1<<m} tilings")
    print("="*60)

    t0 = time.time()
    groups = defaultdict(list)

    for bits in range(1 << m):
        T = tiling_to_adj(n, tiles, bits)
        inv = rich_invariant(T, n)
        groups[inv].append(bits)

        if bits % 2000 == 0 and bits > 0:
            elapsed = time.time() - t0
            rate = bits / elapsed
            eta = ((1 << m) - bits) / rate
            print(f"  {bits}/{1<<m} ({len(groups)} groups, {elapsed:.0f}s, "
                  f"ETA {eta:.0f}s)", flush=True)

    elapsed = time.time() - t0
    nc = len(groups)
    print(f"\n  Groups: {nc} (expected 456)")
    print(f"  Time: {elapsed:.0f}s")

    expected = 456
    if nc == expected:
        print(f"  *** PERFECT MATCH! Rich invariant distinguishes all {expected} classes ***")
    elif nc < expected:
        print(f"  Under-distinguishing: {nc} < {expected}")
        # Check how many groups need splitting
        multi = sum(1 for g in groups.values() if len(g) > 1)
        max_size = max(len(g) for g in groups.values())
        # Check tiling-fibration for violations
        violations = 0
        for inv, members in groups.items():
            h = inv[2]
            size = len(members)
            aut = h / size
            if abs(aut - round(aut)) > 0.001:
                violations += 1
        print(f"  Multi-member groups: {multi}")
        print(f"  Max group size: {max_size}")
        print(f"  Tiling-fibration violations: {violations}")
    else:
        print(f"  Over-distinguishing: {nc} > {expected}")

    # H-value statistics
    h_counts = defaultdict(int)
    for inv in groups:
        h_counts[inv[2]] += 1
    print(f"\n  Distinct H values: {len(h_counts)}")
    print(f"  H range: [{min(h_counts)}, {max(h_counts)}]")

    # Build adjacency for the metagraph
    print(f"\n  Building adjacency...")
    t1 = time.time()

    # Create class map
    class_map = {}
    class_list = list(groups.keys())
    for cid, inv in enumerate(class_list):
        for bits in groups[inv]:
            class_map[bits] = cid

    # Adjacency
    adj = defaultdict(set)
    for bits in range(1 << m):
        c1 = class_map[bits]
        for tile in range(m):
            c2 = class_map[bits ^ (1 << tile)]
            if c1 != c2:
                adj[c1].add(c2)

    num_edges = sum(len(v) for v in adj.values()) // 2
    degrees = [len(adj[c]) for c in range(nc)]

    t2 = time.time()
    print(f"  Edges: {num_edges}")
    print(f"  Degree: min={min(degrees)}, max={max(degrees)}, mean={np.mean(degrees):.1f}")
    print(f"  Adjacency time: {t2-t1:.0f}s")

    # Self-loop statistics
    self_loops = defaultdict(int)
    for bits in range(1 << m):
        c1 = class_map[bits]
        for tile in range(m):
            if class_map[bits ^ (1 << tile)] == c1:
                self_loops[c1] += 1

    total_sl = sum(self_loops.values())
    total_flips = (1 << m) * m
    sl_rate = total_sl / total_flips
    print(f"\n  Self-loop rate: {total_sl}/{total_flips} = {100*sl_rate:.1f}%")

    # H-gradient across edges
    dh_dist = defaultdict(int)
    for c1 in range(nc):
        for c2 in adj[c1]:
            if c2 > c1:
                h1 = class_list[c1][2]
                h2 = class_list[c2][2]
                dh_dist[abs(h1-h2)] += 1
    print(f"\n  ΔH distribution:")
    for dh in sorted(dh_dist.keys())[:15]:
        print(f"    ΔH={dh}: {dh_dist[dh]} edges")

    return groups, adj, class_map

if __name__ == "__main__":
    build_fast_n7()
