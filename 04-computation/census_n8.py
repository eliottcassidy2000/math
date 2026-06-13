#!/usr/bin/env python3
"""
Optimized n=8 census: 2^21 = 2,097,152 tilings.

For n=8, full canonicalization (8! permutations) is infeasible for all tilings.
Strategy: use (sorted_score, H, score_multiset_of_neighbors) as approximate class key.

If two tilings share this key, they are LIKELY isomorphic.
Two tilings with different keys are definitely NOT isomorphic.

Instance: opus-2026-03-06-S3
"""
import sys
import time
from collections import defaultdict, Counter
import array

def count_ham_dp_fast(adj, n):
    """Optimized Hamiltonian path count using array-based DP."""
    full = (1 << n) - 1
    # dp[mask][v] stored as flat array
    size = (1 << n) * n
    dp = array.array('l', [0] * size)
    for v in range(n):
        dp[(1 << v) * n + v] = 1

    for mask in range(1, 1 << n):
        base = mask * n
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            c = dp[base + v]
            if c == 0:
                continue
            adj_v = adj[v]  # pre-fetched row
            nmask_base = mask
            for u in range(n):
                if nmask_base & (1 << u):
                    continue
                if adj_v[u]:
                    dp[(nmask_base | (1 << u)) * n + u] += c

    return sum(dp[full * n + v] for v in range(n))

def build_adj_rows(n, arc_bits, arcs):
    """Build adjacency as list of rows (tuples) for fast DP."""
    adj = [[0] * n for _ in range(n)]
    for i in range(n - 1):
        adj[i][i + 1] = 1
    for k in range(len(arcs)):
        i, j = arcs[k]
        if arc_bits & (1 << k):
            adj[j][i] = 1
        else:
            adj[i][j] = 1
    return adj

def score_seq(adj, n):
    return tuple(sorted(sum(adj[i]) for i in range(n)))

def refined_key(adj, n):
    """Compute a refined class key using neighbor degree structure."""
    out_deg = [sum(adj[i]) for i in range(n)]
    in_deg = [sum(adj[j][i] for j in range(n)) for i in range(n)]

    # For each vertex: (out_deg, in_deg, sorted out-neighbor out-degrees)
    vertex_sigs = []
    for v in range(n):
        out_nbrs = tuple(sorted(out_deg[u] for u in range(n) if adj[v][u]))
        in_nbrs = tuple(sorted(out_deg[u] for u in range(n) if adj[u][v]))
        vertex_sigs.append((out_deg[v], in_deg[v], out_nbrs, in_nbrs))

    return tuple(sorted(vertex_sigs))

def find_pos_arcs(n, arcs):
    return [k for k, (i, j) in enumerate(arcs) if j == 2 * i + 2]

def build_transpose_map(n, arcs):
    arc_idx = {a: k for k, a in enumerate(arcs)}
    trans_map = [0] * len(arcs)
    for k, (i, j) in enumerate(arcs):
        ti = n - 1 - j
        tj = n - 1 - i
        if ti > tj:
            ti, tj = tj, ti
        trans_map[k] = arc_idx.get((ti, tj), k)
    return trans_map

def is_grid_sym(bits, trans_map, m):
    for k in range(m):
        if ((bits >> k) & 1) != ((bits >> trans_map[k]) & 1):
            return False
    return True

def run_n8():
    n = 8
    arcs = [(i, j) for i in range(n) for j in range(i + 2, n)]
    m = len(arcs)  # 21
    total = 1 << m  # 2,097,152
    flip_mask = total - 1
    pos_indices = find_pos_arcs(n, arcs)
    trans_map = build_transpose_map(n, arcs)

    print(f"n=8: {m} arcs, {total} tilings, {len(pos_indices)} POS arcs")
    print(f"POS arcs: {[arcs[k] for k in pos_indices]}")
    print(f"Expected GS count: 2^{(n-1)**2//4} = {2**((n-1)**2//4)}")

    # Phase 1: Compute H, scores, refined key for all tilings
    print(f"\nPhase 1: Computing H and class keys...", flush=True)
    t0 = time.time()

    h_vals = array.array('l', [0] * total)
    keys = [None] * total  # refined class key
    gs_flags = bytearray(total)

    for bits in range(total):
        if bits % 200000 == 0 and bits > 0:
            elapsed = time.time() - t0
            rate = bits / elapsed
            eta = (total - bits) / rate
            print(f"  {bits}/{total} ({100*bits/total:.1f}%), "
                  f"{rate:.0f}/s, ETA {eta:.0f}s", flush=True)

        adj = build_adj_rows(n, bits, arcs)
        h = count_ham_dp_fast(adj, n)
        h_vals[bits] = h

        rkey = refined_key(adj, n)
        keys[bits] = (h, rkey)

        if is_grid_sym(bits, trans_map, m):
            gs_flags[bits] = 1

    elapsed = time.time() - t0
    print(f"  Phase 1 done in {elapsed:.1f}s")

    gs_total = sum(gs_flags)
    print(f"  Grid-symmetric tilings: {gs_total} (expected {2**((n-1)**2//4)})")

    # Phase 2: Group by class key
    print(f"\nPhase 2: Grouping and analysis...", flush=True)
    t1 = time.time()

    class_sizes = Counter()
    class_gs = Counter()
    class_h = {}

    for bits in range(total):
        k = keys[bits]
        class_sizes[k] += 1
        if gs_flags[bits]:
            class_gs[k] += 1
        class_h[k] = h_vals[bits]

    num_classes = len(class_sizes)
    print(f"  Approximate classes: {num_classes}")

    # Phase 3: Self-flip detection
    print(f"\nPhase 3: Self-flip detection...", flush=True)

    sf_count = 0
    blue_count = 0
    black_count = 0
    sf_per_class = Counter()
    blue_per_class = Counter()
    black_per_class = Counter()

    for bits in range(total):
        fb = bits ^ flip_mask
        if keys[bits] == keys[fb]:
            sf_count += 1
            sf_per_class[keys[bits]] += 1
            if gs_flags[bits]:
                blue_count += 1
                blue_per_class[keys[bits]] += 1
            else:
                black_count += 1
                black_per_class[keys[bits]] += 1

    elapsed = time.time() - t1
    print(f"  Phase 2+3 done in {elapsed:.1f}s")

    # Phase 4: POS analysis
    print(f"\nPhase 4: POS analysis...", flush=True)
    pos_global = Counter()
    for bits in range(total):
        pp = tuple((bits >> k) & 1 for k in pos_indices)
        pos_global[pp] += 1

    # ─── Report ───
    print(f"\n{'='*70}")
    print(f"RESULTS FOR n=8")
    print(f"{'='*70}")

    print(f"  Total tilings: {total}")
    print(f"  Approximate classes: {num_classes}")
    print(f"  Grid-symmetric tilings: {gs_total}")
    print(f"  Expected GS: {2**((n-1)**2//4)}")

    print(f"\n  Self-flip tilings: {sf_count}")
    print(f"  Blueself tilings: {blue_count}")
    print(f"  Blackself tilings: {black_count}")

    sf_classes = sum(1 for k in sf_per_class if sf_per_class[k] > 0)
    blue_classes = sum(1 for k in blue_per_class if blue_per_class[k] > 0)
    black_classes = sum(1 for k in black_per_class if black_per_class[k] > 0)
    mixed_classes = sum(1 for k in blue_per_class
                       if blue_per_class.get(k, 0) > 0 and black_per_class.get(k, 0) > 0)

    print(f"\n  Self-flip classes: {sf_classes}")
    print(f"  Blueself classes: {blue_classes}")
    print(f"  Blackself classes: {black_classes}")
    print(f"  Mixed (both blue+black) classes: {mixed_classes}")

    # POS distribution
    print(f"\n  POS orientation distribution ({len(pos_indices)} POS arcs):")
    for pp in sorted(pos_global.keys()):
        print(f"    POS={pp}: {pos_global[pp]} ({100*pos_global[pp]/total:.2f}%)")

    # H distribution
    h_dist = Counter(h_vals[bits] for bits in range(total))
    print(f"\n  H value distribution ({len(h_dist)} distinct values):")
    for h in sorted(h_dist.keys()):
        print(f"    H={h}: {h_dist[h]} tilings")

    # Class size distribution
    size_dist = Counter(class_sizes.values())
    print(f"\n  Class size distribution:")
    for s in sorted(size_dist.keys())[:30]:
        print(f"    Size {s}: {size_dist[s]} classes")
    if len(size_dist) > 30:
        print(f"    ... ({len(size_dist)} distinct sizes)")

    # Self-flip class details
    print(f"\n  Self-flip class details:")
    for k in sorted(sf_per_class.keys(), key=lambda k: class_h[k]):
        h = class_h[k]
        sz = class_sizes[k]
        gs = class_gs.get(k, 0)
        sf = sf_per_class[k]
        bl = blue_per_class.get(k, 0)
        bk = black_per_class.get(k, 0)
        print(f"    H={h:5d} size={sz:6d} GS={gs:4d} SF={sf:4d} "
              f"blue={bl:4d} black={bk:4d}")

    # GS per class (top 20 by GS count)
    print(f"\n  Top 20 classes by grid-symmetric count:")
    gs_sorted = sorted(class_gs.keys(), key=lambda k: -class_gs[k])[:20]
    for k in gs_sorted:
        h = class_h[k]
        sz = class_sizes[k]
        gs = class_gs[k]
        print(f"    H={h:5d} size={sz:6d} GS={gs:4d} ({100*gs/sz:.1f}%)")


if __name__ == '__main__':
    run_n8()
