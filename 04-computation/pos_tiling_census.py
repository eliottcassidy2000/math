#!/usr/bin/env python3
"""
POS (Points of Symmetry) tiling census — OPTIMIZED version.

For each n, compute:
1. POS arcs and sigma structure
2. Tiling count per isomorphism class
3. Blueself/blackself classification
4. Grid-symmetric tiling counts
5. POS orientation patterns per class

Optimizations:
- Flat arrays instead of lists-of-lists
- Nauty-style hash refinement for canonicalization
- Bitmask DP with flat arrays
- Batch processing

Instance: opus-2026-03-06-S3
"""
import sys
import time
from collections import defaultdict, Counter
from itertools import permutations

# ─── Core: Fast Hamiltonian path count ───

def count_ham_dp(adj_flat, n):
    """Count Hamiltonian paths by bitmask DP. adj_flat is n*n flat array."""
    full = (1 << n) - 1
    # Use a flat dict for sparse DP
    dp = {}
    for v in range(n):
        dp[(1 << v) * n + v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            key = mask * n + v
            c = dp.get(key, 0)
            if c == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj_flat[v * n + u]:
                    nkey = (mask | (1 << u)) * n + u
                    dp[nkey] = dp.get(nkey, 0) + c
    return sum(dp.get(full * n + v, 0) for v in range(n))

def build_adj_flat(n, arc_bits, arcs):
    """Build flat adjacency array."""
    adj = bytearray(n * n)
    for i in range(n - 1):
        adj[i * n + (i + 1)] = 1
    for k in range(len(arcs)):
        i, j = arcs[k]
        if arc_bits & (1 << k):
            adj[j * n + i] = 1
        else:
            adj[i * n + j] = 1
    return adj

# ─── Canonicalization ───

def score_sequence(adj_flat, n):
    """Sorted score sequence."""
    scores = [0] * n
    for i in range(n):
        s = 0
        for j in range(n):
            s += adj_flat[i * n + j]
        scores[i] = s
    scores.sort()
    return tuple(scores)

def canon_hash(adj_flat, n):
    """Fast hash-based canonical form.
    For n <= 7: full permutation search (fast enough).
    For n = 8: use degree refinement + partial search.
    """
    if n <= 7:
        return _canon_full(adj_flat, n)
    else:
        return _canon_refined(adj_flat, n)

def _canon_full(adj_flat, n):
    """Full permutation canonical form for small n."""
    best = None
    for perm in permutations(range(n)):
        s = tuple(adj_flat[perm[i] * n + perm[j]] for i in range(n) for j in range(n))
        if best is None or s < best:
            best = s
    return best

def _canon_refined(adj_flat, n):
    """Refined canonical form for n=8 using degree-based partition.
    Groups vertices by (out-degree, in-degree, neighbor-degree-multiset).
    Then tries permutations respecting the partition.
    """
    # Compute out-degree and in-degree
    out_deg = [0] * n
    in_deg = [0] * n
    for i in range(n):
        for j in range(n):
            if adj_flat[i * n + j]:
                out_deg[i] += 1
                in_deg[j] += 1

    # Vertex signature: (out_deg, in_deg, sorted neighbor out-degrees)
    sigs = []
    for v in range(n):
        out_neighbors = tuple(sorted(out_deg[u] for u in range(n) if adj_flat[v * n + u]))
        in_neighbors = tuple(sorted(out_deg[u] for u in range(n) if adj_flat[u * n + v]))
        sigs.append((out_deg[v], in_deg[v], out_neighbors, in_neighbors))

    # Group by signature
    groups = defaultdict(list)
    for v in range(n):
        groups[sigs[v]].append(v)

    # Sort groups by signature for deterministic ordering
    sorted_groups = sorted(groups.values(), key=lambda g: sigs[g[0]])

    # If all groups have size 1, canonical form is determined
    if all(len(g) == 1 for g in sorted_groups):
        perm = [g[0] for g in sorted_groups]
        return tuple(adj_flat[perm[i] * n + perm[j]] for i in range(n) for j in range(n))

    # Generate permutations respecting partition
    from itertools import product as iproduct
    def gen_perms(groups):
        if not groups:
            yield []
            return
        g = groups[0]
        for p in permutations(g):
            for rest in gen_perms(groups[1:]):
                yield list(p) + rest

    best = None
    count = 0
    for perm in gen_perms(sorted_groups):
        s = tuple(adj_flat[perm[i] * n + perm[j]] for i in range(n) for j in range(n))
        if best is None or s < best:
            best = s
        count += 1
        if count > 100000:  # safety limit
            break
    return best

# ─── POS and sigma structure ───

def find_pos_arcs(n, arcs):
    """POS arcs: tile (r,c) with r=c, i.e., arc (i,j) with j=2i+2."""
    return [k for k, (i, j) in enumerate(arcs) if j == 2 * i + 2]

def find_sigma_pairs(n, arcs):
    """Sigma-paired arcs (not POS). Pin-grid sigma: arc (i,j) -> (j-i-2, j)."""
    arc_set = {a: k for k, a in enumerate(arcs)}
    pos_set = set(find_pos_arcs(n, arcs))
    pairs = []
    seen = set()
    for k, (i, j) in enumerate(arcs):
        if k in seen or k in pos_set:
            continue
        sigma_arc = (j - i - 2, j)
        if sigma_arc in arc_set:
            sk = arc_set[sigma_arc]
            pairs.append((k, sk))
            seen.add(k)
            seen.add(sk)
    return pairs

# ─── Grid transpose ───

def build_transpose_map(n, arcs):
    """Grid transpose map on 0-indexed arcs.
    In 1-indexed tiles: (a,b) -> (n+1-b, n+1-a).
    0-indexed arc (i,j) maps to 1-indexed tile (n-i, n-j).
    Transpose tile: (n+1-(n-j), n+1-(n-i)) = (j+1, i+1).
    Back to 0-indexed arc: (n-(j+1), n-(i+1)) = (n-j-1, n-i-1).
    But we need j' > i' and j' >= i'+2.
    """
    arc_idx = {a: k for k, a in enumerate(arcs)}
    trans_map = [0] * len(arcs)
    for k, (i, j) in enumerate(arcs):
        # Transform
        ti = n - 1 - j
        tj = n - 1 - i
        if ti > tj:
            ti, tj = tj, ti
        trans_arc = (ti, tj)
        if trans_arc in arc_idx:
            trans_map[k] = arc_idx[trans_arc]
        else:
            trans_map[k] = k
    return trans_map

def is_grid_sym(bits, trans_map, m):
    """Check if tiling is grid-symmetric."""
    for k in range(m):
        if ((bits >> k) & 1) != ((bits >> trans_map[k]) & 1):
            return False
    return True

# ─── Main census ───

def census(n, verbose=True):
    """Full census at n."""
    arcs = [(i, j) for i in range(n) for j in range(i + 2, n)]
    m = len(arcs)
    pos_indices = find_pos_arcs(n, arcs)
    sigma_pairs = find_sigma_pairs(n, arcs)
    trans_map = build_transpose_map(n, arcs)
    flip_mask = (1 << m) - 1

    if verbose:
        print(f"\n{'='*70}")
        print(f"n={n}: {m} arcs, {len(pos_indices)} POS, {len(sigma_pairs)} sigma-pairs")
        print(f"POS arcs: {[arcs[k] for k in pos_indices]}")
        print(f"2^m = {1 << m} tilings")
        print(f"{'='*70}")

    total = 1 << m

    # Phase 1: Compute H and canon for all tilings
    t0 = time.time()
    h_vals = [0] * total
    canon_vals = [None] * total
    score_vals = [None] * total

    for bits in range(total):
        if bits % 100000 == 0 and bits > 0 and verbose:
            elapsed = time.time() - t0
            rate = bits / elapsed
            eta = (total - bits) / rate
            print(f"  Phase 1: {bits}/{total} ({100*bits/total:.1f}%), "
                  f"{rate:.0f}/s, ETA {eta:.0f}s", flush=True)

        adj = build_adj_flat(n, bits, arcs)
        h = count_ham_dp(adj, n)
        h_vals[bits] = h
        scores = score_sequence(adj, n)
        score_vals[bits] = scores

        if n <= 7:
            canon = canon_hash(adj, n)
        else:
            canon = canon_hash(adj, n)

        canon_vals[bits] = canon

    elapsed = time.time() - t0
    if verbose:
        print(f"  Phase 1 done in {elapsed:.1f}s")

    # Phase 2: Group by class and compute properties
    class_members = defaultdict(list)
    for bits in range(total):
        class_members[canon_vals[bits]].append(bits)

    class_data = []
    gs_total = 0
    pos_global = Counter()

    for canon, members in sorted(class_members.items(), key=lambda kv: h_vals[kv[1][0]]):
        h = h_vals[members[0]]
        scores = score_vals[members[0]]
        size = len(members)

        # Grid-symmetric count
        gs = sum(1 for b in members if is_grid_sym(b, trans_map, m))
        gs_total += gs

        # Self-flip count
        sf = 0
        blue = 0
        black = 0
        for b in members:
            fb = b ^ flip_mask
            if canon_vals[fb] == canon:
                sf += 1
                if is_grid_sym(b, trans_map, m):
                    blue += 1
                else:
                    black += 1

        # POS patterns in this class
        pos_patterns = Counter()
        for b in members:
            pp = tuple((b >> k) & 1 for k in pos_indices)
            pos_patterns[pp] += 1
            pos_global[pp] += 1

        class_data.append({
            'h': h, 'size': size, 'scores': scores,
            'gs': gs, 'sf': sf, 'blue': blue, 'black': black,
            'pos_patterns': pos_patterns,
        })

    num_classes = len(class_data)

    if verbose:
        print(f"\n  Classes: {num_classes}")
        print(f"  Grid-symmetric tilings: {gs_total}/{total}")

        # POS distribution
        print(f"\n  POS orientation distribution ({len(pos_indices)} POS arcs):")
        for pp in sorted(pos_global.keys()):
            print(f"    POS={pp}: {pos_global[pp]} tilings ({100*pos_global[pp]/total:.2f}%)")

        # Class table
        print(f"\n  {'H':>5} {'Size':>6} {'GS':>4} {'SF':>4} {'Blue':>5} {'Black':>5} {'Scores'}")
        print(f"  {'-'*65}")
        for d in class_data:
            if n <= 7 or d['sf'] > 0 or d['gs'] > 0:
                print(f"  {d['h']:5d} {d['size']:6d} {d['gs']:4d} {d['sf']:4d} "
                      f"{d['blue']:5d} {d['black']:5d} {d['scores']}")

        # Summary
        sf_classes = sum(1 for d in class_data if d['sf'] > 0)
        blue_classes = sum(1 for d in class_data if d['blue'] > 0)
        black_classes = sum(1 for d in class_data if d['black'] > 0)
        total_sf = sum(d['sf'] for d in class_data)
        total_blue = sum(d['blue'] for d in class_data)
        total_black = sum(d['black'] for d in class_data)

        print(f"\n  Summary:")
        print(f"    Classes: {num_classes}")
        print(f"    Self-flip classes: {sf_classes} ({total_sf} tilings)")
        print(f"    Blueself classes: {blue_classes} ({total_blue} tilings)")
        print(f"    Blackself classes: {black_classes} ({total_black} tilings)")
        print(f"    Grid-symmetric: {gs_total} = 2^{gs_total.bit_length()-1 if gs_total > 0 else 0}")

        # POS per class: how many classes have each number of POS forward/backward
        print(f"\n  POS consumption per class:")
        for d in class_data:
            if len(d['pos_patterns']) > 1 and verbose and n <= 6:
                print(f"    H={d['h']}, size={d['size']}: POS patterns = {dict(d['pos_patterns'])}")

        # Class size distribution
        size_dist = Counter(d['size'] for d in class_data)
        print(f"\n  Class size distribution:")
        for s in sorted(size_dist.keys()):
            print(f"    Size {s}: {size_dist[s]} classes")

        # H value distribution
        h_dist = Counter(d['h'] for d in class_data)
        print(f"\n  H value distribution ({len(h_dist)} distinct values):")
        for h in sorted(h_dist.keys()):
            total_tilings = sum(d['size'] for d in class_data if d['h'] == h)
            print(f"    H={h}: {h_dist[h]} classes, {total_tilings} tilings")

    # Key formulas
    if verbose:
        print(f"\n  Key formulas:")
        print(f"    Grid-symmetric = 2^floor((n-1)^2/4) = 2^{(n-1)**2//4} = {2**((n-1)**2//4)}")
        print(f"    Actual GS = {gs_total}")
        assert gs_total == 2 ** ((n-1)**2 // 4), "GS count mismatch!"
        print(f"    GS formula VERIFIED")

        # Size = H / |Aut| for each class
        print(f"\n  Size = H / |Aut| verification:")
        all_divide = True
        for d in class_data:
            if d['h'] % d['size'] != 0:
                all_divide = False
                print(f"    FAIL: H={d['h']} not divisible by size={d['size']}")
        if all_divide:
            print(f"    All classes: H divisible by size (|Aut| = H/size)")

    return class_data


if __name__ == '__main__':
    target = int(sys.argv[1]) if len(sys.argv) > 1 else 8

    for n in range(3, min(target + 1, 9)):
        census(n)
