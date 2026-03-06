#!/usr/bin/env python3
"""
Deep census analysis: POS consumption, blueself/blackself structure,
SC maximizer connection, and cross-n patterns.

Runs n=3 through n=7 with full canonicalization.
For n=7 (2^15 = 32768 tilings, 7! = 5040 perms), this is ~165M operations.

Optimizations:
- Flat adjacency arrays
- Hash-based pre-filtering for canonicalization
- Batch POS/GS/SF analysis

Instance: opus-2026-03-06-S3
"""
import sys
import time
from collections import defaultdict, Counter
from itertools import permutations

# ─── Core functions ───

def build_arcs(n):
    return [(i, j) for i in range(n) for j in range(i + 2, n)]

def build_adj_flat(n, bits, arcs):
    adj = bytearray(n * n)
    for i in range(n - 1):
        adj[i * n + (i + 1)] = 1
    for k in range(len(arcs)):
        i, j = arcs[k]
        if bits & (1 << k):
            adj[j * n + i] = 1
        else:
            adj[i * n + j] = 1
    return adj

def count_ham_dp(adj, n):
    full = (1 << n) - 1
    dp = {}
    for v in range(n):
        dp[(1 << v) * n + v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            c = dp.get(mask * n + v, 0)
            if c == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v * n + u]:
                    nk = (mask | (1 << u)) * n + u
                    dp[nk] = dp.get(nk, 0) + c
    return sum(dp.get(full * n + v, 0) for v in range(n))

def score_seq(adj, n):
    return tuple(sorted(sum(adj[i * n + j] for j in range(n)) for i in range(n)))

def canon_full(adj, n):
    best = None
    for perm in permutations(range(n)):
        s = tuple(adj[perm[i] * n + perm[j]] for i in range(n) for j in range(n))
        if best is None or s < best:
            best = s
    return best

def canon_hash_filtered(adj, n, score):
    """Canonical form with score-based filtering for speed."""
    # For n <= 7 just do full search — it's fast enough
    return canon_full(adj, n)

def find_pos_arcs(n, arcs):
    return [k for k, (i, j) in enumerate(arcs) if j == 2 * i + 2]

def build_transpose_map(n, arcs):
    arc_idx = {a: k for k, a in enumerate(arcs)}
    trans_map = [0] * len(arcs)
    for k, (i, j) in enumerate(arcs):
        ti, tj = n - 1 - j, n - 1 - i
        if ti > tj:
            ti, tj = tj, ti
        trans_map[k] = arc_idx.get((ti, tj), k)
    return trans_map

def is_grid_sym(bits, trans_map, m):
    for k in range(m):
        if ((bits >> k) & 1) != ((bits >> trans_map[k]) & 1):
            return False
    return True

# ─── Main census ───

def deep_census(n, verbose=True):
    arcs = build_arcs(n)
    m = len(arcs)
    pos_indices = find_pos_arcs(n, arcs)
    trans_map = build_transpose_map(n, arcs)
    flip_mask = (1 << m) - 1
    total = 1 << m

    if verbose:
        print(f"\n{'=' * 70}")
        print(f"n={n}: {m} arcs, {len(pos_indices)} POS, {total} tilings")
        print(f"POS arcs: {[arcs[k] for k in pos_indices]}")
        print(f"{'=' * 70}")

    t0 = time.time()

    # Phase 1: compute H, scores, canon for all tilings
    h_vals = [0] * total
    canon_vals = [None] * total

    for bits in range(total):
        if bits % 10000 == 0 and bits > 0 and verbose:
            elapsed = time.time() - t0
            rate = bits / elapsed
            eta = (total - bits) / rate
            print(f"  Phase 1: {bits}/{total} ({100*bits/total:.1f}%), "
                  f"{rate:.0f}/s, ETA {eta:.0f}s", flush=True)

        adj = build_adj_flat(n, bits, arcs)
        h = count_ham_dp(adj, n)
        h_vals[bits] = h
        sc = score_seq(adj, n)
        canon_vals[bits] = canon_full(adj, n)

    elapsed = time.time() - t0
    if verbose:
        print(f"  Phase 1 done in {elapsed:.1f}s")

    # Phase 2: group by class
    class_members = defaultdict(list)
    for bits in range(total):
        class_members[canon_vals[bits]].append(bits)

    # Phase 3: per-class analysis
    results = []
    for canon, members in sorted(class_members.items(), key=lambda kv: h_vals[kv[1][0]]):
        h = h_vals[members[0]]
        size = len(members)

        # Grid-symmetric members
        gs_members = [b for b in members if is_grid_sym(b, trans_map, m)]
        gs = len(gs_members)

        # Self-flip analysis
        sf_members = [b for b in members if canon_vals[b ^ flip_mask] == canon]
        sf = len(sf_members)

        blue_members = [b for b in sf_members if is_grid_sym(b, trans_map, m)]
        blue = len(blue_members)
        black = sf - blue

        # POS orientation patterns
        pos_patterns = Counter()
        for b in members:
            pp = tuple((b >> k) & 1 for k in pos_indices)
            pos_patterns[pp] += 1

        # POS pattern in GS members
        gs_pos = Counter()
        for b in gs_members:
            pp = tuple((b >> k) & 1 for k in pos_indices)
            gs_pos[pp] += 1

        # Score sequence (from adjacency)
        adj0 = build_adj_flat(n, members[0], arcs)
        scores = score_seq(adj0, n)

        # Self-converse check: does tournament sigma map class to itself?
        # Tournament sigma maps (i,j) -> (n-1-j, n-1-i)
        b0 = members[0]
        arc_idx = {a: k for k, a in enumerate(arcs)}
        sb = 0
        for k, (i, j) in enumerate(arcs):
            si, sj = n - 1 - j, n - 1 - i
            if si > sj:
                si, sj = sj, si
            sk = arc_idx.get((si, sj), k)
            if (b0 >> k) & 1:
                sb |= (1 << sk)
        is_sc = canon_vals[sb] == canon

        results.append({
            'h': h, 'size': size, 'scores': scores,
            'gs': gs, 'sf': sf, 'blue': blue, 'black': black,
            'is_sc': is_sc,
            'pos_patterns': pos_patterns,
            'gs_pos': gs_pos,
        })

    # ─── Report ───
    if verbose:
        print(f"\n  Classes: {len(results)}")
        gs_total = sum(d['gs'] for d in results)
        sf_total = sum(d['sf'] for d in results)
        blue_total = sum(d['blue'] for d in results)
        black_total = sum(d['black'] for d in results)
        sc_count = sum(1 for d in results if d['is_sc'])

        print(f"  Grid-symmetric tilings: {gs_total} (expected {2**((n-1)**2//4)})")
        print(f"  Self-converse classes: {sc_count}/{len(results)}")
        print(f"  Self-flip tilings: {sf_total} (blue={blue_total}, black={black_total})")

        # ─── POS uniformity verification ───
        print(f"\n  POS UNIFORMITY:")
        pos_global = Counter()
        for d in results:
            for pp, c in d['pos_patterns'].items():
                pos_global[pp] += c
        for pp in sorted(pos_global.keys()):
            expected = total // (2 ** len(pos_indices))
            actual = pos_global[pp]
            match = "UNIFORM" if actual == expected else f"SKEWED ({actual} vs {expected})"
            print(f"    POS={pp}: {actual} ({match})")

        # ─── Class table ───
        print(f"\n  {'H':>5} {'Size':>6} {'GS':>4} {'SF':>4} {'Blue':>5} {'Black':>5} {'SC?':>3} {'Scores'}")
        print(f"  {'-' * 70}")
        for d in results:
            sc_marker = 'Y' if d['is_sc'] else 'N'
            if d['sf'] > 0 or d['gs'] > 0 or d['is_sc']:
                print(f"  {d['h']:5d} {d['size']:6d} {d['gs']:4d} {d['sf']:4d} "
                      f"{d['blue']:5d} {d['black']:5d}   {sc_marker} {d['scores']}")

        # ─── POS consumption per class ───
        print(f"\n  POS CONSUMPTION PER CLASS:")
        for d in results:
            if len(d['pos_patterns']) > 1:
                print(f"    H={d['h']:5d}, size={d['size']:4d}, "
                      f"POS={dict(sorted(d['pos_patterns'].items()))}")

        # ─── SC maximizer within score classes ───
        print(f"\n  SC MAXIMIZER WITHIN SCORE CLASSES:")
        score_groups = defaultdict(list)
        for d in results:
            score_groups[d['scores']].append(d)
        for scores in sorted(score_groups.keys()):
            group = score_groups[scores]
            if len(group) > 1:
                sc_h = [d['h'] for d in group if d['is_sc']]
                nsc_h = [d['h'] for d in group if not d['is_sc']]
                if sc_h and nsc_h:
                    print(f"    scores={scores}: SC H={sorted(sc_h)}, NSC H={sorted(nsc_h)}")
                    print(f"      MAX is {'SC' if max(sc_h) >= max(nsc_h) else 'NSC'}")

        # ─── Blueself position within H ranking ───
        if blue_total > 0:
            print(f"\n  BLUESELF POSITION IN H RANKING:")
            blue_classes = [d for d in results if d['blue'] > 0]
            max_h = max(d['h'] for d in results)
            for d in blue_classes:
                rank = sum(1 for d2 in results if d2['h'] > d['h']) + 1
                print(f"    H={d['h']} (rank {rank}/{len(results)}), "
                      f"scores={d['scores']}, "
                      f"size={d['size']}, blue={d['blue']}, gs={d['gs']}")

        # ─── GS tilings: POS pattern distribution ───
        print(f"\n  GS TILINGS POS PATTERN:")
        gs_pos_global = Counter()
        for d in results:
            for pp, c in d['gs_pos'].items():
                gs_pos_global[pp] += c
        if gs_pos_global:
            for pp in sorted(gs_pos_global.keys()):
                print(f"    POS={pp}: {gs_pos_global[pp]} GS tilings")

        # ─── Blackself class analysis ───
        if black_total > 0:
            print(f"\n  BLACKSELF CLASSES:")
            black_classes = [d for d in results if d['black'] > 0]
            for d in black_classes:
                rank = sum(1 for d2 in results if d2['h'] > d['h']) + 1
                print(f"    H={d['h']} (rank {rank}/{len(results)}), "
                      f"scores={d['scores']}, "
                      f"size={d['size']}, black={d['black']}, "
                      f"SC={'Y' if d['is_sc'] else 'N'}")

        # ─── Key pattern: SF-class implies near-regular scores? ───
        print(f"\n  SELF-FLIP SCORE REGULARITY:")
        for d in results:
            if d['sf'] > 0:
                scores = d['scores']
                mean = sum(scores) / n
                var = sum((s - mean)**2 for s in scores) / n
                print(f"    H={d['h']}, scores={scores}, "
                      f"variance={var:.2f}, SC={'Y' if d['is_sc'] else 'N'}")

    return results


if __name__ == '__main__':
    target = int(sys.argv[1]) if len(sys.argv) > 1 else 7

    all_results = {}
    for n_val in range(3, min(target + 1, 8)):
        all_results[n_val] = deep_census(n_val)

    # ─── Cross-n summary ───
    print(f"\n\n{'=' * 70}")
    print("CROSS-n SUMMARY")
    print(f"{'=' * 70}")

    print(f"\n  {'n':>3} {'Classes':>8} {'SC':>5} {'SF':>5} {'Blue':>5} {'Black':>6} "
          f"{'GS':>6} {'Max H':>7} {'Max SC?':>7}")
    print(f"  {'-' * 65}")
    for n_val in sorted(all_results.keys()):
        res = all_results[n_val]
        n_classes = len(res)
        n_sc = sum(1 for d in res if d['is_sc'])
        n_sf = sum(1 for d in res if d['sf'] > 0)
        n_blue_c = sum(1 for d in res if d['blue'] > 0)
        n_black_c = sum(1 for d in res if d['black'] > 0)
        gs_total = sum(d['gs'] for d in res)
        max_h = max(d['h'] for d in res)
        max_sc = any(d['is_sc'] for d in res if d['h'] == max_h)
        print(f"  {n_val:3d} {n_classes:8d} {n_sc:5d} {n_sf:5d} {n_blue_c:5d} {n_black_c:6d} "
              f"{gs_total:6d} {max_h:7d} {'Y' if max_sc else 'N':>7}")

    # ─── Self-flip tiling fraction ───
    print(f"\n  Self-flip fraction:")
    for n_val in sorted(all_results.keys()):
        res = all_results[n_val]
        m = len(build_arcs(n_val))
        total = 1 << m
        sf_tilings = sum(d['sf'] for d in res)
        blue_tilings = sum(d['blue'] for d in res)
        black_tilings = sum(d['black'] for d in res)
        print(f"    n={n_val}: SF={sf_tilings}/{total} ({100*sf_tilings/total:.2f}%), "
              f"blue={blue_tilings}, black={black_tilings}")

    # ─── Blueself parity pattern ───
    print(f"\n  Blueself by parity:")
    for n_val in sorted(all_results.keys()):
        res = all_results[n_val]
        blue_tilings = sum(d['blue'] for d in res)
        parity = 'even' if n_val % 2 == 0 else 'odd'
        print(f"    n={n_val} ({parity}): {blue_tilings} blueself tilings")
