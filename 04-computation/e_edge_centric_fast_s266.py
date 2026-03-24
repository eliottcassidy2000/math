#!/usr/bin/env python3
"""
e_edge_centric_fast_s266.py — Fast E(G_n) via single-arc edge-centric method
opus-2026-03-23-S266

By kind-pasteur S20dz: E(G_n) = #{distinct class pairs connected through arc (0,1)}.
This means: for each class rep T, flip arc (0,1), canonicalize, check if different.

This is C(n,2)× FASTER than building the full metagraph (which flips ALL arcs).

Expected speedup at n=9: 36× → ~5 min instead of 3 hours.
Expected time at n=10: ~2-3 hours.
"""

import subprocess, sys, time
from math import comb
from collections import defaultdict
from itertools import permutations

sys.stdout.reconfigure(line_buffering=True)

def get_reps(n):
    """Get class reps from gentourng."""
    r = subprocess.run(['gentourng', str(n)], capture_output=True, text=True)
    m = comb(n, 2)
    P = [(i,j) for i in range(n) for j in range(i+1,n)]
    lines = [l.strip() for l in (r.stdout or '').split('\n')
             if len(l.strip()) == m and all(c in '01' for c in l.strip())]
    return lines, P

def bits_to_adj(line, n, P):
    adj = [[0]*n for _ in range(n)]
    for k, (i,j) in enumerate(P):
        if line[k] == '1': adj[i][j] = 1
        else: adj[j][i] = 1
    return adj

def score_and_c3(adj, n):
    """Fast invariant computation."""
    scores = tuple(sorted(sum(adj[i][j] for j in range(n)) for i in range(n)))
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if adj[i][j] and adj[j][k] and adj[k][i]: c3 += 1
                if adj[i][k] and adj[k][j] and adj[j][i]: c3 += 1
    return scores, c3

def canonical_form(adj, n):
    """Canonical form using score-partitioned permutation search."""
    scores = [sum(adj[i][j] for j in range(n)) for i in range(n)]
    groups = defaultdict(list)
    for v in range(n):
        groups[scores[v]].append(v)
    sorted_groups = [groups[s] for s in sorted(groups.keys())]

    best = None
    count = 0
    def gen_perms(gs):
        if not gs:
            yield []
            return
        for p in permutations(gs[0]):
            for rest in gen_perms(gs[1:]):
                yield list(p) + rest

    for perm in gen_perms(sorted_groups):
        sig = tuple(adj[perm[i]][perm[j]] for i in range(n) for j in range(n) if i != j)
        if best is None or sig < best:
            best = sig
        count += 1
        if count > 500000:
            break
    return best

def compute_E_fast(n):
    """Compute E(G_n) using edge-centric single-arc method."""
    print(f"\n  n={n}: ", end='', flush=True)

    t0 = time.time()
    lines, P = get_reps(n)
    V = len(lines)
    m = len(P)
    print(f"{V} classes, ", end='', flush=True)

    # Build lookup: canonical_form → index
    t1 = time.time()
    # Use fast hash for pre-filtering
    hash_map = defaultdict(list)
    canon_map = {}  # canonical_form → index

    for idx, line in enumerate(lines):
        adj = bits_to_adj(line, n, P)
        scores, c3 = score_and_c3(adj, n)
        h = (scores, c3)
        hash_map[h].append(idx)

        cf = canonical_form(adj, n)
        canon_map[cf] = idx

        if idx > 0 and idx % 50000 == 0:
            elapsed = time.time() - t1
            rate = idx / elapsed
            eta = (V - idx) / rate
            print(f"\n    Indexed {idx}/{V} ({rate:.0f}/s, ETA {eta:.0f}s)", end='', flush=True)

    dt_index = time.time() - t1
    print(f"indexed in {dt_index:.1f}s, ", end='', flush=True)

    # For each class, flip arc (0,1), canonicalize, look up
    t2 = time.time()
    edges = set()
    self_loops = 0

    for idx, line in enumerate(lines):
        adj = bits_to_adj(line, n, P)

        # Flip arc (0,1) = position 0 in P
        new_adj = [row[:] for row in adj]
        new_adj[0][1] = 1 - adj[0][1]
        new_adj[1][0] = 1 - adj[1][0]

        # Canonicalize the flipped tournament
        cf = canonical_form(new_adj, n)

        # Look up
        nb = canon_map.get(cf)
        if nb is None:
            # This shouldn't happen
            continue

        if nb == idx:
            self_loops += 1
        else:
            edges.add((min(idx, nb), max(idx, nb)))

        if idx > 0 and idx % 50000 == 0:
            elapsed = time.time() - t2
            rate = idx / elapsed
            eta = (V - idx) / rate
            print(f"\n    Flipped {idx}/{V} ({rate:.0f}/s, ETA {eta:.0f}s, E_so_far={len(edges)})",
                  end='', flush=True)

    dt_flip = time.time() - t2
    E = len(edges)
    total = time.time() - t0

    print(f"E={E}, SL={self_loops}, total={total:.1f}s")
    return E, self_loops, total

if __name__ == '__main__':
    print("=" * 70)
    print("  FAST E(G_n) VIA EDGE-CENTRIC SINGLE-ARC METHOD")
    print("  opus-2026-03-23-S266")
    print("=" * 70)

    known_E = {3:1, 4:5, 5:30, 6:290, 7:4086, 8:91161, 9:3380751}

    # Validate n=3..8
    for n in range(3, 9):
        E, sl, dt = compute_E_fast(n)
        ek = known_E[n]
        match = "✓" if E == ek else f"✗ (expected {ek})"
        print(f"    → {match}")

    # Try n=9 (should be ~5 min instead of 3 hours)
    print(f"\n  ATTEMPTING n=9...")
    E9, sl9, dt9 = compute_E_fast(9)
    ek9 = known_E[9]
    match9 = "✓" if E9 == ek9 else f"✗ (expected {ek9})"
    print(f"    → n=9: E={E9}, {match9}, time={dt9:.1f}s")

    print("\nDONE.")
