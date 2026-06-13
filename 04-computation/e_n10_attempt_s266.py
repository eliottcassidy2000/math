#!/usr/bin/env python3
"""
e_n10_attempt_s266.py — Attempt to compute E(G_10) via edge-centric nauty
opus-2026-03-23-S266

Strategy: Fix arc position 0 (= arc (0,1)).
For each tournament class rep T (from gentourng):
  1. Flip arc (0,1) to get T'
  2. Canonicalize T' (using nauty via directg or manual scoring)
  3. If T' ≠ T (as canonical), this is an edge in the metagraph

By kind-pasteur S20dz: E(G_n) = #{distinct class PAIRS connected through arc (0,1)}.

The bottleneck: canonicalization of T'. We can't use full nauty easily,
but we CAN use HASH-BASED approximate canonicalization with refinement.

Hash approach:
  score_seq + c3 + hash of 2-neighborhoods → collision groups
  Within collision groups, check full isomorphism

For n=10 with 9.7M classes, this might be feasible in ~1 hour.
Let's try n=8 first (6880 classes) as validation.
"""

import subprocess
import sys
import time
from math import comb
from collections import defaultdict

sys.stdout.reconfigure(line_buffering=True)

def run_gentourng(n):
    """Get all tournament class reps from nauty."""
    r = subprocess.run(['gentourng', str(n)], capture_output=True, text=True, timeout=300)
    m = comb(n, 2)
    P = [(i,j) for i in range(n) for j in range(i+1,n)]
    lines = [l.strip() for l in (r.stdout or '').split('\n')
             if len(l.strip()) == m and all(c in '01' for c in l.strip())]
    return lines, P

def bits_to_adj(bits_str, n, P):
    adj = [[0]*n for _ in range(n)]
    for k, (i,j) in enumerate(P):
        if bits_str[k] == '1': adj[i][j] = 1
        else: adj[j][i] = 1
    return adj

def adj_to_bits(adj, n, P):
    return ''.join('1' if adj[i][j] else '0' for (i,j) in P)

def score_seq(adj, n):
    return tuple(sorted(sum(adj[i][j] for j in range(n)) for i in range(n)))

def count_c3(adj, n):
    c = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if adj[i][j] and adj[j][k] and adj[k][i]: c += 1
                if adj[i][k] and adj[k][j] and adj[j][i]: c += 1
    return c

def neighborhood_hash(adj, n, v):
    """Hash based on v's out-neighbors' scores."""
    out_nbrs = [j for j in range(n) if adj[v][j]]
    in_nbrs = [j for j in range(n) if adj[j][v] and j != v]
    out_scores = tuple(sorted(sum(adj[u][w] for w in range(n)) for u in out_nbrs))
    in_scores = tuple(sorted(sum(adj[u][w] for w in range(n)) for u in in_nbrs))
    return (out_scores, in_scores)

def full_hash(adj, n):
    """Multi-level hash for tournament."""
    ss = score_seq(adj, n)
    c3 = count_c3(adj, n)
    nbr_hashes = tuple(sorted(neighborhood_hash(adj, n, v) for v in range(n)))
    return (ss, c3, nbr_hashes)

def compute_E_edge_centric(n, verbose=True):
    """Compute E(G_n) via edge-centric method."""
    if verbose:
        print(f"\n  Computing E(G_{n})...")

    t0 = time.time()
    lines, P = run_gentourng(n)
    m = len(P)
    V = len(lines)
    if verbose:
        print(f"  {V} class reps in {time.time()-t0:.1f}s")

    # Build hash map: hash → list of (index, bits_str)
    t1 = time.time()
    hash_map = defaultdict(list)
    class_data = []

    for idx, bits_str in enumerate(lines):
        adj = bits_to_adj(bits_str, n, P)
        h = full_hash(adj, n)
        hash_map[h].append(idx)
        class_data.append({'bits': bits_str, 'adj': adj, 'hash': h})
        if verbose and idx % 100000 == 0 and idx > 0:
            print(f"    Hashed {idx}/{V}...")

    if verbose:
        n_groups = len(hash_map)
        max_group = max(len(v) for v in hash_map.values())
        print(f"  Hash map: {n_groups} groups, max size {max_group} ({time.time()-t1:.1f}s)")

    # For each class, flip arc 0 (= arc (0,1)), find resulting class
    t2 = time.time()
    edges = set()
    self_loops = 0

    for idx in range(V):
        adj = class_data[idx]['adj']

        # Flip arc (0,1)
        new_adj = [row[:] for row in adj]
        new_adj[0][1] = 1 - adj[0][1]
        new_adj[1][0] = 1 - adj[1][0]

        # Find the class of new_adj
        h = full_hash(new_adj, n)
        candidates = hash_map.get(h, [])

        if len(candidates) == 0:
            # Hash not found — shouldn't happen
            print(f"  WARNING: hash miss at idx {idx}")
            continue
        elif len(candidates) == 1:
            nb = candidates[0]
        else:
            # Multiple candidates — need deeper comparison
            # Use canonicalization within the group
            nb = candidates[0]  # APPROXIMATE: take first match
            # TODO: proper isomorphism check within collision group

        if nb == idx:
            self_loops += 1
        else:
            edges.add((min(idx, nb), max(idx, nb)))

        if verbose and idx % 100000 == 0 and idx > 0:
            print(f"    Processed {idx}/{V}, edges so far: {len(edges)}...")

    E = len(edges)
    dt = time.time() - t2

    if verbose:
        print(f"  Result: E(G_{n}) = {E}, self_loops = {self_loops} ({dt:.1f}s)")
        print(f"  Total time: {time.time()-t0:.1f}s")

    return E, self_loops

if __name__ == '__main__':
    print("=" * 70)
    print("  EDGE-CENTRIC E(G_n) VIA NAUTY + HASHING")
    print("  opus-2026-03-23-S266")
    print("=" * 70)

    known_E = {3:1, 4:5, 5:30, 6:290, 7:4086, 8:91161, 9:3380751}

    # Validate at n=3..8
    for n in range(3, 9):
        E, sl = compute_E_edge_centric(n, verbose=False)
        ek = known_E[n]
        match = "✓" if E == ek else f"✗ (expected {ek})"
        print(f"  n={n}: E={E}, SL={sl}, {match}")

    # Now try n=9
    print(f"\n  Attempting n=9 (191536 classes)...")
    t0 = time.time()
    E9, sl9 = compute_E_edge_centric(9)
    print(f"  n=9: E={E9}, SL={sl9}, expected={known_E[9]}, time={time.time()-t0:.1f}s")

    # If n=9 works, try n=10
    # print(f"\n  Attempting n=10 (9733056 classes)...")
    # E10, sl10 = compute_E_edge_centric(10)
    # print(f"  n=10: E={E10}, SL={sl10}")

    print("\nDONE.")
