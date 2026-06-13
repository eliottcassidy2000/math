#!/usr/bin/env python3
"""
Compute H(T_19): number of Hamiltonian paths in the Paley tournament on 19 vertices.

Paley tournament T_19:
  - Vertices: {0, 1, ..., 18}
  - Arc i->j iff (j - i) mod 19 is a quadratic residue mod 19
  - QR(19) = {1, 4, 5, 6, 7, 9, 11, 16, 17}

Method: Bitmask DP (pure Python with numpy int64 array for storage).
  dp[mask][v] = number of Hamiltonian paths visiting exactly the vertices in `mask`,
                ending at vertex v.

States: 2^19 * 19 = 9,961,472 entries.
"""

import numpy as np
import time
import math

def compute():
    n = 19
    p = 19
    QR19 = {1, 4, 5, 6, 7, 9, 11, 16, 17}

    print(f"=== Computing H(T_{p}) ===")
    print(f"Paley tournament on {p} vertices")
    print(f"QR({p}) = {sorted(QR19)}")

    # Build adjacency bitmasks: adj_out[i] = bitmask of vertices j such that i->j
    adj_out = [0] * n
    for i in range(n):
        for j in range(n):
            if i != j and (j - i) % p in QR19:
                adj_out[i] |= (1 << j)

    # Also build adj_in[j] = bitmask of vertices i such that i->j
    adj_in = [0] * n
    for j in range(n):
        for i in range(n):
            if i != j and (j - i) % p in QR19:
                adj_in[j] |= (1 << i)

    # Verify out-degrees
    for i in range(n):
        od = bin(adj_out[i]).count('1')
        assert od == 9, f"Vertex {i} has out-degree {od}, expected 9"
    print(f"All vertices have out-degree 9 (verified)")

    total_masks = 1 << n  # 524288
    full_mask = total_masks - 1

    print(f"DP table: {total_masks:,} masks x {n} vertices = {total_masks * n:,} entries")
    print(f"Memory: ~{total_masks * n * 8 / 1024 / 1024:.0f} MB (numpy int64)")
    print()

    # Allocate DP table as numpy array: dp[mask, v]
    # Using int64. Max possible value at any cell bounded by 18! ~ 6.4e15 < 2^63.
    dp = np.zeros((total_masks, n), dtype=np.int64)

    # Base case: single vertex
    for v in range(n):
        dp[1 << v, v] = 1

    print("Starting DP computation...")
    start_time = time.time()
    last_report = start_time

    # Main DP loop: iterate over all masks in order.
    # For each (mask, v) with dp[mask,v] > 0, extend to successors of v not in mask.
    #
    # Key optimization: use "push" DP -- from each state, push to successors.
    # This avoids checking all vertices for each mask; we only process nonzero entries.

    # Convert adj_out to a list for fast access
    adj_out_list = adj_out

    # For speed, we'll access the numpy array's underlying buffer via item/itemset
    # or just use direct indexing. Let's try the straightforward approach first.

    for mask in range(1, total_masks):
        # For progress reporting
        if mask & 0xFFFF == 0:
            now = time.time()
            if now - last_report >= 30:
                elapsed = now - start_time
                pct = mask / total_masks * 100
                print(f"  Progress: {pct:.1f}% ({mask:,}/{total_masks:,}), elapsed: {elapsed:.0f}s")
                last_report = now

        # Get the row for this mask
        row = dp[mask]

        # Find vertices in mask that have nonzero dp values
        # Available extension targets: bits NOT in mask
        available_bits = (~mask) & full_mask
        if available_bits == 0:
            continue  # full mask, no extensions possible

        for v in range(n):
            count = row[v]
            if count == 0:
                continue

            # Successors of v that are not yet in mask
            targets = adj_out_list[v] & available_bits
            # Iterate over set bits of targets
            t = targets
            while t:
                w = (t & -t).bit_length() - 1
                new_mask = mask | (1 << w)
                dp[new_mask, w] += count
                t &= t - 1

    elapsed = time.time() - start_time
    print(f"\nDP computation finished in {elapsed:.1f} seconds.")

    # Sum over all ending vertices for the full mask
    H_T19 = int(np.sum(dp[full_mask]))

    # Print per-endpoint counts
    print(f"\nPer-endpoint counts:")
    for v in range(n):
        val = int(dp[full_mask, v])
        if val > 0:
            print(f"  endpoint {v}: {val}")

    print(f"\n{'='*60}")
    print(f"H(T_{p}) = {H_T19}")
    print(f"{'='*60}")

    # Vertex-transitivity check: all endpoints should have equal count
    endpoint_counts = [int(dp[full_mask, v]) for v in range(n)]
    if len(set(endpoint_counts)) == 1:
        print(f"Vertex-transitive check: PASS (all endpoints = {endpoint_counts[0]})")
    else:
        print(f"Vertex-transitive check: FAIL (endpoints vary)")

    # Automorphism group
    # |Aut(T_p)| = p * (p-1)/2 for Paley tournaments (p >= 7)
    # For p=19: 19 * 9 = 171
    aut_size = p * (p - 1) // 2
    print(f"\n|Aut(T_{p})| = {aut_size}")

    quotient, remainder = divmod(H_T19, aut_size)
    print(f"H(T_{p}) / |Aut(T_{p})| = {H_T19} / {aut_size} = {quotient}")
    if remainder == 0:
        print(f"  (Exact division confirmed)")
    else:
        print(f"  WARNING: remainder = {remainder}")

    # Reference values
    fact_19 = math.factorial(19)
    print(f"\n19! = {fact_19}")
    print(f"H(T_{p}) / 19! = {H_T19 / fact_19:.10f}")

    # Parity
    print(f"\nH(T_{p}) mod 2 = {H_T19 % 2}")
    print(f"H(T_{p}) is {'odd' if H_T19 % 2 == 1 else 'even'}")

    # OEIS A038375 note
    print(f"\n=== OEIS A038375 reference ===")
    print(f"A038375 gives the maximum number of Hamiltonian paths over all")
    print(f"tournaments on n vertices. If H(T_19) matches the OEIS value for")
    print(f"n=19, this confirms the Paley Maximizer Conjecture at p=19.")
    print(f"(Look up OEIS A038375 manually to compare.)")

    return H_T19


if __name__ == '__main__':
    compute()
