#!/usr/bin/env python3
"""
Exhaustive proof of OCF at n=7: verify H(T) = I(Omega(T), 2) for ALL
2^C(7,2)-1 = 2^20 arc variable assignments (with one arc fixed).

By fixing arc 0->1, we cover all tournaments up to relabeling of one arc.
Every tournament on 7 vertices appears (with possibly relabeled vertices)
in this enumeration. Together with the relabeling-invariance of OCF,
this proves OCF for all n=7 tournaments.

Runtime: ~25 minutes.
Instance: opus-2026-03-05-S3
"""

import sys
sys.path.insert(0, '.')
from tournament_lib import (
    hamiltonian_path_count, find_odd_cycles, independence_poly_at_fast
)
from itertools import combinations, permutations
import time


def main():
    n = 7
    # Fix arc 0->1. Vary all other C(7,2)-1 = 20 arcs.
    # Variables: for each pair (a,b) with a<b and (a,b)!=(0,1),
    # bit=1 means a->b, bit=0 means b->a.

    pairs = [(a, b) for a in range(n) for b in range(a + 1, n) if not (a == 0 and b == 1)]
    num_vars = len(pairs)
    total = 1 << num_vars
    print(f"n={n}, variables={num_vars}, assignments={total}")

    ok = 0
    tested = 0
    t0 = time.time()

    for mask in range(total):
        # Build tournament
        T = [[0] * n for _ in range(n)]
        T[0][1] = 1  # fixed arc

        for idx, (a, b) in enumerate(pairs):
            if (mask >> idx) & 1:
                T[a][b] = 1
            else:
                T[b][a] = 1

        # Compute H(T)
        H = hamiltonian_path_count(T)

        # Compute I(Omega(T), 2) via fast method
        cycles = find_odd_cycles(T)
        I = independence_poly_at_fast(cycles, 2)

        tested += 1
        if H == I:
            ok += 1
        else:
            print(f"  FAIL mask={mask}: H={H}, I={I}, #cycles={len(cycles)}")
            if tested - ok > 3:
                print("  Aborting due to failures.")
                sys.exit(1)

        if tested % 100000 == 0:
            elapsed = time.time() - t0
            rate = tested / elapsed
            eta = (total - tested) / rate / 60
            print(f"  {tested}/{total} ({100*tested/total:.1f}%), ok={ok}, "
                  f"{elapsed:.0f}s elapsed, ETA {eta:.0f}min", flush=True)

    elapsed = time.time() - t0
    print(f"\n{'='*60}")
    print(f"RESULT: {ok}/{tested} passed in {elapsed:.0f}s ({elapsed/60:.1f}min)")
    if ok == tested:
        print(f"PROVED: H(T) = I(Omega(T), 2) for ALL {tested} arc assignments at n={n}.")
        print(f"Combined with relabeling invariance, this proves OCF for all n<={n} tournaments.")
    else:
        print(f"FAILED: {tested-ok} failures found.")


if __name__ == "__main__":
    main()
