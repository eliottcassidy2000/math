#!/usr/bin/env python3
"""Property (i) of THM-4611: the identity is the only coupling j -> a j + b (a in the coupling group) whose roots all
vanish, i.e. which preserves the multiplier vector.  Exact check for the maps certified in this session."""
from block_balance import ratio_group
MAPS = [(5, [1, 2, 3, 7, 11]), (5, [1, 2, 3, 7, 1]), (5, [1, 1, 2, 3, 7]), (7, [1, 1, 1, 1, 2, 3, 5]), (7, [1, 1, 1, 2, 3, 5, 11]),
        (5, [1, 6, 11, 11, 4]), (5, [1, 1, 6, 39, 11]), (5, [1, 14, 4, 1, 29]), (5, [1, 4, 31, 1, 6]), (5, [1, 4, 1, 11, 34]),
        (5, [1, 6, 4, 31, 1]), (5, [1, 1, 28, 11, 7]), (5, [1, 3, 7, 24, 1])]
for d, m in MAPS:
    G = ratio_group(m, d)
    bad = [(a, b) for a in G for b in range(d) if (a, b) != (1, 0) and all(m[(a * j + b) % d] == m[j] for j in range(d))]
    print(f"Z_{d} m = {m}: coupling group {G}; non-identity couplings preserving m: {bad} -> property (i) {'holds' if not bad else 'FAILS'}")
