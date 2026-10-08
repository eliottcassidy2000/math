#!/usr/bin/env python3
"""Property (i) of THM-4611 for all 27 certified maps: the identity is the only coupling j -> a j + b (a in the coupling group) whose roots all
vanish, i.e. which preserves the multiplier vector.  Exact check for the maps certified in this session."""
from block_balance import ratio_group
MAPS = [
    # THM-4611 (3): fixed-length certificates
    (5, [1, 2, 3, 7, 11]), (5, [1, 2, 3, 7, 1]), (5, [1, 1, 2, 3, 7]), (7, [1, 1, 1, 1, 2, 3, 5]), (7, [1, 1, 1, 2, 3, 5, 11]),
    # THM-4611 (3'): adaptive, coupling group {+-1} (these six are also the census class (3, 2))
    (5, [1, 6, 11, 11, 4]), (5, [1, 1, 6, 39, 11]), (5, [1, 14, 4, 1, 29]), (5, [1, 4, 31, 1, 6]), (5, [1, 4, 1, 11, 34]), (5, [1, 6, 4, 31, 1]),
    # THM-4611 (4): Z_4
    (4, [1, 3, 5, 7]),
    # THM-4611 (5): census class (3, 4)
    (5, [1, 1, 28, 11, 7]), (5, [1, 3, 7, 24, 1]), (5, [1, 6, 8, 21, 3]), (5, [1, 3, 9, 7, 11]), (5, [1, 11, 1, 21, 2]), (5, [1, 2, 3, 14, 8]),
    # THM-4611 (5): census class (4, 4)
    (5, [1, 12, 2, 17, 7]), (5, [1, 23, 2, 3, 21]), (5, [1, 7, 19, 11, 2]), (5, [1, 2, 3, 31, 14]), (5, [1, 7, 2, 22, 6]), (5, [1, 6, 11, 7, 2]),
    # THM-4611 (6): audit F's evidence maps
    (5, [1, 8, 3, 7, 12]), (5, [1, 2, 3, 7, 6]), (7, [1, 2, 3, 5, 1, 1, 1]),
]
assert len(MAPS) == 27 and len({(d, tuple(m)) for d, m in MAPS}) == 27
for d, m in MAPS:
    G = ratio_group(m, d)
    bad = [(a, b) for a in G for b in range(d) if (a, b) != (1, 0) and all(m[(a * j + b) % d] == m[j] for j in range(d))]
    print(f"Z_{d} m = {m}: coupling group {G}; non-identity couplings preserving m: {bad} -> property (i) {'holds' if not bad else 'FAILS'}")
