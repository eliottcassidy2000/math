"""Quick check: is I(Omega,2) = 1 + 2*#cycles + 4*#VD_pairs at n=6?"""

import sys
sys.path.insert(0, '03-artifacts/code')
from tournament_lib import (
    all_tournaments, find_odd_cycles, conflict_graph, independence_poly_at
)

n = 6
mismatches = 0
for t_idx, T in enumerate(all_tournaments(n)):
    if t_idx >= 500:
        break
    cycles = find_odd_cycles(T)
    if not cycles:
        io_gt = 1
    else:
        cg = conflict_graph(cycles)
        io_gt = independence_poly_at(cg, 2)

    # My formula
    vd_pairs = 0
    for a in range(len(cycles)):
        for b in range(a + 1, len(cycles)):
            if not (set(cycles[a]) & set(cycles[b])):
                vd_pairs += 1

    io_formula = 1 + 2 * len(cycles) + 4 * vd_pairs

    if io_gt != io_formula:
        mismatches += 1
        if mismatches <= 3:
            print(f"T#{t_idx}: gt={io_gt}, formula={io_formula}, "
                  f"#cycles={len(cycles)}, #VD_pairs={vd_pairs}")
            # Show the independent sets
            print(f"  Cycles: {cycles}")
            for a in range(len(cycles)):
                for b in range(a + 1, len(cycles)):
                    disjoint = not (set(cycles[a]) & set(cycles[b]))
                    if disjoint:
                        print(f"    VD: {cycles[a]} and {cycles[b]}")

print(f"\nMismatches: {mismatches}")
