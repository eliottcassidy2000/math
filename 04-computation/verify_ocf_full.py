"""
Full exhaustive verification: H(T) = I(Omega(T), 2) for all n <= 6.
This confirms that OCF is a genuine closed-form formula, not just recursive.

Author: opus-2026-03-05-S2
"""
import sys
import os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import (
    all_tournaments, hamiltonian_path_count,
    find_odd_cycles, conflict_graph, independence_poly_at
)

for n in [3, 4, 5, 6]:
    failures = 0
    total = 0
    for T in all_tournaments(n):
        total += 1
        ht = hamiltonian_path_count(T)
        cycles = find_odd_cycles(T)
        if not cycles:
            io = 1
        else:
            cg = conflict_graph(cycles)
            io = independence_poly_at(cg, 2)
        if ht != io:
            failures += 1
            if failures <= 3:
                print(f"  FAILURE: n={n}, H(T)={ht}, I(Omega,2)={io}, #cycles={len(cycles)}")
    print(f"n={n}: {total} tournaments, {failures} failures {'ALL PASSED' if failures == 0 else 'FAILURES FOUND'}")

print("\nIf all passed: H(T) = I(Omega(T), 2) is PROVEN for n <= 6")
print("(conditional on Claim A, which is verified for n <= 6)")
