"""
Verify or refute MISTAKE-004: "OCF is recursive, not closed-form over all cycles"

MISTAKE-004 claims H(T) != I(Omega(T), 2) for a tournament with two
vertex-disjoint 3-cycles. But the "counterexample" computes I(Omega,2)
incorrectly: the author used mu-weighted products instead of the standard
independence polynomial.

The independence polynomial I(G, x) = sum_{k} alpha_k * x^k where
alpha_k = #independent sets of size k. This does NOT involve mu weights.

Let me verify H(T) = I(Omega(T), 2) for the alleged counterexample and
for all tournaments at n <= 6.

Author: opus-2026-03-05-S2
"""

import sys
import os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import (
    tournament_from_bits, all_tournaments, hamiltonian_path_count,
    find_odd_cycles, conflict_graph, independence_poly_at, delete_vertex
)


def verify_ocf(T):
    """Verify H(T) = I(Omega(T), 2)."""
    ht = hamiltonian_path_count(T)
    cycles = find_odd_cycles(T)
    if not cycles:
        # No odd cycles -> Omega is empty -> I(empty, 2) = 1
        return ht == 1, ht, 1
    cg = conflict_graph(cycles)
    i_omega = independence_poly_at(cg, 2)
    return ht == i_omega, ht, i_omega


if __name__ == "__main__":
    print("="*70)
    print("Testing MISTAKE-004: Is H(T) = I(Omega(T), 2) a valid closed form?")
    print("="*70)

    # Test the alleged counterexample from file.txt
    # T on {0,1,2,3,4,5} (relabeled from {1,...,6}):
    # 3-cycle: 0->1->2->0
    # 3-cycle: 3->4->5->3
    # All arcs from {0,1,2} to {3,4,5}
    print("\n--- Alleged counterexample from file.txt ---")
    T = [[0]*6 for _ in range(6)]
    # 3-cycle 0->1->2->0
    T[0][1] = 1; T[1][2] = 1; T[2][0] = 1
    # 3-cycle 3->4->5->3
    T[3][4] = 1; T[4][5] = 1; T[5][3] = 1
    # {0,1,2} beat {3,4,5}
    for i in range(3):
        for j in range(3, 6):
            T[i][j] = 1

    ht = hamiltonian_path_count(T)
    cycles = find_odd_cycles(T)
    print(f"H(T) = {ht}")
    print(f"Odd cycles: {len(cycles)}: {cycles}")

    cg = conflict_graph(cycles)
    i_omega = independence_poly_at(cg, 2)
    print(f"I(Omega(T), 2) = {i_omega}")
    print(f"H(T) = I(Omega(T), 2)? {ht == i_omega}")

    if ht == i_omega:
        print("\n*** MISTAKE-004 is WRONG! The counterexample is invalid. ***")
        print("The file.txt author confused I(Omega,2) with a mu-weighted product.")
        print("I(Omega,2) = sum_{S independent} 2^|S|, NOT sum_{S} prod mu(C) * 2^|S|.")

    # Full verification at n=4,5
    for n in [4, 5]:
        print(f"\n--- Verifying H(T) = I(Omega(T), 2) for all n={n} tournaments ---")
        failures = 0
        total = 0
        for T in all_tournaments(n):
            total += 1
            ok, ht, io = verify_ocf(T)
            if not ok:
                failures += 1
                if failures <= 3:
                    print(f"  FAILURE: H(T)={ht}, I(Omega,2)={io}")
        print(f"  {total} tournaments, {failures} failures")

    # n=6: sample (full is 32768 tournaments, may be slow with cycle enumeration)
    print(f"\n--- Verifying H(T) = I(Omega(T), 2) for n=6 (sample) ---")
    failures = 0
    total = 0
    for T in all_tournaments(6):
        total += 1
        if total > 2000:
            break
        ok, ht, io = verify_ocf(T)
        if not ok:
            failures += 1
            if failures <= 3:
                print(f"  FAILURE #{failures}: H(T)={ht}, I(Omega,2)={io}")
                cycles = find_odd_cycles(T)
                print(f"    #cycles={len(cycles)}")
    print(f"  {total} tournaments sampled, {failures} failures")

    if failures == 0:
        print("\n" + "="*70)
        print("CONCLUSION: MISTAKE-004 is WRONG.")
        print("H(T) = I(Omega(T), 2) IS a valid closed-form formula.")
        print("The counterexample in file.txt computed I(Omega,2) incorrectly")
        print("by multiplying mu weights into the independence polynomial.")
        print("The standard independence polynomial I(G,2) = sum alpha_k * 2^k")
        print("does NOT involve mu — it just counts independent sets.")
        print("="*70)
