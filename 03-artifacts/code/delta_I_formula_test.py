#!/usr/bin/env python3
"""
Test the conjectured formula: delta_I = 2 * (n_gained - n_lost)
under arc flip, where n_gained = # new odd cycles gained, n_lost = # old odd cycles lost.

Also tests: delta_H = 2 * (n_gained - n_lost) [since delta_H = delta_I].

If this holds, it means:
  H(T') - H(T) = 2 * (#{cycles using new arc} - #{cycles using old arc})

This would be a clean combinatorial formula connecting Hamiltonian paths to cycle counts.

Instance: kind-pasteur-2026-03-05-S5
"""

import sys
sys.path.insert(0, r"C:\Users\Eliott\Documents\GitHub\math\03-artifacts\code")
from tournament_lib import *
from collections import Counter
import random


def flip_arc(T, i, j):
    n = len(T)
    Tp = [row[:] for row in T]
    Tp[i][j] = 1 - T[i][j]
    Tp[j][i] = 1 - T[j][i]
    return Tp


def compute_I(T):
    cycles = find_odd_cycles(T)
    if not cycles:
        return 1
    cg = conflict_graph(cycles)
    return independence_poly_at(cg, 2)


def uses_arc(cyc, a, b):
    L = len(cyc)
    return any(cyc[k] == a and cyc[(k+1) % L] == b for k in range(L))


def test_formula(n, exhaustive=True, n_random=1000):
    """Test delta_I = 2*(gained - lost) and delta_H = 2*(gained - lost)."""
    print(f"\n=== Testing n={n} {'(exhaustive)' if exhaustive else f'(random {n_random})'} ===")

    formula_I_ok = 0
    formula_I_fail = 0
    formula_H_ok = 0
    formula_H_fail = 0
    total = 0

    rng = random.Random(42)

    if exhaustive:
        tournament_iter = all_tournaments(n)
    else:
        tournament_iter = (random_tournament(n, rng) for _ in range(n_random))

    fail_examples_I = []
    fail_examples_H = []

    for T in tournament_iter:
        ht = hamiltonian_path_count(T)
        it = compute_I(T)
        cyc_T = find_odd_cycles(T)

        for i in range(n):
            for j in range(i+1, n):
                if not T[i][j]:
                    continue  # only test i->j flips (covers all by symmetry)

                Tp = flip_arc(T, i, j)
                htp = hamiltonian_path_count(Tp)
                itp = compute_I(Tp)
                cyc_Tp = find_odd_cycles(Tp)

                lost = sum(1 for c in cyc_T if uses_arc(c, i, j))
                gained = sum(1 for c in cyc_Tp if uses_arc(c, j, i))

                delta_I = itp - it
                delta_H = htp - ht
                predicted = 2 * (gained - lost)

                total += 1

                if delta_I == predicted:
                    formula_I_ok += 1
                else:
                    formula_I_fail += 1
                    if len(fail_examples_I) < 5:
                        fail_examples_I.append({
                            'i': i, 'j': j,
                            'delta_I': delta_I, 'predicted': predicted,
                            'lost': lost, 'gained': gained,
                            'n_cyc_T': len(cyc_T), 'n_cyc_Tp': len(cyc_Tp),
                        })

                if delta_H == predicted:
                    formula_H_ok += 1
                else:
                    formula_H_fail += 1
                    if len(fail_examples_H) < 5:
                        fail_examples_H.append({
                            'i': i, 'j': j,
                            'delta_H': delta_H, 'predicted': predicted,
                            'lost': lost, 'gained': gained,
                        })

    print(f"Total tests: {total}")
    print(f"delta_I = 2*(gained-lost): {formula_I_ok}/{total} pass, {formula_I_fail} fail")
    print(f"delta_H = 2*(gained-lost): {formula_H_ok}/{total} pass, {formula_H_fail} fail")

    if fail_examples_I:
        print("\nFailed examples (delta_I):")
        for ex in fail_examples_I:
            print(f"  Flip {ex['i']}->{ex['j']}: delta_I={ex['delta_I']}, predicted={ex['predicted']}, "
                  f"lost={ex['lost']}, gained={ex['gained']}, "
                  f"#cyc_T={ex['n_cyc_T']}, #cyc_Tp={ex['n_cyc_Tp']}")

    if fail_examples_H:
        print("\nFailed examples (delta_H):")
        for ex in fail_examples_H:
            print(f"  Flip {ex['i']}->{ex['j']}: delta_H={ex['delta_H']}, predicted={ex['predicted']}, "
                  f"lost={ex['lost']}, gained={ex['gained']}")

    return formula_I_fail == 0 and formula_H_fail == 0


def main():
    print("=== Testing: delta_I = delta_H = 2*(#gained_cycles - #lost_cycles) ===")

    # Exhaustive for n <= 5
    for n in range(3, 6):
        test_formula(n, exhaustive=True)

    # Random for n=6
    test_formula(6, exhaustive=False, n_random=500)

    # Random for n=7
    test_formula(7, exhaustive=False, n_random=50)

    print("\n=== Summary ===")
    print("If the formula holds, then for any arc flip i->j to j->i:")
    print("  H(T') - H(T) = 2 * (#{odd cycles using j->i in T'} - #{odd cycles using i->j in T})")
    print("This directly connects Hamiltonian path changes to cycle count changes!")


if __name__ == "__main__":
    main()
