#!/usr/bin/env python3
"""
Verify the correct algebraic formula for delta_I under arc flip.

When flipping arc i->j to j->i:
  Lost cycles: those using arc i->j (all contain both i and j, pairwise conflicting)
  Gained cycles: those using arc j->i (same structure)
  Common cycles: all others (unchanged)

Since all lost/gained cycles are pairwise adjacent (share i and j),
at most one can appear in any independent set. By inclusion-exclusion:

  I(Omega(T), 2) = I(Omega_common, 2) + sum_{lost c} 2 * I(Omega_common \ N(c), 2)
  I(Omega(T'), 2) = I(Omega_common, 2) + sum_{gained c} 2 * I(Omega_common \ N(c), 2)

  delta_I = 2 * [sum_{gained c} I(R_c, 2) - sum_{lost c} I(R_c, 2)]

where R_c = conflict graph of common cycles NOT adjacent to c.

At n<=5, R_c is always empty (I=1), giving the simple formula.
At n>=6, R_c can be non-empty.

ALSO: verify that delta_H = delta_I using this decomposition.

Instance: kind-pasteur-2026-03-05-S5
"""

import sys
sys.path.insert(0, r"C:\Users\Eliott\Documents\GitHub\math\03-artifacts\code")
from tournament_lib import *
import random


def flip_arc(T, i, j):
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


def compute_delta_I_formula(T, i, j):
    """Compute delta_I using the algebraic decomposition.
    T has arc i->j. We compute I(Omega(T')) - I(Omega(T)) where T' flips i->j.

    Returns (delta_I_formula, delta_I_direct, details).
    """
    Tp = flip_arc(T, i, j)
    cyc_T = find_odd_cycles(T)
    cyc_Tp = find_odd_cycles(Tp)

    # Classify cycles
    lost = [c for c in cyc_T if uses_arc(c, i, j)]     # in T only
    gained = [c for c in cyc_Tp if uses_arc(c, j, i)]   # in T' only
    common = [c for c in cyc_T if not uses_arc(c, i, j)] # in both

    # Sanity: common should match
    common_set_T = {frozenset(c) for c in common}
    common_set_Tp = {frozenset(c) for c in cyc_Tp if not uses_arc(c, j, i)}
    assert common_set_T == common_set_Tp, "Common cycle mismatch"

    # Build conflict graph on common cycles
    common_vsets = [set(c) for c in common]

    # For each lost/gained cycle c, find R_c = common cycles not conflicting with c
    def compute_Rc_value(c):
        """I(Omega_common \ N(c), 2) where N(c) = common cycles sharing a vertex with c."""
        c_vset = set(c)
        # Indices of common cycles NOT conflicting with c
        free_idx = [idx for idx, vs in enumerate(common_vsets) if not (vs & c_vset)]
        if not free_idx:
            return 1  # empty graph

        # Build conflict graph on free cycles
        free_cycles = [common[idx] for idx in free_idx]
        cg = conflict_graph(free_cycles)
        return independence_poly_at(cg, 2)

    # Compute formula
    gained_sum = sum(compute_Rc_value(c) for c in gained)
    lost_sum = sum(compute_Rc_value(c) for c in lost)
    delta_I_formula = 2 * (gained_sum - lost_sum)

    # Direct computation
    it = compute_I(T)
    itp = compute_I(Tp)
    delta_I_direct = itp - it

    # Also compute delta_H
    ht = hamiltonian_path_count(T)
    htp = hamiltonian_path_count(Tp)
    delta_H = htp - ht

    # Compute individual R_c values for debugging
    lost_Rc = [compute_Rc_value(c) for c in lost]
    gained_Rc = [compute_Rc_value(c) for c in gained]

    return {
        'delta_I_formula': delta_I_formula,
        'delta_I_direct': delta_I_direct,
        'delta_H': delta_H,
        'n_lost': len(lost),
        'n_gained': len(gained),
        'n_common': len(common),
        'lost_Rc': lost_Rc,
        'gained_Rc': gained_Rc,
        'formula_ok': delta_I_formula == delta_I_direct,
        'H_eq_I': delta_H == delta_I_direct,
    }


def main():
    print("=== Correct delta_I formula verification ===\n")

    rng = random.Random(42)

    for n in [5, 6]:
        print(f"--- n={n} ---")
        formula_ok = 0
        H_eq_I_ok = 0
        total = 0
        n_trials = 200 if n <= 5 else 100

        Rc_values_seen = set()

        for trial in range(n_trials):
            if n <= 5:
                T = list(all_tournaments(n))[trial % (1 << (n*(n-1)//2))]
            else:
                T = random_tournament(n, rng)

            # Pick random arc
            i, j = rng.sample(range(n), 2)
            if not T[i][j]:
                i, j = j, i

            r = compute_delta_I_formula(T, i, j)
            total += 1

            if r['formula_ok']:
                formula_ok += 1
            else:
                print(f"  FORMULA FAIL: dI_formula={r['delta_I_formula']}, "
                      f"dI_direct={r['delta_I_direct']}")

            if r['H_eq_I']:
                H_eq_I_ok += 1

            # Track R_c values
            for v in r['lost_Rc'] + r['gained_Rc']:
                Rc_values_seen.add(v)

        print(f"  formula delta_I: {formula_ok}/{total} match")
        print(f"  delta_H = delta_I: {H_eq_I_ok}/{total}")
        print(f"  R_c values seen: {sorted(Rc_values_seen)}")

    # Detailed analysis at n=6
    print("\n--- Detailed n=6 analysis ---")
    rng = random.Random(123)
    for trial in range(20):
        T = random_tournament(6, rng)
        i, j = rng.sample(range(6), 2)
        if not T[i][j]:
            i, j = j, i

        r = compute_delta_I_formula(T, i, j)
        rc_str = f"lost_Rc={r['lost_Rc']}, gained_Rc={r['gained_Rc']}"
        match = "OK" if r['formula_ok'] and r['H_eq_I'] else "FAIL"
        print(f"  flip {i}->{j}: dH={r['delta_H']:+3d}, dI={r['delta_I_direct']:+3d}, "
              f"formula={r['delta_I_formula']:+3d}, "
              f"lost={r['n_lost']}, gained={r['n_gained']}, common={r['n_common']} "
              f"{rc_str} [{match}]")


    # Key question: can we decompose delta_H similarly?
    print("\n=== Key Question ===")
    print("The algebraic formula correctly gives delta_I.")
    print("Does delta_H = delta_I provide a path to proof?")
    print("If H(T) = I(Omega(T), 2), then delta_H = delta_I is a tautology.")
    print("But if we can independently prove delta_H = delta_I for arc flips,")
    print("then since both equal 0 for transitive tournaments, H = I(Omega, 2) follows.")


if __name__ == "__main__":
    main()
