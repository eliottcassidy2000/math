#!/usr/bin/env python3
"""
OCF Arc-Flip Invariance Study
================================
The OCF formula says H(T) = I(Omega(T), 2).
Define E(T) = H(T) - I(Omega(T), 2).

Claim A is equivalent to E(T) = 0 for all T (given Claim B is proved).

Strategy: show E(T) is invariant under arc flips.
  - E(transitive) = 1 - 1 = 0 (trivially)
  - Any tournament is reachable from transitive by arc flips
  - If E is flip-invariant, E = 0 for all T.

This script computes delta_E = delta_H - delta_I for every arc flip,
and studies the STRUCTURE of delta_H and delta_I individually.

Instance: kind-pasteur-2026-03-05-S5
"""

import sys
sys.path.insert(0, r"C:\Users\Eliott\Documents\GitHub\math\03-artifacts\code")
from tournament_lib import *
from collections import Counter, defaultdict
import random

def flip_arc(T, i, j):
    """Return T' obtained by reversing arc between i and j."""
    n = len(T)
    Tp = [row[:] for row in T]
    Tp[i][j] = 1 - T[i][j]
    Tp[j][i] = 1 - T[j][i]
    return Tp

def compute_I(T):
    """Compute I(Omega(T), 2)."""
    cycles = find_odd_cycles(T)
    if not cycles:
        return 1
    cg = conflict_graph(cycles)
    return independence_poly_at(cg, 2)

def compute_E(T):
    """Compute E(T) = H(T) - I(Omega(T), 2)."""
    return hamiltonian_path_count(T) - compute_I(T)


def main():
    print("=== OCF Arc-Flip Invariance Study ===\n")

    # Phase 1: Verify E(T) = 0 for all n<=5
    print("Phase 1: Verify E(T) = 0 exhaustively for n <= 5")
    print("-" * 60)
    for n in range(1, 6):
        count = 0
        failures = 0
        for T in all_tournaments(n):
            e = compute_E(T)
            count += 1
            if e != 0:
                failures += 1
        print(f"  n={n}: {count} tournaments, {failures} failures, E=0 for all: {'YES' if failures == 0 else 'NO'}")

    # Phase 2: Compute delta_H and delta_I for every arc flip at n=5
    print("\nPhase 2: delta_H and delta_I under arc flips at n=5")
    print("-" * 60)

    n = 5
    delta_H_vals = []
    delta_I_vals = []
    delta_E_vals = []
    triple_counts = Counter()

    for T in all_tournaments(n):
        for i in range(n):
            for j in range(i+1, n):
                # Flip arc between i and j (whichever direction it goes)
                if T[i][j]:
                    # flip i->j to j->i
                    Tp = flip_arc(T, i, j)
                else:
                    # flip j->i to i->j
                    Tp = flip_arc(T, j, i)

                ht = hamiltonian_path_count(T)
                htp = hamiltonian_path_count(Tp)
                it = compute_I(T)
                itp = compute_I(Tp)

                dh = htp - ht
                di = itp - it
                de = dh - di  # should be 0

                delta_H_vals.append(dh)
                delta_I_vals.append(di)
                delta_E_vals.append(de)
                triple_counts[(dh, di)] += 1

    print(f"Tested {len(delta_E_vals)} (T, arc) pairs")
    print(f"delta_E = 0 for all: {all(x == 0 for x in delta_E_vals)}")
    print(f"delta_H range: [{min(delta_H_vals)}, {max(delta_H_vals)}]")
    print(f"delta_I range: [{min(delta_I_vals)}, {max(delta_I_vals)}]")
    print(f"delta_E range: [{min(delta_E_vals)}, {max(delta_E_vals)}]")

    print("\n(delta_H, delta_I) distribution:")
    for (dh, di), cnt in sorted(triple_counts.items()):
        match = "[=]" if dh == di else "[!]"
        print(f"  ({dh:+4d}, {di:+4d}) -> {cnt:6d}  {match}")

    # Phase 3: n=6 spot check
    print("\nPhase 3: n=6 random spot check (E=0 and delta_E=0)")
    print("-" * 60)

    random.seed(42)
    n = 6
    n_E_nonzero = 0
    n_dE_nonzero = 0
    n_tests = 0

    for _ in range(500):
        T = random_tournament(n, random)
        e = compute_E(T)
        if e != 0:
            n_E_nonzero += 1

        i, j = random.sample(range(n), 2)
        if T[i][j]:
            Tp = flip_arc(T, i, j)
        else:
            Tp = flip_arc(T, j, i)

        ep = compute_E(Tp)
        if e != ep:
            n_dE_nonzero += 1
        n_tests += 1

    print(f"E != 0: {n_E_nonzero}/{n_tests}")
    print(f"delta_E != 0: {n_dE_nonzero}/{n_tests}")

    # Phase 4: Understanding delta_H for a flip
    # When we flip i->j to j->i:
    # delta_H = (paths using j->i in T') - (paths using i->j in T)
    # Can we relate this to local structure around i and j?
    print("\nPhase 4: Local structure of delta_H and delta_I")
    print("-" * 60)

    # For each flip, compute:
    # - in-degree/out-degree of i and j
    # - number of common predecessors/successors
    # - delta_H and delta_I
    n = 5
    local_data = []
    for T in all_tournaments(n):
        for i in range(n):
            for j in range(i+1, n):
                if not T[i][j]:
                    continue  # only flip i->j

                Tp = flip_arc(T, i, j)
                dh = hamiltonian_path_count(Tp) - hamiltonian_path_count(T)
                di = compute_I(Tp) - compute_I(T)

                # Local structure
                # Vertices that beat both i and j
                beat_both = sum(1 for k in range(n) if k != i and k != j
                                and T[k][i] and T[k][j])
                # Vertices beaten by both
                lose_both = sum(1 for k in range(n) if k != i and k != j
                                and T[i][k] and T[j][k])
                # k beats i, j beats k (the "through" vertices for paths ...k->i->j->...)
                # After flip: ...k->j->i->... needs k->j and i->...
                pass_through = sum(1 for k in range(n) if k != i and k != j
                                   and T[k][i] and T[j][k])

                local_data.append({
                    'dh': dh, 'di': di,
                    'beat_both': beat_both,
                    'lose_both': lose_both,
                    'pass_through': pass_through,
                })

    # Check if delta_H is determined by local structure
    print("Correlating delta_H with local vertex structure...")
    # Group by (beat_both, lose_both, pass_through) and see if dh is determined
    by_local = defaultdict(list)
    for d in local_data:
        key = (d['beat_both'], d['lose_both'], d['pass_through'])
        by_local[key].append((d['dh'], d['di']))

    print(f"{'(bb,lb,pt)':>12s} | {'dh values':>30s} | {'di values':>30s}")
    for key in sorted(by_local.keys()):
        dh_vals = sorted(set(v[0] for v in by_local[key]))
        di_vals = sorted(set(v[1] for v in by_local[key]))
        det = "DET" if len(dh_vals) == 1 else "---"
        print(f"  {str(key):>10s} | {str(dh_vals):>30s} | {str(di_vals):>30s}  {det}")

    # Phase 5: Study delta_I in terms of cycle structure
    print("\nPhase 5: delta_I decomposition by cycle changes")
    print("-" * 60)

    # For a flip i->j to j->i, I(Omega(T')) - I(Omega(T)):
    # Cycles that existed in T and not in T' are those using arc i->j
    # Cycles that exist in T' but not T are those using arc j->i
    # Cycles that exist in both are unchanged
    # But the conflict graph changes because the cycle set changes

    n = 5
    cycle_change_data = []
    for T_idx, T in enumerate(all_tournaments(n)):
        if T_idx > 100:  # limit for speed
            break
        for i in range(n):
            for j in range(i+1, n):
                if not T[i][j]:
                    continue

                Tp = flip_arc(T, i, j)
                cyc_T = find_odd_cycles(T)
                cyc_Tp = find_odd_cycles(Tp)

                # Cycles using arc i->j in T
                def uses_arc(cyc, a, b):
                    L = len(cyc)
                    return any(cyc[k] == a and cyc[(k+1)%L] == b for k in range(L))

                lost = [c for c in cyc_T if uses_arc(c, i, j)]
                gained = [c for c in cyc_Tp if uses_arc(c, j, i)]
                common_T = [c for c in cyc_T if not uses_arc(c, i, j)]
                common_Tp = [c for c in cyc_Tp if not uses_arc(c, j, i)]

                # Sanity: common cycles should be the same set
                assert set(frozenset(c) for c in common_T) == set(frozenset(c) for c in common_Tp), \
                    "Common cycle sets should match"

                di = compute_I(Tp) - compute_I(T)
                n_lost = len(lost)
                n_gained = len(gained)
                lost_lengths = sorted([len(c) for c in lost])
                gained_lengths = sorted([len(c) for c in gained])

                cycle_change_data.append({
                    'di': di, 'n_lost': n_lost, 'n_gained': n_gained,
                    'lost_lengths': tuple(lost_lengths),
                    'gained_lengths': tuple(gained_lengths),
                })

    by_change = defaultdict(list)
    for d in cycle_change_data:
        key = (d['n_lost'], d['n_gained'], d['lost_lengths'], d['gained_lengths'])
        by_change[key].append(d['di'])

    print(f"{'(lost, gained, lost_lens, gained_lens)':>50s} | delta_I values")
    for key in sorted(by_change.keys()):
        di_vals = sorted(set(by_change[key]))
        cnt = len(by_change[key])
        det = "DET" if len(di_vals) == 1 else "---"
        print(f"  {str(key):>48s} | {str(di_vals):>20s} (n={cnt}) {det}")

    print("\n=== Key Question ===")
    print("Is delta_I determined by the set of lost/gained cycles alone,")
    print("or does it depend on the structure of the common cycles?")

if __name__ == "__main__":
    main()
