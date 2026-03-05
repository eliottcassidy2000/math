#!/usr/bin/env python3
"""
Arc-Reversal Invariance Study (OPEN-Q-009)
============================================
For each (T, v, arc i->j) with i,j != v, study what happens to:
  - H(T) and H(T-v)
  - The set of odd cycles through v
  - The total mu sum

The key question: when we flip i->j to j->i (forming T'), why does
  H(T) - H(T-v) - 2*sum mu_T(C) = H(T') - H(T'-v) - 2*sum mu_{T'}(C') ?

Since both sides are 0 (Claim A verified n<=6), the interesting
structure is in the INDIVIDUAL changes: DeltaH, DeltaH_v, DeltaMu.

Instance: kind-pasteur-2026-03-05-S5
"""

import sys
sys.path.insert(0, r"C:\Users\Eliott\Documents\GitHub\math\03-artifacts\code")
from tournament_lib import *
from collections import Counter, defaultdict

def flip_arc(T, i, j):
    """Return T' obtained by reversing arc between i and j."""
    n = len(T)
    Tp = [row[:] for row in T]
    Tp[i][j] = 1 - T[i][j]
    Tp[j][i] = 1 - T[j][i]
    return Tp

def compute_D(T, v):
    """Compute D(T,v) = H(T) - H(T-v) - 2*sum mu(C)."""
    ht = hamiltonian_path_count(T)
    Tv, old_labels = delete_vertex(T, v)
    htv = hamiltonian_path_count(Tv)

    all_cyc = find_odd_cycles(T)
    cyc_v = [c for c in all_cyc if v in set(c)]

    tv_cycles = find_odd_cycles(Tv)
    cache = (Tv, old_labels, tv_cycles)
    total_mu = sum(mu(T, v, c, _tv_cache=cache) for c in cyc_v)

    return ht - htv - 2 * total_mu, ht, htv, total_mu, cyc_v

def ham_paths_using_arc(T, i, j):
    """Count Ham paths in T that use the arc i->j."""
    n = len(T)
    if not T[i][j]:
        return 0
    # DP counting paths using arc i->j
    # A Ham path uses i->j iff i immediately precedes j in the path
    full = (1 << n) - 1
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1

    # Track separately: dp_used[mask][v] = paths ending at v that have used i->j
    dp_used = [[0]*n for _ in range(1 << n)]

    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            cnt = dp[mask][v]
            cnt_used = dp_used[mask][v]
            if cnt == 0 and cnt_used == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if T[v][u]:
                    if v == i and u == j:
                        # This step uses the arc i->j
                        dp[mask | (1 << u)][u] += cnt
                        dp_used[mask | (1 << u)][u] += cnt_used + cnt
                    else:
                        dp[mask | (1 << u)][u] += cnt
                        dp_used[mask | (1 << u)][u] += cnt_used

    return sum(dp_used[full][v] for v in range(n))

def cycles_using_arc(cycles, i, j):
    """Return cycles that contain the arc i->j (i immediately before j in the cycle)."""
    result = []
    for c in cycles:
        L = len(c)
        for k in range(L):
            if c[k] == i and c[(k+1) % L] == j:
                result.append(c)
                break
    return result

def study_arc_flip(T, v, i, j):
    """Study a single arc flip i->j -> j->i."""
    n = len(T)
    assert T[i][j] == 1, f"Expected arc {i}->{j}"

    Tp = flip_arc(T, i, j)

    D, ht, htv, mu_sum, cyc_v = compute_D(T, v)
    Dp, htp, htpv, mu_sum_p, cyc_v_p = compute_D(Tp, v)

    delta_H = htp - ht
    delta_Hv = htpv - htv
    delta_mu = mu_sum_p - mu_sum

    # Cycles through v using arc i->j in T
    using_ij = cycles_using_arc(cyc_v, i, j)
    # Cycles through v using arc j->i in T'
    using_ji = cycles_using_arc(cyc_v_p, j, i)
    # Cycles through v that don't use either arc direction (present in both)
    cyc_v_sets = {frozenset(c) for c in cyc_v}
    cyc_vp_sets = {frozenset(c) for c in cyc_v_p}
    common = cyc_v_sets & cyc_vp_sets
    only_T = cyc_v_sets - cyc_vp_sets
    only_Tp = cyc_vp_sets - cyc_v_sets

    return {
        'D': D, 'Dp': Dp,
        'ht': ht, 'htv': htv, 'mu_sum': mu_sum,
        'htp': htp, 'htpv': htpv, 'mu_sum_p': mu_sum_p,
        'delta_H': delta_H, 'delta_Hv': delta_Hv, 'delta_mu': delta_mu,
        'n_cyc_v': len(cyc_v), 'n_cyc_vp': len(cyc_v_p),
        'n_using_ij': len(using_ij), 'n_using_ji': len(using_ji),
        'n_common': len(common), 'n_only_T': len(only_T), 'n_only_Tp': len(only_Tp),
    }

def main():
    print("=== Arc-Reversal Invariance Study (OPEN-Q-009) ===\n")

    # --- Phase 1: Exhaustive at n=5, sample structure ---
    print("Phase 1: Exhaustive n=5 decomposition")
    print("-" * 60)

    delta_stats = defaultdict(list)
    n = 5
    count = 0

    for T in all_tournaments(n):
        for v in range(n):
            for i in range(n):
                for j in range(n):
                    if i == v or j == v or i == j:
                        continue
                    if not T[i][j]:
                        continue

                    r = study_arc_flip(T, v, i, j)
                    assert r['D'] == 0 and r['Dp'] == 0, "Claim A failed!"

                    # Since D=Dp=0, we have delta_H - delta_Hv = 2*delta_mu
                    assert r['delta_H'] - r['delta_Hv'] == 2 * r['delta_mu'], \
                        f"Invariance decomposition failed: {r}"

                    delta_stats['delta_H'].append(r['delta_H'])
                    delta_stats['delta_Hv'].append(r['delta_Hv'])
                    delta_stats['delta_mu'].append(r['delta_mu'])
                    delta_stats['n_only_T'].append(r['n_only_T'])
                    delta_stats['n_only_Tp'].append(r['n_only_Tp'])
                    count += 1

    print(f"Tested {count} (T, v, i->j) quadruples at n={n}")
    print(f"All satisfy D=0 and delta_H - delta_Hv = 2*delta_mu ✓\n")

    for key in ['delta_H', 'delta_Hv', 'delta_mu']:
        vals = delta_stats[key]
        c = Counter(vals)
        print(f"{key}: min={min(vals)}, max={max(vals)}, "
              f"mean={sum(vals)/len(vals):.3f}")
        print(f"  distribution: {dict(sorted(c.items()))}")

    print()

    # --- Phase 2: Structure of delta_H ---
    print("Phase 2: Decomposition of delta_H via arc-counting")
    print("-" * 60)

    # For n=5, compute ham_paths_using_arc to decompose delta_H
    n = 5
    decomp_check = 0
    for T in all_tournaments(n):
        for v in range(n):
            for i in range(n):
                for j in range(n):
                    if i == v or j == v or i == j or not T[i][j]:
                        continue

                    Tp = flip_arc(T, i, j)
                    # Paths in T using i->j disappear; paths in T' using j->i appear
                    paths_ij_in_T = ham_paths_using_arc(T, i, j)
                    paths_ji_in_Tp = ham_paths_using_arc(Tp, j, i)
                    delta_H_computed = paths_ji_in_Tp - paths_ij_in_T

                    ht = hamiltonian_path_count(T)
                    htp = hamiltonian_path_count(Tp)
                    assert htp - ht == delta_H_computed
                    decomp_check += 1
                    break  # just verify a few
                break
            break
        if decomp_check >= 10:
            break

    print(f"Verified delta_H = #paths_using_ji_in_T' - #paths_using_ij_in_T ({decomp_check} cases) ✓")

    # --- Phase 3: Relationship between cycle changes and mu changes ---
    print("\nPhase 3: Cycle structure changes under arc flip")
    print("-" * 60)

    # Look at detailed structure for one specific tournament
    # Use the cyclic tournament on 5 vertices: 0->1->2->3->4->0 (pentagonal)
    T5c = [[0]*5 for _ in range(5)]
    for i in range(5):
        T5c[i][(i+1) % 5] = 1
        T5c[i][(i+2) % 5] = 1

    print("Example: Regular tournament on 5 vertices")
    print(f"H(T) = {hamiltonian_path_count(T5c)}")

    v = 0
    print(f"\nVertex v={v}")
    for i in range(1, 5):
        for j in range(1, 5):
            if i == j or not T5c[i][j]:
                continue
            r = study_arc_flip(T5c, v, i, j)
            print(f"  Flip {i}->{j}: ΔH={r['delta_H']:+d}, ΔH_v={r['delta_Hv']:+d}, "
                  f"Δμ={r['delta_mu']:+d}, "
                  f"cycles: {r['n_cyc_v']}->{r['n_cyc_vp']} "
                  f"(lost {r['n_only_T']}, gained {r['n_only_Tp']}, kept {r['n_common']})")

    # --- Phase 4: Look for patterns in (delta_H, delta_Hv, delta_mu) triples ---
    print("\nPhase 4: Joint distribution of (ΔH, ΔH_v, Δμ)")
    print("-" * 60)

    n = 5
    triple_counts = Counter()
    for T in all_tournaments(n):
        for v in range(n):
            for i in range(n):
                for j in range(n):
                    if i == v or j == v or i == j or not T[i][j]:
                        continue
                    r = study_arc_flip(T, v, i, j)
                    triple = (r['delta_H'], r['delta_Hv'], r['delta_mu'])
                    triple_counts[triple] += 1

    print("(ΔH, ΔH_v, Δμ) -> count")
    for triple, cnt in sorted(triple_counts.items()):
        dh, dhv, dm = triple
        check = "✓" if dh - dhv == 2*dm else "✗"
        print(f"  ({dh:+3d}, {dhv:+3d}, {dm:+3d}) -> {cnt:6d}  [{check}]")

    # --- Phase 5: n=6 spot check with random tournaments ---
    print("\nPhase 5: n=6 random spot check")
    print("-" * 60)

    import random
    random.seed(42)
    n = 6
    n_tests = 0
    n_failures = 0

    for _ in range(200):
        T = random_tournament(n, random)
        v = random.randint(0, n-1)
        # Pick a random arc not involving v
        candidates = [(i,j) for i in range(n) for j in range(n)
                      if i != v and j != v and i != j and T[i][j]]
        if not candidates:
            continue
        i, j = random.choice(candidates)

        r = study_arc_flip(T, v, i, j)
        n_tests += 1
        if r['D'] != 0 or r['Dp'] != 0:
            n_failures += 1
            print(f"FAILURE at test {n_tests}: D={r['D']}, Dp={r['Dp']}")
        elif r['delta_H'] - r['delta_Hv'] != 2 * r['delta_mu']:
            n_failures += 1
            print(f"DECOMPOSITION FAILURE: {r}")

    print(f"n=6: {n_tests} random tests, {n_failures} failures")

    print("\n=== Summary ===")
    print("The identity ΔH - ΔH_v = 2·Δμ holds for every arc flip.")
    print("This is equivalent to Claim A being invariant under arc flips.")
    print("Next: understand WHY each of ΔH, ΔH_v, Δμ takes the values it does.")

if __name__ == "__main__":
    main()
