"""
Investigate OPEN-Q-011: Near-cancellation of two error effects at n=6.

Define for each (T, v):
  A = sum_P' #TypeII(P')               (unweighted per-path count)
  B = sum_P' sum_{TypeII at j} mu(v, P'[j], P'[j+1])  (weighted by mu)
  D = sum_{all odd C through v} mu(C)  (Claim A RHS / 2)

Claim A says: (H(T) - H(T-v)) / 2 = D.
The per-path sum gives: (H(T) - H(T-v)) / 2 = A (if per-path identity held).

So Claim A <=> A = D. But at n=6, A != D for many pairs.
However: A - B ≈ -(B - D), meaning the two error effects nearly cancel.

Can we prove A - B = -(B - D), i.e., A = 2B - D? Or more precisely,
does A - D = 0 always hold (which IS Claim A)?

Let me investigate the exact relationship.

Author: opus-2026-03-05-S2
"""

import sys
sys.path.insert(0, '03-artifacts/code')
from tournament_lib import (
    all_tournaments, hamiltonian_path_count, delete_vertex,
    find_odd_cycles, conflict_graph, independence_poly_at, mu
)
from itertools import permutations


def compute_abd(T, v, n):
    """Compute A, B, D for a (T, v) pair at size n."""
    Tv, old_labels = delete_vertex(T, v)
    tv_cycles = find_odd_cycles(Tv)
    cache = (Tv, old_labels, tv_cycles)

    # Find Ham paths of T-v
    tmv_verts = list(range(n - 1))  # vertices of Tv

    A = 0  # sum of #TypeII over all paths
    B = 0  # sum of mu-weighted TypeII over all paths

    # Precompute mu for all 3-cycles (v, old_a, old_b) where v->old_a->old_b->v
    # We need to map back to original labels
    old_to_new = {old: new for new, old in enumerate(old_labels)}

    for perm in permutations(tmv_verts):
        # Check if this is a valid Ham path of Tv
        valid = True
        for i in range(len(perm) - 1):
            if Tv[perm[i]][perm[i + 1]] != 1:
                valid = False
                break
        if not valid:
            continue

        # Compute signature with respect to v (using original labels)
        sig = []
        for j in range(len(perm)):
            old_j = old_labels[perm[j]]
            sig.append(1 if T[v][old_j] == 1 else 0)

        # Count Type-II positions and their mu weights
        for j in range(len(sig) - 1):
            if sig[j] == 1 and sig[j + 1] == 0:
                # Type-II position
                A += 1
                # The 3-cycle is v -> old_labels[perm[j]] -> old_labels[perm[j+1]] -> v
                old_a = old_labels[perm[j]]
                old_b = old_labels[perm[j + 1]]
                # Verify this IS a 3-cycle
                assert T[v][old_a] == 1 and T[old_a][old_b] == 1 and T[old_b][v] == 1
                cycle = (v, old_a, old_b)  # directed cycle as tuple
                mu_val = mu(T, v, cycle, _tv_cache=cache)
                B += mu_val

    # D = sum of mu(C) over ALL odd cycles through v
    all_cycles = find_odd_cycles(T)
    cycles_v = [c for c in all_cycles if v in set(c)]
    D = sum(mu(T, v, c, _tv_cache=cache) for c in cycles_v)

    return A, B, D


if __name__ == "__main__":
    n = 6
    print(f"OPEN-Q-011: Near-cancellation analysis at n={n}")
    print("="*70)

    results = []
    total = 0

    # Sample tournaments (full n=6 would be slow with per-path enumeration)
    for t_idx, T in enumerate(all_tournaments(n)):
        if t_idx >= 500:  # sample first 500
            break
        for v in range(n):
            total += 1
            A, B, D = compute_abd(T, v, n)
            ht = hamiltonian_path_count(T)
            Tv, _ = delete_vertex(T, v)
            htv = hamiltonian_path_count(Tv)
            half_lhs = (ht - htv) // 2

            results.append({
                'A': A, 'B': B, 'D': D,
                'half_lhs': half_lhs,
                'A-B': A - B, 'B-D': B - D, 'A-D': A - D,
                'claim_a_ok': half_lhs == D,
            })

    print(f"Sampled {total} (T, v) pairs\n")

    # Check: does A = half_lhs always? (per-path sum = insertion count)
    a_eq_lhs = sum(1 for r in results if r['A'] == r['half_lhs'])
    print(f"A = (H(T)-H(T-v))/2: {a_eq_lhs}/{total}")

    # Check: does D = half_lhs always? (Claim A)
    claim_a_ok = sum(1 for r in results if r['claim_a_ok'])
    print(f"D = (H(T)-H(T-v))/2 (Claim A): {claim_a_ok}/{total}")

    # Statistics on A-B, B-D, A-D
    ab = [r['A-B'] for r in results]
    bd = [r['B-D'] for r in results]
    ad = [r['A-D'] for r in results]

    print(f"\nA-B (mu>1 effect): mean={sum(ab)/len(ab):.4f}, "
          f"min={min(ab)}, max={max(ab)}")
    print(f"B-D (5-cycle effect): mean={sum(bd)/len(bd):.4f}, "
          f"min={min(bd)}, max={max(bd)}")
    print(f"A-D (total error): mean={sum(ad)/len(ad):.4f}, "
          f"min={min(ad)}, max={max(ad)}")

    # Key test: does A-B = -(B-D) exactly? i.e., A-D = 0?
    ad_zero = sum(1 for r in results if r['A-D'] == 0)
    print(f"\nA = D exactly: {ad_zero}/{total}")
    print(f"A-B = -(B-D) exactly (A=D): {ad_zero}/{total}")

    # If A = D always, then the per-path sum equals the full cycle sum
    # This would mean: sum_P' #TypeII = sum_{all odd C through v} mu(C)
    # which is a different way of stating Claim A!
    if ad_zero == total:
        print("\n*** A = D for ALL pairs! ***")
        print("This means: sum_P' #TypeII(P') = sum_{all C through v} mu(C)")
        print("Which IS Claim A (since sum_P' #TypeII = (H(T)-H(T-v))/2)")
    else:
        # Show the distribution of A-D values
        ad_dist = {}
        for v in ad:
            ad_dist[v] = ad_dist.get(v, 0) + 1
        print(f"\nA-D distribution: {dict(sorted(ad_dist.items()))}")

        # Check if A-B + B-D = 0 exactly (should be, since A-D = A-B + B-D)
        check = sum(1 for r in results if r['A-B'] + r['B-D'] == r['A-D'])
        print(f"A-B + B-D = A-D (sanity): {check}/{total}")

    # Decompose: when A != D, look at the structure
    print("\n--- Cases where A != D (if any) ---")
    shown = 0
    for i, r in enumerate(results):
        if r['A-D'] != 0 and shown < 5:
            print(f"  Pair {i}: A={r['A']}, B={r['B']}, D={r['D']}, "
                  f"A-B={r['A-B']}, B-D={r['B-D']}, A-D={r['A-D']}")
            shown += 1
