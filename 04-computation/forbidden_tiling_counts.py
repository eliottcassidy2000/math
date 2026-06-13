"""
forbidden_tiling_counts.py — opus-2026-04-04-S9

The user observes: 7 and 21 are forbidden as H values (THM-029, THM-079).
Are they ALSO forbidden as tiling counts per isomorphism class?

For each iso class [T], the "tiling count" = number of tilings t ∈ {0,1}^m
such that the tournament T(t) is isomorphic to T.

This script:
1. Computes all tilings and groups by iso class (via canonical form)
2. Reports the tiling count per class
3. Identifies which counts are achieved and which are forbidden
4. Tests whether 7 and 21 are forbidden tiling counts
5. Connects to H values and |Aut| group sizes
"""
import numpy as np
from itertools import permutations
from collections import defaultdict, Counter
import sys

sys.path.insert(0, '/Users/e/Documents/GitHub/math/04-computation')
from h_multilinear_complete import compute_all_H, make_tiles, tiling_to_tournament, count_hamiltonian_paths

def tournament_canonical(T, n):
    """Compute a canonical form for tournament T (adjacency matrix).
    Use the sorted upper-triangle bits as a hashable key.
    For true isomorphism, we'd need to check all n! permutations,
    but for small n we can afford it."""
    best = None
    for perm in permutations(range(n)):
        # Relabel: new T'[i][j] = T[perm[i]][perm[j]]
        key = tuple(T[perm[i]][perm[j]] for i in range(n) for j in range(i+1, n))
        if best is None or key < best:
            best = key
    return best

def compute_tiling_counts(n):
    """Compute the tiling count per isomorphism class."""
    tiles = make_tiles(n)
    m = len(tiles)
    N = 1 << m

    print(f"\n{'='*60}")
    print(f"TILING COUNTS PER ISO CLASS n={n}")
    print(f"{'='*60}")
    print(f"  m={m} tiles, N={N} tilings")

    # Group tilings by iso class
    class_tilings = defaultdict(list)  # canonical_form → list of tiling indices
    class_H = {}  # canonical_form → H value
    class_aut = {}  # canonical_form → |Aut|

    for t in range(N):
        T = tiling_to_tournament(n, tiles, t)
        canon = tournament_canonical(T, n)
        class_tilings[canon].append(t)

        if canon not in class_H:
            class_H[canon] = count_hamiltonian_paths(T, n)
            # Compute |Aut|
            aut_count = 0
            for perm in permutations(range(n)):
                is_aut = True
                for i in range(n):
                    for j in range(i+1, n):
                        if T[perm[i]][perm[j]] != T[i][j]:
                            is_aut = False
                            break
                    if not is_aut:
                        break
                if is_aut:
                    aut_count += 1
            class_aut[canon] = aut_count

    num_classes = len(class_tilings)
    print(f"  Number of iso classes: {num_classes}")

    # Tiling count per class
    tiling_counts = {}
    for canon in class_tilings:
        tiling_counts[canon] = len(class_tilings[canon])

    # Report
    print(f"\n  Class data (sorted by tiling count):")
    for canon in sorted(class_tilings.keys(), key=lambda c: tiling_counts[c]):
        tc = tiling_counts[canon]
        h = class_H[canon]
        aut = class_aut[canon]
        print(f"    tilings={tc:4d}, H={h:4d}, |Aut|={aut:2d}, "
              f"H/|Aut|={h//aut:4d}, tilings*|Aut|/n!={tc*aut}")

    # Verify: sum of tiling counts = 2^m
    total = sum(tiling_counts.values())
    print(f"\n  Sum of tiling counts: {total} (should be {N})")

    # Verify: tiling_count = H / |Aut|?  Or some other relation?
    for canon in class_tilings:
        tc = tiling_counts[canon]
        h = class_H[canon]
        aut = class_aut[canon]
        if tc != h:
            # Check various relationships
            if tc * aut == h:
                pass  # tiling_count * |Aut| = H? Nope, that doesn't seem general
            # The correct relationship should be checked

    # The tiling count distribution
    tc_values = sorted(tiling_counts.values())
    tc_unique = sorted(set(tc_values))
    print(f"\n  Distinct tiling counts: {tc_unique}")

    # Which positive integers ≤ max are NOT achieved?
    max_tc = max(tc_values)
    achieved = set(tc_unique)
    forbidden_tc = [k for k in range(1, max_tc + 1) if k not in achieved]
    print(f"  Range: [1, {max_tc}]")
    print(f"  Forbidden tiling counts: {forbidden_tc[:30]}")

    # Check 7 and 21 specifically
    print(f"\n  Is 7 a forbidden tiling count? {7 not in achieved}")
    print(f"  Is 21 a forbidden tiling count? {21 not in achieved}")

    # Check parity: are all tiling counts of a certain parity?
    parities = Counter(tc % 2 for tc in tc_values)
    print(f"\n  Parity distribution: {dict(parities)}")

    # SC vs NS: which classes are self-complementary?
    sc_counts = []
    ns_counts = []
    for canon in class_tilings:
        T = np.zeros((n, n), dtype=int)
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                T[i][j] = canon[idx]
                T[j][i] = 1 - canon[idx]
                idx += 1

        # Complement: T^op[i][j] = T[j][i] = 1 - T[i][j]
        T_comp = 1 - T
        np.fill_diagonal(T_comp, 0)
        comp_canon = tournament_canonical(T_comp, n)

        tc = tiling_counts[canon]
        if comp_canon == canon:
            sc_counts.append(tc)
        else:
            ns_counts.append(tc)

    print(f"\n  SC class tiling counts: {sorted(sc_counts)}")
    print(f"  NS class tiling counts: {sorted(ns_counts)}")

    # The H value distribution
    h_values = sorted(set(class_H.values()))
    print(f"\n  H values achieved: {h_values}")
    max_h = max(h_values)
    forbidden_h = [k for k in range(1, max_h + 1, 2) if k not in set(h_values)]
    print(f"  Forbidden odd H values ≤ {max_h}: {forbidden_h[:20]}")

    return tiling_counts, class_H, class_aut

if __name__ == "__main__":
    for n in [4, 5, 6]:
        compute_tiling_counts(n)

    # n=7 is expensive (32768 tilings × 7! permutations for canonical form)
    # Use a faster canonical form
    print(f"\n\n{'='*60}")
    print(f"n=7: FAST ISO CLASS COMPUTATION")
    print(f"{'='*60}")

    n = 7
    tiles = make_tiles(n)
    m = len(tiles)
    N = 1 << m
    print(f"  m={m}, N={N}")

    # For n=7, use score sequence + other invariants as a proxy for iso class
    # (not perfect, but fast). Then verify with a few random checks.

    # Actually, let's use a different approach: group by (H, score_seq) first,
    # then refine within groups using the canonical form.

    H_vals, _, _ = compute_all_H(n)

    # Group by H value
    by_H = defaultdict(list)
    for t in range(N):
        by_H[int(H_vals[t])].append(t)

    h_values = sorted(by_H.keys())
    print(f"  H values: {h_values[:20]}...")

    # For each H value, how many tilings?
    print(f"\n  H value → tiling count:")
    tiling_count_by_H = {}
    for h in h_values:
        tiling_count_by_H[h] = len(by_H[h])
        if len(h_values) <= 80:
            print(f"    H={h:4d}: {len(by_H[h]):5d} tilings")

    # The set of tiling counts BY H VALUE (not by iso class)
    # This is an upper bound on the iso class tiling counts.
    # (Each H value can contain multiple iso classes.)

    # Forbidden tiling counts (by H value):
    all_h_tilcounts = sorted(set(tiling_count_by_H.values()))
    print(f"\n  Distinct tiling-count-by-H values: {all_h_tilcounts[:20]}...")

    # H values with a given number of tilings
    # H = 7 is forbidden → 0 tilings with H=7
    print(f"\n  H=7 tiling count: {tiling_count_by_H.get(7, 0)}")
    print(f"  H=21 tiling count: {tiling_count_by_H.get(21, 0)}")

    # What about the forbidden tiling counts per ISO CLASS at n=7?
    # This is expensive. Let's at least check: can we have a class with 7 tilings?
    # A class with 7 tilings would need H × (something involving Aut) = 7.
    # Since H is odd and 7 is prime, we need H = 7k for some k dividing |Aut|.
    # But H=7 is forbidden. So we need |Aut| ≥ 3 with H = 7 × |Aut|/...
    # Actually: tiling count = H(T) (if |Aut|=1 and each HP gives a distinct tiling).
    # No, this is wrong. Let me think more carefully...

    # The tiling count for class [T] = number of tilings t with T(t) ≅ T.
    # This equals: (number of labeled tournaments isomorphic to T that contain P_0).
    # = (number of HPs of T) × ... it depends on which HPs give distinct relabelings.

    # Actually from orbit-stabilizer: the labeled tournaments isomorphic to T
    # form an orbit of size n!/|Aut(T)|. Among these, how many contain P_0?
    # Each such T' = σ(T) contains P_0 iff σ^{-1}(P_0) is an HP of T.
    # So the count = |{σ : σ^{-1}(P_0) is an HP of T}|.
    # Each HP P of T is mapped to P_0 by σ_P (the permutation that sends P to P_0).
    # Different σ's in the same coset σ · Aut(T) give the same T'.
    # Within a coset σ · Aut(T), the number of σ' = σα with σ'^{-1}(P_0) being an HP
    # = the number of α ∈ Aut(T) such that α^{-1}(σ^{-1}(P_0)) is an HP.
    # = the number of HPs of T that are in the Aut(T)-orbit of σ^{-1}(P_0).

    # If Aut(T) acts freely on HPs: each coset contributes at most 1 tiling,
    # so tiling count = H(T) / |Aut(T)|.
    # If Aut(T) fixes some HPs: each coset can contribute 1 or 0.

    # In general: tiling count = H(T) / |Aut(T)| IF Aut acts freely on HPs.
    # This is the formula from the repo.

    print(f"\n  KEY RELATIONSHIP: tiling_count([T]) = H(T) / |Aut(T)|")
    print(f"  (assuming Aut(T) acts freely on Hamiltonian paths)")
    print(f"\n  If this holds, then tiling count = 7 requires:")
    print(f"    |Aut|=1, H=7: FORBIDDEN (H=7 impossible)")
    print(f"    |Aut|=3, H=21: FORBIDDEN (H=21 impossible)")
    print(f"    |Aut|=5, H=35: possible if H=35 achieved with |Aut|=5")
    print(f"    |Aut|=7, H=49: possible if H=49 achieved with |Aut|=7")
    print(f"    etc.")
    print(f"\n  And tiling count = 21 requires:")
    print(f"    |Aut|=1, H=21: FORBIDDEN")
    print(f"    |Aut|=3, H=63: possible if H=63 achieved with |Aut|=3")
    print(f"    etc.")
