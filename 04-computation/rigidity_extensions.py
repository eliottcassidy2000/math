"""
rigidity_extensions.py — opus-2026-04-04-S11

PRACTICAL EXTENSIONS of the rigidity/forbidden-value framework.

1. Full alpha_1 gap landscape (n=3..7, sample n=8)
2. (alpha_1, alpha_2) achievable pairs — the 2D gap structure
3. H-gap census with "first appearance" table
4. Alpha_k distributions — shape, moments
5. The 21 mechanism: which (alpha_1, alpha_2) decompositions of H=21 fail?
"""
import numpy as np
from itertools import combinations, permutations
from collections import Counter, defaultdict
import sys, time

sys.path.insert(0, '/Users/e/Documents/GitHub/math/04-computation')
from h_multilinear_complete import compute_all_H, make_tiles, tiling_to_tournament, count_hamiltonian_paths

def find_odd_cycles_and_independence(T, n):
    """Find all directed odd cycles and compute independence polynomial data.
    Returns (alpha_1, alpha_2, cycle_list)."""
    cycles = []

    # 3-cycles
    for a in range(n):
        for b in range(a+1, n):
            for c in range(b+1, n):
                if T[a][b] + T[b][c] + T[c][a] == 3:
                    cycles.append(frozenset({a, b, c}))
                elif T[a][c] + T[c][b] + T[b][a] == 3:
                    cycles.append(frozenset({a, b, c}))

    # 5-cycles (only for n >= 5)
    if n >= 5:
        for verts in combinations(range(n), 5):
            for perm in permutations(verts):
                rotations = [perm[i:]+perm[:i] for i in range(5)]
                if perm != min(rotations):
                    continue
                if all(T[perm[i]][perm[(i+1)%5]] for i in range(5)):
                    cycles.append(frozenset(verts))

    # 7-cycles (only for n >= 7, expensive)
    if n >= 7 and n <= 7:
        for verts in combinations(range(n), 7):
            for perm in permutations(verts):
                rotations = [perm[i:]+perm[:i] for i in range(7)]
                if perm != min(rotations):
                    continue
                if all(T[perm[i]][perm[(i+1)%7]] for i in range(7)):
                    cycles.append(frozenset(verts))
                    break  # one orientation found

    alpha_1 = len(cycles)

    # Compute alpha_2: count vertex-disjoint pairs
    alpha_2 = 0
    for i in range(len(cycles)):
        for j in range(i+1, len(cycles)):
            if not (cycles[i] & cycles[j]):
                alpha_2 += 1

    return alpha_1, alpha_2, cycles

# ==============================================================
# 1. ALPHA_1 GAP LANDSCAPE
# ==============================================================
def alpha1_gap_landscape():
    """Map achievable alpha_1 values at each n."""
    print(f"{'='*60}")
    print(f"ALPHA_1 GAP LANDSCAPE")
    print(f"{'='*60}")

    for n in range(3, 8):
        tiles = make_tiles(n)
        m = len(tiles)
        N = 1 << m

        alpha1_set = set()
        alpha1_dist = Counter()

        t0 = time.time()
        for t in range(N):
            T = tiling_to_tournament(n, tiles, t)
            a1, _, _ = find_odd_cycles_and_independence(T, n)
            alpha1_set.add(a1)
            alpha1_dist[a1] += 1

        elapsed = time.time() - t0
        alpha1_sorted = sorted(alpha1_set)
        gaps = [a for a in range(max(alpha1_sorted)+1) if a not in alpha1_set]

        print(f"\n  n={n} ({N} tournaments, {elapsed:.1f}s):")
        print(f"    Achievable alpha_1: {alpha1_sorted}")
        print(f"    Gaps: {gaps}")
        print(f"    Distribution: {dict(sorted(alpha1_dist.items()))}")

        # The H values implied by the gaps
        # H = 1 + 2*alpha_1 + 4*alpha_2 + ...
        # If alpha_1=k is a gap, then H = 1+2k is AT LEAST gap-adjacent
        gap_H_implications = [1 + 2*g for g in gaps]
        print(f"    Implied H gap candidates: {gap_H_implications}")

# ==============================================================
# 2. (ALPHA_1, ALPHA_2) ACHIEVABLE PAIRS
# ==============================================================
def alpha12_pairs():
    """Compute all achievable (alpha_1, alpha_2) pairs."""
    print(f"\n{'='*60}")
    print(f"(ALPHA_1, ALPHA_2) ACHIEVABLE PAIRS")
    print(f"{'='*60}")

    for n in [5, 6]:
        tiles = make_tiles(n)
        m = len(tiles)
        N = 1 << m

        pairs = set()
        pair_to_H = defaultdict(set)

        for t in range(N):
            T = tiling_to_tournament(n, tiles, t)
            a1, a2, _ = find_odd_cycles_and_independence(T, n)
            h = count_hamiltonian_paths(T, n)
            pairs.add((a1, a2))
            pair_to_H[(a1, a2)].add(int(h))

        print(f"\n  n={n}:")
        print(f"    Achievable (alpha_1, alpha_2) pairs: {len(pairs)}")

        # Print as a grid
        max_a1 = max(a1 for a1, _ in pairs)
        max_a2 = max(a2 for _, a2 in pairs)

        print(f"    Grid (alpha_1 = rows, alpha_2 = cols):")
        header = "    a1\\a2 " + " ".join(f"{a2:4d}" for a2 in range(max_a2+1))
        print(header)
        for a1 in range(max_a1+1):
            row = f"    {a1:4d}  "
            for a2 in range(max_a2+1):
                if (a1, a2) in pairs:
                    h_vals = pair_to_H[(a1, a2)]
                    h_str = str(min(h_vals))
                    row += f" {h_str:>4s}"
                else:
                    row += "    ."
            print(row)

        # Which (alpha_1, alpha_2) would give H=7?
        # H = 1 + 2*a1 + 4*a2 + 8*a3 + ...
        # H = 7 → 2*a1 + 4*a2 + ... = 6 → a1 + 2*a2 + ... = 3
        # Possibilities: (3,0), (1,1)
        print(f"\n    H=7 requires (a1,a2,...) with a1+2*a2+...=3:")
        print(f"      (3,0): achieved? {(3,0) in pairs}")
        print(f"      (1,1): achieved? {(1,1) in pairs}")
        if (1,1) in pairs:
            print(f"        H at (1,1): {pair_to_H[(1,1)]}")

        # H=21 requires a1+2*a2+...=10
        print(f"\n    H=21 requires a1+2*a2+...=10:")
        target_combos = [(a1, a2) for a1 in range(11) for a2 in range((10-a1)//2+1)
                        if a1 + 2*a2 == 10]
        for a1, a2 in target_combos:
            achieved = (a1, a2) in pairs
            h_at = pair_to_H.get((a1, a2), set())
            status = f"H={h_at}" if achieved else "NOT ACHIEVED"
            # Even if achieved, H might not be 21 (higher alpha_3 could contribute)
            has_21 = 21 in h_at if achieved else False
            print(f"      ({a1},{a2}): {status}, H=21? {has_21}")

# ==============================================================
# 3. H-GAP CENSUS: FIRST APPEARANCE TABLE
# ==============================================================
def h_gap_census():
    """For each odd number up to ~200, find the smallest n where it's first achieved as H."""
    print(f"\n{'='*60}")
    print(f"H-GAP CENSUS: FIRST APPEARANCE")
    print(f"{'='*60}")

    first_appearance = {}  # odd value → smallest n

    for n in range(3, 8):
        tiles = make_tiles(n)
        m = len(tiles)
        H, _, _ = compute_all_H(n)
        achieved = set(int(h) for h in H)

        for h in achieved:
            if h not in first_appearance:
                first_appearance[h] = n

    # Print table
    max_h = 200
    print(f"\n  First appearance of odd H values (up to {max_h}):")
    print(f"  {'H':>5s} {'n':>3s}  {'H':>5s} {'n':>3s}  {'H':>5s} {'n':>3s}  {'H':>5s} {'n':>3s}")

    odd_vals = list(range(1, max_h+1, 2))
    for i in range(0, len(odd_vals), 4):
        row = "  "
        for j in range(4):
            if i+j < len(odd_vals):
                h = odd_vals[i+j]
                n_val = first_appearance.get(h, "?")
                if n_val == "?":
                    if h in [7, 21]:
                        n_val = "∞"  # permanently forbidden
                    else:
                        n_val = "≥8"
                row += f"{h:5d} {str(n_val):>3s}  "
        print(row)

    # Summary statistics
    by_first_n = defaultdict(list)
    for h, n_val in first_appearance.items():
        if h <= max_h:
            by_first_n[n_val].append(h)

    print(f"\n  Summary:")
    for n_val in sorted(by_first_n.keys()):
        vals = sorted(by_first_n[n_val])
        print(f"    First at n={n_val}: {len(vals)} values, range [{min(vals)}, {max(vals)}]")

    not_yet = [h for h in range(1, max_h+1, 2) if h not in first_appearance]
    print(f"    Not yet achieved (≤{max_h}): {not_yet}")
    print(f"    Permanently forbidden: {[h for h in not_yet if h in [7, 21]]}")
    print(f"    Likely transient: {[h for h in not_yet if h not in [7, 21]]}")

# ==============================================================
# 4. ALPHA_K DISTRIBUTIONS
# ==============================================================
def alpha_distributions():
    """Distribution of alpha_1, alpha_2 over uniform random tournaments."""
    print(f"\n{'='*60}")
    print(f"ALPHA_K DISTRIBUTIONS")
    print(f"{'='*60}")

    for n in [5, 6]:
        tiles = make_tiles(n)
        m = len(tiles)
        N = 1 << m

        a1_vals = []
        a2_vals = []
        h_vals = []

        for t in range(N):
            T = tiling_to_tournament(n, tiles, t)
            a1, a2, _ = find_odd_cycles_and_independence(T, n)
            h = count_hamiltonian_paths(T, n)
            a1_vals.append(a1)
            a2_vals.append(a2)
            h_vals.append(int(h))

        a1_arr = np.array(a1_vals, dtype=float)
        a2_arr = np.array(a2_vals, dtype=float)
        h_arr = np.array(h_vals, dtype=float)

        print(f"\n  n={n} ({N} tournaments):")
        print(f"    alpha_1: mean={np.mean(a1_arr):.2f}, std={np.std(a1_arr):.2f}, "
              f"min={int(np.min(a1_arr))}, max={int(np.max(a1_arr))}")
        print(f"    alpha_2: mean={np.mean(a2_arr):.2f}, std={np.std(a2_arr):.2f}, "
              f"min={int(np.min(a2_arr))}, max={int(np.max(a2_arr))}")

        # Skewness and kurtosis of alpha_1
        a1_centered = a1_arr - np.mean(a1_arr)
        skew = np.mean(a1_centered**3) / np.std(a1_arr)**3
        kurt = np.mean(a1_centered**4) / np.std(a1_arr)**4 - 3
        print(f"    alpha_1 skewness: {skew:.3f}, excess kurtosis: {kurt:.3f}")

        # Correlation between alpha_1 and alpha_2
        if np.std(a2_arr) > 0:
            corr_12 = np.corrcoef(a1_arr, a2_arr)[0,1]
            print(f"    Corr(alpha_1, alpha_2) = {corr_12:.4f}")

        # How well does H = 1 + 2*alpha_1 + 4*alpha_2 predict actual H?
        h_pred = 1 + 2*a1_arr + 4*a2_arr
        r2 = 1 - np.sum((h_arr - h_pred)**2) / np.sum((h_arr - np.mean(h_arr))**2)
        print(f"    H ≈ 1 + 2*a1 + 4*a2: R² = {r2:.6f}")

        max_err = np.max(np.abs(h_arr - h_pred))
        print(f"    Max |H - (1+2a1+4a2)| = {max_err}")
        if max_err > 0:
            # The excess is 8*alpha_3 + ...
            excess = h_arr - h_pred
            print(f"    Excess (= 8*a3 + 16*a4 + ...): "
                  f"nonzero for {np.sum(excess != 0)}/{N} tournaments")
            if np.sum(excess != 0) > 0:
                print(f"    Excess values: {sorted(set(int(e) for e in excess if e != 0))}")

        # PRACTICAL: can we predict H from just alpha_1?
        from numpy.linalg import lstsq
        A = np.column_stack([a1_arr, np.ones(N)])
        beta, _, _, _ = lstsq(A, h_arr, rcond=None)
        h_pred_a1 = A @ beta
        r2_a1 = 1 - np.sum((h_arr - h_pred_a1)**2) / np.sum((h_arr - np.mean(h_arr))**2)
        print(f"    H ≈ {beta[0]:.3f}*a1 + {beta[1]:.3f}: R² = {r2_a1:.4f}")

# ==============================================================
# 5. THE H=21 MECHANISM IN DETAIL
# ==============================================================
def h21_mechanism():
    """Why is H=21 forbidden? Analyze the (alpha_1, alpha_2) landscape."""
    print(f"\n{'='*60}")
    print(f"THE H=21 MECHANISM")
    print(f"{'='*60}")

    # H=21 = 1 + 2*a1 + 4*a2 + 8*a3 + ...
    # So 2*a1 + 4*a2 + 8*a3 + ... = 20
    # a1 + 2*a2 + 4*a3 + ... = 10

    # The "level-2 decompositions":
    # (a1, a2, a3) with a1 + 2*a2 + 4*a3 = 10 and a3 = 0 (simplest):
    # (10,0), (8,1), (6,2), (4,3), (2,4), (0,5)

    print(f"  H=21 requires a1 + 2*a2 + 4*a3 + ... = 10")
    print(f"  Level-2 decompositions (a3=0): (10,0), (8,1), (6,2), (4,3), (2,4), (0,5)")

    # For each n, check which decompositions are achievable
    for n in [5, 6, 7]:
        tiles = make_tiles(n)
        m = len(tiles)
        N = 1 << m

        print(f"\n  n={n}:")

        # Compute (a1, a2) for all tournaments and check if any give H=21
        pair_counts = Counter()
        h21_decomps = set()

        sample_size = min(N, 5000)
        indices = range(N) if N <= 5000 else np.random.RandomState(42).choice(N, 5000, replace=False)

        for t in indices:
            T = tiling_to_tournament(n, tiles, t)
            h = count_hamiltonian_paths(T, n)
            a1, a2, _ = find_odd_cycles_and_independence(T, n)
            pair_counts[(a1, a2)] += 1

            if int(h) == 21:
                h21_decomps.add((a1, a2))

        print(f"    Sampled {sample_size} tournaments")
        print(f"    H=21 found: {len(h21_decomps) > 0}")

        # For the target decompositions of H=21, are they achieved?
        for a1_target in [10, 8, 6, 4, 2, 0]:
            a2_target = (10 - a1_target) // 2
            if (10 - a1_target) % 2 != 0:
                continue
            achieved = any(a1 == a1_target and a2 == a2_target
                          for (a1, a2) in pair_counts)
            count = pair_counts.get((a1_target, a2_target), 0)
            print(f"    (a1={a1_target}, a2={a2_target}): "
                  f"{'ACHIEVED' if achieved else 'NOT FOUND'} ({count} times)")

            # If achieved, what H values do these tournaments actually have?
            if achieved:
                actual_H = set()
                for t in indices:
                    T = tiling_to_tournament(n, tiles, t)
                    a1, a2, _ = find_odd_cycles_and_independence(T, n)
                    if a1 == a1_target and a2 == a2_target:
                        h = count_hamiltonian_paths(T, n)
                        actual_H.add(int(h))
                print(f"      Actual H values at this (a1,a2): {sorted(actual_H)}")

# ==============================================================
# 6. PRACTICAL: H-VALUE DENSITY AND GAP PREDICTIONS
# ==============================================================
def h_density_and_predictions():
    """Compute the density of achievable H values and predict future gaps."""
    print(f"\n{'='*60}")
    print(f"H-VALUE DENSITY AND GAP PREDICTIONS")
    print(f"{'='*60}")

    achieved_by_n = {}
    for n in range(3, 8):
        tiles = make_tiles(n)
        m = len(tiles)
        H, _, _ = compute_all_H(n)
        achieved_by_n[n] = sorted(set(int(h) for h in H))

    for n in sorted(achieved_by_n.keys()):
        vals = achieved_by_n[n]
        max_h = max(vals)
        total_odd = (max_h + 1) // 2
        density = len(vals) / total_odd
        gaps_odd = [k for k in range(1, max_h+1, 2) if k not in set(vals)]
        print(f"\n  n={n}: {len(vals)} achievable H values, max={max_h}, "
              f"density={density:.3f}")
        print(f"    Odd gaps: {gaps_odd[:20]}{'...' if len(gaps_odd)>20 else ''}")

    # Gap evolution: which gaps CLOSE at each n?
    print(f"\n  Gap evolution:")
    all_gaps = {}
    for n in range(3, 8):
        s = set(achieved_by_n[n])
        max_h = max(achieved_by_n[n])
        all_gaps[n] = set(k for k in range(1, max_h+1, 2) if k not in s)

    for n in range(4, 8):
        closed = all_gaps[n-1] - all_gaps[n]
        opened = all_gaps[n] - all_gaps[n-1]
        persistent = all_gaps[n-1] & all_gaps[n]
        if closed:
            print(f"    n={n-1}→{n}: CLOSED gaps: {sorted(closed)[:15]}")
        if opened:
            print(f"    n={n-1}→{n}: NEW gaps: {sorted(opened)[:15]}")
        print(f"    n={n}: persistent gaps ≤50: "
              f"{sorted(g for g in persistent if g <= 50)}")

    # PREDICTION: gaps that persist through ALL tested n
    persistent_all = all_gaps[3]
    for n in range(4, 8):
        persistent_all &= all_gaps[n]
    # But only look at values achievable at some larger n
    permanent_candidates = sorted(g for g in persistent_all if g <= 100)
    print(f"\n  Gaps persistent through n=3..7 (≤100): {permanent_candidates}")
    print(f"  Of these, PROVED permanent: {[g for g in permanent_candidates if g in [7, 21]]}")
    print(f"  Unproved but persistent: {[g for g in permanent_candidates if g not in [7, 21]]}")

if __name__ == "__main__":
    alpha1_gap_landscape()
    alpha12_pairs()
    h_gap_census()
    alpha_distributions()
    h21_mechanism()
    h_density_and_predictions()
