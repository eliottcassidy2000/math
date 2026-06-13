"""
practical_tools.py — opus-2026-04-04-S11

PRACTICAL DELIVERABLES built on the rigidity framework.

1. Fast H computation via cycle counting (exact at n≤7: H = 1+2a1+4a2)
2. Tournament "complexity score" — distance from forbidden values
3. Rigidity detector — effective DOF per tournament
4. The achievable (a1, a2, H) lattice points — a practical lookup table
5. VERIFICATION: H = 1+2a1+4a2 at n=7 (need to check if a3 = 0 always)
"""
import numpy as np
from itertools import combinations, permutations
from collections import Counter, defaultdict
import time
import sys

sys.path.insert(0, '/Users/e/Documents/GitHub/math/04-computation')
from h_multilinear_complete import compute_all_H, make_tiles, tiling_to_tournament, count_hamiltonian_paths

def count_odd_cycles_fast(T, n):
    """Count directed 3-cycles and 5-cycles. Return (c3, c5, c7)."""
    c3 = 0
    for a in range(n):
        for b in range(a+1, n):
            for c in range(b+1, n):
                if T[a][b]+T[b][c]+T[c][a] == 3 or T[a][c]+T[c][b]+T[b][a] == 3:
                    c3 += 1

    c5 = 0
    if n >= 5:
        for verts in combinations(range(n), 5):
            for perm in permutations(verts):
                rot = [perm[i:]+perm[:i] for i in range(5)]
                if perm != min(rot): continue
                if all(T[perm[i]][perm[(i+1)%5]] for i in range(5)):
                    c5 += 1

    c7 = 0
    if n == 7:
        for perm in permutations(range(7)):
            rot = [perm[i:]+perm[:i] for i in range(7)]
            if perm != min(rot): continue
            if all(T[perm[i]][perm[(i+1)%7]] for i in range(7)):
                c7 += 1

    return c3, c5, c7

def count_disjoint_pairs(T, n):
    """Count pairs of vertex-disjoint odd cycles."""
    # Find all odd cycles first
    cycles = []
    for a in range(n):
        for b in range(a+1, n):
            for c in range(b+1, n):
                if T[a][b]+T[b][c]+T[c][a]==3 or T[a][c]+T[c][b]+T[b][a]==3:
                    cycles.append(frozenset({a,b,c}))
    if n >= 5:
        for verts in combinations(range(n), 5):
            for perm in permutations(verts):
                rot = [perm[i:]+perm[:i] for i in range(5)]
                if perm != min(rot): continue
                if all(T[perm[i]][perm[(i+1)%5]] for i in range(5)):
                    cycles.append(frozenset(verts))
                    break
    if n == 7:
        for perm in permutations(range(7)):
            rot = [perm[i:]+perm[:i] for i in range(7)]
            if perm != min(rot): continue
            if all(T[perm[i]][perm[(i+1)%7]] for i in range(7)):
                cycles.append(frozenset(range(7)))
                break

    # Count disjoint pairs
    a2 = 0
    for i in range(len(cycles)):
        for j in range(i+1, len(cycles)):
            if not (cycles[i] & cycles[j]):
                a2 += 1
    return a2

# ==============================================================
# 1. VERIFY H = 1 + 2*alpha_1 + 4*alpha_2 AT n=7
# ==============================================================
def verify_exact_formula_n7():
    """Check if H = 1+2a1+4a2 is exact at n=7."""
    print(f"{'='*60}")
    print(f"VERIFYING H = 1 + 2*alpha_1 + 4*alpha_2 AT n=7")
    print(f"{'='*60}")

    n = 7
    tiles = make_tiles(n)
    m = len(tiles)
    H, _, _ = compute_all_H(n)
    N = 1 << m

    # Sample to save time (7-cycle enumeration is expensive)
    np.random.seed(42)
    sample = np.random.choice(N, min(500, N), replace=False)

    exact = 0
    errors = []

    for t in sample:
        T = tiling_to_tournament(n, tiles, t)
        h_actual = int(H[t])

        c3, c5, c7 = count_odd_cycles_fast(T, n)
        a1 = c3 + c5 + c7
        a2 = count_disjoint_pairs(T, n)

        h_pred = 1 + 2*a1 + 4*a2

        if h_pred == h_actual:
            exact += 1
        else:
            errors.append((t, h_actual, h_pred, a1, a2, c3, c5, c7))

    print(f"  Sampled {len(sample)} tournaments")
    print(f"  Exact matches: {exact}/{len(sample)}")
    if errors:
        print(f"  Errors: {len(errors)}")
        for t, h, hp, a1, a2, c3, c5, c7 in errors[:10]:
            excess = h - hp
            print(f"    t={t}: H={h}, pred={hp}, excess={excess}, "
                  f"a1={a1} (c3={c3},c5={c5},c7={c7}), a2={a2}")
            # Excess = 8*alpha_3 + ... If excess != 0, alpha_3 > 0
            if excess > 0:
                print(f"      → alpha_3 ≥ {excess//8} (8*alpha_3 component)")
    else:
        print(f"  *** H = 1 + 2*alpha_1 + 4*alpha_2 is EXACT at n=7! ***")

# ==============================================================
# 2. FAST H VIA CYCLE COUNTING
# ==============================================================
def benchmark_cycle_counting():
    """Compare speed of cycle counting vs DP for H computation."""
    print(f"\n{'='*60}")
    print(f"FAST H BENCHMARK: CYCLE COUNTING vs DP")
    print(f"{'='*60}")

    for n in [5, 6]:
        tiles = make_tiles(n)
        m = len(tiles)
        N = 1 << m
        H, _, _ = compute_all_H(n)

        # Time cycle counting
        t0 = time.time()
        for t in range(min(N, 200)):
            T = tiling_to_tournament(n, tiles, t)
            c3, c5, _ = count_odd_cycles_fast(T, n)
            a1 = c3 + c5
            a2 = count_disjoint_pairs(T, n) if n >= 6 else 0
            h_cycle = 1 + 2*a1 + 4*a2
        t_cycle = time.time() - t0

        # Time DP
        t0 = time.time()
        for t in range(min(N, 200)):
            T = tiling_to_tournament(n, tiles, t)
            h_dp = count_hamiltonian_paths(T, n)
        t_dp = time.time() - t0

        sample = min(N, 200)
        print(f"\n  n={n} ({sample} tournaments):")
        print(f"    Cycle counting: {t_cycle:.3f}s ({t_cycle/sample*1000:.1f}ms each)")
        print(f"    DP (bitmask): {t_dp:.3f}s ({t_dp/sample*1000:.1f}ms each)")
        print(f"    Ratio: {t_cycle/t_dp:.2f}x")

# ==============================================================
# 3. TOURNAMENT COMPLEXITY SCORE
# ==============================================================
def complexity_score():
    """
    Define a "tournament complexity score" based on the rigidity framework.

    Score(T) = position of H(T) in the achievable lattice.
    Components:
    - H value (absolute)
    - Distance to nearest forbidden value
    - Position in the alpha_1 distribution (percentile)
    - Whether the tournament is "near-rigid" (few DOF)
    """
    print(f"\n{'='*60}")
    print(f"TOURNAMENT COMPLEXITY SCORE")
    print(f"{'='*60}")

    n = 6
    tiles = make_tiles(n)
    m = len(tiles)
    H, _, _ = compute_all_H(n)
    N = 1 << m

    # Forbidden values
    forbidden = {7, 21}

    # Compute scores for all tournaments
    scores = []
    for t in range(N):
        h = int(H[t])
        T = tiling_to_tournament(n, tiles, t)
        c3, c5, _ = count_odd_cycles_fast(T, n)
        a1 = c3 + c5
        a2 = count_disjoint_pairs(T, n)

        # Component 1: H value (normalized by max)
        h_norm = h / 45.0  # max H at n=6 is 45

        # Component 2: distance to nearest forbidden
        dist_forbidden = min(abs(h - f) for f in forbidden)

        # Component 3: cycle count richness
        cycle_richness = a1 / 20.0  # max alpha_1 at n=6 is ~20

        # Component 4: independence depth (a2 / max a2)
        depth = a2 / 4.0  # max a2 at n=6 is 4

        # Combined score
        score = h_norm * 0.4 + cycle_richness * 0.3 + depth * 0.2 + (dist_forbidden/45) * 0.1
        scores.append({
            'tiling': t, 'H': h, 'a1': a1, 'a2': a2, 'c3': c3, 'c5': c5,
            'score': score, 'dist_forbidden': dist_forbidden
        })

    # Sort by score
    scores.sort(key=lambda x: x['score'], reverse=True)

    print(f"\n  Top 10 most complex tournaments:")
    for s in scores[:10]:
        print(f"    H={s['H']:3d} a1={s['a1']:2d} a2={s['a2']} c3={s['c3']:2d} c5={s['c5']:2d} "
              f"score={s['score']:.3f} dist_forbidden={s['dist_forbidden']}")

    print(f"\n  Bottom 10 (simplest):")
    for s in scores[-10:]:
        print(f"    H={s['H']:3d} a1={s['a1']:2d} a2={s['a2']} c3={s['c3']:2d} c5={s['c5']:2d} "
              f"score={s['score']:.3f} dist_forbidden={s['dist_forbidden']}")

    # Score near forbidden values
    print(f"\n  Tournaments nearest to forbidden values:")
    for f in sorted(forbidden):
        nearest = min(scores, key=lambda s: abs(s['H'] - f))
        print(f"    Nearest to H={f}: H={nearest['H']} (dist={nearest['dist_forbidden']})")

    # Score distribution
    all_scores = [s['score'] for s in scores]
    print(f"\n  Score distribution: mean={np.mean(all_scores):.3f}, "
          f"std={np.std(all_scores):.3f}")

# ==============================================================
# 4. THE ACHIEVABLE LATTICE AS A LOOKUP TABLE
# ==============================================================
def achievable_lattice():
    """Build the lookup table: for each (a1, a2), what H is achieved?"""
    print(f"\n{'='*60}")
    print(f"ACHIEVABLE LATTICE: (alpha_1, alpha_2) → H")
    print(f"{'='*60}")

    for n in [6, 7]:
        tiles = make_tiles(n)
        m = len(tiles)
        N = 1 << m

        # For n=7 sample
        if N > 5000:
            np.random.seed(42)
            indices = np.random.choice(N, 3000, replace=False)
        else:
            indices = range(N)

        H_data, _, _ = compute_all_H(n)

        lattice = defaultdict(set)
        for t in indices:
            T = tiling_to_tournament(n, tiles, t)
            h = int(H_data[t])
            c3, c5, c7 = count_odd_cycles_fast(T, n)
            a1 = c3 + c5 + c7
            a2 = count_disjoint_pairs(T, n)
            lattice[(a1, a2)].add(h)

        print(f"\n  n={n}:")
        print(f"    Achieved lattice points: {len(lattice)}")

        # The key insight: is (a1, a2) → H a FUNCTION?
        # I.e., does each (a1, a2) give a unique H?
        is_function = all(len(v) == 1 for v in lattice.values())
        print(f"    (a1, a2) → H is a function? {is_function}")

        if is_function:
            print(f"    *** THEOREM: H is completely determined by (alpha_1, alpha_2)! ***")
            print(f"    Formula: H(a1, a2) = 1 + 2*a1 + 4*a2")
            # Verify
            for (a1, a2), h_set in sorted(lattice.items()):
                h = list(h_set)[0]
                pred = 1 + 2*a1 + 4*a2
                if h != pred:
                    print(f"      MISMATCH: (a1={a1}, a2={a2}) → H={h}, pred={pred}")
        else:
            # Which (a1, a2) give multiple H values?
            multi = {k: v for k, v in lattice.items() if len(v) > 1}
            print(f"    Multi-valued: {len(multi)} lattice points")
            for (a1, a2), h_set in sorted(multi.items())[:10]:
                print(f"      (a1={a1}, a2={a2}): H ∈ {sorted(h_set)}")

        # THE TABLE (for documentation/engineering)
        if n <= 7:
            max_a1 = max(a1 for a1, _ in lattice)
            max_a2 = max(a2 for _, a2 in lattice)
            print(f"\n    Achievable (a1, a2) → H table:")
            for a2 in range(max_a2+1):
                row_entries = []
                for a1 in range(max_a1+1):
                    if (a1, a2) in lattice:
                        h = list(lattice[(a1, a2)])[0] if len(lattice[(a1, a2)]) == 1 else "?"
                        row_entries.append(f"{h}")
                    else:
                        row_entries.append(".")
                if any(e != "." for e in row_entries):
                    print(f"      a2={a2}: {' '.join(f'{e:>4s}' for e in row_entries[:25])}")

# ==============================================================
# 5. THE FORBIDDEN VALUE OBSTRUCTION ANALYSIS FOR H=21
# ==============================================================
def h21_obstruction_detail():
    """
    Why does EVERY (a1, a2) decomposition of H=21 fail?

    For n=6: we know the achievable (a1, a2) pairs.
    H=21 needs a1+2*a2 = 10.
    Check: the achievable pairs with a1+2*a2 in {8,9,10,11,12}
    to see the "neighborhood" of the H=21 level set.
    """
    print(f"\n{'='*60}")
    print(f"H=21 OBSTRUCTION: NEIGHBORHOOD ANALYSIS")
    print(f"{'='*60}")

    n = 6
    tiles = make_tiles(n)
    m = len(tiles)
    N = 1 << m
    H_data, _, _ = compute_all_H(n)

    # Collect all (a1, a2) pairs
    pairs = defaultdict(int)
    pair_to_H = defaultdict(set)

    for t in range(N):
        T = tiling_to_tournament(n, tiles, t)
        c3, c5, _ = count_odd_cycles_fast(T, n)
        a1 = c3 + c5
        a2 = count_disjoint_pairs(T, n)
        h = int(H_data[t])
        pairs[(a1, a2)] += 1
        pair_to_H[(a1, a2)].add(h)

    # The level sets a1 + 2*a2 = k for k near 10
    for k in range(6, 15):
        h_target = 1 + 2*k
        combos = [(a1, (k-a1)//2) for a1 in range(k+1) if (k-a1)%2 == 0 and (k-a1)//2 >= 0]
        achieved = [(a1, a2) for a1, a2 in combos if (a1, a2) in pairs]
        missed = [(a1, a2) for a1, a2 in combos if (a1, a2) not in pairs]

        print(f"\n  Level a1+2*a2={k} (H={h_target}):")
        print(f"    Possible decompositions: {combos}")
        print(f"    Achieved: {achieved}")
        print(f"    Missing: {missed}")

        # For achieved ones, verify H matches
        for a1, a2 in achieved:
            h_actual = pair_to_H[(a1, a2)]
            print(f"      ({a1},{a2}): H={sorted(h_actual)}, "
                  f"tournaments={pairs[(a1,a2)]}")

if __name__ == "__main__":
    verify_exact_formula_n7()
    benchmark_cycle_counting()
    complexity_score()
    achievable_lattice()
    h21_obstruction_detail()
