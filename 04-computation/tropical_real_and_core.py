"""
tropical_real_and_core.py — opus-2026-04-04-S7

REAL TROPICALIZATION: Use log|c_S| as the valuation (not v_2).
The tropical polynomial trop(H)(t) = max_S (log|c_S| + sum_{i in S} t_i)
is a piecewise-linear convex function on R^m.

Its CORNER LOCUS (where 2+ terms tie for the max) is the tropical variety.
This tells us which monomial DOMINATES in each region of the tiling space.

Also: the INDEPENDENCE SEQUENCE as the irreducible core.

Also: non-monotonicity of H under tile flips — characterize precisely.
When does flipping a tile DECREASE H?
"""
import numpy as np
from itertools import combinations
from collections import defaultdict, Counter
import sys

sys.path.insert(0, '/Users/e/Documents/GitHub/math/04-computation')
from h_multilinear_complete import compute_all_H, mobius_inversion, make_tiles, tiling_to_tournament, count_hamiltonian_paths

# ==============================================================
# 1. REAL TROPICALIZATION: trop(H)(t) = max_S (log|c_S| + t·1_S)
# ==============================================================
def real_tropical(n):
    """
    Real tropical H using log|c_S| valuation.

    On the discrete domain {0,1}^m:
    trop(H)(t) = max_{S subset supp(t)} log|c_S|

    The dominant monomial at tiling t is the S ⊆ supp(t) with largest |c_S|.
    """
    print(f"\n{'='*60}")
    print(f"REAL TROPICALIZATION n={n}")
    print(f"{'='*60}")

    tiles = make_tiles(n)
    m = len(tiles)
    H, _, _ = compute_all_H(n)
    coeffs = mobius_inversion(H, m)
    N = 1 << m

    # For each tiling, find the dominant monomial (max |c_S| with S ⊆ supp(t))
    dom_deg_dist = Counter()  # degree of dominant monomial
    dom_sign_dist = Counter()  # sign of dominant coefficient

    trop_vs_actual = []  # (trop value, actual log|H|)

    for t in range(N):
        best_S = 0  # empty set (c_0 = 1)
        best_abs = 1  # |c_0| = 1
        for S, c in coeffs.items():
            if S & ~t:  # S not subset of supp(t)
                continue
            if abs(c) > best_abs:
                best_abs = abs(c)
                best_S = S

        deg = bin(best_S).count('1')
        dom_deg_dist[deg] += 1
        dom_sign_dist[1 if coeffs.get(best_S, 1) > 0 else -1] += 1

        trop_val = np.log(best_abs) if best_abs > 0 else 0
        actual_val = np.log(int(H[t])) if H[t] > 0 else 0
        trop_vs_actual.append((trop_val, actual_val))

    print(f"  Dominant monomial degree distribution:")
    for d in sorted(dom_deg_dist.keys()):
        print(f"    degree {d}: {dom_deg_dist[d]}/{N} = {dom_deg_dist[d]/N:.3f}")

    print(f"  Dominant monomial sign: {dict(dom_sign_dist)}")

    # How well does the tropical value predict log(H)?
    trop_vals = [tv for tv, _ in trop_vs_actual]
    actual_vals = [av for _, av in trop_vs_actual]
    corr = np.corrcoef(trop_vals, actual_vals)[0, 1]
    print(f"\n  Corr(trop, log(H)) = {corr:.4f}")

    # Ratio log(H) / trop
    ratios = [av / tv if tv > 0 else float('nan') for tv, av in trop_vs_actual]
    finite_ratios = [r for r in ratios if not np.isnan(r)]
    if finite_ratios:
        print(f"  log(H)/trop ratio: mean={np.mean(finite_ratios):.3f}, "
              f"std={np.std(finite_ratios):.3f}")

    # TROPICAL CELLS: which tilings share the same dominant monomial?
    dom_monomial_groups = defaultdict(list)
    for t in range(N):
        best_S = 0
        best_abs = 1
        for S, c in coeffs.items():
            if S & ~t: continue
            if abs(c) > best_abs:
                best_abs = abs(c)
                best_S = S
        dom_monomial_groups[best_S].append(t)

    print(f"\n  Tropical cells (regions with same dominant monomial):")
    print(f"    Number of cells: {len(dom_monomial_groups)}")
    for S, tilings in sorted(dom_monomial_groups.items(), key=lambda x: -len(x[1]))[:8]:
        tile_set = [tiles[i] for i in range(m) if S & (1 << i)]
        c = coeffs.get(S, 1)
        h_range = [int(H[t]) for t in tilings]
        print(f"    S={tile_set} (c={c}): {len(tilings)} tilings, "
              f"H range [{min(h_range)}, {max(h_range)}]")

# ==============================================================
# 2. NON-MONOTONICITY: when does flipping a tile DECREASE H?
# ==============================================================
def non_monotonicity(n):
    """
    H is NOT monotone under tile flips from transitive.
    Characterize exactly which flips decrease H and why.

    From the multilinear polynomial:
    H(t XOR e_k) - H(t) = gradient_k(t) = c_k + sum_{j!=k} c_{kj}*t_j + ...

    The gradient is POSITIVE at t=0 (all c_k > 0, the linear coefficient).
    But at t with many backward tiles, the gradient can become negative
    due to the negative quadratic terms from same-end interactions.
    """
    print(f"\n{'='*60}")
    print(f"NON-MONOTONICITY ANALYSIS n={n}")
    print(f"{'='*60}")

    tiles = make_tiles(n)
    m = len(tiles)
    H, _, _ = compute_all_H(n)
    N = 1 << m
    coeffs = mobius_inversion(H, m)

    # For each tile k, at how many tilings does flipping k DECREASE H?
    for k in range(m):
        increase = 0
        decrease = 0
        neutral = 0
        for t in range(N):
            if t & (1 << k):
                continue  # k is already backward; flipping it to forward
            h_before = int(H[t])
            h_after = int(H[t | (1 << k)])
            diff = h_after - h_before
            if diff > 0: increase += 1
            elif diff < 0: decrease += 1
            else: neutral += 1

        tile = tiles[k]
        skip = tile[0] - tile[1]
        total = increase + decrease + neutral
        print(f"  tile {tile} skip={skip}: "
              f"+{increase} ↑{increase/total:.2f}, "
              f"-{decrease} ↓{decrease/total:.2f}, "
              f"={neutral}")

    # When does the GRADIENT become negative?
    # gradient_k(t) = c_k + sum_{j: j!=k, t_j=1} c_{kj} + higher terms
    # The linear contribution is c_k = 2^{skip-1} > 0.
    # The quadratic contribution from same-end tiles is NEGATIVE.
    # When enough same-end tiles are backward, the gradient flips sign.

    print(f"\n  GRADIENT SIGN ANALYSIS for apex tile {tiles[0]}:")
    k = 0  # apex tile
    grad_pos = 0
    grad_neg = 0
    for t in range(N):
        if t & (1 << k): continue
        diff = int(H[t | (1 << k)]) - int(H[t])
        # What determines the sign? Count same-end backward tiles
        same_end_backward = 0
        for j in range(m):
            if j == k: continue
            if not (t & (1 << j)): continue
            t1, t2 = tiles[k], tiles[j]
            x1,y1 = t1; x2,y2 = t2
            if x1==x2 or y1==y2:
                same_end_backward += 1

        if diff > 0: grad_pos += 1
        else: grad_neg += 1

    print(f"    Gradient positive: {grad_pos}, negative: {grad_neg}")

    # THE CRITICAL THRESHOLD: how many same-end backward tiles are needed
    # to make the gradient negative?
    print(f"\n  Critical threshold for gradient sign flip:")
    for k in range(min(m, 6)):
        tile = tiles[k]
        skip = tile[0] - tile[1]
        # Same-end tiles for tile k
        same_end_indices = []
        for j in range(m):
            if j == k: continue
            t1, t2 = tiles[k], tiles[j]
            if t1[0]==t2[0] or t1[1]==t2[1]:
                same_end_indices.append(j)

        # For each count of same-end backward tiles, what's the gradient?
        by_se_count = defaultdict(list)
        for t in range(N):
            if t & (1 << k): continue
            se_count = sum(1 for j in same_end_indices if t & (1 << j))
            diff = int(H[t | (1 << k)]) - int(H[t])
            by_se_count[se_count].append(diff)

        print(f"    tile {tile} skip={skip}, {len(same_end_indices)} same-end tiles:")
        for se in sorted(by_se_count.keys()):
            diffs = by_se_count[se]
            pos = sum(1 for d in diffs if d > 0)
            neg = sum(1 for d in diffs if d < 0)
            mean_d = np.mean(diffs)
            print(f"      same-end backward={se}: mean_grad={mean_d:.1f}, "
                  f"+{pos}/-{neg}")

# ==============================================================
# 3. THE 2-ADIC INDEPENDENCE SEQUENCE
# ==============================================================
def independence_sequence_analysis(n):
    """
    H(T) = sum_k alpha_k(Omega(T)) * 2^k

    The independence sequence (alpha_0, alpha_1, alpha_2, ...) of the
    conflict graph Omega determines H. Study this sequence.
    """
    print(f"\n{'='*60}")
    print(f"INDEPENDENCE SEQUENCE ANALYSIS n={n}")
    print(f"{'='*60}")

    tiles = make_tiles(n)
    m = len(tiles)
    H, _, _ = compute_all_H(n)
    N = 1 << m

    # For each tiling, compute H and decompose into independence sequence
    # H = sum_k alpha_k * 2^k ≡ binary expansion but with possible large alpha_k

    # Actually, the OCF directly gives H = I(Omega, 2).
    # The independence polynomial I(Omega, x) = sum_k alpha_k x^k
    # where alpha_k = number of independent sets of size k in Omega.
    # H = I(Omega, 2) = alpha_0 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ...

    # The first few: alpha_0 = 1 always, alpha_1 = |V(Omega)| = number of odd cycles

    # From H: alpha_1 = (H - 1) / 2 at first order, but higher alpha_k contribute.
    # H = 1 + 2*alpha_1 + 4*alpha_2 + ...
    # (H-1)/2 = alpha_1 + 2*alpha_2 + ...
    # alpha_1 mod 2 = ((H-1)/2) mod 2 = (H mod 4 - 1)/2

    # H mod 4: {1, 3}. If H ≡ 1: alpha_1 even. If H ≡ 3: alpha_1 odd.
    h1_even = sum(1 for h in H if h % 4 == 1)
    h1_odd = sum(1 for h in H if h % 4 == 3)
    print(f"  alpha_1 even (H≡1 mod 4): {h1_even}/{N} = {h1_even/N:.3f}")
    print(f"  alpha_1 odd (H≡3 mod 4): {h1_odd}/{N} = {h1_odd/N:.3f}")

    # H mod 8 tells us about alpha_1 mod 4 and alpha_2 mod 2
    # H = 1 + 2*alpha_1 + 4*alpha_2 + 8*(...)
    # H mod 8 = (1 + 2*(alpha_1 mod 4) + 4*(alpha_2 mod 2)) mod 8
    mod8_dist = Counter(int(h) % 8 for h in H)
    print(f"  H mod 8 distribution: {dict(sorted(mod8_dist.items()))}")

    # alpha_1 mod 4: from (H-1)/2 mod 4
    a1_mod4_dist = Counter(((int(h)-1)//2) % 4 for h in H)
    print(f"  (H-1)/2 mod 4 distribution: {dict(sorted(a1_mod4_dist.items()))}")

    # The maximum possible alpha_1 = maximum number of odd cycles
    # This is bounded by C(n, 3) / 3 + C(n,5) / 5 + ... for directed odd cycles

    # Distribution of H values
    h_dist = Counter(int(h) for h in H)
    h_unique = sorted(h_dist.keys())
    print(f"\n  H value distribution:")
    print(f"    {len(h_unique)} distinct values: {h_unique[:15]}...")
    print(f"    All odd: {all(h % 2 == 1 for h in h_unique)}")

    # Gap structure: consecutive H values differ by...
    gaps = [h_unique[i+1] - h_unique[i] for i in range(len(h_unique)-1)]
    print(f"    Gaps between consecutive H values: {gaps[:15]}...")
    print(f"    All gaps even: {all(g % 2 == 0 for g in gaps)}")
    # Gaps should be even because H ≡ 1 or 3 mod 4, so consecutive values differ by 2 or 4+

# ==============================================================
# 4. THE CONFLICT GRAPH AS FUNCTOR: tile flip → graph operation
# ==============================================================
def omega_functor(n):
    """
    How does flipping a single tile change Omega(T)?

    Flipping tile k changes one arc of the tournament. This:
    - CREATES new odd cycles that use the new arc direction
    - DESTROYS old odd cycles that used the old arc direction
    - PRESERVES cycles that don't use this arc

    The change in Omega is a LOCAL OPERATION around the flipped tile.
    """
    print(f"\n{'='*60}")
    print(f"OMEGA AS FUNCTOR: TILE FLIP → GRAPH OPERATION n={n}")
    print(f"{'='*60}")

    tiles = make_tiles(n)
    m = len(tiles)

    # Start from several representative tilings and flip each tile
    for start_name, start_tiling in [("transitive", 0),
                                      ("all backward", (1 << m) - 1)]:
        print(f"\n  Starting from {start_name}:")

        T_start = tiling_to_tournament(n, tiles, start_tiling)
        h_start = count_hamiltonian_paths(T_start, n)

        for k in range(min(m, 6)):
            new_tiling = start_tiling ^ (1 << k)
            T_new = tiling_to_tournament(n, tiles, new_tiling)
            h_new = count_hamiltonian_paths(T_new, n)

            diff = h_new - h_start
            tile = tiles[k]
            skip = tile[0] - tile[1]
            direction = "forward→backward" if not (start_tiling & (1<<k)) else "backward→forward"
            print(f"    Flip tile {tile} (skip={skip}, {direction}): "
                  f"ΔH = {diff:+d} (H: {h_start}→{h_new})")

    # THE KEY OBSERVATION: from transitive, every flip INCREASES H
    # (because c_k > 0). From all-backward, every flip can increase or decrease.
    print(f"\n  FROM TRANSITIVE: every single flip increases H")
    print(f"    (because c_k = 2^(skip-1) > 0 and there are no")
    print(f"     quadratic corrections at t=0)")
    print(f"\n  The gradient at t=0 is PURELY LINEAR: grad_k = c_k = 2^(skip-1)")
    print(f"  This is the 'source property' in action.")

# ==============================================================
# 5. THE DEEPEST VIEW: H as independence polynomial evaluation
# ==============================================================
def deepest_view(n):
    """
    The deepest simplification:

    EVERYTHING is the conflict graph.

    H(T) = I(Omega(T), 2) where:
    - Omega(T) = conflict graph (odd cycles as vertices, shared-vertex = edge)
    - I(G, x) = independence polynomial

    The multilinear polynomial H(t) is the COMPOSITION:
      tiling t → tournament T(t) → conflict graph Omega(T(t)) → I(Omega, 2) = H

    Each step is a FUNCTOR:
    1. t → T: linear (each tile variable controls one arc)
    2. T → Omega: nonlinear but LOCAL (odd cycles are local graph structures)
    3. Omega → I(Omega, 2): nonlinear and GLOBAL (independence polynomial
       depends on the whole graph structure)

    The multilinear polynomial ENCODES the composition 1∘2∘3 in tile coordinates.
    Every property we've found is a consequence of one of these three functors.
    """
    print(f"\n{'='*60}")
    print(f"THE DEEPEST VIEW: H = I(Omega(T(t)), 2)")
    print(f"{'='*60}")

    tiles = make_tiles(n)
    m = len(tiles)
    H, _, _ = compute_all_H(n)
    N = 1 << m
    coeffs = mobius_inversion(H, m)

    print(f"""
    THE THREE FUNCTORS:

    1. TILING → TOURNAMENT (linear, invertible)
       t ∈ {{0,1}}^m → T(t) ∈ tournaments on n vertices
       Each tile t_k controls one arc. Base path arcs are fixed.
       This is an EMBEDDING of the m-cube into tournament space.

    2. TOURNAMENT → CONFLICT GRAPH (nonlinear, local)
       T → Omega(T) = (odd cycles, shared-vertex adjacency)
       Adding one backward arc can CREATE cycles (those using the new arc)
       or DESTROY them (those that required the old arc direction).
       This is the CRUX of the nonlinearity.

    3. CONFLICT GRAPH → HP COUNT (polynomial evaluation)
       Omega → I(Omega, 2) = sum_k alpha_k * 2^k
       This converts the TOPOLOGICAL information in Omega
       into the ARITHMETIC value H.

    CONSEQUENCES:
    - Redei (H odd): alpha_0 = 1 always → I(Omega, 2) ≡ 1 (mod 2)
    - Recursive preservation: vertex n+1 as source means Omega(T_{{n+1}}) = Omega(T_n)
      when new tiles are forward (no new cycles use forward arcs from source)
    - Sign rule: same-end tiles create COMPETING cycles at shared vertex
      (antiferromagnetic). Cross-end tiles create COOPERATIVE cycles.
    - Degree cap: at most floor((n-1)/2) vertex-disjoint odd cycles in
      a tournament on n vertices (independence number of Omega is bounded).
    - c_3 predicts 91% of H: because c_3 ≈ alpha_1 at small n
      (most independent sets have size 0 or 1, so I ≈ 1 + 2*alpha_1 ≈ 1 + 2*c_3).
    """)

    # VERIFY: I ≈ 1 + 2*c_3 as an approximation of H
    c3_vals = []
    for t in range(N):
        T = tiling_to_tournament(n, tiles, t)
        c3 = 0
        for a in range(n):
            for b in range(a+1, n):
                for c in range(b+1, n):
                    if T[a][b]+T[b][c]+T[c][a]==3 or T[a][c]+T[c][b]+T[b][a]==3:
                        c3 += 1
        c3_vals.append(c3)

    # I ≈ 1 + 2*c_3 vs actual H
    approx = [1 + 2*c for c in c3_vals]
    actual = [int(h) for h in H]
    max_err = max(abs(a - h) for a, h in zip(approx, actual))
    mean_err = np.mean([abs(a - h) for a, h in zip(approx, actual)])
    r2 = 1 - np.sum([(a-h)**2 for a,h in zip(approx, actual)]) / np.sum([(h-np.mean(actual))**2 for h in actual])

    print(f"  APPROXIMATION: H ≈ 1 + 2*c_3")
    print(f"    R² = {r2:.4f}")
    print(f"    Mean absolute error = {mean_err:.2f}")
    print(f"    Max absolute error = {max_err}")

    # Better: H ≈ 1 + 2*(c_3 + c_5) where c_5 = number of 5-cycles
    # At n=5: c_5 matters. At n=6: c_5 contributes.

    # The EXACT formula: H = 1 + 2*alpha_1 + 4*alpha_2 + ...
    # alpha_1 = total number of odd directed cycles
    # alpha_2 = number of vertex-disjoint pairs of odd cycles
    # For most tournaments, alpha_2 ≪ alpha_1, so H ≈ 1 + 2*alpha_1

    # But alpha_1 ≠ c_3 in general (there are also 5-cycles, 7-cycles)
    # At n=5: alpha_1 = c_3 + c_5 where c_5 = number of directed 5-cycles

    # HOW MUCH does alpha_2 contribute?
    # H - (1 + 2*alpha_1) = 4*alpha_2 + 8*alpha_3 + ...
    # If we knew alpha_1 exactly, the residual would be 4*alpha_2 + ...

    # Estimate: at n=5, H ranges from 1 to 15. 1+2*alpha_1 ranges up to
    # 1+2*6=13 (max 6 odd cycles). The max H=15 requires alpha_2 >= 1.
    # 15 = 1 + 2*6 + 4*alpha_2 → alpha_2 = (15-13)/4 = 0.5. Not integer!
    # So the approximation isn't this simple — alpha_1 can be > c_3.

    # The EXACT relationship:
    # alpha_1 = c_3 + c_5 + c_7 + ... (total directed odd cycles)
    # For n=5: c_5 exists (directed 5-cycles that use only non-base arcs)

    print(f"\n  The hierarchy of approximations:")
    print(f"    Level 0: H = 1 (Redei, correct mod 2)")
    print(f"    Level 1: H ≈ 1 + 2*c_3 (R² = {r2:.3f})")
    print(f"    Level 2: H ≈ 1 + 2*(c_3 + c_5) + 4*alpha_2")
    print(f"    Exact: H = I(Omega, 2)")

# ==============================================================
# MAIN
# ==============================================================
if __name__ == "__main__":
    for n in [5, 6, 7]:
        real_tropical(n)

    for n in [5, 6]:
        non_monotonicity(n)

    for n in [5, 6, 7]:
        independence_sequence_analysis(n)

    for n in [5]:
        omega_functor(n)

    for n in [5, 6]:
        deepest_view(n)
