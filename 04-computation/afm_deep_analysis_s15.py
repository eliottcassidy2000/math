#!/usr/bin/env python3
"""
afm_deep_analysis_s15.py — opus-2026-04-04-S15

DEEPER ANTIFERROMAGNETIC ANALYSIS

1. Magnon dispersion at n=6 (where isotropy breaks)
2. Boltzmann-weighted spin correlations: <σ_i σ_j>_β
3. The "true" AFM energy: E_frustration based on vertex triples
4. Susceptibility and specific heat from H-weighted ensemble
5. The antiferromagnetic GROUND STATE: which tournaments minimize E_frustration?
6. H=7 and H=21 gap analysis via AFM quantization
"""

import numpy as np
from itertools import combinations, permutations
from math import comb, factorial
from collections import defaultdict, Counter
import random
import sys
import time

random.seed(42)
np.random.seed(42)

# ─────────────────────────────────────────────
#  Core tournament utilities
# ─────────────────────────────────────────────

def all_tournaments(n):
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i, j) in enumerate(edges):
            if (mask >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        yield mask, A

def compute_H(A, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if not (mask & (1 << u)) and A[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[(1 << n) - 1])

def scores(A, n):
    return tuple(sorted(sum(A[i][j] for j in range(n) if j != i) for i in range(n)))

def count_3cycles(A, n):
    c = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j] and A[j][k] and A[k][i]:
                    c += 1
                elif A[i][k] and A[k][j] and A[j][i]:
                    c += 1
    return c

def staircase_tiles(n):
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            if x - y >= 2:
                tiles.append((x, y))
    return tiles

def tournament_to_spins(A, n, tiles):
    spins = []
    for (a, b) in tiles:
        if A[a-1][b-1]:
            spins.append(1)
        else:
            spins.append(-1)
    return spins

def tile_neighbors(tiles, n):
    adj = defaultdict(set)
    for i, (a1, b1) in enumerate(tiles):
        for j, (a2, b2) in enumerate(tiles):
            if i < j and (a1 == a2 or a1 == b2 or b1 == a2 or b1 == b2):
                adj[i].add(j)
                adj[j].add(i)
    return adj


# ─────────────────────────────────────────────
#  1. MAGNON DISPERSION AT n=6
# ─────────────────────────────────────────────

def magnon_dispersion_n6():
    """Compute magnon spectrum at n=6 where isotropy is expected to break."""
    n = 6
    tiles = staircase_tiles(n)
    m = len(tiles)
    print(f"\n{'='*80}")
    print(f"  MAGNON DISPERSION AT n=6 (10 tiles, 32768 tournaments)")
    print(f"{'='*80}")
    print(f"  Tiles: {tiles}")

    edges = [(i, j) for i in range(n) for j in range(i+1, n)]

    delta_H_by_tile = defaultdict(list)
    delta_H_abs_by_tile = defaultdict(list)

    count = 0
    for mask, A in all_tournaments(n):
        H0 = compute_H(A, n)
        for t_idx, (a, b) in enumerate(tiles):
            A2 = [row[:] for row in A]
            A2[a-1][b-1] = 1 - A2[a-1][b-1]
            A2[b-1][a-1] = 1 - A2[b-1][a-1]
            H1 = compute_H(A2, n)
            dH = H1 - H0
            delta_H_by_tile[t_idx].append(dH)
            delta_H_abs_by_tile[t_idx].append(abs(dH))

        count += 1
        if count % 5000 == 0:
            print(f"  ... {count}/{1 << comb(n,2)} done", file=sys.stderr)

    print(f"\n  Tile-by-tile magnon spectrum:")
    print(f"  {'Tile':>8}  {'skip':>4}  {'mean|ΔH|':>10}  {'std|ΔH|':>10}  "
          f"{'mean ΔH':>10}  {'ΔH=0 %':>8}  {'median|ΔH|':>12}")

    results = []
    for t_idx, (a, b) in enumerate(tiles):
        skip = a - b
        abs_vals = np.array(delta_H_abs_by_tile[t_idx])
        signed_vals = np.array(delta_H_by_tile[t_idx])
        mean_abs = np.mean(abs_vals)
        std_abs = np.std(abs_vals)
        mean_signed = np.mean(signed_vals)
        zero_pct = 100 * np.sum(abs_vals == 0) / len(abs_vals)
        median_abs = np.median(abs_vals)

        results.append({
            'tile': (a, b), 'skip': skip,
            'mean_abs': mean_abs, 'std_abs': std_abs,
            'mean_signed': mean_signed, 'zero_pct': zero_pct,
            'median_abs': median_abs
        })
        print(f"  ({a},{b})   {skip:4d}  {mean_abs:10.4f}  {std_abs:10.4f}  "
              f"{mean_signed:10.4f}  {zero_pct:8.2f}  {median_abs:12.1f}")

    # Group by skip
    skip_groups = defaultdict(list)
    for r in results:
        skip_groups[r['skip']].append(r)

    print(f"\n  Mean |ΔH| by skip distance:")
    for skip in sorted(skip_groups.keys()):
        items = skip_groups[skip]
        mean_abs_avg = np.mean([r['mean_abs'] for r in items])
        zero_avg = np.mean([r['zero_pct'] for r in items])
        print(f"    skip={skip}: mean|ΔH| = {mean_abs_avg:.4f}, ΔH=0 rate = {zero_avg:.2f}%"
              f"  ({len(items)} tiles)")

    # Check for ANISOTROPY: do different tiles with same skip have different spectra?
    print(f"\n  ANISOTROPY CHECK (tiles with same skip):")
    for skip in sorted(skip_groups.keys()):
        items = skip_groups[skip]
        if len(items) > 1:
            abs_vals = [r['mean_abs'] for r in items]
            spread = max(abs_vals) - min(abs_vals)
            print(f"    skip={skip}: mean|ΔH| range = [{min(abs_vals):.4f}, {max(abs_vals):.4f}], "
                  f"spread = {spread:.4f}")
            for r in items:
                print(f"      tile {r['tile']}: mean|ΔH| = {r['mean_abs']:.4f}, "
                      f"ΔH=0 = {r['zero_pct']:.2f}%")

    return results


# ─────────────────────────────────────────────
#  2. BOLTZMANN-WEIGHTED CORRELATIONS
# ─────────────────────────────────────────────

def boltzmann_correlations(n, beta_values=[0, 0.1, 0.5, 1.0, 2.0]):
    """Compute <σ_i σ_j>_β for different inverse temperatures.
    At β=0: uniform measure (paramagnetic, all correlations 0).
    At β>0: H-weighted measure (should show AFM correlations).
    At β<0: anti-H-weighted (ferromagnetic preference)."""

    print(f"\n{'='*80}")
    print(f"  BOLTZMANN-WEIGHTED SPIN CORRELATIONS — n = {n}")
    print(f"{'='*80}")

    tiles = staircase_tiles(n)
    m = len(tiles)
    adj = tile_neighbors(tiles, n)

    # Compute lattice distances
    from collections import deque
    dist_matrix = [[float('inf')]*m for _ in range(m)]
    for i in range(m):
        dist_matrix[i][i] = 0
        q = deque([i])
        while q:
            v = q.popleft()
            for w in adj[v]:
                if dist_matrix[i][w] > dist_matrix[i][v] + 1:
                    dist_matrix[i][w] = dist_matrix[i][v] + 1
                    q.append(w)

    # Collect all tournaments with H and spins
    tournaments = []
    for mask, A in all_tournaments(n):
        H = compute_H(A, n)
        spins = tournament_to_spins(A, n, tiles)
        tournaments.append((H, spins))

    N = len(tournaments)
    print(f"  {N} tournaments, {m} tiles")

    for beta in beta_values:
        print(f"\n  β = {beta:.2f}:")

        # Compute Boltzmann weights
        H_arr = np.array([t[0] for t in tournaments], dtype=np.float64)
        if beta != 0:
            log_weights = beta * H_arr
            log_weights -= np.max(log_weights)  # Prevent overflow
            weights = np.exp(log_weights)
        else:
            weights = np.ones(N)
        weights /= weights.sum()

        # Compute weighted <σ_i σ_j>
        corr_matrix = np.zeros((m, m))
        for idx, (H, spins) in enumerate(tournaments):
            w = weights[idx]
            for i in range(m):
                for j in range(m):
                    corr_matrix[i][j] += w * spins[i] * spins[j]

        # Average by distance
        max_d = max(dist_matrix[i][j] for i in range(m) for j in range(m)
                     if dist_matrix[i][j] < float('inf'))
        corr_by_dist = defaultdict(list)
        for i in range(m):
            for j in range(m):
                d = int(dist_matrix[i][j])
                if d < float('inf'):
                    corr_by_dist[d].append(corr_matrix[i][j])

        print(f"    {'dist':>5}  {'<σσ>_β':>12}  {'sign':>6}")
        for d in range(int(max_d) + 1):
            vals = corr_by_dist[d]
            mean_c = np.mean(vals)
            sign = "+" if mean_c > 0.005 else ("-" if mean_c < -0.005 else "~0")
            print(f"    {d:5d}  {mean_c:12.6f}  {sign:>6}")

        # Also compute average magnetization and staggered magnetization
        mean_M = sum(weights[idx] * sum(tournaments[idx][1]) for idx in range(N))
        mean_M2 = sum(weights[idx] * sum(tournaments[idx][1])**2 for idx in range(N))
        mean_H = np.sum(weights * H_arr)
        mean_H2 = np.sum(weights * H_arr**2)
        chi_H = beta * (mean_H2 - mean_H**2)  # Susceptibility
        print(f"    <M> = {mean_M:.6f}, <M²> = {mean_M2:.4f}, <H> = {mean_H:.4f}")
        print(f"    χ_H = β·Var(H) = {chi_H:.4f}")


# ─────────────────────────────────────────────
#  3. VERTEX-TRIPLE FRUSTRATION ENERGY
# ─────────────────────────────────────────────

def frustration_energy_analysis(n):
    """Define a 'frustration energy' based on vertex triples, not tile pairs.

    E_frust(T) = number of TRANSITIVE triples (unfrustrated)
                = C(n,3) - c₃(T)

    This is the 'true' AFM energy: it counts the number of 'satisfied' (ordered)
    triples. The AFM ground state MAXIMIZES c₃ (frustration) = MINIMIZES E_frust.

    Also define the ALTERNATING energy:
    E_alt(T) = Σ_{(a,b) tile} (-1)^{a+b} σ_{a,b}
    This measures the staggered magnetization on the staircase.
    """
    print(f"\n{'='*80}")
    print(f"  VERTEX-TRIPLE FRUSTRATION ENERGY — n = {n}")
    print(f"{'='*80}")

    tiles = staircase_tiles(n)
    m = len(tiles)

    data = []
    for mask, A in all_tournaments(n):
        H = compute_H(A, n)
        c3 = count_3cycles(A, n)
        sc = scores(A, n)
        sv = sum((s - (n-1)/2)**2 for s in sc) / n
        spins = tournament_to_spins(A, n, tiles)

        # Frustration energy
        E_frust = comb(n, 3) - c3  # Number of transitive triples

        # Staggered magnetization on staircase
        M_stag = sum((-1)**(a+b) * s for (a, b), s in zip(tiles, spins))

        # "Topological charge" — sum of spins weighted by position
        Q = sum(s * (a - b) for (a, b), s in zip(tiles, spins))  # Range-weighted

        data.append({
            'H': H, 'c3': c3, 'sv': sv, 'E_frust': E_frust,
            'M_stag': M_stag, 'Q': Q, 'scores': sc
        })

    N = len(data)
    H_arr = np.array([d['H'] for d in data])
    c3_arr = np.array([d['c3'] for d in data])
    Ef_arr = np.array([d['E_frust'] for d in data])
    Ms_arr = np.array([d['M_stag'] for d in data])
    Q_arr = np.array([d['Q'] for d in data])
    sv_arr = np.array([d['sv'] for d in data])

    print(f"\n  Correlations:")
    print(f"    corr(E_frust, H) = {np.corrcoef(Ef_arr, H_arr)[0,1]:.6f}")
    print(f"    corr(c₃, H)     = {np.corrcoef(c3_arr, H_arr)[0,1]:.6f}")
    print(f"    corr(M_stag, H) = {np.corrcoef(Ms_arr, H_arr)[0,1]:.6f}")
    print(f"    corr(Q, H)      = {np.corrcoef(Q_arr, H_arr)[0,1]:.6f}")
    print(f"    corr(|M_stag|, H) = {np.corrcoef(np.abs(Ms_arr), H_arr)[0,1]:.6f}")
    print(f"    corr(|Q|, H)    = {np.corrcoef(np.abs(Q_arr), H_arr)[0,1]:.6f}")

    # The KEY question: is there an energy functional E(T) such that
    # H(T) = f(E(T)) for some simple function f?

    # Test: H = a * c₃² + b * c₃ + c ?
    X = np.column_stack([c3_arr**2, c3_arr, np.ones(N)])
    coeffs, _, _, _ = np.linalg.lstsq(X, H_arr, rcond=None)
    H_pred = X @ coeffs
    R2_quad = 1 - np.sum((H_arr - H_pred)**2) / np.sum((H_arr - H_arr.mean())**2)
    print(f"\n  Quadratic model: H ≈ {coeffs[0]:.6f}*c₃² + {coeffs[1]:.6f}*c₃ + {coeffs[2]:.6f}")
    print(f"  R²(quadratic in c₃) = {R2_quad:.6f}")

    # Test: does adding 5-cycle count help?
    # At n=5: count 5-cycles
    if n <= 6:
        c5_arr = np.zeros(N)
        for idx in range(N):
            A = [[0]*n for _ in range(n)]
            # Reconstruct from data... need the actual tournament
            pass
        # Actually let's compute c5 inline
        c5_data = []
        for mask, A in all_tournaments(n):
            c5 = 0
            for combo in combinations(range(n), 5):
                verts = list(combo)
                for perm in permutations(verts):
                    if all(A[perm[i]][perm[(i+1) % 5]] for i in range(5)):
                        c5 += 1
                # Each 5-cycle counted 10 times (5 rotations × 2 directions... no, only 5 rotations for directed)
            c5 //= 5  # 5 rotations
            c5_data.append(c5)
        # This is too slow for n=6 (32768 × C(6,5) × 5! = too much)
        # Skip for now

    # The OCF decomposition: H = 1 + 2α₁ + 4α₂ + 8α₃ + ...
    # At n=4: α₁ = c₃, α₂ = 0, so H = 1 + 2c₃ EXACTLY
    # At n=5: need to compute α₁ = c₃ + c₅ and α₂
    print(f"\n  OCF decomposition analysis:")
    print(f"  H = 1 + 2α₁ + 4α₂ + 8α₃ + ...")
    print(f"  (H-1)/2 = α₁ + 2α₂ + 4α₃ + ...")

    # Group by (c₃, H) to see how much variance α₂ explains
    c3_H_groups = defaultdict(list)
    for d in data:
        c3_H_groups[(d['c3'], d['H'])].append(d)

    print(f"\n  (c₃, H) table:")
    print(f"  {'c₃':>4}  {'H':>4}  {'count':>6}  {'(H-1)/2':>8}  {'α₁_min':>7}")
    for (c3, H) in sorted(c3_H_groups.keys()):
        count = len(c3_H_groups[(c3, H)])
        h_half = (H - 1) // 2
        # α₁ ≥ c₃ (since every 3-cycle is an odd cycle)
        print(f"  {c3:4d}  {H:4d}  {count:6d}  {h_half:8d}  {c3:7d}")

    return data


# ─────────────────────────────────────────────
#  4. THE ANTIFERROMAGNETIC GROUND STATE
# ─────────────────────────────────────────────

def afm_ground_state_analysis(n):
    """Identify the AFM ground state: the tournament(s) that maximize frustration (c₃)
    while minimizing score variance.

    In condensed matter: AFM ground state = maximum Néel order.
    In tournaments: regular tournaments with maximum c₃ = regular with max H.

    Question: Is the H-maximizer always the AFM ground state?"""

    print(f"\n{'='*80}")
    print(f"  AFM GROUND STATE ANALYSIS — n = {n}")
    print(f"{'='*80}")

    max_c3 = 0
    max_H = 0
    min_sv = float('inf')

    data_by_H = defaultdict(list)
    data_by_c3 = defaultdict(list)

    for mask, A in all_tournaments(n):
        H = compute_H(A, n)
        c3 = count_3cycles(A, n)
        sc = scores(A, n)
        sv = sum((s - (n-1)/2)**2 for s in sc) / n
        max_c3 = max(max_c3, c3)
        max_H = max(max_H, H)
        min_sv = min(min_sv, sv)
        data_by_H[H].append({'c3': c3, 'sv': sv, 'scores': sc, 'mask': mask})
        data_by_c3[c3].append({'H': H, 'sv': sv, 'scores': sc, 'mask': mask})

    print(f"\n  Maximum c₃ = {max_c3}")
    print(f"  Maximum H = {max_H}")
    print(f"  Minimum score variance = {min_sv:.4f}")

    # H-maximizers
    H_max_data = data_by_H[max_H]
    print(f"\n  H-MAXIMIZERS (H = {max_H}): {len(H_max_data)} tournaments")
    c3_vals = [d['c3'] for d in H_max_data]
    sv_vals = [d['sv'] for d in H_max_data]
    score_seqs = set(d['scores'] for d in H_max_data)
    print(f"    c₃ values: {sorted(set(c3_vals))}")
    print(f"    Score variances: {sorted(set(round(v, 4) for v in sv_vals))}")
    print(f"    Score sequences: {sorted(score_seqs)}")

    # c₃-maximizers
    c3_max_data = data_by_c3[max_c3]
    print(f"\n  FRUSTRATION-MAXIMIZERS (c₃ = {max_c3}): {len(c3_max_data)} tournaments")
    H_vals = [d['H'] for d in c3_max_data]
    sv_vals = [d['sv'] for d in c3_max_data]
    score_seqs = set(d['scores'] for d in c3_max_data)
    print(f"    H values: {sorted(set(H_vals))}")
    print(f"    Score variances: {sorted(set(round(v, 4) for v in sv_vals))}")
    print(f"    Score sequences: {sorted(score_seqs)}")

    # Key question: do H-max and c₃-max overlap?
    H_max_masks = set(d['mask'] for d in H_max_data)
    c3_max_masks = set(d['mask'] for d in c3_max_data)
    overlap = H_max_masks & c3_max_masks

    print(f"\n  OVERLAP between H-max and c₃-max: {len(overlap)} tournaments")
    if len(overlap) > 0:
        print(f"    → The H-maximizer IS (sometimes) the maximally frustrated tournament!")
    else:
        print(f"    → The H-maximizer is NOT the maximally frustrated tournament.")
        # Which is more frustrated?
        H_max_c3 = max(c3_vals)
        c3_max_H = max(H_vals)
        print(f"    H-max has c₃ = {H_max_c3} (out of max {max_c3})")
        print(f"    c₃-max has H = {c3_max_H} (out of max {max_H})")

    # Score-variance minimizers (most regular)
    min_sv_data = [(c3, d) for c3 in data_by_c3 for d in data_by_c3[c3] if abs(d['sv'] - min_sv) < 0.001]
    print(f"\n  MOST REGULAR tournaments (sv = {min_sv:.4f}): {len(min_sv_data)}")
    H_vals_reg = [d['H'] for _, d in min_sv_data]
    c3_vals_reg = [c3 for c3, _ in min_sv_data]
    print(f"    H range: {min(H_vals_reg)} to {max(H_vals_reg)}")
    print(f"    c₃ range: {min(c3_vals_reg)} to {max(c3_vals_reg)}")

    return data_by_H, data_by_c3


# ─────────────────────────────────────────────
#  5. SPECIFIC HEAT AND SUSCEPTIBILITY
# ─────────────────────────────────────────────

def thermodynamic_functions(n, beta_range=np.linspace(-2, 3, 51)):
    """Compute thermodynamic functions of the tournament ensemble
    with Boltzmann weight exp(β H(T)).

    The tournament 'partition function' Z(β) = Σ_T exp(β H(T)).
    Free energy: F(β) = -log(Z) / β
    Internal energy: U(β) = <H>_β
    Specific heat: C(β) = β² Var(H)_β
    Susceptibility: χ(β) = β Var(c₃)_β

    Phase transition: where C(β) peaks."""

    print(f"\n{'='*80}")
    print(f"  THERMODYNAMIC FUNCTIONS — n = {n}")
    print(f"{'='*80}")

    # Collect all data
    H_list = []
    c3_list = []
    sv_list = []
    for mask, A in all_tournaments(n):
        H_list.append(compute_H(A, n))
        c3_list.append(count_3cycles(A, n))
        sv_list.append(score_variance(A, n))

    H_arr = np.array(H_list, dtype=np.float64)
    c3_arr = np.array(c3_list, dtype=np.float64)
    sv_arr = np.array(sv_list, dtype=np.float64)
    N = len(H_arr)

    print(f"  {N} tournaments")
    print(f"  H range: [{H_arr.min():.0f}, {H_arr.max():.0f}], mean = {H_arr.mean():.2f}")
    print(f"  c₃ range: [{c3_arr.min():.0f}, {c3_arr.max():.0f}], mean = {c3_arr.mean():.2f}")

    print(f"\n  {'β':>8}  {'<H>':>10}  {'Var(H)':>10}  {'C(β)':>10}  "
          f"{'<c₃>':>8}  {'<sv>':>8}  {'F':>12}")

    max_C = 0
    beta_max_C = 0

    for beta in beta_range:
        log_w = beta * H_arr
        log_w -= np.max(log_w)
        w = np.exp(log_w)
        Z = w.sum()
        w /= Z

        mean_H = np.sum(w * H_arr)
        var_H = np.sum(w * (H_arr - mean_H)**2)
        C = beta**2 * var_H if abs(beta) > 1e-10 else 0
        mean_c3 = np.sum(w * c3_arr)
        mean_sv = np.sum(w * sv_arr)
        F = -np.log(Z) / beta if abs(beta) > 1e-10 else 0

        if C > max_C:
            max_C = C
            beta_max_C = beta

        if abs(beta) < 0.01 or abs(beta - 1.0) < 0.06 or abs(beta - 2.0) < 0.06 or \
           abs(beta - beta_range[0]) < 0.01 or abs(beta - beta_range[-1]) < 0.01 or \
           C > max_C * 0.8:
            print(f"  {beta:8.3f}  {mean_H:10.4f}  {var_H:10.4f}  {C:10.4f}  "
                  f"{mean_c3:8.3f}  {mean_sv:8.4f}  {F:12.4f}")

    print(f"\n  SPECIFIC HEAT PEAK: C_max = {max_C:.4f} at β = {beta_max_C:.3f}")
    print(f"  This is the PHASE TRANSITION TEMPERATURE of the tournament ensemble.")
    print(f"  At β > β_c: system in 'ordered' (high H, high c₃, low sv) phase")
    print(f"  At β < β_c: system in 'disordered' (low H, low c₃, high sv) phase")

    # Negative beta: anti-H weighting (prefers transitive/ferromagnetic)
    print(f"\n  At β = -2 (ferromagnetic preference):")
    log_w = -2 * H_arr
    log_w -= np.max(log_w)
    w = np.exp(log_w)
    w /= w.sum()
    print(f"    <H> = {np.sum(w * H_arr):.4f}, <c₃> = {np.sum(w * c3_arr):.4f}, <sv> = {np.sum(w * sv_arr):.4f}")

    print(f"\n  At β = +2 (antiferromagnetic preference):")
    log_w = 2 * H_arr
    log_w -= np.max(log_w)
    w = np.exp(log_w)
    w /= w.sum()
    print(f"    <H> = {np.sum(w * H_arr):.4f}, <c₃> = {np.sum(w * c3_arr):.4f}, <sv> = {np.sum(w * sv_arr):.4f}")


# ─────────────────────────────────────────────
#  6. FORBIDDEN VALUES AND SPECTRAL GAPS
# ─────────────────────────────────────────────

def forbidden_analysis_deep(n):
    """Deep analysis of H=7 gap through the AFM lens.

    The transfer matrix M = [[1,2,0],[0,0,1],[1,1,0]] has
    ζ_M(-3) = 7 and ζ_M(-5) = 21.

    In the AFM picture: forbidden values correspond to energy levels
    that CANNOT be achieved by any spin configuration on the staircase.
    This is a GEOMETRIC constraint, not an energetic one.

    Key observation: H = 1 + 2Σ_{k≥1} 2^{k-1} α_k
    H = 7 requires α₁ = 3, α_k = 0 for k ≥ 2 (since 7 = 1 + 2·3)
    This means EXACTLY 3 odd cycles, NO pair vertex-disjoint.

    In AFM terms: 3 frustrated triangles that all share vertices.
    On a small lattice, this is a TOPOLOGICAL CONSTRAINT."""

    print(f"\n{'='*80}")
    print(f"  H=7 FORBIDDEN VALUE: AFM TOPOLOGICAL ANALYSIS — n = {n}")
    print(f"{'='*80}")

    # For each tournament, compute the OCF decomposition H = 1 + 2α₁ + 4α₂ + ...
    ocf_data = []

    for mask, A in all_tournaments(n):
        H = compute_H(A, n)
        c3 = count_3cycles(A, n)

        # Count independent sets in Ω(T) at various sizes
        # First, find all directed odd cycles
        odd_cycles = find_odd_cycles(A, n)
        alpha_1 = len(odd_cycles)

        # Build conflict graph (share vertex = adjacent)
        conflict_adj = build_conflict_graph(odd_cycles)

        # Count independent sets by size
        alphas = count_independent_sets(odd_cycles, conflict_adj)

        # Verify OCF
        H_ocf = sum(2**k * alphas.get(k, 0) for k in range(len(odd_cycles)+1))
        assert H_ocf == H, f"OCF failed: H={H}, OCF={H_ocf}, alphas={alphas}"

        ocf_data.append({
            'H': H, 'c3': c3, 'alpha_1': alpha_1,
            'alpha_2': alphas.get(2, 0),
            'alpha_3': alphas.get(3, 0),
            'alphas': alphas, 'scores': scores(A, n)
        })

    # Group by H
    by_H = defaultdict(list)
    for d in ocf_data:
        by_H[d['H']].append(d)

    print(f"\n  OCF decomposition for each achievable H:")
    print(f"  {'H':>4}  {'count':>6}  {'α₁':>4}  {'α₂':>4}  {'α₃':>4}  {'formula':>20}")
    for H in sorted(by_H.keys()):
        items = by_H[H]
        a1s = sorted(set(d['alpha_1'] for d in items))
        a2s = sorted(set(d['alpha_2'] for d in items))
        a3s = sorted(set(d['alpha_3'] for d in items))

        formula = f"1+2·{a1s[0]}+4·{a2s[0]}"
        if a3s[0] > 0:
            formula += f"+8·{a3s[0]}"
        print(f"  {H:4d}  {len(items):6d}  {a1s}  {a2s}  {a3s}  {formula}")

    # Why H=7 is impossible: need α₁=3, α₂=0
    print(f"\n  WHY H=7 IS IMPOSSIBLE (AFM analysis):")
    print(f"  H=7 requires: 1 + 2α₁ + 4α₂ + ... = 7")
    print(f"  → α₁ = 3, α₂ = 0 (since 4α₂ would push H ≥ 11)")
    print(f"  → Exactly 3 odd cycles, ALL sharing vertices (α₂=0)")

    # Check: is α₁=3 achievable?
    a1_3_data = [d for d in ocf_data if d['alpha_1'] == 3]
    if a1_3_data:
        print(f"  α₁=3 IS achievable: {len(a1_3_data)} tournaments")
        a2_vals = set(d['alpha_2'] for d in a1_3_data)
        print(f"  But when α₁=3: α₂ values = {sorted(a2_vals)}")
        if 0 in a2_vals:
            # α₁=3, α₂=0 achievable → H would be 1+6=7
            a1_3_a2_0 = [d for d in a1_3_data if d['alpha_2'] == 0]
            print(f"  α₁=3, α₂=0 combinations: {len(a1_3_a2_0)}")
            print(f"  Their H values: {sorted(set(d['H'] for d in a1_3_a2_0))}")
        else:
            print(f"  α₁=3, α₂=0 is NEVER achieved! → H=7 impossible")
            print(f"  WHY: when you have 3 odd cycles, at least 2 must be vertex-disjoint")
            print(f"  This is the TOPOLOGICAL CONSTRAINT on the cycle packing")
    else:
        print(f"  α₁=3 is NOT achievable at n={n}")
        print(f"  Achievable α₁ values: {sorted(set(d['alpha_1'] for d in ocf_data))}")

    return ocf_data


def find_odd_cycles(A, n):
    """Find all directed odd cycles in tournament A."""
    cycles = []
    for length in range(3, n+1, 2):  # 3, 5, 7, ...
        for combo in combinations(range(n), length):
            verts = list(combo)
            # Check all cyclic orderings (directed)
            # For efficiency, try only one canonical rotation
            for perm in permutations(verts):
                if perm[0] == min(perm):  # Canonical: start with smallest vertex
                    if all(A[perm[i]][perm[(i+1) % length]] for i in range(length)):
                        cycles.append(frozenset(perm))
    # Remove duplicates (same vertex set, different starting point)
    unique = set()
    result = []
    for c in cycles:
        if c not in unique:
            unique.add(c)
            result.append(c)
    return result


def build_conflict_graph(cycles):
    """Build conflict graph: two cycles are adjacent iff they share a vertex."""
    adj = defaultdict(set)
    for i in range(len(cycles)):
        for j in range(i+1, len(cycles)):
            if cycles[i] & cycles[j]:  # Share a vertex
                adj[i].add(j)
                adj[j].add(i)
    return adj


def count_independent_sets(cycles, adj):
    """Count independent sets by size using brute force (small n only)."""
    m = len(cycles)
    alphas = defaultdict(int)
    for mask in range(1 << m):
        # Check if mask is an independent set
        subset = [i for i in range(m) if mask & (1 << i)]
        is_indep = True
        for i in range(len(subset)):
            for j in range(i+1, len(subset)):
                if subset[j] in adj.get(subset[i], set()):
                    is_indep = False
                    break
            if not is_indep:
                break
        if is_indep:
            alphas[len(subset)] += 1
    return dict(alphas)


# ─────────────────────────────────────────────
#  MAIN
# ─────────────────────────────────────────────

if __name__ == "__main__":
    print("╔══════════════════════════════════════════════════════════════════╗")
    print("║  THE ANTIFERROMAGNETIC TOURNAMENT — DEEP ANALYSIS              ║")
    print("║  opus-2026-04-04-S15                                           ║")
    print("╚══════════════════════════════════════════════════════════════════╝")

    t0 = time.time()

    # 1. Magnon dispersion at n=6
    magnon_dispersion_n6()

    t1 = time.time()
    print(f"\n  [Magnon dispersion took {t1-t0:.1f}s]")

    # 2. Boltzmann correlations at n=5
    boltzmann_correlations(5, beta_values=[0, 0.1, 0.5, 1.0, 2.0, 5.0])

    t2 = time.time()
    print(f"\n  [Boltzmann correlations took {t2-t1:.1f}s]")

    # 3. Frustration energy at n=5
    frustration_energy_analysis(5)

    # 4. AFM ground state at n=5, n=6
    for nn in [5, 6]:
        afm_ground_state_analysis(nn)

    t3 = time.time()
    print(f"\n  [Ground state analysis took {t3-t2:.1f}s]")

    # 5. Thermodynamic functions at n=5
    thermodynamic_functions(5)

    t4 = time.time()
    print(f"\n  [Thermodynamics took {t4-t3:.1f}s]")

    # 6. Forbidden value analysis at n=5
    forbidden_analysis_deep(5)

    t5 = time.time()
    print(f"\n  [Forbidden analysis took {t5-t4:.1f}s]")

    print(f"\n\n{'='*80}")
    print(f"  TOTAL TIME: {t5-t0:.1f}s")
    print(f"{'='*80}")
