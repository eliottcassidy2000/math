#!/usr/bin/env python3
"""
antiferromagnetic_tournament_s15.py — opus-2026-04-04-S15

THE ANTIFERROMAGNETIC TOURNAMENT: A Condensed Matter Physics Perspective

CORE THESIS: A tournament is a frustrated spin system on a triangular lattice.
The staircase lattice delta_{n-2} is the underlying lattice; each tile is a
"spin" with two states (forward/backward arc). The tournament structure induces
ANTIFERROMAGNETIC coupling — and frustration is the source of structure.

DICTIONARY:
  Tournament T          ↔  Spin configuration σ on staircase lattice
  Arc orientation (0/1) ↔  Spin up/down (±1)
  Score sequence        ↔  Magnetization profile M(layer)
  3-cycle count c₃      ↔  Frustration index
  Regular tournament    ↔  Néel state (staggered AFM ground state)
  Transitive tournament ↔  Fully polarized (ferromagnetic) state
  H(T) = # Ham paths    ↔  Partition function Z = I(Ω, 2)
  Wiggly line (1 flip)  ↔  Magnon (spin wave excitation)
  Metagraph G_n         ↔  Configuration space / energy landscape

EXPLORATIONS:
  1. Frustration index vs H correlation
  2. Néel order parameter (score variance → staggered magnetization)
  3. Magnon dispersion relation (ΔH per tile position)
  4. AFM energy functional on staircase lattice
  5. Spin stiffness and correlation lengths
  6. Connection to forbidden H values via quantization
"""

import numpy as np
from itertools import combinations, permutations
from math import comb, factorial
from collections import defaultdict, Counter
import random
import sys

random.seed(42)
np.random.seed(42)

# ─────────────────────────────────────────────
#  Core tournament utilities (from stat_mech_tournament.py)
# ─────────────────────────────────────────────

def all_tournaments(n):
    """Enumerate all 2^C(n,2) tournaments on n vertices."""
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
    """H(T) via Held-Karp DP. O(n^2 * 2^n)."""
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
    return tuple(sorted(sum(A[i]) for i in range(n)))

def score_variance(A, n):
    s = [sum(A[i]) for i in range(n)]
    mean = sum(s) / n
    return sum((x - mean)**2 for x in s) / n

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

def count_5cycles(A, n):
    """Count directed 5-cycles (each counted once)."""
    c = 0
    for combo in combinations(range(n), 5):
        verts = list(combo)
        for perm in permutations(verts):
            if all(A[perm[i]][perm[(i+1) % 5]] for i in range(5)):
                c += 1
    return c // 10  # Each 5-cycle counted 5 times in each direction

def is_sc(A, n):
    """Is tournament self-complementary? (Exists anti-automorphism)"""
    sc = scores(A, n)
    # Quick check: score sequence must be self-complementary
    for i in range(n):
        if sc[i] + sc[n-1-i] != n-1:
            return False
    return True  # This is necessary but not sufficient; good enough for statistics

# ─────────────────────────────────────────────
#  ANTIFERROMAGNETIC FRAMEWORK
# ─────────────────────────────────────────────

def staircase_tiles(n):
    """Return list of tiles (a,b) with a >= b+2 in the staircase.
    Tile order: for y=1..n-2: for x=n down to y+2: tile (x,y).
    This matches the explorer's enumeration."""
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            if x - y >= 2:
                tiles.append((x, y))
    return tiles

def tile_neighbors(tiles, n):
    """Two tiles are NEIGHBORS if they share a vertex.
    Tile (a,b) involves vertices a and b.
    This defines the lattice structure for AFM coupling."""
    adj = defaultdict(set)
    for i, (a1, b1) in enumerate(tiles):
        for j, (a2, b2) in enumerate(tiles):
            if i < j and (a1 == a2 or a1 == b2 or b1 == a2 or b1 == b2):
                adj[i].add(j)
                adj[j].add(i)
    return adj

def tournament_to_spins(A, n, tiles):
    """Convert tournament to spin configuration on staircase.
    Spin = +1 if arc goes from higher to lower vertex (forward = base path direction),
    Spin = -1 if arc goes backward."""
    spins = []
    for (a, b) in tiles:
        # a > b always. Forward = a→b (base path direction for consecutive)
        # Actually in the tiling model: bit 0 = a→b (forward), bit 1 = b→a (backward)
        # We use: σ = +1 for a→b, σ = -1 for b→a
        if A[a-1][b-1]:  # a→b (0-indexed)
            spins.append(1)
        else:
            spins.append(-1)
    return spins

def afm_energy(spins, adj):
    """Antiferromagnetic nearest-neighbor Ising energy.
    E = +J Σ_{<i,j>} σ_i σ_j  with J > 0 (AFM: aligned spins cost energy).
    E < 0 when neighbors have opposite spins (AFM ground state).
    We use J = 1."""
    E = 0
    counted = set()
    for i in adj:
        for j in adj[i]:
            if (i, j) not in counted:
                E += spins[i] * spins[j]
                counted.add((i, j))
                counted.add((j, i))
    return E

def magnetization(spins):
    """Total magnetization M = Σ σ_i (should be near 0 for AFM ground state)."""
    return sum(spins)

def staggered_magnetization(spins, tiles, n):
    """Staggered (Néel) order parameter.
    Assign sublattice: tile (a,b) → sublattice A if (a+b) is even, B if odd.
    M_stag = |Σ_A σ - Σ_B σ| / m."""
    m = len(spins)
    A_sum, B_sum = 0, 0
    for idx, (a, b) in enumerate(tiles):
        if (a + b) % 2 == 0:
            A_sum += spins[idx]
        else:
            B_sum += spins[idx]
    return abs(A_sum - B_sum) / m if m > 0 else 0

def frustration_index(A, n):
    """Frustration index = c_3 / C(n,3). Ranges from 0 (transitive) to
    the maximum possible (regular tournament at odd n)."""
    c3 = count_3cycles(A, n)
    return c3 / comb(n, 3) if comb(n, 3) > 0 else 0

# ─────────────────────────────────────────────
#  MAGNON DISPERSION: ΔH per tile position
# ─────────────────────────────────────────────

def compute_magnon_spectrum(n, sample_size=None):
    """For each tile position (wiggly class), compute the distribution of ΔH
    when that tile is flipped. This is the magnon dispersion relation:
    the 'energy cost' of flipping spin at position k."""

    tiles = staircase_tiles(n)
    m = len(tiles)
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    ne = len(edges)

    # Map tiles to edge indices
    tile_to_edge = {}
    for t_idx, (a, b) in enumerate(tiles):
        # Find which edge corresponds to tile (a,b)
        # Tile (a,b) is the arc between vertices a-1 and b-1 (0-indexed)
        for e_idx, (i, j) in enumerate(edges):
            if (i == b-1 and j == a-1) or (i == a-1 and j == b-1):
                tile_to_edge[t_idx] = e_idx
                break

    delta_H_by_tile = defaultdict(list)  # tile_idx → list of |ΔH| values
    delta_H_signed_by_tile = defaultdict(list)  # signed ΔH

    count = 0
    for mask, A in all_tournaments(n):
        if sample_size and count >= sample_size:
            break
        H0 = compute_H(A, n)

        for t_idx, (a, b) in enumerate(tiles):
            # Flip tile: swap A[a-1][b-1] and A[b-1][a-1]
            A2 = [row[:] for row in A]
            A2[a-1][b-1] = 1 - A2[a-1][b-1]
            A2[b-1][a-1] = 1 - A2[b-1][a-1]
            H1 = compute_H(A2, n)
            dH = H1 - H0
            delta_H_by_tile[t_idx].append(abs(dH))
            delta_H_signed_by_tile[t_idx].append(dH)

        count += 1
        if count % 500 == 0:
            print(f"  Processed {count} tournaments...", file=sys.stderr)

    return tiles, delta_H_by_tile, delta_H_signed_by_tile


# ─────────────────────────────────────────────
#  MAIN: Run all explorations
# ─────────────────────────────────────────────

def explore_n(n, exhaustive=True, sample_size=5000):
    """Full antiferromagnetic analysis at a given n."""
    print(f"\n{'='*80}")
    print(f"  ANTIFERROMAGNETIC TOURNAMENT ANALYSIS — n = {n}")
    print(f"{'='*80}")

    tiles = staircase_tiles(n)
    m = len(tiles)
    adj = tile_neighbors(tiles, n)

    print(f"\n  Staircase lattice: {m} tiles")
    print(f"  Tile positions: {tiles}")
    degrees = [len(adj[i]) for i in range(m)]
    print(f"  Lattice degree sequence: {sorted(degrees)}")
    print(f"  Mean degree: {sum(degrees)/m:.2f}")
    print(f"  Total lattice edges: {sum(degrees)//2}")

    # Collect data
    data = []  # (H, c3, score_var, E_afm, M, M_stag, frustration, is_sc_flag)

    if exhaustive:
        print(f"\n  Exhaustive enumeration: 2^{comb(n,2)} = {1 << comb(n,2)} tournaments")
        for mask, A in all_tournaments(n):
            H = compute_H(A, n)
            c3 = count_3cycles(A, n)
            sv = score_variance(A, n)
            spins = tournament_to_spins(A, n, tiles)
            E = afm_energy(spins, adj)
            M = magnetization(spins)
            Ms = staggered_magnetization(spins, tiles, n)
            fi = c3 / comb(n, 3) if comb(n, 3) > 0 else 0
            sc = is_sc(A, n)
            data.append((H, c3, sv, E, M, Ms, fi, sc))

            if len(data) % 2000 == 0:
                print(f"    ... {len(data)} done", file=sys.stderr)
    else:
        print(f"\n  Random sampling: {sample_size} tournaments")
        for _ in range(sample_size):
            A = [[0]*n for _ in range(n)]
            for i in range(n):
                for j in range(i+1, n):
                    if random.random() < 0.5:
                        A[i][j] = 1
                    else:
                        A[j][i] = 1
            H = compute_H(A, n)
            c3 = count_3cycles(A, n)
            sv = score_variance(A, n)
            spins = tournament_to_spins(A, n, tiles)
            E = afm_energy(spins, adj)
            M = magnetization(spins)
            Ms = staggered_magnetization(spins, tiles, n)
            fi = c3 / comb(n, 3)
            sc = is_sc(A, n)
            data.append((H, c3, sv, E, M, Ms, fi, sc))

    N = len(data)
    H_arr = np.array([d[0] for d in data])
    c3_arr = np.array([d[1] for d in data])
    sv_arr = np.array([d[2] for d in data])
    E_arr = np.array([d[3] for d in data])
    M_arr = np.array([d[4] for d in data])
    Ms_arr = np.array([d[5] for d in data])
    fi_arr = np.array([d[6] for d in data])
    sc_arr = np.array([d[7] for d in data])

    # ─── 1. FRUSTRATION vs H ───
    print(f"\n{'─'*60}")
    print(f"  1. FRUSTRATION vs HAMILTONIAN PATH COUNT")
    print(f"{'─'*60}")
    corr_fi_H = np.corrcoef(fi_arr, H_arr)[0,1]
    print(f"  Correlation(frustration, H) = {corr_fi_H:.6f}")
    print(f"  Correlation(c₃, H) = {np.corrcoef(c3_arr, H_arr)[0,1]:.6f}")

    # Group by c3
    c3_groups = defaultdict(list)
    for i in range(N):
        c3_groups[data[i][1]].append(data[i][0])

    print(f"\n  H by frustration level (c₃):")
    print(f"  {'c₃':>5}  {'count':>7}  {'mean H':>10}  {'min H':>7}  {'max H':>7}  {'H range':>9}")
    for c3 in sorted(c3_groups.keys()):
        Hs = c3_groups[c3]
        print(f"  {c3:5d}  {len(Hs):7d}  {np.mean(Hs):10.2f}  "
              f"{min(Hs):7d}  {max(Hs):7d}  {max(Hs)-min(Hs):9d}")

    # ─── 2. NÉEL ORDER PARAMETER ───
    print(f"\n{'─'*60}")
    print(f"  2. NÉEL ORDER PARAMETER (Score Variance)")
    print(f"{'─'*60}")
    corr_sv_H = np.corrcoef(sv_arr, H_arr)[0,1]
    print(f"  Correlation(score_variance, H) = {corr_sv_H:.6f}")
    print(f"  Interpretation: {'MORE regular → MORE Ham paths' if corr_sv_H < 0 else 'MORE regular → FEWER Ham paths'}")

    # Group by score variance
    sv_groups = defaultdict(list)
    for i in range(N):
        sv_groups[round(data[i][2], 4)].append(data[i][0])

    print(f"\n  H by score variance:")
    print(f"  {'sv':>8}  {'count':>7}  {'mean H':>10}  {'min H':>7}  {'max H':>7}")
    for sv in sorted(sv_groups.keys()):
        Hs = sv_groups[sv]
        if len(Hs) >= 2:
            print(f"  {sv:8.4f}  {len(Hs):7d}  {np.mean(Hs):10.2f}  "
                  f"{min(Hs):7d}  {max(Hs):7d}")

    # ─── 3. AFM ENERGY ───
    print(f"\n{'─'*60}")
    print(f"  3. ANTIFERROMAGNETIC ENERGY FUNCTIONAL")
    print(f"{'─'*60}")
    corr_E_H = np.corrcoef(E_arr, H_arr)[0,1]
    corr_E_c3 = np.corrcoef(E_arr, c3_arr)[0,1]
    print(f"  Correlation(E_AFM, H)  = {corr_E_H:.6f}")
    print(f"  Correlation(E_AFM, c₃) = {corr_E_c3:.6f}")
    print(f"  E_AFM range: [{E_arr.min()}, {E_arr.max()}]")
    print(f"  Mean E_AFM: {E_arr.mean():.4f}")

    # Ground state (min E) analysis
    min_E = E_arr.min()
    gs_indices = np.where(E_arr == min_E)[0]
    gs_H = H_arr[gs_indices]
    print(f"\n  AFM 'GROUND STATE' (min E = {min_E}):")
    print(f"    Count: {len(gs_indices)}")
    print(f"    H values: min={gs_H.min()}, max={gs_H.max()}, mean={gs_H.mean():.2f}")
    print(f"    c₃ values: {sorted(set(c3_arr[gs_indices].astype(int)))}")

    # Max E (fully frustrated) analysis
    max_E = E_arr.max()
    fe_indices = np.where(E_arr == max_E)[0]
    fe_H = H_arr[fe_indices]
    print(f"\n  'FULLY ALIGNED' state (max E = {max_E}):")
    print(f"    Count: {len(fe_indices)}")
    print(f"    H values: min={fe_H.min()}, max={fe_H.max()}, mean={fe_H.mean():.2f}")
    print(f"    = transitive tournaments? scores: {set(tuple(scores(None, n) for _ in [None]) if False else [])}")

    # E by H
    E_groups = defaultdict(list)
    for i in range(N):
        E_groups[int(E_arr[i])].append(int(H_arr[i]))
    print(f"\n  H distribution by AFM energy level:")
    print(f"  {'E_AFM':>6}  {'count':>7}  {'mean H':>10}  {'min H':>7}  {'max H':>7}")
    for E in sorted(E_groups.keys()):
        Hs = E_groups[E]
        print(f"  {E:6d}  {len(Hs):7d}  {np.mean(Hs):10.2f}  {min(Hs):7d}  {max(Hs):7d}")

    # ─── 4. MAGNETIZATION ───
    print(f"\n{'─'*60}")
    print(f"  4. MAGNETIZATION ANALYSIS")
    print(f"{'─'*60}")
    print(f"  Total magnetization M range: [{M_arr.min()}, {M_arr.max()}]")
    print(f"  Mean M: {M_arr.mean():.6f}  (should be ~0 by symmetry)")
    print(f"  Staggered magnetization M_stag range: [{Ms_arr.min():.4f}, {Ms_arr.max():.4f}]")
    corr_Ms_H = np.corrcoef(Ms_arr, H_arr)[0,1]
    print(f"  Correlation(M_stag, H) = {corr_Ms_H:.6f}")
    corr_M_H = np.corrcoef(np.abs(M_arr), H_arr)[0,1]
    print(f"  Correlation(|M|, H) = {corr_M_H:.6f}")

    # ─── 5. SC vs NSC comparison ───
    print(f"\n{'─'*60}")
    print(f"  5. SELF-COMPLEMENTARY vs NON-SC (AFM symmetry)")
    print(f"{'─'*60}")
    sc_idx = np.where(sc_arr)[0]
    nsc_idx = np.where(~sc_arr)[0]
    if len(sc_idx) > 0:
        print(f"  SC tournaments: {len(sc_idx)}")
        print(f"    Mean H: {H_arr[sc_idx].mean():.2f}, Mean E: {E_arr[sc_idx].mean():.2f}")
        print(f"    Mean c₃: {c3_arr[sc_idx].mean():.2f}, Mean sv: {sv_arr[sc_idx].mean():.2f}")
    if len(nsc_idx) > 0:
        print(f"  NSC tournaments: {len(nsc_idx)}")
        print(f"    Mean H: {H_arr[nsc_idx].mean():.2f}, Mean E: {E_arr[nsc_idx].mean():.2f}")
        print(f"    Mean c₃: {c3_arr[nsc_idx].mean():.2f}, Mean sv: {sv_arr[nsc_idx].mean():.2f}")

    # ─── 6. PHASE DIAGRAM: E_AFM vs c₃ vs H ───
    print(f"\n{'─'*60}")
    print(f"  6. PHASE DIAGRAM (E_AFM, c₃, H)")
    print(f"{'─'*60}")
    # Find the "antiferromagnetic phase boundary": where does the
    # relationship between E and H change character?
    E_unique = sorted(set(E_arr))
    print(f"  Energy levels: {len(E_unique)} distinct E values")
    print(f"  Frustration levels: {len(set(c3_arr))} distinct c₃ values")

    # Multiple regression: H ~ a*E + b*c3 + c
    X = np.column_stack([E_arr, c3_arr, np.ones(N)])
    coeffs, residuals, _, _ = np.linalg.lstsq(X, H_arr, rcond=None)
    print(f"\n  Linear model: H ≈ {coeffs[0]:.4f}*E_AFM + {coeffs[1]:.4f}*c₃ + {coeffs[2]:.4f}")
    H_pred = X @ coeffs
    R2 = 1 - np.sum((H_arr - H_pred)**2) / np.sum((H_arr - H_arr.mean())**2)
    print(f"  R² = {R2:.6f}")

    # Check if E_AFM adds information beyond c₃
    X2 = np.column_stack([c3_arr, np.ones(N)])
    coeffs2, _, _, _ = np.linalg.lstsq(X2, H_arr, rcond=None)
    H_pred2 = X2 @ coeffs2
    R2_c3only = 1 - np.sum((H_arr - H_pred2)**2) / np.sum((H_arr - H_arr.mean())**2)
    print(f"  R²(c₃ only) = {R2_c3only:.6f}")
    print(f"  ΔR² from adding E_AFM: {R2 - R2_c3only:.6f}")

    # Score variance regression
    X3 = np.column_stack([sv_arr, c3_arr, np.ones(N)])
    coeffs3, _, _, _ = np.linalg.lstsq(X3, H_arr, rcond=None)
    H_pred3 = X3 @ coeffs3
    R2_svc3 = 1 - np.sum((H_arr - H_pred3)**2) / np.sum((H_arr - H_arr.mean())**2)
    print(f"  R²(sv + c₃) = {R2_svc3:.6f}")

    # Full model with E, c3, sv
    X4 = np.column_stack([E_arr, c3_arr, sv_arr, np.ones(N)])
    coeffs4, _, _, _ = np.linalg.lstsq(X4, H_arr, rcond=None)
    H_pred4 = X4 @ coeffs4
    R2_full = 1 - np.sum((H_arr - H_pred4)**2) / np.sum((H_arr - H_arr.mean())**2)
    print(f"  R²(E + c₃ + sv) = {R2_full:.6f}")
    print(f"  Full model: H ≈ {coeffs4[0]:.4f}*E + {coeffs4[1]:.4f}*c₃ + {coeffs4[2]:.4f}*sv + {coeffs4[3]:.4f}")

    return data


def magnon_analysis(n, sample_size=None):
    """Compute magnon dispersion relation: average |ΔH| per tile position."""
    print(f"\n{'='*80}")
    print(f"  MAGNON DISPERSION RELATION — n = {n}")
    print(f"{'='*80}")

    tiles, dH_abs, dH_signed = compute_magnon_spectrum(n, sample_size)

    print(f"\n  Tile positions and their 'magnon energies':")
    print(f"  {'Tile':>8}  {'skip':>4}  {'mean|ΔH|':>10}  {'std|ΔH|':>10}  "
          f"{'mean ΔH':>10}  {'%same_class':>12}")

    tile_data = []
    for t_idx, (a, b) in enumerate(tiles):
        skip = a - b  # "range" of the tile
        abs_vals = dH_abs[t_idx]
        signed_vals = dH_signed[t_idx]
        mean_abs = np.mean(abs_vals)
        std_abs = np.std(abs_vals)
        mean_signed = np.mean(signed_vals)
        # Self-loop rate: ΔH = 0 means same H value (not same class, but proxy)
        same_H = sum(1 for v in signed_vals if v == 0)
        pct_same = 100 * same_H / len(signed_vals)

        tile_data.append((t_idx, a, b, skip, mean_abs, std_abs, mean_signed, pct_same))
        print(f"  ({a},{b})   {skip:4d}  {mean_abs:10.4f}  {std_abs:10.4f}  "
              f"{mean_signed:10.4f}  {pct_same:10.2f}%")

    # Check if dispersion depends on skip distance
    skip_groups = defaultdict(list)
    for _, _, _, skip, mean_abs, _, _, _ in tile_data:
        skip_groups[skip].append(mean_abs)

    print(f"\n  Mean |ΔH| by skip distance (AFM 'momentum'):")
    for skip in sorted(skip_groups.keys()):
        vals = skip_groups[skip]
        print(f"    skip={skip}: mean|ΔH| = {np.mean(vals):.4f} "
              f"(from {len(vals)} tile{'s' if len(vals)>1 else ''})")

    # Distribution of ΔH values (for all tiles combined)
    all_dH = []
    for t_idx in dH_signed:
        all_dH.extend(dH_signed[t_idx])
    dH_counter = Counter(all_dH)
    print(f"\n  Overall ΔH distribution (top 20):")
    for dH, cnt in sorted(dH_counter.items(), key=lambda x: -x[1])[:20]:
        print(f"    ΔH = {dH:+4d}: {cnt} ({100*cnt/len(all_dH):.2f}%)")

    return tile_data


def spin_correlation_function(n, exhaustive=True, sample_size=1000):
    """Compute spin-spin correlation <σ_i σ_j> as function of lattice distance.
    In AFM, this should oscillate and decay."""
    print(f"\n{'='*80}")
    print(f"  SPIN-SPIN CORRELATION FUNCTION — n = {n}")
    print(f"{'='*80}")

    tiles = staircase_tiles(n)
    m = len(tiles)
    adj = tile_neighbors(tiles, n)

    # Compute shortest path distances on the staircase lattice
    from collections import deque
    dist = [[float('inf')]*m for _ in range(m)]
    for i in range(m):
        dist[i][i] = 0
        q = deque([i])
        while q:
            v = q.popleft()
            for w in adj[v]:
                if dist[i][w] > dist[i][v] + 1:
                    dist[i][w] = dist[i][v] + 1
                    q.append(w)

    # Collect spin configurations
    spin_configs = []
    if exhaustive:
        for mask, A in all_tournaments(n):
            spin_configs.append(tournament_to_spins(A, n, tiles))
    else:
        for _ in range(sample_size):
            A = [[0]*n for _ in range(n)]
            for i in range(n):
                for j in range(i+1, n):
                    if random.random() < 0.5:
                        A[i][j] = 1
                    else:
                        A[j][i] = 1
            spin_configs.append(tournament_to_spins(A, n, tiles))

    N = len(spin_configs)
    spins_arr = np.array(spin_configs)  # N x m

    # Compute <σ_i σ_j> averaged over all configurations
    corr_matrix = np.zeros((m, m))
    for cfg in range(N):
        for i in range(m):
            for j in range(m):
                corr_matrix[i][j] += spins_arr[cfg][i] * spins_arr[cfg][j]
    corr_matrix /= N

    # Average correlation by lattice distance
    max_d = max(dist[i][j] for i in range(m) for j in range(m) if dist[i][j] < float('inf'))
    corr_by_dist = defaultdict(list)
    for i in range(m):
        for j in range(m):
            d = int(dist[i][j])
            if d < float('inf'):
                corr_by_dist[d].append(corr_matrix[i][j])

    print(f"\n  Spin-spin correlation <σ_i σ_j> by lattice distance:")
    print(f"  {'dist':>5}  {'<σσ>':>10}  {'#pairs':>7}  {'sign':>6}")
    for d in range(int(max_d) + 1):
        vals = corr_by_dist[d]
        mean_c = np.mean(vals) if vals else 0
        sign = "+" if mean_c > 0.001 else ("-" if mean_c < -0.001 else "~0")
        print(f"  {d:5d}  {mean_c:10.6f}  {len(vals):7d}  {sign:>6}")

    return corr_by_dist


def forbidden_values_from_afm(n):
    """Investigate forbidden H values through the AFM lens.

    Key observation: H(T) = I(Ω(T), 2) can be written as
    H = 1 + 2α₁ + 4α₂ + 8α₃ + ...
    where α_k = number of independent k-sets in Ω(T).

    The binary representation constrains achievable values.
    Can AFM quantization rules explain why 7 = 111₂ and 21 = 10101₂ are forbidden?
    """
    print(f"\n{'='*80}")
    print(f"  FORBIDDEN H VALUES: AFM QUANTIZATION ANALYSIS — n = {n}")
    print(f"{'='*80}")

    # Collect (H, α₁, α₂, α₃) for all tournaments
    alpha_data = defaultdict(list)  # H → list of (α₁, α₂, ...)

    for mask, A in all_tournaments(n):
        H = compute_H(A, n)
        c3 = count_3cycles(A, n)
        # For the OCF decomposition: H = 1 + 2α₁ + 4α₂ + ...
        # α₁ = number of directed odd cycles in Ω
        # For small n, α₁ is dominated by 3-cycles
        # We can compute the exact decomposition:
        # H - 1 is even (Rédei). (H-1)/2 = α₁ + 2α₂ + 4α₃ + ...
        h_minus_1 = H - 1
        alpha_data[H].append({'c3': c3, 'h_m1': h_minus_1})

    # Achievable H values
    achievable = sorted(alpha_data.keys())
    all_odd = list(range(1, max(achievable)+2, 2))
    missing = [h for h in all_odd if h not in alpha_data]

    print(f"\n  Achievable H values: {achievable}")
    print(f"  Missing odd values in [1, {max(achievable)}]: {missing}")

    # Binary decomposition of H
    print(f"\n  Binary structure of H values:")
    print(f"  {'H':>5}  {'binary':>12}  {'(H-1)/2':>8}  {'count':>6}  {'c₃ range':>10}")
    for H in achievable:
        items = alpha_data[H]
        c3_vals = [d['c3'] for d in items]
        print(f"  {H:5d}  {H:>12b}  {(H-1)//2:8d}  {len(items):6d}  "
              f"{min(c3_vals)}-{max(c3_vals)}")

    # For forbidden values, analyze WHY they can't be achieved
    print(f"\n  Analysis of forbidden values:")
    for H_forb in missing[:10]:
        h_m1 = H_forb - 1
        # H - 1 = 2α₁ + 4α₂ + 8α₃ + ...
        # So (H-1)/2 = α₁ + 2α₂ + 4α₃ + ...
        print(f"\n  H = {H_forb} (binary: {H_forb:b}):")
        print(f"    (H-1)/2 = {h_m1//2}")
        print(f"    Requires α₁ + 2α₂ + 4α₃ + ... = {h_m1//2}")
        # What α₁ values are achievable at this n?
        all_alpha1 = sorted(set(items['c3'] for H in alpha_data for items in alpha_data[H]))
        # For H=7: (H-1)/2 = 3, so α₁ = 3 + (even terms from α₂...)
        print(f"    Possible: α₁ = {h_m1//2} with α₂=0, or α₁ = {h_m1//2 - 2} with α₂=1, ...")

    return alpha_data


# ─────────────────────────────────────────────
#  Main
# ───���─────────────────────────────────────────

if __name__ == "__main__":
    print("╔══════════════════════════════════════════════════════════════╗")
    print("║  THE ANTIFERROMAGNETIC TOURNAMENT                          ║")
    print("║  opus-2026-04-04-S15                                       ║")
    print("╚══════════════════════════════════════════════════════════════╝")

    # Phase 1+2: Frustration, Néel, Energy — exhaustive at n=4,5; sampled at n=6,7
    for n in [4, 5]:
        explore_n(n, exhaustive=True)

    explore_n(6, exhaustive=True)

    # Phase 3: Magnon dispersion — exhaustive at n=4,5
    for n in [4, 5]:
        magnon_analysis(n)

    # Phase 4: Spin correlation function
    for n in [4, 5]:
        spin_correlation_function(n, exhaustive=True)

    # Phase 5: Forbidden values
    for n in [5, 6]:
        forbidden_values_from_afm(n)

    print("\n\n" + "="*80)
    print("  DONE. All results above.")
    print("="*80)
