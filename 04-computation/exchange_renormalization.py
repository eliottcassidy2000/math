"""
exchange_renormalization.py — opus-2026-04-04-S20

INTEGRATE all recent ideas into a RENORMALIZATION GROUP framework.

THE PIECES:
1. F(n) = 2^n - 4(n-1): total frustration (S19)
2. ΔH ≈ a·Δc₃·H_sub + b·Δc₃ + c·H_sub + d: exchange coupling (S17b)
3. E[Δc₃|s] = s(n-1-s)/2: parabolic law (S17b)
4. H = 1 + 2α₁ + 4α₂: OCF truncation (S11)
5. Fiber bundle V(n) from V(n-1) (S12)

THE GOAL: A RG equation dJ/dn = f(J,n) for the effective coupling.

PRACTICAL: Predict H(T_n) from H(T_{n-1}) and the extension signature σ,
without computing Hamiltonian paths at all.
"""
import numpy as np
from itertools import combinations, permutations
from collections import defaultdict
from math import comb
import sys, time

sys.path.insert(0, '/Users/e/Documents/GitHub/math/04-computation')
from h_multilinear_complete import make_tiles, tiling_to_tournament, count_hamiltonian_paths

# ==============================================================
# 1. VERIFY the exchange coupling at n=5→6→7
# ==============================================================
def exchange_coupling_analysis(n_from, n_to):
    """Measure how H propagates through the fiber bundle."""
    print(f"\n{'='*60}")
    print(f"EXCHANGE COUPLING n={n_from} → n={n_to}")
    print(f"{'='*60}")

    tiles_sub = make_tiles(n_from)
    m_sub = len(tiles_sub)
    tiles_full = make_tiles(n_to)
    m_full = len(tiles_full)

    # Compute H for all sub-tournaments
    H_sub_all, _, _ = compute_all_H_fast(n_from)
    H_full_all, _, _ = compute_all_H_fast(n_to)

    N_sub = 1 << m_sub
    N_full = 1 << m_full

    # For each full tiling, decompose into (sub_tiling, extension σ)
    # Extension tiles are those with x = n_to (0-indexed: x-1 = n_to-1)
    ext_indices = [i for i, (x,y) in enumerate(tiles_full) if x == n_to]
    sub_indices = [i for i, (x,y) in enumerate(tiles_full) if x < n_to]

    # Map sub tiles at n_to to tiles at n_from
    sub_map = {}
    for j, idx in enumerate(sub_indices):
        tile = tiles_full[idx]
        for k, t_sub in enumerate(tiles_sub):
            if t_sub == tile:
                sub_map[j] = k
                break

    data = []
    for t_full in range(N_full):
        # Extract sub-tiling
        sub_bits = 0
        for j, idx in enumerate(sub_indices):
            if t_full & (1 << idx):
                sub_bits |= (1 << sub_map[j])

        # Extension signature
        sigma = tuple((t_full >> idx) & 1 for idx in ext_indices)
        score_new = sum(sigma)

        h_sub = int(H_sub_all[sub_bits])
        h_full = int(H_full_all[t_full])
        delta_h = h_full - h_sub

        # Compute Δc₃ via the tournament
        T_full = tiling_to_tournament(n_to, tiles_full, t_full)
        T_sub_mat = tiling_to_tournament(n_from, tiles_sub, sub_bits)

        c3_full = sum(1 for a in range(n_to) for b in range(a+1,n_to) for c in range(b+1,n_to)
                     if T_full[a][b]+T_full[b][c]+T_full[c][a]==3
                     or T_full[a][c]+T_full[c][b]+T_full[b][a]==3)
        c3_sub = sum(1 for a in range(n_from) for b in range(a+1,n_from) for c in range(b+1,n_from)
                    if T_sub_mat[a][b]+T_sub_mat[b][c]+T_sub_mat[c][a]==3
                    or T_sub_mat[a][c]+T_sub_mat[c][b]+T_sub_mat[b][a]==3)
        delta_c3 = c3_full - c3_sub

        data.append({
            'h_sub': h_sub, 'h_full': h_full, 'delta_h': delta_h,
            'score': score_new, 'delta_c3': delta_c3,
        })

    # Fit: ΔH = a·Δc₃·H_sub + b·Δc₃ + c·H_sub + d
    X = np.array([[d['delta_c3']*d['h_sub'], d['delta_c3'], d['h_sub'], 1] for d in data])
    y = np.array([d['delta_h'] for d in data])

    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    pred = X @ beta
    r2 = 1 - np.sum((y-pred)**2) / np.sum((y-np.mean(y))**2)

    print(f"  Fit: ΔH = {beta[0]:.4f}·Δc₃·H_sub + {beta[1]:.4f}·Δc₃ + "
          f"{beta[2]:.4f}·H_sub + {beta[3]:.4f}")
    print(f"  R² = {r2:.6f}")

    # Also: simple model ΔH = a·Δc₃ + b
    X2 = np.array([[d['delta_c3'], 1] for d in data])
    beta2 = np.linalg.lstsq(X2, y, rcond=None)[0]
    pred2 = X2 @ beta2
    r2_simple = 1 - np.sum((y-pred2)**2) / np.sum((y-np.mean(y))**2)
    print(f"  Simple: ΔH = {beta2[0]:.4f}·Δc₃ + {beta2[1]:.4f} (R²={r2_simple:.4f})")

    # The exchange coupling J_eff = coefficient of Δc₃ in the full model
    # depends on H_sub: J_eff(H_sub) = a·H_sub + b
    j_eff = lambda h: beta[0]*h + beta[1]
    print(f"\n  Exchange coupling: J_eff(H_sub) = {beta[0]:.4f}·H_sub + {beta[1]:.4f}")
    for h in sorted(set(d['h_sub'] for d in data))[:8]:
        print(f"    H_sub={h:4d}: J_eff = {j_eff(h):.3f}")

    # The RENORMALIZATION: J(n) evolves with n
    # From n=4→5: J_eff at various H_sub
    # From n=5→6: J_eff at various H_sub
    # The ratio gives the RG flow

    return beta, r2

def compute_all_H_fast(n):
    """Compute H for all tilings."""
    tiles = make_tiles(n)
    m = len(tiles)
    N = 1 << m
    H = np.zeros(N, dtype=np.int64)
    for t in range(N):
        T = tiling_to_tournament(n, tiles, t)
        H[t] = count_hamiltonian_paths(T, n)
    return H, tiles, m

# ==============================================================
# 2. THE RG EQUATION
# ==============================================================
def rg_equation():
    """Extract the RG flow: how J_eff changes with n."""
    print(f"\n{'='*60}")
    print(f"RENORMALIZATION GROUP FLOW")
    print(f"{'='*60}")

    # Compute exchange coupling at each transition
    transitions = {}
    for n_from, n_to in [(4,5), (5,6)]:
        beta, r2 = exchange_coupling_analysis(n_from, n_to)
        transitions[(n_from, n_to)] = {
            'a': beta[0], 'b': beta[1], 'c': beta[2], 'd': beta[3], 'r2': r2
        }

    print(f"\n  RG FLOW TABLE:")
    print(f"  {'transition':>12s} {'a (J*H)':>10s} {'b (J)':>10s} {'c (H)':>10s} {'d':>10s} {'R²':>8s}")
    for (nf, nt), d in transitions.items():
        print(f"  {nf}→{nt:>3d}      {d['a']:10.4f} {d['b']:10.4f} {d['c']:10.4f} {d['d']:10.4f} {d['r2']:8.4f}")

    # The coupling coefficient a grows: does it follow a pattern?
    a_vals = [d['a'] for d in transitions.values()]
    print(f"\n  Coupling coefficient a: {a_vals}")
    if len(a_vals) >= 2:
        ratio = a_vals[1] / a_vals[0] if a_vals[0] != 0 else 0
        print(f"  Ratio a(5→6)/a(4→5) = {ratio:.3f}")

# ==============================================================
# 3. PRACTICAL: PREDICT H from fiber bundle data
# ==============================================================
def predict_H_from_fiber(n_target):
    """
    Predict H(T_n) from H(T_{n-1}) and Δc₃ using the exchange coupling model.
    This avoids computing Hamiltonian paths at n — only needs c₃ counting.
    """
    print(f"\n{'='*60}")
    print(f"PRACTICAL H PREDICTOR n={n_target}")
    print(f"{'='*60}")

    n_sub = n_target - 1
    tiles_sub = make_tiles(n_sub)
    tiles_full = make_tiles(n_target)
    m_sub = len(tiles_sub)
    m_full = len(tiles_full)

    H_sub_all, _, _ = compute_all_H_fast(n_sub)
    H_full_all, _, _ = compute_all_H_fast(n_target)

    ext_indices = [i for i, (x,y) in enumerate(tiles_full) if x == n_target]
    sub_indices = [i for i, (x,y) in enumerate(tiles_full) if x < n_target]
    sub_map = {}
    for j, idx in enumerate(sub_indices):
        tile = tiles_full[idx]
        for k, t_sub in enumerate(tiles_sub):
            if t_sub == tile:
                sub_map[j] = k
                break

    # First: fit the model on all data
    data_all = []
    for t_full in range(1 << m_full):
        sub_bits = 0
        for j, idx in enumerate(sub_indices):
            if t_full & (1 << idx):
                sub_bits |= (1 << sub_map[j])

        sigma = tuple((t_full >> idx) & 1 for idx in ext_indices)
        score_new = sum(sigma)
        h_sub = int(H_sub_all[sub_bits])
        h_full = int(H_full_all[t_full])

        # Compute Δc₃ efficiently: count 3-cycles involving vertex n_target-1 (0-indexed)
        T_full = tiling_to_tournament(n_target, tiles_full, t_full)
        v = n_target - 1  # new vertex (0-indexed)
        delta_c3 = 0
        for a in range(n_target):
            if a == v: continue
            for b in range(a+1, n_target):
                if b == v: continue
                # 3-cycle {v, a, b}
                if T_full[v][a]+T_full[a][b]+T_full[b][v]==3 or \
                   T_full[v][b]+T_full[b][a]+T_full[a][v]==3:
                    delta_c3 += 1

        data_all.append((h_sub, delta_c3, score_new, h_full))

    # Fit
    X = np.array([[hs*dc, dc, hs, 1] for hs, dc, _, _ in data_all])
    y = np.array([hf for _, _, _, hf in data_all])
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    pred = X @ beta
    r2 = 1 - np.sum((y-pred)**2) / np.sum((y-np.mean(y))**2)

    print(f"  Full model: H_n = {beta[0]:.4f}·Δc₃·H_sub + {beta[1]:.4f}·Δc₃ + "
          f"{beta[2]:.4f}·H_sub + {beta[3]:.4f}")
    print(f"  R² = {r2:.6f}")

    # Error distribution
    errors = np.abs(y - pred)
    print(f"  Mean |error| = {np.mean(errors):.3f}")
    print(f"  Max |error| = {np.max(errors):.3f}")
    print(f"  Fraction |error| ≤ 1: {np.mean(errors <= 1):.4f}")
    print(f"  Fraction |error| ≤ 2: {np.mean(errors <= 2):.4f}")

    # IMPROVED: Use H_sub + 2·(exact Δα₁) + 4·(exact Δα₂)
    # This should be EXACT by THM-304!
    print(f"\n  EXACT PREDICTOR (OCF truncation):")
    exact_matches = 0
    for hs, dc, s, hf in data_all:
        # At n≤8: H = 1 + 2α₁ + 4α₂ exactly
        # So: H_n = H_sub + 2·Δα₁ + 4·Δα₂
        # where Δα₁ = new odd cycles involving vertex n
        # Δα₂ = new vertex-disjoint pairs involving vertex n
        # This IS exact (not an approximation) but requires computing Δα₁ and Δα₂
        # which is basically the same as computing H from OCF...
        pass

    # So the exchange coupling model is a STATISTICAL approximation
    # while the OCF truncation is exact.
    # The exchange coupling's value: it PREDICTS H from CHEAPER features
    # (just Δc₃ and H_sub, not full cycle enumeration).

    # HOW GOOD is just Δc₃ + H_sub?
    X_cheap = np.array([[dc, hs, 1] for hs, dc, _, _ in data_all])
    beta_cheap = np.linalg.lstsq(X_cheap, y, rcond=None)[0]
    pred_cheap = X_cheap @ beta_cheap
    r2_cheap = 1 - np.sum((y-pred_cheap)**2) / np.sum((y-np.mean(y))**2)
    print(f"\n  Cheap model: H_n ≈ {beta_cheap[0]:.4f}·Δc₃ + {beta_cheap[1]:.4f}·H_sub + "
          f"{beta_cheap[2]:.4f} (R²={r2_cheap:.4f})")

    # And: just H_sub?
    X_trivial = np.array([[hs, 1] for hs, _, _, _ in data_all])
    beta_triv = np.linalg.lstsq(X_trivial, y, rcond=None)[0]
    pred_triv = X_trivial @ beta_triv
    r2_triv = 1 - np.sum((y-pred_triv)**2) / np.sum((y-np.mean(y))**2)
    print(f"  Trivial: H_n ≈ {beta_triv[0]:.4f}·H_sub + {beta_triv[1]:.4f} "
          f"(R²={r2_triv:.4f})")

# ==============================================================
# 4. FRUSTRATION PROPAGATION meets EXCHANGE COUPLING
# ==============================================================
def unified_propagation():
    """
    The frustration F(n) = 2^n - 4(n-1) gives the TOTAL negative coupling.
    The exchange coupling J_eff gives how EACH frustration unit translates to ΔH.
    Together: predict mean_H(n) from mean_H(n-1).

    mean_H(n) ≈ mean_H(n-1) + J_eff(mean_H(n-1)) · E[Δc₃]
    where E[Δc₃] ≈ (n-1)/4 · (n-2)/4 ... the parabolic average.
    """
    print(f"\n{'='*60}")
    print(f"UNIFIED PROPAGATION: FRUSTRATION × EXCHANGE")
    print(f"{'='*60}")

    # Known mean_H from S13: mean_H(n) = W(n)/2^{n-1}
    W = {3:8, 4:32, 5:158, 6:928, 7:6350}
    mean_H = {n: W[n] / 2**(n-1) for n in W}

    # The parabolic law gives E[Δc₃]:
    # For uniform random extension at step n: s ~ Binomial(n-1, 1/2)
    # E[Δc₃] = E[s(n-1-s)/2] = (n-1)/4 · (n-2) ... let me compute exactly
    # E[s(n-1-s)] = E[s](n-1) - E[s²] = (n-1)/2 · (n-1) - ((n-1)/4 + (n-1)²/4)
    # = (n-1)²/2 - (n-1)/4 - (n-1)²/4 = (n-1)²/4 - (n-1)/4 = (n-1)(n-2)/4
    # So E[Δc₃] = (n-1)(n-2)/8

    print(f"  {'n':>3s} {'mean_H':>10s} {'ΔH':>8s} {'E[Δc₃]':>8s} {'J=ΔH/E[Δc₃]':>14s} {'F(n)/m':>8s}")
    prev_mH = None
    for n in sorted(mean_H):
        mH = mean_H[n]
        m = comb(n-1, 2)
        F = 2**n - 4*(n-1)
        E_dc3 = (n-1)*(n-2)/8

        if prev_mH is not None:
            delta_mH = mH - prev_mH
            J = delta_mH / E_dc3 if E_dc3 > 0 else 0
            print(f"  {n:3d} {mH:10.4f} {delta_mH:8.4f} {E_dc3:8.3f} {J:14.4f} {F/m:8.2f}")
        else:
            print(f"  {n:3d} {mH:10.4f}      ---      ---           --- {F/m:8.2f}")
        prev_mH = mH

    # The RG flow: J(n) vs n
    print(f"\n  EFFECTIVE EXCHANGE COUPLING GROWTH:")
    J_vals = []
    prev_mH = mean_H[3]
    for n in range(4, 8):
        mH = mean_H[n]
        E_dc3 = (n-1)*(n-2)/8
        J = (mH - prev_mH) / E_dc3 if E_dc3 > 0 else 0
        J_vals.append((n, J))
        prev_mH = mH

    for n, J in J_vals:
        print(f"    n={n}: J = {J:.4f}")

    if len(J_vals) >= 2:
        Js = [j for _, j in J_vals]
        ratios = [Js[i+1]/Js[i] for i in range(len(Js)-1)]
        print(f"  Ratios J(n+1)/J(n): {[f'{r:.3f}' for r in ratios]}")

# ==============================================================
# MAIN
# ==============================================================
if __name__ == "__main__":
    rg_equation()
    for n in [6]:
        predict_H_from_fiber(n)
    unified_propagation()
