"""
sign_proof_disjoint_and_embedding.py — opus-2026-04-04-S6

THREAD 1: Complete the sign rule proof for disjoint pairs.
  At n=6, (6,1)-(4,2) has chi values {1:3, -1:1}. The -1 comes from
  a cycle where the tiles are used in opposite directions despite no
  shared vertex. We need: net A(i,j) ≥ 0, and A+2B > 0.

  KEY IDEA: The chi=-1 cycles through disjoint tiles must be counted
  against the chi=+1 cycles. We need to show that +1 always wins.
  The argument: the number of +1 cycles exceeds -1 cycles because
  of the base-path constraint (most cycle orientations require
  same-direction tile usage for disjoint tiles).

THREAD 2: Practical tournament embedding via Q eigenvectors.

THREAD 3: The inverse covariance (precision) matrix and its relation to Q.

THREAD 4: Random tournament statistics — how do coefficients distribute?
"""
import numpy as np
from itertools import combinations
from collections import defaultdict
from fractions import Fraction
import sys

sys.path.insert(0, '/Users/e/Documents/GitHub/math/04-computation')
from h_multilinear_complete import compute_all_H, mobius_inversion, make_tiles, tiling_to_tournament, count_hamiltonian_paths

def tile_relation(t1, t2):
    x1, y1 = t1; x2, y2 = t2
    if x1 == x2 or y1 == y2: return "same-end"
    if x1 == y2 or y1 == x2: return "cross-end"
    return "disjoint"

# ==============================================================
# THREAD 1: Disjoint chi=-1 cycles — what makes them negative?
# ==============================================================
def analyze_disjoint_chi_minus1(n):
    """Deep analysis of the chi=-1 cycles for disjoint pairs."""
    print(f"\n{'='*60}")
    print(f"DISJOINT CHI=-1 CYCLES DEEP ANALYSIS n={n}")
    print(f"{'='*60}")

    tiles = make_tiles(n)
    m = len(tiles)
    tile_to_idx = {t: i for i, t in enumerate(tiles)}

    from h_ocf_multilinear_v2 import enumerate_directed_cycles, cycle_indicator

    all_cycles = []
    for k in range(3, n+1, 2):
        raw = enumerate_directed_cycles(n, k)
        for cyc in raw:
            chi = cycle_indicator(cyc, n, tiles, tile_to_idx)
            if chi is not None:
                all_cycles.append((cyc, chi))

    for i in range(m):
        for j in range(i+1, m):
            t1, t2 = tiles[i], tiles[j]
            if tile_relation(t1, t2) != "disjoint":
                continue

            ij_set = frozenset({i, j})
            pos_cycles = []
            neg_cycles = []
            for cyc, chi in all_cycles:
                if ij_set in chi:
                    coeff = chi[ij_set]
                    if coeff > 0:
                        pos_cycles.append((cyc, coeff))
                    elif coeff < 0:
                        neg_cycles.append((cyc, coeff))

            if neg_cycles:
                print(f"\n  {t1}-{t2} s=({t1[0]-t1[1]},{t2[0]-t2[1]}): "
                      f"+1 cycles: {len(pos_cycles)}, -1 cycles: {len(neg_cycles)}")

                # Analyze the negative cycles
                for cyc, coeff in neg_cycles:
                    # What direction does each tile arc go in this cycle?
                    # Tile i = (x1,y1): arc between vertices x1-1 and y1-1 (0-indexed)
                    x1, y1 = t1; x2, y2 = t2
                    a1, b1 = x1-1, y1-1  # 0-indexed endpoints of tile 1
                    a2, b2 = x2-1, y2-1

                    # Find the arc direction for each tile in the cycle
                    L = len(cyc)
                    dir1 = None
                    dir2 = None
                    for ci in range(L):
                        u = cyc[ci]; v = cyc[(ci+1)%L]
                        if (u == a1 and v == b1): dir1 = "forward"  # high→low
                        elif (u == b1 and v == a1): dir1 = "backward"  # low→high
                        if (u == a2 and v == b2): dir2 = "forward"
                        elif (u == b2 and v == a2): dir2 = "backward"

                    print(f"    Chi=-1 cycle {cyc}: tile1 dir={dir1}, tile2 dir={dir2}")

                # Compare with positive cycles
                if len(pos_cycles) <= 5:
                    for cyc, coeff in pos_cycles:
                        x1, y1 = t1; x2, y2 = t2
                        a1, b1 = x1-1, y1-1; a2, b2 = x2-1, y2-1
                        L = len(cyc)
                        dir1 = dir2 = None
                        for ci in range(L):
                            u = cyc[ci]; v = cyc[(ci+1)%L]
                            if (u==a1 and v==b1): dir1="forward"
                            elif (u==b1 and v==a1): dir1="backward"
                            if (u==a2 and v==b2): dir2="forward"
                            elif (u==b2 and v==a2): dir2="backward"
                        print(f"    Chi=+1 cycle {cyc}: tile1 dir={dir1}, tile2 dir={dir2}")

# ==============================================================
# THREAD 2: Tournament embedding via Q eigenvectors
# ==============================================================
def tournament_embedding(n):
    """Embed tournaments into Q's eigenspace. Practical application."""
    print(f"\n{'='*60}")
    print(f"TOURNAMENT EMBEDDING via Q EIGENVECTORS n={n}")
    print(f"{'='*60}")

    tiles = make_tiles(n)
    m = len(tiles)
    H, _, _ = compute_all_H(n)
    coeffs = mobius_inversion(H, m)
    N = 1 << m

    Q = np.zeros((m, m))
    for S, c in coeffs.items():
        if bin(S).count('1') != 2: continue
        bits = [i for i in range(m) if S & (1 << i)]
        Q[bits[0]][bits[1]] = c
        Q[bits[1]][bits[0]] = c

    eigvals, eigvecs = np.linalg.eigh(Q)

    # Project each tiling into top-k eigenvectors
    # Tiling t → vector x(t) = (t_0-0.5, t_1-0.5, ..., t_{m-1}-0.5) centered
    # Projection: y_k(t) = v_k^T x(t)

    # Sort by |eigenvalue|
    order = sorted(range(m), key=lambda i: abs(eigvals[i]), reverse=True)

    # Compute embeddings
    top_k = min(4, m)
    embeddings = np.zeros((N, top_k))
    for t in range(N):
        x = np.array([(t >> i) & 1 for i in range(m)], dtype=float) - 0.5
        for k_idx in range(top_k):
            embeddings[t, k_idx] = eigvecs[:, order[k_idx]].dot(x)

    # Group by H value
    h_to_tilings = defaultdict(list)
    for t in range(N):
        h_to_tilings[int(H[t])].append(t)

    print(f"  Top {top_k} eigenvalues (by |lambda|): "
          f"{[f'{eigvals[order[k]]:.2f}' for k in range(top_k)]}")
    print(f"  Corresponding eigenvectors project tilings into R^{top_k}")

    # How well does the top-1 projection predict H?
    proj1 = embeddings[:, 0]
    corr_H_proj1 = np.corrcoef(proj1, H)[0, 1]
    print(f"\n  Corr(top-1 projection, H) = {corr_H_proj1:.4f}")

    # Top-2 projection
    if top_k >= 2:
        from numpy.linalg import lstsq
        X = np.column_stack([embeddings[:, :2], np.ones(N)])
        beta, res, _, _ = lstsq(X, H.astype(float), rcond=None)
        pred = X @ beta
        r2 = 1 - np.sum((H - pred)**2) / np.sum((H - H.mean())**2)
        print(f"  R² (linear model with top-2 projections) = {r2:.4f}")

    # Full quadratic model
    X_full = np.column_stack([embeddings, np.ones(N)])
    beta_f, _, _, _ = lstsq(X_full, H.astype(float), rcond=None)
    pred_f = X_full @ beta_f
    r2_full = 1 - np.sum((H - pred_f)**2) / np.sum((H - H.mean())**2)
    print(f"  R² (linear model with top-{top_k} projections) = {r2_full:.4f}")

    # Quadratic model
    features = [embeddings[:, k] for k in range(top_k)]
    for k1 in range(top_k):
        for k2 in range(k1, top_k):
            features.append(embeddings[:, k1] * embeddings[:, k2])
    features.append(np.ones(N))
    X_quad = np.column_stack(features)
    beta_q, _, _, _ = lstsq(X_quad, H.astype(float), rcond=None)
    pred_q = X_quad @ beta_q
    r2_quad = 1 - np.sum((H - pred_q)**2) / np.sum((H - H.mean())**2)
    print(f"  R² (quadratic model with top-{top_k}) = {r2_quad:.4f}")

    # Can we recover iso class from embedding?
    # At n=5, there are 12 iso classes. Do they form clusters?
    if n <= 6:
        print(f"\n  Iso-class separation in embedding space:")
        # Two tilings are in the same iso class iff they produce isomorphic tournaments
        # Quick proxy: same H value (not perfect but correlated)
        h_vals = sorted(set(int(H[t]) for t in range(N)))
        for h in h_vals[:8]:
            tilings = h_to_tilings[h]
            if len(tilings) >= 2:
                coords = embeddings[tilings, :2]
                spread = np.std(coords, axis=0)
                print(f"    H={h}: {len(tilings)} tilings, "
                      f"spread in top-2: ({spread[0]:.3f}, {spread[1]:.3f})")

# ==============================================================
# THREAD 3: Precision matrix (inverse covariance) vs Q
# ==============================================================
def precision_matrix_analysis(n):
    """Compare the precision matrix (inverse of covariance) with Q."""
    print(f"\n{'='*60}")
    print(f"PRECISION MATRIX vs Q n={n}")
    print(f"{'='*60}")

    tiles = make_tiles(n)
    m = len(tiles)
    H, _, _ = compute_all_H(n)
    N = 1 << m
    coeffs = mobius_inversion(H, m)

    Z = float(sum(H))
    P = np.array(H, dtype=float) / Z

    # Compute full covariance matrix
    means = np.zeros(m)
    for t in range(N):
        x = np.array([(t >> k) & 1 for k in range(m)], dtype=float)
        means += P[t] * x

    cov = np.zeros((m, m))
    for t in range(N):
        x = np.array([(t >> k) & 1 for k in range(m)], dtype=float)
        diff = x - means
        cov += P[t] * np.outer(diff, diff)

    # Precision matrix = inverse covariance
    try:
        prec = np.linalg.inv(cov)
    except np.linalg.LinAlgError:
        print("  Covariance matrix singular, using pseudoinverse")
        prec = np.linalg.pinv(cov)

    # Build Q
    Q = np.zeros((m, m))
    for S, c in coeffs.items():
        if bin(S).count('1') != 2: continue
        bits = [i for i in range(m) if S & (1 << i)]
        Q[bits[0]][bits[1]] = c
        Q[bits[1]][bits[0]] = c

    # Compare Q and precision matrix
    # Normalize both to have unit Frobenius norm for comparison
    Q_norm = Q / np.linalg.norm(Q, 'fro')
    prec_offdiag = prec.copy()
    np.fill_diagonal(prec_offdiag, 0)
    prec_norm = prec_offdiag / np.linalg.norm(prec_offdiag, 'fro')

    fro_diff = np.linalg.norm(Q_norm - prec_norm, 'fro')
    cos_sim = np.sum(Q_norm * prec_norm)
    print(f"  Frobenius distance (normalized): {fro_diff:.4f}")
    print(f"  Cosine similarity: {cos_sim:.4f}")

    # Sign agreement of off-diagonal entries
    sign_agree = 0
    sign_total = 0
    for i in range(m):
        for j in range(i+1, m):
            if Q[i][j] != 0:
                sign_total += 1
                if (Q[i][j] > 0) == (prec[i][j] > 0):
                    sign_agree += 1

    print(f"  Sign agreement (Q vs precision): {sign_agree}/{sign_total}")

    # In Ising models, the precision matrix equals the coupling matrix J
    # For H(t) = exp(sum h_i t_i + sum J_ij t_i t_j + ...):
    # J_ij = lambda_ij (log-interaction) from THM-290
    # The precision of the Boltzmann distribution approximates J to first order

    # Compute log-interactions for comparison
    log_Q = np.zeros((m, m))
    for i in range(m):
        for j in range(i+1, m):
            hi = int(H[1 << i])
            hj = int(H[1 << j])
            hij = int(H[(1 << i) | (1 << j)])
            h0 = 1  # H(transitive)
            if hi > 0 and hj > 0 and h0 > 0 and hij > 0:
                log_ij = np.log(hij * h0 / (hi * hj))
                log_Q[i][j] = log_ij
                log_Q[j][i] = log_ij

    # Compare sign of log-interaction with Q
    log_sign_agree = 0
    log_total = 0
    for i in range(m):
        for j in range(i+1, m):
            if Q[i][j] != 0 and abs(log_Q[i][j]) > 1e-15:
                log_total += 1
                if (Q[i][j] > 0) == (log_Q[i][j] > 0):
                    log_sign_agree += 1

    print(f"  Sign agreement (Q vs log-interaction): {log_sign_agree}/{log_total}")
    if log_total > 0 and log_sign_agree < log_total:
        for i in range(m):
            for j in range(i+1, m):
                if Q[i][j] != 0 and abs(log_Q[i][j]) > 1e-15:
                    if (Q[i][j] > 0) != (log_Q[i][j] > 0):
                        print(f"    DISAGREE: {tiles[i]}-{tiles[j]}: "
                              f"Q={Q[i][j]:.1f}, log={log_Q[i][j]:.4f}")

    # Eigenvalue comparison
    Q_eig = sorted(np.linalg.eigvalsh(Q), reverse=True)
    prec_eig = sorted(np.linalg.eigvalsh(prec), reverse=True)
    log_eig = sorted(np.linalg.eigvalsh(log_Q), reverse=True)

    neg_Q = sum(1 for e in Q_eig if e < -1e-10)
    neg_prec = sum(1 for e in prec_eig if e < -1e-10)
    neg_log = sum(1 for e in log_eig if e < -1e-10)
    print(f"\n  Negative eigenvalue count: Q={neg_Q}, Prec={neg_prec}, Log={neg_log}")

# ==============================================================
# THREAD 4: Random tournament statistics
# ==============================================================
def random_tournament_stats(n):
    """Statistics of H, quadratic interactions, etc. over random tournaments."""
    print(f"\n{'='*60}")
    print(f"RANDOM TOURNAMENT STATISTICS n={n}")
    print(f"{'='*60}")

    tiles = make_tiles(n)
    m = len(tiles)
    H, _, _ = compute_all_H(n)
    N = 1 << m

    h_vals = H.astype(float)
    print(f"  H statistics over all {N} tilings:")
    print(f"    mean = {np.mean(h_vals):.3f}")
    print(f"    std  = {np.std(h_vals):.3f}")
    print(f"    min  = {np.min(h_vals)}, max = {np.max(h_vals)}")
    print(f"    median = {np.median(h_vals)}")
    print(f"    skewness = {float(np.mean((h_vals - np.mean(h_vals))**3) / np.std(h_vals)**3):.3f}")
    print(f"    kurtosis = {float(np.mean((h_vals - np.mean(h_vals))**4) / np.std(h_vals)**4 - 3):.3f}")

    # H mod 2: always 1 (Redei)
    # H mod 4: distribution
    mod4 = defaultdict(int)
    mod8 = defaultdict(int)
    for h in H:
        mod4[int(h) % 4] += 1
        mod8[int(h) % 8] += 1
    print(f"    H mod 4: {dict(sorted(mod4.items()))}")
    print(f"    H mod 8: {dict(sorted(mod8.items()))}")

    # Distribution of pairwise "flip effects"
    # For each tile pair (i,j), the "quadratic effect" is
    # delta_{ij} = H(ij) - H(i) - H(j) + H(0) = c_{ij}
    # But THIS IS THE SAME FOR ALL TILINGS (it's a multilinear coefficient).
    # The TILING-DEPENDENT effect is: H(t XOR e_i XOR e_j) - H(t XOR e_i) - H(t XOR e_j) + H(t)
    # which is the "discrete second derivative" at point t.

    print(f"\n  Discrete second derivatives D²H / Dt_i Dt_j:")
    # This equals ∂²H/∂t_i∂t_j for multilinear H, which is just c_{ij}!
    # Wait, for MULTILINEAR H, the second partial ∂²H/∂t_i∂t_j = c_{ij} + higher-order terms
    # involving other variables. So the discrete Hessian at a point t is:
    # H_ij(t) = sum_{S containing i,j} c_S * prod_{k in S\{i,j}} t_k

    # Compute Hessian at different points
    coeffs = mobius_inversion(H, m)
    for point_name, point in [("transitive (0)", 0),
                               ("anti-transitive (all 1)", N-1),
                               ("random", N//3)]:
        print(f"\n    At {point_name}:")
        hess = np.zeros((m, m))
        for S, c in coeffs.items():
            if bin(S).count('1') < 2: continue
            bits = [i for i in range(m) if S & (1 << i)]
            for i_idx, bi in enumerate(bits):
                for j_idx, bj in enumerate(bits):
                    if bi >= bj: continue
                    # Contribution: c * prod_{k in S\{bi,bj}} t_k
                    other_bits = [b for b in bits if b != bi and b != bj]
                    prod_val = 1
                    for b in other_bits:
                        prod_val *= (point >> b) & 1
                    hess[bi][bj] += c * prod_val
                    hess[bj][bi] += c * prod_val

        eigvals = sorted(np.linalg.eigvalsh(hess))
        neg = sum(1 for e in eigvals if e < -1e-10)
        pos = sum(1 for e in eigvals if e > 1e-10)
        print(f"      Hessian signature: ({pos}+, {neg}-, {m-pos-neg} zero)")
        print(f"      Hessian eigenvalues: {[f'{e:.2f}' for e in eigvals[:5]]}...")

# ==============================================================
# THREAD 5: The sign rule via log-interactions
# ==============================================================
def sign_rule_via_log(n):
    """
    THM-290 says all log-interactions lambda_{ij} < 0 (antiferromagnetic).
    The sign of c_{ij} and sign of lambda_{ij} should agree.
    If they always agree, the sign rule follows from antiferromagnetism.
    """
    print(f"\n{'='*60}")
    print(f"SIGN RULE via LOG-INTERACTIONS n={n}")
    print(f"{'='*60}")

    tiles = make_tiles(n)
    m = len(tiles)
    H, _, _ = compute_all_H(n)
    coeffs = mobius_inversion(H, m)

    agree = 0
    disagree = 0
    for i in range(m):
        for j in range(i+1, m):
            S = (1 << i) | (1 << j)
            c_ij = coeffs.get(S, 0)
            if c_ij == 0: continue

            hi = int(H[1 << i])
            hj = int(H[1 << j])
            hij = int(H[S])
            h0 = 1

            log_ij = np.log(hij * h0 / (hi * hj))

            if (c_ij > 0) == (log_ij > 0):
                agree += 1
            else:
                disagree += 1
                rel = tile_relation(tiles[i], tiles[j])
                print(f"  DISAGREE: {tiles[i]}-{tiles[j]} ({rel}): "
                      f"c={c_ij}, log={log_ij:.4f}")

    print(f"  Sign agreement: {agree}, disagreement: {disagree}")
    if disagree == 0:
        print(f"  PERFECT: c_ij and log-interaction ALWAYS have same sign!")
    else:
        print(f"  IMPERFECT: the sign rule does NOT simply follow from log-interaction sign")

if __name__ == "__main__":
    for n in [6]:
        analyze_disjoint_chi_minus1(n)

    for n in [5, 6]:
        sign_rule_via_log(n)

    for n in [5, 6, 7]:
        precision_matrix_analysis(n)

    for n in [5, 6]:
        tournament_embedding(n)

    for n in [5, 6, 7]:
        random_tournament_stats(n)
