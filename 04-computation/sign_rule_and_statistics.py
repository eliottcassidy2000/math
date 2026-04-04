"""
sign_rule_and_statistics.py — opus-2026-04-04-S6

GOALS:
1. PROVE sign rule (THM-301) from OCF cycle geometry
2. Decompose Q = Q_same + Q_cross + Q_disj, analyze each part's spectrum
3. Statistical: tile correlation matrix, mutual information under H(t)/Z
4. Information-theoretic: entropy of H-distribution, KL from uniform
5. Practical: Q eigenvectors as tournament embedding

KEY PROOF IDEA for sign rule:
  Same-end tiles share vertex v → EVERY cycle through both tiles visits v
  → ANY disjoint pair of cycles through tiles i,j separately share v
  → B(i,j)=0 (no independent pair contribution)
  → c_{ij} = 2*A(i,j) where A counts single-cycle contributions
  → each single cycle through same-end tiles uses them in OPPOSITE directions
     (one enters v, the other exits v) → chi = -1 → A < 0 → c_{ij} < 0.

  Cross-end / disjoint: chi = +1 for single cycles, B(i,j) ≥ 0, so c > 0.
"""
import numpy as np
from itertools import combinations, permutations
from collections import defaultdict
from fractions import Fraction
import sys

sys.path.insert(0, '/Users/e/Documents/GitHub/math/04-computation')
from h_multilinear_complete import compute_all_H, mobius_inversion, make_tiles

def tile_relation(t1, t2):
    x1, y1 = t1; x2, y2 = t2
    if x1 == x2 or y1 == y2: return "same-end"
    if x1 == y2 or y1 == x2: return "cross-end"
    return "disjoint"

def shared_vertex(t1, t2):
    """Return the shared vertex (1-indexed) for same-end or cross-end tiles, or None."""
    x1, y1 = t1; x2, y2 = t2
    shared = set()
    if x1 == x2: shared.add(x1)
    if y1 == y2: shared.add(y1)
    if x1 == y2: shared.add(x1)
    if y1 == x2: shared.add(y1)
    return shared

# ==============================================================
# PART 1: Verify the B(i,j)=0 argument for same-end pairs
# ==============================================================
def verify_B_zero_same_end(n):
    """
    Verify: for same-end tile pairs, the |I|=2 OCF contribution B(i,j) is zero.
    Reason: any two cycles through the two tiles both visit the shared vertex,
    so they cannot be vertex-disjoint (independent).
    """
    print(f"\n{'='*60}")
    print(f"B(i,j)=0 VERIFICATION FOR SAME-END PAIRS n={n}")
    print(f"{'='*60}")

    tiles = make_tiles(n)
    m = len(tiles)
    tile_to_idx = {t: i for i, t in enumerate(tiles)}

    # Enumerate all directed odd cycles (3 and 5)
    from h_ocf_multilinear_v2 import enumerate_directed_cycles, cycle_indicator

    all_cycles = []
    for k in range(3, n+1, 2):
        raw = enumerate_directed_cycles(n, k)
        for cyc in raw:
            chi = cycle_indicator(cyc, n, tiles, tile_to_idx)
            if chi is not None:
                all_cycles.append((cyc, chi))

    # For each same-end pair, check B = 0
    se_tested = 0
    se_B_nonzero = 0

    for i in range(m):
        for j in range(i+1, m):
            t1, t2 = tiles[i], tiles[j]
            if tile_relation(t1, t2) != "same-end":
                continue
            se_tested += 1

            sv = shared_vertex(t1, t2)

            # Find cycles through tile i
            cyc_i = []
            for idx, (cyc, chi) in enumerate(all_cycles):
                for S, c in chi.items():
                    if i in S:
                        cyc_i.append(idx)
                        break

            # Find cycles through tile j
            cyc_j = []
            for idx, (cyc, chi) in enumerate(all_cycles):
                for S, c in chi.items():
                    if j in S:
                        cyc_j.append(idx)
                        break

            # Check: can any pair from cyc_i x cyc_j be vertex-disjoint?
            has_disjoint_pair = False
            for ci in cyc_i:
                for cj in cyc_j:
                    if ci == cj:
                        continue
                    v_ci = set(all_cycles[ci][0])
                    v_cj = set(all_cycles[cj][0])
                    if not (v_ci & v_cj):
                        has_disjoint_pair = True
                        break
                if has_disjoint_pair:
                    break

            if has_disjoint_pair:
                se_B_nonzero += 1
                print(f"  UNEXPECTED: same-end {t1}-{t2} HAS disjoint cycle pair!")

            # Also verify that every cycle through BOTH tiles visits the shared vertex
            for idx, (cyc, chi) in enumerate(all_cycles):
                uses_i = any(i in S for S in chi)
                uses_j = any(j in S for S in chi)
                if uses_i and uses_j:
                    # Does this cycle visit the shared vertex?
                    verts = set(cyc)
                    # Shared vertex in 0-indexed
                    sv_0 = {v-1 for v in sv}
                    if not (verts & sv_0):
                        print(f"  UNEXPECTED: cycle {cyc} through {t1},{t2} "
                              f"doesn't visit shared vertex {sv}!")

    print(f"  Tested {se_tested} same-end pairs")
    print(f"  B>0 violations: {se_B_nonzero}")
    if se_B_nonzero == 0:
        print(f"  CONFIRMED: B(i,j)=0 for all same-end pairs at n={n}")

# ==============================================================
# PART 2: Verify chi = -1 for ALL single cycles through same-end pairs
# ==============================================================
def verify_chi_negative_same_end(n):
    """
    For same-end tiles: every single cycle through both tiles has chi = -1
    (contributes negatively to A(i,j)).
    """
    print(f"\n{'='*60}")
    print(f"CHI = -1 VERIFICATION FOR SAME-END CYCLES n={n}")
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

    violations = 0
    total = 0

    for i in range(m):
        for j in range(i+1, m):
            t1, t2 = tiles[i], tiles[j]
            if tile_relation(t1, t2) != "same-end":
                continue

            # For each single cycle through BOTH tiles i and j,
            # the coefficient of the {i,j} monomial in chi should be -1
            for cyc, chi in all_cycles:
                ij_set = frozenset({i, j})
                if ij_set in chi:
                    coeff = chi[ij_set]
                    total += 1
                    if coeff != -1:
                        violations += 1
                        if violations <= 5:
                            print(f"  UNEXPECTED: {t1}-{t2}: cycle {cyc} "
                                  f"has chi[{{i,j}}] = {coeff} (expected -1)")

    print(f"  Tested {total} (cycle, same-end pair) combinations")
    print(f"  Violations: {violations}")
    if violations == 0:
        print(f"  CONFIRMED: chi = -1 for ALL single cycles through same-end pairs")

# ==============================================================
# PART 3: Verify chi = +1 for cross-end cycles
# ==============================================================
def verify_chi_positive_cross_end(n):
    """For cross-end tiles: chi = +1 for all single cycles through both."""
    print(f"\n{'='*60}")
    print(f"CHI = +1 VERIFICATION FOR CROSS-END CYCLES n={n}")
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

    violations = 0
    total = 0

    for i in range(m):
        for j in range(i+1, m):
            t1, t2 = tiles[i], tiles[j]
            if tile_relation(t1, t2) != "cross-end":
                continue

            for cyc, chi in all_cycles:
                ij_set = frozenset({i, j})
                if ij_set in chi:
                    coeff = chi[ij_set]
                    total += 1
                    if coeff != 1:
                        violations += 1
                        if violations <= 5:
                            print(f"  UNEXPECTED: {t1}-{t2}: cycle {cyc} "
                                  f"has chi[{{i,j}}] = {coeff} (expected +1)")

    print(f"  Tested {total} (cycle, cross-end pair) combinations")
    print(f"  Violations: {violations}")
    if violations == 0:
        print(f"  CONFIRMED: chi = +1 for ALL single cycles through cross-end pairs")

# ==============================================================
# PART 4: Disjoint pair chi analysis
# ==============================================================
def analyze_disjoint_chi(n):
    """For disjoint tiles: what are the chi coefficients?"""
    print(f"\n{'='*60}")
    print(f"DISJOINT PAIR CHI ANALYSIS n={n}")
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

    chi_vals = defaultdict(lambda: defaultdict(int))

    for i in range(m):
        for j in range(i+1, m):
            t1, t2 = tiles[i], tiles[j]
            if tile_relation(t1, t2) != "disjoint":
                continue

            for cyc, chi in all_cycles:
                ij_set = frozenset({i, j})
                if ij_set in chi:
                    coeff = chi[ij_set]
                    chi_vals[(i,j)][coeff] += 1

    print(f"  Chi coefficient distribution for disjoint pairs:")
    all_pos = True
    for (i,j), dist in sorted(chi_vals.items()):
        t1, t2 = tiles[i], tiles[j]
        s1, s2 = t1[0]-t1[1], t2[0]-t2[1]
        total_chi = sum(coeff * count for coeff, count in dist.items())
        print(f"    {t1}-{t2} s=({s1},{s2}): chi values = {dict(dist)}, net = {total_chi}")
        if total_chi < 0:
            all_pos = False

    print(f"  All net chi positive? {all_pos}")

# ==============================================================
# PART 5: Q decomposition into same-end / cross-end / disjoint
# ==============================================================
def decompose_Q(n):
    """Decompose Q = Q_same + Q_cross + Q_disj. Analyze each spectrum."""
    print(f"\n{'='*60}")
    print(f"Q DECOMPOSITION n={n}")
    print(f"{'='*60}")

    tiles = make_tiles(n)
    m = len(tiles)
    H, _, _ = compute_all_H(n)
    coeffs = mobius_inversion(H, m)

    Q = np.zeros((m, m))
    Q_same = np.zeros((m, m))
    Q_cross = np.zeros((m, m))
    Q_disj = np.zeros((m, m))

    for S, c in coeffs.items():
        if bin(S).count('1') != 2:
            continue
        bits = [i for i in range(m) if S & (1 << i)]
        i, j = bits[0], bits[1]
        Q[i][j] = Q[j][i] = c

        rel = tile_relation(tiles[i], tiles[j])
        if rel == "same-end":
            Q_same[i][j] = Q_same[j][i] = c
        elif rel == "cross-end":
            Q_cross[i][j] = Q_cross[j][i] = c
        else:
            Q_disj[i][j] = Q_disj[j][i] = c

    for name, M in [("Q_same", Q_same), ("Q_cross", Q_cross), ("Q_disj", Q_disj), ("Q", Q)]:
        eigvals = sorted(np.linalg.eigvalsh(M), reverse=True)
        pos = sum(1 for e in eigvals if e > 1e-10)
        neg = sum(1 for e in eigvals if e < -1e-10)
        zero = m - pos - neg
        print(f"\n  {name}: signature ({pos}+, {neg}-, {zero} zero)")
        print(f"    Eigenvalues: {[f'{e:.3f}' for e in eigvals[:8]]}{'...' if len(eigvals)>8 else ''}")
        if neg > 0:
            print(f"    Neg eigenvalues: {[f'{e:.3f}' for e in eigvals if e < -1e-10]}")

    # KEY TEST: Is neg(Q_same) = n-2?
    same_neg = sum(1 for e in np.linalg.eigvalsh(Q_same) if e < -1e-10)
    print(f"\n  neg(Q_same) = {same_neg}, n-2 = {n-2}")
    print(f"  MATCH: {same_neg == n-2}")

    # Decompose Q_same further by shared vertex
    print(f"\n  Q_same decomposition by shared vertex:")
    vertex_matrices = {}
    for v in range(1, n+1):
        M_v = np.zeros((m, m))
        count = 0
        for i in range(m):
            for j in range(i+1, m):
                if Q_same[i][j] == 0:
                    continue
                t1, t2 = tiles[i], tiles[j]
                if t1[0] == t2[0] == v or t1[1] == t2[1] == v:
                    M_v[i][j] = M_v[j][i] = Q_same[i][j]
                    count += 1
        if count > 0:
            vertex_matrices[v] = M_v
            eigv = sorted(np.linalg.eigvalsh(M_v))
            neg_v = sum(1 for e in eigv if e < -1e-10)
            print(f"    vertex {v}: {count} pairs, neg eigenvalues: {neg_v}, "
                  f"smallest: {eigv[0]:.3f}")

    # Check: is the sum of vertex matrices = Q_same?
    # (They overlap! A same-end pair can share BOTH a column and row vertex.)
    sum_v = sum(vertex_matrices.values())
    diff = np.max(np.abs(Q_same - sum_v))
    print(f"  Sum of vertex matrices = Q_same? max diff = {diff:.6f}")

# ==============================================================
# PART 6: Statistical analysis — Boltzmann distribution
# ==============================================================
def statistical_analysis(n):
    """Treat P(t) = H(t) / sum(H) as a distribution. Compute statistics."""
    print(f"\n{'='*60}")
    print(f"STATISTICAL / INFORMATION ANALYSIS n={n}")
    print(f"{'='*60}")

    tiles = make_tiles(n)
    m = len(tiles)
    H, _, _ = compute_all_H(n)
    N = 1 << m

    Z = float(sum(H))  # partition function
    P = np.array(H, dtype=float) / Z  # probability distribution

    # 1. Entropy of the distribution
    entropy = -sum(p * np.log2(p) for p in P if p > 0)
    uniform_entropy = m  # log2(2^m) = m
    print(f"\n  Distribution H(t)/Z:")
    print(f"    Z = sum(H) = {int(Z)}")
    print(f"    Entropy = {entropy:.4f} bits (uniform = {uniform_entropy})")
    print(f"    Redundancy = {uniform_entropy - entropy:.4f} bits")
    print(f"    Efficiency = {entropy / uniform_entropy:.4f}")

    # 2. KL divergence from uniform
    kl = sum(p * np.log2(p * N) for p in P if p > 0)
    print(f"    KL(P || uniform) = {kl:.4f} bits")

    # 3. Per-tile marginals
    print(f"\n  Per-tile marginals E[t_k] under H/Z:")
    marginals = []
    for k in range(m):
        # E[t_k] = sum_{t: t_k=1} P(t)
        prob_1 = sum(P[t] for t in range(N) if t & (1 << k))
        marginals.append(prob_1)
        tile = tiles[k]
        skip = tile[0] - tile[1]
        print(f"    tile {tile} skip={skip}: E[t_k] = {prob_1:.4f} "
              f"(bias = {prob_1 - 0.5:.4f})")

    # 4. Pairwise correlations
    print(f"\n  Pairwise correlations Corr(t_i, t_j):")
    corr_matrix = np.zeros((m, m))
    for i in range(m):
        for j in range(i+1, m):
            E_ij = sum(P[t] for t in range(N) if (t & (1<<i)) and (t & (1<<j)))
            cov = E_ij - marginals[i] * marginals[j]
            std_i = np.sqrt(marginals[i] * (1 - marginals[i]))
            std_j = np.sqrt(marginals[j] * (1 - marginals[j]))
            corr = cov / (std_i * std_j) if std_i > 0 and std_j > 0 else 0
            corr_matrix[i][j] = corr_matrix[j][i] = corr

    # Compare correlation sign with quadratic coefficient sign
    coeffs = mobius_inversion(H, m)
    sign_match = 0
    sign_total = 0
    for i in range(m):
        for j in range(i+1, m):
            S = (1 << i) | (1 << j)
            c_ij = coeffs.get(S, 0)
            if c_ij == 0:
                continue
            sign_total += 1
            if (c_ij > 0) == (corr_matrix[i][j] > 0):
                sign_match += 1

    print(f"    Sign(Corr) = Sign(c_ij)? {sign_match}/{sign_total}")

    # Correlation eigenvalues
    corr_eigvals = sorted(np.linalg.eigvalsh(corr_matrix), reverse=True)
    pos_corr = sum(1 for e in corr_eigvals if e > 1e-10)
    neg_corr = sum(1 for e in corr_eigvals if e < -1e-10)
    print(f"    Correlation matrix: ({pos_corr}+, {neg_corr}-)")
    print(f"    Top eigenvalues: {[f'{e:.4f}' for e in corr_eigvals[:5]]}")

    # 5. Mutual information between tiles
    print(f"\n  Pairwise mutual information I(t_i; t_j):")
    mi_matrix = np.zeros((m, m))
    for i in range(m):
        for j in range(i+1, m):
            # I(X;Y) = sum p(x,y) log(p(x,y) / (p(x)*p(y)))
            p00 = sum(P[t] for t in range(N) if not(t&(1<<i)) and not(t&(1<<j)))
            p01 = sum(P[t] for t in range(N) if not(t&(1<<i)) and (t&(1<<j)))
            p10 = sum(P[t] for t in range(N) if (t&(1<<i)) and not(t&(1<<j)))
            p11 = sum(P[t] for t in range(N) if (t&(1<<i)) and (t&(1<<j)))
            pi = marginals[i]
            pj = marginals[j]
            mi = 0
            for px, py, pxy in [(1-pi, 1-pj, p00), (1-pi, pj, p01),
                                 (pi, 1-pj, p10), (pi, pj, p11)]:
                if pxy > 0 and px > 0 and py > 0:
                    mi += pxy * np.log2(pxy / (px * py))
            mi_matrix[i][j] = mi_matrix[j][i] = mi

    # Mutual information vs |c_ij|
    print(f"  MI vs |c_ij| correlation:")
    mi_vals = []
    c_vals = []
    for i in range(m):
        for j in range(i+1, m):
            S = (1 << i) | (1 << j)
            c_ij = coeffs.get(S, 0)
            if c_ij != 0:
                mi_vals.append(mi_matrix[i][j])
                c_vals.append(abs(c_ij))
    if mi_vals:
        corr_mi_c = np.corrcoef(mi_vals, c_vals)[0, 1]
        print(f"    Pearson correlation(MI, |c_ij|) = {corr_mi_c:.4f}")

    # MI by tile relation type
    for rel in ["same-end", "cross-end", "disjoint"]:
        mi_by_rel = []
        for i in range(m):
            for j in range(i+1, m):
                if tile_relation(tiles[i], tiles[j]) == rel:
                    mi_by_rel.append(mi_matrix[i][j])
        if mi_by_rel:
            print(f"    {rel}: mean MI = {np.mean(mi_by_rel):.6f}, "
                  f"max = {max(mi_by_rel):.6f}")

    # 6. Total correlation (multivariate generalization)
    total_corr = sum(-p * np.log2(p) for p in marginals) + \
                 sum(-(1-p) * np.log2(1-p) for p in marginals) - entropy
    print(f"\n  Total correlation C(t_1,...,t_m) = {total_corr:.4f} bits")
    print(f"    (= sum H(t_k) - H(t_1,...,t_m))")

    return corr_matrix, mi_matrix, marginals

# ==============================================================
# PART 7: Correlation matrix vs Q eigenstructure comparison
# ==============================================================
def compare_corr_and_Q(n, corr_matrix):
    """Compare the Boltzmann correlation matrix with the quadratic form Q."""
    print(f"\n{'='*60}")
    print(f"CORRELATION vs Q EIGENSTRUCTURE n={n}")
    print(f"{'='*60}")

    tiles = make_tiles(n)
    m = len(tiles)
    H, _, _ = compute_all_H(n)
    coeffs = mobius_inversion(H, m)

    Q = np.zeros((m, m))
    for S, c in coeffs.items():
        if bin(S).count('1') != 2:
            continue
        bits = [i for i in range(m) if S & (1 << i)]
        Q[bits[0]][bits[1]] = c
        Q[bits[1]][bits[0]] = c

    # Eigenvalue comparison
    Q_eigvals, Q_eigvecs = np.linalg.eigh(Q)
    C_eigvals, C_eigvecs = np.linalg.eigh(corr_matrix)

    # Do Q and Corr share eigenvectors?
    # Compute overlap: |<v_i^Q, v_j^C>|
    overlap = np.abs(Q_eigvecs.T @ C_eigvecs)
    max_overlap_per_Q = np.max(overlap, axis=1)
    print(f"  Max eigenvector overlap per Q eigenvector:")
    for i in range(m):
        print(f"    Q[{i}] (lambda={Q_eigvals[i]:.3f}): "
              f"max |<v,w>| = {max_overlap_per_Q[i]:.4f}")

    # Negative eigenvalue count comparison
    neg_Q = sum(1 for e in Q_eigvals if e < -1e-10)
    neg_C = sum(1 for e in C_eigvals if e < -1e-10)
    print(f"\n  Negative eigenvalue count: Q={neg_Q}, Corr={neg_C}")

    # Subspace alignment: angle between negative eigenspaces
    if neg_Q > 0 and neg_C > 0:
        Q_neg_space = Q_eigvecs[:, Q_eigvals < -1e-10]
        C_neg_space = C_eigvecs[:, C_eigvals < -1e-10]
        # Principal angles
        svd_overlap = np.linalg.svd(Q_neg_space.T @ C_neg_space, compute_uv=False)
        svd_overlap = np.clip(svd_overlap, 0, 1)
        angles = np.arccos(svd_overlap) * 180 / np.pi
        print(f"  Principal angles between negative eigenspaces: "
              f"{[f'{a:.1f}' for a in angles[:min(neg_Q, neg_C)]]}")

# ==============================================================
# MAIN
# ==============================================================
if __name__ == "__main__":
    # Sign rule proof components
    for n in [5, 6]:
        verify_B_zero_same_end(n)
    for n in [5, 6]:
        verify_chi_negative_same_end(n)
    for n in [5, 6]:
        verify_chi_positive_cross_end(n)
    for n in [5, 6]:
        analyze_disjoint_chi(n)

    # Q decomposition
    for n in [5, 6, 7]:
        decompose_Q(n)

    # Statistical analysis
    for n in [5, 6]:
        corr, mi, marg = statistical_analysis(n)
        compare_corr_and_Q(n, corr)
