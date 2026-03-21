#!/usr/bin/env python3
"""
spectral_kurtosis_deep_s97.py — opus-2026-03-21-S97

DEEP INVESTIGATION: SPECTRAL KURTOSIS, THE TWO-PARAMETER H FORMULA,
AND THE CARTAN BRIDGE

Building on:
- kind-pasteur S13: Spectral flatness principle (min tr(S^4) <=> max H among regulars)
- kind-pasteur S14: Two-parameter formula H ~ f(S_2, tr(S^4)), 96% via S_2
- opus S95c: Cartan commutator ||[A,S]||^2 = n*S_2/2

Key questions for tonight:
1. Does the two-parameter formula extend to ALL tournaments (not just regular)?
2. What is the EXACT relationship between tr(S^4) and c3, alpha_1, alpha_2?
3. Can we PROVE that min kurtosis => max H?
4. Higher moments: tr(S^6), tr(S^8) — do they refine H prediction?
5. Unified Cartan picture: H ~ g(coupling, kurtosis)

Author: opus-2026-03-21-S97
"""

import sys
import numpy as np
from itertools import combinations, permutations
from collections import defaultdict
from scipy import stats
sys.stdout.reconfigure(line_buffering=True)

def all_tournaments(n):
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(pairs)
    for bits in range(2**m):
        A = np.zeros((n, n), dtype=int)
        for k, (i, j) in enumerate(pairs):
            if (bits >> k) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        yield A

def count_hp(A, n):
    dp = defaultdict(int)
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[(mask, v)] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    dp[(mask | (1 << w), w)] += dp[(mask, v)]
    full = (1 << n) - 1
    return sum(dp[(full, v)] for v in range(n))

def count_3cycles(A, n):
    c3 = 0
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]:
            c3 += 1
        if A[i][k] and A[k][j] and A[j][i]:
            c3 += 1
    return c3

def count_all_odd_cycles(A, n):
    """Count all directed odd cycles (alpha_1 in OCF)."""
    total = 0
    for length in range(3, n + 1, 2):
        for verts in combinations(range(n), length):
            first = verts[0]
            rest = list(verts[1:])
            for perm in permutations(rest):
                path = [first] + list(perm)
                ok = True
                for k in range(length):
                    if not A[path[k]][path[(k + 1) % length]]:
                        ok = False
                        break
                if ok:
                    total += 1
    return total

def score_sequence(A, n):
    return tuple(sorted(A.sum(axis=1)))

def spectral_data(A, n):
    """Compute all spectral invariants of the skew-adjacency B = A - A^T."""
    B = (A - A.T).astype(float)
    eigs_B = np.linalg.eigvals(B)
    # B is skew => eigenvalues are purely imaginary (±iθ_k) or 0
    imag_parts = np.sort(np.abs(eigs_B.imag))[::-1]
    # Remove near-zero (the 0 eigenvalue for odd n)
    thetas = imag_parts[imag_parts > 1e-10]

    # Spectral moments: tr(S^{2k}) = (-1)^k * Σ θ_j^{2k}
    # Actually tr(B^{2k}) = Σ (iθ)^{2k} = (-1)^k Σ θ^{2k}
    tr_B2 = round(np.real(np.trace(B @ B)))  # = -Σθ² = -n(n-1) always
    B4 = B @ B @ B @ B
    tr_B4 = round(np.real(np.trace(B4)))     # = Σθ⁴
    B6 = B4 @ B @ B
    tr_B6 = round(np.real(np.trace(B6)))     # = -Σθ⁶
    B8 = B6 @ B @ B
    tr_B8 = round(np.real(np.trace(B8)))     # = Σθ⁸

    # Score sequence and Landau irregularity
    scores = A.sum(axis=1)
    S2 = np.sum((scores - (n - 1) / 2) ** 2)

    # Spectral kurtosis: a normalized measure
    # κ_spec = tr(B^4) / tr(B^2)^2 * n (normalized by dimension)
    kurtosis = tr_B4 / (tr_B2 ** 2) * n if tr_B2 != 0 else 0

    return {
        'thetas': thetas,
        'tr_B2': tr_B2,
        'tr_B4': tr_B4,
        'tr_B6': tr_B6,
        'tr_B8': tr_B8,
        'S2': S2,
        'scores': tuple(sorted(scores)),
        'kurtosis': kurtosis,
    }

print("=" * 72)
print("  SPECTRAL KURTOSIS DEEP INVESTIGATION")
print("  opus-2026-03-21-S97")
print("=" * 72)

# ========================================================================
# PART 1: Two-parameter formula for ALL tournaments
# ========================================================================
print("\n" + "=" * 72)
print("  PART 1: TWO-PARAMETER H FORMULA — ALL TOURNAMENTS")
print("=" * 72)

for n in [5, 6, 7]:
    print(f"\n--- n = {n} ---")
    data = []
    count = 0

    if n <= 6:
        for A in all_tournaments(n):
            H = count_hp(A, n)
            c3 = count_3cycles(A, n)
            sd = spectral_data(A, n)
            data.append({
                'H': H, 'c3': c3,
                'S2': sd['S2'],
                'tr_B4': sd['tr_B4'],
                'tr_B6': sd['tr_B6'],
                'tr_B8': sd['tr_B8'],
                'kurtosis': sd['kurtosis'],
                'scores': sd['scores'],
            })
    else:
        # Sample at n=7
        import random
        random.seed(42)
        for _ in range(5000):
            A = np.zeros((n, n), dtype=int)
            for i in range(n):
                for j in range(i + 1, n):
                    if random.random() < 0.5:
                        A[i][j] = 1
                    else:
                        A[j][i] = 1
            H = count_hp(A, n)
            c3 = count_3cycles(A, n)
            sd = spectral_data(A, n)
            data.append({
                'H': H, 'c3': c3,
                'S2': sd['S2'],
                'tr_B4': sd['tr_B4'],
                'tr_B6': sd['tr_B6'],
                'tr_B8': sd['tr_B8'],
                'kurtosis': sd['kurtosis'],
                'scores': sd['scores'],
            })

    # Correlations
    Hs = np.array([d['H'] for d in data])
    S2s = np.array([d['S2'] for d in data])
    trB4s = np.array([d['tr_B4'] for d in data])
    trB6s = np.array([d['tr_B6'] for d in data])
    kurtoses = np.array([d['kurtosis'] for d in data])
    c3s = np.array([d['c3'] for d in data])

    corr_S2 = np.corrcoef(Hs, -S2s)[0, 1]
    corr_trB4 = np.corrcoef(Hs, -trB4s)[0, 1]
    corr_kurt = np.corrcoef(Hs, -kurtoses)[0, 1]
    corr_c3 = np.corrcoef(Hs, c3s)[0, 1]

    # Spearman rank correlations
    spear_S2 = stats.spearmanr(Hs, -S2s)[0]
    spear_trB4 = stats.spearmanr(Hs, -trB4s)[0]
    spear_kurt = stats.spearmanr(Hs, -kurtoses)[0]

    print(f"  {len(data)} tournaments analyzed")
    print(f"\n  Pearson correlations with H:")
    print(f"    corr(H, -S₂)     = {corr_S2:+.4f}  (Landau irregularity)")
    print(f"    corr(H, -tr(B⁴)) = {corr_trB4:+.4f}  (spectral kurtosis)")
    print(f"    corr(H, -κ_spec) = {corr_kurt:+.4f}  (normalized kurtosis)")
    print(f"    corr(H, c₃)      = {corr_c3:+.4f}  (3-cycle count)")

    print(f"\n  Spearman rank correlations:")
    print(f"    ρ(H, -S₂)     = {spear_S2:+.4f}")
    print(f"    ρ(H, -tr(B⁴)) = {spear_trB4:+.4f}")
    print(f"    ρ(H, -κ_spec) = {spear_kurt:+.4f}")

    # Two-parameter regression: H ~ a*S2 + b*tr(B^4) + c
    X = np.column_stack([S2s, trB4s, np.ones(len(data))])
    beta, residuals, rank, sv = np.linalg.lstsq(X, Hs, rcond=None)
    H_pred = X @ beta
    R2 = 1 - np.sum((Hs - H_pred) ** 2) / np.sum((Hs - np.mean(Hs)) ** 2)

    print(f"\n  Two-parameter regression: H ~ a·S₂ + b·tr(B⁴) + c")
    print(f"    a = {beta[0]:.6f}")
    print(f"    b = {beta[1]:.6f}")
    print(f"    c = {beta[2]:.4f}")
    print(f"    R² = {R2:.6f}")

    # Three-parameter: add tr(B^6)
    X3 = np.column_stack([S2s, trB4s, trB6s, np.ones(len(data))])
    beta3, _, _, _ = np.linalg.lstsq(X3, Hs, rcond=None)
    H_pred3 = X3 @ beta3
    R2_3 = 1 - np.sum((Hs - H_pred3) ** 2) / np.sum((Hs - np.mean(Hs)) ** 2)
    print(f"\n  Three-parameter: H ~ a·S₂ + b·tr(B⁴) + c·tr(B⁶) + d")
    print(f"    R² = {R2_3:.6f}  (improvement: {R2_3 - R2:.6f})")

    # S2-only regression
    X1 = np.column_stack([S2s, np.ones(len(data))])
    beta1, _, _, _ = np.linalg.lstsq(X1, Hs, rcond=None)
    H_pred1 = X1 @ beta1
    R2_1 = 1 - np.sum((Hs - H_pred1) ** 2) / np.sum((Hs - np.mean(Hs)) ** 2)
    print(f"\n  S₂-only: H ~ a·S₂ + b")
    print(f"    R² = {R2_1:.6f}")
    print(f"    Adding tr(B⁴) improves R² by {R2 - R2_1:.6f}")

    # c3-only for comparison
    Xc3 = np.column_stack([c3s, np.ones(len(data))])
    betac3, _, _, _ = np.linalg.lstsq(Xc3, Hs, rcond=None)
    H_predc3 = Xc3 @ betac3
    R2_c3 = 1 - np.sum((Hs - H_predc3) ** 2) / np.sum((Hs - np.mean(Hs)) ** 2)
    print(f"\n  c₃-only: R² = {R2_c3:.6f}")

    # Conditional correlation within score classes
    if n <= 6:
        score_groups = defaultdict(list)
        for d in data:
            score_groups[d['scores']].append(d)

        cond_corrs = []
        for sc, group in score_groups.items():
            if len(group) < 5:
                continue
            Hs_g = [d['H'] for d in group]
            trB4_g = [d['tr_B4'] for d in group]
            if np.std(Hs_g) > 0 and np.std(trB4_g) > 0:
                cc = np.corrcoef(Hs_g, [-t for t in trB4_g])[0, 1]
                cond_corrs.append(cc)

        if cond_corrs:
            print(f"\n  Conditional corr(H, -tr(B⁴) | score seq fixed):")
            print(f"    Mean: {np.mean(cond_corrs):.4f}")
            print(f"    Median: {np.median(cond_corrs):.4f}")
            print(f"    Range: [{min(cond_corrs):.4f}, {max(cond_corrs):.4f}]")

# ========================================================================
# PART 2: EXACT RELATIONSHIP — tr(B^4), c3, and H
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 2: WHAT DETERMINES tr(B⁴)?")
print("=" * 72)
print("""
  tr(B⁴) = ||B²||_F² = Σ_{i,j} (B²[i,j])²

  B²[i,j] for i≠j counts (with sign) the number of 2-paths between i,j.
  Specifically: B²[i,j] = #{k: both agree on k} - #{k: they disagree on k}
  where "agree" means both beat k or both lose to k.

  For the DIAGONAL: B²[i,i] = -(n-1) always.

  So tr(B⁴) = n(n-1)² + Σ_{i≠j} (B²[i,j])².

  The off-diagonal part Σ_{i≠j}(B²[i,j])² depends on the tournament.

  QUESTION: Is tr(B⁴) determined by (c3, score_sequence)?
""")

for n in [5, 6]:
    print(f"\n--- n = {n}: what determines tr(B⁴)? ---")
    by_key = defaultdict(set)

    for A in all_tournaments(n):
        c3 = count_3cycles(A, n)
        sc = score_sequence(A, n)
        sd = spectral_data(A, n)
        by_key[(c3, sc)].add(sd['tr_B4'])

    all_unique = True
    for key, vals in sorted(by_key.items()):
        if len(vals) > 1:
            all_unique = False
            c3, sc = key
            if n <= 5:
                print(f"  c3={c3}, scores={sc}: tr(B⁴) ∈ {sorted(vals)}")

    if all_unique:
        print(f"  tr(B⁴) is determined by (c3, score_sequence): YES")
    else:
        print(f"  tr(B⁴) is determined by (c3, score_sequence): NO")

    # Check if determined by score sequence alone
    by_score = defaultdict(set)
    for A in all_tournaments(n):
        sc = score_sequence(A, n)
        sd = spectral_data(A, n)
        by_score[sc].add(sd['tr_B4'])

    score_unique = all(len(v) == 1 for v in by_score.values())
    print(f"  tr(B⁴) determined by score_sequence alone: {score_unique}")

    # Check if adding c3 helps
    by_c3 = defaultdict(set)
    for A in all_tournaments(n):
        c3 = count_3cycles(A, n)
        sd = spectral_data(A, n)
        by_c3[c3].add(sd['tr_B4'])

    c3_unique = all(len(v) == 1 for v in by_c3.values())
    print(f"  tr(B⁴) determined by c3 alone: {c3_unique}")

# ========================================================================
# PART 3: EXACT FORMULA for tr(B⁴) in terms of tournament invariants
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 3: DERIVING tr(B⁴) FROM TOURNAMENT INVARIANTS")
print("=" * 72)
print("""
  tr(B⁴) = Σ_i (B²)²[i,i] + Σ_{i≠j} (B²[i,j])²
  But (B²)[i,i] = -(n-1) for all i, so diagonal contributes n(n-1)².

  Off-diagonal: B²[i,j] = Σ_{k≠i,j} B[i,k]·B[k,j] for i≠j.

  B[i,k]·B[k,j] =
    +1 if (i→k AND k→j) OR (k→i AND j→k) — "aligned through k"
    -1 if (i→k AND j→k) OR (k→i AND k→j) — "anti-aligned through k"

  So B²[i,j] = #{aligned k} - #{anti-aligned k}

  Define a(i,j) = #{aligned} and d(i,j) = #{anti-aligned}.
  Then B²[i,j] = a(i,j) - d(i,j), and a(i,j) + d(i,j) = n-2.
  So B²[i,j] = 2·a(i,j) - (n-2).

  The off-diagonal sum is:
  Σ_{i≠j}(B²[i,j])² = Σ_{i≠j}(2a(i,j)-(n-2))²

  Now a(i,j) = #{k≠i,j: i→k→j or j→k→i}
             = w(i,j) + w(j,i) [from S95c notation]

  And we showed in S95c that w(i,j)+w(j,i) is NOT determined by
  score sequence alone.

  But Σ_{i≠j} a(i,j)² might have a cleaner form.

  Let's check: is Σ_{i≠j} a(i,j)² determined by (c3, scores)?
""")

for n in [5]:
    print(f"\n--- n = {n}: sum of a(i,j)² ---")
    by_key = defaultdict(set)

    for A in all_tournaments(n):
        c3 = count_3cycles(A, n)
        H = count_hp(A, n)
        sc = score_sequence(A, n)
        B = (A - A.T).astype(float)
        B2 = B @ B

        # a(i,j) = (B²[i,j] + (n-2)) / 2
        sum_a2 = 0
        for i in range(n):
            for j in range(n):
                if i != j:
                    a_ij = (B2[i, j] + (n - 2)) / 2
                    sum_a2 += a_ij ** 2

        by_key[(H, c3, sc)].add(round(sum_a2))

    print(f"  {'H':>4s} {'c3':>4s} {'score':>20s} {'Σa²':>10s}")
    for key in sorted(by_key.keys()):
        H, c3, sc = key
        vals = sorted(by_key[key])
        print(f"  {H:4d} {c3:4d} {str(sc):>20s} {str(vals):>10s}")

# ========================================================================
# PART 4: SPECTRAL KURTOSIS AND CYCLE PACKING — WHY DOES IT WORK?
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 4: WHY DOES SPECTRAL FLATNESS MAXIMIZE H?")
print("=" * 72)
print("""
  THEOREM ATTEMPT:

  H(T) = I(Ω(T), 2) where Ω is the conflict graph on odd cycles.
  I(Ω, 2) = Σ_k α_k · 2^k where α_k = #{independent sets of size k}.

  To maximize H, we need to maximize Σ α_k · 2^k.
  This is maximized when Ω has the most independent sets.

  Ω has many independent sets when:
  (a) Ω has many vertices (many odd cycles) — this is c3 + c5 + ...
  (b) Ω has few edges (cycles rarely share vertices)

  For REGULAR tournaments:
  - All vertices have equal out-degree (n-1)/2
  - This distributes 3-cycles as uniformly as possible
  - The conflict graph Ω should have the most uniform structure

  For SPECTRALLY FLAT regular tournaments (Paley):
  - Eigenvalues of B are all ±i√n with equal multiplicity
  - This means B² = -nI + J (doubly regular condition)
  - Equivalently: for ALL pairs (i,j), the number of common successors
    is exactly (n-3)/4 (when n ≡ 3 mod 4)

  This DOUBLE REGULARITY means:
  - Every vertex is in exactly the same number of 3-cycles
  - Every PAIR of 3-cycles shares vertices in a predictable pattern
  - The conflict graph Ω is itself highly regular

  A regular conflict graph has the most independent sets among all
  graphs with the same number of vertices and edges (this is related
  to the Kahn-Galvin-Tetali theorem for regular bipartite graphs).

  So: spectral flatness → doubly regular → regular Ω → max indep. sets → max H.

  Let's verify the chain computationally.
""")

if True:  # n=7 regular tournaments
    n = 7
    print(f"\n--- n = {n}: verifying the chain ---")

    # Generate Paley T_7 and non-Paley regulars
    QR7 = set()
    for a in range(1, 7):
        QR7.add((a * a) % 7)

    A_paley = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(n):
            if i != j and ((j - i) % 7) in QR7:
                A_paley[i][j] = 1

    # Count odd cycles and build conflict graph for Paley
    H_paley = count_hp(A_paley, n)
    c3_paley = count_3cycles(A_paley, n)

    # Check doubly regular condition: B^2 = -nI + J?
    B_paley = (A_paley - A_paley.T).astype(float)
    B2_paley = B_paley @ B_paley
    expected = -n * np.eye(n) + np.ones((n, n))
    is_drt = np.allclose(B2_paley, expected)

    print(f"  Paley T_7: H={H_paley}, c3={c3_paley}")
    print(f"  Doubly regular (B² = -nI + J): {is_drt}")

    # Off-diagonal of B² for Paley
    off_diag = B2_paley[0, 1]
    print(f"  B²[0,1] = {off_diag} (should be 1 for DRT)")
    print(f"  All off-diagonal = 1: {np.all(B2_paley[np.eye(n) == 0] == 1)}")

    # For each pair, count aligned vertices
    for i in range(min(3, n)):
        for j in range(i + 1, min(3, n)):
            a_ij = (B2_paley[i, j] + (n - 2)) / 2
            print(f"  a({i},{j}) = {a_ij} (#{'{'}aligned k{'}'})")

    # Now check: how many 3-cycles per vertex?
    per_vertex_c3 = [0] * n
    for i, j, k in combinations(range(n), 3):
        if A_paley[i][j] and A_paley[j][k] and A_paley[k][i]:
            per_vertex_c3[i] += 1
            per_vertex_c3[j] += 1
            per_vertex_c3[k] += 1
        if A_paley[i][k] and A_paley[k][j] and A_paley[j][i]:
            per_vertex_c3[i] += 1
            per_vertex_c3[j] += 1
            per_vertex_c3[k] += 1

    print(f"  3-cycles per vertex: {per_vertex_c3}")
    print(f"  All equal: {len(set(per_vertex_c3)) == 1}")

# ========================================================================
# PART 5: HIGHER SPECTRAL MOMENTS — HOW MANY TO DETERMINE H?
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 5: HOW MANY SPECTRAL MOMENTS DETERMINE H?")
print("=" * 72)

for n in [5, 6]:
    print(f"\n--- n = {n} ---")
    by_moments = defaultdict(set)

    for A in all_tournaments(n):
        H = count_hp(A, n)
        sd = spectral_data(A, n)

        # Use moments up to order 2k
        for k_max in [2, 3, 4, 5]:
            moments = (sd['tr_B4'],) if k_max >= 2 else ()
            if k_max >= 3:
                moments += (sd['tr_B6'],)
            if k_max >= 4:
                moments += (sd['tr_B8'],)
            by_moments[(k_max, moments)].add(H)

    # Check discrimination power at each moment order
    for k_max in [2, 3, 4, 5]:
        relevant = {k: v for k, v in by_moments.items() if k[0] == k_max}
        n_classes = len(relevant)
        n_ambiguous = sum(1 for v in relevant.values() if len(v) > 1)
        max_ambig = max(len(v) for v in relevant.values())
        print(f"  Moments up to tr(B^{2*k_max}): "
              f"{n_classes} classes, "
              f"{n_ambiguous} ambiguous, "
              f"max H-values per class: {max_ambig}")

# ========================================================================
# PART 6: UNIFIED CARTAN-KURTOSIS FORMULA
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 6: UNIFIED CARTAN-KURTOSIS FORMULA")
print("=" * 72)
print("""
  The two-parameter model has:
    Parameter 1: S₂ = ||[A_anti, A_sym]||² × 2/n (Cartan COUPLING)
    Parameter 2: tr(B⁴) (spectral KURTOSIS within tournament sector)

  S₂ measures the BETWEEN-SECTOR entanglement (how much the tournament
  mode couples to the cooperation mode).

  tr(B⁴) measures the WITHIN-SECTOR non-uniformity (how unevenly
  distributed the tournament's eigenvalues are).

  UNIFIED FORMULA:
    H ≈ H_max(n) - α·(between-sector entanglement) - β·(within-sector kurtosis)

  where H_max(n) is the maximum H for tournaments on n vertices.

  Let's fit this and compute the coefficients.
""")

for n in [5, 6, 7]:
    print(f"\n--- n = {n} ---")
    data = []

    if n <= 6:
        for A in all_tournaments(n):
            H = count_hp(A, n)
            sd = spectral_data(A, n)
            data.append((H, sd['S2'], sd['tr_B4']))
    else:
        import random
        random.seed(42)
        for _ in range(5000):
            A = np.zeros((n, n), dtype=int)
            for i in range(n):
                for j in range(i + 1, n):
                    if random.random() < 0.5:
                        A[i][j] = 1
                    else:
                        A[j][i] = 1
            H = count_hp(A, n)
            sd = spectral_data(A, n)
            data.append((H, sd['S2'], sd['tr_B4']))

    Hs = np.array([d[0] for d in data])
    S2s = np.array([d[1] for d in data])
    trB4s = np.array([d[2] for d in data])
    H_max = max(Hs)
    H_mean = np.mean(Hs)

    # Normalized variables
    S2_max = max(S2s)
    trB4_max = max(trB4s)
    trB4_min = min(trB4s)

    # Fit: H = c0 - c1*(S2/S2_max) - c2*((trB4-trB4_min)/(trB4_max-trB4_min))
    if S2_max > 0 and trB4_max > trB4_min:
        x1 = S2s / S2_max
        x2 = (trB4s - trB4_min) / (trB4_max - trB4_min)
        X = np.column_stack([x1, x2, np.ones(len(data))])
        beta, _, _, _ = np.linalg.lstsq(X, Hs, rcond=None)
        H_pred = X @ beta
        R2 = 1 - np.sum((Hs - H_pred) ** 2) / np.sum((Hs - np.mean(Hs)) ** 2)

        print(f"  H ≈ {beta[2]:.2f} - {-beta[0]:.2f}·(S₂/S₂_max) "
              f"- {-beta[1]:.2f}·(trB⁴_norm)")
        print(f"  R² = {R2:.6f}")
        print(f"  H_max = {H_max}, H_mean = {H_mean:.1f}")
        print(f"  Intercept ≈ H_max = {beta[2]:.2f} vs actual {H_max}")

        # Variance decomposition
        # How much variance does S2 alone explain?
        X1 = np.column_stack([x1, np.ones(len(data))])
        b1, _, _, _ = np.linalg.lstsq(X1, Hs, rcond=None)
        R2_S2 = 1 - np.sum((Hs - X1 @ b1) ** 2) / np.sum((Hs - H_mean) ** 2)

        X2 = np.column_stack([x2, np.ones(len(data))])
        b2, _, _, _ = np.linalg.lstsq(X2, Hs, rcond=None)
        R2_trB4 = 1 - np.sum((Hs - X2 @ b2) ** 2) / np.sum((Hs - H_mean) ** 2)

        print(f"\n  Variance decomposition:")
        print(f"    S₂ alone:    R² = {R2_S2:.6f}")
        print(f"    tr(B⁴) alone: R² = {R2_trB4:.6f}")
        print(f"    Combined:    R² = {R2:.6f}")
        print(f"    S₂ share:    {R2_S2/R2*100:.1f}%")

# ========================================================================
# PART 7: SPECTRAL KURTOSIS AS ENGINEERING DIAGNOSTIC
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 7: SPECTRAL KURTOSIS ENGINEERING APPLICATIONS")
print("=" * 72)
print("""
  APPLICATION: Given pairwise comparison data, the spectral kurtosis
  of the resulting tournament tells us about the QUALITY of the ranking.

  Low kurtosis (flat spectrum) → well-balanced competition, many paths
  High kurtosis (peaked spectrum) → dominated competition, few paths

  This applies to:
  - Sports leagues: low kurtosis = competitive league, high = dominated
  - Election preferences: low = many viable orderings, high = clear winner
  - LLM attention: low = balanced reasoning, high = focused attention
  - Consumer choice: low = product parity, high = brand dominance
""")

# Demo: create tournament ranking scenarios and compute kurtosis
n = 7
print(f"\n  Demo: n={n} tournament scenarios")
print(f"  {'Scenario':>25s} {'H':>6s} {'S₂':>6s} {'tr(B⁴)':>8s} {'κ':>8s}")

scenarios = {
    "Transitive (total order)": np.array([[1 if i < j else 0
                                           for j in range(n)] for i in range(n)]),
    "Paley T_7 (balanced)": None,  # will build below
    "Near-transitive": None,
    "Random": None,
}

# Build Paley
QR7 = set()
for a in range(1, 7):
    QR7.add((a * a) % 7)
A_pal = np.zeros((n, n), dtype=int)
for i in range(n):
    for j in range(n):
        if i != j and ((j - i) % 7) in QR7:
            A_pal[i][j] = 1
scenarios["Paley T_7 (balanced)"] = A_pal

# Near-transitive: flip one arc from transitive
A_nt = np.zeros((n, n), dtype=int)
for i in range(n):
    for j in range(i + 1, n):
        A_nt[i][j] = 1
A_nt[0][1] = 0
A_nt[1][0] = 1  # flip 0→1 to 1→0
scenarios["Near-transitive"] = A_nt

# Random
np.random.seed(42)
A_rand = np.zeros((n, n), dtype=int)
for i in range(n):
    for j in range(i + 1, n):
        if np.random.random() < 0.5:
            A_rand[i][j] = 1
        else:
            A_rand[j][i] = 1
scenarios["Random"] = A_rand

for name, A in scenarios.items():
    H = count_hp(A, n)
    sd = spectral_data(A, n)
    print(f"  {name:>25s} {H:6d} {sd['S2']:6.1f} {sd['tr_B4']:8d} "
          f"{sd['kurtosis']:8.4f}")

print("""
  INTERPRETATION:
  - Transitive has MAXIMUM kurtosis (one total order, H=1)
  - Paley has MINIMUM kurtosis (maximally balanced, H=189)
  - Near-transitive perturbs slightly (still mostly ordered)
  - Random falls between

  For a RANKING QUALITY metric:
    Quality(T) = 1 - κ(T)/κ_max
  where κ_max = kurtosis of transitive tournament.

  Quality = 1 means maximally balanced (Paley-like)
  Quality = 0 means total order (transitive)
""")

for name, A in scenarios.items():
    sd = spectral_data(A, n)
    # Compute kurtosis relative to transitive
    kappa_max = spectral_data(scenarios["Transitive (total order)"], n)['kurtosis']
    quality = 1 - sd['kurtosis'] / kappa_max if kappa_max > 0 else 0
    print(f"  {name:>25s}: Quality = {quality:.4f}")

print("\n\nDone. opus-2026-03-21-S97")
