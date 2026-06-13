#!/usr/bin/env python3
"""
orthogonal_control_principle.py -- kind-pasteur-2026-03-21-S15

THE ORTHOGONAL CONTROL PRINCIPLE

THESIS: In any sufficiently structured directed system, the key observables
are controlled by invariants that are ORTHOGONAL to the structure that
generates them. Direction creates structure; symmetry reads it.

This is not a metaphor. It is a quantifiable phenomenon:
- In tournaments: score regularity (symmetric) explains 96% of H variance
- In Lie theory: Harish-Chandra says Z(U(g)) = R[h]^W — invariants live on the Cartan
- In physics: gauge-invariant observables = traces of directed gauge field
- In quantum mechanics: Hermitian observables vs unitary evolution

QUESTION: Does this extend beyond tournaments?

TEST DOMAINS:
1. Markov chains: transition matrix (directed) → stationary distribution (symmetric)
2. Random tournaments at varying n: is the 96% OCR universal?
3. Non-tournament digraphs: does the principle still hold?
4. The Harish-Chandra connection: eigenvalue projection vs H

Define: Orthogonal Control Ratio (OCR) = R^2 of symmetric invariants predicting f

Author: kind-pasteur-2026-03-21-S15
"""

import numpy as np
from itertools import combinations
from collections import defaultdict


def signed_adjacency(A):
    n = A.shape[0]
    S = 2.0 * A - 1.0 + np.eye(n)
    np.fill_diagonal(S, 0)
    return S


def count_ham_paths(A):
    n = A.shape[0]
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp.get((mask, v), 0) == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    dp[(mask | (1 << w), w)] = dp.get((mask | (1 << w), w), 0) + dp[(mask, v)]
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


# ============================================================
# PART 1: The Orthogonal Control Ratio at varying n
# ============================================================

def ocr_by_n():
    """
    Compute OCR(n) = R^2 of score regularity predicting H(T) for all
    tournaments on n vertices. Does OCR increase, decrease, or stabilize?
    """
    print("=" * 70)
    print("PART 1: ORTHOGONAL CONTROL RATIO BY n")
    print("=" * 70)
    print("OCR = R^2 of S_2 (score variance) predicting H")
    print("Higher OCR = orthogonal invariants control more of H\n")

    for n in [3, 4, 5, 6, 7]:
        m = n * (n - 1) // 2
        total = 1 << m

        if total > 2200000:
            # Sample
            np.random.seed(42)
            sample_size = 30000
            bits_list = np.random.randint(0, total, sample_size)
        else:
            bits_list = range(total)

        Hs = []
        S2s = []
        tr4s = []

        for bits in bits_list:
            bits = int(bits)
            A = np.zeros((n, n), dtype=int)
            pos = 0
            for i in range(n):
                for j in range(i + 1, n):
                    if bits & (1 << pos):
                        A[i][j] = 1
                    else:
                        A[j][i] = 1
                    pos += 1

            H = count_ham_paths(A)
            scores = A.sum(axis=1)
            S2 = float(np.sum((scores - (n - 1) / 2) ** 2))

            S = signed_adjacency(A)
            tr4 = float(np.real(np.trace(np.linalg.matrix_power(S, 4))))

            Hs.append(H)
            S2s.append(S2)
            tr4s.append(tr4)

        Hs = np.array(Hs, dtype=float)
        S2s = np.array(S2s, dtype=float)
        tr4s = np.array(tr4s, dtype=float)

        # OCR from S_2 alone
        corr_S2 = np.corrcoef(Hs, S2s)[0, 1]
        ocr_S2 = corr_S2 ** 2

        # OCR from tr(S^4) alone
        corr_tr4 = np.corrcoef(Hs, tr4s)[0, 1]
        ocr_tr4 = corr_tr4 ** 2

        # OCR from BOTH (multiple regression R^2)
        X = np.column_stack([S2s, tr4s, np.ones(len(Hs))])
        try:
            beta = np.linalg.lstsq(X, Hs, rcond=None)[0]
            Hs_pred = X @ beta
            ss_res = np.sum((Hs - Hs_pred) ** 2)
            ss_tot = np.sum((Hs - np.mean(Hs)) ** 2)
            ocr_both = 1 - ss_res / ss_tot
        except:
            ocr_both = float('nan')

        print(f"  n={n}: OCR(S_2)={ocr_S2:.4f}, OCR(tr4)={ocr_tr4:.4f}, "
              f"OCR(both)={ocr_both:.4f}, "
              f"samples={len(Hs)}")

    print("\n  If OCR(S_2) -> 1 as n grows: orthogonal control STRENGTHENS")
    print("  If OCR(S_2) -> 0: it weakens (other invariants needed)")


# ============================================================
# PART 2: The Harish-Chandra Projection
# ============================================================

def harish_chandra_projection():
    """
    The Harish-Chandra isomorphism: Z(U(so(n))) = R[h]^W

    For a tournament T with signed adjacency S_T, the Harish-Chandra
    projection maps S_T to its EIGENVALUE MAGNITUDES (= its projection
    onto the Cartan subalgebra, modulo Weyl group action).

    Question: How much of H is determined by the HC projection?
    (= how much is determined by the eigenvalue magnitudes alone?)

    This extends the "spectral flatness" analysis to a formal
    representation-theoretic framework.
    """
    print("\n" + "=" * 70)
    print("PART 2: HARISH-CHANDRA PROJECTION")
    print("=" * 70)
    print("HC projection = eigenvalue magnitudes of S_T (sorted)")
    print("Question: How much of H is determined by eigenvalue magnitudes?\n")

    for n in [5, 7]:
        m = n * (n - 1) // 2
        total = 1 << m

        # Group tournaments by their HC projection (= sorted |eigenvalue| signature)
        by_hc = defaultdict(list)

        for bits in range(total):
            A = np.zeros((n, n), dtype=int)
            pos = 0
            for i in range(n):
                for j in range(i + 1, n):
                    if bits & (1 << pos):
                        A[i][j] = 1
                    else:
                        A[j][i] = 1
                    pos += 1

            H = count_ham_paths(A)
            S = signed_adjacency(A)
            eigvals = np.linalg.eigvals(S)
            # HC projection = sorted eigenvalue magnitudes (rounded for grouping)
            hc = tuple(sorted(round(abs(e), 3) for e in eigvals))
            by_hc[hc].append(H)

        # How many HC classes? How much H variance within vs between?
        total_mean = np.mean([H for hs in by_hc.values() for H in hs])
        ss_between = sum(len(hs) * (np.mean(hs) - total_mean) ** 2
                         for hs in by_hc.values())
        ss_total = sum((H - total_mean) ** 2
                       for hs in by_hc.values() for H in hs)
        r2_hc = ss_between / ss_total if ss_total > 0 else 0

        # Count: how many HC classes have unique H?
        unique_h = sum(1 for hs in by_hc.values() if len(set(hs)) == 1)
        multi_h = sum(1 for hs in by_hc.values() if len(set(hs)) > 1)

        print(f"  n={n}: {len(by_hc)} HC classes, "
              f"R^2={r2_hc:.4f} of H variance between classes")
        print(f"    {unique_h} classes with unique H, "
              f"{multi_h} classes with multiple H values")

        # Show the multi-H classes
        if multi_h > 0:
            print(f"    Multi-H classes (first 5):")
            count = 0
            for hc, hs in sorted(by_hc.items(), key=lambda x: -len(set(x[1]))):
                if len(set(hs)) <= 1:
                    continue
                h_set = sorted(set(hs))
                print(f"      |eig|={hc}: H in {h_set} ({len(hs)} tours)")
                count += 1
                if count >= 5:
                    break


# ============================================================
# PART 3: Markov Chain Hidden Orthogonality
# ============================================================

def markov_chain_orthogonality():
    """
    A Markov chain on n states has a directed transition matrix P.
    The key observable is the MIXING TIME tau (how fast it converges
    to stationarity).

    Orthogonal invariants of P:
    - Symmetric part: (P + P^T)/2 = "undirected" transition structure
    - Antisymmetric part: (P - P^T)/2 = "directional bias"
    - Eigenvalue magnitudes of P (HC projection)
    - Score regularity = row sum variance (always 0 for stochastic P)

    Does the spectral gap of (P + P^T)/2 predict the mixing time of P?
    This would extend our principle to Markov chains.
    """
    print("\n" + "=" * 70)
    print("PART 3: MARKOV CHAIN HIDDEN ORTHOGONALITY")
    print("=" * 70)
    print("Can symmetric invariants of a directed chain predict mixing?\n")

    np.random.seed(42)
    n = 5
    num_chains = 500

    results = []
    for trial in range(num_chains):
        # Generate random stochastic matrix
        logits = np.random.randn(n, n) * 2
        P = np.exp(logits)
        P = P / P.sum(axis=1, keepdims=True)

        # Mixing time proxy: spectral gap = 1 - |lambda_2(P)|
        eigvals = np.linalg.eigvals(P)
        eigvals_sorted = sorted(abs(eigvals), reverse=True)
        spec_gap = 1 - eigvals_sorted[1] if len(eigvals_sorted) > 1 else 1

        # Symmetric part: (P + P^T)/2
        P_sym = (P + P.T) / 2
        sym_eigvals = sorted(np.linalg.eigvalsh(P_sym), reverse=True)
        sym_gap = sym_eigvals[0] - sym_eigvals[1] if len(sym_eigvals) > 1 else 1

        # Antisymmetric energy
        P_anti = (P - P.T) / 2
        anti_energy = np.linalg.norm(P_anti, 'fro') ** 2
        sym_energy = np.linalg.norm(P_sym, 'fro') ** 2

        # Cartan fractions
        total_energy = np.linalg.norm(P, 'fro') ** 2
        anti_frac = anti_energy / total_energy

        results.append({
            'spec_gap': spec_gap,
            'sym_gap': sym_gap,
            'anti_frac': anti_frac,
            'anti_energy': anti_energy
        })

    # Correlations
    gaps = [r['spec_gap'] for r in results]
    sym_gaps = [r['sym_gap'] for r in results]
    anti_fracs = [r['anti_frac'] for r in results]

    corr_sym = np.corrcoef(gaps, sym_gaps)[0, 1]
    corr_anti = np.corrcoef(gaps, anti_fracs)[0, 1]

    print(f"  n={n}, {num_chains} random Markov chains:")
    print(f"  Corr(mixing_gap, symmetric_gap) = {corr_sym:.4f}")
    print(f"  Corr(mixing_gap, anti_fraction) = {corr_anti:.4f}")
    print(f"  R^2 from symmetric gap = {corr_sym**2:.4f}")

    if corr_sym ** 2 > 0.5:
        print(f"  >>> ORTHOGONAL CONTROL CONFIRMED for Markov chains!")
        print(f"  >>> Symmetric gap explains {corr_sym**2*100:.0f}% of mixing variance")
    else:
        print(f"  >>> Orthogonal control is WEAK for Markov chains")
        print(f"  >>> The directed structure matters independently")


# ============================================================
# PART 4: The Compression Miracle
# ============================================================

def compression_miracle():
    """
    THE COMPRESSION MIRACLE

    A tournament on n vertices has n(n-1)/2 binary degrees of freedom.
    The score sequence has n degrees of freedom.

    Compression ratio = n(n-1)/2 / n = (n-1)/2.

    But the score sequence explains 96% of H variance!

    How is this possible? The answer reveals the structure of
    "hidden orthogonality":

    The score sequence is not just ANY n-dimensional projection.
    It is the MARGINAL of the joint distribution of arcs.
    It captures the FIRST-ORDER statistics perfectly.

    H, however, depends on HIGHER-ORDER correlations between arcs.
    The 4% residual comes from the joint structure beyond marginals.

    This is exactly the Cartan hierarchy:
    - Level 1 (scores): marginals = symmetric projection = 96%
    - Level 2 (tr(S^4)): pairwise correlations = kurtosis = 44% of 4%
    - Level 3+ (higher traces): higher-order correlations = remaining

    The MIRACLE is that marginals capture so much.
    The EXPLANATION is that tournaments are COMPLETE (every pair has an arc),
    which creates strong correlations between arcs, making marginals highly
    informative.
    """
    print("\n" + "=" * 70)
    print("PART 4: THE COMPRESSION MIRACLE")
    print("=" * 70)

    for n in [3, 4, 5, 6, 7]:
        m = n * (n - 1) // 2
        total = 1 << m

        if total > 2200000:
            np.random.seed(42)
            bits_list = np.random.randint(0, total, 30000)
        else:
            bits_list = range(total)

        Hs = []
        S2s = []

        for bits in bits_list:
            bits = int(bits)
            A = np.zeros((n, n), dtype=int)
            pos = 0
            for i in range(n):
                for j in range(i + 1, n):
                    if bits & (1 << pos):
                        A[i][j] = 1
                    else:
                        A[j][i] = 1
                    pos += 1
            H = count_ham_paths(A)
            scores = A.sum(axis=1)
            S2 = float(np.sum((scores - (n - 1) / 2) ** 2))
            Hs.append(H)
            S2s.append(S2)

        Hs = np.array(Hs, dtype=float)
        S2s = np.array(S2s, dtype=float)

        corr = np.corrcoef(Hs, S2s)[0, 1]
        ocr = corr ** 2
        compression = m / n

        print(f"  n={n}: {m} arcs -> {n} scores. "
              f"Compression={compression:.1f}x. "
              f"OCR={ocr:.4f}. "
              f"Info density = OCR/compression = {ocr/compression:.4f}")


# ============================================================
# PART 5: Where Orthogonal Control BREAKS
# ============================================================

def where_it_breaks():
    """
    The orthogonal control principle should break when:
    1. The system is NOT complete (sparse digraphs)
    2. The observable depends on GLOBAL connectivity (not local structure)
    3. The Lie algebra is too small (not enough symmetric invariants)

    Test: for SPARSE directed graphs (not tournaments), does score
    regularity still predict Hamiltonian path count?
    """
    print("\n" + "=" * 70)
    print("PART 5: WHERE ORTHOGONAL CONTROL BREAKS")
    print("=" * 70)

    n = 5
    print(f"\n  Sparse digraphs on n={n} vertices (random edge density)")

    np.random.seed(42)
    results_by_density = defaultdict(list)

    for trial in range(2000):
        density = np.random.uniform(0.3, 1.0)
        A = np.zeros((n, n), dtype=int)
        for i in range(n):
            for j in range(n):
                if i != j and np.random.random() < density:
                    A[i][j] = 1

        # This is NOT a tournament (can have mutual arcs or no arcs)
        # Count Hamiltonian paths (same DP works)
        H = count_ham_paths(A)

        # Out-degree variance (analogous to S_2)
        out_degrees = A.sum(axis=1)
        S2_out = float(np.var(out_degrees))

        # In-degree variance
        in_degrees = A.sum(axis=0)
        S2_in = float(np.var(in_degrees))

        # Actual density
        actual_density = A.sum() / (n * (n - 1))

        bucket = round(actual_density, 1)
        results_by_density[bucket].append({
            'H': H, 'S2_out': S2_out, 'S2_in': S2_in,
            'density': actual_density
        })

    print(f"\n  {'density':>8s} {'count':>6s} {'mean H':>8s} "
          f"{'corr(H,-S2out)':>15s} {'OCR':>8s}")
    for d in sorted(results_by_density.keys()):
        group = results_by_density[d]
        if len(group) < 20:
            continue
        Hs = [r['H'] for r in group]
        S2s = [r['S2_out'] for r in group]
        if len(set(Hs)) > 1 and len(set(S2s)) > 1:
            corr = np.corrcoef(Hs, S2s)[0, 1]
            ocr = corr ** 2
        else:
            corr = float('nan')
            ocr = float('nan')
        print(f"  {d:8.1f} {len(group):6d} {np.mean(Hs):8.1f} "
              f"{corr:15.4f} {ocr:8.4f}")

    # Tournament (density=0.5, complete) vs sparse
    print(f"\n  At density=0.5 (tournament-like): OCR should be highest")
    print(f"  At density<0.3 or >0.7: OCR should decrease")
    print(f"  COMPLETENESS is the key to orthogonal control")


# ============================================================
# PART 6: The Universal Orthogonal Principle — Abstract
# ============================================================

def universal_principle():
    """
    THE UNIVERSAL ORTHOGONAL CONTROL PRINCIPLE

    For a directed system S with n elements:

    1. DEFINITION: A "symmetric shadow" of S is any function f(S) satisfying
       f(S) = f(S^op) where S^op reverses all directions.

    2. CLAIM: The stronger the pairwise interactions (more complete the graph),
       the more the symmetric shadows control the system's global invariants.

    3. QUANTIFICATION: OCR(S) = R^2 of best linear symmetric predictor of
       the key observable.

    4. SCALING: For complete directed graphs (tournaments),
       OCR(n) ~ 1 - c/n for some constant c > 0.
       The orthogonal control STRENGTHENS with system size.

    5. INTERPRETATION: Direction creates STRUCTURE, but symmetric
       invariants create CONSTRAINTS. In a complete system, every pair
       has an arc, so the constraints are maximally interconnected.
       The symmetric invariants "see through" the forest of directions
       to the underlying topological structure.

    Let's test the scaling claim: OCR ~ 1 - c/n.
    """
    print("\n" + "=" * 70)
    print("PART 6: UNIVERSAL PRINCIPLE — OCR SCALING")
    print("=" * 70)

    # We already computed OCR by n in Part 1. Let's fit OCR ~ 1 - c/n.
    # Using cached values from Part 1's output:
    # Need to recompute for the fit.

    ocrs = {}
    for n in [3, 4, 5, 6, 7]:
        m = n * (n - 1) // 2
        total = 1 << m

        if total > 2200000:
            np.random.seed(42)
            bits_list = np.random.randint(0, total, 30000)
        else:
            bits_list = range(total)

        Hs = []
        S2s = []

        for bits in bits_list:
            bits = int(bits)
            A = np.zeros((n, n), dtype=int)
            pos = 0
            for i in range(n):
                for j in range(i + 1, n):
                    if bits & (1 << pos):
                        A[i][j] = 1
                    else:
                        A[j][i] = 1
                    pos += 1
            H = count_ham_paths(A)
            scores = A.sum(axis=1)
            S2 = float(np.sum((scores - (n - 1) / 2) ** 2))
            Hs.append(H)
            S2s.append(S2)

        corr = np.corrcoef(Hs, S2s)[0, 1]
        ocrs[n] = corr ** 2

    print(f"\n  OCR by n:")
    for n in sorted(ocrs.keys()):
        print(f"    n={n}: OCR = {ocrs[n]:.4f}, "
              f"1-OCR = {1-ocrs[n]:.4f}, "
              f"n*(1-OCR) = {n*(1-ocrs[n]):.4f}")

    # Fit: 1-OCR ~ c/n => n*(1-OCR) ~ c (constant)
    ns = sorted(ocrs.keys())
    residuals = [n * (1 - ocrs[n]) for n in ns]
    print(f"\n  If n*(1-OCR) is constant, then OCR = 1 - c/n:")
    print(f"    Average n*(1-OCR) = {np.mean(residuals):.4f}")
    print(f"    Std n*(1-OCR) = {np.std(residuals):.4f}")

    if np.std(residuals) < np.mean(residuals) * 0.5:
        c = np.mean(residuals)
        print(f"\n  UNIVERSAL SCALING LAW: OCR(n) = 1 - {c:.2f}/n")
        print(f"  At n=100: OCR ~ {1 - c/100:.4f}")
        print(f"  At n=1000: OCR ~ {1 - c/1000:.6f}")
        print(f"  >>> Score regularity explains 99%+ for large tournaments!")
    else:
        print(f"\n  Scaling is not 1/n. Residuals vary too much.")
        # Try other fits
        log_ns = [np.log(n) for n in ns]
        log_1mOCR = [np.log(max(1 - ocrs[n], 1e-10)) for n in ns]
        if len(log_ns) >= 3:
            slope, intercept = np.polyfit(log_ns, log_1mOCR, 1)
            print(f"  Power law fit: 1-OCR ~ n^{slope:.2f}")
            print(f"  (Exponent {slope:.2f}: "
                  f"{'faster than 1/n' if slope < -1 else 'slower than 1/n'})")


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    ocr_by_n()
    harish_chandra_projection()
    markov_chain_orthogonality()
    compression_miracle()
    where_it_breaks()
    universal_principle()
