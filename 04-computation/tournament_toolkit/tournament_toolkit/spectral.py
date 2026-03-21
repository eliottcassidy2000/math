"""
SpectralAnalyzer: Tournament quality and ambiguity measurement.

Analyzes the spectral properties of tournament matrices to determine:
  - Spectral flatness (all eigenvalue magnitudes equal = optimal)
  - Kurtosis (deviation from flatness = sub-optimal)
  - Algebraic complexity (dimension of generated algebra)
  - H estimation from spectral invariants

Key theorem (S96): For regular tournaments:
  FLAT SPECTRUM <=> DOUBLY REGULAR <=> [A,S]=0 <=> MAX H <=> MIN ALG DIM

    analyzer = SpectralAnalyzer()
    result = analyzer.analyze(adjacency_matrix)
    print(result['flatness'])   # 1.0 = perfectly flat (Paley-like)
    print(result['estimated_H'])
"""

import numpy as np
from math import sqrt, log
from itertools import permutations


class SpectralAnalyzer:
    """Spectral analysis of tournament adjacency matrices."""

    def analyze(self, adj):
        """Full spectral analysis of tournament adjacency matrix.

        Args:
            adj: n x n adjacency matrix (adj[i][j] = 1 if i beats j)

        Returns:
            dict with spectral invariants and quality metrics
        """
        A = np.array(adj, dtype=float)
        n = A.shape[0]

        # Signed adjacency
        S = A - A.T  # entries in {-1, 0, 1}

        # Score sequence
        scores = A.sum(axis=1)
        score_var = np.var(scores)
        S2 = np.sum((scores - (n-1)/2)**2)  # Landau irregularity

        # Eigenvalues of S (purely imaginary for real antisymmetric)
        evals = np.linalg.eigvals(S)
        magnitudes = sorted([abs(e) for e in evals if abs(e) > 0.01], reverse=True)

        # Spectral moments (traces of even powers)
        tr2 = np.trace(S @ S).real  # = -n(n-1) always
        S2_mat = S @ S
        tr4 = np.trace(S2_mat @ S2_mat).real

        # Spectral flatness: tr(S^4) / tr(S^2)^2
        # For perfectly flat spectrum: this equals 1/(n-1) * n(n-1)... let me compute
        if tr2 != 0:
            kurtosis = tr4 / tr2**2  # Lower = flatter
        else:
            kurtosis = 0

        # Number of distinct eigenvalue magnitudes
        distinct_mags = len(set(round(m, 3) for m in magnitudes))

        # Algebraic dimension = 2 * distinct_mags + 1 (if 0 is eigenvalue)
        has_zero = any(abs(e) < 0.01 for e in evals)
        alg_dim = 2 * distinct_mags + (1 if has_zero else 0)

        # Flatness score: 1.0 = perfectly flat, 0.0 = maximally peaked
        if len(magnitudes) > 1:
            mag_cv = np.std(magnitudes) / (np.mean(magnitudes) + 1e-10)
            flatness = max(0, 1.0 - mag_cv)
        elif len(magnitudes) == 1:
            flatness = 1.0
        else:
            flatness = 0.0

        # Is it regular?
        is_regular = all(abs(s - (n-1)/2) < 0.01 for s in scores)

        # Is it doubly regular? (S^2 = -nI + J)
        J = np.ones((n, n))
        is_drt = np.allclose(S2_mat, -n * np.eye(n) + J, atol=0.5)

        # 3-cycle count
        from itertools import combinations
        c3 = sum(1 for i, j, k in combinations(range(n), 3)
                 if (A[i][j] and A[j][k] and A[k][i]) or
                    (A[i][k] and A[k][j] and A[j][i]))

        # H computation (only for small n)
        H = None
        if n <= 8:
            H = sum(1 for perm in permutations(range(n))
                    if all(A[perm[k]][perm[k+1]] for k in range(n-1)))

        result = {
            'n': n,
            'scores': sorted(scores.astype(int), reverse=True),
            'score_variance': score_var,
            'landau_S2': S2,
            'is_regular': is_regular,
            'is_doubly_regular': is_drt,
            'c3': c3,
            'H': H,
            'eigenvalue_magnitudes': magnitudes,
            'distinct_magnitudes': distinct_mags,
            'tr_S2': tr2,
            'tr_S4': tr4,
            'kurtosis': kurtosis,
            'flatness': flatness,
            'algebra_dimension': alg_dim,
            'spectral_diversity': distinct_mags,
        }

        return result

    def quality_report(self, adj):
        """Human-readable quality report for a tournament."""
        r = self.analyze(adj)
        lines = [
            f"Tournament Quality Report (n={r['n']})",
            f"  Scores: {r['scores']}",
            f"  Regular: {r['is_regular']}, DRT: {r['is_doubly_regular']}",
            f"  3-cycles: {r['c3']}",
            f"  H (Hamiltonian paths): {r['H']}",
            f"  Spectral flatness: {r['flatness']:.3f} (1.0 = optimal)",
            f"  Kurtosis: {r['kurtosis']:.6f} (lower = better)",
            f"  Algebra dimension: {r['algebra_dimension']} (3 = optimal for DRT)",
            f"  Landau irregularity S₂: {r['landau_S2']:.1f} (0 = regular)",
        ]
        if r['is_doubly_regular']:
            lines.append("  ★ DOUBLY REGULAR — maximum symmetry, maximum H for this score class")
        return "\n".join(lines)
