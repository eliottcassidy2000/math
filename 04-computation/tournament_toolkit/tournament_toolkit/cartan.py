"""
CartanProbe: AI confidence diagnostics via Cartan decomposition.

Any n x n matrix M decomposes as:
  M = A + S_tl + sigma * I

where:
  A = (M - M^T)/2           antisymmetric (TOURNAMENT: who beats whom)
  S_tl = (M + M^T)/2 - sI   symmetric traceless (COOPERATION: self-knowledge)
  sigma = Tr(M)/n            scalar (SCALE: overall confidence)

For LLM attention matrices, this decomposition reveals:
  - tournament_fraction: how competitive/decisive the attention is
  - cooperation_fraction: how much self-knowledge is present
  - entanglement: how coupled the competition and self-knowledge are

Key insight (proved S95): [A, S_tl] is SYMMETRIC, meaning
"interaction of competition and self-knowledge produces self-knowledge."

    probe = CartanProbe()
    report = probe.analyze(attention_matrix)
    print(report['verdict'])  # CONFIDENT / UNCERTAIN / RISKY
"""

import numpy as np
from math import tanh, atanh, sqrt


class CartanProbe:
    """Zero-parameter AI confidence diagnostic.

    Decomposes any matrix into tournament (competition) and cooperation
    (self-knowledge) sectors. The balance between them diagnoses whether
    a system is confident or uncertain.
    """

    def analyze(self, M, logits=None):
        """Full Cartan analysis of matrix M.

        Args:
            M: n x n numpy array (e.g., attention matrix)
            logits: optional logit vector for combined diagnosis

        Returns:
            dict with Cartan decomposition and diagnostic report
        """
        M = np.array(M, dtype=float)
        n = M.shape[0]

        # Cartan decomposition
        A = (M - M.T) / 2
        S_full = (M + M.T) / 2
        sigma = np.trace(S_full) / n
        S_tl = S_full - sigma * np.eye(n)

        # Frobenius norms
        norm2_M = np.sum(M**2)
        norm2_A = np.sum(A**2)
        norm2_S = np.sum(S_tl**2)
        norm2_sigma = n * sigma**2

        if norm2_M < 1e-15:
            return {'verdict': 'EMPTY', 'tournament_fraction': 0,
                    'cooperation_fraction': 0, 'entanglement': 0}

        tournament_frac = norm2_A / norm2_M
        cooperation_frac = norm2_S / norm2_M
        scalar_frac = norm2_sigma / norm2_M

        # Commutator [A, S_tl] — measures sector coupling
        comm = A @ S_tl - S_tl @ A
        norm2_comm = np.sum(comm**2)
        denom = sqrt(norm2_A * norm2_S) if norm2_A > 0 and norm2_S > 0 else 1.0
        entanglement = sqrt(norm2_comm) / denom

        # Decisive vs deliberating
        total_directed = tournament_frac + cooperation_frac
        decisiveness = tournament_frac / total_directed if total_directed > 0 else 0.5

        # Logit analysis if provided
        logit_risk = None
        if logits is not None:
            logits = np.array(logits).flatten()
            sorted_l = np.sort(logits)[::-1]
            if len(sorted_l) >= 2:
                gap = sorted_l[0] - sorted_l[1]
                logit_risk = 1.0 - tanh(gap)

        # Combined risk using formal group F(a,b) = (a+b)/(1+ab)
        cartan_risk = max(0, min(1, 1.0 - tournament_frac * 2))
        if logit_risk is not None:
            combined = (logit_risk + cartan_risk) / (1 + logit_risk * cartan_risk)
        else:
            combined = cartan_risk

        if combined < 0.2:
            verdict = 'CONFIDENT'
        elif combined < 0.5:
            verdict = 'UNCERTAIN'
        elif combined < 0.8:
            verdict = 'RISKY'
        else:
            verdict = 'HALLUCINATION_LIKELY'

        return {
            'tournament_fraction': tournament_frac,
            'cooperation_fraction': cooperation_frac,
            'scalar_fraction': scalar_frac,
            'scalar': sigma,
            'entanglement': entanglement,
            'decisiveness': decisiveness,
            'logit_risk': logit_risk,
            'cartan_risk': cartan_risk,
            'combined_risk': combined,
            'verdict': verdict,
            'antisymmetric': A,
            'symmetric_traceless': S_tl,
        }

    def quick_check(self, M):
        """Fast check: just return the verdict string."""
        return self.analyze(M)['verdict']

    def compare(self, M1, M2):
        """Compare two matrices (e.g., attention at different layers).

        Returns how the Cartan balance shifts between them.
        """
        r1 = self.analyze(M1)
        r2 = self.analyze(M2)
        return {
            'tournament_shift': r2['tournament_fraction'] - r1['tournament_fraction'],
            'cooperation_shift': r2['cooperation_fraction'] - r1['cooperation_fraction'],
            'entanglement_shift': r2['entanglement'] - r1['entanglement'],
            'verdict_change': f"{r1['verdict']} -> {r2['verdict']}",
        }
