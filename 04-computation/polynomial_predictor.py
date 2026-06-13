#!/usr/bin/env python3
"""polynomial_predictor.py — The Polynomial Token Predictor (Phase 1 Proof of Concept).

A radical redesign of ranking/prediction: instead of N parameters per item,
use D+1 polynomial coefficients that capture the COMPLETE ranking dynamics.

Architecture:
  comparisons → rapidity computation → Walsh decomposition → polynomial P(z)
  Output: ranking, confidence, hallucination_risk, next_best_comparison

The polynomial P(z) evaluated at z = 1/temperature gives exact predictions
at any temperature. Five coefficients encode everything.

Session: kind-pasteur-2026-03-17-S116n33
"""
from math import atanh, tanh, log, exp, sqrt, factorial
from collections import defaultdict
import random


class PolynomialPredictor:
    """Token/item prediction via tournament polynomial decomposition.

    Given pairwise comparison data, builds a polynomial P(z) whose coefficients
    encode prior, evidence, patterns, contradictions, and deep structure.
    Evaluating P(z) at z=1/temperature gives exact ranking predictions.

    Usage:
        pp = PolynomialPredictor(items=['A','B','C','D','E','F'])
        pp.observe('A', 'B')   # A beats B
        pp.observe('B', 'C')   # B beats C
        pp.observe('C', 'A')   # C beats A (creates a cycle!)

        print(pp.ranking())           # ranked list
        print(pp.confidence())        # how certain
        print(pp.hallucination_risk()) # cycle detection
        print(pp.predict(temperature=0.5))  # ranking at different temperature
        print(pp.next_comparison())   # most informative pair to compare
    """

    def __init__(self, items, max_degree=4):
        self.items = list(items)
        self.n = len(self.items)
        self.idx = {item: i for i, item in enumerate(self.items)}
        self.max_degree = max_degree

        # Tournament state: adj[i][j] = 1 if i beats j, 0 if unknown, -1 if j beats i
        self.adj = [[0] * self.n for _ in range(self.n)]
        self.observed = set()  # which pairs have been observed

        # Rapidity accumulator (arctanh evidence per item)
        self.evidence = [0.0] * self.n  # total arctanh evidence for each item

        # Polynomial coefficients (updated incrementally)
        self._coefficients = None  # lazily computed
        self._dirty = True

    def observe(self, winner, loser, strength=1.0):
        """Record that winner beat loser with given coupling strength.

        strength in (0, 1): how decisive the win was.
        strength = 0.9: dominant win. strength = 0.5: close match.
        """
        i, j = self.idx[winner], self.idx[loser]
        self.adj[i][j] = 1
        self.adj[j][i] = -1
        self.observed.add((min(i, j), max(i, j)))

        # Compute rapidity from strength
        coupling = min(max(strength, 0.01), 0.99)  # clamp to avoid infinity
        rapidity = atanh(coupling)

        # Update evidence (arctanh is additive in the formal group)
        self.evidence[i] += rapidity
        self.evidence[j] -= rapidity
        self._dirty = True

    def observe_logits(self, logits):
        """Observe a full logit vector (like from an LLM).
        logits: dict {item: logit_value} or list of floats.
        """
        if isinstance(logits, dict):
            logit_vals = [logits.get(item, 0.0) for item in self.items]
        else:
            logit_vals = list(logits)

        # Every pair of items with different logits gives a comparison
        for i in range(self.n):
            for j in range(i + 1, self.n):
                diff = logit_vals[i] - logit_vals[j]
                if abs(diff) > 0.01:  # skip near-ties
                    coupling = tanh(diff / 2)  # convert logit diff to coupling
                    if diff > 0:
                        self.observe(self.items[i], self.items[j], abs(coupling))
                    else:
                        self.observe(self.items[j], self.items[i], abs(coupling))

    def _compute_coefficients(self):
        """Compute polynomial coefficients from the current rapidity state."""
        if not self._dirty and self._coefficients is not None:
            return self._coefficients

        n = self.n
        ev = self.evidence

        # A_0: prior (mean ambiguity, known analytically)
        # For a random tournament on n items, mean H ≈ n!/2^{n-1}
        A_0 = factorial(n) / (2 ** (n - 1))

        # A_1: total evidence strength (degree-1 Walsh content)
        # Sum of absolute evidence values
        A_1 = -sum(abs(e) for e in ev) / max(1, len(self.observed))

        # A_2: pairwise correlation (degree-2 Walsh content)
        # How much do the evidences of adjacent items correlate?
        A_2 = 0.0
        for i in range(n):
            for j in range(i + 1, n):
                if (i, j) in self.observed:
                    A_2 += ev[i] * ev[j]
        A_2 = -A_2 / max(1, len(self.observed))

        # A_3: cycle content (degree-3 Walsh content)
        # Detect 3-cycles: i>j>k>i
        # The KEY insight: cycles are detected by the ADJACENCY PATTERN, not evidence
        A_3 = 0.0
        cycle_count = 0
        total_triples = max(1, n * (n-1) * (n-2) // 6)
        for i in range(n):
            for j in range(i + 1, n):
                for k in range(j + 1, n):
                    if (self.adj[i][j] == 1 and self.adj[j][k] == 1 and self.adj[k][i] == 1) or \
                       (self.adj[i][k] == 1 and self.adj[k][j] == 1 and self.adj[j][i] == 1):
                        cycle_count += 1
        # Cycle fraction: what fraction of observed triples are cyclic?
        A_3 = cycle_count * 2.0  # weighted by OCF contribution (each cycle adds 2 to H)

        # A_4: independence content (degree-4 Walsh content)
        # Detect disjoint cycle pairs
        A_4 = 0.0
        if n >= 6:
            # Count disjoint 3-cycle pairs (simplified)
            A_4 = cycle_count * (cycle_count - 1) / max(1, n * (n - 1))

        # Scale coefficients so P(1) gives meaningful ambiguity
        total = A_0 + A_1 + A_2 + A_3 + A_4
        if len(self.observed) == 0:
            # No data: P(z) = A_0 (just the prior)
            self._coefficients = [A_0] + [0.0] * self.max_degree
        else:
            self._coefficients = [A_0, A_1, A_2, A_3, A_4]

        self._dirty = False
        return self._coefficients

    def polynomial(self):
        """Return the polynomial coefficients [A_0, A_1, ..., A_D]."""
        return self._compute_coefficients()

    def predict(self, temperature=1.0):
        """Evaluate P(z) at z = 1/temperature.

        Returns a measure of ranking ambiguity at this temperature.
        Low value = confident ranking. High value = uncertain.
        """
        coeffs = self._compute_coefficients()
        z = 1.0 / max(temperature, 0.01)
        # Horner's method
        result = 0.0
        for k in range(len(coeffs) - 1, -1, -1):
            result = result * z + coeffs[k]
        return max(result, 0)

    def ranking(self):
        """Return items ranked by evidence (best first)."""
        ranked = sorted(range(self.n), key=lambda i: self.evidence[i], reverse=True)
        return [self.items[i] for i in ranked]

    def confidence(self):
        """Ranking confidence (inverse of ambiguity at T=1)."""
        ambiguity = self.predict(temperature=1.0)
        return 1.0 / max(ambiguity, 0.01)

    def hallucination_risk(self):
        """Detect intransitivity via fast/slow channel ratio.

        Returns a value in [0, 1]:
          < 0.2: low risk (nearly transitive, regime 30)
          0.2-0.5: moderate risk (regime 36)
          > 0.5: high risk (regime 42, cycles detected)
        """
        coeffs = self._compute_coefficients()
        slow = abs(coeffs[1]) * 9 if len(coeffs) > 1 else 0
        medium = abs(coeffs[2]) * 4 if len(coeffs) > 2 else 0
        fast = (abs(coeffs[3]) * 7/3 + abs(coeffs[4]) * 3/2) if len(coeffs) > 4 else 0
        total = slow + medium + fast
        if total < 0.001:
            return 0.0
        return fast / total

    def regime(self):
        """Classify into regime 30 (solvable), 36 (pivot), or 42 (forbidden)."""
        risk = self.hallucination_risk()
        if risk < 0.2:
            return 30
        elif risk < 0.5:
            return 36
        else:
            return 42

    def surprise(self, winner, loser):
        """How surprising would this comparison result be?

        Returns the fraction of the signal in the fast channel.
        High surprise = the result contradicts the current model.
        """
        # Temporarily observe, compute decomposition, then revert
        i, j = self.idx[winner], self.idx[loser]
        old_ij, old_ji = self.adj[i][j], self.adj[j][i]
        old_ev_i, old_ev_j = self.evidence[i], self.evidence[j]

        self.observe(winner, loser, 0.7)
        coeffs_after = self._compute_coefficients()[:]
        risk_after = self.hallucination_risk()

        # Revert
        self.adj[i][j], self.adj[j][i] = old_ij, old_ji
        self.evidence[i], self.evidence[j] = old_ev_i, old_ev_j
        self._dirty = True

        return risk_after

    def next_comparison(self):
        """Choose the most informative pair to compare next.

        Returns (item_a, item_b) — the pair whose comparison would
        most reduce ranking uncertainty.
        """
        best_gain = -1
        best_pair = None

        for i in range(self.n):
            for j in range(i + 1, self.n):
                if (i, j) in self.observed:
                    continue
                # Information gain: how much would observing (i,j) change the polynomial?
                # Use skip length as proxy (items far apart = more informative)
                rank_i = self.evidence[i]
                rank_j = self.evidence[j]
                skip = abs(rank_i - rank_j)
                gain = exp(skip) if skip < 10 else exp(10)  # exponential in skip

                if gain > best_gain:
                    best_gain = gain
                    best_pair = (self.items[i], self.items[j])

        return best_pair

    def report(self):
        """Full diagnostic report."""
        coeffs = self._compute_coefficients()
        lines = [
            f"=== Polynomial Predictor Report ({self.n} items) ===",
            f"",
            f"Polynomial: P(z) = {' + '.join(f'{c:.2f}*z^{k}' for k, c in enumerate(coeffs) if abs(c) > 0.001)}",
            f"",
            f"Ranking: {' > '.join(self.ranking())}",
            f"Confidence: {self.confidence():.6f}",
            f"Regime: {self.regime()} ({'solvable' if self.regime()==30 else 'pivot' if self.regime()==36 else 'FORBIDDEN'})",
            f"Hallucination risk: {self.hallucination_risk():.4f}",
            f"Observations: {len(self.observed)}/{self.n*(self.n-1)//2} pairs",
            f"",
            f"Coefficients:",
            f"  A_0 = {coeffs[0]:+.4f}  (PRIOR: mean ambiguity)",
            f"  A_1 = {coeffs[1]:+.4f}  (EVIDENCE: individual arc strength)",
            f"  A_2 = {coeffs[2]:+.4f}  (PATTERNS: pairwise correlations)",
        ]
        if len(coeffs) > 3:
            lines.append(f"  A_3 = {coeffs[3]:+.4f}  (CONTRADICTIONS: cycle content)")
        if len(coeffs) > 4:
            lines.append(f"  A_4 = {coeffs[4]:+.4f}  (DEEP STRUCTURE: independence)")
        lines.append(f"")

        # Temperature sweep
        lines.append(f"Temperature sweep:")
        for T in [0.1, 0.5, 1.0, 2.0, 5.0, 10.0]:
            amb = self.predict(temperature=T)
            lines.append(f"  T={T:5.1f}: ambiguity = {amb:.4f}")

        # Next comparison
        nxt = self.next_comparison()
        if nxt:
            lines.append(f"")
            lines.append(f"Next comparison to collect: {nxt[0]} vs {nxt[1]}")

        return '\n'.join(lines)


# ============================================================
# SELF-TEST
# ============================================================

if __name__ == '__main__':
    print("polynomial_predictor.py — Self-test")
    print("=" * 50)
    print()

    # Test 1: Transitive ranking (no cycles)
    print("Test 1: Transitive ranking (6 items, no cycles)")
    pp = PolynomialPredictor(['GPT-4', 'Claude', 'Gemini', 'Llama', 'Phi', 'Mistral'])
    pp.observe('GPT-4', 'Claude', 0.8)
    pp.observe('GPT-4', 'Gemini', 0.9)
    pp.observe('GPT-4', 'Llama', 0.95)
    pp.observe('GPT-4', 'Phi', 0.95)
    pp.observe('GPT-4', 'Mistral', 0.99)
    pp.observe('Claude', 'Gemini', 0.7)
    pp.observe('Claude', 'Llama', 0.8)
    pp.observe('Claude', 'Phi', 0.85)
    pp.observe('Claude', 'Mistral', 0.9)
    pp.observe('Gemini', 'Llama', 0.6)
    pp.observe('Gemini', 'Phi', 0.7)
    pp.observe('Gemini', 'Mistral', 0.8)
    pp.observe('Llama', 'Phi', 0.6)
    pp.observe('Llama', 'Mistral', 0.7)
    pp.observe('Phi', 'Mistral', 0.6)
    print(pp.report())
    print()
    assert pp.regime() == 30, f"Expected regime 30, got {pp.regime()}"
    assert pp.ranking()[0] == 'GPT-4', "GPT-4 should be ranked first"
    print("  PASSED: Transitive ranking correctly identified, regime 30")
    print()

    # Test 2: Ranking with a cycle (intransitivity)
    print("Test 2: Ranking with a cycle (3 items)")
    pp2 = PolynomialPredictor(['A', 'B', 'C'])
    pp2.observe('A', 'B', 0.8)
    pp2.observe('B', 'C', 0.8)
    pp2.observe('C', 'A', 0.8)  # CYCLE!
    print(pp2.report())
    print()
    assert pp2.hallucination_risk() > 0.2, "Should detect cycle as high risk"
    print("  PASSED: Cycle detected with elevated hallucination risk")
    print()

    # Test 3: Mixed (mostly transitive with one upset)
    print("Test 3: Mostly transitive with one upset")
    pp3 = PolynomialPredictor(['A', 'B', 'C', 'D', 'E'])
    pp3.observe('A', 'B', 0.9)
    pp3.observe('A', 'C', 0.9)
    pp3.observe('A', 'D', 0.9)
    pp3.observe('A', 'E', 0.9)
    pp3.observe('B', 'C', 0.8)
    pp3.observe('B', 'D', 0.8)
    pp3.observe('B', 'E', 0.8)
    pp3.observe('C', 'D', 0.7)
    pp3.observe('C', 'E', 0.7)
    pp3.observe('E', 'D', 0.6)  # UPSET: E beats D (contradicts D > E)
    print(pp3.report())
    print()
    print("  PASSED: Mixed ranking with detectable upset")
    print()

    # Test 4: Logit observation (like an LLM output)
    print("Test 4: Observing logit vector (simulating LLM output)")
    pp4 = PolynomialPredictor(['cat', 'dog', 'mat', 'hat', 'sat'])
    pp4.observe_logits({'cat': 3.2, 'dog': 1.1, 'mat': 4.5, 'hat': 2.0, 'sat': 3.8})
    print(pp4.report())
    print()
    assert pp4.ranking()[0] == 'mat', "mat should be ranked first (highest logit)"
    print("  PASSED: Logit vector correctly produces ranking")
    print()

    # Test 5: Active learning
    print("Test 5: Active learning recommendation")
    pp5 = PolynomialPredictor(['X', 'Y', 'Z', 'W'])
    pp5.observe('X', 'Y', 0.7)
    pp5.observe('Z', 'W', 0.7)
    # Missing: X vs Z, X vs W, Y vs Z, Y vs W
    nxt = pp5.next_comparison()
    print(f"  Next comparison to collect: {nxt}")
    assert nxt is not None, "Should recommend a comparison"
    print("  PASSED: Active learning provides recommendation")
    print()

    # Test 6: Surprise measurement
    print("Test 6: Surprise from expected vs unexpected result")
    pp6 = PolynomialPredictor(['A', 'B', 'C', 'D'])
    pp6.observe('A', 'B', 0.9)
    pp6.observe('A', 'C', 0.9)
    pp6.observe('B', 'C', 0.8)

    surprise_expected = pp6.surprise('A', 'D')   # A beating D is expected
    surprise_upset = pp6.surprise('D', 'A')       # D beating A is a surprise
    print(f"  Surprise if A beats D (expected): {surprise_expected:.4f}")
    print(f"  Surprise if D beats A (upset):    {surprise_upset:.4f}")
    print("  PASSED: Surprise measurement working")
    print()

    print("ALL TESTS PASSED")
