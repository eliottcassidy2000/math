"""
FormalRank: One-pass ranking from pairwise comparisons.

Uses the formal group logarithm arctanh to convert multiplicative
evidence into additive scores. O(m) time, O(n) space.

No iteration. No convergence criterion. No hyperparameters.
Equivalent to Bradley-Terry for direct evidence.

    ranker = FormalRank()
    ranker.add("Alice", "Bob")        # Alice beats Bob
    ranker.add("Alice", "Carol")
    ranker.add("Bob", "Carol")
    print(ranker.summary())
"""

from math import atanh, tanh, sqrt
from collections import defaultdict


class FormalRank:
    """One-pass ranking via formal group logarithm.

    For items i, j with w_ij wins for i and w_ji for j:
      coupling = (w_ij - w_ji) / (w_ij + w_ji)
      evidence = arctanh(coupling) = (1/2) * ln(w_ij / w_ji)

    Total score for item i = sum of evidence over all opponents.
    Ranking = sort by total score (descending).

    This is exact, one-pass, and O(m) for m comparisons.
    """

    def __init__(self):
        self._evidence = defaultdict(float)
        self._opponents = defaultdict(set)
        self._n_comparisons = defaultdict(int)
        self._wins = defaultdict(lambda: defaultdict(int))
        self._items = set()

    def add(self, winner, loser, wins=1, losses=0):
        """Record comparison(s). winner beat loser `wins` times, lost `losses` times."""
        self._items.add(winner)
        self._items.add(loser)
        self._wins[winner][loser] += wins
        self._wins[loser][winner] += losses
        total = wins + losses
        if total > 0:
            coupling = max(-0.9999, min(0.9999, (wins - losses) / total))
            ev = atanh(coupling)
            self._evidence[winner] += ev
            self._evidence[loser] -= ev
            self._opponents[winner].add(loser)
            self._opponents[loser].add(winner)
            self._n_comparisons[winner] += 1
            self._n_comparisons[loser] += 1

    def rank(self):
        """Return rankings as list of (item, score, confidence, n_opponents)."""
        results = []
        for item in sorted(self._items):
            theta = self._evidence[item]
            nc = max(self._n_comparisons[item], 1)
            conf = tanh(abs(theta) / sqrt(nc))
            results.append({
                'item': item,
                'score': theta,
                'confidence': conf,
                'opponents': len(self._opponents[item]),
                'comparisons': self._n_comparisons[item],
            })
        results.sort(key=lambda x: -x['score'])
        return results

    def summary(self, top_n=None):
        """Pretty-print the ranking."""
        results = self.rank()
        if top_n:
            results = results[:top_n]
        lines = [f"{'Rank':>4} | {'Item':>20} | {'Score':>8} | {'Conf':>6} | {'Opps':>4}"]
        lines.append("-" * 52)
        for i, r in enumerate(results, 1):
            lines.append(f"{i:4d} | {str(r['item']):>20} | {r['score']:+8.3f} | "
                        f"{r['confidence']:5.1%} | {r['opponents']:4d}")
        return "\n".join(lines)

    def ambiguity(self):
        """Compute the tournament's H-value (ambiguity measure).

        H = number of Hamiltonian paths = number of consistent total orderings.
        H = 1: unambiguous. H >> 1: highly ambiguous. H = 7 or 21: impossible.

        Only computable for small n (≤ 10). Returns None for larger tournaments.
        """
        items = sorted(self._items)
        n = len(items)
        if n > 10:
            return None

        # Build adjacency matrix from win counts
        idx = {item: i for i, item in enumerate(items)}
        adj = [[0] * n for _ in range(n)]
        for a in items:
            for b in items:
                if a != b and self._wins[a].get(b, 0) > self._wins[b].get(a, 0):
                    adj[idx[a]][idx[b]] = 1

        # Count Hamiltonian paths via DP
        from itertools import permutations
        return sum(1 for perm in permutations(range(n))
                   if all(adj[perm[k]][perm[k+1]] for k in range(n-1)))

    def most_informative_comparison(self):
        """Find the pair whose next comparison would most reduce ambiguity.

        Returns the pair (a, b) with the smallest |evidence|, meaning
        the least-decided head-to-head matchup.
        """
        best_pair = None
        min_evidence = float('inf')
        items = sorted(self._items)
        for i, a in enumerate(items):
            for b in items[i+1:]:
                w_ab = self._wins[a].get(b, 0)
                w_ba = self._wins[b].get(a, 0)
                total = w_ab + w_ba
                if total == 0:
                    return (a, b)  # Never compared = most informative
                coupling = abs(w_ab - w_ba) / total
                if coupling < min_evidence:
                    min_evidence = coupling
                    best_pair = (a, b)
        return best_pair

    @property
    def n_items(self):
        return len(self._items)

    @property
    def n_total_comparisons(self):
        return sum(self._n_comparisons.values()) // 2
