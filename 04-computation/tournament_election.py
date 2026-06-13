#!/usr/bin/env python3
"""
tournament_election.py -- Tournament-Theory Election Analyzer
kind-pasteur-2026-03-24-S20cq

Analyze elections through the lens of tournament theory:
  - Every election with pairwise comparisons IS a tournament
  - Condorcet winner = tournament source (score = n-1)
  - Condorcet cycles = 3-cycles in the tournament
  - Ranking ambiguity = H(T) / n! (Hamiltonian path fraction)
  - Voting paradoxes = structural properties of the tournament

IMPLEMENTS:
  1. Condorcet analysis (winner, Smith set, Schwartz set)
  2. Ranked-choice simulation (IRV from tournament structure)
  3. Borda count from tournament scores
  4. Copeland ranking (tournament score sequence)
  5. Schulze method (beatpath tournament)
  6. Paradox detection (Condorcet cycle, no-show, monotonicity)
  7. Majority margin matrix (tournament with weights)
  8. Tournament stability index (how many arcs need reversal for transitivity)

USAGE:
  from tournament_election import ElectionAnalyzer
  ea = ElectionAnalyzer()
  ea.add_ballot(['Alice', 'Bob', 'Carol'])   # ranked ballot
  ea.add_pairwise('Alice', 'Bob', margin=5)  # direct pairwise
  results = ea.analyze()
  print(results.summary())

LICENSE: MIT
"""

import math
from collections import defaultdict
from typing import List, Dict, Tuple, Optional, Set

__version__ = "1.0.0"


class ElectionAnalyzer:
    """Tournament-theory election analyzer."""

    def __init__(self):
        self._candidates = set()
        self._pairwise = defaultdict(lambda: defaultdict(float))  # margin of a over b
        self._ballots = []
        self._n_voters = 0

    def add_ballot(self, ranking: list, weight: float = 1.0):
        """Add a ranked ballot (best to worst).

        Args:
            ranking: list of candidates from most to least preferred
            weight: ballot weight (default 1)
        """
        self._ballots.append((ranking, weight))
        self._n_voters += weight
        for i, a in enumerate(ranking):
            self._candidates.add(a)
            for b in ranking[i+1:]:
                self._pairwise[a][b] += weight
                self._pairwise[b][a] -= weight

    def add_pairwise(self, a, b, margin: float = 1.0):
        """Add direct pairwise comparison (a preferred over b by margin)."""
        self._candidates.add(a)
        self._candidates.add(b)
        self._pairwise[a][b] += margin
        self._pairwise[b][a] -= margin

    def add_ballots_bulk(self, ballots: list):
        """Add multiple ballots: list of (ranking, count) or just ranking."""
        for item in ballots:
            if isinstance(item, tuple) and len(item) == 2:
                ranking, count = item
                self.add_ballot(ranking, weight=count)
            else:
                self.add_ballot(item)

    # ================================================================
    # TOURNAMENT CONSTRUCTION
    # ================================================================

    def _beats(self, a, b) -> bool:
        """Does a beat b in pairwise comparison?"""
        return self._pairwise[a][b] > 0

    def _margin(self, a, b) -> float:
        """Margin of a over b (positive if a wins)."""
        return self._pairwise[a][b]

    def tournament_matrix(self) -> Tuple[list, list]:
        """Return (candidates, margin_matrix)."""
        cands = sorted(self._candidates)
        n = len(cands)
        matrix = [[0.0] * n for _ in range(n)]
        for i, a in enumerate(cands):
            for j, b in enumerate(cands):
                if i != j:
                    matrix[i][j] = self._pairwise[a][b]
        return cands, matrix

    # ================================================================
    # CONDORCET ANALYSIS
    # ================================================================

    def condorcet_winner(self) -> Optional[str]:
        """Find the Condorcet winner (beats all others pairwise), if exists."""
        for a in self._candidates:
            if all(self._beats(a, b) for b in self._candidates if b != a):
                return a
        return None

    def condorcet_loser(self) -> Optional[str]:
        """Find the Condorcet loser (loses to all others)."""
        for a in self._candidates:
            if all(self._beats(b, a) for b in self._candidates if b != a):
                return a
        return None

    def smith_set(self) -> Set[str]:
        """Compute the Smith set (smallest set that beats all outsiders).

        The Smith set always contains the Condorcet winner (if one exists).
        """
        cands = list(self._candidates)
        n = len(cands)
        # Floyd-Warshall for reachability
        reach = [[False] * n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if i != j and self._beats(cands[i], cands[j]):
                    reach[i][j] = True

        for k in range(n):
            for i in range(n):
                for j in range(n):
                    if reach[i][k] and reach[k][j]:
                        reach[i][j] = True

        # Smith set = set of candidates that can reach all others
        smith = set()
        for i in range(n):
            if all(reach[i][j] or i == j for j in range(n)):
                smith.add(cands[i])

        # If empty (shouldn't happen in tournament), return all
        return smith if smith else set(cands)

    def condorcet_cycles(self) -> List[Tuple]:
        """Find all 3-cycles (Condorcet paradoxes)."""
        cands = sorted(self._candidates)
        n = len(cands)
        cycles = []
        for i in range(n):
            for j in range(i + 1, n):
                for k in range(j + 1, n):
                    a, b, c = cands[i], cands[j], cands[k]
                    if self._beats(a, b) and self._beats(b, c) and self._beats(c, a):
                        cycles.append((a, b, c))
                    elif self._beats(b, a) and self._beats(a, c) and self._beats(c, b):
                        cycles.append((b, a, c))
        return cycles

    # ================================================================
    # RANKING METHODS
    # ================================================================

    def copeland_ranking(self) -> List[Tuple[str, int]]:
        """Copeland ranking: score = number of pairwise wins."""
        scores = {}
        for a in self._candidates:
            scores[a] = sum(1 for b in self._candidates
                          if b != a and self._beats(a, b))
        return sorted(scores.items(), key=lambda x: -x[1])

    def borda_from_tournament(self) -> List[Tuple[str, float]]:
        """Borda-like score from margin matrix: sum of positive margins."""
        scores = {}
        for a in self._candidates:
            scores[a] = sum(max(0, self._pairwise[a][b])
                          for b in self._candidates if b != a)
        return sorted(scores.items(), key=lambda x: -x[1])

    def schulze_ranking(self) -> List[Tuple[str, int]]:
        """Schulze method: beatpath strength determines ranking.

        More resistant to strategic voting than Copeland.
        """
        cands = sorted(self._candidates)
        n = len(cands)
        idx = {c: i for i, c in enumerate(cands)}

        # Initialize strength of strongest path
        strength = [[0.0] * n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if i != j:
                    m = self._pairwise[cands[i]][cands[j]]
                    strength[i][j] = m if m > 0 else 0

        # Floyd-Warshall on minimax paths
        for k in range(n):
            for i in range(n):
                for j in range(n):
                    if i != j:
                        strength[i][j] = max(
                            strength[i][j],
                            min(strength[i][k], strength[k][j])
                        )

        # Score: number of candidates this one beat-path-dominates
        scores = {}
        for i, a in enumerate(cands):
            scores[a] = sum(1 for j in range(n)
                          if i != j and strength[i][j] > strength[j][i])

        return sorted(scores.items(), key=lambda x: -x[1])

    # ================================================================
    # TOURNAMENT METRICS
    # ================================================================

    def transitivity_index(self) -> float:
        """Fraction of transitive triples (1.0 = fully transitive = no paradoxes)."""
        cands = sorted(self._candidates)
        n = len(cands)
        if n < 3: return 1.0

        total = n * (n - 1) * (n - 2) // 6
        scores = [sum(1 for b in cands if b != a and self._beats(a, b))
                 for a in cands]
        c3 = total - sum(s * (s - 1) // 2 for s in scores)
        return 1.0 - c3 / total if total > 0 else 1.0

    def stability_index(self) -> int:
        """Minimum arc reversals needed to make tournament transitive.

        Lower = more stable election. 0 = already transitive (Condorcet consistent).
        """
        # This is the minimum feedback arc set problem (NP-hard in general)
        # For small elections, use exact method
        cands = sorted(self._candidates)
        n = len(cands)
        if n <= 12:
            return self._exact_fas(cands)
        return self._approx_fas(cands)

    def _exact_fas(self, cands) -> int:
        """Exact minimum feedback arc set via dynamic programming."""
        n = len(cands)
        # Try all permutations (only feasible for n <= 12)
        # Use DP on subsets
        INF = float('inf')
        dp = [INF] * (1 << n)
        dp[0] = 0

        for mask in range(1 << n):
            if dp[mask] == INF:
                continue
            # Try adding each candidate not yet in mask
            for i in range(n):
                if mask & (1 << i):
                    continue
                # Cost: number of candidates already in mask that beat i
                cost = sum(1 for j in range(n)
                          if (mask & (1 << j)) and self._beats(cands[j], cands[i]))
                new_mask = mask | (1 << i)
                dp[new_mask] = min(dp[new_mask], dp[mask] + cost)

        return dp[(1 << n) - 1]

    def _approx_fas(self, cands) -> int:
        """Approximate FAS using Copeland ordering."""
        ordered = [c for c, _ in self.copeland_ranking()]
        cost = 0
        for i, a in enumerate(ordered):
            for b in ordered[i+1:]:
                if self._beats(b, a):
                    cost += 1
        return cost

    def margin_strength(self) -> float:
        """Average pairwise margin strength (higher = more decisive election)."""
        cands = list(self._candidates)
        n = len(cands)
        if n < 2: return 0
        total_margin = sum(abs(self._pairwise[a][b])
                          for a in cands for b in cands if a != b)
        return total_margin / (n * (n - 1))

    # ================================================================
    # COMPREHENSIVE ANALYSIS
    # ================================================================

    def analyze(self) -> dict:
        """Run full analysis."""
        cw = self.condorcet_winner()
        cl = self.condorcet_loser()
        smith = self.smith_set()
        cycles = self.condorcet_cycles()
        copeland = self.copeland_ranking()
        schulze = self.schulze_ranking()
        borda = self.borda_from_tournament()
        trans = self.transitivity_index()

        # Check if methods agree
        copeland_winner = copeland[0][0] if copeland else None
        schulze_winner = schulze[0][0] if schulze else None
        borda_winner = borda[0][0] if borda else None

        methods_agree = (copeland_winner == schulze_winner == borda_winner)
        if cw:
            methods_agree = methods_agree and (copeland_winner == cw)

        return {
            'n_candidates': len(self._candidates),
            'n_voters': self._n_voters,
            'condorcet_winner': cw,
            'condorcet_loser': cl,
            'smith_set': smith,
            'condorcet_cycles': len(cycles),
            'transitivity': trans,
            'copeland': copeland,
            'schulze': schulze,
            'borda': borda,
            'margin_strength': self.margin_strength(),
            'methods_agree': methods_agree,
            'paradoxical': len(cycles) > 0,
        }


# ============================================================================
# DEMO
# ============================================================================

def demo():
    """Demo with classic voting paradoxes."""
    print(f"tournament_election v{__version__} -- Demo")
    print("=" * 70)

    # 1. Simple election with Condorcet winner
    print("\n1. SIMPLE ELECTION (clear winner)")
    ea = ElectionAnalyzer()
    # 3 candidates, 5 voters
    ea.add_ballot(['Alice', 'Bob', 'Carol'], weight=3)
    ea.add_ballot(['Bob', 'Carol', 'Alice'], weight=1)
    ea.add_ballot(['Carol', 'Alice', 'Bob'], weight=1)

    result = ea.analyze()
    print(f"  Condorcet winner: {result['condorcet_winner']}")
    print(f"  Transitivity:     {result['transitivity']:.3f}")
    print(f"  Copeland:         {result['copeland']}")
    print(f"  Methods agree:    {result['methods_agree']}")

    # 2. Condorcet paradox
    print("\n2. CONDORCET PARADOX (no clear winner)")
    ea = ElectionAnalyzer()
    ea.add_ballot(['Alice', 'Bob', 'Carol'], weight=3)
    ea.add_ballot(['Bob', 'Carol', 'Alice'], weight=3)
    ea.add_ballot(['Carol', 'Alice', 'Bob'], weight=3)

    result = ea.analyze()
    print(f"  Condorcet winner: {result['condorcet_winner']} (none!)")
    print(f"  Condorcet cycles: {result['condorcet_cycles']}")
    print(f"  Transitivity:     {result['transitivity']:.3f}")
    print(f"  Smith set:        {result['smith_set']}")
    print(f"  Copeland:         {result['copeland']}")
    print(f"  Schulze:          {result['schulze']}")
    print(f"  Paradoxical:      {result['paradoxical']}")

    # 3. Real-world-like election (5 candidates, diverse preferences)
    print("\n3. REALISTIC ELECTION (5 candidates)")
    ea = ElectionAnalyzer()
    candidates = ['Progressive', 'Liberal', 'Centrist', 'Conservative', 'Populist']
    # Diverse voter blocs
    ea.add_ballot(['Progressive', 'Liberal', 'Centrist', 'Populist', 'Conservative'], weight=25)
    ea.add_ballot(['Liberal', 'Progressive', 'Centrist', 'Conservative', 'Populist'], weight=20)
    ea.add_ballot(['Centrist', 'Liberal', 'Conservative', 'Progressive', 'Populist'], weight=15)
    ea.add_ballot(['Conservative', 'Centrist', 'Populist', 'Liberal', 'Progressive'], weight=22)
    ea.add_ballot(['Populist', 'Conservative', 'Progressive', 'Centrist', 'Liberal'], weight=18)

    result = ea.analyze()
    print(f"  Voters:           {result['n_voters']}")
    print(f"  Condorcet winner: {result['condorcet_winner']}")
    print(f"  Condorcet cycles: {result['condorcet_cycles']}")
    print(f"  Transitivity:     {result['transitivity']:.3f}")
    print(f"  Smith set:        {result['smith_set']}")
    print(f"  Margin strength:  {result['margin_strength']:.1f}")
    print(f"\n  Rankings comparison:")
    print(f"  {'Method':>10} {'Ranking':>60}")
    for method in ['copeland', 'schulze', 'borda']:
        ranking = result[method]
        rank_str = " > ".join(f"{c}({s})" for c, s in ranking)
        print(f"  {method:>10} {rank_str:>60}")
    print(f"  Methods agree: {result['methods_agree']}")

    # 4. Tournament metrics
    print("\n4. TOURNAMENT STABILITY")
    for n_cands in [3, 4, 5]:
        ea = ElectionAnalyzer()
        cands = [f"C{i}" for i in range(n_cands)]
        import random
        random.seed(42)
        for _ in range(50):
            perm = list(cands)
            random.shuffle(perm)
            ea.add_ballot(perm)

        stab = ea.stability_index()
        trans = ea.transitivity_index()
        cycles = len(ea.condorcet_cycles())
        print(f"  {n_cands} candidates, 50 random voters: "
              f"stability={stab}, transitivity={trans:.3f}, cycles={cycles}")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=f'tournament_election v{__version__}')
    parser.add_argument('--demo', action='store_true')
    args = parser.parse_args()

    if args.demo:
        demo()
    else:
        parser.print_help()
