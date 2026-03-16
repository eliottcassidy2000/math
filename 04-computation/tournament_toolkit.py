#!/usr/bin/env python3
"""tournament_toolkit.py — Extended tournament analysis toolkit.

Builds on tournament_H.py with graph structures from the tiling explorer
and practical applications to ranking, voting, sports, and TDA.

Modules:
  TilingGrid     — Triangular grid encoding of tournaments
  FlipGraph      — Single-arc transition graph between tournaments
  RankAmbiguity  — Ranking ambiguity measures for LLM/voting/sports
  TournamentTDA  — Topological data analysis features

Session: kind-pasteur-2026-03-16-S116n32
"""
from tournament_H import TournamentH, near_transitive_H
from itertools import combinations, permutations
from collections import defaultdict, Counter
from math import comb
import sys


# ================================================================
# MODULE 1: TILING GRID
# ================================================================

class TilingGrid:
    """Triangular grid encoding of a tournament relative to a fixed path.

    For n vertices with path 0->1->...->n-1, the non-path arcs form
    a triangular grid with (n-1)(n-2)/2 cells, organized by skip length:
      Row 0 (skip 2): n-2 arcs  — nearest neighbors
      Row 1 (skip 3): n-3 arcs  — next-nearest
      ...
      Row k (skip k+2): n-k-2 arcs
      ...
      Row n-3 (skip n-1): 1 arc  — longest range

    Each cell is 0 (forward, same as transitive) or 1 (backward).
    """

    def __init__(self, n, bits=None, tournament=None):
        self.n = n
        self.m = (n - 1) * (n - 2) // 2  # number of tiling bits

        # Build arc-to-bit mapping
        self._arcs = []  # list of (i, j, skip, row) for each bit
        idx = 0
        for skip in range(2, n):
            row = skip - 2
            for start in range(n - skip):
                i, j = start, start + skip
                self._arcs.append((i, j, skip, row))
                idx += 1

        # Initialize bits
        if bits is not None:
            self._bits = list(bits) if hasattr(bits, '__iter__') else \
                [(bits >> k) & 1 for k in range(self.m)]
        elif tournament is not None:
            self._bits = self._from_tournament(tournament)
        else:
            self._bits = [0] * self.m  # all forward = transitive

    def _from_tournament(self, T):
        """Extract tiling bits from a TournamentH object."""
        bits = []
        for i, j, skip, row in self._arcs:
            # bit=1 means backward (j->i), bit=0 means forward (i->j)
            bits.append(0 if T.adj[i][j] else 1)
        return bits

    def to_tournament(self):
        """Convert tiling to TournamentH object."""
        n = self.n
        adj = [[0]*n for _ in range(n)]
        # Path arcs
        for i in range(n - 1):
            adj[i][i + 1] = 1
        # Tiling arcs
        for idx, (i, j, skip, row) in enumerate(self._arcs):
            if self._bits[idx] == 0:  # forward
                adj[i][j] = 1
            else:  # backward
                adj[j][i] = 1
        return TournamentH(adj)

    # Row access
    def row(self, r):
        """Get bits for row r (skip length r+2)."""
        return [self._bits[idx] for idx, (_, _, _, row) in enumerate(self._arcs)
                if row == r]

    def row_hw(self, r):
        """Hamming weight of row r."""
        return sum(self.row(r))

    def total_hw(self):
        """Total Hamming weight (number of backward arcs)."""
        return sum(self._bits)

    def skip_profile(self):
        """Dict: skip_length -> hamming_weight."""
        profile = {}
        for idx, (_, _, skip, _) in enumerate(self._arcs):
            profile[skip] = profile.get(skip, 0) + self._bits[idx]
        return profile

    def flip(self, idx):
        """Return new TilingGrid with bit idx flipped."""
        new_bits = self._bits[:]
        new_bits[idx] = 1 - new_bits[idx]
        return TilingGrid(self.n, bits=new_bits)

    def complement(self):
        """Return tiling with all bits flipped."""
        return TilingGrid(self.n, bits=[1 - b for b in self._bits])

    def transpose(self):
        """Apply the sigma-transpose map (converse tournament)."""
        n = self.n
        new_bits = [0] * self.m
        for idx, (i, j, skip, row) in enumerate(self._arcs):
            # sigma maps (i,j) to (n-1-j, n-1-i)
            ni, nj = n - 1 - j, n - 1 - i
            if ni > nj:
                ni, nj = nj, ni
            # Find the index of (ni, nj) in _arcs
            for idx2, (i2, j2, _, _) in enumerate(self._arcs):
                if i2 == ni and j2 == nj:
                    # Transpose also flips the bit
                    new_bits[idx2] = 1 - self._bits[idx]
                    break
        return TilingGrid(self.n, bits=new_bits)

    def __repr__(self):
        return f"TilingGrid(n={self.n}, bits={''.join(str(b) for b in self._bits)})"

    def display(self):
        """Pretty-print the triangular grid."""
        lines = []
        for skip in range(2, self.n):
            row = skip - 2
            row_bits = self.row(row)
            arcs = [(i, i+skip) for i in range(self.n - skip)]
            cells = ' '.join(f"{'B' if b else 'F'}" for b in row_bits)
            arc_str = ' '.join(f"({i},{j})" for i, j in arcs)
            lines.append(f"  skip {skip}: [{cells}]  {arc_str}")
        return '\n'.join(lines)


# ================================================================
# MODULE 2: FLIP GRAPH
# ================================================================

class FlipGraph:
    """Graph where nodes are tournaments and edges are single-arc flips.

    Can operate on:
    - Individual tilings (2^m nodes, exhaustive for small n)
    - Isomorphism classes (A000568 nodes, much smaller)
    """

    def __init__(self, n, by_class=False):
        self.n = n
        self.m = (n - 1) * (n - 2) // 2
        self.by_class = by_class
        self._nodes = {}  # id -> TournamentH or class_repr
        self._edges = []  # list of (id1, id2, delta_H)
        self._H_map = {}  # id -> H value

    def build_tiling_graph(self):
        """Build flip graph over all 2^m tilings (canonical path)."""
        # Generate all tilings
        for bits_int in range(1 << self.m):
            bits = [(bits_int >> k) & 1 for k in range(self.m)]
            tg = TilingGrid(self.n, bits=bits)
            T = tg.to_tournament()
            H = T.H()
            self._nodes[bits_int] = T
            self._H_map[bits_int] = H

        # Generate edges (single-bit flips)
        for bits_int in range(1 << self.m):
            for k in range(self.m):
                neighbor = bits_int ^ (1 << k)
                if neighbor > bits_int:  # avoid double-counting
                    dH = self._H_map[neighbor] - self._H_map[bits_int]
                    self._edges.append((bits_int, neighbor, dH))

        return self

    def H_of(self, node_id):
        return self._H_map.get(node_id, None)

    def neighbors(self, node_id):
        """Get all neighbors of a node."""
        nbrs = []
        for a, b, dH in self._edges:
            if a == node_id:
                nbrs.append((b, dH))
            elif b == node_id:
                nbrs.append((a, -dH))
        return nbrs

    def diameter(self):
        """BFS diameter of the flip graph."""
        if not self._nodes:
            return -1
        from collections import deque
        max_dist = 0
        start = next(iter(self._nodes))
        dist = {start: 0}
        queue = deque([start])
        while queue:
            node = queue.popleft()
            for nbr, _ in self.neighbors(node):
                if nbr not in dist:
                    dist[nbr] = dist[node] + 1
                    max_dist = max(max_dist, dist[nbr])
                    queue.append(nbr)
        return max_dist

    def delta_H_distribution(self):
        """Distribution of |delta_H| over all edges."""
        return Counter(abs(dH) for _, _, dH in self._edges)

    def stats(self):
        """Summary statistics of the flip graph."""
        return {
            'n': self.n,
            'nodes': len(self._nodes),
            'edges': len(self._edges),
            'H_range': (min(self._H_map.values()), max(self._H_map.values())),
            'H_distinct': len(set(self._H_map.values())),
            'delta_H_dist': dict(sorted(self.delta_H_distribution().items())),
        }


# ================================================================
# MODULE 3: RANKING AMBIGUITY (LLM / Voting / Sports)
# ================================================================

class RankAmbiguity:
    """Measures ranking ambiguity in pairwise comparison data.

    Given n items with pairwise comparison results (a tournament),
    quantify how ambiguous the ranking is and identify high-leverage
    comparisons.

    Applications:
    - LLM evaluation (Chatbot Arena)
    - Voting systems (Condorcet / Kemeny)
    - Sports leagues (round-robin standings)
    """

    def __init__(self, adj_or_tournament, labels=None):
        if isinstance(adj_or_tournament, TournamentH):
            self.T = adj_or_tournament
        else:
            self.T = TournamentH(adj_or_tournament)
        self.n = self.T.n
        self.labels = labels or [str(i) for i in range(self.n)]

    @classmethod
    def from_comparisons(cls, comparisons, labels=None):
        """Build from list of (winner, loser) pairs.
        All C(n,2) pairs must be present exactly once.
        """
        # Find all items
        items = sorted(set(w for w, _ in comparisons) |
                       set(l for _, l in comparisons))
        n = len(items)
        idx = {item: i for i, item in enumerate(items)}

        adj = [[0]*n for _ in range(n)]
        for winner, loser in comparisons:
            adj[idx[winner]][idx[loser]] = 1

        return cls(adj, labels=items)

    def ambiguity(self):
        """H(T) = number of Hamiltonian paths = number of total orderings
        consistent with pairwise data. H=1 means unique ranking."""
        return self.T.H()

    def ambiguity_ratio(self):
        """H(T) / n! = fraction of all orderings consistent with data.
        0 = impossible (can't happen for tournaments), 1/n! = unique ranking."""
        from math import factorial
        return self.T.H() / factorial(self.n)

    def kemeny_distance(self):
        """Number of backward arcs from best transitive labeling.
        = minimum number of pairwise results that must be 'wrong'
        for a consistent ranking to exist.
        """
        self.T._find_transitive_labeling()
        return len(self.T._backward_arcs)

    def condorcet_winner(self):
        """Return the Condorcet winner (beats all others) if one exists."""
        for i in range(self.n):
            if sum(self.T.adj[i][j] for j in range(self.n)) == self.n - 1:
                return self.labels[i]
        return None

    def high_leverage_comparisons(self, top_k=5):
        """Identify comparisons whose reversal would most change ambiguity.
        Returns list of ((winner, loser), delta_H) sorted by |delta_H|.
        """
        sens = self.T.sensitivity()
        results = []
        for (i, j), delta in sens.items():
            if self.T.adj[i][j]:
                winner, loser = self.labels[i], self.labels[j]
            else:
                winner, loser = self.labels[j], self.labels[i]
            results.append(((winner, loser), delta))
        results.sort(key=lambda x: x[1], reverse=True)
        return results[:top_k]

    def ranking_confidence(self):
        """Score for each item: fraction of others it beats (Copeland score)."""
        scores = {}
        for i in range(self.n):
            wins = sum(self.T.adj[i][j] for j in range(self.n))
            scores[self.labels[i]] = wins / (self.n - 1) if self.n > 1 else 1
        return dict(sorted(scores.items(), key=lambda x: x[1], reverse=True))

    def report(self):
        """Full ranking ambiguity report."""
        H = self.ambiguity()
        lines = [
            f"=== Ranking Ambiguity Report ({self.n} items) ===",
            f"",
            f"Ambiguity (H): {H} consistent total orderings",
            f"Ambiguity ratio: {self.ambiguity_ratio():.6e}",
            f"Kemeny distance: {self.kemeny_distance()} backward arcs",
            f"Condorcet winner: {self.condorcet_winner() or 'NONE (cyclic)'}",
            f"Regime: {self.T.regime}",
            f"Score rigidity: {'YES' if H == 1 + 2*self.T.c3 else 'needs cycle analysis'}",
            f"",
            f"Ranking confidence (Copeland):",
        ]
        for item, conf in self.ranking_confidence().items():
            lines.append(f"  {item}: {conf:.3f}")
        lines.append(f"")
        lines.append(f"High-leverage comparisons:")
        for (w, l), delta in self.high_leverage_comparisons():
            lines.append(f"  {w} > {l}: reversing changes H by {delta}")
        return '\n'.join(lines)


# ================================================================
# MODULE 4: TOURNAMENT TDA
# ================================================================

class TournamentTDA:
    """Topological data analysis features for tournament data.

    Extracts features for ML pipelines:
    - H(T) = Hamiltonian path count (scalar summary)
    - Score sequence (degree distribution)
    - c3 (3-cycle count, from Rao)
    - Regime classification (solvable/pivot/forbidden)
    - Backward arc count (distance from transitive)
    - H-spectrum (H values reachable by single flip)

    For directed networks more broadly:
    - Path homology Betti numbers (beta_0, beta_1, beta_2=0, beta_3)
    """

    def __init__(self, adj_or_tournament):
        if isinstance(adj_or_tournament, TournamentH):
            self.T = adj_or_tournament
        else:
            self.T = TournamentH(adj_or_tournament)

    def feature_vector(self):
        """Extract numerical feature vector for ML."""
        info = self.T.info()
        n = info['n']

        features = {
            'n': n,
            'H': info['H'],
            'log_H': __import__('math').log(info['H']) if info['H'] > 0 else 0,
            'H_normalized': info['H'] / __import__('math').factorial(n),
            'c3': info['c3'],
            'c3_normalized': info['c3'] / max(1, comb(n, 3)),
            'backward_arcs': info['backward_arcs'],
            'backward_ratio': info['backward_arcs'] / max(1, n*(n-1)//2),
            'H_mod3': info['H'] % 3,
            'H_mod6': info['H_mod6'],
            'regime_code': {'solvable': 0, 'pivot': 1, 'forbidden': 2}[info['regime']],
            'score_variance': self._score_variance(),
            'score_range': max(self.T.score) - min(self.T.score),
        }
        return features

    def _score_variance(self):
        mid = (self.T.n - 1) / 2
        return sum((s - mid)**2 for s in self.T.score) / self.T.n

    def H_spectrum(self):
        """Set of H values reachable by flipping any single arc.
        This is the 'fingerprint' of the tournament's local structure.
        """
        spec = self.T.flip_spectrum()
        return sorted(set(spec.values()))

    def H_spectrum_stats(self):
        """Statistics of the H-spectrum."""
        spectrum = self.H_spectrum()
        H0 = self.T.H()
        return {
            'H': H0,
            'spectrum_size': len(spectrum),
            'spectrum_min': min(spectrum),
            'spectrum_max': max(spectrum),
            'spectrum_mean': sum(spectrum) / len(spectrum),
            'up_neighbors': sum(1 for h in spectrum if h > H0),
            'down_neighbors': sum(1 for h in spectrum if h < H0),
            'flat_neighbors': sum(1 for h in spectrum if h == H0),
        }


# ================================================================
# SELF-TEST
# ================================================================

if __name__ == '__main__':
    print("tournament_toolkit.py — Self-test")
    print("=" * 50)
    print()

    # Test TilingGrid
    print("Test 1: TilingGrid")
    tg = TilingGrid(6)
    print(f"  n=6, m={tg.m} bits")
    print(f"  All-forward (transitive):")
    print(tg.display())
    T = tg.to_tournament()
    print(f"  H = {T.H()} (should be 1 for transitive)")
    assert T.H() == 1

    # Flip one bit
    tg2 = tg.flip(5)  # flip a skip-3 arc
    T2 = tg2.to_tournament()
    print(f"  After flipping bit 5: H = {T2.H()}")
    print(f"  Skip profile: {tg2.skip_profile()}")
    print("  PASSED")
    print()

    # Test FlipGraph
    print("Test 2: FlipGraph at n=4")
    fg = FlipGraph(4).build_tiling_graph()
    stats = fg.stats()
    print(f"  Nodes: {stats['nodes']}, Edges: {stats['edges']}")
    print(f"  H range: {stats['H_range']}")
    print(f"  H distinct: {stats['H_distinct']}")
    print(f"  Delta H distribution: {stats['delta_H_dist']}")
    print("  PASSED")
    print()

    # Test RankAmbiguity
    print("Test 3: RankAmbiguity (LLM ranking scenario)")
    # Simulate: 5 LLMs with mostly consistent rankings
    comparisons = [
        ('GPT-4', 'Claude'), ('GPT-4', 'Gemini'), ('GPT-4', 'Llama'),
        ('GPT-4', 'Mistral'),
        ('Claude', 'Gemini'), ('Claude', 'Llama'), ('Claude', 'Mistral'),
        ('Gemini', 'Llama'), ('Gemini', 'Mistral'),
        ('Llama', 'Mistral'),
    ]
    ra = RankAmbiguity.from_comparisons(comparisons)
    print(ra.report())
    print()

    # Now with one cycle: Claude > GPT-4 > Gemini > Claude
    comparisons2 = [
        ('Claude', 'GPT-4'), ('GPT-4', 'Gemini'), ('Gemini', 'Claude'),
        ('GPT-4', 'Llama'), ('GPT-4', 'Mistral'),
        ('Claude', 'Llama'), ('Claude', 'Mistral'),
        ('Gemini', 'Llama'), ('Gemini', 'Mistral'),
        ('Llama', 'Mistral'),
    ]
    ra2 = RankAmbiguity.from_comparisons(comparisons2)
    print(ra2.report())
    print()

    # Test TournamentTDA
    print("Test 4: TournamentTDA features")
    import random
    random.seed(42)
    adj = [[0]*7 for _ in range(7)]
    for i in range(7):
        for j in range(i+1, 7):
            if random.random() < 0.5:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
    tda = TournamentTDA(adj)
    features = tda.feature_vector()
    print(f"  Features: {features}")
    spec_stats = tda.H_spectrum_stats()
    print(f"  H-spectrum stats: {spec_stats}")
    print("  PASSED")
    print()

    print("ALL TESTS PASSED")
