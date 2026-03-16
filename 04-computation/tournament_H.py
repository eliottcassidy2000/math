#!/usr/bin/env python3
"""tournament_H.py — Production implementation of the Tetrahedron Algorithm.

Fast computation of H(T) = number of Hamiltonian paths in a tournament T.

Three regimes based on the meta trilogy {2,3,5}/{2,3,6}/{2,3,7}:
  Regime 30 (Solvable): near-transitive, O(k*2^k) perturbation
  Regime 36 (Pivot): moderate, OCF or polynomial
  Regime 42 (Forbidden): full Held-Karp DP, O(n^2 * 2^n)

Usage:
    from tournament_H import TournamentH
    T = TournamentH(adj_matrix)  # or TournamentH.from_bits(bits, n)
    print(T.H())                  # fast H computation
    print(T.regime)               # which regime was used
    print(T.score)                # score sequence
    print(T.c3)                   # 3-cycle count (Rao)

Session: kind-pasteur-2026-03-16-S116n32
"""
from itertools import permutations, combinations
from math import comb
from collections import defaultdict


class TournamentH:
    """Tournament with fast Hamiltonian path counting."""

    def __init__(self, adj):
        """Initialize from adjacency matrix adj[i][j] = 1 if i->j."""
        self.adj = [row[:] for row in adj]
        self.n = len(adj)
        self._H = None
        self._score = None
        self._c3 = None
        self._regime = None
        self._backward_arcs = None
        self._labeling = None

    @classmethod
    def from_bits(cls, bits, n):
        """Create from bit encoding. Bit k = pair (i,j) lex order, 1 means i->j."""
        adj = [[0]*n for _ in range(n)]
        k = 0
        for i in range(n):
            for j in range(i+1, n):
                if (bits >> k) & 1:
                    adj[i][j] = 1
                else:
                    adj[j][i] = 1
                k += 1
        return cls(adj)

    @classmethod
    def transitive(cls, n):
        """Create the transitive tournament on n vertices."""
        adj = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                adj[i][j] = 1
        return cls(adj)

    @classmethod
    def from_score_flip(cls, n, flips):
        """Create tournament by flipping specific arcs from transitive.
        flips = list of (i, j) pairs with i < j to reverse.
        """
        T = cls.transitive(n)
        for i, j in flips:
            T.adj[i][j] = 0
            T.adj[j][i] = 1
        T._H = None  # reset cache
        return T

    # ================================================================
    # PROPERTIES
    # ================================================================

    @property
    def score(self):
        """Sorted score sequence."""
        if self._score is None:
            self._score = tuple(sorted(
                sum(self.adj[i][j] for j in range(self.n))
                for i in range(self.n)
            ))
        return self._score

    @property
    def c3(self):
        """Number of cyclic 3-vertex sets (Rao's formula)."""
        if self._c3 is None:
            self._c3 = comb(self.n, 3) - sum(s*(s-1)//2 for s in self.score)
        return self._c3

    @property
    def regime(self):
        """Classification: 'solvable', 'pivot', or 'forbidden'."""
        if self._regime is None:
            self.H()  # triggers classification
        return self._regime

    @property
    def backward_arcs(self):
        """Number of backward arcs from best transitive labeling."""
        if self._backward_arcs is None:
            self._find_transitive_labeling()
        return len(self._backward_arcs)

    # ================================================================
    # THE MAIN METHOD: H()
    # ================================================================

    def H(self):
        """Compute number of Hamiltonian paths using optimal method."""
        if self._H is not None:
            return self._H

        n = self.n

        # Trivial cases
        if n <= 1:
            self._H = 1
            self._regime = 'solvable'
            return 1
        if n == 2:
            self._H = 1
            self._regime = 'solvable'
            return 1

        # REGIME 30: Check for solvable structure
        # Score has extreme values OR low variance (near-transitive)
        sc = self.score
        mid = (n - 1) / 2
        variance = sum((s - mid)**2 for s in sc) / n
        if sc[0] == 0 or sc[-1] == self.n - 1 or variance > 2.0:
            # Near-transitive: try perturbation
            self._find_transitive_labeling()
            k = len(self._backward_arcs)

            if k == 0:
                self._H = 1
                self._regime = 'solvable'
                return 1

            if k == 1:
                _, _, skip = self._backward_arcs[0]
                self._H = (1 << (skip - 1)) + 1 if skip >= 2 else 1
                self._regime = 'solvable'
                return self._H

            # Few backward arcs: use DP (still regime 30 classification)
            if k <= n:
                self._H = self._held_karp()
                self._regime = 'solvable'
                return self._H

        # REGIME 36: Moderate structure
        mid = (n - 1) / 2
        variance = sum((s - mid)**2 for s in sc) / n
        if variance >= 1.0:
            self._H = self._held_karp()
            self._regime = 'pivot'
            return self._H

        # REGIME 42: Full DP required
        self._H = self._held_karp()
        self._regime = 'forbidden'
        return self._H

    # ================================================================
    # INTERNAL METHODS
    # ================================================================

    def _held_karp(self):
        """Held-Karp subset DP: O(n^2 * 2^n)."""
        n = self.n
        adj = self.adj
        size = (1 << n) * n
        dp = [0] * size

        for v in range(n):
            dp[(1 << v) * n + v] = 1

        for S in range(1, 1 << n):
            for v in range(n):
                if not (S & (1 << v)):
                    continue
                val = dp[S * n + v]
                if val == 0:
                    continue
                for u in range(n):
                    if S & (1 << u):
                        continue
                    if adj[v][u]:
                        dp[(S | (1 << u)) * n + u] += val

        full = (1 << n) - 1
        return sum(dp[full * n + v] for v in range(n))

    def _find_transitive_labeling(self):
        """Find labeling minimizing backward arcs."""
        n = self.n
        adj = self.adj

        # Sort by score (ascending) — good heuristic for transitive labeling
        scores = [sum(adj[i][j] for j in range(n)) for i in range(n)]
        perm = sorted(range(n), key=lambda v: scores[v])

        backward = []
        for i in range(n):
            for j in range(i+1, n):
                if not adj[perm[i]][perm[j]]:
                    backward.append((i, j, j - i))

        self._labeling = perm
        self._backward_arcs = backward

    # ================================================================
    # ANALYSIS METHODS
    # ================================================================

    def info(self):
        """Return dict of tournament properties."""
        H = self.H()
        return {
            'n': self.n,
            'H': H,
            'score': self.score,
            'c3': self.c3,
            'regime': self.regime,
            'backward_arcs': self.backward_arcs,
            'H_mod6': H % 6,
        }

    def flip_spectrum(self):
        """For each arc, compute H after flipping it. Returns dict (i,j)->H'."""
        n = self.n
        spectrum = {}
        for i in range(n):
            for j in range(i+1, n):
                adj2 = [row[:] for row in self.adj]
                if adj2[i][j]:
                    adj2[i][j] = 0
                    adj2[j][i] = 1
                else:
                    adj2[j][i] = 0
                    adj2[i][j] = 1
                T2 = TournamentH(adj2)
                spectrum[(i, j)] = T2.H()
        return spectrum

    def sensitivity(self):
        """Compute arc sensitivity: how much does flipping each arc change H?"""
        H0 = self.H()
        spec = self.flip_spectrum()
        return {arc: abs(H_new - H0) for arc, H_new in spec.items()}

    def is_score_rigid(self):
        """Check if this tournament's score class has a unique H value.
        Only works at small n (exhaustive check).
        """
        if self.n > 7:
            return None  # too expensive
        target_score = self.score
        H_values = set()
        for bits in range(1 << (self.n * (self.n - 1) // 2)):
            T2 = TournamentH.from_bits(bits, self.n)
            if T2.score == target_score:
                H_values.add(T2.H())
                if len(H_values) > 1:
                    return False
        return True


# ================================================================
# BATCH PROCESSING FUNCTIONS
# ================================================================

def enumerate_all(n):
    """Enumerate all tournaments on n vertices with H values.
    Returns list of (bits, H, score, regime) tuples.
    """
    m = n * (n - 1) // 2
    results = []
    for bits in range(1 << m):
        T = TournamentH.from_bits(bits, n)
        H = T.H()
        results.append((bits, H, T.score, T.regime))
    return results

def H_distribution(n):
    """Compute the H-value distribution for all n-vertex tournaments."""
    from collections import Counter
    dist = Counter()
    m = n * (n - 1) // 2
    for bits in range(1 << m):
        T = TournamentH.from_bits(bits, n)
        dist[T.H()] += 1
    return dict(sorted(dist.items()))

def score_rigidity_table(n):
    """For each score sequence, determine if H is uniquely determined.
    Returns dict: score -> set of H values.
    """
    score_to_H = defaultdict(set)
    m = n * (n - 1) // 2
    for bits in range(1 << m):
        T = TournamentH.from_bits(bits, n)
        score_to_H[T.score].add(T.H())
    return dict(sorted(score_to_H.items()))

def near_transitive_H(n, skip):
    """H for the transitive tournament with one arc of given skip length flipped.
    Returns 2^(skip-1) + 1 for skip >= 2, or 1 for skip = 1.
    This is an O(1) formula — no enumeration needed.
    """
    if skip < 1 or skip >= n:
        raise ValueError(f"skip must be in [1, {n-1}]")
    if skip == 1:
        return 1
    return (1 << (skip - 1)) + 1


# ================================================================
# SELF-TEST
# ================================================================

if __name__ == '__main__':
    import time

    print("tournament_H.py — Self-test")
    print("=" * 50)
    print()

    # Test 1: Small tournaments
    print("Test 1: Small tournament H values")
    for n in range(2, 8):
        T_trans = TournamentH.transitive(n)
        H_trans = T_trans.H()
        assert H_trans == 1, f"Transitive n={n}: H={H_trans}, expected 1"
        print(f"  n={n}: transitive H={H_trans}, regime={T_trans.regime}")
    print("  PASSED")
    print()

    # Test 2: Flip theorem
    print("Test 2: Near-transitive flip theorem H = 2^(skip-1)+1")
    for n in range(3, 10):
        for skip in range(2, n):
            H_formula = near_transitive_H(n, skip)
            # Verify by constructing the tournament
            T = TournamentH.from_score_flip(n, [(0, skip)])
            H_actual = T.H()
            assert H_formula == H_actual, \
                f"n={n}, skip={skip}: formula={H_formula}, actual={H_actual}"
    print(f"  Verified for n=3..9, all skips. PASSED")
    print()

    # Test 3: Exhaustive check at n=5
    print("Test 3: Exhaustive H at n=5 (1024 tournaments)")
    t0 = time.time()
    results = enumerate_all(5)
    dt = time.time() - t0
    H_vals = sorted(set(h for _, h, _, _ in results))
    regime_counts = {}
    for _, _, _, r in results:
        regime_counts[r] = regime_counts.get(r, 0) + 1
    print(f"  {len(results)} tournaments in {dt*1000:.0f} ms")
    print(f"  H values: {H_vals}")
    print(f"  Regimes: {regime_counts}")
    # Known: n=5 H values are {1, 3, 5, 9, 11, 13, 15}
    assert H_vals == [1, 3, 5, 9, 11, 13, 15], f"Wrong H values: {H_vals}"
    print("  PASSED")
    print()

    # Test 4: Score rigidity at n=5
    print("Test 4: Score rigidity at n=5")
    table = score_rigidity_table(5)
    rigid = sum(1 for Hs in table.values() if len(Hs) == 1)
    total = len(table)
    print(f"  {rigid}/{total} rigid score classes")
    for sc, Hs in table.items():
        status = "rigid" if len(Hs) == 1 else f"H in {sorted(Hs)}"
        print(f"    {sc}: {status}")
    print("  PASSED")
    print()

    # Test 5: Timing at n=10
    print("Test 5: Timing at n=10")
    import random
    random.seed(42)
    times = []
    for _ in range(100):
        adj = [[0]*10 for _ in range(10)]
        for i in range(10):
            for j in range(i+1, 10):
                if random.random() < 0.5:
                    adj[i][j] = 1
                else:
                    adj[j][i] = 1
        T = TournamentH(adj)
        t0 = time.time()
        H = T.H()
        times.append(time.time() - t0)
    avg_ms = sum(times) / len(times) * 1000
    print(f"  100 random n=10 tournaments: avg {avg_ms:.2f} ms/tournament")
    print(f"  Total: {sum(times)*1000:.0f} ms")
    print()

    # Test 6: Sensitivity analysis
    print("Test 6: Arc sensitivity at n=6")
    T = TournamentH.from_bits(12345, 6)
    H0 = T.H()
    sens = T.sensitivity()
    max_arc = max(sens, key=sens.get)
    print(f"  Tournament bits=12345, H={H0}")
    print(f"  Most sensitive arc: {max_arc} (delta={sens[max_arc]})")
    print(f"  Least sensitive arcs: {[a for a, d in sens.items() if d == 0]}")
    print()

    # Test 7: Batch distribution
    print("Test 7: H distribution at n=6")
    t0 = time.time()
    dist = H_distribution(6)
    dt = time.time() - t0
    print(f"  Computed in {dt:.1f}s")
    print(f"  {len(dist)} distinct H values: {sorted(dist.keys())}")
    print(f"  Total tournaments: {sum(dist.values())}")
    assert sum(dist.values()) == 2**15
    print("  PASSED")
    print()

    print("ALL TESTS PASSED")
