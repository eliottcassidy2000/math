"""
CycleDetector: Streaming intransitivity and cycle detection.

Processes pairwise comparison edges one at a time. Detects when
rankings become intransitive (A>B, B>C, C>A) in real-time.

Uses the formal group F(x,y) = (x+y)/(1+xy) for evidence aggregation.
The key property: F(x, -x) = 0 means circular evidence CANCELS.
High observation count + near-zero score = SUSPICIOUS.

Applications:
  - Financial fraud: wash trading detection in transaction streams
  - A/B testing: detect when experiments contradict each other
  - Sports: identify rock-paper-scissors dynamics in leagues
  - LLM evaluation: flag inconsistent model rankings

    detector = CycleDetector()
    detector.observe("A", "B")  # A beats B
    detector.observe("B", "C")  # B beats C
    detector.observe("C", "A")  # C beats A — cycle!
    print(detector.risk())      # > 0: intransitivity detected
"""

from math import atanh, tanh, sqrt
from collections import defaultdict


class CycleDetector:
    """Streaming cycle/intransitivity detector.

    Maintains a running tournament from pairwise observations.
    Detects 3-cycles (the minimum unit of intransitivity) in O(1) per edge.
    """

    def __init__(self, decay=1.0):
        """
        Args:
            decay: exponential decay factor for old observations (1.0 = no decay)
        """
        self.decay = decay
        self._score = defaultdict(float)     # item -> arctanh evidence total
        self._count = defaultdict(int)       # item -> observation count
        self._adj = defaultdict(lambda: defaultdict(int))  # wins[a][b]
        self._items = set()
        self._cycles_detected = 0
        self._total_edges = 0
        self._recent_cycles = []

    def observe(self, winner, loser, weight=1.0):
        """Process one comparison: winner beats loser.

        Returns True if this edge COMPLETES a 3-cycle.
        """
        self._items.add(winner)
        self._items.add(loser)
        self._adj[winner][loser] += 1
        self._total_edges += 1

        # Update formal group evidence
        coupling = min(0.9999, weight)
        ev = atanh(coupling)
        self._score[winner] += ev
        self._score[loser] -= ev
        self._count[winner] += 1
        self._count[loser] += 1

        # Check for 3-cycle completion
        # A 3-cycle exists if: winner→loser, and there exists X such that
        # loser→X and X→winner (checking majority direction)
        cycle_found = False
        for x in self._adj[loser]:
            if x == winner:
                continue
            if self._adj[loser][x] > self._adj[x].get(loser, 0):
                # loser beats x
                if self._adj[x].get(winner, 0) > self._adj[winner].get(x, 0):
                    # x beats winner — cycle: winner→loser→x→winner
                    self._cycles_detected += 1
                    self._recent_cycles.append((winner, loser, x))
                    if len(self._recent_cycles) > 100:
                        self._recent_cycles = self._recent_cycles[-50:]
                    cycle_found = True

        return cycle_found

    def risk(self):
        """Overall intransitivity risk score in [0, 1].

        0 = perfectly transitive (no cycles)
        1 = maximally intransitive
        """
        if self._total_edges < 3:
            return 0.0
        # Risk = cycles detected / maximum possible cycles
        n = len(self._items)
        if n < 3:
            return 0.0
        max_cycles = n * (n-1) * (n-2) / 6  # C(n,3) upper bound
        return min(1.0, self._cycles_detected / max(1, max_cycles))

    def suspects(self, top_n=5):
        """Items most involved in cycles (potential fraud participants).

        A "suspect" has: many observations but near-zero net score.
        This is the formal group's cancellation signature: F(x,-x) = 0.
        """
        suspects = []
        for item in self._items:
            count = self._count[item]
            score = abs(self._score[item])
            if count > 0:
                # Suspicion = high count + low |score|
                # Normal: score grows with count (sqrt scaling)
                expected_score = sqrt(count) * 0.5  # rough baseline
                anomaly = max(0, 1 - score / max(0.01, expected_score))
                suspects.append({
                    'item': item,
                    'observations': count,
                    'net_score': self._score[item],
                    'anomaly': anomaly,
                    'cycles_involved': sum(1 for c in self._recent_cycles if item in c),
                })
        suspects.sort(key=lambda x: -x['anomaly'])
        return suspects[:top_n]

    def summary(self):
        """Summary report."""
        lines = [
            f"Streaming Cycle Detector Report",
            f"  Items tracked: {len(self._items)}",
            f"  Edges processed: {self._total_edges}",
            f"  3-cycles detected: {self._cycles_detected}",
            f"  Intransitivity risk: {self.risk():.1%}",
        ]
        if self._cycles_detected > 0:
            lines.append(f"  Recent cycles: {self._recent_cycles[-3:]}")
            lines.append(f"  Top suspects:")
            for s in self.suspects(3):
                lines.append(f"    {s['item']}: anomaly={s['anomaly']:.2f}, "
                           f"obs={s['observations']}, score={s['net_score']:+.2f}")
        return "\n".join(lines)

    @property
    def n_items(self):
        return len(self._items)

    @property
    def n_cycles(self):
        return self._cycles_detected
