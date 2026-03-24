#!/usr/bin/env python3
"""
tournament_predict.py -- Tournament-Based Sequence Predictor & Anomaly Detector
kind-pasteur-2026-03-24-S20cq

THE IDEA: Model a time series as a tournament, where each time step "beats"
future steps when the value increases. The tournament structure (3-cycle count,
score sequence, independence polynomial) captures temporal patterns that
simple statistics miss.

FEATURES:
  1. Sliding-window tournament construction from any numeric sequence
  2. Tournament-derived features: score entropy, 3-cycle density, transitivity
  3. Anomaly detection: sudden changes in tournament structure = anomaly
  4. Prediction: tournament structure determines likely next-value direction
  5. Regime detection: tournament topology changes = regime change

APPLICATIONS:
  - Financial time series (detect regime changes in stock/crypto)
  - Server monitoring (detect anomalies in latency/CPU/memory)
  - Sensor data (detect equipment degradation patterns)
  - Sports analytics (detect momentum shifts)

USAGE:
  from tournament_predict import TournamentPredictor
  tp = TournamentPredictor(window=16)
  tp.fit(time_series)
  anomalies = tp.detect_anomalies(threshold=2.0)
  prediction = tp.predict_direction()

LICENSE: MIT
"""

import math
import sys
from collections import deque
from dataclasses import dataclass, field
from typing import List, Tuple, Optional, Sequence

__version__ = "1.0.0"


@dataclass
class TournamentFeatures:
    """Features extracted from a tournament on a time window."""
    # Basic
    n: int = 0
    score_seq: list = field(default_factory=list)

    # Derived
    score_entropy: float = 0.0      # entropy of score sequence
    score_variance: float = 0.0     # variance of scores
    transitivity: float = 0.0       # fraction of transitive triples
    c3_density: float = 0.0         # 3-cycle density (among all triples)
    dom_ratio: float = 0.0          # fraction dominated by leader
    landau_deficit: float = 0.0     # sum of score deficits from max possible
    upset_rate: float = 0.0         # fraction of "upsets" (later beats earlier)
    momentum: float = 0.0           # recent trend strength

    def to_vector(self) -> list:
        return [
            self.score_entropy,
            self.score_variance,
            self.transitivity,
            self.c3_density,
            self.dom_ratio,
            self.landau_deficit,
            self.upset_rate,
            self.momentum,
        ]


def _build_tournament(values: Sequence[float]) -> list:
    """Build tournament adjacency from time series values.

    Rule: i beats j if values[i] > values[j], or if equal and i < j (tie-break).
    Returns list of lists: adj[i] = set of vertices that i beats.
    """
    n = len(values)
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if values[i] > values[j]:
                adj[i].append(j)
            elif values[i] < values[j]:
                adj[j].append(i)
            else:
                # Tie-break: earlier beats later
                adj[i].append(j)
    return adj


def _scores(adj: list, n: int) -> list:
    """Score sequence (out-degrees), sorted."""
    return sorted(len(adj[i]) for i in range(n))


def _count_c3(adj: list, n: int) -> int:
    """Count directed 3-cycles."""
    # Use score-based formula: c3 = C(n,3) - sum(C(s_i, 2))
    total_triples = n * (n - 1) * (n - 2) // 6
    score_contribution = sum(len(adj[i]) * (len(adj[i]) - 1) // 2 for i in range(n))
    return total_triples - score_contribution


def _entropy(values: list) -> float:
    """Shannon entropy of a distribution."""
    if not values: return 0.0
    total = sum(values)
    if total == 0: return 0.0
    return -sum(v/total * math.log2(v/total) for v in values if v > 0)


def extract_features(values: Sequence[float]) -> TournamentFeatures:
    """Extract tournament features from a sequence of values."""
    n = len(values)
    if n < 3:
        return TournamentFeatures(n=n)

    adj = _build_tournament(values)
    scores = _scores(adj, n)
    c3 = _count_c3(adj, n)
    total_triples = n * (n - 1) * (n - 2) // 6

    feat = TournamentFeatures()
    feat.n = n
    feat.score_seq = scores

    # Score entropy (normalized)
    if sum(scores) > 0:
        feat.score_entropy = _entropy(scores) / math.log2(n) if n > 1 else 0

    # Score variance
    mean_s = sum(scores) / n
    feat.score_variance = sum((s - mean_s)**2 for s in scores) / n

    # Transitivity
    feat.transitivity = 1 - c3 / total_triples if total_triples > 0 else 1.0
    feat.c3_density = c3 / total_triples if total_triples > 0 else 0.0

    # Domination ratio (max score / (n-1))
    feat.dom_ratio = max(scores) / (n - 1) if n > 1 else 0

    # Landau deficit
    max_possible = sum(range(n))  # maximum total score = sum(0..n-1) = n(n-1)/2
    feat.landau_deficit = 1.0 - sum(abs(s - (n-1)/2) for s in scores) / (max_possible / 2 + 0.01)

    # Upset rate: among pairs (i, j) with i < j, how often does i beat j?
    # "Upset" = earlier time step has HIGHER value (decreasing direction).
    # In a monotone increasing series: 0 upsets. Random = 50%.
    upsets = 0
    for i in range(n):
        for j in adj[i]:
            if j > i:
                upsets += 1  # i < j but i beats j => value decreased
    feat.upset_rate = upsets / (n * (n - 1) / 2) if n > 1 else 0

    # Momentum: score of the last element vs average
    last_score = len(adj[n - 1])
    feat.momentum = (last_score - (n - 1) / 2) / ((n - 1) / 2) if n > 1 else 0

    return feat


class TournamentPredictor:
    """Sliding-window tournament predictor and anomaly detector.

    Usage:
        tp = TournamentPredictor(window=16)
        for value in stream:
            tp.update(value)
            if tp.is_anomaly():
                print(f"ANOMALY at {tp.position}")
            direction = tp.predict_direction()
    """

    def __init__(self, window: int = 16, history: int = 100):
        """
        Args:
            window: size of tournament window
            history: number of feature snapshots to keep for anomaly detection
        """
        self.window = window
        self.history = history
        self._buffer = deque(maxlen=window)
        self._feature_history = deque(maxlen=history)
        self._value_history = deque(maxlen=history * window)
        self.position = 0
        self._baseline_mean = None
        self._baseline_std = None

    def update(self, value: float):
        """Add a new value to the stream."""
        self._buffer.append(value)
        self._value_history.append(value)
        self.position += 1

        if len(self._buffer) >= self.window:
            feat = extract_features(list(self._buffer))
            self._feature_history.append(feat)
            self._update_baseline()

    def _update_baseline(self):
        """Update rolling baseline statistics."""
        if len(self._feature_history) < 5:
            return
        vectors = [f.to_vector() for f in self._feature_history]
        n = len(vectors)
        dim = len(vectors[0])
        self._baseline_mean = [sum(v[d] for v in vectors) / n for d in range(dim)]
        self._baseline_std = [
            math.sqrt(sum((v[d] - self._baseline_mean[d])**2 for v in vectors) / n + 1e-10)
            for d in range(dim)
        ]

    def current_features(self) -> Optional[TournamentFeatures]:
        """Get features of current window."""
        if self._feature_history:
            return self._feature_history[-1]
        return None

    def anomaly_score(self) -> float:
        """Z-score based anomaly score of current window.

        Returns: float >= 0. Higher = more anomalous. >2.0 is suspicious, >3.0 is anomaly.
        """
        if self._baseline_mean is None or not self._feature_history:
            return 0.0

        current = self._feature_history[-1].to_vector()
        z_scores = [
            abs(current[d] - self._baseline_mean[d]) / self._baseline_std[d]
            for d in range(len(current))
        ]
        return max(z_scores)  # max z-score across all features

    def is_anomaly(self, threshold: float = 3.0) -> bool:
        """Check if current window is anomalous."""
        return self.anomaly_score() > threshold

    def predict_direction(self) -> float:
        """Predict likely direction of next value.

        Returns: float in [-1, 1]. Positive = likely increase, negative = likely decrease.
        Based on tournament momentum and transitivity.
        """
        if not self._feature_history:
            return 0.0

        feat = self._feature_history[-1]
        # Momentum from tournament structure
        direction = feat.momentum * 0.6

        # High transitivity = strong trend, amplify direction
        if feat.transitivity > 0.8:
            direction *= 1.5
        # High c3 density = choppy/mean-reverting, dampen direction
        elif feat.c3_density > 0.3:
            direction *= 0.3

        # Upset rate correction: many upsets = reversal likely
        if feat.upset_rate > 0.6:
            direction *= -0.5

        return max(-1.0, min(1.0, direction))

    def regime_change_score(self) -> float:
        """Detect regime changes via tournament topology shift.

        Returns: float >= 0. Higher = more likely regime change.
        """
        if len(self._feature_history) < 10:
            return 0.0

        # Compare recent features vs older features
        recent = [f.to_vector() for f in list(self._feature_history)[-5:]]
        older = [f.to_vector() for f in list(self._feature_history)[-10:-5]]

        dim = len(recent[0])
        recent_mean = [sum(v[d] for v in recent) / 5 for d in range(dim)]
        older_mean = [sum(v[d] for v in older) / 5 for d in range(dim)]
        older_std = [
            math.sqrt(sum((v[d] - older_mean[d])**2 for v in older) / 5 + 1e-10)
            for d in range(dim)
        ]

        shifts = [abs(recent_mean[d] - older_mean[d]) / older_std[d] for d in range(dim)]
        return max(shifts)

    def fit(self, series: Sequence[float]) -> 'TournamentPredictor':
        """Fit on a full series (batch mode)."""
        for value in series:
            self.update(value)
        return self

    def detect_anomalies(self, series: Optional[Sequence[float]] = None,
                         threshold: float = 3.0) -> List[Tuple[int, float]]:
        """Detect all anomalies in a series.

        Returns: list of (position, score) tuples for anomalous points.
        """
        if series is not None:
            self.__init__(self.window, self.history)
            self.fit(series)

        anomalies = []
        # Re-scan with anomaly detection
        tp = TournamentPredictor(self.window, self.history)
        for i, v in enumerate(self._value_history):
            tp.update(v)
            score = tp.anomaly_score()
            if score > threshold:
                anomalies.append((i, score))
        return anomalies

    def summary(self) -> str:
        """Summary of current state."""
        feat = self.current_features()
        if feat is None:
            return f"TournamentPredictor(window={self.window}, buffered={len(self._buffer)})"

        lines = [
            f"TournamentPredictor v{__version__}",
            f"  Window:        {self.window}",
            f"  Position:      {self.position}",
            f"  Score entropy: {feat.score_entropy:.3f}",
            f"  Transitivity:  {feat.transitivity:.3f}",
            f"  C3 density:    {feat.c3_density:.3f}",
            f"  Upset rate:    {feat.upset_rate:.3f}",
            f"  Momentum:      {feat.momentum:+.3f}",
            f"  Direction:     {self.predict_direction():+.3f}",
            f"  Anomaly score: {self.anomaly_score():.2f}",
            f"  Regime shift:  {self.regime_change_score():.2f}",
        ]
        return "\n".join(lines)


# ============================================================================
# DEMO
# ============================================================================

def demo():
    """Demonstrate tournament prediction on synthetic data."""
    import random
    random.seed(42)

    print(f"tournament_predict v{__version__} -- Demo")
    print("=" * 70)

    # 1. Trending series
    print("\n1. TRENDING SERIES (should predict positive direction)")
    series = [10 + i * 0.5 + random.gauss(0, 1) for i in range(100)]
    tp = TournamentPredictor(window=16)
    tp.fit(series)
    print(tp.summary())

    # 2. Mean-reverting series
    print("\n2. MEAN-REVERTING SERIES (should show high c3, low momentum)")
    series = [50 + 10 * math.sin(i / 5) + random.gauss(0, 2) for i in range(100)]
    tp = TournamentPredictor(window=16)
    tp.fit(series)
    print(tp.summary())

    # 3. Regime change
    print("\n3. REGIME CHANGE (two phases with different dynamics)")
    phase1 = [50 + random.gauss(0, 1) for _ in range(50)]
    phase2 = [80 + random.gauss(0, 5) for _ in range(50)]
    series = phase1 + phase2
    tp = TournamentPredictor(window=16)
    tp.fit(series)
    print(tp.summary())
    print(f"  Regime change score: {tp.regime_change_score():.2f}")

    # 4. Random walk
    print("\n4. RANDOM WALK (no predictable direction)")
    series = [0.0]
    for _ in range(99):
        series.append(series[-1] + random.gauss(0, 1))
    tp = TournamentPredictor(window=16)
    tp.fit(series)
    print(tp.summary())

    # 5. Anomaly detection
    print("\n5. ANOMALY DETECTION")
    normal = [50 + random.gauss(0, 2) for _ in range(80)]
    anomalous = [50 + random.gauss(0, 2) for _ in range(10)]
    anomalous[5] = 200  # spike!
    anomalous[7] = -100  # crash!
    series = normal + anomalous + [50 + random.gauss(0, 2) for _ in range(10)]

    tp = TournamentPredictor(window=16, history=50)
    anomalies = tp.detect_anomalies(series, threshold=2.0)
    print(f"  Found {len(anomalies)} anomalies:")
    for pos, score in anomalies[:10]:
        print(f"    Position {pos}: score={score:.2f} value={series[pos]:.1f}")

    # 6. Streaming demo
    print("\n6. STREAMING MODE")
    tp = TournamentPredictor(window=8)
    print(f"  {'Step':>4} {'Value':>7} {'Dir':>6} {'Anom':>6} {'Trans':>6} {'C3':>6}")
    for i in range(30):
        v = 50 + 10 * math.sin(i / 3) + random.gauss(0, 1)
        tp.update(v)
        feat = tp.current_features()
        if feat and feat.n >= 3:
            print(f"  {i:4d} {v:7.1f} {tp.predict_direction():+5.2f} "
                  f"{tp.anomaly_score():5.2f} {feat.transitivity:5.3f} {feat.c3_density:5.3f}")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=f'tournament_predict v{__version__}')
    parser.add_argument('--demo', action='store_true', help='Run demo')
    parser.add_argument('--file', '-f', help='CSV file with numeric values (one per line)')
    parser.add_argument('--window', '-w', type=int, default=16, help='Window size')
    parser.add_argument('--threshold', '-t', type=float, default=3.0, help='Anomaly threshold')
    args = parser.parse_args()

    if args.demo:
        demo()
    elif args.file:
        with open(args.file) as f:
            values = [float(line.strip()) for line in f if line.strip()]
        tp = TournamentPredictor(window=args.window)
        tp.fit(values)
        print(tp.summary())
        anomalies = tp.detect_anomalies(values, threshold=args.threshold)
        if anomalies:
            print(f"\n  Anomalies ({len(anomalies)}):")
            for pos, score in anomalies:
                print(f"    [{pos}] value={values[pos]:.2f}, score={score:.2f}")
    else:
        parser.print_help()
