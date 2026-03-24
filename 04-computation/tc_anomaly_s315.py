#!/usr/bin/env python3
"""
tc_anomaly_s315.py — Tournament Anomaly Detector
opus-2026-03-24-S315

INSIGHT: Compressibility IS normalcy. If a data chunk compresses well
under the learned model, it's "normal." If it doesn't compress well,
it's "anomalous" — it deviates from the expected score distribution.

This turns the tournament codec into an ANOMALY DETECTOR:
  - Train: stream normal data, learn score distributions
  - Detect: for each new chunk, compute compression ratio
  - Alert: if ratio drops below threshold → anomaly

APPLICATIONS:
  - Network intrusion detection (packet structure changes)
  - Sensor fault detection (reading pattern changes)
  - Manufacturing quality control (signal deviates from norm)
  - Financial fraud (transaction pattern changes)

ALSO: A tournament-based DATA SKETCHING tool that compresses
a dataset into a compact "sketch" (score distribution) that
can answer approximate queries.
"""

import sys, os, struct, zlib, time
import numpy as np
from math import comb, ceil, log2
from collections import deque, Counter

sys.stdout.reconfigure(line_buffering=True)

ROW_WIDTH = 32

class TournamentAnomalyDetector:
    """Detect anomalies via compressibility shift."""

    def __init__(self, window_size=50, threshold_sigma=3.0):
        self.window_size = window_size
        self.threshold_sigma = threshold_sigma
        self.score_model = np.ones(ROW_WIDTH + 1, dtype=np.float64)
        self.score_total = ROW_WIDTH + 1
        self.ratio_history = deque(maxlen=window_size)
        self.baseline_mean = None
        self.baseline_std = None
        self.trained = False
        self.chunk_count = 0

    def _to_matrix(self, data):
        n = len(data)
        arr = np.frombuffer(data, dtype=np.uint8)
        gray = arr ^ (arr >> 1)
        H = (n + ROW_WIDTH - 1) // ROW_WIDTH
        padded = np.zeros(H * ROW_WIDTH, dtype=np.uint8)
        padded[:n] = gray
        return padded.reshape(H, ROW_WIDTH)

    def _score_cost(self, matrix):
        """Compute score-conditioned encoding cost."""
        H, W = matrix.shape
        total = 0.0
        for b in range(8):
            plane = ((matrix >> (7 - b)) & 1).astype(np.uint8)
            # Delta-row
            delta = np.zeros_like(plane)
            delta[0] = plane[0]
            for i in range(1, H): delta[i] = plane[i] ^ plane[i - 1]
            for i in range(H):
                weight = int(np.sum(delta[i]))
                # Score entropy under learned model
                p = self.score_model[weight] / self.score_total
                total += -log2(max(p, 1e-10))
                # Combinadic cost
                if 0 < weight < W:
                    total += log2(max(comb(W, weight), 1))
        return total

    def _update_model(self, matrix):
        """Update score distribution from observed data."""
        H, W = matrix.shape
        for b in range(8):
            plane = ((matrix >> (7 - b)) & 1).astype(np.uint8)
            delta = np.zeros_like(plane)
            delta[0] = plane[0]
            for i in range(1, H): delta[i] = plane[i] ^ plane[i - 1]
            for i in range(H):
                weight = int(np.sum(delta[i]))
                self.score_model[weight] += 1
                self.score_total += 1

    def train(self, data_chunks):
        """Train on a sequence of normal data chunks."""
        for chunk in data_chunks:
            matrix = self._to_matrix(chunk)
            cost = self._score_cost(matrix)
            raw = len(chunk) * 8
            ratio = raw / max(cost, 1)
            self.ratio_history.append(ratio)
            self._update_model(matrix)
            self.chunk_count += 1

        ratios = list(self.ratio_history)
        self.baseline_mean = np.mean(ratios)
        self.baseline_std = max(np.std(ratios), 0.01)
        self.trained = True

    def detect(self, chunk):
        """Check if a chunk is anomalous. Returns (is_anomaly, score, details)."""
        matrix = self._to_matrix(chunk)
        cost = self._score_cost(matrix)
        raw = len(chunk) * 8
        ratio = raw / max(cost, 1)

        self._update_model(matrix)  # keep learning
        self.ratio_history.append(ratio)
        self.chunk_count += 1

        if not self.trained:
            return False, 0.0, {'ratio': ratio, 'message': 'not trained yet'}

        # Z-score: how many standard deviations from baseline?
        z = (ratio - self.baseline_mean) / self.baseline_std
        is_anomaly = abs(z) > self.threshold_sigma

        return is_anomaly, z, {
            'ratio': ratio,
            'baseline_mean': self.baseline_mean,
            'baseline_std': self.baseline_std,
            'z_score': z,
        }


class TournamentSketch:
    """Compact data sketch using tournament score distributions.

    Compresses a dataset into a small summary (the score histogram
    of its bit-plane delta-row decomposition) that can answer:
    - Is this new data similar to the training data?
    - What's the approximate entropy?
    - What's the approximate compression ratio?
    """

    def __init__(self, width=ROW_WIDTH):
        self.width = width
        self.score_hist = np.zeros((8, width + 1), dtype=np.int64)  # per plane
        self.n_rows = 0
        self.n_chunks = 0

    def add(self, data):
        """Add data to the sketch."""
        arr = np.frombuffer(data, dtype=np.uint8)
        gray = arr ^ (arr >> 1)
        H = (len(data) + self.width - 1) // self.width
        padded = np.zeros(H * self.width, dtype=np.uint8)
        padded[:len(data)] = gray
        matrix = padded.reshape(H, self.width)

        for b in range(8):
            plane = ((matrix >> (7 - b)) & 1).astype(np.uint8)
            delta = np.zeros_like(plane)
            delta[0] = plane[0]
            for i in range(1, H): delta[i] = plane[i] ^ plane[i - 1]
            for i in range(H):
                weight = int(np.sum(delta[i]))
                self.score_hist[b, weight] += 1

        self.n_rows += H
        self.n_chunks += 1

    def similarity(self, other):
        """Cosine similarity between two sketches."""
        a = self.score_hist.flatten().astype(np.float64)
        b = other.score_hist.flatten().astype(np.float64)
        dot = np.dot(a, b)
        norm = np.sqrt(np.dot(a, a) * np.dot(b, b))
        return dot / max(norm, 1e-10)

    def estimated_ratio(self):
        """Estimated compression ratio from the score distribution."""
        total_bits = 0
        raw_bits = self.n_rows * self.width * 8
        for b in range(8):
            for w in range(self.width + 1):
                count = self.score_hist[b, w]
                if count > 0 and 0 < w < self.width:
                    total_bits += count * log2(max(comb(self.width, w), 1))
                    # Score encoding cost
                    total_bits += count * log2(self.width + 1)
        return raw_bits / max(total_bits, 1)

    def size_bytes(self):
        """Size of the sketch in bytes."""
        return self.score_hist.nbytes + 16  # + metadata

    def __repr__(self):
        return f"TournamentSketch(rows={self.n_rows}, chunks={self.n_chunks}, " \
               f"size={self.size_bytes()}B)"


# ============================================================================
# DEMO
# ============================================================================

def main():
    print("=" * 72)
    print("  TOURNAMENT ANOMALY DETECTOR + DATA SKETCH")
    print("  opus-2026-03-24-S315")
    print("=" * 72)

    np.random.seed(42)

    # ============================================================
    # 1. ANOMALY DETECTION
    # ============================================================
    print(f"\n{'='*72}")
    print(f"  1. ANOMALY DETECTION ON SENSOR DATA")
    print(f"{'='*72}")

    # Generate normal sensor data: slow ramp + small noise
    def gen_normal(n=512):
        t = np.arange(n)
        signal = (t * 0.1 + np.random.normal(0, 2, n)) % 256
        return np.clip(signal, 0, 255).astype(np.uint8).tobytes()

    # Generate anomalous data: sudden spike / mode change
    def gen_anomaly_spike(n=512):
        signal = np.ones(n) * 128
        signal[200:300] = 255  # sudden spike
        signal += np.random.normal(0, 2, n)
        return np.clip(signal, 0, 255).astype(np.uint8).tobytes()

    def gen_anomaly_noise(n=512):
        return bytes(np.random.randint(0, 256, n).astype(np.uint8))

    def gen_anomaly_constant(n=512):
        return b'\x80' * n

    # Train detector
    detector = TournamentAnomalyDetector(window_size=30, threshold_sigma=2.5)
    train_chunks = [gen_normal() for _ in range(50)]
    detector.train(train_chunks)
    print(f"  Trained on {len(train_chunks)} normal chunks")
    print(f"  Baseline ratio: {detector.baseline_mean:.3f} ± {detector.baseline_std:.3f}")

    # Test
    print(f"\n  {'Chunk type':>20} {'Ratio':>7} {'Z-score':>8} {'Anomaly':>8}")
    print(f"  {'-'*20} {'-'*7} {'-'*8} {'-'*8}")

    for label, gen_fn, expected in [
        ("normal_1",       gen_normal,           False),
        ("normal_2",       gen_normal,           False),
        ("normal_3",       gen_normal,           False),
        ("spike",          gen_anomaly_spike,    True),
        ("random_noise",   gen_anomaly_noise,    True),
        ("constant",       gen_anomaly_constant, True),
        ("normal_4",       gen_normal,           False),
        ("normal_5",       gen_normal,           False),
    ]:
        chunk = gen_fn()
        is_anom, z, details = detector.detect(chunk)
        status = "⚠ ALERT" if is_anom else "  ok"
        correct = "✓" if (is_anom == expected) else "✗"
        print(f"  {label:>20} {details['ratio']:>7.3f} {z:>+8.2f} {status:>8} {correct}")

    # ============================================================
    # 2. DATA SKETCHING
    # ============================================================
    print(f"\n{'='*72}")
    print(f"  2. DATA SKETCHING — Compact similarity fingerprints")
    print(f"{'='*72}")

    # Create sketches of different data types
    sketch_sensor = TournamentSketch()
    for _ in range(100):
        sketch_sensor.add(gen_normal())

    sketch_text = TournamentSketch()
    import random
    random.seed(42)
    for _ in range(100):
        text = ''.join(random.choices('etaoinshrdlu ETAOIN.,\n', k=512))
        sketch_text.add(text.encode('ascii'))

    sketch_random = TournamentSketch()
    for _ in range(100):
        sketch_random.add(gen_anomaly_noise())

    sketch_gradient = TournamentSketch()
    for _ in range(100):
        offset = np.random.randint(256)
        data = bytes((i + offset) % 256 for i in range(512))
        sketch_gradient.add(data)

    sketches = {
        'sensor':   sketch_sensor,
        'text':     sketch_text,
        'random':   sketch_random,
        'gradient': sketch_gradient,
    }

    print(f"\n  Sketch sizes: {sketch_sensor.size_bytes()}B each ({sketch_sensor})")
    print(f"  vs raw data: {100 * 512}B per dataset")
    print(f"  Compression: {100 * 512 / sketch_sensor.size_bytes():.0f}x")

    print(f"\n  SIMILARITY MATRIX (cosine of score distributions):")
    names = list(sketches.keys())
    print(f"  {'':>10}", end='')
    for n in names: print(f" {n:>10}", end='')
    print()
    for n1 in names:
        print(f"  {n1:>10}", end='')
        for n2 in names:
            sim = sketches[n1].similarity(sketches[n2])
            print(f" {sim:>10.4f}", end='')
        print()

    # Test: which sketch is a new chunk most similar to?
    print(f"\n  CLASSIFICATION TEST (new chunks → nearest sketch):")
    test_chunks = [
        ("sensor_new",  gen_normal()),
        ("text_new",    ''.join(random.choices('etaoinshrdlu ', k=512)).encode()),
        ("random_new",  gen_anomaly_noise()),
        ("gradient_new", bytes(i % 256 for i in range(512))),
    ]

    for label, chunk in test_chunks:
        query = TournamentSketch()
        query.add(chunk)
        sims = {name: query.similarity(sk) for name, sk in sketches.items()}
        best = max(sims.keys(), key=lambda k: sims[k])
        print(f"    {label:>15} → {best:>10} "
              f"(sim={sims[best]:.4f}, 2nd={sorted(sims.values())[-2]:.4f})")

    # ============================================================
    # 3. COMPRESSION-BASED CLUSTERING
    # ============================================================
    print(f"\n{'='*72}")
    print(f"  3. COMPRESSION-BASED CLUSTERING")
    print(f"{'='*72}")

    # NCD (Normalized Compression Distance) using our codec
    def ncd(x, y):
        """Normalized Compression Distance using zlib."""
        cx = len(zlib.compress(x, 9))
        cy = len(zlib.compress(y, 9))
        cxy = len(zlib.compress(x + y, 9))
        return (cxy - min(cx, cy)) / max(cx, cy)

    # Generate clusters
    clusters = {
        'A_gradient': [bytes((i + k * 10) % 256 for i in range(256)) for k in range(5)],
        'B_text': [''.join(random.choices('hello world test ', k=256)).encode() for _ in range(5)],
        'C_sparse': [bytes(np.random.RandomState(k).choice([0]*95+[1]*5, 256).astype(np.uint8))
                    for k in range(5)],
    }

    print(f"\n  NCD matrix (lower = more similar):")
    all_items = []
    for cname, items in clusters.items():
        for i, item in enumerate(items):
            all_items.append((f"{cname[:1]}{i}", item))

    # Just show intra-cluster vs inter-cluster distances
    intra_dists = []
    inter_dists = []
    for cname, items in clusters.items():
        for i in range(len(items)):
            for j in range(i + 1, len(items)):
                intra_dists.append(ncd(items[i], items[j]))
        for other_name, other_items in clusters.items():
            if other_name == cname: continue
            for item1 in items[:2]:
                for item2 in other_items[:2]:
                    inter_dists.append(ncd(item1, item2))

    print(f"    Intra-cluster NCD: mean={np.mean(intra_dists):.3f}, std={np.std(intra_dists):.3f}")
    print(f"    Inter-cluster NCD: mean={np.mean(inter_dists):.3f}, std={np.std(inter_dists):.3f}")
    print(f"    Separation: {np.mean(inter_dists) - np.mean(intra_dists):.3f}")
    print(f"    → {'GOOD separation ✓' if np.mean(inter_dists) > np.mean(intra_dists) + np.std(intra_dists) else 'Weak separation'}")

    print(f"\n{'='*72}")
    print(f"  APPLICATIONS")
    print(f"{'='*72}")
    print("""
    1. ANOMALY DETECTION:
       - Train on normal data → learn score distribution
       - Anomalous data compresses differently → z-score alert
       - Use cases: intrusion detection, sensor faults, fraud

    2. DATA SKETCHING:
       - Compress dataset into 280-byte score histogram
       - 180x compression of the distribution
       - Cosine similarity answers "is this data like that data?"
       - Use cases: data deduplication, approximate search, monitoring

    3. COMPRESSION-BASED CLUSTERING:
       - NCD (Normalized Compression Distance) groups similar data
       - No feature engineering needed — the compressor IS the model
       - Use cases: malware classification, document clustering, music genre
    """)

    print("DONE.")

if __name__ == '__main__':
    main()
