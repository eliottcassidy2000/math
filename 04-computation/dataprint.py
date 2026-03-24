#!/usr/bin/env python3
"""
dataprint.py -- Universal Data Fingerprinter & Classifier
kind-pasteur-2026-03-24-S20cq

Fast structural fingerprint of any byte sequence. Uses entropy profiling,
stride correlation analysis, and transform responsiveness to classify
data type and recommend optimal compression.

USAGE:
  python dataprint.py FILE                    # fingerprint a file
  python dataprint.py --classify FILE         # classify data type
  python dataprint.py --recommend FILE        # recommend compression
  python dataprint.py --compare FILE1 FILE2   # structural similarity

OUTPUT: A DataPrint object with:
  - entropy:       Shannon entropy (bits/byte)
  - byte_hist:     256-bin histogram (normalized)
  - stride_corr:   correlation at strides 1,2,3,4,8,16
  - delta_entropy: entropy after delta encoding
  - run_density:   fraction of bytes in runs >= 4
  - unique_ratio:  unique_bytes / 256
  - ascii_ratio:   fraction of printable ASCII bytes
  - null_ratio:    fraction of zero bytes
  - bigram_entropy: entropy of byte pairs
  - class_label:   classified type (text, binary, structured, random, sparse)

APPLICATIONS:
  1. Choose compression algorithm without trial-and-error
  2. Detect data corruption (fingerprint changes)
  3. Cluster similar files by structural similarity
  4. Pre-screen data before ML pipeline
  5. Network traffic classification

LICENSE: MIT
"""

import sys
import os
import math
import argparse
from collections import Counter
from dataclasses import dataclass, field
from typing import List, Tuple, Optional

__version__ = "1.0.0"


@dataclass
class DataPrint:
    """Structural fingerprint of a byte sequence."""
    size: int = 0
    entropy: float = 0.0
    delta_entropy: float = 0.0
    xor_entropy: float = 0.0
    bigram_entropy: float = 0.0
    unique_bytes: int = 0
    unique_ratio: float = 0.0
    ascii_ratio: float = 0.0
    null_ratio: float = 0.0
    run_density: float = 0.0
    max_run: int = 0
    stride_corrs: dict = field(default_factory=dict)
    best_stride: int = 1
    byte_hist: list = field(default_factory=list)
    top_bytes: list = field(default_factory=list)
    class_label: str = "unknown"
    class_confidence: float = 0.0
    recommended_pipeline: str = ""
    recommended_backend: str = ""

    def summary(self) -> str:
        lines = [
            f"DataPrint v{__version__}",
            f"  Size:          {self.size:,} bytes",
            f"  Entropy:       {self.entropy:.3f} bits/byte ({self.entropy/8*100:.1f}%)",
            f"  Delta entropy: {self.delta_entropy:.3f} bits/byte",
            f"  XOR entropy:   {self.xor_entropy:.3f} bits/byte",
            f"  Bigram ent:    {self.bigram_entropy:.3f} bits/pair",
            f"  Unique bytes:  {self.unique_bytes}/256 ({self.unique_ratio:.1%})",
            f"  ASCII ratio:   {self.ascii_ratio:.1%}",
            f"  Null ratio:    {self.null_ratio:.1%}",
            f"  Run density:   {self.run_density:.1%} (max run: {self.max_run})",
            f"  Best stride:   {self.best_stride} (corr={self.stride_corrs.get(self.best_stride, 0):.3f})",
            f"  Stride corrs:  " + ", ".join(f"s{k}={v:.3f}" for k, v in sorted(self.stride_corrs.items())),
            f"  Top bytes:     {self.top_bytes[:5]}",
            f"  Class:         {self.class_label} ({self.class_confidence:.0%})",
            f"  Recommended:   {self.recommended_pipeline} + {self.recommended_backend}",
        ]
        return "\n".join(lines)

    def to_vector(self) -> list:
        """Convert to fixed-size feature vector (for ML)."""
        return [
            self.entropy,
            self.delta_entropy,
            self.xor_entropy,
            self.bigram_entropy,
            self.unique_ratio,
            self.ascii_ratio,
            self.null_ratio,
            self.run_density,
            self.stride_corrs.get(1, 0),
            self.stride_corrs.get(2, 0),
            self.stride_corrs.get(3, 0),
            self.stride_corrs.get(4, 0),
            self.stride_corrs.get(8, 0),
            self.stride_corrs.get(16, 0),
        ]


def _entropy(data: bytes) -> float:
    """Shannon entropy in bits/byte."""
    if not data: return 0.0
    counts = Counter(data)
    n = len(data)
    return -sum(c/n * math.log2(c/n) for c in counts.values() if c > 0)


def _delta(data: bytes) -> bytes:
    if len(data) < 2: return data
    out = bytearray(len(data))
    out[0] = data[0]
    for i in range(1, len(data)):
        out[i] = (data[i] - data[i-1]) & 0xFF
    return bytes(out)


def _xor(data: bytes) -> bytes:
    if len(data) < 2: return data
    out = bytearray(len(data))
    out[0] = data[0]
    for i in range(1, len(data)):
        out[i] = data[i] ^ data[i-1]
    return bytes(out)


def _stride_correlation(data: bytes, stride: int) -> float:
    """Pearson correlation between byte[i] and byte[i+stride]."""
    if len(data) <= stride: return 0.0
    n = len(data) - stride
    if n < 10: return 0.0
    # Use a sample for speed
    sample_size = min(n, 4096)
    step = max(1, n // sample_size)

    sx = sy = sxx = syy = sxy = 0.0
    count = 0
    for i in range(0, n, step):
        x = data[i]
        y = data[i + stride]
        sx += x
        sy += y
        sxx += x * x
        syy += y * y
        sxy += x * y
        count += 1

    if count < 2: return 0.0
    mx = sx / count
    my = sy / count
    vx = sxx / count - mx * mx
    vy = syy / count - my * my
    if vx <= 0 or vy <= 0: return 0.0
    cov = sxy / count - mx * my
    return cov / math.sqrt(vx * vy)


def _run_analysis(data: bytes) -> Tuple[float, int]:
    """Return (run_density, max_run_length)."""
    if not data: return 0.0, 0
    run_bytes = 0
    max_run = 1
    current_run = 1
    for i in range(1, len(data)):
        if data[i] == data[i-1]:
            current_run += 1
        else:
            if current_run >= 4:
                run_bytes += current_run
            max_run = max(max_run, current_run)
            current_run = 1
    if current_run >= 4:
        run_bytes += current_run
    max_run = max(max_run, current_run)
    return run_bytes / len(data), max_run


def _bigram_entropy(data: bytes) -> float:
    """Entropy of consecutive byte pairs."""
    if len(data) < 2: return 0.0
    # Sample for speed
    sample_size = min(len(data) - 1, 8192)
    step = max(1, (len(data) - 1) // sample_size)
    pairs = Counter()
    count = 0
    for i in range(0, len(data) - 1, step):
        pairs[(data[i], data[i+1])] += 1
        count += 1
    if count == 0: return 0.0
    return -sum(c/count * math.log2(c/count) for c in pairs.values() if c > 0)


def fingerprint(data: bytes) -> DataPrint:
    """Compute full structural fingerprint."""
    dp = DataPrint()
    dp.size = len(data)
    if not data:
        return dp

    # Basic entropy
    dp.entropy = _entropy(data)
    dp.delta_entropy = _entropy(_delta(data))
    dp.xor_entropy = _entropy(_xor(data))
    dp.bigram_entropy = _bigram_entropy(data)

    # Byte statistics
    counts = Counter(data)
    dp.unique_bytes = len(counts)
    dp.unique_ratio = dp.unique_bytes / 256
    dp.ascii_ratio = sum(1 for b in data if 32 <= b <= 126 or b in (9, 10, 13)) / len(data)
    dp.null_ratio = counts.get(0, 0) / len(data)
    dp.top_bytes = [(b, c) for b, c in counts.most_common(10)]

    # Run analysis
    dp.run_density, dp.max_run = _run_analysis(data)

    # Stride correlations
    for s in [1, 2, 3, 4, 8, 16]:
        dp.stride_corrs[s] = abs(_stride_correlation(data, s))
    # Best stride = highest correlation beyond stride 1
    best_s = max([s for s in [2, 3, 4, 8, 16] if s in dp.stride_corrs],
                 key=lambda s: dp.stride_corrs[s], default=1)
    dp.best_stride = best_s if dp.stride_corrs.get(best_s, 0) > 0.1 else 1

    # Histogram
    dp.byte_hist = [counts.get(i, 0) / len(data) for i in range(256)]

    # Classify
    dp.class_label, dp.class_confidence = _classify(dp)

    # Recommend compression
    dp.recommended_pipeline, dp.recommended_backend = _recommend(dp)

    return dp


def _classify(dp: DataPrint) -> Tuple[str, float]:
    """Classify data type from fingerprint."""
    scores = {}

    # Text: high ASCII ratio, moderate entropy, low null ratio
    text_score = dp.ascii_ratio * 0.5 + (1 - dp.null_ratio) * 0.2
    if dp.entropy > 3.0 and dp.entropy < 6.5:
        text_score += 0.3
    scores['text'] = text_score

    # Structured: strong stride correlation, delta reduces entropy significantly
    struct_score = max(dp.stride_corrs.get(s, 0) for s in [2, 3, 4, 8, 16]) * 0.4
    entropy_reduction = (dp.entropy - dp.delta_entropy) / dp.entropy if dp.entropy > 0 else 0
    struct_score += max(0, entropy_reduction) * 0.4
    if dp.unique_ratio < 0.5:
        struct_score += 0.2
    scores['structured'] = struct_score

    # Random: high entropy, near-uniform distribution, no patterns
    random_score = dp.entropy / 8.0 * 0.5
    if dp.unique_ratio > 0.9:
        random_score += 0.2
    if dp.delta_entropy >= dp.entropy * 0.98:
        random_score += 0.3
    scores['random'] = random_score

    # Sparse: high null ratio, low entropy
    sparse_score = dp.null_ratio * 0.5 + dp.run_density * 0.3
    if dp.entropy < 3.0:
        sparse_score += 0.2
    scores['sparse'] = sparse_score

    # Binary executable: moderate entropy, specific byte patterns
    binary_score = 0.0
    if 5.0 < dp.entropy < 7.5 and dp.ascii_ratio < 0.5:
        binary_score = 0.5
    if dp.stride_corrs.get(4, 0) > 0.1 or dp.stride_corrs.get(8, 0) > 0.1:
        binary_score += 0.2
    scores['binary'] = binary_score

    # Low-entropy: few unique values, very low entropy
    low_ent_score = 0.0
    if dp.unique_bytes <= 16:
        low_ent_score = 0.5 + (1 - dp.unique_bytes / 16) * 0.5
    scores['low-entropy'] = low_ent_score

    best = max(scores, key=scores.get)
    confidence = scores[best] / (sum(scores.values()) + 0.01)
    return best, min(confidence, 0.99)


def _recommend(dp: DataPrint) -> Tuple[str, str]:
    """Recommend optimal compression pipeline."""
    # Based on class and features
    c = dp.class_label

    if c == 'sparse' or dp.run_density > 0.3:
        return 'rle+delta', 'zlib'

    if c == 'structured':
        if dp.best_stride > 1:
            return f'stride-{dp.best_stride}+delta', 'zlib'
        if dp.delta_entropy < dp.entropy * 0.5:
            return 'delta', 'zlib'
        return 'delta', 'lzma'

    if c == 'text':
        if dp.size < 8192:
            return 'bwt+mtf', 'zlib'
        return 'stride-3', 'bz2'

    if c == 'random':
        return 'raw', 'zlib-1'

    if c == 'low-entropy':
        return 'nibble-split', 'bz2'

    if c == 'binary':
        if dp.stride_corrs.get(4, 0) > dp.stride_corrs.get(8, 0):
            return 'stride-4+delta', 'zlib'
        return 'stride-8+delta', 'zlib'

    return 'raw', 'lzma'


def similarity(dp1: DataPrint, dp2: DataPrint) -> float:
    """Structural similarity between two fingerprints (0 to 1)."""
    v1 = dp1.to_vector()
    v2 = dp2.to_vector()
    # Cosine similarity
    dot = sum(a * b for a, b in zip(v1, v2))
    mag1 = math.sqrt(sum(a * a for a in v1))
    mag2 = math.sqrt(sum(b * b for b in v2))
    if mag1 == 0 or mag2 == 0: return 0.0
    return dot / (mag1 * mag2)


# ============================================================================
# CLI
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description=f'dataprint v{__version__} -- Universal Data Fingerprinter',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('files', nargs='*', help='Files to fingerprint')
    parser.add_argument('--classify', action='store_true', help='Show classification only')
    parser.add_argument('--recommend', action='store_true', help='Show compression recommendation')
    parser.add_argument('--compare', action='store_true', help='Compare two files')
    parser.add_argument('--vector', action='store_true', help='Output feature vector (for ML)')
    parser.add_argument('--demo', action='store_true', help='Run demo on synthetic data')

    args = parser.parse_args()

    if args.demo:
        _run_demo()
        return

    if not args.files:
        parser.print_help()
        return

    if args.compare and len(args.files) >= 2:
        with open(args.files[0], 'rb') as f:
            dp1 = fingerprint(f.read())
        with open(args.files[1], 'rb') as f:
            dp2 = fingerprint(f.read())
        sim = similarity(dp1, dp2)
        print(f"Similarity: {sim:.4f}")
        print(f"\n  {args.files[0]}: {dp1.class_label} ({dp1.entropy:.2f} b/B)")
        print(f"  {args.files[1]}: {dp2.class_label} ({dp2.entropy:.2f} b/B)")
        return

    for filepath in args.files:
        with open(filepath, 'rb') as f:
            data = f.read()
        dp = fingerprint(data)

        if args.classify:
            print(f"{filepath}: {dp.class_label} ({dp.class_confidence:.0%})")
        elif args.recommend:
            print(f"{filepath}: {dp.recommended_pipeline} + {dp.recommended_backend}")
        elif args.vector:
            print(f"{filepath}: {dp.to_vector()}")
        else:
            print(f"\n{'='*60}")
            print(f"  File: {filepath}")
            print(f"{'='*60}")
            print(dp.summary())


def _run_demo():
    """Demo on synthetic data types."""
    import time

    print(f"dataprint v{__version__} -- Demo")
    print("=" * 70)

    demos = {}

    # Text
    demos['english'] = (b"the quick brown fox jumps over the lazy dog. " * 100)[:4096]
    demos['python'] = (b"def f(x):\n    return x*x + 2*x + 1\n" * 120)[:4096]
    demos['json'] = (b'{"id":1,"name":"test","values":[1,2,3]}\n' * 105)[:4096]

    # Structured
    demos['counter'] = bytes([i % 256 for i in range(4096)])
    demos['quadratic'] = bytes([int((i*i/16) % 256) for i in range(4096)])

    # Binary
    demos['binary_exe'] = bytes([(i*7+13) % 256 for i in range(4096)])

    # Random
    import random
    random.seed(42)
    demos['random'] = bytes(random.randint(0, 255) for _ in range(4096))

    # Sparse
    demos['mostly_zeros'] = bytes([0]*3800 + list(range(256)) + [0]*40)
    demos['sparse'] = bytes(random.choice([0]*10 + [42]) for _ in range(4096))

    # Low entropy
    demos['low_ent_4'] = bytes(random.choice([0, 1, 2, 3]) for _ in range(4096))

    print(f"\n  {'Name':>15} {'Size':>6} {'Ent':>6} {'dEnt':>6} {'Uniq':>5} {'ASCII':>6} "
          f"{'Null':>6} {'Runs':>6} {'BestS':>5} {'Class':>12} {'Recommend':>25}")
    print("  " + "-" * 120)

    for name, data in demos.items():
        t0 = time.time()
        dp = fingerprint(data)
        elapsed = (time.time() - t0) * 1000

        print(f"  {name:>15} {dp.size:5d}B {dp.entropy:5.2f} {dp.delta_entropy:5.2f} "
              f"{dp.unique_bytes:4d} {dp.ascii_ratio:5.1%} {dp.null_ratio:5.1%} "
              f"{dp.run_density:5.1%} {dp.best_stride:4d} {dp.class_label:>12} "
              f"{dp.recommended_pipeline}+{dp.recommended_backend:>10}")


if __name__ == "__main__":
    main()
