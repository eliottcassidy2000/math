#!/usr/bin/env python3
"""
streaming_tournament_format.py — Streaming tournament data format using the chain
A PRACTICAL PRODUCT.

THE CHAIN: delta_0 -> square_1 -> delta_1 -> square_2 -> ... -> delta_k

A tournament grows vertex-by-vertex. At each step, k new bits arrive.
The STREAMING FORMAT encodes this naturally:

  HEADER: n (number of vertices)
  LAYER 1: 1 bit   (vertex 2's choice)
  LAYER 2: 2 bits  (vertex 3's choices)
  ...
  LAYER k: k bits  (vertex k+1's choices)

FEATURES:
  1. PROGRESSIVE: can decode after any prefix (partial tournament)
  2. STREAMABLE: process observations as they arrive
  3. COMPRESSIBLE: each layer can be entropy-coded independently
  4. RANDOM ACCESS: jump to any layer in O(1) with offset table
  5. DELTA-ENCODABLE: difference between two tournaments = XOR per layer

APPLICATIONS:
  - Live sports: encode round-robin results as matches happen
  - A/B testing: stream comparison results to analytics
  - Database: store tournament collection with shared prefixes (trie)
  - Network: send partial rankings with progressive refinement
"""

import sys
import struct
import time
import random
from math import comb
from collections import defaultdict

__version__ = "1.0.0"


class StreamingTournamentFormat:
    """Efficient streaming format for tournament data."""

    def __init__(self, n):
        self.n = n
        self.layers = []
        self.layer_offsets = [0]  # bit offsets for random access
        offset = 0
        for k in range(n - 2):
            offset += k + 1
            self.layer_offsets.append(offset)
        self.total_bits = offset

    def encode_tournament(self, layers):
        """Encode a tournament (list of layer ints) to bytes."""
        # Pack layers into a byte stream
        bits = 0
        offset = 0
        for k, layer in enumerate(layers):
            width = k + 1
            bits |= (layer & ((1 << width) - 1)) << offset
            offset += width

        # Convert to bytes
        n_bytes = (offset + 7) // 8
        return bits.to_bytes(n_bytes, 'little')

    def decode_tournament(self, data):
        """Decode bytes back to layer list."""
        bits = int.from_bytes(data, 'little')
        layers = []
        offset = 0
        for k in range(self.n - 2):
            width = k + 1
            layer = (bits >> offset) & ((1 << width) - 1)
            layers.append(layer)
            offset += width
        return layers

    def encode_progressive(self, layers, up_to_layer=None):
        """Encode partial tournament (first k layers only)."""
        if up_to_layer is None:
            up_to_layer = len(layers)
        return self.encode_tournament(layers[:up_to_layer])

    def delta_encode(self, layers1, layers2):
        """Delta between two tournaments: XOR per layer."""
        assert len(layers1) == len(layers2)
        return [l1 ^ l2 for l1, l2 in zip(layers1, layers2)]

    def score_from_layers(self, layers):
        """Compute score sequence from layers without full adjacency."""
        n = len(layers) + 2
        scores = [0] * n
        for i in range(n - 1):
            scores[i] += 1  # base path
        for k, layer in enumerate(layers):
            v = k + 2
            for j in range(k + 1):
                if layer & (1 << j):
                    scores[j] += 1
                else:
                    scores[v] += 1
        return tuple(sorted(scores))

    def hamming_weight_from_layers(self, layers):
        """Total backward arcs from layers."""
        return sum(bin(l).count('1') for l in layers)

    def shared_prefix_length(self, layers1, layers2):
        """How many layers are identical (for trie compression)."""
        for k in range(min(len(layers1), len(layers2))):
            if layers1[k] != layers2[k]:
                return k
        return min(len(layers1), len(layers2))


def demo():
    print("=" * 70)
    print("  STREAMING TOURNAMENT FORMAT")
    print("=" * 70)

    # Demo 1: Encode/decode
    print("\n  1. ENCODE/DECODE:")
    for n in [5, 10, 20, 50, 100]:
        fmt = StreamingTournamentFormat(n)
        layers = [random.getrandbits(k+1) & ((1 << (k+1)) - 1) for k in range(n-2)]

        encoded = fmt.encode_tournament(layers)
        decoded = fmt.decode_tournament(encoded)

        assert layers == decoded, "Roundtrip failed!"
        raw_bytes = (comb(n, 2) + 7) // 8
        fmt_bytes = len(encoded)
        savings = (1 - fmt_bytes / raw_bytes) * 100 if raw_bytes > 0 else 0

        print(f"    n={n:4d}: {fmt.total_bits} bits = {fmt_bytes} bytes (raw adj: {raw_bytes} bytes, saves {savings:.0f}%)")

    # Demo 2: Progressive decoding
    print("\n  2. PROGRESSIVE DECODING:")
    n = 10
    fmt = StreamingTournamentFormat(n)
    layers = [random.getrandbits(k+1) & ((1 << (k+1)) - 1) for k in range(n-2)]

    for up_to in [1, 2, 4, 8]:
        if up_to > n - 2: break
        partial = fmt.encode_progressive(layers, up_to)
        partial_score = fmt.score_from_layers(layers[:up_to])
        vertices = up_to + 2
        print(f"    After {up_to} layers ({vertices} vertices): {len(partial)} bytes, score={partial_score}")

    # Demo 3: Delta encoding
    print("\n  3. DELTA ENCODING:")
    n = 20
    fmt = StreamingTournamentFormat(n)
    t1 = [random.getrandbits(k+1) & ((1 << (k+1)) - 1) for k in range(n-2)]
    t2 = list(t1)  # copy
    # Flip 3 random tiles
    for _ in range(3):
        layer = random.randint(0, n-3)
        bit = random.randint(0, layer)
        t2[layer] ^= (1 << bit)

    delta = fmt.delta_encode(t1, t2)
    delta_bits = sum(bin(d).count('1') for d in delta)
    delta_bytes = len(fmt.encode_tournament(delta))
    full_bytes = len(fmt.encode_tournament(t2))
    print(f"    Full tournament: {full_bytes} bytes")
    print(f"    Delta (3 flips): {delta_bits} bits changed, {delta_bytes} bytes")
    print(f"    Delta compression: {full_bytes/max(delta_bytes,1):.1f}x")

    # Demo 4: Trie compression (shared prefix)
    print("\n  4. TRIE COMPRESSION (shared prefixes):")
    n = 50
    fmt = StreamingTournamentFormat(n)
    base = [random.getrandbits(k+1) & ((1 << (k+1)) - 1) for k in range(n-2)]
    variants = []
    for _ in range(100):
        v = list(base)
        # Mutate a few high layers (late vertices)
        for _ in range(random.randint(1, 5)):
            layer = random.randint(max(0, n-10), n-3)
            bit = random.randint(0, layer)
            v[layer] ^= (1 << bit)
        variants.append(v)

    prefix_lengths = [fmt.shared_prefix_length(base, v) for v in variants]
    avg_prefix = sum(prefix_lengths) / len(prefix_lengths)
    prefix_bits = sum(sum(range(1, p+1)) for p in prefix_lengths) / len(prefix_lengths)
    total_bits = fmt.total_bits
    print(f"    n={n}: avg shared prefix = {avg_prefix:.1f} layers out of {n-2}")
    print(f"    Avg shared bits: {prefix_bits:.0f} / {total_bits} = {prefix_bits/total_bits*100:.0f}%")
    print(f"    Trie savings: {prefix_bits/total_bits*100:.0f}% of bits shared in trie nodes")

    # Demo 5: Speed
    print("\n  5. SPEED BENCHMARK:")
    for n in [10, 100, 1000]:
        fmt = StreamingTournamentFormat(n)
        layers = [random.getrandbits(k+1) & ((1 << (k+1)) - 1) for k in range(n-2)]

        reps = 10000 if n <= 100 else 1000
        t0 = time.time()
        for _ in range(reps):
            fmt.encode_tournament(layers)
        enc_time = (time.time() - t0) / reps * 1e6

        encoded = fmt.encode_tournament(layers)
        t0 = time.time()
        for _ in range(reps):
            fmt.decode_tournament(encoded)
        dec_time = (time.time() - t0) / reps * 1e6

        print(f"    n={n:5d}: encode={enc_time:.0f}us decode={dec_time:.0f}us ({len(encoded)} bytes)")


if __name__ == "__main__":
    demo()
    print("\nDONE.")
