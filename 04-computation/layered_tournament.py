#!/usr/bin/env python3
"""
layered_tournament.py — Tournament as a sequence of bit layers
A PRACTICAL PRODUCT exploiting the recursive tiling structure.

Layer k (0-indexed): k bits encoding vertex (k+2)'s arcs to vertices 0..k-1.
Total bits: 0 + 1 + 2 + ... + (n-3) = C(n-2,2)...

Wait — let me get the indexing right from the staircase.

In the tiling model (base path n->n-1->...->1):
  Tile (x,y) where x-y >= 2.
  At SKIP s = x-y-1: tiles are (s+2,1), (s+3,2), ..., (n, n-s-1).
  Count at skip s: n-s-2 tiles.

But for the VERTEX ADDITION model:
  Add vertex v (v=2,3,...,n-1 in 0-indexed adjacency).
  Vertex v connects to vertices 0,...,v-1.
  Base arc: v -> v-1 (always forward).
  Free arcs: v to 0,...,v-2. That's v-1 free arcs.

  Layer v-2 (0-indexed): v-1 bits for vertex v.
  Layer 0: vertex 2, 1 bit.
  Layer 1: vertex 3, 2 bits.
  Layer k: vertex k+2, k+1 bits.
  Total: 1+2+...+(n-2) = C(n-1,2). Correct!

FAST OPERATIONS:
  - Random: generate each layer independently
  - Score: popcount per layer
  - Comparison: layer-by-layer, early exit
  - Sub-tournament: truncate layers
  - Tile flip: XOR one bit in one layer
"""

import time
import random
from math import comb


class LayeredTournament:
    """Tournament stored as sequence of bit layers for maximum efficiency."""

    __slots__ = ['n', 'layers']

    def __init__(self, n, layers=None):
        self.n = n
        if layers is None:
            self.layers = tuple(0 for _ in range(n-2))
        else:
            self.layers = tuple(layers)

    @classmethod
    def random(cls, n):
        """O(n) random tournament generation."""
        return cls(n, [random.getrandbits(k+1) & ((1 << (k+1)) - 1) for k in range(n-2)])

    @classmethod
    def transitive(cls, n):
        """The transitive tournament (all layers = 0)."""
        return cls(n)

    def score_sequence(self):
        """O(n^2) score computation via layer popcount."""
        n = self.n
        scores = [0] * n
        # Base path: vertex i beats i+1
        for i in range(n-1):
            scores[i] += 1
        # Layer k: vertex k+2's free arcs to 0,...,k
        for k in range(n-2):
            v = k + 2
            layer = self.layers[k]
            for j in range(k+1):
                if layer & (1 << j):
                    scores[j] += 1  # j beats v (backward)
                else:
                    scores[v] += 1  # v beats j (forward)
        return tuple(sorted(scores))

    def hamming_weight(self):
        """O(n) total backward arcs."""
        return sum(bin(l).count('1') for l in self.layers)

    def vertex_score(self, v):
        """O(n) score of a single vertex."""
        n = self.n
        s = 0
        # Base path contributions
        if v < n-1: s += 1  # v beats v+1
        # Layer contributions
        # As the adder (layer v-2): count forward arcs
        if v >= 2:
            layer = self.layers[v-2]
            s += (v-1) - bin(layer).count('1')  # forward = not backward
        # As target in later layers
        for k in range(max(v-1, 0), n-2):
            layer_v = k + 2
            if v <= k:
                if not (self.layers[k] & (1 << v)):
                    pass  # layer_v beats v (doesn't add to v's score)
                else:
                    s += 1  # v beats layer_v
        return s

    def flip(self, layer_idx, bit_idx):
        """Return NEW tournament with one tile flipped."""
        new_layers = list(self.layers)
        new_layers[layer_idx] ^= (1 << bit_idx)
        return LayeredTournament(self.n, new_layers)

    def progressive_equal(self, other):
        """Fast equality via early exit. Returns (equal, first_diff_layer)."""
        for k in range(self.n - 2):
            if self.layers[k] != other.layers[k]:
                return False, k
        return True, -1

    def sub(self, up_to_vertex):
        """Sub-tournament on vertices {0,...,up_to_vertex}. O(1)."""
        return LayeredTournament(up_to_vertex + 1, self.layers[:up_to_vertex - 1])

    def fingerprint(self):
        """Fast fingerprint for isomorphism pre-screening."""
        return (self.score_sequence(), self.hamming_weight())

    def layer_weights(self):
        """Hamming weight per layer (diagnostic)."""
        return tuple(bin(l).count('1') for l in self.layers)

    def to_int(self):
        """Pack all layers into a single integer."""
        result = 0
        offset = 0
        for k in range(self.n - 2):
            result |= (self.layers[k] << offset)
            offset += k + 1
        return result

    @classmethod
    def from_int(cls, n, bits):
        """Unpack from a single integer."""
        layers = []
        offset = 0
        for k in range(n - 2):
            width = k + 1
            layer = (bits >> offset) & ((1 << width) - 1)
            layers.append(layer)
            offset += width
        return cls(n, layers)


def benchmark():
    """Comprehensive speed benchmark."""
    print("=" * 70)
    print("  LAYERED TOURNAMENT ENGINE — BENCHMARK")
    print("=" * 70)

    print(f"\n  {'n':>6} {'random':>10} {'score':>10} {'hw':>10} {'fp':>10} {'cmp':>10} {'flip':>10}")
    print(f"  {'':>6} {'(us)':>10} {'(us)':>10} {'(us)':>10} {'(us)':>10} {'(us)':>10} {'(us)':>10}")

    for n in [10, 50, 100, 500, 1000, 5000]:
        reps = 10000 if n <= 100 else 1000 if n <= 1000 else 100

        # Random generation
        t0 = time.time()
        ts = [LayeredTournament.random(n) for _ in range(reps)]
        gen = (time.time() - t0) / reps * 1e6

        # Score
        t0 = time.time()
        for t in ts[:min(reps, 1000)]:
            t.score_sequence()
        sc = (time.time() - t0) / min(reps, 1000) * 1e6

        # Hamming weight
        t0 = time.time()
        for t in ts[:min(reps, 5000)]:
            t.hamming_weight()
        hw = (time.time() - t0) / min(reps, 5000) * 1e6

        # Fingerprint
        t0 = time.time()
        for t in ts[:min(reps, 1000)]:
            t.fingerprint()
        fp = (time.time() - t0) / min(reps, 1000) * 1e6

        # Comparison
        t0 = time.time()
        for i in range(min(reps-1, 4999)):
            ts[i].progressive_equal(ts[(i+1) % len(ts)])
        cm = (time.time() - t0) / min(reps-1, 4999) * 1e6

        # Flip
        t0 = time.time()
        for t in ts[:min(reps, 5000)]:
            t.flip(min(n-3, 5), 0)
        fl = (time.time() - t0) / min(reps, 5000) * 1e6

        print(f"  {n:6d} {gen:10.1f} {sc:10.1f} {hw:10.1f} {fp:10.1f} {cm:10.1f} {fl:10.1f}")

    # Progressive comparison advantage
    print(f"\n  PROGRESSIVE COMPARISON EARLY EXIT:")
    for n in [100, 1000, 10000]:
        t1 = LayeredTournament.random(n)
        t2 = LayeredTournament.random(n)
        eq, diff_layer = t1.progressive_equal(t2)
        bits_checked = sum(range(1, diff_layer + 2)) if diff_layer >= 0 else comb(n-1, 2)
        total_bits = comb(n-1, 2)
        print(f"    n={n}: diff at layer {diff_layer}, checked {bits_checked}/{total_bits} bits = {bits_checked/total_bits*100:.2f}%")

    # Pack/unpack roundtrip
    print(f"\n  PACK/UNPACK:")
    for n in [10, 50, 100]:
        t = LayeredTournament.random(n)
        packed = t.to_int()
        t2 = LayeredTournament.from_int(n, packed)
        assert t.layers == t2.layers
        print(f"    n={n}: {comb(n-1,2)} bits packed to {packed.bit_length()}-bit int. Roundtrip OK.")


if __name__ == "__main__":
    benchmark()
    print("\nDONE.")
