#!/usr/bin/env python3
"""
tournament_hash.py -- Locality-Sensitive Hashing via Tournament Structure
kind-pasteur-2026-03-24-S20cq

THE INSIGHT: A tournament on n vertices is determined by C(n,2) arc directions.
If we map data to tournaments via pairwise comparisons of features, then:
  - SIMILAR data -> FEW arc reversals -> CLOSE in tournament space
  - DIFFERENT data -> MANY arc reversals -> FAR in tournament space

The "tournament distance" (Kendall tau distance = number of arc reversals)
is a well-studied metric. By projecting data onto random pairwise comparisons,
we get a locality-sensitive hash family.

KEY PROPERTY: Tournament hashes are ORDER-PRESERVING.
  If data is sorted by feature k, tournaments on feature k are TRANSITIVE.
  Partially sorted data -> partially transitive tournament -> low c3 count.
  This connects to our tournament theory: H(T) encodes structural regularity.

APPLICATIONS:
  1. Near-duplicate detection (documents, images, records)
  2. Approximate nearest neighbors (replacement for MinHash/SimHash)
  3. Clustering via tournament structure (c3 density = cluster quality)
  4. Data quality: tournament transitivity = consistency of data
  5. Anomaly detection: items that reverse many arcs are outliers

USAGE:
  from tournament_hash import TournamentHasher
  th = TournamentHasher(n_features=8, n_projections=16)
  hash1 = th.hash(vector1)
  hash2 = th.hash(vector2)
  dist = th.distance(hash1, hash2)    # 0 = identical, 1 = maximally different
  sim = th.similarity(hash1, hash2)   # 1 = identical, 0 = maximally different

  # Batch: build index and query
  th.build_index(vectors)
  neighbors = th.query(query_vector, k=10)

LICENSE: MIT
"""

import math
import sys
import struct
import hashlib
from collections import defaultdict
from dataclasses import dataclass, field
from typing import List, Tuple, Optional, Sequence

__version__ = "1.0.0"


class TournamentHash:
    """A tournament hash is a compact binary representation of a tournament.

    For n vertices, we need C(n,2) bits. Stored as an integer.
    Bit (i,j) with i < j is 1 if i beats j, 0 if j beats i.
    """

    __slots__ = ['n', 'bits', 'n_bits']

    def __init__(self, n: int, bits: int = 0):
        self.n = n
        self.n_bits = n * (n - 1) // 2
        self.bits = bits & ((1 << self.n_bits) - 1)

    def _bit_index(self, i: int, j: int) -> int:
        """Get bit index for pair (i, j) with i < j."""
        if i > j: i, j = j, i
        return i * self.n - i * (i + 1) // 2 + (j - i - 1)

    def set_arc(self, i: int, j: int, i_beats_j: bool):
        """Set arc direction."""
        if i == j: return
        if i > j:
            i, j = j, i
            i_beats_j = not i_beats_j
        idx = self._bit_index(i, j)
        if i_beats_j:
            self.bits |= (1 << idx)
        else:
            self.bits &= ~(1 << idx)

    def get_arc(self, i: int, j: int) -> bool:
        """Returns True if i beats j."""
        if i == j: return False
        if i > j:
            return not self.get_arc(j, i)
        idx = self._bit_index(i, j)
        return bool(self.bits & (1 << idx))

    def score(self, v: int) -> int:
        """Out-degree (score) of vertex v."""
        s = 0
        for j in range(self.n):
            if j != v and self.get_arc(v, j):
                s += 1
        return s

    def score_sequence(self) -> list:
        """Sorted score sequence."""
        return sorted(self.score(v) for v in range(self.n))

    def c3_count(self) -> int:
        """Number of directed 3-cycles."""
        n = self.n
        total = n * (n - 1) * (n - 2) // 6
        return total - sum(self.score(v) * (self.score(v) - 1) // 2 for v in range(n))

    def hamming_distance(self, other: 'TournamentHash') -> int:
        """Number of arc reversals (Kendall tau distance)."""
        assert self.n == other.n
        return bin(self.bits ^ other.bits).count('1')

    def normalized_distance(self, other: 'TournamentHash') -> float:
        """Normalized distance in [0, 1]."""
        return self.hamming_distance(other) / self.n_bits if self.n_bits > 0 else 0

    def to_bytes(self) -> bytes:
        """Serialize to bytes."""
        n_bytes = (self.n_bits + 7) // 8
        return struct.pack('!B', self.n) + self.bits.to_bytes(n_bytes, 'big')

    @classmethod
    def from_bytes(cls, data: bytes) -> 'TournamentHash':
        n = data[0]
        n_bits = n * (n - 1) // 2
        n_bytes = (n_bits + 7) // 8
        bits = int.from_bytes(data[1:1+n_bytes], 'big')
        return cls(n, bits)

    def hex(self) -> str:
        """Hex representation."""
        n_bytes = (self.n_bits + 7) // 8
        return f"T{self.n}:{self.bits.to_bytes(n_bytes, 'big').hex()}"

    def __eq__(self, other):
        return self.n == other.n and self.bits == other.bits

    def __hash__(self):
        return hash((self.n, self.bits))

    def __repr__(self):
        return f"TournamentHash(n={self.n}, bits=0x{self.bits:x}, c3={self.c3_count()})"


class TournamentHasher:
    """Locality-sensitive hasher using tournament structure.

    Maps vectors to tournaments via random pairwise feature comparisons.
    """

    def __init__(self, n_vertices: int = 8, n_projections: int = 0, seed: int = 42):
        """
        Args:
            n_vertices: tournament size (8 = 28 bits, 12 = 66 bits, 16 = 120 bits)
            n_projections: if > 0, use random projections to n_projections dimensions first
            seed: random seed for reproducibility
        """
        self.n = n_vertices
        self.n_projections = n_projections
        self.seed = seed
        self._projections = None
        self._index = {}  # hash -> list of (id, vector)
        self._hashes = {}  # id -> TournamentHash

    def _ensure_projections(self, dim: int):
        """Initialize random projections if needed."""
        if self._projections is not None and len(self._projections[0]) == dim:
            return
        import random
        rng = random.Random(self.seed)
        n_proj = self.n_projections if self.n_projections > 0 else self.n
        self._projections = [
            [rng.gauss(0, 1) for _ in range(dim)]
            for _ in range(n_proj)
        ]

    def hash(self, vector: Sequence[float]) -> TournamentHash:
        """Hash a vector to a tournament.

        The tournament is constructed by comparing projected values:
        vertex i beats vertex j if projection[i] . vector > projection[j] . vector
        """
        dim = len(vector)

        if self.n_projections > 0:
            # Use random projections
            self._ensure_projections(dim)
            projected = [sum(p[k] * vector[k] for k in range(dim))
                        for p in self._projections[:self.n]]
        elif dim >= self.n:
            # Use first n features directly
            projected = list(vector[:self.n])
        else:
            # Use features cyclically
            projected = [vector[i % dim] + i * 1e-10 for i in range(self.n)]

        # Build tournament: i beats j if projected[i] > projected[j]
        th = TournamentHash(self.n)
        for i in range(self.n):
            for j in range(i + 1, self.n):
                if projected[i] > projected[j]:
                    th.set_arc(i, j, True)
                elif projected[i] < projected[j]:
                    th.set_arc(i, j, False)
                else:
                    # Tie-break by index
                    th.set_arc(i, j, True)

        return th

    def hash_text(self, text: str, ngram: int = 3) -> TournamentHash:
        """Hash text via character n-gram frequency comparison.

        Each vertex corresponds to an n-gram bucket. Arc direction
        determined by frequency comparison.
        """
        # Count n-grams
        ngrams = defaultdict(int)
        for i in range(len(text) - ngram + 1):
            ngrams[text[i:i+ngram]] += 1

        # Hash n-grams to buckets
        n = self.n
        buckets = [0.0] * n
        for ng, count in ngrams.items():
            bucket_id = int(hashlib.md5(ng.encode()).hexdigest(), 16) % n
            buckets[bucket_id] += count

        # Build tournament from bucket values
        th = TournamentHash(n)
        for i in range(n):
            for j in range(i + 1, n):
                th.set_arc(i, j, buckets[i] > buckets[j])
        return th

    def hash_bytes(self, data: bytes) -> TournamentHash:
        """Hash raw bytes via byte frequency at different offsets."""
        n = self.n
        # Split data into n segments, count byte frequencies per segment
        seg_size = max(1, len(data) // n)
        segments = []
        for i in range(n):
            start = i * seg_size
            end = start + seg_size if i < n - 1 else len(data)
            segment = data[start:end]
            # Entropy-like feature
            if not segment:
                segments.append(0.0)
            else:
                counts = defaultdict(int)
                for b in segment:
                    counts[b] += 1
                total = len(segment)
                ent = -sum(c/total * math.log2(c/total) for c in counts.values() if c > 0)
                segments.append(ent)

        # Build tournament
        th = TournamentHash(n)
        for i in range(n):
            for j in range(i + 1, n):
                th.set_arc(i, j, segments[i] > segments[j])
        return th

    def distance(self, h1: TournamentHash, h2: TournamentHash) -> float:
        """Normalized distance between two hashes (0 = same, 1 = opposite)."""
        return h1.normalized_distance(h2)

    def similarity(self, h1: TournamentHash, h2: TournamentHash) -> float:
        """Similarity between two hashes (1 = same, 0 = opposite)."""
        return 1.0 - self.distance(h1, h2)

    def build_index(self, items: list, ids: Optional[list] = None):
        """Build search index from a list of items.

        Items can be vectors, strings, or bytes.
        """
        if ids is None:
            ids = list(range(len(items)))

        self._index.clear()
        self._hashes.clear()

        for item_id, item in zip(ids, items):
            if isinstance(item, str):
                h = self.hash_text(item)
            elif isinstance(item, bytes):
                h = self.hash_bytes(item)
            else:
                h = self.hash(item)

            self._hashes[item_id] = h
            # Store by exact hash for O(1) lookup
            key = h.bits
            if key not in self._index:
                self._index[key] = []
            self._index[key].append(item_id)

    def query(self, item, k: int = 10) -> List[Tuple[int, float]]:
        """Find k nearest neighbors by tournament distance.

        Returns: list of (item_id, distance) sorted by distance.
        """
        if isinstance(item, str):
            h = self.hash_text(item)
        elif isinstance(item, bytes):
            h = self.hash_bytes(item)
        elif isinstance(item, TournamentHash):
            h = item
        else:
            h = self.hash(item)

        # Brute force for now (fine for < 100K items)
        results = []
        for item_id, item_hash in self._hashes.items():
            dist = h.normalized_distance(item_hash)
            results.append((item_id, dist))

        results.sort(key=lambda x: x[1])
        return results[:k]

    def cluster_quality(self, labels: dict) -> float:
        """Measure cluster quality: items with same label should have low tournament distance.

        Args:
            labels: dict mapping item_id -> cluster_label

        Returns: float in [0, 1]. Higher = better clustering (same-cluster items are close).
        """
        intra = []
        inter = []
        ids = list(labels.keys())
        for i in range(len(ids)):
            for j in range(i + 1, len(ids)):
                id_i, id_j = ids[i], ids[j]
                if id_i not in self._hashes or id_j not in self._hashes:
                    continue
                dist = self._hashes[id_i].normalized_distance(self._hashes[id_j])
                if labels[id_i] == labels[id_j]:
                    intra.append(dist)
                else:
                    inter.append(dist)

        if not intra or not inter:
            return 0.0

        avg_intra = sum(intra) / len(intra)
        avg_inter = sum(inter) / len(inter)
        return (avg_inter - avg_intra) / (avg_inter + avg_intra + 1e-10)


# ============================================================================
# DEMO
# ============================================================================

def demo():
    """Comprehensive demo of tournament hashing."""
    import random
    random.seed(42)

    print(f"tournament_hash v{__version__} -- Demo")
    print("=" * 70)

    # 1. Vector hashing
    print("\n1. VECTOR HASHING (similar vectors -> similar tournaments)")
    th = TournamentHasher(n_vertices=8)

    v1 = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
    v2 = [1.1, 2.1, 3.1, 4.1, 5.1, 6.1, 7.1, 8.1]  # very similar
    v3 = [8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0]  # reversed
    v4 = [random.random() * 10 for _ in range(8)]     # random

    h1 = th.hash(v1)
    h2 = th.hash(v2)
    h3 = th.hash(v3)
    h4 = th.hash(v4)

    print(f"  v1 = [1..8]:    {h1.hex()}  scores={h1.score_sequence()}")
    print(f"  v2 = [1.1..8.1]: {h2.hex()}  scores={h2.score_sequence()}")
    print(f"  v3 = [8..1]:    {h3.hex()}  scores={h3.score_sequence()}")
    print(f"  v4 = random:    {h4.hex()}  scores={h4.score_sequence()}")
    print(f"  dist(v1,v2) = {th.distance(h1,h2):.3f} (should be ~0)")
    print(f"  dist(v1,v3) = {th.distance(h1,h3):.3f} (should be ~1)")
    print(f"  dist(v1,v4) = {th.distance(h1,h4):.3f} (should be ~0.5)")

    # 2. Text hashing
    print("\n2. TEXT HASHING (similar texts -> similar tournaments)")
    t1 = "the quick brown fox jumps over the lazy dog"
    t2 = "the fast brown fox leaps over the lazy dog"
    t3 = "completely unrelated text about quantum physics"
    t4 = "THE QUICK BROWN FOX JUMPS OVER THE LAZY DOG"

    ht1 = th.hash_text(t1)
    ht2 = th.hash_text(t2)
    ht3 = th.hash_text(t3)
    ht4 = th.hash_text(t4)

    print(f"  t1 = '{t1[:30]}...'")
    print(f"  t2 = '{t2[:30]}...' (similar)")
    print(f"  t3 = '{t3[:30]}...' (different)")
    print(f"  t4 = '{t4[:30]}...' (uppercase)")
    print(f"  dist(t1,t2) = {th.distance(ht1,ht2):.3f}")
    print(f"  dist(t1,t3) = {th.distance(ht1,ht3):.3f}")
    print(f"  dist(t1,t4) = {th.distance(ht1,ht4):.3f}")

    # 3. Near-duplicate detection
    print("\n3. NEAR-DUPLICATE DETECTION")
    documents = [
        "machine learning is a subset of artificial intelligence",
        "machine learning is a part of artificial intelligence",  # near-dup
        "deep learning uses neural networks for pattern recognition",
        "neural networks are used in deep learning for patterns",  # near-dup
        "quantum computing leverages quantum mechanics for computation",
        "the weather today is sunny with clear skies",
        "the weather today is sunny with some clouds",  # near-dup
    ]

    th.build_index(documents)
    print(f"  Indexed {len(documents)} documents")
    print(f"\n  Pairwise distances:")
    for i in range(len(documents)):
        for j in range(i + 1, len(documents)):
            dist = th.distance(
                th.hash_text(documents[i]),
                th.hash_text(documents[j])
            )
            if dist < 0.3:
                print(f"    [{i}] vs [{j}]: dist={dist:.3f} NEAR-DUPLICATE")

    # 4. Byte-level hashing
    print("\n4. BYTE-LEVEL HASHING (structural similarity)")
    data1 = bytes([i % 256 for i in range(1024)])
    data2 = bytes([(i + 1) % 256 for i in range(1024)])  # shifted
    data3 = bytes([random.randint(0, 255) for _ in range(1024)])  # random
    data4 = bytes([0] * 1024)  # all zeros

    hb1 = th.hash_bytes(data1)
    hb2 = th.hash_bytes(data2)
    hb3 = th.hash_bytes(data3)
    hb4 = th.hash_bytes(data4)

    print(f"  counter vs shifted:  dist={th.distance(hb1, hb2):.3f}")
    print(f"  counter vs random:   dist={th.distance(hb1, hb3):.3f}")
    print(f"  counter vs zeros:    dist={th.distance(hb1, hb4):.3f}")
    print(f"  random vs zeros:     dist={th.distance(hb3, hb4):.3f}")

    # 5. Nearest neighbor search
    print("\n5. NEAREST NEIGHBOR SEARCH")
    # Generate clustered data
    vectors = []
    labels = {}
    for cluster in range(3):
        center = [cluster * 10 + random.gauss(0, 1) for _ in range(8)]
        for i in range(20):
            v = [c + random.gauss(0, 2) for c in center]
            item_id = cluster * 20 + i
            vectors.append(v)
            labels[item_id] = cluster

    th2 = TournamentHasher(n_vertices=8)
    th2.build_index(vectors)

    # Query
    query = vectors[0]
    neighbors = th2.query(query, k=5)
    print(f"  Query: item 0 (cluster {labels[0]})")
    print(f"  Nearest neighbors:")
    for item_id, dist in neighbors:
        print(f"    item {item_id:3d} (cluster {labels[item_id]}): dist={dist:.3f}")

    # Cluster quality
    quality = th2.cluster_quality(labels)
    print(f"  Cluster quality: {quality:.3f} (higher = better separation)")

    # 6. Hash statistics
    print("\n6. TOURNAMENT HASH PROPERTIES")
    n = 8
    n_bits = n * (n - 1) // 2
    print(f"  Vertices: {n}")
    print(f"  Bits per hash: {n_bits}")
    print(f"  Possible tournaments: 2^{n_bits} = {2**n_bits:,}")
    print(f"  Hash size: {(n_bits + 7) // 8 + 1} bytes")
    print(f"  Distance range: [0, {n_bits}] arcs")

    # Show that transitive tournaments are the "sorted" case
    sorted_v = list(range(n))
    h_sorted = th.hash(sorted_v)
    print(f"\n  Sorted vector hash:  c3={h_sorted.c3_count()}, scores={h_sorted.score_sequence()}")
    print(f"  (Transitive = perfectly sorted data)")

    # c3 distribution for random data
    c3_counts = []
    for _ in range(1000):
        rv = [random.random() for _ in range(n)]
        rh = th.hash(rv)
        c3_counts.append(rh.c3_count())
    avg_c3 = sum(c3_counts) / len(c3_counts)
    max_c3 = max(c3_counts)
    print(f"  Random data avg c3: {avg_c3:.1f} (max possible = {n*(n-1)*(n-2)//6 - sum(k*(k-1)//2 for k in range(n))})")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=f'tournament_hash v{__version__}')
    parser.add_argument('--demo', action='store_true')
    parser.add_argument('--hash-file', help='Hash a file')
    parser.add_argument('--hash-text', help='Hash a text string')
    parser.add_argument('--vertices', '-n', type=int, default=8, help='Tournament size')
    args = parser.parse_args()

    if args.demo:
        demo()
    elif args.hash_file:
        th = TournamentHasher(n_vertices=args.vertices)
        with open(args.hash_file, 'rb') as f:
            h = th.hash_bytes(f.read())
        print(f"{args.hash_file}: {h.hex()}")
        print(f"  scores={h.score_sequence()}, c3={h.c3_count()}")
    elif args.hash_text:
        th = TournamentHasher(n_vertices=args.vertices)
        h = th.hash_text(args.hash_text)
        print(f"text: {h.hex()}")
        print(f"  scores={h.score_sequence()}, c3={h.c3_count()}")
    else:
        parser.print_help()
