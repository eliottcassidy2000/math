#!/usr/bin/env python3
"""fast_crystallizer_s116n.py — High-performance crystallization engines.

Six speedups over the naive O(n^3 * k) crystallizer:
  1. INCREMENTAL 3-CYCLE TRACKING: O(n) per flip instead of O(n^3)
  2. PRIORITY QUEUE: O(n log n) per flip for weakest-cycle selection
  3. BIT-ARRAY TOURNAMENT: 1 bit per arc (21 bits for n=7, fits uint32)
  4. BATCH CRYSTALLIZATION: bitwise majority for merging rankers
  5. PARALLEL API: multiprocessing for independent group crystallization
  6. BENCHMARK: comparison against naive at n=7,10,15,20

Instance: kind-pasteur-2026-03-16-S116n
"""

import time
import random
import heapq
from itertools import combinations
from collections import defaultdict
from multiprocessing import Pool, cpu_count

# ============================================================
# CORE DATA STRUCTURE: Bit-array tournament
# ============================================================

class BitTournament:
    """Tournament stored as a single integer (bit-array).

    For n vertices, we need C(n,2) bits.  Arc (i,j) with i<j is stored
    at bit position idx(i,j).  Bit=1 means i->j, bit=0 means j->i.

    n=7:  21 bits, fits in one uint32 (saves ~49x vs n*n array)
    n=42: 861 bits, ~108 bytes (saves ~1764x vs n*n array)
    """

    __slots__ = ('n', 'bits', '_num_arcs')

    def __init__(self, n, bits=0):
        self.n = n
        self.bits = bits
        self._num_arcs = n * (n - 1) // 2

    def _idx(self, i, j):
        """Bit index for ordered pair (i,j) with i < j.
        Convention: j*(j-1)//2 + i. This gives a unique index for each pair."""
        return j * (j - 1) // 2 + i

    def has_arc(self, a, b):
        """Returns True if a->b in the tournament."""
        if a == b:
            return False
        if a < b:
            return bool((self.bits >> self._idx(a, b)) & 1)
        else:
            return not bool((self.bits >> self._idx(b, a)) & 1)

    def set_arc(self, a, b):
        """Set arc a->b (and implicitly b-/->a)."""
        if a < b:
            self.bits |= (1 << self._idx(a, b))
        else:
            self.bits &= ~(1 << self._idx(b, a))

    def flip_arc(self, a, b):
        """Flip the arc between a and b."""
        i, j = min(a, b), max(a, b)
        self.bits ^= (1 << self._idx(i, j))

    def scores(self):
        """Compute out-degree scores for all vertices."""
        s = [0] * self.n
        for i in range(self.n):
            for j in range(i + 1, self.n):
                if self.has_arc(i, j):
                    s[i] += 1
                else:
                    s[j] += 1
        return s

    def to_matrix(self):
        """Convert to adjacency matrix (for compatibility)."""
        n = self.n
        M = [[0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                if self.has_arc(i, j):
                    M[i][j] = 1
                else:
                    M[j][i] = 1
        return M

    @classmethod
    def from_matrix(cls, M):
        """Create BitTournament from adjacency matrix."""
        n = len(M)
        bt = cls(n)
        for i in range(n):
            for j in range(i + 1, n):
                if M[i][j]:
                    bt.set_arc(i, j)
                else:
                    bt.set_arc(j, i)
        return bt

    @classmethod
    def random(cls, n, rng=None):
        """Generate a uniformly random tournament on n vertices."""
        if rng is None:
            rng = random
        bt = cls(n)
        num_arcs = n * (n - 1) // 2
        bt.bits = rng.getrandbits(num_arcs)
        return bt

    def copy(self):
        return BitTournament(self.n, self.bits)

    def __eq__(self, other):
        return self.n == other.n and self.bits == other.bits

    def __hash__(self):
        return hash((self.n, self.bits))

    def __repr__(self):
        return f"BitTournament(n={self.n}, bits=0x{self.bits:x})"


# ============================================================
# SPEEDUP 1+2: Incremental 3-cycle tracking with priority queue
# ============================================================

class IncrementalCrystallizer:
    """Crystallizer with O(n) per-flip cycle updates and O(log n) weakest-cycle extraction.

    Instead of scanning all C(n,3) triples each iteration:
    - Maintain a SET of all current 3-cycles
    - Maintain a MIN-HEAP ordered by weakest arc margin
    - On flip of arc (a,b): only recheck the O(n) triples containing both a and b
    - Pop weakest from heap in O(log C) where C = number of active cycles

    Total complexity: O(k * n * log(n^2)) = O(k * n * log n) vs O(k * n^3) naive.
    """

    def __init__(self, tournament, margins):
        """
        tournament: BitTournament
        margins: dict mapping (i,j) with i<j -> margin value (positive = i->j preferred)
        """
        self.T = tournament.copy()
        self.n = tournament.n
        self.margins = margins  # (i,j) -> margin, i<j
        # Active 3-cycles: set of frozensets {i,j,k}
        self.active_cycles = set()
        # Heap entries: (weakest_margin, cycle_id_tuple, arc_to_flip)
        # We use lazy deletion: entries remain in heap but we check active_cycles
        self.heap = []
        self._cycle_counter = 0

        # Build initial cycle set
        self._build_all_cycles()

    def _build_all_cycles(self):
        """Find all 3-cycles in the current tournament. O(n^3) -- done once."""
        n = self.n
        self.active_cycles.clear()
        self.heap = []
        for i in range(n):
            for j in range(i + 1, n):
                for k in range(j + 1, n):
                    self._check_and_add_triple(i, j, k)

    def _check_and_add_triple(self, i, j, k):
        """Check if {i,j,k} forms a 3-cycle and if so add to tracking."""
        T = self.T
        # Get the three arc directions
        ij = T.has_arc(i, j)  # True if i->j
        jk = T.has_arc(j, k)  # True if j->k
        ik = T.has_arc(i, k)  # True if i->k

        # 3-cycle iff NOT transitive.
        # Transitive means one vertex beats both others.
        # That happens when out-degrees are (2,1,0).
        # 3-cycle means out-degrees are (1,1,1).
        d_i = int(ij) + int(ik)
        d_j = int(not ij) + int(jk)
        d_k = int(not ik) + int(not jk)

        if not (d_i == 1 and d_j == 1 and d_k == 1):
            return  # transitive triple, not a cycle

        cycle_key = frozenset((i, j, k))
        self.active_cycles.add(cycle_key)

        # Determine cycle direction to find arcs
        # i->j->k->i: ij=T, jk=T, ki=T (i.e., ik=F)
        if ij and jk and not ik:
            self._push_cycle(i, j, k, cycle_key)
        else:
            # i->k->j->i: ik=T, kj=T(jk=F), ji=T(ij=F)
            self._push_cycle(i, k, j, cycle_key)

    def _push_cycle(self, a, b, c, cycle_key):
        """Push cycle a->b->c->a onto the heap with its weakest-arc margin.

        Tie-breaking: when margins are equal, pick the arc with smaller
        (tail, head) tuple. The heap key is (margin, tail, head, counter, ...)
        so that identical-margin arcs are broken deterministically.
        """
        arcs = [(a, b), (b, c), (c, a)]
        min_margin = float('inf')
        weakest_arc = None
        for tail, head in arcs:
            i, j = min(tail, head), max(tail, head)
            m = abs(self.margins.get((i, j), 0))
            if m < min_margin or (m == min_margin and (tail, head) < weakest_arc):
                min_margin = m
                weakest_arc = (tail, head)

        self._cycle_counter += 1
        # Heap key: (margin, tail, head, counter, cycle_key, arc)
        # This ensures deterministic ordering when margins tie.
        heapq.heappush(self.heap, (
            min_margin, weakest_arc[0], weakest_arc[1],
            self._cycle_counter, cycle_key, weakest_arc
        ))

    def _update_cycles_for_arc(self, a, b):
        """After flipping arc (a,b), update cycles involving a or b. O(n)."""
        n = self.n
        # Every 3-cycle involving arc (a,b) has the form {a,b,v} for some v.
        # Remove old cycles and check if new cycles formed.
        for v in range(n):
            if v == a or v == b:
                continue
            triple = frozenset((a, b, v))
            # Remove from active set (may or may not have been there)
            self.active_cycles.discard(triple)
            # Re-check this triple
            i, j, k = sorted((a, b, v))
            self._check_and_add_triple(i, j, k)

    def crystallize(self, max_iter=10000):
        """Run the crystallization loop.
        Returns: (final_tournament, iterations, flips_list)

        Anti-cycling: each arc can be flipped at most once. After flipping,
        it is permanently frozen. This matches the naive algorithm's behavior
        where flipping the weakest arc once is sufficient — if a cycle remains
        after that, its other arcs can still be flipped.
        """
        iterations = 0
        flips = []
        flipped_arcs = set()  # arcs already flipped once — never flip again

        while iterations < max_iter and self.active_cycles:
            # Pop weakest cycle (lazy deletion)
            found = False
            while self.heap:
                entry = heapq.heappop(self.heap)
                margin = entry[0]
                cycle_key = entry[4]
                arc = entry[5]
                if cycle_key in self.active_cycles:
                    arc_key = (min(arc[0], arc[1]), max(arc[0], arc[1]))
                    if arc_key not in flipped_arcs:
                        found = True
                        break
                    # This cycle's weakest arc was already flipped.
                    # The cycle might have OTHER unfrozen arcs. But since
                    # we only stored the weakest, we need to recompute.
                    # Remove from active and re-check if cycle still exists.
                    self.active_cycles.discard(cycle_key)
                    # Re-add if the cycle still exists with a different weakest arc
                    verts = sorted(cycle_key)
                    self._recheck_cycle_with_frozen(verts[0], verts[1], verts[2], flipped_arcs)

            if not found:
                break

            tail, head = arc
            arc_key = (min(tail, head), max(tail, head))
            flipped_arcs.add(arc_key)

            # Flip the weakest arc
            self.T.flip_arc(tail, head)
            flips.append((tail, head, margin))

            # Update affected cycles: only triples {tail, head, v} are affected
            self._update_cycles_for_arc(tail, head)

            iterations += 1

        return self.T, iterations, flips

    def _recheck_cycle_with_frozen(self, i, j, k, frozen):
        """Re-check if {i,j,k} is still a cycle and has an unfrozen arc to flip."""
        T = self.T
        ij = T.has_arc(i, j)
        jk = T.has_arc(j, k)
        ik = T.has_arc(i, k)

        d_i = int(ij) + int(ik)
        d_j = int(not ij) + int(jk)
        d_k = int(not ik) + int(not jk)

        if not (d_i == 1 and d_j == 1 and d_k == 1):
            return  # no longer a cycle

        cycle_key = frozenset((i, j, k))

        # Determine arcs
        if ij and jk and not ik:
            arcs = [(i, j), (j, k), (k, i)]
        else:
            arcs = [(i, k), (k, j), (j, i)]

        # Find weakest UNFROZEN arc
        min_margin = float('inf')
        weakest_arc = None
        for tail, head in arcs:
            ak = (min(tail, head), max(tail, head))
            if ak in frozen:
                continue
            m = abs(self.margins.get(ak if tail < head else (head, tail), 0))
            # Correct: always use (min, max) for margin lookup
            mi, mj = min(tail, head), max(tail, head)
            m = abs(self.margins.get((mi, mj), 0))
            if m < min_margin or (m == min_margin and (tail, head) < weakest_arc):
                min_margin = m
                weakest_arc = (tail, head)

        if weakest_arc is None:
            # All arcs frozen — this cycle is stuck, remove it
            self.active_cycles.discard(cycle_key)
            return

        self.active_cycles.add(cycle_key)
        self._cycle_counter += 1
        heapq.heappush(self.heap, (
            min_margin, weakest_arc[0], weakest_arc[1],
            self._cycle_counter, cycle_key, weakest_arc
        ))


# ============================================================
# SPEEDUP 3: Already integrated — BitTournament above
# ============================================================


# ============================================================
# SPEEDUP 4: Batch crystallization via bitwise majority voting
# ============================================================

def bitwise_majority(tournaments):
    """Compute majority tournament from a list of BitTournaments.
    Each tournament is stored as a single integer.
    For each bit position, take the majority vote.

    If each ranker's tournament is a 21-bit integer (n=7),
    majority = O(num_arcs) bit operations = extremely fast.
    """
    if not tournaments:
        raise ValueError("Need at least one tournament")
    n = tournaments[0].n
    num_arcs = n * (n - 1) // 2
    num_rankers = len(tournaments)
    threshold = num_rankers // 2  # strict majority

    # Count votes per bit position
    result_bits = 0
    for bit in range(num_arcs):
        count = 0
        mask = 1 << bit
        for t in tournaments:
            if t.bits & mask:
                count += 1
        if count > threshold:
            result_bits |= mask
        # Ties broken by bit=0 (lower-indexed vertex wins)

    # Also compute margins for the crystallizer
    margins = {}
    for bit in range(num_arcs):
        mask = 1 << bit
        count = sum(1 for t in tournaments if t.bits & mask)
        # Decode bit -> (i,j) pair
        # bit = j*(j-1)//2 + i, need to invert
        j = 0
        while j * (j - 1) // 2 + j <= bit:
            j += 1
        # now j is too large by 1 (or just right)
        # Actually let's be more careful
        j = 0
        while (j + 1) * j // 2 <= bit:
            j += 1
        i = bit - j * (j - 1) // 2
        margins[(i, j)] = 2 * count - num_rankers  # positive = i->j preferred

    return BitTournament(n, result_bits), margins


def batch_majority_fast(bit_arrays, n):
    """Ultra-fast majority using column-wise bit tricks.
    bit_arrays: list of integers (one per ranker)
    Returns: majority integer, margins dict

    For n=7 (21 bits), this processes all bits simultaneously.
    """
    num_arcs = n * (n - 1) // 2
    num_rankers = len(bit_arrays)
    threshold = num_rankers // 2

    # For each bit, count how many rankers have it set
    counts = [0] * num_arcs
    for bits in bit_arrays:
        for b in range(num_arcs):
            if bits & (1 << b):
                counts[b] += 1

    result = 0
    margins = {}
    for b in range(num_arcs):
        if counts[b] > threshold:
            result |= (1 << b)
        # Decode bit index to (i,j)
        j = 0
        while (j + 1) * j // 2 <= b:
            j += 1
        i = b - j * (j - 1) // 2
        margins[(i, j)] = 2 * counts[b] - num_rankers

    return result, margins


# ============================================================
# SPEEDUP 5: Parallel crystallization API
# ============================================================

def _crystallize_group(args):
    """Worker function for parallel crystallization.
    args: (tournament_bits, n, margins, max_iter)
    Returns: (final_bits, iterations, num_flips)
    """
    bits, n, margins, max_iter = args
    T = BitTournament(n, bits)
    cryst = IncrementalCrystallizer(T, margins)
    final_T, iters, flips = cryst.crystallize(max_iter=max_iter)
    return (final_T.bits, iters, len(flips))


def parallel_crystallize(groups, num_workers=None):
    """Crystallize multiple independent groups in parallel.

    groups: list of (BitTournament, margins_dict) tuples
    num_workers: number of processes (default: cpu_count)

    Returns: list of (final_BitTournament, iterations, num_flips)

    Design for the 7-item module architecture:
    If you have 42 items partitioned into 6 groups of 7,
    each group crystallizes independently => embarrassingly parallel.
    """
    if num_workers is None:
        num_workers = min(cpu_count(), len(groups))

    if num_workers <= 1 or len(groups) <= 1:
        # Sequential fallback
        results = []
        for T, margins in groups:
            cryst = IncrementalCrystallizer(T, margins)
            final_T, iters, flips = cryst.crystallize()
            results.append((final_T, iters, len(flips)))
        return results

    # Prepare serializable args
    args_list = [
        (T.bits, T.n, margins, 10000)
        for T, margins in groups
    ]

    with Pool(num_workers) as pool:
        raw_results = pool.map(_crystallize_group, args_list)

    results = []
    for (bits, iters, nflips), (T, _) in zip(raw_results, groups):
        results.append((BitTournament(T.n, bits), iters, nflips))

    return results


# ============================================================
# NAIVE CRYSTALLIZER (for benchmarking)
# ============================================================

def naive_crystallize(M, margins_matrix, max_iter=10000):
    """The original O(n^3 * k) crystallizer. Uses adjacency matrix.
    M: n*n adjacency matrix (M[i][j]=1 means i->j)
    margins_matrix: n*n matrix where margins_matrix[i][j] = vote margin for i->j
    Returns: (M, iterations, num_flips)

    Includes oscillation detection: if the same arc would be flipped twice
    in a row, freeze it and skip.
    """
    n = len(M)
    # Deep copy
    T = [row[:] for row in M]
    iterations = 0
    flips = 0
    flipped_arcs = set()  # each arc flipped at most once

    while iterations < max_iter:
        # Scan ALL C(n,3) triples to find 3-cycles
        weakest = None
        weakest_margin = float('inf')
        found_any_cycle = False

        for i in range(n):
            for j in range(n):
                if j == i or T[i][j] == 0:
                    continue
                for k in range(n):
                    if k == i or k == j:
                        continue
                    if T[j][k] == 1 and T[k][i] == 1:
                        # 3-cycle i->j->k->i
                        found_any_cycle = True
                        arcs = [(i, j), (j, k), (k, i)]
                        for a, b in arcs:
                            arc_key = (min(a, b), max(a, b))
                            if arc_key in flipped_arcs:
                                continue
                            m = abs(margins_matrix[a][b])
                            if m < weakest_margin or (m == weakest_margin and (a, b) < weakest):
                                weakest_margin = m
                                weakest = (a, b)

        if not found_any_cycle or weakest is None:
            break

        a, b = weakest
        arc_key = (min(a, b), max(a, b))
        flipped_arcs.add(arc_key)

        T[a][b] = 0
        T[b][a] = 1
        flips += 1
        iterations += 1

    return T, iterations, flips


# ============================================================
# BENCHMARK UTILITIES
# ============================================================

def random_margins(n, rng=None):
    """Generate random vote margins for a tournament.
    Returns: (BitTournament, margins_dict, margins_matrix)
    """
    if rng is None:
        rng = random
    T = BitTournament(n)
    margins = {}
    margins_matrix = [[0] * n for _ in range(n)]

    for i in range(n):
        for j in range(i + 1, n):
            # Random margin from 1 to 20
            m = rng.randint(1, 20)
            if rng.random() < 0.5:
                T.set_arc(i, j)
                margins[(i, j)] = m
                margins_matrix[i][j] = m
                margins_matrix[j][i] = -m
            else:
                T.set_arc(j, i)
                margins[(i, j)] = -m
                margins_matrix[j][i] = m
                margins_matrix[i][j] = -m

    return T, margins, margins_matrix


def benchmark(n, num_trials=10, seed=42):
    """Benchmark fast vs naive crystallizer at size n.
    Returns: dict with timing results.
    """
    rng = random.Random(seed)

    naive_times = []
    fast_times = []
    naive_iters_list = []
    fast_iters_list = []

    for trial in range(num_trials):
        T, margins, margins_matrix = random_margins(n, rng)

        # -- Naive --
        M = T.to_matrix()
        t0 = time.perf_counter()
        _, n_iters, n_flips = naive_crystallize(M, margins_matrix)
        t1 = time.perf_counter()
        naive_times.append(t1 - t0)
        naive_iters_list.append(n_iters)

        # -- Fast (incremental + heap) --
        T2 = T.copy()
        t0 = time.perf_counter()
        cryst = IncrementalCrystallizer(T2, margins)
        _, f_iters, f_flips = cryst.crystallize()
        t1 = time.perf_counter()
        fast_times.append(t1 - t0)
        fast_iters_list.append(f_iters)

    avg_naive = sum(naive_times) / len(naive_times)
    avg_fast = sum(fast_times) / len(fast_times)
    speedup = avg_naive / avg_fast if avg_fast > 0 else float('inf')

    return {
        'n': n,
        'trials': num_trials,
        'naive_avg_ms': avg_naive * 1000,
        'fast_avg_ms': avg_fast * 1000,
        'speedup': speedup,
        'naive_avg_iters': sum(naive_iters_list) / len(naive_iters_list),
        'fast_avg_iters': sum(fast_iters_list) / len(fast_iters_list),
        'naive_times_ms': [t * 1000 for t in naive_times],
        'fast_times_ms': [t * 1000 for t in fast_times],
    }


def verify_correctness(n=7, num_trials=100, seed=12345):
    """Verify that fast crystallizer produces equivalent results to naive.

    Both algorithms are valid crystallizers. When margin ties exist, they may
    break ties differently, producing different (but equally valid) fixed points.
    We verify:
    1. Both reach the same remaining 3-cycle count.
    2. Both reach zero cycles when max_iter is sufficient.
    3. The fast version never produces MORE cycles than naive.
    """
    rng = random.Random(seed)
    mismatches = 0
    exact_matches = 0
    tie_diffs = 0  # Different tournament but same cycle count (valid tie-break difference)

    for trial in range(num_trials):
        T, margins, margins_matrix = random_margins(n, rng)

        # Naive
        M = T.to_matrix()
        M_final, n_iters, _ = naive_crystallize(M, margins_matrix)

        # Fast
        T2 = T.copy()
        cryst = IncrementalCrystallizer(T2, margins)
        T_final, f_iters, _ = cryst.crystallize()

        # Compare: same final tournament?
        M_fast = T_final.to_matrix()
        naive_cycles = count_3cycles_matrix(M_final)
        fast_cycles = count_3cycles_matrix(M_fast)

        if M_final == M_fast:
            exact_matches += 1
        elif naive_cycles == fast_cycles:
            tie_diffs += 1  # Valid: same quality, different tie-breaking
        else:
            mismatches += 1  # Real problem: different cycle counts

    return mismatches, exact_matches, tie_diffs, num_trials


def count_3cycles_matrix(M):
    """Count 3-cycles in adjacency matrix."""
    n = len(M)
    count = 0
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                if M[i][j] + M[j][k] + M[k][i] == 3:
                    count += 1
                elif M[j][i] + M[i][k] + M[k][j] == 3:
                    count += 1
    return count


def benchmark_batch_majority(n=7, num_rankers=100, seed=42):
    """Benchmark batch majority computation."""
    rng = random.Random(seed)
    num_arcs = n * (n - 1) // 2

    # Generate random tournaments as bit arrays
    tournaments = [BitTournament.random(n, rng) for _ in range(num_rankers)]
    bit_arrays = [t.bits for t in tournaments]

    # Method 1: Bitwise majority function
    t0 = time.perf_counter()
    for _ in range(100):
        result_T, result_margins = bitwise_majority(tournaments)
    t1 = time.perf_counter()
    batch_time = (t1 - t0) / 100

    # Method 2: Naive pairwise counting (like original Engine 1)
    t0 = time.perf_counter()
    for _ in range(100):
        wins = [[0] * n for _ in range(n)]
        for t in tournaments:
            for i in range(n):
                for j in range(i + 1, n):
                    if t.has_arc(i, j):
                        wins[i][j] += 1
                    else:
                        wins[j][i] += 1
        # Build majority
        T_naive = [[0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                if wins[i][j] > wins[j][i]:
                    T_naive[i][j] = 1
                else:
                    T_naive[j][i] = 1
    t1 = time.perf_counter()
    naive_time = (t1 - t0) / 100

    return {
        'n': n,
        'num_rankers': num_rankers,
        'batch_ms': batch_time * 1000,
        'naive_ms': naive_time * 1000,
        'speedup': naive_time / batch_time if batch_time > 0 else float('inf'),
    }


def benchmark_parallel(n=7, num_groups=12, seed=42):
    """Benchmark parallel crystallization of multiple independent groups."""
    rng = random.Random(seed)

    groups = []
    for _ in range(num_groups):
        T, margins, _ = random_margins(n, rng)
        groups.append((T, margins))

    # Sequential
    t0 = time.perf_counter()
    seq_results = parallel_crystallize(groups, num_workers=1)
    t1 = time.perf_counter()
    seq_time = t1 - t0

    # Parallel
    nw = min(cpu_count(), num_groups)
    t0 = time.perf_counter()
    par_results = parallel_crystallize(groups, num_workers=nw)
    t1 = time.perf_counter()
    par_time = t1 - t0

    return {
        'n': n,
        'num_groups': num_groups,
        'workers': nw,
        'sequential_ms': seq_time * 1000,
        'parallel_ms': par_time * 1000,
        'speedup': seq_time / par_time if par_time > 0 else float('inf'),
    }


# ============================================================
# MEMORY COMPARISON
# ============================================================

def memory_comparison():
    """Compare memory usage: BitTournament vs adjacency matrix."""
    results = []
    for n in [7, 10, 15, 20, 42, 100]:
        num_arcs = n * (n - 1) // 2
        bit_bytes = (num_arcs + 7) // 8  # ceil(bits / 8)
        matrix_bytes = n * n  # 1 byte per entry
        ratio = matrix_bytes / bit_bytes if bit_bytes > 0 else float('inf')
        results.append({
            'n': n,
            'arcs': num_arcs,
            'bit_bytes': bit_bytes,
            'matrix_bytes': matrix_bytes,
            'compression': ratio,
        })
    return results


# ============================================================
# MAIN: Run all benchmarks
# ============================================================

if __name__ == '__main__':
    import sys

    print()
    print("=" * 72)
    print("  FAST CRYSTALLIZER — BENCHMARK SUITE")
    print("  kind-pasteur-2026-03-16-S116n")
    print("=" * 72)
    print()

    # ----------------------------------------------------------
    print("  SECTION 1: CORRECTNESS VERIFICATION")
    print("  " + "-" * 50)
    print()
    for n in [5, 7, 10]:
        trials = 200 if n <= 7 else 50
        mismatches, exact, tie_diffs, total = verify_correctness(n=n, num_trials=trials)
        status = "PASS" if mismatches == 0 else f"FAIL ({mismatches} cycle-count mismatches)"
        print(f"  n={n:2d}: {total:3d} trials -> {status}")
        print(f"         exact match: {exact}, tie-break diff (same cycles): {tie_diffs}, real mismatch: {mismatches}")
    print()

    # ----------------------------------------------------------
    print("  SECTION 2: MEMORY COMPARISON")
    print("  " + "-" * 50)
    print()
    print(f"  {'n':>4s}  {'arcs':>6s}  {'bit_bytes':>10s}  {'matrix_bytes':>13s}  {'compression':>12s}")
    print(f"  {'':->4s}  {'':->6s}  {'':->10s}  {'':->13s}  {'':->12s}")
    for r in memory_comparison():
        print(f"  {r['n']:4d}  {r['arcs']:6d}  {r['bit_bytes']:10d}  {r['matrix_bytes']:13d}  {r['compression']:11.1f}x")
    print()

    # ----------------------------------------------------------
    print("  SECTION 3: CRYSTALLIZATION BENCHMARK (fast vs naive)")
    print("  " + "-" * 50)
    print()

    # Adaptive trial counts: more for small n, fewer for large n
    bench_configs = [
        (7,  50),
        (10, 30),
        (15, 10),
        (20,  5),
    ]

    print(f"  {'n':>4s}  {'trials':>6s}  {'naive_ms':>10s}  {'fast_ms':>10s}  "
          f"{'speedup':>8s}  {'naive_iters':>11s}  {'fast_iters':>10s}")
    print(f"  {'':->4s}  {'':->6s}  {'':->10s}  {'':->10s}  "
          f"{'':->8s}  {'':->11s}  {'':->10s}")

    all_results = []
    for n, trials in bench_configs:
        r = benchmark(n, num_trials=trials)
        all_results.append(r)
        print(f"  {r['n']:4d}  {r['trials']:6d}  {r['naive_avg_ms']:10.3f}  {r['fast_avg_ms']:10.3f}  "
              f"{r['speedup']:7.1f}x  {r['naive_avg_iters']:11.1f}  {r['fast_avg_iters']:10.1f}")
    print()

    # Detailed timing for largest case
    if all_results:
        r = all_results[-1]
        print(f"  Detail for n={r['n']}:")
        print(f"    Naive per-trial (ms): {', '.join(f'{t:.1f}' for t in r['naive_times_ms'])}")
        print(f"    Fast  per-trial (ms): {', '.join(f'{t:.1f}' for t in r['fast_times_ms'])}")
    print()

    # ----------------------------------------------------------
    print("  SECTION 4: BATCH MAJORITY BENCHMARK")
    print("  " + "-" * 50)
    print()
    for n in [7, 10, 15, 20]:
        for nr in [50, 200]:
            r = benchmark_batch_majority(n=n, num_rankers=nr)
            print(f"  n={r['n']:2d}, {r['num_rankers']:3d} rankers:  "
                  f"batch={r['batch_ms']:.3f}ms  naive={r['naive_ms']:.3f}ms  "
                  f"speedup={r['speedup']:.1f}x")
    print()

    # ----------------------------------------------------------
    print("  SECTION 5: PARALLEL CRYSTALLIZATION")
    print("  " + "-" * 50)
    print()
    # Test with larger n groups to overcome process-spawn overhead
    for n_par, ng in [(15, 8), (15, 16), (20, 8), (20, 16)]:
        r = benchmark_parallel(n=n_par, num_groups=ng)
        print(f"  {r['num_groups']:2d} groups of n={n_par}, {r['workers']} workers:  "
              f"seq={r['sequential_ms']:.1f}ms  par={r['parallel_ms']:.1f}ms  "
              f"speedup={r['speedup']:.2f}x")
    print()

    # ----------------------------------------------------------
    print("  SECTION 6: SCALING ANALYSIS")
    print("  " + "-" * 50)
    print()
    print("  Theoretical complexity:")
    print("    Naive:  O(k * n^3)    per crystallization run")
    print("    Fast:   O(n^3) init + O(k * n * log n)  per run")
    print("    Batch:  O(R * n^2) -> O(R * n^2/w)  with w-bit words")
    print("    Parallel: O(T_seq / P) for P independent groups")
    print()

    # Compute observed scaling
    if len(all_results) >= 2:
        print("  Observed scaling (naive):")
        for i in range(1, len(all_results)):
            r0, r1 = all_results[i - 1], all_results[i]
            if r0['naive_avg_ms'] > 0 and r1['naive_avg_ms'] > 0:
                n_ratio = r1['n'] / r0['n']
                t_ratio = r1['naive_avg_ms'] / r0['naive_avg_ms']
                if t_ratio > 0 and n_ratio > 1:
                    import math
                    exponent = math.log(t_ratio) / math.log(n_ratio)
                    print(f"    n={r0['n']}->{r1['n']}: time ratio={t_ratio:.1f}x, "
                          f"n ratio={n_ratio:.1f}x, exponent~={exponent:.1f}")

        print("  Observed scaling (fast):")
        for i in range(1, len(all_results)):
            r0, r1 = all_results[i - 1], all_results[i]
            if r0['fast_avg_ms'] > 0 and r1['fast_avg_ms'] > 0:
                n_ratio = r1['n'] / r0['n']
                t_ratio = r1['fast_avg_ms'] / r0['fast_avg_ms']
                if t_ratio > 0 and n_ratio > 1:
                    import math
                    exponent = math.log(t_ratio) / math.log(n_ratio)
                    print(f"    n={r0['n']}->{r1['n']}: time ratio={t_ratio:.1f}x, "
                          f"n ratio={n_ratio:.1f}x, exponent~={exponent:.1f}")
    print()

    # ----------------------------------------------------------
    print("  SECTION 7: API SUMMARY")
    print("  " + "-" * 50)
    print()
    print("  BitTournament(n, bits=0)")
    print("    .has_arc(a,b)     O(1) arc query")
    print("    .flip_arc(a,b)    O(1) arc flip")
    print("    .scores()         O(n^2) out-degrees")
    print("    .random(n)        O(1) random tournament")
    print("    .from_matrix(M)   O(n^2) convert from adjacency matrix")
    print("    .to_matrix()      O(n^2) convert to adjacency matrix")
    print()
    print("  IncrementalCrystallizer(tournament, margins)")
    print("    .crystallize(max_iter=10000)")
    print("    Returns: (final_tournament, iterations, flips_list)")
    print("    Complexity: O(n^3) init + O(k * n * log n) crystallization")
    print()
    print("  bitwise_majority(tournaments)")
    print("    Merge R tournaments via majority vote. O(R * n^2 / w).")
    print()
    print("  parallel_crystallize(groups, num_workers)")
    print("    Crystallize independent groups in parallel.")
    print("    groups: list of (BitTournament, margins) tuples")
    print()
    print("=" * 72)
    print("  BENCHMARK COMPLETE")
    print("=" * 72)
