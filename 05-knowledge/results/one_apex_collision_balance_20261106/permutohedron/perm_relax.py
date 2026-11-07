"""S84: sink spectra of partially coherent orientations of the permutohedron P_(n-1).

Edges of the permutohedron graph on S_n: pi -- pi s_i (swap positions i, i+1); the edge swaps the values {a,b}.
Orientation convention: the edge points INTO the arrangement where the 'first' element of the class precedes.
A vertex pi is a sink iff every adjacent pair (pi_i, pi_(i+1)) is 'forward' for its class.
Coherence levels (each defines a partition of the edges into classes oriented as one unit):
  L0 (parallel-coherent = tournaments): class = pair {a,b}.
  L1 (square-coherent):  class = (pair {a,b}, set S of elements before the pair).
  L2 (hexagon-coherent): class = (pair {a,b}, relative order of the other n-2 elements).
  FREE: class = edge.
A class orientation is a bit: bit 1 means min(a,b) before max(a,b) is forward.
"""
import itertools, sys
import numpy as np


def edges_and_classes(n, level):
    perms = list(itertools.permutations(range(n)))
    index = {p: k for k, p in enumerate(perms)}
    # each vertex/position gives a literal (class id, required bit)
    cls = {}
    lits = []  # per vertex: list of (class, bit)
    for p in perms:
        L = []
        for i in range(n - 1):
            a, b = p[i], p[i + 1]
            pair = (min(a, b), max(a, b))
            if level == 'L0':
                key = pair
            elif level == 'L1':
                key = (pair, frozenset(p[:i]))
            elif level == 'L2':
                key = (pair, tuple(x for x in p if x not in pair))
            elif level == 'FREE':
                q = list(p); q[i], q[i + 1] = q[i + 1], q[i]
                key = (pair, min(p, tuple(q)))
            else:
                raise ValueError(level)
            if key not in cls:
                cls[key] = len(cls)
            # p is a sink for this edge iff its order a-before-b is forward: bit 1 iff a < b
            L.append((cls[key], 1 if a < b else 0))
        lits.append(L)
    return perms, cls, lits


def spectrum_exhaustive(n, level, check_acyclic=False):
    perms, cls, lits = edges_and_classes(n, level)
    nc = len(cls)
    if nc > 26:
        raise ValueError('too many classes for exhaustion: %d' % nc)
    M = np.arange(1 << nc, dtype=np.int64)
    sinks = np.zeros(1 << nc, dtype=np.int16)
    for L in lits:
        ok = np.ones(1 << nc, dtype=bool)
        for c, bit in L:
            ok &= (((M >> c) & 1) == bit)
        sinks += ok
    vals, counts = np.unique(sinks, return_counts=True)
    return nc, dict(zip(vals.tolist(), counts.tolist())), sinks, perms, lits


if __name__ == '__main__':
    for n in (3, 4):
        for level in ('L0', 'L2', 'L1') if n == 4 else ('L0', 'L2', 'L1', 'FREE'):
            nc, spec, *_ = spectrum_exhaustive(n, level)
            print('n=%d %-4s classes=%2d  sink spectrum {count: #orientations} = %s' % (n, level, nc, spec))
