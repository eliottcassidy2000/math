#!/usr/bin/env python3
"""cycle_counts_spectral_boundary_kps_S128c136.py -- kind-pasteur-2026-07-21-S128c136

RECONNECTING a forgotten thread to the live frame.  THM-172 (opus, era-2, forgotten: cited by
1 file) proved the total directed 5-CYCLE count c5 is LAMBDA-DETERMINED (a function of the
spectrum).  THM-1780 (mine, current) proved H = Hamiltonian-path count is NOT spectral -- it
splits within a co-spectral class at n = 6.  Between them is a boundary: for which cycle length
k is the simple directed k-cycle count c_k spectral (constant on co-spectral classes) and where
does it split?  On the moment-nullcone / binary-form ladder (THM-1775/1810), 'spectral' = 'a
function of the trace moments = the char_A coefficients', so this maps exactly where the cycle
statistics sit on the ladder.

Computes c_k (simple directed k-cycles, k=3..n) over all tournaments, groups by the trace-moment
vector (= char poly), and reports the first (n, k) where c_k SPLITS -- the spectral boundary.
"""
import sys
import numpy as np
from itertools import permutations, combinations
from collections import defaultdict


def from_bits(bits, n):
    A = np.zeros((n, n), dtype=np.int64)
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if bits >> idx & 1:
                A[i, j] = 1
            else:
                A[j, i] = 1
            idx += 1
    return A


def moment_vec(A, n):
    v = []
    P = np.eye(n, dtype=np.int64)
    for k in range(1, n + 1):
        P = P @ A
        v.append(int(np.trace(P)))
    return tuple(v)


def simple_cycle_counts(A, n, kmax):
    """c_k = number of simple directed k-cycles, for k = 3..kmax.  Enumerate k-subsets and
    cyclic orderings (each directed cycle counted once)."""
    c = {k: 0 for k in range(3, kmax + 1)}
    for k in range(3, kmax + 1):
        for S in combinations(range(n), k):
            # fix S[0] as start to avoid k-fold rotation overcount; count directed cycles
            rest = S[1:]
            for perm in permutations(rest):
                seq = (S[0],) + perm
                ok = all(A[seq[i], seq[(i + 1) % k]] for i in range(k))
                if ok:
                    c[k] += 1
        c[k] //= 1   # each cycle counted twice? start fixed, direction free -> counted once per direction
        # with start fixed and permuting rest, each directed cycle appears exactly once
    return c


print("=" * 84)
print("Is c_k (simple directed k-cycle count) SPECTRAL? split within co-spectral classes?")
print("=" * 84)
for n in (4, 5, 6):
    m = n * (n - 1) // 2
    kmax = n
    groups = defaultdict(lambda: defaultdict(set))   # mv -> k -> set of c_k values
    for bits in range(1 << m):
        A = from_bits(bits, n)
        mv = moment_vec(A, n)
        cc = simple_cycle_counts(A, n, kmax)
        for k in range(3, kmax + 1):
            groups[mv][k].add(cc[k])
    split_k = {}
    for mv, kd in groups.items():
        for k, vals in kd.items():
            if len(vals) > 1:
                split_k.setdefault(k, mv)
    line = []
    for k in range(3, kmax + 1):
        line.append("c%d:%s" % (k, "SPLIT" if k in split_k else "spectral"))
    print("  n=%d : %s" % (n, "  ".join(line)))
    for k in sorted(split_k):
        print("       c%d first splits at moment-vector %s" % (k, split_k[k]))
    sys.stdout.flush()

print()
print("  READING vs the ladder:")
print("  c3 = tr(A^3)/3 is always spectral (a trace moment).  THM-172: c5 spectral.  H (the")
print("  n-cycle/path level) splits at n=6 (THM-1780).  Any c_k that SPLITS above is the first")
print("  cycle statistic to leave the rational/spectral floor -- the boundary between the")
print("  forgotten cycle-count thread (THM-171/172/173) and the #P permanent (H).")
