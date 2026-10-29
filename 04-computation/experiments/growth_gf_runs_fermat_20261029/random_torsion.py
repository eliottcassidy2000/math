"""HYP-9445 (1) beyond the census: random and low-cycle strong tournaments on 9 and 10 vertices, two-run observer.

Samples (seeded): uniformly random labelled tournaments, and 'low-cycle' tournaments found by a random walk of arc
flips from a transitive tournament with one back arc, kept when strong with 2 or 3 Hamiltonian cycles.
For each sample: number of Hamiltonian cycles, and whether coker M_2 has torsion (PARI matsnf).
usage: python3 random_torsion.py n n_random n_lowcycle seed
"""
import sys
import random
import itertools
import numpy as np
from fast_runs import m2_fast, ham_cycles, tuples2
from runs import gp_matrix, gp_batch, parse_vec


def random_tournament(n, rng):
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


def hc_count(A):
    _, hp = m2_fast(A) if False else (None, None)
    return None


def main():
    n, nr, nl, seed = map(int, sys.argv[1:5])
    rng = random.Random(seed)
    samples = []
    for _ in range(nr):
        samples.append(('random', random_tournament(n, rng)))
    # low-cycle walk
    A = [[1 if j > i else 0 for j in range(n)] for i in range(n)]
    A[0][n - 1], A[n - 1][0] = 0, 1          # one back arc: a unique Hamiltonian cycle
    found = 0
    steps = 0
    while found < nl and steps < 200000:
        steps += 1
        i, j = rng.sample(range(n), 2)
        A[i][j], A[j][i] = A[j][i], A[i][j]
        T, hp = m2_fast(A)
        hc = ham_cycles(A, hp)
        if hc == 0 or hc > 6:
            A[i][j], A[j][i] = A[j][i], A[i][j]    # reject: stay near the low-cycle region
            continue
        if hc in (2, 3):
            samples.append(('lowcycle', [row[:] for row in A]))
            found += 1
    tup = tuples2(n)
    cmds, meta = [], []
    for kind, A in samples:
        T, hp = m2_fast(A)
        hc = ham_cycles(A, hp)
        if hc == 0:
            continue
        M = [[int(T[s[0], s[1], t[0], t[1]]) for t in tup] for s in tup]
        cmds.append('matsnf(%s)' % gp_matrix(M))
        meta.append((kind, hc))
    out = []
    for i in range(0, len(cmds), 50):
        out += gp_batch(cmds[i:i + 50])
    stats = {}
    bad = []
    for (kind, hc), line in zip(meta, out):
        d = [abs(x) for x in parse_vec(line)]
        torsion = any(x > 1 for x in d)
        key = (kind, 'HC=1' if hc == 1 else 'HC>=2', 'torsion' if torsion else 'torsion-free')
        stats[key] = stats.get(key, 0) + 1
        if hc >= 2 and not torsion:
            bad.append((kind, hc, line[:120]))
    print('n=%d seed=%d samples=%d (strong): %s' % (n, seed, len(meta), sorted(stats.items())))
    print('counterexamples to HYP-9445 (1) at k=2:', bad if bad else 'none')


if __name__ == '__main__':
    main()
