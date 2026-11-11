#!/usr/bin/env python3
"""H mod 4 of all circulant tournaments on Z_n by Theorem Q (folding only, no Hamiltonian-path count).
For each odd n: number of connection sets with H = 1 and H = 3 (mod 4), and the distribution of |D_1|.
Usage: circulant_fold_table.py NMAX"""
import sys, itertools
from collections import Counter
from fold import fold, paths_ending_in
nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 21
for n in range(3, nmax + 1, 2):
    half = list(range(1, (n + 1) // 2)); res = Counter(); byt = Counter()
    for signs in itertools.product((0, 1), repeat=len(half)):
        S = [h if s == 0 else n - h for h, s in zip(half, signs)]
        out = [0] * n
        for x in range(n):
            for s in S: out[x] |= 1 << ((x + s) % n)
        Z, D1, plus = fold(out, n, [(-x) % n for x in range(n)])
        r = (1 + 2 * paths_ending_in(Z, D1)) % 4
        res[r] += 1; byt[(bin(D1).count('1'), r)] += 1
    print("n = %2d: H = 1 (mod 4): %5d, H = 3 (mod 4): %5d;  by |D_1|: %s" % (n, res[1], res[3], " ".join("%d:%d/%d" % (t, byt[(t, 1)], byt[(t, 3)]) for t in range(len(half) + 1) if byt[(t, 1)] + byt[(t, 3)])))
    sys.stdout.flush()
