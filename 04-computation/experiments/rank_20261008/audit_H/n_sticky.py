#!/usr/bin/env python3
"""audit H, task E: sticky reflections.
(a) every subset of Z/5 is invariant under some reflection j -> b - j (so every two-valued residue pattern is symmetric);
    and for d = 7, 11 the same statement FAILS (counts), so the "every {+-1} map has a sticky reflection" claim is Z_5-specific;
(b) for each of the six {+-1} maps and the Z_5 (1,2,3,7,1) family: from states with M = -1, e = b (mod 5) (b sticky), every
    digit leads to M' = -1 mod 5 (exact pair-chain steps on random exact states);
(c) on Z_5 every reflection covariance has rank <= 2, and under property (i) rank >= 1."""
import random, itertools
from fractions import Fraction as Fr
import numpy as np
from hcore import Map, residue
for d in (5, 7, 11):
    bad = [S for r in range(d + 1) for S in itertools.combinations(range(d), r)
           if not any({(b - j) % d for j in S} == set(S) for b in range(d))]
    print(f"d = {d}: subsets of Z/{d} not invariant under any reflection j -> b - j: {len(bad)} of {2 ** d}")
rnd = random.Random(3)
for m in ([1, 4, 1, 11, 34], [1, 6, 11, 11, 4], [1, 1, 6, 39, 11], [1, 14, 4, 1, 29], [1, 4, 31, 1, 6], [1, 6, 4, 31, 1], [1, 2, 3, 7, 1]):
    mp = Map(5, m); mb = [x % 5 for x in m]
    sticky = [b for b in range(5) if all(mb[(b - j) % 5] == mb[j] for j in range(5))]
    ok = True; ranks = []
    for b in range(5):
        ranks.append(int(np.linalg.matrix_rank(np.array(mp.Dmat(4, b), dtype=float))))
    for b in sticky:
        for _ in range(200):
            # random exact state with M = -1 mod 5, e = b mod 5
            while True:
                M = Fr(1)
                for i in range(1, 5): M *= Fr(m[i], m[0]) ** rnd.randint(-3, 3)
                if residue(M, 5) == 4: break
            e = Fr(5 * rnd.randint(-10 ** 5, 10 ** 5) + b, rnd.choice([1, 3, 7]))
            e = e + (b - residue(e, 5))   # force e = b mod 5
            assert residue(e, 5) == b
            for j in range(5):
                M2, e2, i = mp.step(M, e, j)
                ok &= residue(M2, 5) == 4
    print(f"m = {m}: residues {mb}, sticky b {sticky}, all digits keep a = -1 from (M,e) = (-1, b): {ok}; reflection covariance ranks by b: {ranks}")
