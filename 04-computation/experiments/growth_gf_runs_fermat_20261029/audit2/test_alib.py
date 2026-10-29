"""Self-test of alib against brute force (ordered set partitions x permutations inside blocks), n <= 6."""
import random
import numpy as np
from alib import (gentourng, adj_from_ascii, path_table, m2_matrix, m3_matrix, brute_mk, ham_cycles, is_strong,
                  snf_diag, snf_summary, fmt_snf, rank_mod_p, det_exact)
import itertools

random.seed(1)
# 1. path table vs brute force
for n in (3, 4, 5, 6):
    for s in gentourng(n)[:20]:
        A = adj_from_ascii(n, s)
        H = path_table(A)
        for U in range(1, 1 << n):
            verts = [v for v in range(n) if (U >> v) & 1]
            cnt = np.zeros((n, n), dtype=np.int64)
            for p in itertools.permutations(verts):
                if all(A[p[i], p[i + 1]] for i in range(len(p) - 1)):
                    cnt[p[0], p[-1]] += 1
            assert (cnt == H[:, :, U]).all(), (n, s, U)
print('path_table OK vs brute force (n=3..6, 20 classes each)')

# 2. M2, M3 vs brute force
for n in (3, 4, 5, 6):
    cls = gentourng(n)
    for s in cls[: (12 if n < 6 else 6)]:
        A = adj_from_ascii(n, s)
        M2, pr = m2_matrix(A)
        b2 = brute_mk(A, 2)
        for i, a in enumerate(pr):
            for j, b in enumerate(pr):
                assert M2[i, j] == b2.get((a, b), 0)
        assert sum(b2.values()) == M2.sum()
        if n >= 3:
            M3, bases, tup = m3_matrix(A)
            b3 = brute_mk(A, 3)
            for i, a in enumerate(tup):
                for j, b in enumerate(tup):
                    assert M3[i, j] == b3.get((a, b), 0), (n, s, a, b)
            assert sum(b3.values()) == M3.sum()
print('M2 and M3 OK vs brute force (n=3..6)')

# 3. HC count and strong vs brute force
for n in (3, 4, 5, 6, 7):
    for s in gentourng(n)[:60]:
        A = adj_from_ascii(n, s)
        H = path_table(A)
        hc = ham_cycles(A, H)
        bf = 0
        for p in itertools.permutations(range(1, n)):
            cyc = (0,) + p
            if all(A[cyc[i], cyc[(i + 1) % n]] for i in range(n)):
                bf += 1
        assert hc == bf
        assert is_strong(A) == (hc > 0)
print('HC counts OK vs brute force; strong <=> HC>0 (Camion) on samples n=3..7')

# 4. class counts
for n in range(3, 9):
    print('n=%d classes=%d' % (n, len(gentourng(n))))

# 5. SNF sanity
print(fmt_snf(snf_summary(snf_diag([[2, 4], [6, 8]]))), det_exact([[2, 4], [6, 8]]))
print(fmt_snf(snf_summary(snf_diag([[0, 0], [0, 0]]))))
print(fmt_snf(snf_summary(snf_diag([[1, 2, 3], [4, 5, 6], [7, 8, 9]]))))
