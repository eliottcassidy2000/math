#!/usr/bin/env python3
# audit J1 (THM-4911 (3a)): the signed sub-Pfaffian vector z of an odd-order tournament's skew matrix M is in ker M, adj M = z z^T, every z_i is odd, sum z odd; plus u not in range(M) for all classes n <= 9.
import random, subprocess, sys
from functools import lru_cache
import numpy as np
import sympy as sp

GT = "/opt/homebrew/bin/gentourng"
def classes(n):
    m = n * (n - 1) // 2; iu = np.triu_indices(n, 1)
    out = subprocess.run([GT, "-q", str(n)], capture_output=True).stdout
    a = np.frombuffer(out, dtype=np.uint8).reshape(-1, m + 1)[:, :m] - 48
    A = np.zeros((len(a), n, n), dtype=np.int64); A[:, iu[0], iu[1]] = a; A[:, iu[1], iu[0]] = 1 - a
    return A

def pfaffian(M):
    Mt = tuple(tuple(int(v) for v in r) for r in M)
    @lru_cache(maxsize=None)
    def pf(idx):
        if not idx: return 1
        i0 = idx[0]; tot = 0
        for pos in range(1, len(idx)):
            j = idx[pos]; rest = idx[1:pos] + idx[pos + 1:]
            tot += (-1) ** (pos + 1) * Mt[i0][j] * pf(rest)
        return tot
    return pf(tuple(range(len(Mt))))

def check(M, full_adj):
    n = len(M)
    z = np.array([(-1) ** i * pfaffian(np.delete(np.delete(M, i, 0), i, 1)) for i in range(n)], dtype=object)  # 0-based: (-1)^i = (-1)^{(i+1)+1}
    ok_ker = all(v == 0 for v in M.astype(object).dot(z))
    ok_odd = all(int(v) % 2 == 1 for v in z) and int(sum(z)) % 2 == 1
    ok_adj = True
    if full_adj:
        adj = sp.Matrix(M.tolist()).adjugate()
        ok_adj = all(adj[i, j] == z[i] * z[j] for i in range(n) for j in range(n))
    return ok_ker, ok_odd, ok_adj

for n in (3, 5, 7):
    A = classes(n); M = A - A.transpose(0, 2, 1)
    res = [check(M[b], True) for b in range(len(A))]
    print(f"n={n}: all {len(A)} classes: z in ker M {sum(r[0] for r in res)}, z_i all odd & sum z odd {sum(r[1] for r in res)}, adj M = z z^T {sum(r[2] for r in res)}")
A = classes(9); M = A - A.transpose(0, 2, 1)
random.seed(9); sample = random.sample(range(len(A)), 400)
res = [check(M[b], False) for b in sample]
print(f"n=9: 400 random classes of {len(A)}: z in ker M {sum(r[0] for r in res)}, z_i all odd & sum z odd {sum(r[1] for r in res)}")
rng = np.random.default_rng(11); good = 0
for _ in range(60):
    n = 11; U = np.triu(2 * rng.integers(0, 2, (n, n)) - 1, 1); Mx = U - U.T      # random tournament, M_ij = +-1
    r = check(Mx, False); good += r[0] and r[1]
print(f"n=11: 60 random tournaments: kernel + parity OK in {good}/60")
# u not in range(M) for every class of odd order n <= 9 (rank test in floats; M has rank n-1, [M|u] rank n)
for n in (3, 5, 7, 9):
    A = classes(n); M = (A - A.transpose(0, 2, 1)).astype(np.float64)
    rM = np.linalg.matrix_rank(M); rMu = np.linalg.matrix_rank(np.concatenate([M, np.ones((len(A), n, 1))], axis=2))
    print(f"n={n}: rank M = n-1 for all classes: {(rM == n - 1).all()}; u not in range(M) for all classes: {(rMu == n).all()}")
print("(n-2)!! perfect matchings of K_(n-1) is odd for odd n: ", all(int(sp.factorial2(n - 2)) % 2 == 1 for n in range(3, 40, 2)))
