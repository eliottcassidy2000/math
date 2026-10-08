#!/usr/bin/env python3
"""audit H, task D: THM-4611 statement 4.  For every b != 0 mod d and k <= KMAX, the block at the hidden state
(M, e) = (1, b d^(k-1)) mod d^k equals sum over digit words (j_1..j_(k-1)) of D(j -> j + b prod m_(j_s) mod d)
(memoized residue recursion, itself checked against brute force).  On Z_4 (1,3,5,7), b = 2: A_k = 4^(k-1) D(j -> j+2), rank 2.
Also: the minimal block rank over ALL hidden states for Z_4 at k <= 6 (full tables), and the rank of D(j -> j+2)."""
import itertools, numpy as np
from hcore import Map
from hdp import Tables
for d, m, KMAX in [(4, [1, 3, 5, 7], 8), (5, [1, 2, 3, 7, 1], 6), (5, [1, 4, 1, 11, 34], 6), (7, [1, 1, 1, 1, 2, 3, 5], 4)]:
    mp = Map(d, m); bad = 0; n = 0
    for k in range(1, KMAX + 1):
        memo = {}
        for b in range(1, d):
            A = np.array(mp.A_rec(k, 1, (b * d ** (k - 1)) % d ** k, memo))
            ref = np.zeros_like(A)
            # sum over words of length k-1 of D(translation by b * prod m_j mod d)
            prods = {1 % d: 1}
            for _ in range(k - 1):
                new = {}
                for p, c in prods.items():
                    for j in range(d):
                        q = p * m[j] % d; new[q] = new.get(q, 0) + c
                prods = new
            for p, c in prods.items():
                ref += c * np.array(mp.Dmat(1, b * p % d))
            n += 1
            if not np.array_equal(A, ref): bad += 1
    print(f"Z_{d} m = {m}: identity-run blocks A_k(1, b d^(k-1)) = sum_paths D(translation by b prod m) for k <= {KMAX}, all b: "
          f"{n} cases, mismatches {bad}")
mp = Map(4, [1, 3, 5, 7])
D2 = np.array(mp.Dmat(1, 2))
print("Z_4: D(j -> j+2) =", D2.tolist(), "rank", np.linalg.matrix_rank(D2.astype(float)))
for k in range(1, 9):
    A = np.array(mp.A_rec(k, 1, 2 * 4 ** (k - 1) % 4 ** k, {}))
    print(f"   k = {k}: A_k(1, 2*4^(k-1)) == 4^(k-1) D(j -> j+2): {np.array_equal(A, 4 ** (k - 1) * D2)}")
tb = Tables(4, [1, 3, 5, 7])
for k in range(1, 7):
    T = tb.table(k).reshape(-1, 3, 3)
    nz = np.any(T.reshape(len(T), -1) != 0, axis=1)
    U = np.unique(T[nz].reshape(-1, 9), axis=0).reshape(-1, 3, 3)
    rk = np.linalg.matrix_rank(U.astype(float))
    print(f"   Z_4 full table k = {k}: {len(T)} states, {int((~nz).sum())} zero, {len(U)} distinct nonzero, min rank {rk.min()}, #rank<=2 distinct blocks {(rk <= 2).sum()}")
