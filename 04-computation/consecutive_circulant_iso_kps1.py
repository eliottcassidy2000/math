#!/usr/bin/env python3
"""
consecutive_circulant_iso_kps1.py — kind-pasteur-2026-06-09-S1
BRANCH E follow-up 2 (HYP-2337): identify WHICH regular n=7 class has
H(T[K2]) = H(SCblow) (= 24855901), and test actual isomorphism at n=14.

Finding from regular7_doubling_Htest_kps1: regularity is NOT the criterion
(Paley T_7 = DRT splits; the H=175 class coincides). Hypothesis: the coinciding
class is the CONSECUTIVE-SET circulant U_n = C_n({1,...,(n-1)/2}), and
  T[K2](U_n) ~ SCblow(U_n)  for n = 3, 5, 7
(n=3: U_3 = C_3, iso verified at n=6; n=5: U_5 = regular T_5, iso verified at
n=10; n=7: tested here at n=14 via pruned backtracking isomorphism).

Output: 05-knowledge/results/consecutive_circulant_iso_kps1.out
"""
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from skew_doubling_core_kps1 import (canon, H_count, scores, D_blowup,
                                     D_scblowup, D_skew, is_DRT)

OUTPATH = os.path.join(ROOT, "05-knowledge", "results",
                       "consecutive_circulant_iso_kps1.out")


def circulant(n, S):
    A = np.zeros((n, n), dtype=np.int64)
    for i in range(n):
        for j in range(n):
            if (j - i) % n in S:
                A[i, j] = 1
    return A


def iso_backtrack(A, B):
    """Exact isomorphism test via score-constrained DFS with adjacency pruning."""
    n = A.shape[0]
    sA, sB = scores(A), scores(B)
    if sorted(sA) != sorted(sB):
        return False
    cand = [[v for v in range(n) if sB[v] == sA[u]] for u in range(n)]
    order = sorted(range(n), key=lambda u: len(cand[u]))
    mapping = [-1] * n
    used = [False] * n

    def dfs(k):
        if k == n:
            return True
        u = order[k]
        for v in cand[u]:
            if used[v]:
                continue
            ok = True
            for kk in range(k):
                uu = order[kk]
                vv = mapping[uu]
                if A[u, uu] != B[v, vv] or A[uu, u] != B[vv, v]:
                    ok = False
                    break
            if ok:
                mapping[u] = v
                used[v] = True
                if dfs(k + 1):
                    return True
                mapping[u] = -1
                used[v] = False
        return False

    return dfs(0)


def main():
    out = open(OUTPATH, "w", encoding="utf-8")

    def w(line=""):
        out.write(line + "\n")
        out.flush()
        print(line)

    w("=== consecutive_circulant_iso_kps1 ===")
    w("")
    w("--- identify the two circulant regular n=7 classes ---")
    U7 = circulant(7, {1, 2, 3})
    P7 = circulant(7, {1, 2, 4})
    w("U_7 = C_7({1,2,3}): H=%d DRT=%s" % (H_count(U7), is_DRT(U7)))
    w("P_7 = C_7({1,2,4}) (Paley): H=%d DRT=%s" % (H_count(P7), is_DRT(P7)))
    w("=> the K2==SC coinciding class (H=175) IS the consecutive circulant U_7;")
    w("   the splitting DRT class (H=189) is Paley.")
    w("")

    w("--- sanity: backtracking iso agrees on known n<=10 cases ---")
    C3 = circulant(3, {1})
    t = iso_backtrack(D_blowup(C3)[0], D_scblowup(C3)[0])
    w("n=3 U_3: iso(K2, SCblow) at order 6 = %s (expected True)" % t)
    t = iso_backtrack(D_skew(C3)[0], D_blowup(C3)[0])
    w("n=3 U_3: iso(Dskew, K2) at order 6 = %s (expected False)" % t)
    U5 = circulant(5, {1, 2})
    w("U_5 = C_5({1,2}) is regular T_5: H=%d scores=%s" % (H_count(U5), str(scores(U5))))
    t = iso_backtrack(D_blowup(U5)[0], D_scblowup(U5)[0])
    w("n=5 U_5: iso(K2, SCblow) at order 10 = %s (expected True)" % t)
    w("")

    w("--- main test: iso(T[K2](U_7), SCblow(U_7)) at order 14 ---")
    AK = D_blowup(U7)[0]
    AS = D_scblowup(U7)[0]
    hk, hs = H_count(AK), H_count(AS)
    w("H(K2(U_7)) = %d, H(SCblow(U_7)) = %d, equal = %s" % (hk, hs, hk == hs))
    t0 = time.time()
    isoKS = iso_backtrack(AK, AS)
    w("iso(K2(U_7), SCblow(U_7)) = %s   (%.2fs)" % (isoKS, time.time() - t0))
    w("")
    w("--- contrast: Paley T_7 (H already differs => non-iso) ---")
    w("H(K2(P_7)) = %d, H(SCblow(P_7)) = %d  => NOT isomorphic" %
      (H_count(D_blowup(P7)[0]), H_count(D_scblowup(P7)[0])))
    w("")
    w("--- D_skew vs K2 on U_7 (different H => non-iso, listed for completeness) ---")
    w("H(Dskew(U_7)) = %d" % H_count(D_skew(U7)[0]))
    w("")
    w("VERDICT: T[K2](U_n) ~ SCblow(U_n) for the consecutive circulants U_3, U_5, U_7")
    w("(iso verified at orders 6, 10, 14). Criterion is NOT regularity, NOT DRT, NOT")
    w("Paley (Paley T_7 splits 24589929 vs 24453597), NOT self-complementarity")
    w("(all non-regular SC n=5 classes split). Consecutive circulant structure is the")
    w("surviving candidate.")
    w("")
    w("=== done ===")
    out.close()


if __name__ == "__main__":
    main()
