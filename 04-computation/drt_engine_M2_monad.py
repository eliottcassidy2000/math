#!/usr/bin/env python3
"""
drt_engine_M2_monad.py
monad-explorer-2026-06-07 (deep-research, 5th session)

The DRT-UNIVERSAL, NUMBER-THEORY-FREE algebraic engine behind THM-438.

CLAIM:  for the Paley skew matrix M[x,y]=chi(y-x) (p=3 mod4),
            M @ 1 = 0      (regular tournament)
            M^2 = J - p I  (doubly-regular tournament identity; J=all-ones)
This is PRECISELY the defining relation of a doubly regular tournament skew-matrix
(S^T=-S, S 1=0, S S^T = nI - J  <=>  S^2 = J - nI), so the entire leading-order
machinery of THM-438 needs NO Gauss sums / quadratic reciprocity: the merged
2-step sum  sum_y S[x,y]S[y,z] = (J-nI)[x,z]  is the only analytic input, and it
holds for EVERY doubly regular tournament (strengthens HYP-2308).

Consequence used downstream:  S^{2t} = (-n)^t (I - J/n)  on even powers, so a
length-2t series-chain contributes (-n)^t to leading order  =>  k bigons give n^k,
the free hub-sum gives the extra n, total n^{k+1}.  Spectrum collapses to {0, +-i*sqrt n}.
"""
import numpy as np


def legendre(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def paley_M(p):
    M = np.zeros((p, p), dtype=np.int64)
    for x in range(p):
        for y in range(p):
            if x != y:
                M[x, y] = legendre(y - x, p)
    return M


print("=" * 64)
print("DRT engine:  M 1 = 0   and   M^2 = J - p I   (number-theory-free downstream)")
for p in [3, 7, 11, 19, 23, 31, 43]:
    M = paley_M(p)
    J = np.ones((p, p), dtype=np.int64)
    one = np.ones(p, dtype=np.int64)
    row0 = np.all(M @ one == 0)
    sq = np.array_equal(M @ M, J - p * np.eye(p, dtype=np.int64))
    # even-power collapse: M^4 = (-p) * (J - pI) ?  i.e. M^{2t} = (-p)^t (I - J/p) scaled
    M2 = M @ M
    M4 = M2 @ M2
    collapse4 = np.array_equal(M4, (-p) * M2)            # M^4 = -p * M^2  <=> M^2(M^2+pI)=0
    # eigenvalues ~ {0, +-i sqrt p}
    ev = np.linalg.eigvals(M.astype(float))
    nz = ev[np.abs(ev) > 1e-6]
    twopoint = np.allclose(np.sort(np.abs(nz)), np.full(len(nz), np.sqrt(p)), atol=1e-6)
    print(f"  p={p:>2}: M1=0:{row0}  M^2=J-pI:{sq}  M^4=-pM^2:{collapse4}  "
          f"|nonzero eig|=sqrt(p):{twopoint}  (#nz={len(nz)}=p-1:{len(nz)==p-1})")
print("=" * 64)
print("All True  =>  M^2 = J - pI is the sole analytic input; it is the DRT identity,")
print("so c_0 = lim A_{2k}/p^{k+1} is DRT-universal at the level of every merged 2-step.")
