"""SIGNED-LRC collisions via the DST reflection matrix R_sigma = S^{-1} diag(sigma) S.
monad-explorer-2026-06-06-S707.

Phi_eps = S eps,  S[t,i] = sin(2 pi t i / C),  t,i in {1,..,n-1}=H (half-system).
Collision eps~eps' <=> Phi_eps'(t) = sigma_t Phi_eps(t) for a per-freq sign sigma.
=> eps' = R_sigma eps,  R_sigma = S^{-1} diag(sigma) S.

For the SUBGROUP REFLECTION sigma_q: sigma_t = +1 if q|t else -1 (q prime, q|C),
we test:
  (1) Is R_{sigma_q} an exact integer signed-permutation-like matrix?  Print it.
  (2) Equivalent runner-space description of the silent move.
  (3) The exact silent condition on eps and the resulting collision count.
This pins the mechanism so the converse (composite => collision) becomes constructive.
"""
import numpy as np
import math
from itertools import product
from collections import Counter

np.set_printoptions(precision=3, suppress=True, linewidth=200)


def dst_matrix(C):
    n1 = (C - 1) // 2
    H = list(range(1, n1 + 1))
    S = np.array([[math.sin(2 * math.pi * t * i / C) for i in H] for t in H])
    return S, H


def subgroup_reflection_sigma(C, q):
    # frequencies t = 1..(C-1)/2 ; sigma_t = +1 if q|t else -1
    n1 = (C - 1) // 2
    return np.array([1.0 if (t % q == 0) else -1.0 for t in range(1, n1 + 1)])


def analyze_C(C, q):
    S, H = dst_matrix(C)
    Sinv = np.linalg.inv(S)
    sigma = subgroup_reflection_sigma(C, q)
    R = Sinv @ np.diag(sigma) @ S
    Rr = np.round(R)
    err = np.max(np.abs(R - Rr))
    n1 = len(H)
    print(f"\n=== C={C}  q={q}  (subgroup order q, points multiples of C/q) ===")
    print(f"  runners H = {H}")
    print(f"  sigma (freq signs, +1 iff q|t) = {[int(x) for x in sigma]}")
    print(f"  R_sigma rounded? max|R - round(R)| = {err:.3e}")
    if err < 1e-6:
        print("  R_sigma is INTEGER. matrix (rows=output runner index, cols=input):")
        print(Rr.astype(int))
        # is it a signed permutation? each row/col exactly one nonzero +-1?
        is_signed_perm = True
        for v in [Rr, Rr.T]:
            for row in v:
                nz = np.nonzero(row)[0]
                if len(nz) != 1 or abs(row[nz[0]]) != 1:
                    is_signed_perm = False
        print(f"  signed permutation? {is_signed_perm}")
        # If not signed perm, find which eps in {+-1}^n1 map to {+-1}^n1 under R (the silent cuts).
        if n1 <= 16:
            cnt = 0
            for bits in product([1, -1], repeat=n1):
                eps = np.array(bits)
                epsp = Rr @ eps
                if np.all(np.abs(epsp) == 1) and not np.array_equal(epsp, eps) and not np.array_equal(epsp, -eps):
                    cnt += 1
            print(f"  # eps in cube with R eps in cube and != +-eps (gross, counts each pair twice & both signs): {cnt}")
    return R


# smallest composite cases, each prime factor q
for C, qs in [(9, [3]), (15, [3, 5]), (21, [3, 7]), (25, [5]), (27, [3]), (33, [3, 11]), (35, [5, 7])]:
    for q in qs:
        analyze_C(C, q)
