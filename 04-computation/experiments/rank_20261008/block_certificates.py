#!/usr/bin/env python3
"""THM-4611 statement 3: re-verify the five block certificates from scratch.  For each map: check the branch condition,
independence and rank, the coupling group, the minimal one-step coupling-covariance rank (<= 2 means no one-step form),
property (i) (only the identity coupling has all roots zero), then compute every k-step block matrix exactly, count the
zero blocks, verify the stated integer form exactly (block_lmi2.exact_check_vec), and report that form's own margin
min (tr - 2 lambda_max)/tr and alpha_max = min tr/lambda_max - 2 (float, informational).  Usage: python3 block_certificates.py"""
import math, time
import numpy as np
from block_balance import int_coords, ratio_group, coupling_D, std_r
from block_lmi2 import last_level_unique, exact_check_vec

CERTS = [
    (5, [1, 2, 3, 7, 11], 2, [[4, -1, -1, -1], [-1, 4, -1, -1], [-1, -1, 4, -1], [-1, -1, -1, 4]]),
    (5, [1, 2, 3, 7, 1], 4, [[20, -4, -8], [-4, 15, -3], [-8, -3, 20]]),
    (5, [1, 1, 2, 3, 7], 5, [[8, -2, -3], [-2, 6, -1], [-3, -1, 8]]),
    (7, [1, 1, 1, 1, 2, 3, 5], 4, [[6, -1, -1], [-1, 5, -1], [-1, -1, 6]]),
    (7, [1, 1, 1, 2, 3, 5, 11], 2, [[5, -1, -2, -1], [-1, 5, 0, -1], [-2, 0, 6, -1], [-1, -1, -1, 5]]),
]

def form_margin(Qi, U):
    Q = np.array(Qi, dtype=float); L = np.linalg.cholesky(Q); Li = np.linalg.inv(L)
    worst = 1.0; amax = np.inf
    for s in range(0, len(U), 400000):
        W = Li @ U[s:s + 400000].astype(float) @ Li.T
        ev = np.linalg.eigvalsh(W); tr = ev.sum(axis=1)
        worst = min(worst, float(((tr - 2 * ev[:, -1]) / tr).min())); amax = min(amax, float((tr / ev[:, -1]).min() - 2))
    return worst, amax

if __name__ == '__main__':
    for d, m, k, Qi in CERTS:
        t0 = time.time()
        r = std_r(m)
        assert all((m[i] * i + r[i]) % d == 0 for i in range(d)) and all(math.gcd(x, d) == 1 for x in m)
        rho, v = int_coords(m)
        Lam = sum(math.log(x / d) for x in m) / d
        G = ratio_group(m, d)
        ranks = {}; zero_noid = []
        for a in G:
            for b in range(d):
                D = np.array(coupling_D(v, a, b, d, rho)); rk = int(np.linalg.matrix_rank(D)) if D.any() else 0
                ranks[(a, b)] = rk
                if rk == 0 and (a, b) != (1, 0): zero_noid.append((a, b))
        min1 = min(rk for key, rk in ranks.items() if key != (1, 0))
        U, nstates, frozen = last_level_unique(d, m, r, v, rho, k)
        ok, unsure = exact_check_vec(Qi, U)
        mg, amax = form_margin(Qi, U)
        print(f"Z_{d} m = {m} r = {r}: rank {rho}, Lambda = {Lam:+.4f}, coupling group {G}, min one-step coupling rank {min1}"
              f" ({'no one-step form' if min1 <= 2 else 'one-step form possible'}), zero-root non-identity couplings {zero_noid};"
              f" k = {k}: {nstates} states, {frozen} zero blocks, {len(U)} distinct nonzero blocks; Q = {Qi}: exact {ok}"
              f" (integer fallbacks {unsure}), margin {mg:+.4f}, alpha_max {amax:.4f}  [{time.time() - t0:.0f}s]", flush=True)
