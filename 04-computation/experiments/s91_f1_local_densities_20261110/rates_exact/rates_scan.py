#!/usr/bin/env python3
"""S91 rates_exact: spectral radii of the twisted transfer operators T_k modulo every odd m in a range.

Reduction (PROVED in the report; checked numerically in complement_symmetry.py): T_k commutes with the complement
operator (C u)(eta, beta) = w^(-(eta - beta)) u(-1 - eta, 1 - beta), w = zeta^(k/2), C^2 = zeta^k, so the spectrum of
the 2m x 2m matrix T_k is the union of the spectra of the two m x m matrices

        (T_k^{+-} f)(eta) = f(2 eta) +- w^(-(eta + 1)) f(-2 eta - 2).

k and m - k give complex-conjugate spectra, so only 1 <= k <= (m - 1)/2 is scanned.

Method (all floating point; status VERIFIED, not exact):
  * m <= DENSE_ALL: dense eigenvalues (LAPACK) of every T_k^{+-}; the power-iteration estimate below is also computed
    and compared, to calibrate the filter (the largest error is printed and must be below MARGIN / 3);
  * m > DENSE_ALL: a power-iteration estimate of every spectral radius (geometric mean growth of the norm over WINDOW
    steps after BURN steps, random complex start vector), then dense eigenvalues of every (k, sign) whose estimate is
    within MARGIN of the largest estimate for that m, or is at least GLOBAL (candidates for the global list).
    The per-m maximum and every value >= GLOBAL + MARGIN reported downstream are dense values.

Usage: python3 rates_scan.py M_LO M_HI OUT.npz [DENSE_ALL]
Output arrays: m, k, sign, est (estimate / 2), rate (dense rho / 2, or nan when not computed).
"""
import sys
import time

import numpy as np
from numba import njit

MARGIN = 0.012
GLOBAL = 0.835
BURN = 500
WINDOW = 1000


@njit(cache=False)
def power_estimates(m, kks, burn, window, seed):
    np.random.seed(seed)
    nk = kks.shape[0]
    out = np.zeros((nk, 2))
    ph = np.zeros(m, np.complex128)
    for t in range(m):
        ph[t] = np.exp(-2j * np.pi * t / m)
    idx0 = np.zeros(m, np.int64)
    idx1 = np.zeros(m, np.int64)
    for eta in range(m):
        idx0[eta] = (2 * eta) % m
        idx1[eta] = (-2 * eta - 2) % m
    v = np.zeros(m, np.complex128)
    nv = np.zeros(m, np.complex128)
    for j in range(nk):
        kk = kks[j]
        for si in range(2):
            sgn = 1.0 if si == 0 else -1.0
            for eta in range(m):
                v[eta] = np.random.randn() + 1j * np.random.randn()
            acc = 0.0
            for step in range(burn + window):
                s2 = 0.0
                for eta in range(m):
                    val = v[idx0[eta]] + sgn * ph[(kk * (eta + 1)) % m] * v[idx1[eta]]
                    nv[eta] = val
                    s2 += val.real * val.real + val.imag * val.imag
                nrm = np.sqrt(s2)
                if step >= burn:
                    acc += np.log(nrm)
                for eta in range(m):
                    v[eta] = nv[eta] / nrm
            out[j, si] = np.exp(acc / window)
    return out


def reduced_dense(m, k, sign):
    inv2 = pow(2, -1, m)
    kk = (k * inv2) % m
    eta = np.arange(m)
    M = np.zeros((m, m), dtype=np.complex128)
    np.add.at(M, (eta, (2 * eta) % m), 1.0)
    np.add.at(M, (eta, (-2 * eta - 2) % m), sign * np.exp(-2j * np.pi * ((kk * (eta + 1)) % m) / m))
    return M


def rho_dense(m, k, sign):
    return float(np.abs(np.linalg.eigvals(reduced_dense(m, k, sign))).max())


def main():
    lo, hi, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    dense_all = int(sys.argv[4]) if len(sys.argv) > 4 else 201
    rows = []
    t0 = time.time()
    ok = True
    worst_filter = 0.0
    ndense = 0
    for m in range(lo | 1, hi + 1, 2):
        ks = np.arange(1, (m - 1) // 2 + 1)
        inv2 = pow(2, -1, m)
        kks = (ks * inv2) % m
        est = power_estimates(m, kks.astype(np.int64), BURN, WINDOW, 1000003 * m + 7) / 2.0
        rate = np.full(est.shape, np.nan)
        if m <= dense_all:
            for j, k in enumerate(ks):
                for si, sign in enumerate((+1, -1)):
                    rate[j, si] = rho_dense(m, int(k), sign) / 2.0
                    ndense += 1
            err = float(np.abs(est - rate).max())
            worst_filter = max(worst_filter, err)
            # the filter must be accurate enough for the margin used above DENSE_ALL
            if err > MARGIN / 3:
                ok = False
                print(f"FILTER ERROR {err:.5f} at m = {m}", flush=True)
        else:
            top = est.max()
            for j, k in enumerate(ks):
                for si, sign in enumerate((+1, -1)):
                    if est[j, si] >= top - MARGIN or est[j, si] >= GLOBAL:
                        rate[j, si] = rho_dense(m, int(k), sign) / 2.0
                        ndense += 1
                        worst_filter = max(worst_filter, abs(rate[j, si] - est[j, si]))
                        if abs(rate[j, si] - est[j, si]) > MARGIN / 3:
                            ok = False
                            print(f"FILTER ERROR {abs(rate[j, si] - est[j, si]):.5f} at m = {m} k = {k} sign = {sign}", flush=True)
        best = np.nanmax(rate)
        if not (best < 1.0 - 1e-9):
            ok = False
            print(f"rate >= 1 at m = {m}", flush=True)
        for j, k in enumerate(ks):
            for si, sign in enumerate((+1, -1)):
                rows.append((m, int(k), sign, est[j, si], rate[j, si]))
        if (m // 2) % 10 == 0 or m == hi:
            jj, ss = np.unravel_index(np.nanargmax(rate), rate.shape)
            print(f"m = {m:4d}: max_k rho/2 = {best:.6f} at k = {ks[jj]} sign {'+' if ss == 0 else '-'}   [{time.time() - t0:.0f} s]", flush=True)
    arr = np.array(rows)
    np.savez(out, m=arr[:, 0].astype(np.int64), k=arr[:, 1].astype(np.int64), sign=arr[:, 2].astype(np.int64),
             est=arr[:, 3], rate=arr[:, 4])
    print(f"range [{lo}, {hi}]: {len(rows)} operators; dense eigenvalue computations: {ndense}; "
          f"largest |power estimate - dense| = {worst_filter:.5f} (margin {MARGIN}); time {time.time() - t0:.0f} s")
    print("ALL CHECKS PASSED" if ok else "CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
