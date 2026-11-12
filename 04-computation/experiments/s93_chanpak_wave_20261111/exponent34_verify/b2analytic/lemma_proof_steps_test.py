#!/usr/bin/env python3
"""
Numerical test of the two CITED inequalities, in the exact forms used in the proof of the Discrepancy lemma of
THM-5141 (B2), on the sequences u_n = a/(2^ds k'), k' odd in (L, 2L]   (verifier: b2analytic).

(V) van der Corput estimate, Kuipers-Niederreiter ch. 1 Thm 2.7:  for integers A < Bq and f twice differentiable
    on [A, Bq] with f'' >= rho > 0 or f'' <= -rho < 0:
        |sum_{n=A}^{Bq} e(f(n))| <= (|f'(Bq) - f'(A)| + 2) (4/sqrt(rho) + 3).
    Here f(l) = F_h(l) = h a/(2^ds (2l+1)), rho = min |F_h''| = 8 h |a| / (2^ds kmax^3),
    |F_h'(lmax) - F_h'(lmin)| = (2 h |a| / 2^ds) (1/kmin^2 - 1/kmax^2).
    Tested: (V1) the cited inequality with these exact rho and derivative difference;
            (V2) the draft's consequence  |S_h| <= 56 L lam^(1/2) + 8 lam^(-1/2) + 6,  lam = h |a|/(2^ds L^3).

(E) Erdos-Turan, Kuipers-Niederreiter ch. 2 Thm 2.5 / eq. (2.40) (discrepancy over arcs mod 1):
        N D_N <= 6 N/(m+1) + (4/pi) sum_{h=1}^{m} (1/h - 1/(m+1)) |S_h|,     every integer m >= 1.
    Tested for m = 1, 2, 4, ..., and m = H = ceil(L |a|^(-1/3)) (the draft's choice) when H <= HMAX.

(C) the chain of the proof with explicit constants: for L >= |a|^(1/3)
        N D_N <= 6 |a|^(1/3) + (4/pi) [112 sqrt2 |a|^(1/3) + 8 zeta(3/2) sqrt2 L^(3/2)|a|^(-1/2) + 6 (1 + ln 2 + ln L)]
    (so the lemma holds with C = 208); the ratio LHS/RHS is reported.

All phases are computed with exact integer reduction ((h mod m)(a mod m) mod m)/m.
usage: python lemma_proof_steps_test.py [n_inst] [seed]
"""
import math
import sys
import time

import numpy as np
from numba import njit, prange

HMAX = 128
ZETA32 = 2.612375348685488


@njit(cache=False)
def one_instance(a, ds, L, res):
    """res: [max ratio V1, max ratio V2, max ratio E, ratio C, N, H, hV1, hV2, mE, ratioE_at_H]"""
    lo = int(math.floor(L)) + 1
    hi = int(math.floor(2.0 * L))
    if lo % 2 == 0:
        lo += 1
    if hi % 2 == 0:
        hi -= 1
    if hi < lo:
        res[4] = 0
        return
    n = (hi - lo) // 2 + 1
    A = abs(float(a))
    pts = np.empty(n, dtype=np.float64)
    ks = np.empty(n, dtype=np.int64)
    for i in range(n):
        k = lo + 2 * i
        m = k << ds
        r = a % m
        if r < 0:
            r += m
        pts[i] = r / m
        ks[i] = k
    p = np.sort(pts)
    # N D_N over arcs mod 1
    best = 0.0
    gmin = 1e300
    gmax = -1e300
    for t in range(2 * n):
        if t < n:
            g = t - n * p[t]
        else:
            g = t - n * (p[t - n] + 1.0)
        if g < gmin:
            gmin = g
        v = g - gmin + 1.0
        if v > best:
            best = v
        if gmax > -1e299:
            v = gmax - g + 1.0
            if v > best:
                best = v
        if g > gmax:
            gmax = g
    ND = best
    H = int(math.ceil(L / A ** (1.0 / 3.0)))
    if H < 1:
        H = 1
    hmax = HMAX
    absS = np.empty(hmax + 1, dtype=np.float64)
    twopi = 2.0 * math.pi
    kmin = float(lo)
    kmax = float(hi)
    rV1 = 0.0
    rV2 = 0.0
    hV1 = 0
    hV2 = 0
    for h in range(1, hmax + 1):
        sr = 0.0
        si = 0.0
        for i in range(n):
            m = ks[i] << ds
            r = ((h % m) * (a % m)) % m
            if r < 0:
                r += m
            ang = twopi * r / m
            sr += math.cos(ang)
            si += math.sin(ang)
        s = math.sqrt(sr * sr + si * si)
        absS[h] = s
        lam = h * A / ((1 << ds) * L ** 3)
        b2 = 56.0 * L * math.sqrt(lam) + 8.0 / math.sqrt(lam) + 6.0
        q = s / b2
        if q > rV2:
            rV2 = q
            hV2 = h
        if n >= 2:
            rho = 8.0 * h * A / ((1 << ds) * kmax ** 3)
            dF = 2.0 * h * A / (1 << ds) * (1.0 / kmin ** 2 - 1.0 / kmax ** 2)
            b1 = (dF + 2.0) * (4.0 / math.sqrt(rho) + 3.0)
            q = s / b1
            if q > rV1:
                rV1 = q
                hV1 = h
    # Erdos-Turan for m = 1, 2, 4, ..., hmax and m = H
    rE = 0.0
    mE = 0
    rEH = -1.0
    m = 1
    while m <= hmax:
        acc = 0.0
        for h in range(1, m + 1):
            acc += (1.0 / h - 1.0 / (m + 1)) * absS[h]
        rhs = 6.0 * n / (m + 1) + 4.0 / math.pi * acc
        q = ND / rhs
        if q > rE:
            rE = q
            mE = m
        m *= 2
    if H <= hmax:
        acc = 0.0
        for h in range(1, H + 1):
            acc += (1.0 / h - 1.0 / (H + 1)) * absS[h]
        rhs = 6.0 * n / (H + 1) + 4.0 / math.pi * acc
        rEH = ND / rhs
        if rEH > rE:
            rE = rEH
            mE = H
    rC = -1.0
    if L ** 3 >= A:
        rhsC = 6.0 * A ** (1.0 / 3.0) + 4.0 / math.pi * (
            112.0 * math.sqrt(2.0) * A ** (1.0 / 3.0) + 8.0 * ZETA32 * math.sqrt(2.0) * L ** 1.5 / math.sqrt(A)
            + 6.0 * (1.0 + math.log(2.0) + math.log(L)))
        rC = ND / rhsC
    res[0] = rV1
    res[1] = rV2
    res[2] = rE
    res[3] = rC
    res[4] = n
    res[5] = H
    res[6] = hV1
    res[7] = hV2
    res[8] = mE
    res[9] = rEH


@njit(parallel=True, cache=False)
def run(a_arr, ds_arr, L_arr, out):
    for t in prange(a_arr.shape[0]):
        one_instance(a_arr[t], ds_arr[t], L_arr[t], out[t])


def main():
    n_inst = int(sys.argv[1]) if len(sys.argv) > 1 else 12000
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 5141
    rng = np.random.default_rng(seed)
    t0 = time.time()
    blocks = []
    # (i) generic: |a| log-uniform <= 2^40, L log-uniform <= 2^11
    la = rng.uniform(0, 40, n_inst)
    a1 = np.maximum(1, np.floor(2.0 ** la)).astype(np.int64) * (rng.integers(0, 2, n_inst) * 2 - 1)
    L1 = 2.0 ** rng.uniform(0.6, 11, n_inst)
    blocks.append(("generic |a| <= 2^40, L <= 2^11", a1, rng.integers(0, 2, n_inst).astype(np.int64), L1))
    # (ii) the regime of the lemma: |a|^(1/3) <= L <= |a|
    L2 = 2.0 ** rng.uniform(1, 11, n_inst)
    gam = rng.uniform(1.0, 3.0, n_inst)
    a2 = np.maximum(1, np.floor(L2 ** gam)).astype(np.int64)
    blocks.append(("regime |a|^(1/3) <= L <= |a|", a2, rng.integers(0, 2, n_inst).astype(np.int64), L2))
    # (iii) small a, large L (lam tiny) and large a, small L (lam >= 1: the 'trivial bound' regime)
    a3 = rng.integers(1, 50, n_inst).astype(np.int64)
    L3 = 2.0 ** rng.uniform(1, 11, n_inst)
    blocks.append(("small |a| < 50", a3, rng.integers(0, 2, n_inst).astype(np.int64), L3))
    a4 = (2.0 ** rng.uniform(20, 46, n_inst)).astype(np.int64)
    L4 = 2.0 ** rng.uniform(0.6, 6, n_inst)
    blocks.append(("large |a| >= 2^20, L <= 64 (lam >= 1 mostly)", a4, rng.integers(0, 2, n_inst).astype(np.int64), L4))
    for (name, a, ds, L) in blocks:
        out = np.zeros((len(a), 10))
        run(a, ds, L, out)
        ok = out[:, 4] >= 1
        print(f"[{name}] instances = {int(ok.sum())} (h <= {HMAX})   ({time.time()-t0:.1f}s)")
        for col, label in ((0, "(V1) |S_h| / K-N Thm 2.7 bound (exact rho, exact f' difference)"),
                           (1, "(V2) |S_h| / (56 L lam^(1/2) + 8 lam^(-1/2) + 6)"),
                           (2, "(E)  N D_N / Erdos-Turan bound (m = 1,2,4,...,128 and m = H)"),
                           (3, "(C)  N D_N / explicit chain bound (L >= |a|^(1/3) only)")):
            i = int(np.argmax(np.where(ok, out[:, col], -1)))
            extra = ""
            if col == 0:
                extra = f", h = {int(out[i,6])}"
            if col == 1:
                extra = f", h = {int(out[i,7])}"
            if col == 2:
                extra = f", m = {int(out[i,8])}, H = {int(out[i,5])}"
            print(f"   max {label}: {out[i,col]:.6f}   at a = {a[i]}, ds = {ds[i]}, L = {L[i]:.4f}, N = {int(out[i,4])}{extra}")
        sel = ok & (out[:, 9] >= 0)
        if sel.any():
            print(f"   (E at the draft's m = H, {int(sel.sum())} instances with H <= {HMAX}): max ratio {out[sel, 9].max():.6f}")
        sys.stdout.flush()
    print(f"total time {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
