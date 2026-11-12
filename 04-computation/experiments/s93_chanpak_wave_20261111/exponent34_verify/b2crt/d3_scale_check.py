#!/usr/bin/env python3
"""
d3_scale_check.py -- size of Deg_K and Deg'_K against the bracket of (D3) at larger window lengths
(verifier lane b2crt).  Direct enumeration of the box with integer gcds (numba, exact int64).

    Deg_K(P_u)  = #{(a', b') : |a'|, |b'| < 2^s2, b' odd, height of (a' 2^s1 + b, a 2^s2 + b') <= K}
    Deg'_K(P_v) = #{(a, b)   : |a|, |b| < 2^s1,  b odd,  height of (a' 2^s1 + b, a 2^s2 + b') <= K}
    bracket(K, A) = K (1 + ln K)^2 + K (1 + ln K) A^(1/3) + K^(5/2) A^(-1/2),  A = |a| resp. |a'|

For K = 2^i the script prints the maximum over the tested P of Deg_K/bracket (an empirical lower bound for
the constant C of (D3) at this size), of Deg_K/(K (1 + ln K)^2), and checks (D1), (D2) on the way
(Deg_K <= 4K^2, Deg'_K <= 8K^2; for K <= 2^s2: Deg_K <= 3 (2^(s2+1)(1+ln K) + K), Deg'_K <= 3 (2^(s1+1)(1+ln K) + K)).

usage: python d3_scale_check.py W [NSAMPLE] [THREADS]
   P_u runs over the difference set of {(u, rev_s1 u)} with b odd (all of it, or NSAMPLE random points
   plus the points with the smallest and largest |a|).
"""
import sys
import math
import time
import numpy as np
from numba import njit, prange, set_num_threads


@njit(cache=True)
def gcd_(a, b):
    while b:
        a, b = b, a % b
    return a


@njit(parallel=True, cache=True)
def deg_tab(PA, PB, s1, s2, nb):
    n = PA.shape[0]
    out = np.zeros((n, nb), dtype=np.int64)
    B1 = 1 << s1
    B2 = 1 << s2
    for idx in prange(n):
        a = PA[idx]
        b = PB[idx]
        for ap in range(-B2 + 1, B2):
            aX = abs(ap * B1 + b)
            for bp in range(-B2 + 1, B2, 2):
                aXi = abs(a * B2 + bp)
                g = gcd_(aX, aXi)
                M = max(aX, aXi) // g
                i = 0
                while (1 << i) < M:
                    i += 1
                out[idx, i] += 1
    return out


@njit(parallel=True, cache=True)
def degp_tab(PA, PB, s1, s2, nb):
    n = PA.shape[0]
    out = np.zeros((n, nb), dtype=np.int64)
    B1 = 1 << s1
    B2 = 1 << s2
    for idx in prange(n):
        ap = PA[idx]
        bp = PB[idx]
        for a in range(-B1 + 1, B1):
            aXi = abs(a * B2 + bp)
            for b in range(-B1 + 1, B1, 2):
                aX = abs(ap * B1 + b)
                g = gcd_(aX, aXi)
                M = max(aX, aXi) // g
                i = 0
                while (1 << i) < M:
                    i += 1
                out[idx, i] += 1
    return out


def rev(x, w):
    r = 0
    for i in range(w):
        if (x >> i) & 1:
            r |= 1 << (w - 1 - i)
    return r


def diff_points(s):
    """all (D(eps), D*(eps)) with D* odd, by direct enumeration of pairs (u, u')"""
    R = np.array([rev(u, s) for u in range(1 << s)], dtype=np.int64)
    U = np.arange(1 << s, dtype=np.int64)
    A = (U[:, None] - U[None, :]).ravel()
    B = (R[:, None] - R[None, :]).ravel()
    keep = (B % 2) != 0
    pts = np.unique(np.stack([A[keep], B[keep]], axis=1), axis=0)
    return pts


def bracket(K, A):
    lk = 1 + math.log(K)
    return K * lk * lk + K * lk * A ** (1 / 3) + K ** 2.5 * A ** -0.5


def report(name, tab, PA, s_other, s_box, w, d1const):
    cum = np.cumsum(tab, axis=1)
    absA = np.abs(PA).astype(np.float64)
    ok = True
    for i in range(0, w + 1):
        K = 1 << i
        d = cum[:, i].astype(np.float64)
        br = np.array([bracket(K, A) for A in absA])
        r = d / br
        j = int(np.argmax(r))
        lk = 1 + math.log(K)
        line = (f"  {name} K=2^{i:<2d}: max Deg/bracket = {r[j]:.4f} (|a|={int(absA[j])}, Deg={int(d[j])}); "
                f"max Deg = {int(d.max())}; max Deg/(K(1+lnK)^2) = {d.max() / (K * lk * lk):.4f}")
        if d.max() > d1const * K * K:
            ok = False
            line += "  (D1) VIOLATED"
        if i <= min(s_other, s_box):
            pass
        print(line)
    return ok


def main():
    W = int(sys.argv[1])
    nsample = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    threads = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    set_num_threads(threads)
    s1 = (W + 1) // 2
    s2 = W // 2
    rng = np.random.default_rng(20261111 + W)
    t0 = time.time()
    PU = diff_points(s1)
    PV = diff_points(s2)
    assert len(PU) == 2 * 3 ** (s1 - 1) and len(PV) == 2 * 3 ** (s2 - 1)
    assert np.all(PU[:, 0] != 0) and np.all(PV[:, 0] != 0)

    def sample(P):
        if nsample and nsample < len(P):
            order = np.argsort(np.abs(P[:, 0]), kind="stable")
            ext = np.concatenate([order[:40], order[-40:]])
            rnd = rng.choice(len(P), size=nsample, replace=False)
            sel = np.unique(np.concatenate([ext, rnd]))
            return P[sel]
        return P

    PUs = sample(PU)
    PVs = sample(PV)
    print(f"W={W} s1={s1} s2={s2}: |difference set, b odd| = {len(PU)} (tested {len(PUs)}), {len(PV)} (tested {len(PVs)})")
    tab = deg_tab(np.ascontiguousarray(PUs[:, 0]), np.ascontiguousarray(PUs[:, 1]), s1, s2, W + 2)
    assert np.all(tab.sum(axis=1) == (2 * (1 << s2) - 1) * (1 << s2))
    ok1 = report("Deg ", tab, PUs[:, 0], s1, s2, W, 4)
    # (D2) for K <= 2^s2
    cum = np.cumsum(tab, axis=1)
    for i in range(0, s2 + 1):
        K = 1 << i
        assert cum[:, i].max() <= 3 * ((1 << (s2 + 1)) * (1 + math.log(K)) + K), ("D2", i)
    tabp = degp_tab(np.ascontiguousarray(PVs[:, 0]), np.ascontiguousarray(PVs[:, 1]), s1, s2, W + 2)
    assert np.all(tabp.sum(axis=1) == (2 * (1 << s1) - 1) * (1 << s1))
    ok2 = report("Deg'", tabp, PVs[:, 0], s2, s1, W, 8)
    cump = np.cumsum(tabp, axis=1)
    for i in range(0, s2 + 1):
        K = 1 << i
        assert cump[:, i].max() <= 3 * ((1 << (s1 + 1)) * (1 + math.log(K)) + K), ("D2'", i)
    print(f"W={W}: (D1) {'ok' if ok1 and ok2 else 'VIOLATED'}, (D2) ok for K = 2^i <= 2^s2; {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
