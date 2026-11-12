#!/usr/bin/env python3
"""
Independent test of the "Discrepancy lemma" in the proof of THM-5141 (B2)  (verifier: b2analytic).

Lemma (draft).  a != 0 integer, ds in {0,1}, L >= 1, theta real, dl > 0:
    #{odd k' in (L, 2L] : ||a/(2^ds k') + theta|| <= dl}
        <= 2 dl L + C (|a|^(1/3) + L^(3/2) |a|^(-1/2) + 1 + ln L).

For one instance (a, ds, L) the supremum over ALL (theta, dl) of  count - 2 dl L  is computed exactly from the sorted
points p_n = frac(a/(2^ds k')):  a closed arc with first point p_i and last point p_j (cyclic order) has minimal
length p_j - p_i and contains j - i + 1 points, so
    E_L(a, ds, L) = max over i <= j of (j - i + 1) - L (p_j - p_i)            (doubled array p_(i+N) = p_i + 1).
The same computation with N (the number of points) in place of L gives N D_N^+ ; the open arcs give N D_N^- ;
N D_N = max of the two is the discrepancy (over arcs mod 1) that the Erdos-Turan inequality bounds.

Bound function:  B(a, L) = |a|^(1/3) + L^(3/2) |a|^(-1/2) + 1 + ln L.
We report  max E_L / B  (the constant the lemma needs) and  max N D_N / B  (the constant the proof's inner
statement "N D_N <= C B for L >= |a|^(1/3)" needs).

Points are computed with exact integer reduction (a mod m)/m, m = 2^ds k' <= 2^18; two distinct such fractions
differ by >= 2^-36, equal fractions give equal doubles, so the float sort reproduces the exact order.
The worst instances are re-verified with fractions.Fraction (exact rational arithmetic) at the end.

usage: python lemma_discrepancy_test.py [n_random] [seed]
"""
import math
import sys
import time
from fractions import Fraction

import numpy as np
from numba import njit, prange


@njit(cache=False)
def fill_points(a, ds, L, buf):
    lo = int(math.floor(L)) + 1
    hi = int(math.floor(2.0 * L))
    if lo % 2 == 0:
        lo += 1
    n = 0
    k = lo
    while k <= hi:
        m = k << ds
        r = a % m
        if r < 0:
            r += m
        buf[n] = r / m
        n += 1
        k += 2
    return n


@njit(cache=False)
def excess(p, n, c):
    """max over i<=j (cyclic, doubled) of (j-i+1) - c (p_j - p_i); p sorted, n points, c >= n."""
    best = 0.0
    gmin = 1e300
    for t in range(2 * n):
        if t < n:
            g = t - c * p[t]
        else:
            g = t - c * (p[t - n] + 1.0)
        if g < gmin:
            gmin = g
        v = g - gmin + 1.0
        if v > best:
            best = v
    return best


@njit(cache=False)
def deficiency(p, n):
    """N D_N^- : max over i<j<=i+n (cyclic) of n (p_j - p_i) - (j - i - 1)."""
    best = 0.0
    gmax = -1e300
    for t in range(2 * n):
        if t < n:
            g = t - n * p[t]
        else:
            g = t - n * (p[t - n] + 1.0)
        # value for the pair (i, t), i < t:  g_i - g_t + 1
        if gmax > -1e299:
            v = gmax - g + 1.0
            if v > best:
                best = v
        if g > gmax:
            gmax = g
    return best


@njit(cache=False)
def bound_B(a, L):
    A = abs(float(a))
    return A ** (1.0 / 3.0) + L ** 1.5 / math.sqrt(A) + 1.0 + math.log(L)


@njit(parallel=True, cache=False)
def run_batch(a_arr, ds_arr, L_arr, out):
    n_inst = a_arr.shape[0]
    for t in prange(n_inst):
        a = a_arr[t]
        ds = ds_arr[t]
        L = L_arr[t]
        buf = np.empty(int(L) + 3, dtype=np.float64)
        n = fill_points(a, ds, L, buf)
        if n == 0:
            out[t, 0] = 0.0
            out[t, 1] = 0.0
            out[t, 2] = 0.0
            out[t, 3] = bound_B(a, L)
            continue
        p = np.sort(buf[:n])
        eL = excess(p, n, L)
        dplus = excess(p, n, float(n))
        dminus = deficiency(p, n)
        out[t, 0] = eL
        out[t, 1] = max(dplus, dminus)
        out[t, 2] = n
        out[t, 3] = bound_B(a, L)


@njit(parallel=True, cache=False)
def run_exhaustive(amax, jmax, res):
    """all a in [1, amax], L = j/2 for j in [2, jmax], ds in {0,1}.  res[a] = (max ratio E_L/B, argmax j, ds,
    max ratio ND/B restricted to L >= a^(1/3), argmax j, ds)."""
    for a in prange(1, amax + 1):
        buf = np.empty(jmax + 3, dtype=np.float64)
        best1 = 0.0
        best2 = 0.0
        arg1 = 0
        arg2 = 0
        for j in range(2, jmax + 1):
            L = 0.5 * j
            for ds in range(2):
                n = fill_points(a, ds, L, buf)
                if n == 0:
                    continue
                p = np.sort(buf[:n])
                B = bound_B(a, L)
                r1 = excess(p, n, L) / B
                if r1 > best1:
                    best1 = r1
                    arg1 = 2 * j + ds
                if L * L * L >= a:
                    dplus = excess(p, n, float(n))
                    dminus = deficiency(p, n)
                    r2 = max(dplus, dminus) / B
                    if r2 > best2:
                        best2 = r2
                        arg2 = 2 * j + ds
        res[a, 0] = best1
        res[a, 1] = arg1
        res[a, 2] = best2
        res[a, 3] = arg2


def exact_instance(a, ds, L):
    """Exact (Fraction) recomputation of E_L and N D_N for one instance; L given as a Fraction."""
    L = Fraction(L)
    lo = math.floor(L) + 1
    hi = math.floor(2 * L)
    pts = []
    for k in range(lo, hi + 1):
        if k % 2 == 1:
            m = k << ds
            pts.append(Fraction(a % m, m))
    pts.sort()
    n = len(pts)
    if n == 0:
        return Fraction(0), Fraction(0), 0
    ext = pts + [x + 1 for x in pts]
    bestL = Fraction(0)
    bestN = Fraction(0)
    gminL = None
    gminN = None
    gmaxN = None
    bestD = Fraction(0)
    for t in range(2 * n):
        gL = t - L * ext[t]
        gN = t - n * ext[t]
        if gminL is None or gL < gminL:
            gminL = gL
        if gminN is None or gN < gminN:
            gminN = gN
        bestL = max(bestL, gL - gminL + 1)
        bestN = max(bestN, gN - gminN + 1)
        if gmaxN is not None:
            bestD = max(bestD, gmaxN - gN + 1)
        if gmaxN is None or gN > gmaxN:
            gmaxN = gN
    return bestL, max(bestN, bestD), n


def report(tag, a_arr, ds_arr, L_arr, out, top=6):
    B = out[:, 3]
    r1 = out[:, 0] / B
    mask = (L_arr ** 3 >= np.abs(a_arr.astype(np.float64)))
    r2 = np.where(mask, out[:, 1] / B, 0.0)
    print(f"[{tag}] instances = {len(a_arr)}")
    i1 = np.argsort(-r1)[:top]
    print(f"  max (count - 2 dl L)/B over all theta, dl:  {r1.max():.6f}")
    for i in i1:
        print(f"     a = {a_arr[i]}, ds = {ds_arr[i]}, L = {L_arr[i]!r}: N = {int(out[i,2])}, E_L = {out[i,0]:.6f}, "
              f"B = {B[i]:.6f}, ratio = {r1[i]:.6f}")
    i2 = np.argsort(-r2)[:top]
    print(f"  max N D_N / B over instances with L >= |a|^(1/3):  {r2.max():.6f}")
    for i in i2:
        print(f"     a = {a_arr[i]}, ds = {ds_arr[i]}, L = {L_arr[i]!r}: N = {int(out[i,2])}, N D_N = {out[i,1]:.6f}, "
              f"B = {B[i]:.6f}, ratio = {r2[i]:.6f}")
    # the lemma with the proof's explicit constant 208 (see the report): must never fail
    viol = np.sum(out[:, 0] > 208.0 * B)
    viol1 = np.sum(out[:, 0] > 1.0 * B + 1e-9)
    print(f"  violations of the lemma with C = 208: {viol};  with C = 1: {viol1}")
    return i1, i2


def main():
    n_random = int(sys.argv[1]) if len(sys.argv) > 1 else 200000
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 20261111
    rng = np.random.default_rng(seed)
    t0 = time.time()

    # ---------------- 1. exhaustive small range
    amax, jmax = 4096, 1024
    res = np.zeros((amax + 1, 4))
    run_exhaustive(amax, jmax, res)
    ia = int(np.argmax(res[:, 0]))
    ib = int(np.argmax(res[:, 2]))
    print(f"[exhaustive] a in [1, {amax}], L = j/2 (2 <= j <= {jmax}), ds in {{0,1}}: "
          f"{amax * (jmax - 1) * 2} instances  ({time.time()-t0:.1f}s)")
    c1 = int(res[ia, 1])
    print(f"  max (count - 2 dl L)/B = {res[ia,0]:.6f} at a = {ia}, L = {(c1//2)/2}, ds = {c1 % 2}")
    c2 = int(res[ib, 3])
    print(f"  max N D_N / B (L >= a^(1/3)) = {res[ib,2]:.6f} at a = {ib}, L = {(c2//2)/2}, ds = {c2 % 2}")
    eL, nd, n = exact_instance(ia, c1 % 2, Fraction(c1 // 2, 2))
    print(f"  exact recheck of the first: N = {n}, E_L = {eL} = {float(eL):.6f}")
    eL, nd, n = exact_instance(ib, c2 % 2, Fraction(c2 // 2, 2))
    print(f"  exact recheck of the second: N = {n}, N D_N = {nd} = {float(nd):.6f}")
    # top few a for each
    order = np.argsort(-res[:, 0])[:5]
    print("  top a for E_L/B:", [(int(a), round(float(res[a, 0]), 4), (int(res[a, 1]) // 2) / 2, int(res[a, 1]) % 2)
                                   for a in order])
    order = np.argsort(-res[:, 2])[:5]
    print("  top a for N D_N/B:", [(int(a), round(float(res[a, 2]), 4), (int(res[a, 3]) // 2) / 2, int(res[a, 3]) % 2)
                                     for a in order])
    sys.stdout.flush()

    # ---------------- 2. random instances: a log-uniform up to 2^40, L log-uniform up to 2^16
    la = rng.uniform(0, 40, n_random)
    a_arr = np.maximum(1, np.floor(2.0 ** la)).astype(np.int64)
    sign = rng.integers(0, 2, n_random) * 2 - 1
    a_arr = a_arr * sign
    L_arr = 2.0 ** rng.uniform(0, 16, n_random)
    ds_arr = rng.integers(0, 2, n_random).astype(np.int64)
    out = np.zeros((n_random, 4))
    run_batch(a_arr, ds_arr, L_arr, out)
    i1, i2 = report("random: |a| log-uniform <= 2^40, L log-uniform <= 2^16", a_arr, ds_arr, L_arr, out)
    print(f"  ({time.time()-t0:.1f}s)")
    sys.stdout.flush()

    # ---------------- 3. random instances in the non-trivial regime |a|^(1/3) <= L <= |a|, dyadic L
    m = n_random
    jj = rng.integers(1, 17, m)
    L3 = 2.0 ** jj
    gam = rng.uniform(1.0, 3.0, m)
    a3 = np.minimum(np.floor(L3 ** gam), 2.0 ** 46).astype(np.int64)
    a3 = np.maximum(a3 + rng.integers(-3, 4, m), 1)
    ds3 = rng.integers(0, 2, m).astype(np.int64)
    out3 = np.zeros((m, 4))
    run_batch(a3, ds3, L3, out3)
    report("random: L = 2^j (j <= 16), a ~ L^gamma, gamma in [1, 3]", a3, ds3, L3, out3)
    print(f"  ({time.time()-t0:.1f}s)")
    sys.stdout.flush()

    # ---------------- 4. structured instances
    A_list, D_list, L_list = [], [], []
    # 4a. numbers with many odd divisors (clusters at 0 mod 1)
    primes = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43]
    smooth = []
    def gen(i, v):
        if v > 2 ** 40:
            return
        if v > 2 ** 20:
            smooth.append(v)
        if i == len(primes):
            return
        e = 0
        x = v
        while x <= 2 ** 40 and e <= 4:
            gen(i + 1, x)
            x *= primes[i]
            e += 1
    gen(0, 1)
    smooth = sorted(set(smooth))
    # keep those with the largest number of divisors
    def ndiv(v):
        c = 1
        for p in primes:
            e = 0
            while v % p == 0:
                v //= p
                e += 1
            c *= e + 1
        return c
    smooth.sort(key=lambda v: -ndiv(v))
    for v in smooth[:400]:
        for j2 in range(4, 33):
            for sh in (0, 1, 2):
                A_list.append(v << sh)
                D_list.append(sh % 2)
                L_list.append(2.0 ** (j2 / 2.0))
    # 4b. stationary points: 2a/(2^ds k0^2) = p/q with small q  (a = p 2^ds k0^2 / (2 q) when integral)
    for _ in range(40000):
        j = int(rng.integers(3, 17))
        L = 2.0 ** j
        k0 = int(rng.integers(int(L) + 1, int(2 * L) + 1)) | 1
        p = int(rng.integers(1, 64))
        q = int(rng.choice([1, 1, 1, 2, 3, 4, 5, 6, 8]))
        ds = int(rng.integers(0, 2))
        num = p * (1 << ds) * k0 * k0
        a = num // (2 * q)
        if a < 1 or a > 2 ** 46:
            continue
        A_list.append(a + int(rng.integers(-1, 2)) if a > 1 else a)
        D_list.append(ds)
        L_list.append(L)
    # 4c. small |a| with all dyadic and half-dyadic L (the regime L >= |a|)
    for a in list(range(1, 200)) + [2 ** e for e in range(8, 41)] + [3 ** e for e in range(5, 26)]:
        for j2 in range(0, 33):
            for ds in (0, 1):
                A_list.append(a)
                D_list.append(ds)
                L_list.append(2.0 ** (j2 / 2.0))
    # 4d. perfect cubes / squares neighbourhoods (boundary L = |a|^(1/3), L = |a|^(1/2))
    for _ in range(20000):
        j = int(rng.integers(2, 17))
        L = float(2 ** j)
        e = rng.choice([1.0, 1.5, 2.0, 2.5, 3.0])
        a = int(round(L ** e)) + int(rng.integers(-2, 3))
        if a < 1 or a > 2 ** 48:
            continue
        A_list.append(a)
        D_list.append(int(rng.integers(0, 2)))
        L_list.append(L)
    a4 = np.array(A_list, dtype=np.int64)
    ds4 = np.array(D_list, dtype=np.int64)
    L4 = np.array(L_list, dtype=np.float64)
    out4 = np.zeros((len(a4), 4))
    run_batch(a4, ds4, L4, out4)
    j1, j2_ = report("structured: smooth a, stationary points, small a, boundary exponents", a4, ds4, L4, out4)
    print(f"  ({time.time()-t0:.1f}s)")

    # exact recheck of the worst structured / random instances (dyadic or exactly representable L only)
    print("[exact rechecks with fractions.Fraction]")
    for (arr_a, arr_d, arr_L, o, idxs) in ((a_arr, ds_arr, L_arr, out, list(i1[:2]) + list(i2[:2])),
                                             (a4, ds4, L4, out4, list(j1[:3]) + list(j2_[:3]))):
        for i in idxs:
            Lf = Fraction(float(arr_L[i]))   # exact value of the double
            if arr_L[i] > 5000:
                continue
            eL, nd, n = exact_instance(int(arr_a[i]), int(arr_d[i]), Lf)
            print(f"  a = {arr_a[i]}, ds = {arr_d[i]}, L = {arr_L[i]!r}: exact E_L = {float(eL):.9f} (float run "
                  f"{o[i,0]:.9f}), exact N D_N = {float(nd):.9f} (float run {o[i,1]:.9f}), N = {n}")
    print(f"total time {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
