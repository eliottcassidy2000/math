#!/usr/bin/env python3
"""S93 coordinator: finite tests of the lemmas of EXPONENT_SKETCH.txt (tests, not proofs).

(T1) block law 2q(A 2^m + B) = 2q(A) + 2q(B) + A rev_m(B) + (A mod 2) B_(m-1), random and exhaustive small cases.
(T2) every colliding pair of n-bit words: collision equation c X - P Xi = gamma, |gamma| <= 2^w, X and Xi odd;
     for run-free words (no run of more than R equal bits): |s - t'| <= R + 3.
(T3) height histogram of all ordered pairs with X, Xi odd: N(M in [2^i, 2^(i+1))) for w <= WMAX, with the
     comparison values 6^(w/2), 8^(w/2); S_w^odd.
(T4) half-split bookkeeping: with e = 2^ds k' a' - k a:  k b' - k' b = 2^s2 e and |e| < |k| + 2^ds k'.
(T5) negative moment: sum r^2 |a|^(-1/2) against 42 * 6^s 2^(-s/2); weighted share of |D(eps)| < 2^(s-j).
Usage: python3 exponent_lemma_tests.py [WMAX] [NMAX]
"""
import math
import random
import sys

import numpy as np
from numba import njit, prange


def q_of(h):
    q = 0
    b = 0
    while (h >> b) > 0:
        if (h >> b) & 1:
            q += ((h >> (b + 1)) + 1) >> 1
        b += 1
    return q


def rev(v, w):
    r = 0
    for i in range(w):
        if (v >> i) & 1:
            r |= 1 << (w - 1 - i)
    return r


def t1():
    ok = True
    for A in range(0, 64):
        for m in range(1, 7):
            for B in range(1 << m):
                lhs = 2 * q_of((A << m) + B)
                rhs = 2 * q_of(A) + 2 * q_of(B) + A * rev(B, m) + (A & 1) * ((B >> (m - 1)) & 1)
                ok &= lhs == rhs
    rng = random.Random(93)
    for _ in range(20000):
        m = rng.randint(1, 40)
        A = rng.getrandbits(rng.randint(1, 40))
        B = rng.getrandbits(m)
        lhs = 2 * q_of((A << m) + B)
        rhs = 2 * q_of(A) + 2 * q_of(B) + A * rev(B, m) + (A & 1) * ((B >> (m - 1)) & 1)
        ok &= lhs == rhs
    print("(T1) block law, exhaustive A < 64, m <= 6 and 20000 random cases:", ok, flush=True)
    return ok


@njit(cache=True)
def q_nb(h):
    q = 0
    b = 0
    while (h >> b) > 0:
        if (h >> b) & 1:
            q += ((h >> (b + 1)) + 1) >> 1
        b += 1
    return q


@njit(cache=True)
def rev_nb(v, w):
    r = 0
    for i in range(w):
        if (v >> i) & 1:
            r |= 1 << (w - 1 - i)
    return r


@njit(cache=True)
def maxrun(h, n):
    best = 1
    cur = 1
    for i in range(1, n):
        if ((h >> i) & 1) == ((h >> (i - 1)) & 1):
            cur += 1
            if cur > best:
                best = cur
        else:
            cur = 1
    return best


@njit(cache=True)
def t2_scan(n, keys, runs, out):
    """keys sorted: value << 32 | index.  out[0] = pairs, out[1] = equation failures, out[2] = parity failures,
    out[3 + R] = max |s - t'| over pairs of words whose runs are <= R, for R = 1..12."""
    N = keys.shape[0]
    lo = 1 << (n - 1)
    mask = (1 << 32) - 1
    i = 0
    while i < N:
        v = keys[i] >> 32
        j = i + 1
        while j < N and (keys[j] >> 32) == v:
            j += 1
        for a in range(i, j):
            ha = lo + (keys[a] & mask)
            for b in range(i, a):
                hb = lo + (keys[b] & mask)
                d = ha ^ hb
                top = 0
                dd = d
                while dd > 1:
                    dd >>= 1
                    top += 1
                s = 0
                while ((d >> s) & 1) == 0:
                    s += 1
                tp = n - 1 - top
                w = top - s + 1
                xa = (ha >> s) & ((1 << w) - 1)
                xb = (hb >> s) & ((1 << w) - 1)
                P = ha >> (top + 1)
                Z = ha & ((1 << s) - 1)
                c = (1 << (s + 1)) - rev_nb(Z, s)
                X = xa - xb
                Xi = rev_nb(xa, w) - rev_nb(xb, w)
                gamma = (2 * (q_nb(xa) - q_nb(xb)) + (P & 1) * (((xa >> (w - 1)) & 1) - ((xb >> (w - 1)) & 1))
                         + (c & 1) * ((xa & 1) - (xb & 1)))
                out[0] += 1
                if c * X - P * Xi != gamma or abs(gamma) > (1 << w):
                    out[1] += 1
                if (X & 1) == 0 or (Xi & 1) == 0:
                    out[2] += 1
                r = max(runs[keys[a] & mask], runs[keys[b] & mask])
                diff = abs(s - tp)
                for R in range(r, 13):
                    if diff > out[3 + R]:
                        out[3 + R] = diff
        i = j


@njit(parallel=True, cache=True)
def t2_fill(n, keys, runs):
    lo = 1 << (n - 1)
    for i in prange(keys.shape[0]):
        h = lo + i
        keys[i] = ((2 * h + 1 - 2 * q_nb(h)) << 32) | i
        runs[i] = maxrun(h, n)


def t2(nmax):
    ok = True
    for n in range(8, nmax + 1, 2):
        N = 1 << (n - 1)
        keys = np.empty(N, np.int64)
        runs = np.empty(N, np.int64)
        t2_fill(n, keys, runs)
        keys.sort()
        out = np.zeros(16, np.int64)
        t2_scan(n, keys, runs, out)
        gaps = [int(out[3 + R]) for R in range(1, 9)]
        good = out[1] == 0 and out[2] == 0 and all(gaps[R - 1] <= R + 3 for R in range(1, 9))
        ok &= bool(good)
        print(f"(T2) n = {n}: {int(out[0])} unordered colliding pairs; equation failures {int(out[1])}, parity failures "
              f"{int(out[2])}; max |s - t'| among pairs with runs <= R, R = 1..8: {gaps} (bound R + 3): {bool(good)}",
              flush=True)
    return ok


@njit(cache=True)
def gcd_nb(a, b):
    while b:
        a, b = b, a % b
    return a


@njit(parallel=True, cache=True)
def t3_hist(w, hist, sums):
    size = 1 << w
    rv = np.zeros(size, np.int64)
    for x in range(size):
        rv[x] = rev_nb(x, w)
    for x in prange(size):
        row = x                      # one row per x: no two threads write the same row
        for y in range(size):
            X = x - y
            if (X & 1) == 0:
                continue
            Xi = rv[x] - rv[y]
            if (Xi & 1) == 0:
                continue
            ax = abs(X)
            ai = abs(Xi)
            g = gcd_nb(ax, ai)
            m = (ax if ax > ai else ai) // g
            lg = 0
            while (m >> (lg + 1)) > 0:
                lg += 1
            hist[row, lg] += 1
            sums[row] += 1.0 / m


def t3(wmax):
    for w in range(4, wmax + 1):
        hist = np.zeros((1 << w, w + 2), np.int64)
        sums = np.zeros(1 << w, np.float64)
        t3_hist(w, hist, sums)
        h = hist.sum(axis=0)
        S = sums.sum()
        six = 6 ** (w / 2)
        ratios = [h[i] / (2 ** i) / six for i in range(w + 1)]
        worst = max(ratios)
        print(f"(T3) w = {w:2d}: S_odd = {S:.1f} = {S / six:.4f} * 6^(w/2) = {S / 8 ** (w / 2):.5f} * 8^(w/2); "
              f"max over i of N(M in [2^i, 2^(i+1)))/2^i = {worst:.4f} * 6^(w/2) at i = {ratios.index(worst)}",
              flush=True)
        if w == wmax:
            print("     N(M ~ 2^i)/2^i / 6^(w/2) for i = 0..w:", " ".join(f"{r:.3f}" for r in ratios))


def t4():
    ok = True
    rng = random.Random(7)
    for w in (9, 10, 11, 12, 13):
        s1 = (w + 1) // 2
        s2 = w // 2
        ds = s1 - s2
        bad = 0
        cnt = 0
        for _ in range(200000):
            x = rng.getrandbits(w)
            xp = rng.getrandbits(w)
            X = x - xp
            Xi = rev(x, w) - rev(xp, w)
            if X % 2 == 0 or Xi % 2 == 0:
                continue
            v, y = x >> s1, x & ((1 << s1) - 1)
            vp, yp = xp >> s1, xp & ((1 << s1) - 1)
            a = rev(y, s1) - rev(yp, s1)
            b = y - yp
            ap = v - vp
            bp = rev(v, s2) - rev(vp, s2)
            if X != ap * (1 << s1) + b or Xi != a * (1 << s2) + bp:
                bad += 1
                continue
            g = math.gcd(X, Xi)
            k, kp = X // g, Xi // g
            if kp < 0:
                k, kp = -k, -kp
            e = (1 << ds) * kp * ap - k * a
            if k * bp - kp * b != (1 << s2) * e or abs(e) >= abs(k) + (1 << ds) * kp:
                bad += 1
            cnt += 1
        ok &= bad == 0
        print(f"(T4) w = {w} (s1, s2) = ({s1}, {s2}): {cnt} random odd pairs, failures of the split identities: {bad}",
              flush=True)
    return ok


def t5():
    ok = True
    for s in range(2, 13):
        n = 3 ** s
        idx = np.arange(n, dtype=np.int64)
        D = np.zeros(n, np.int64)
        z = np.zeros(n, np.int64)
        t = idx.copy()
        for i in range(s):
            e = (t % 3) - 1
            t //= 3
            D += e << i
            z += (e == 0)
        wt = 4.0 ** z
        nz = D != 0
        mom = float((wt[nz] / np.sqrt(np.abs(D[nz]))).sum())
        bound = 42 * 6.0 ** s * 2 ** (-s / 2)
        shares = []
        good = mom <= bound
        for j in range(0, s + 1):
            share = float(wt[np.abs(D) < (1 << (s - j))].sum()) / 6.0 ** s
            shares.append(share)
            good &= share <= (5 / 3) * (2 / 3) ** j + 1e-12
        ok &= bool(good)
        print(f"(T5) s = {s:2d}: sum r^2 |a|^(-1/2) = {mom / (6.0 ** s * 2 ** (-s / 2)):.3f} * 6^s 2^(-s/2) (bound 42); "
              f"shares of |D| < 2^(s-j), j = 0..3: {' '.join(f'{v:.4f}' for v in shares[:4])} "
              f"(bounds {' '.join(f'{(5 / 3) * (2 / 3) ** j:.4f}' for j in range(4))}): {bool(good)}", flush=True)
    return ok


def main():
    wmax = int(sys.argv[1]) if len(sys.argv) > 1 else 13
    nmax = int(sys.argv[2]) if len(sys.argv) > 2 else 22
    if len(sys.argv) > 3 and sys.argv[3] == "t3only":
        t3(wmax)
        return
    ok = t1()
    ok &= t4()
    ok &= t5()
    ok &= t2(nmax)
    t3(wmax)
    print("ALL CHECKS PASSED" if ok else "SOME CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
