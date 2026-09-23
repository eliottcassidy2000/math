#!/usr/bin/env python3
"""FINITE-EXACT engines for the transversality foundry (collatz-procgen-20260922).

All arithmetic is exact (Python/gmpy2 integers).  Each function returns a dict
of results and raises AssertionError if an identity it checks fails.

  erdos_scan        ternary digit 2 in 2^K (Erdos), 0 <= K <= Kmax
  ternary_profile   nonzero digits, digit-2 counts and runs of 2^K in base 3
  low_digit_classes #{K mod 2*3^(j-1): lowest j ternary digits of 2^K in {0,1}}
  kappa_identity    v_3(2^K - r) = [K = e mod 2] (1 + v_3(K - kappa(r))), kappa in Z_3
  v3_scan           max_K v_3(2^K - r) against log_3 K
  pow3_binary       nonzero bits and zero runs of 3^a
  base5_scan        base-5 analogue of the Erdos question
  beatty_landing    2^(K0(s)) w mod 3^D along the loop halving counts K0(s)=floor((s+1)log2 3)
"""
import math
import re

import gmpy2
from gmpy2 import mpz


# ------------------------------------------------------------------ helpers
def v_p(n, p):
    n = mpz(n)
    if n == 0:
        return math.inf
    return int(gmpy2.remove(n, p)[1])


def digits(n, b):
    """Digit string of n >= 0 in base b, most significant first."""
    return gmpy2.digits(mpz(n), b)


def longest_run(s, ch=None):
    """Longest run of the character ch (or of any character if ch is None) in s."""
    chars = [ch] if ch is not None else sorted(set(s))
    best = 0
    for c in chars:
        others = "".join(sorted(set(s) - {c}))
        if not others:
            best = max(best, len(s))
            continue
        parts = re.split("[" + re.escape(others) + "]", s)
        best = max(best, max(len(x) for x in parts))
    return best


# ------------------------------------------------------------------ Erdos
def erdos_scan(Kmax, low=60):
    """For 0 <= K <= Kmax: does 2^K contain the ternary digit 2?  Low digits are
    read from 2^K mod 3^low; the full expansion is converted only when needed."""
    M = mpz(3) ** low
    r = mpz(1)
    no2 = []
    deepest = (0, 0)          # (depth of the lowest digit 2, K) over K with >= low digits
    full_checks = 0
    for K in range(Kmax + 1):
        s = digits(r, 3)[::-1]  # least significant first (length <= low)
        pos = s.find("2")
        ndig = int(K * math.log(2) / math.log(3)) + 2
        if pos < 0 and ndig > len(s):
            # either 2^K < 3^low (then s is the full expansion) or look higher
            full = digits(mpz(2) ** K, 3)[::-1]
            full_checks += 1
            pos = full.find("2")
        if pos < 0:
            no2.append(K)
        elif pos + 1 > deepest[0]:
            deepest = (pos + 1, K)
        r = (2 * r) % M
    return dict(Kmax=Kmax, no2=no2, deepest=deepest, full_checks=full_checks)


def ternary_profile(Kmax, K0=50):
    """Nonzero-digit and digit-2 proportions and runs in the full ternary expansion of 2^K."""
    worst_nnz = (2.0, None)
    worst_two = (2.0, None)
    worst_nnz_abs = (10 ** 9, None)
    run0 = (0, None)
    run_any = (0, None)
    top0 = (0, None)       # zeros right after the leading digit (Ren-Roettger: O(log K))
    x = mpz(1) << K0
    for K in range(K0, Kmax + 1):
        s = digits(x, 3)
        L = len(s)
        t0 = len(s[1:]) - len(s[1:].lstrip("0"))
        if t0 > top0[0]:
            top0 = (t0, K)
        # bottom run (LTE): K even -> 2^K = 1 + 3^(1+v_3(K)) u: zeros above the unit digit = v_3(K);
        # K odd -> 2^K = -1 mod 3^(1+v_3(K)): the lowest 1+v_3(K) digits are 2
        low = s[::-1]
        if K % 2 == 0:
            z = len(low[1:]) - len(low[1:].lstrip("0"))
            assert low[0] == "1" and z == v_p(K, 3), (K, z)
        else:
            t2 = len(low) - len(low.lstrip("2"))
            assert t2 == 1 + v_p(K, 3), (K, t2)
        nnz = L - s.count("0")
        two = s.count("2")
        if nnz / L < worst_nnz[0]:
            worst_nnz = (nnz / L, K)
        if two / L < worst_two[0]:
            worst_two = (two / L, K)
        if nnz < worst_nnz_abs[0]:
            worst_nnz_abs = (nnz, K)
        r0 = longest_run(s, "0")
        if r0 > run0[0]:
            run0 = (r0, K)
        ra = longest_run(s)
        if ra > run_any[0]:
            run_any = (ra, K)
        x <<= 1
    return dict(K0=K0, Kmax=Kmax, worst_nnz=worst_nnz, worst_two=worst_two, worst_nnz_abs=worst_nnz_abs,
                run0=run0, run_any=run_any, top0=top0)


def low_digit_classes(jmax, brute_max=8):
    """N_j = #{K mod 2*3^(j-1) : lowest j ternary digits of 2^K all in {0,1}}.
    Computed by exact lifting (each class lifts to 3 classes, digit j of 2^K then
    takes each value once), and by brute force for j <= brute_max.  PROVED: N_j = 2^(j-1)."""
    S = [0]  # j = 1: K even  (2^K = 1 mod 3)
    out = {1: 1}
    for j in range(1, jmax):
        m = 2 * 3 ** (j - 1)          # ord of 2 mod 3^j
        mod = 3 ** (j + 1)
        S_new = []
        for k in S:
            for t in range(3):
                K = k + t * m
                d = (pow(2, K, mod) // 3 ** j) % 3
                if d != 2:
                    S_new.append(K)
        S = S_new
        out[j + 1] = len(S)
        assert len(S) == 2 ** j, (j + 1, len(S))
    for j in range(1, brute_max + 1):
        m = 2 * 3 ** (j - 1)
        cnt = 0
        for K in range(m):
            s = digits(pow(2, K, 3 ** j), 3).zfill(j)
            if "2" not in s:
                cnt += 1
        assert cnt == out[j] == 2 ** (j - 1), (j, cnt, out[j])
    return out


# ------------------------------------------------------------------ 3-adic logarithm and kappa
def log3adic(x_num, x_den, D):
    """3-adic log of the unit x = x_num/x_den = 1 mod 3, modulo 3^D (returned as int)."""
    M = 3 ** (D + 8)
    x = (x_num * pow(x_den, -1, M)) % M
    assert x % 3 == 1
    t = (x - 1) % M         # t = x - 1, v_3(t) >= 1
    res = 0
    tn = 1
    n = 1
    while True:
        tn = (tn * t) % M
        k = v_p(n, 3)
        if n - k >= D + 2 and n > 2 * D + 10:
            break
        m = n // 3 ** k
        # t^n / n = (t^n / 3^k) / m ; t^n divisible by 3^n >= 3^k
        term = (tn // 3 ** k) * pow(m, -1, M) % M
        res = (res + (term if n % 2 == 1 else -term)) % M
        n += 1
    return res % (3 ** D)


def kappa(r_num, r_den, D=40):
    """Returns (e, kappa mod 3^D) with 2^K - r = 0 mod 3 iff K = e mod 2, and
    v_3(2^K - r) = 1 + v_3(K - kappa) for such K (kappa in Z_3)."""
    Mx = 3 ** (D + 10)
    r = (r_num * pow(r_den, -1, Mx)) % Mx
    assert r % 3 != 0
    e = 0 if r % 3 == 1 else 1
    r1 = r if e == 0 else (-r) % Mx        # r1 = 1 mod 3
    Lr = log3adic(r1, 1, D + 2)
    Lm2 = log3adic((-2) % Mx, 1, D + 2)    # log(-2), valuation exactly 1
    assert Lm2 % 3 == 0 and (Lm2 // 3) % 3 != 0
    assert Lr % 3 == 0
    k = (Lr // 3) * pow(Lm2 // 3, -1, 3 ** (D + 1)) % (3 ** D)
    return e, k


def kappa_identity(rs, Kmax=2000, D=40):
    """Check v_3(2^K b - a) against the kappa formula for 0 <= K <= Kmax."""
    res = {}
    for (a, b) in rs:
        e, k = kappa(a, b, D)
        worst = 0
        for K in range(Kmax + 1):
            direct = v_p(mpz(2) ** K * b - a, 3)
            if K % 2 != e:
                pred = 0
            else:
                pred = 1 + min(v_p(K - k, 3), D) if K != k else math.inf
            assert direct == pred, (a, b, K, direct, pred)
            worst = max(worst, direct)
        res[(a, b)] = dict(e=e, kappa_low_digits=digits(k, 3)[::-1][:12], max_v3=worst)
    return res


def v3_scan(rs, Kmax, D=80):
    """max over 1 <= K <= Kmax of v_3(2^K b - a) and of v_3(...) - log_3 K."""
    M = mpz(3) ** D
    out = {}
    for (a, b) in rs:
        p2 = mpz(1)
        best = (0, None)
        excess = (-1e9, None)
        for K in range(1, Kmax + 1):
            p2 = (2 * p2) % M
            n = (p2 * b - a) % M
            v = D if n == 0 else v_p(n, 3)
            assert v < D, "precision exhausted"
            if v > best[0]:
                best = (v, K)
            ex = v - math.log(K, 3)
            if ex > excess[0]:
                excess = (ex, K)
        out[(a, b)] = dict(max_v3=best, max_excess_over_log3K=excess)
    return out


# ------------------------------------------------------------------ powers of 3 in binary
def pow3_binary(amax, a0=20):
    worst = (2.0, None)
    worst_abs = (10 ** 9, None)
    run0 = (0, None)
    low_run = (0, None)
    x = mpz(3) ** a0
    for a in range(a0, amax + 1):
        L = x.bit_length()
        pc = gmpy2.popcount(x)
        if pc / L < worst[0]:
            worst = (pc / L, a)
        if pc < worst_abs[0]:
            worst_abs = (pc, a)
        s = digits(x, 2)
        r0 = longest_run(s, "0")
        if r0 > run0[0]:
            run0 = (r0, a)
        lr = v_p(x - 1, 2) - 1  # zeros just above the lowest 1 bit
        if lr > low_run[0]:
            low_run = (lr, a)
        x *= 3
    # LTE check: v_2(3^a - 1) = 1 (a odd), 2 + v_2(a) (a even)
    for a in range(1, 3000):
        v = v_p(mpz(3) ** a - 1, 2)
        assert v == (1 if a % 2 else 2 + v_p(a, 2))
    return dict(a0=a0, amax=amax, worst_ratio=worst, worst_abs=worst_abs, run0=run0, low_run=low_run)


def small_popcount_pow3(amax, thresh=22):
    """Exponents a <= amax with popcount(3^a) <= thresh (compare Dimitrov-Howe)."""
    out = []
    x = mpz(1)
    for a in range(amax + 1):
        if gmpy2.popcount(x) <= thresh:
            out.append((a, int(gmpy2.popcount(x))))
        x *= 3
    return out


# ------------------------------------------------------------------ base 5
def base5_scan(Kmax, low=60):
    """K <= Kmax with all base-5 digits of 2^K in {0,1} (the 5-analogue of Erdos)."""
    M = mpz(5) ** low
    r = mpz(1)
    hits = []
    for K in range(Kmax + 1):
        s = digits(r, 5)
        ndig = int(K * math.log(2) / math.log(5)) + 2
        bad = any(c in "234" for c in s)
        if not bad and ndig > len(s):
            s2 = digits(mpz(2) ** K, 5)
            bad = any(c in "234" for c in s2)
        if not bad:
            hits.append(K)
        r = (2 * r) % M
    return dict(Kmax=Kmax, only01=hits)


# ------------------------------------------------------------------ Q2 landing along loop clocks
def loop_clocks(smax):
    """K0(s) = floor((s+1) log2 3) for 0 <= s <= smax, exactly: log2 3 is replaced by a
    600-bit rational L with |L - log2 3| < 2^-600, and every (s+1)L is checked to lie
    farther than 2^-500 from an integer (so the floor is that of the true value)."""
    import mpmath
    mpmath.mp.prec = 700
    L = int(mpmath.nint(mpmath.log(3, 2) * mpmath.mpf(2) ** 600))
    one = 1 << 600
    margin = 1 << 100
    out = []
    acc = 0
    for s in range(smax + 1):
        acc += L
        q, rem = divmod(acc, one)
        assert margin < rem < one - margin, s
        out.append(q)
    # spot check against exact integer comparison 2^K < 3^(s+1) < 2^(K+1)
    for s in list(range(0, 60)) + [smax // 2, smax]:
        K = out[s]
        assert (1 << K) < 3 ** (s + 1) < (1 << (K + 1)), s
    return out


def beatty_landing(ws, smax, D=13):
    """Classes 2^(K0(s)) w mod 3^D along the loop halving counts K0(s) = floor((s+1)log2 3)
    (loops_and_escapes sec. 1).  Counts landings in the 3-adic balls of radius 3^-D around the
    hostile points 1 and 1/2, i.e. v_3(2^K0(s) w - h) >= D."""
    M = 3 ** D
    K0 = loop_clocks(smax)
    half = pow(2, -1, M)
    out = {}
    for (a, b) in ws:
        w = a * pow(b, -1, M) % M
        y = pow(2, K0[0], M) * w % M
        hits = {"1": 0, "1/2": 0}
        first = {"1": None, "1/2": None}
        for s in range(smax + 1):
            if s > 0:
                y = (y << (K0[s] - K0[s - 1])) % M
            if y == 1:
                hits["1"] += 1
                if first["1"] is None:
                    first["1"] = s
            if y == half:
                hits["1/2"] += 1
                if first["1/2"] is None:
                    first["1/2"] = s
        out[(a, b)] = dict(hits=hits, first=first)
    return dict(D=D, smax=smax, res=out, K0_head=K0[:12])


if __name__ == "__main__":
    print(erdos_scan(2000))
    print(low_digit_classes(12))
    print(kappa_identity([(1, 1), (1, 2), (5, 1), (7, 1), (1, 5), (43, 32)], 500))
