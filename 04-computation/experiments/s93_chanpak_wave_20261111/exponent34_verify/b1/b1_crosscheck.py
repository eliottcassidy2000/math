#!/usr/bin/env python3
"""Second, deliberately naive implementation (pure Python, math.gcd, Fractions, bit lists) of the objects of the
half split of THM-5141, used to cross-check b1_deg_checks.py, plus sampled tests of (D1), (D2) at larger w.

Part 1 (w = 4..WSLOW): from the raw definitions, with words handled as bit lists:
   - every ordered pair (x, x'): X, Xi, oddness, height M (Fraction arithmetic for 1/M);
   - Deg_K, Deg'_K for all box points and all integer K, by direct enumeration of the box;
   - N(K), S_w, right side of (1), right side of (2), the two displayed lines of (B1).
Part 2 (sampled, numba): Deg_K(P_u), Deg'_K(P_v) for random and structured P_u, P_v at w = 16, 18, 20, 21, 22, 24.

Usage: python b1_crosscheck.py WSLOW [sample]
"""
import os
import sys
import math
import random
import time
from fractions import Fraction
from math import gcd

os.environ.setdefault("NUMBA_NUM_THREADS", "4")


def bits(x, w):
    return [(x >> i) & 1 for i in range(w)]          # bits[i] = bit i


def val(bl):
    return sum(b << i for i, b in enumerate(bl))


def rev_bits(x, w):
    return val(list(reversed(bits(x, w))))


def slow_check(w):
    s1 = -(-w // 2)
    s2 = w // 2
    n = 1 << w
    n1 = 1 << s1
    n2 = 1 << s2
    # r tables as dictionaries, from the definition
    r1 = {}
    for u in range(n1):
        for u2 in range(n1):
            key = (u - u2, rev_bits(u, s1) - rev_bits(u2, s1))
            r1[key] = r1.get(key, 0) + 1
    r2 = {}
    for v in range(n2):
        for v2 in range(n2):
            key = (v - v2, rev_bits(v, s2) - rev_bits(v2, s2))
            r2[key] = r2.get(key, 0) + 1
    assert sum(c * c for c in r1.values()) == 6 ** s1
    assert sum(c * c for c in r2.values()) == 6 ** s2

    def height(X, Xi):
        return max(abs(X), abs(Xi)) // gcd(X, Xi)

    # direct enumeration of pairs
    cntM = {}
    cntF = {}
    npairs = 0
    R = [rev_bits(x, w) for x in range(n)]
    for x in range(n):
        for x2 in range(n):
            X = x - x2
            Xi = R[x] - R[x2]
            if X % 2 and Xi % 2:
                M = height(X, Xi)
                assert M % 2 == 1
                cntM[M] = cntM.get(M, 0) + 1
                g = gcd(X, Xi)
                mx = max(abs(X), abs(Xi))
                gg = gcd(g, mx)
                kf = (g // gg, mx // gg)                 # the summand gcd/max of the definition, reduced
                cntF[kf] = cntF.get(kf, 0) + 1
                npairs += 1
    assert npairs == 4 ** (w - 1)
    S = sum((Fraction(c * num, den) for (num, den), c in cntF.items()), Fraction(0))
    assert S == sum((Fraction(c, M) for M, c in cntM.items()), Fraction(0))

    def N(K):
        return sum(c for M, c in cntM.items() if M <= K)

    # degrees from the definition
    boxU = [(a, b) for a in range(-(n1 - 1), n1) for b in range(-(n1 - 1), n1) if b % 2]
    boxV = [(a_, b_) for a_ in range(-(n2 - 1), n2) for b_ in range(-(n2 - 1), n2) if b_ % 2]
    Mcap = n
    degU = {}
    degV = {P: [0] * (Mcap + 1) for P in boxV}
    Nsplit = [0] * (Mcap + 1)
    for (a, b) in boxU:
        h = [0] * (Mcap + 1)
        ru = r1.get((a, b), 0)
        for (a_, b_) in boxV:
            X = a_ * n1 + b
            Xi = a * n2 + b_
            M = height(X, Xi)
            h[M] += 1
            degV[(a_, b_)][M] += 1
            Nsplit[M] += ru * r2.get((a_, b_), 0)
        degU[(a, b)] = h
    # cumulative
    for P in degU:
        h = degU[P]
        for K in range(1, Mcap + 1):
            h[K] += h[K - 1]
    for P in degV:
        h = degV[P]
        for K in range(1, Mcap + 1):
            h[K] += h[K - 1]
    for K in range(1, Mcap + 1):
        Nsplit[K] += Nsplit[K - 1]
    bad = 0
    worst = {"D1": 0.0, "D1p": 0.0, "D2": 0.0, "D2p": 0.0, "one": 0.0, "B1N": 0.0}
    out_rows = []
    for K in range(1, Mcap):
        NK = N(K)
        if NK != Nsplit[K]:
            bad += 1
        mU = max(h[K] for h in degU.values())
        mV = max(h[K] for h in degV.values())
        if mU > 4 * K * K or mV > 8 * K * K:
            bad += 1
        worst["D1"] = max(worst["D1"], mU / (4 * K * K))
        worst["D1p"] = max(worst["D1p"], mV / (8 * K * K))
        rhs1 = sum(r1[P] ** 2 * degU[P][K] for P in r1 if P[1] % 2) + sum(r2[P] ** 2 * degV[P][K] for P in r2 if P[1] % 2)
        if 2 * NK > rhs1:
            bad += 1
        worst["one"] = max(worst["one"], 2 * NK / rhs1)
        if NK > 8 * K * K * 6 ** s1:
            bad += 1
        if K <= n2:
            bnd = 3 * n1 * (3 + 2 * math.log(K))
            if mU > bnd or mV > bnd:
                bad += 1
            worst["D2"] = max(worst["D2"], mU / bnd)
            worst["D2p"] = max(worst["D2p"], mV / bnd)
            bb = 6 ** s1 * min(8 * K * K, bnd)
            if NK > bb:
                bad += 1
            worst["B1N"] = max(worst["B1N"], NK / bb)
        if K & (K - 1) == 0:
            out_rows.append((K, mU, mV, NK, rhs1 // 2 if rhs1 % 2 == 0 else rhs1 / 2))
    # (2)
    rhs2 = Fraction(0)
    for i in range(0, 6 * w):
        rhs2 += Fraction(N(1 << (i + 1)), 1 << i)
    rhs2 += Fraction(npairs, 1 << (6 * w - 1))
    if S > rhs2:
        bad += 1
    i0 = s1 // 2
    line1 = 6 ** s1 * (sum(32 * 2 ** i for i in range(i0 + 1)) + Fraction(3 * n1 * (3 + 2 * s2), 2 ** i0)) + Fraction(2 * 4 ** w, 2 ** s2)
    # the derivation of line 1 from (2): i <= i0 with 8K^2, i0 < i <= s2-1 with (D2), i >= s2 with 4^w
    derived = Fraction(0)
    for i in range(0, s2):
        K = 1 << (i + 1)
        if i <= i0:
            derived += Fraction(6 ** s1 * 8 * K * K, 1 << i)
        else:
            derived += Fraction(6 ** s1 * 3 * n1 * (3 + 2 * s2), 1 << i)
    derived += Fraction(2 * 4 ** w, 2 ** s2)
    line2 = (6 * math.sqrt(2)) ** s1 * (82 + 12 * s2) + 2 ** (1.5 * w + 2)
    if not (S <= rhs2 <= derived <= line1) or float(line1) > line2:
        bad += 1
    print(f"w = {w}: S_w = {S} = {float(S):.6f}; N(K) split == direct; worst ratios: "
          f"D1 {worst['D1']:.4f}, D1' {worst['D1p']:.4f}, D2 {worst['D2']:.4f}, D2' {worst['D2p']:.4f}, "
          f"(1) {worst['one']:.4f}, N(K)/bound {worst['B1N']:.5f}; RHS(2) = {float(rhs2):.3f} <= derived {float(derived):.1f} "
          f"<= line1 {float(line1):.1f} <= line2 {line2:.1f}; failures: {bad}")
    for row in out_rows:
        print(f"      K = {row[0]:5d}: max Deg = {row[1]:5d}, max Deg' = {row[2]:5d}, N(K) = {row[3]:8d}, RHS(1) = {row[4]}")
    return bad


def sampled(seed=20261111):
    import numpy as np
    from numba import njit, prange

    @njit(cache=True)
    def g_(a, b):
        while b != 0:
            a, b = b, a % b
        return a

    @njit(parallel=True, cache=True)
    def deg_points(s_out, s_in, P, Q, Kmax):
        """cumulative Deg_K (K = 0..Kmax) for the outer points (P[i], Q[i])."""
        n_out = 1 << s_out
        n_in = 1 << s_in
        m = P.shape[0]
        out = np.zeros((m, Kmax + 1), dtype=np.int64)
        for i in prange(m):
            p = P[i]
            q = Q[i]
            for pp in range(-(n_in - 1), n_in):
                A = pp * n_out + q
                aA = abs(A)
                for qq in range(-(n_in - 1), n_in, 2):      # n_in is even, so these are the odd q'
                    B = p * n_in + qq
                    aB = abs(B)
                    g = g_(aA, aB)
                    if aA > aB:
                        M = aA // g
                    else:
                        M = aB // g
                    if M <= Kmax:
                        out[i, M] += 1
            for K in range(1, Kmax + 1):
                out[i, K] += out[i, K - 1]
        return out

    rng = random.Random(seed)
    total_bad = 0
    for w, nrand in ((16, 3000), (17, 1500), (18, 1500), (20, 300), (21, 200), (22, 120), (24, 24)):
        t0 = time.time()
        s1 = -(-w // 2)
        s2 = w // 2
        n1 = 1 << s1
        n2 = 1 << s2
        for name, s_out, s_in, cD1 in (("Deg", s1, s2, 4), ("Deg'", s2, s1, 8)):
            no = 1 << s_out
            pts = set()
            specials_p = [0, 1, -1, 2, -2, 3, -3, 4, 5, 6, 7, 8, 12, 15, 16, 24, 30, 60, no - 1, -(no - 1), no // 2, no // 2 - 1,
                          no // 2 + 1, no // 3, no // 4, 105, 315, 945]
            specials_q = [1, -1, 3, -3, 5, 7, 9, 15, 45, 105, 315, no - 1, -(no - 1), no // 2 - 1, no // 2 + 1]
            for p in specials_p:
                for q in specials_q:
                    if abs(p) < no and abs(q) < no and q % 2:
                        pts.add((p, q))
            pts = sorted(pts)
            if w >= 22:
                pts = pts[::7]
            if w >= 24:
                pts = pts[::3]
            for _ in range(nrand):
                p = rng.randrange(-(no - 1), no)
                q = rng.randrange(-(no - 1), no) | 1
                if abs(q) >= no:
                    q = no - 1
                pts.append((p, q))
            P = np.array([t[0] for t in pts], dtype=np.int64)
            Q = np.array([t[1] for t in pts], dtype=np.int64)
            Kmax = 4 * n2
            D = deg_points(s_out, s_in, P, Q, Kmax)
            mx = D.max(axis=0)
            bad = 0
            w1 = 0.0
            w2 = (0.0, 0)
            for K in range(1, Kmax + 1):
                if mx[K] > cD1 * K * K:
                    bad += 1
                w1 = max(w1, mx[K] / (cD1 * K * K))
                if K <= n2:
                    bnd = 3 * n1 * (3 + 2 * math.log(K))
                    if mx[K] > bnd:
                        bad += 1
                    if mx[K] / bnd > w2[0]:
                        w2 = (mx[K] / bnd, K)
            total_bad += bad
            print(f"w = {w} ({name}, {len(pts)} points: structured + {nrand} random): (D1) worst ratio {w1:.4f}; "
                  f"(D2) worst ratio {w2[0]:.4f} at K = {w2[1]}; max Deg at K = 2^s2: {int(mx[n2])}; violations {bad}", flush=True)
        print(f"   [{time.time() - t0:.1f} s]", flush=True)
    return total_bad


def main():
    wslow = int(sys.argv[1])
    bad = 0
    for w in range(2, wslow + 1):
        bad += slow_check(w)
        sys.stdout.flush()
    if len(sys.argv) > 2 and sys.argv[2] == "sample":
        bad += sampled()
    print("TOTAL FAILURES:", bad)


if __name__ == "__main__":
    main()
