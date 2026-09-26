#!/usr/bin/env python3
"""collatz_thin_20260925_controls2.py -- further FINITE-EXACT controls for THM-4476.

(C1) Periodic points of T_b for several odd b: all cycles with minimum <= N are found
     by the standard descent search, and the union of their elements is counted
     against X^(h*) (Corollary 4: the union of all cycles is thin).
(C2) A contracting Conway/Matthews-Watts map with m = 3: g(x) = x/3, (2x+1)/3, (4x+1)/3
     on x = 0, 1, 2 mod 3 (prod p_i = 8 < 27, max p_i = 4 < 9). Its Chernoff rate I(g)
     and the direct no-dip count F_g(X, theta) for X = 3^t, compared with X^(1 - I_theta).
     Also the residue-word bijection modulo 3^k is checked exhaustively for k <= 7.
(C3) The same bijection check for T_(+-1) modulo 2^k, and its FAILURE for the sign
     strategy -chi_{-4} (the withdrawn corollary's witness): the number of distinct
     parity words among the 2^k residues is far below 2^k.
(C4) Corollary 9 illustration: the log-drift word d_j = floor(j log_2 3 - a log_2 j)
     with a = 1.03 (excluded by Corollary 9 since a < 1/h* = 1.0527). Every finite
     prefix is realised by some odd integer (Terras), and the least such integer is
     printed for each prefix length: it grows like 2^(d_k), i.e. finite realisability
     at every length with the least representative tending to infinity.
"""
import math, sys


def h(p):
    return 0.0 if p <= 0 or p >= 1 else -p * math.log2(p) - (1 - p) * math.log2(1 - p)


# ------------------------------------------------------------ (C1)
def T(x, b):
    return x // 2 if x % 2 == 0 else (3 * x + b) // 2


def cycles_below(N, b, cap=1 << 62, maxsteps=100000):
    """all cycles of T_b (on positive integers) whose minimum is <= N; returns list of sets."""
    cycles = []
    for n in range(1, N + 1):
        x = n
        steps = 0
        while True:
            if x < n:
                break
            if x == n and steps > 0:
                cyc = set()
                y = n
                while True:
                    cyc.add(y)
                    y = T(y, b)
                    if y == n:
                        break
                cycles.append(cyc)
                break
            x = T(x, b)
            steps += 1
            if x > cap or steps > maxsteps or x <= 0:
                break
    return cycles


def control_c1(N=200000):
    print("== (C1) union of all cycles of T_b with minimum <= %d ==" % N)
    hs = h(1 / math.log2(3))
    for b in (1, -1, 5, -5, 7, -7, 11, 13, -11, -13, 17, 23):
        cs = cycles_below(N, b)
        pts = set().union(*cs) if cs else set()
        mins = sorted(min(c) for c in cs)
        X = max(pts) if pts else 1
        print("   b=%+3d: %d cycles (minima %s), %d periodic points, max %d, count/X^h* = %.4f"
              % (b, len(cs), mins[:8], len(pts), X, len(pts) / X ** hs))


# ------------------------------------------------------------ (C2)
P3 = {0: (1, 0), 1: (2, 1), 2: (4, 1)}   # g(x) = (p x + q)/3 on x = i mod 3


def g3(x):
    p, q = P3[x % 3]
    return (p * x + q) // 3


def chernoff_rate(m, ps):
    best = 10.0
    for i in range(1, 20001):
        lam = i / 2000
        v = sum((p / m) ** lam for p in ps) / m
        best = min(best, v)
    return -math.log(best) / math.log(m)


def bijection_check(mapf, m, kmax):
    ok = True
    for k in range(1, kmax + 1):
        M = m ** k
        words = set()
        for r in range(M):
            x = r
            w = []
            for _ in range(k):
                w.append(x % m)
                x = mapf(x)
            words.add(tuple(w))
        if len(words) != M:
            ok = False
            print("      k=%d: %d distinct words among %d residues" % (k, len(words), M))
    return ok


def F_count_general(X, theta, mapf, m):
    k = int(math.log(X) / math.log(m))
    thr = X ** (-theta)
    cnt = 0
    for y in range(1, X + 1):
        x = y
        fl = y * thr
        ok = True
        for _ in range(k):
            x = mapf(x)
            if x < fl:
                ok = False
                break
        if ok:
            cnt += 1
    return cnt, k


def control_c2():
    print("== (C2) contracting Conway map m=3: g = x/3, (2x+1)/3, (4x+1)/3 ==")
    ps = [1, 2, 4]
    I0 = chernoff_rate(3, ps)
    print("   prod p_i = %d < 27, max p_i = 4 < 9; Chernoff rate I(g) = %.5f (base 3), i.e. exponent 1 - I = %.5f"
          % (math.prod(ps), I0, 1 - I0))
    print("   residue-word bijection mod 3^k, k <= 7: %s" % bijection_check(g3, 3, 7))
    for theta in (0.03, 0.1):
        best = 10.0
        for i in range(1, 20001):
            lam = i / 2000
            v = sum((p / 3) ** lam for p in ps) / 3 * 3 ** (theta * lam)
            best = min(best, v)
        I_th = -math.log(best) / math.log(3)
        row = []
        for t in (8, 10, 12):
            X = 3 ** t
            F, k = F_count_general(X, theta, g3, 3)
            row.append((t, F, math.log(F) / math.log(X)))
        print("   theta=%.2f: bound exponent 1 - I_theta = %.4f; counts " % (theta, 1 - I_th) +
              ", ".join("3^%d: F=%d (exp %.4f)" % r for r in row))


# ------------------------------------------------------------ (C3)
def control_c3():
    print("== (C3) parity-word bijection mod 2^k: constant shifts vs the strategy -chi_{-4} ==")
    for b in (1, -1):
        print("   T_%+d: bijection for k <= 12: %s" % (b, bijection_check(lambda x, b=b: T(x, b), 2, 12)))

    def Tneg(x):
        if x % 2 == 0:
            return x // 2
        s = -1 if x % 4 == 1 else 1
        return (3 * x + s) // 2
    print("   -chi_{-4} strategy (sigma(1 mod 4) = -1, sigma(3 mod 4) = +1): bijection for k <= 10: %s"
          % bijection_check(Tneg, 2, 10))


# ------------------------------------------------------------ (C4)
def log_drift_word(a, K):
    d = [0]
    for j in range(1, K + 1):
        nxt = 2 if j == 1 else math.floor(j * math.log2(3) - a * math.log2(j))
        d.append(max(d[-1] + 1, nxt))
    return [d[j] - d[j - 1] for j in range(1, K + 1)]


def first_k(n, k):
    out = []
    m = n
    for _ in range(k):
        x = 3 * m + 1
        v = (x & -x).bit_length() - 1
        out.append(v)
        m = x >> v
    return out


def control_c4(a=1.03, K=40, cap=1 << 24):
    print("== (C4) Corollary 9 illustration: least odd n realising each prefix of the log-drift word a=%.2f ==" % a)
    vs = log_drift_word(a, K)
    print("   halving counts v_1..v_20: %s" % vs[:20])
    n = 1
    for k in range(2, K + 1):
        target = vs[:k]
        while n < cap and first_k(n, k) != target:
            n += 2
        if n >= cap:
            print("   k=%d: no odd n < 2^24 has this prefix (d_k = %d)" % (k, sum(target)))
            break
        print("   k=%2d: least odd n = %d (log2 n = %.2f, d_k = %d)" % (k, n, math.log2(n), sum(target)))


if __name__ == '__main__':
    control_c1()
    control_c2()
    control_c3()
    control_c4()
