#!/usr/bin/env python3
"""collatz_thin_20260925_counts.py -- FINITE-EXACT controls for the thin-orbit theorem
(session collatz-squares-doubles-20260925, opus; THM-4476 candidate).

Theorem (thin divergence): every T_b-orbit (b = +-1) of a positive integer that is
not eventually periodic has  N(X) := #{j : m_j <= X} = O_eps(X^(h* + eps)),
h* = h(log_3 2) = 0.94995..., h the binary entropy.  Hence sum 1/m_j < infinity
(HYP-9160), the growth constant c = lim m_L 2^(d_L)/3^L is finite and positive, and
on the minus sheet R(d) < n strictly.

The proof has two counting ingredients, both checked here:
  (L1) class count: #{words of length k with o_k >= rho k - c} <= poly(k) 2^(k h(rho)),
       rho = (1-theta)/log_2 3 > 1/2  [exact binomial sums vs the entropy bound];
  (L2) integer count: F(X,theta) = #{m <= X : T^i(m) >= m X^(-theta) for all i <= floor(log_2 X)}
       is O(X^(h(rho)) log^2 X)  [direct enumeration for X = 2^t, t <= 20, both sheets].
It also prints the constants h*, 1/h*, the corollary thresholds, and the optimal theta
of the un-bootstrapped version (1 - theta = h(rho(theta))).
"""
import math, sys

ALPHA = math.log2(3.0)
RHO0 = 1.0 / ALPHA                     # log_3 2
THETA0 = 1.0 - ALPHA / 2.0             # 0.20752: mean drift per T-step (bits)


def h(p):
    if p <= 0 or p >= 1:
        return 0.0
    return -p * math.log2(p) - (1 - p) * math.log2(1 - p)


def rho(theta):
    return (1.0 - theta) / ALPHA


def constants():
    hs = h(RHO0)
    print("== constants ==")
    print("   alpha = log_2 3 = %.6f, rho0 = log_3 2 = %.6f, h* = h(rho0) = %.6f, 1/h* = %.6f" % (ALPHA, RHO0, hs, 1 / hs))
    print("   theta0 = 1 - alpha/2 = %.6f (the largest theta with rho(theta) > 1/2)" % THETA0)
    # thresholds of the corollaries
    print("   Corollary: no non-periodic positive orbit has m_j <= K j^a for all j with a < 1/h* = %.5f" % (1 / hs))
    print("   plus sheet log-band |Delta_j| <= C log_2 j excluded for C < (1/h* - 1)/2 = %.5f" % ((1 / hs - 1) / 2))
    print("   minus sheet log-band |delta_j| <= C log_2 j excluded for C < 1/h* = %.5f" % (1 / hs))
    # un-bootstrapped optimum: 1 - theta = h(rho(theta))
    lo, hi = 0.0, THETA0
    for _ in range(100):
        mid = (lo + hi) / 2
        if 1 - mid > h(rho(mid)):
            lo = mid
        else:
            hi = mid
    print("   un-bootstrapped exponent: theta* = %.5f, exponent = %.5f" % (lo, 1 - lo))
    # bootstrapped: exponent -> h* as theta -> 0; show h(rho(theta)) for small theta
    for th in (0.1, 0.05, 0.02, 0.01, 0.001):
        print("   theta = %.3f: rho = %.5f, h(rho) = %.5f" % (th, rho(th), h(rho(th))))


def class_count_check():
    print("== (L1) class counts: exact sum_{o >= rho k - c} C(k,o) versus 2^(k h(rho)) ==")
    for theta in (0.0, 0.03, 0.1):
        r = rho(theta)
        for k in (20, 40, 80, 160, 320):
            c = 2
            lo = max(0, math.ceil(r * k - c))
            s = sum(math.comb(k, o) for o in range(lo, k + 1))
            ratio = s / 2 ** (k * h(r))
            print("   theta=%.2f k=%3d: count = 2^%.3f, entropy bound 2^%.3f, ratio/(k+1) = %.4f"
                  % (theta, k, math.log2(s), k * h(r), ratio / (k + 1)))


def F_count(X, theta, b):
    """direct count of m <= X (odd and even) with T_b^i(m) >= m X^{-theta} for all i <= k."""
    k = X.bit_length() - 1
    thr = X ** (-theta)
    cnt = 0
    for m in range(1, X + 1):
        x = m
        ok = True
        floor_v = m * thr
        for _ in range(k):
            x = x // 2 if x % 2 == 0 else (3 * x + b) // 2
            if x < floor_v:
                ok = False
                break
        if ok:
            cnt += 1
    return cnt


def integer_count_check(tmax=20):
    print("== (L2) integer counts F(X,theta) for X = 2^t, both sheets, versus X^(h(rho)) ==")
    for theta in (0.03, 0.1):
        r = rho(theta)
        for b in (1, -1):
            row = []
            for t in range(10, tmax + 1, 2):
                X = 1 << t
                F = F_count(X, theta, b)
                row.append((t, F, math.log(F) / math.log(X)))
            print("   theta=%.2f sheet 3n%+d (h(rho)=%.4f): " % (theta, b, h(r)) +
                  ", ".join("2^%d: F=%d (exp %.4f)" % (t, F, e) for t, F, e in row))


if __name__ == '__main__':
    constants()
    class_count_check()
    tmax = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    integer_count_check(tmax)
