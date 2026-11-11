#!/usr/bin/env python3
"""S91 (THM-5128 (D)): the mass constants of the F1 family and the limiting sheet shares.

With rho = x/y uniform on [1/2, 1] and first state (x, y), the number of second states with N <= X on the padding s is
(X / x) ln((rho + 2^s) / (rho/2 + 2^s)) asymptotically.  Telescoping over s:
        sum_{s >= 0} ln((rho + 2^s)/(rho/2 + 2^s)) = ln(1 + rho),          sum_{s >= 1} = ln(1 + rho/2),
        sum_{s >= 0} (s + 1) ln((rho + 2^s)/(rho/2 + 2^s)) = sum_{s >= 0} ln(1 + rho 2^-s)      (word pairs, weight s + 1).
Averages over rho (density 2 on [1/2, 1]) of these divided by rho are dilogarithm values:
        triples, all sheets : 2 (pi^2/12 + Li2(-1/2))                = g_bar
        triples, sheet s = 0: 2 (pi^2/12 + 2 Li2(-1/2) - Li2(-1/4))
        triples, sheet s >= 1: 2 (Li2(-1/4) - Li2(-1/2))
        word pairs, all     : pi^2/6
Checks: the telescoping identities numerically; the closed forms against numerical integration; and the first
coefficient against the measured means of r on the three windows (slope of mean r against mean ln N).
"""
import math
import sys


def li2(z, terms=200000):
    """dilogarithm for |z| <= 1 by its series (enough digits for |z| <= 1/2; for z = -1 use the closed form)."""
    if z == -1:
        return -math.pi ** 2 / 12
    s, zk = 0.0, 1.0
    for k in range(1, 400):
        zk *= z
        s += zk / (k * k)
    return s


def integrate(f, a, b, n=200000):
    h = (b - a) / n
    return h * sum(f(a + (i + 0.5) * h) for i in range(n))


def main():
    ok = True
    for rho in (0.5, 0.6180339887, 0.75, 1.0):
        tot = sum(math.log((rho + 2.0 ** s) / (rho / 2 + 2.0 ** s)) for s in range(0, 200))
        tot1 = sum(math.log((rho + 2.0 ** s) / (rho / 2 + 2.0 ** s)) for s in range(1, 200))
        totp = sum((s + 1) * math.log((rho + 2.0 ** s) / (rho / 2 + 2.0 ** s)) for s in range(0, 200))
        ref = sum(math.log1p(rho * 2.0 ** (-s)) for s in range(0, 200))
        ok &= abs(tot - math.log(1 + rho)) < 1e-12 and abs(tot1 - math.log(1 + rho / 2)) < 1e-12 and abs(totp - ref) < 1e-9
    print("telescoping identities:", ok)
    L2, L4 = li2(-0.5), li2(-0.25)
    g_all = 2 * (math.pi ** 2 / 12 + L2)
    g0 = 2 * (math.pi ** 2 / 12 + 2 * L2 - L4)
    g1 = 2 * (L4 - L2)
    num_all = integrate(lambda r: 2 * math.log(1 + r) / r, 0.5, 1.0)
    num0 = integrate(lambda r: 2 * math.log(2 * (1 + r) / (2 + r)) / r, 0.5, 1.0)
    num1 = integrate(lambda r: 2 * math.log(1 + r / 2) / r, 0.5, 1.0)
    nump = integrate(lambda r: 2 * sum(math.log1p(r * 2.0 ** (-s)) for s in range(0, 60)) / r, 0.5, 1.0, n=20000)
    ok &= abs(g_all - num_all) < 1e-9 and abs(g0 - num0) < 1e-9 and abs(g1 - num1) < 1e-9 and abs(nump - math.pi ** 2 / 6) < 1e-7
    ok &= abs(g0 + g1 - g_all) < 1e-12
    print(f"Li2(-1/2) = {L2:.12f}   Li2(-1/4) = {L4:.12f}")
    print(f"g_bar = pi^2/6 + 2 Li2(-1/2) = {g_all:.10f}   (mean of r over odd N is g_bar ln N + O(1))")
    print(f"sheet s = 0 : {g0:.10f}   share w_0 = {g0 / g_all:.8f}")
    print(f"sheet s >= 1: {g1:.10f}   share w_1 = {g1 / g_all:.8f}")
    print(f"word pairs (weight s + 1): total pi^2/6 = {math.pi ** 2 / 6:.10f}; sheet s = 0 share = {g0 / (math.pi ** 2 / 6):.8f}, "
          f"sheet s >= 1 share = {1 - g0 / (math.pi ** 2 / 6):.8f}   (letter-count model of THM-4916: 1/4, 3/4)")
    # measured means of r (window_split outputs): (mean ln N of the window, measured triples per odd N)
    def mean_ln(lo, hi):
        a, b = 2.0 ** lo, 2.0 ** hi
        return (b * math.log(b) - a * math.log(a)) / (b - a) - 1
    data = [(mean_ln(18, 21), 11.0291), (mean_ln(21, 23), 12.1898), (mean_ln(24, 27), 14.1405)]
    slopes = [(data[i + 1][1] - data[i][1]) / (data[i + 1][0] - data[i][0]) for i in range(2)]
    print("measured slope of mean r against mean ln N between consecutive windows:", [round(s, 4) for s in slopes],
          " (g_bar = %.4f)" % g_all)
    ok &= all(abs(s - g_all) < 0.01 for s in slopes)
    print("ALL CHECKS PASSED" if ok else "CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
