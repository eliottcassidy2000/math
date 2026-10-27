#!/usr/bin/env python3
"""
fast_set_exponent_model_20261027.py  (opus-2026-10-27-S64)  -- HEURISTIC model, no theorem.

The independent-3-adic (Lagarias-Weiss / Applegate-Lagarias tree-count) model of the fast set
F_C = { n : tau(n) <= C ln n }, tau the ordinary total stopping time.

In the T-convention (one-division map) an integer n in [X/2, X) that reaches 1 after k T-steps of which d are odd
satisfies 3^d n + c_w = 2^k with 0 < c_w < 2^k, so k ln 2 - d ln 3 = ln X + O(1); its ordinary stopping time is
tau = k + d. Among the C(k, d) parity words of this shape, the model takes the integrality condition 3^d | 2^k - c_w
to hold with probability 3^(-d), independently. With kappa = k / ln X and delta = d / ln X this gives

    #{n in [X/2, X) : tau(n) <= C ln n} = X^(E(C) + o(1)),
    E(C) = max { kappa h(delta/kappa) - delta ln 3 : delta = (kappa ln 2 - 1)/ln 3 >= 0, kappa + delta <= C },

h the entropy in nats. The unconstrained maximum is attained at delta/kappa = 1/2 (one odd step in two), i.e. at
kappa = 1/(ln 2 - (ln 3)/2) = 6.952 and tau/ln X = (3/2) kappa = 3/ln(4/3) = 10.428 = C*, where E = 1 exactly: the
model puts the whole population at the typical-word constant, so E(C) < 1 (density zero) for C < C* and E(C) = 1 for
C >= C*. Proposition F of the note proves O(X^(1 - c)) unconditionally only for C < 3.307. Output: a table of E(C)
and the optimising word shape; compared in the note with the measured two-block exponents of the 2^28 run.
"""
import math

LN2, LN3 = math.log(2.0), math.log(3.0)
CSTAR = 3.0 / math.log(4.0 / 3.0)


def h(p):
    return 0.0 if p <= 0.0 or p >= 1.0 else -(p * math.log(p) + (1.0 - p) * math.log(1.0 - p))


def exponent(C, N=200000):
    kmin = 1.0 / LN2
    kmax = (C + 1.0 / LN3) / (1.0 + LN2 / LN3)
    if kmax < kmin:
        return 0.0, kmin, 0.0
    best = (-1.0, kmin, 0.0)
    for i in range(N + 1):
        kappa = kmin + (kmax - kmin) * i / N
        delta = (kappa * LN2 - 1.0) / LN3
        val = kappa * h(delta / kappa) - delta * LN3
        if val > best[0]:
            best = (val, kappa, delta)
    return best


def main():
    print("C* = 3/ln(4/3) = %.6f; typical word kappa = %.4f, delta = %.4f (odd share 1/2)" % (CSTAR, 1.0 / (LN2 - LN3 / 2), 0.5 / (LN2 - LN3 / 2)))
    print("     C     E(C)    kappa   delta  odd share   tau/lnX at the optimum")
    for C in [1.5, 2.0, 2.5, 3.0, 3.307, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, CSTAR, 11.0, 12.0]:
        e, k, d = exponent(C)
        print("%8.4f  %.4f  %7.3f  %6.3f   %.3f      %.3f" % (C, e, k, d, d / k, k + d))
    e, k, d = exponent(100.0)
    print("unconstrained: E = %.6f at kappa = %.4f, delta = %.4f, tau/lnX = %.4f" % (e, k, d, k + d))


if __name__ == "__main__":
    main()
