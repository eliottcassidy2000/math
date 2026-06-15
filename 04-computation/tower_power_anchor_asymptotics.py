"""
tower_power_anchor_asymptotics.py

Asymptotics for the real balancing anchor a_p(n) defined by

    sum_{j=0}^n (a+j)^p = sum_{j=1}^n (a+n+j)^p.

The natural midpoint coordinate is c = a + n. Writing u = n(n+1), the balance
equation becomes

    c^p - 2 * sum_{r odd} binom(p,r) c^(p-r) S_r(n) = 0,

where S_r(n) = sum_{j=1}^n j^r. Expanding around c = p*u gives

    a_p(n) = p*n^2 + (p-1)*n + alpha_p + beta_p/u + O(u^-2),

with

    alpha_p = (p-1)(p-2)/(12p),
    beta_p  = -(p-1)(p-2)(2p^2-4p-1)/(180 p^3).
"""

from __future__ import annotations

from decimal import Decimal, getcontext


def alpha(p: int) -> Decimal:
    return Decimal((p - 1) * (p - 2)) / Decimal(12 * p)


def beta(p: int) -> Decimal:
    num = -(p - 1) * (p - 2) * (2 * p * p - 4 * p - 1)
    den = 180 * p * p * p
    return Decimal(num) / Decimal(den)


def balance_difference(p: int, n: int, a: Decimal) -> Decimal:
    left = sum((a + Decimal(j)) ** p for j in range(n + 1))
    right = sum((a + Decimal(n + 1 + j)) ** p for j in range(n))
    return left - right


def approx_anchor(p: int, n: int) -> Decimal:
    u = Decimal(n * (n + 1))
    main = Decimal(p * n * n + (p - 1) * n)
    return main + alpha(p) + beta(p) / u


def real_anchor(p: int, n: int) -> Decimal:
    # The asymptotic location is already tight enough to bracket the root.
    lo = approx_anchor(p, n) - Decimal(1)
    hi = approx_anchor(p, n) + Decimal(1)
    while balance_difference(p, n, lo) > 0:
        lo -= 1
    while balance_difference(p, n, hi) < 0:
        hi += 1
    for _ in range(120):
        mid = (lo + hi) / 2
        if balance_difference(p, n, mid) >= 0:
            hi = mid
        else:
            lo = mid
    return hi


def main() -> None:
    getcontext().prec = 100
    print("Asymptotic balancing anchor for power towers")
    print()
    print("a_p(n) = p*n^2 + (p-1)*n + alpha_p + beta_p/(n(n+1)) + O(n^-4)")
    print()
    print("p   alpha_p                  beta_p")
    for p in range(1, 11):
        print(f"{p:2d}  {alpha(p):<24} {beta(p)}")
    print()
    print("verification: error * (n(n+1))^2 should stay bounded")
    print()
    for p in range(3, 10):
        print(f"p={p}")
        for n in (10, 20, 50, 100):
            root = real_anchor(p, n)
            approx = approx_anchor(p, n)
            u = Decimal(n * (n + 1))
            err = root - approx
            scaled = err * u * u
            print(
                f"  n={n:3d}  root={root}  approx={approx}  "
                f"err={err}  err*u^2={scaled}"
            )
        print()


if __name__ == "__main__":
    main()
