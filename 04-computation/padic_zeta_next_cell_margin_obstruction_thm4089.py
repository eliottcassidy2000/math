#!/usr/bin/env python3
"""Independent exact-interval audit for THM-4089.

This script imports neither the external next-case verifier nor its pinned
certificate.  It uses exact Fraction interval arithmetic, independent
alternating-series enclosures, and exact integer-square-root bounds.
"""

from __future__ import annotations

from decimal import Decimal, localcontext
from fractions import Fraction
from math import isqrt


Interval = tuple[Fraction, Fraction]
TOL = Fraction(1, 10**32)
SQRT_SCALE = 10**42


CASES = {
    (2, 31): (Fraction(5_137_457, 250_000_000), Fraction(20_549_829, 1_000_000_000)),
    (3, 13): (Fraction(18_268_941, 500_000_000), Fraction(36_537_883, 1_000_000_000)),
    (5, 7): (Fraction(10_794_307, 250_000_000), Fraction(43_177_229, 1_000_000_000)),
    (7, 5): (Fraction(46_640_751, 1_000_000_000), Fraction(2_915_047, 62_500_000)),
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def iadd(a: Interval, b: Interval) -> Interval:
    return a[0] + b[0], a[1] + b[1]


def isub(a: Interval, b: Interval) -> Interval:
    return a[0] - b[1], a[1] - b[0]


def imul(a: Interval, b: Interval) -> Interval:
    products = (a[0] * b[0], a[0] * b[1], a[1] * b[0], a[1] * b[1])
    return min(products), max(products)


def iscale(a: Interval, q: Fraction | int) -> Interval:
    q = Fraction(q)
    if q >= 0:
        return q * a[0], q * a[1]
    return q * a[1], q * a[0]


def sqrt_bounds(q: Fraction) -> Interval:
    require(q >= 0, "negative square-root argument")
    target = q.numerator * SQRT_SCALE * SQRT_SCALE
    n = isqrt(target // q.denominator)
    while (n + 1) ** 2 * q.denominator <= target:
        n += 1
    while n * n * q.denominator > target:
        n -= 1
    lo = Fraction(n, SQRT_SCALE)
    if n * n * q.denominator == target:
        return lo, lo
    return lo, Fraction(n + 1, SQRT_SCALE)


def atan_small(q: Fraction) -> Interval:
    """Alternating-series enclosure for 0 <= q <= 1/2."""

    require(0 <= q <= Fraction(1, 2), "atan_small argument outside [0,1/2]")
    total = Fraction(0)
    k = 0
    while True:
        term = q ** (2 * k + 1) / (2 * k + 1)
        total = total + term if k % 2 == 0 else total - term
        nxt = q ** (2 * k + 3) / (2 * k + 3)
        if nxt <= TOL:
            if k % 2 == 0:
                return total - nxt, total
            return total, total + nxt
        k += 1


def pi_bounds() -> Interval:
    a = atan_small(Fraction(1, 5))
    b = atan_small(Fraction(1, 239))
    return isub(iscale(a, 16), iscale(b, 4))


PI = pi_bounds()


def atan_unit(q: Fraction) -> Interval:
    """Enclose atan(q), 0 <= q < 1, with independent range reduction."""

    require(0 <= q < 1, "atan_unit argument outside [0,1)")
    if q <= Fraction(1, 2):
        return atan_small(q)
    w = (1 - q) / (1 + q)
    aw = atan_small(w)
    return PI[0] / 4 - aw[1], PI[1] / 4 - aw[0]


def log_integer_bounds(n: int) -> Interval:
    """Positive atanh-series enclosure of log(n)."""

    z = Fraction(n - 1, n + 1)
    total = Fraction(0)
    k = 0
    while True:
        total += z ** (2 * k + 1) / (2 * k + 1)
        nxt = z ** (2 * k + 3) / (2 * k + 3)
        tail = nxt / (1 - z * z)
        if tail <= TOL:
            return 2 * total, 2 * (total + tail)
        k += 1


def acos_bounds(x: Fraction) -> Interval:
    require(0 < x < 1, "acos argument outside (0,1)")
    t = sqrt_bounds((1 - x) / (1 + x))
    lo = atan_unit(t[0])
    hi = atan_unit(t[1])
    return 2 * lo[0], 2 * hi[1]


def euler_phi(n: int) -> int:
    result = n
    d = 2
    x = n
    while d * d <= x:
        if x % d == 0:
            while x % d == 0:
                x //= d
            result -= result // d
        d += 1
    if x > 1:
        result -= result // x
    return result


def harmonic(k: int) -> Fraction:
    return sum((Fraction(1, j) for j in range(1, k + 1)), Fraction(0))


def xi_star(p: int, s: int) -> Fraction:
    c = Fraction(p + 1, 12)
    return Fraction(s * s + s - 1, 1) / (s * s + (s - 1) * c)


def prime_window(s: int, xi: Fraction) -> Fraction:
    k = s // xi
    delta = max(Fraction(0), s + 1 - (k + 1) * xi)
    return Fraction(2 * s + 1, 2) * harmonic(k) - k * xi + delta * delta / (2 * (k + 1))


def tau(p: int, s: int, xi: Fraction) -> Fraction:
    c = Fraction(p + 1, 12)
    require(c * xi < 1, "unsaturated Hasse formula used past its wall")
    j = xi - c * xi * xi / 2
    small = s * s * xi - (s - 1) * j
    return Fraction(2, (s + 1) ** 2) * (small + s * prime_window(s, xi))


def derivative(p: int, s: int, y: Fraction) -> Interval:
    total: Interval = (Fraction(0), Fraction(0))
    q = p
    while q * y < 1:
        root = sqrt_bounds(1 - q * q * y * y)
        total = iadd(total, iscale(root, Fraction(euler_phi(q), q)))
        q += p
    bracket = isub(iscale(total, 4), (Fraction(s), Fraction(s)))
    return iscale(imul(PI, bracket), 2)


def margin(p: int, s: int, xi: Fraction, y: Fraction) -> Interval:
    logp = log_integer_bounds(p)
    lam = isub(iscale(logp, Fraction(12, p - 1)), iscale(PI, 2 * y))

    collision: Interval = (Fraction(0), Fraction(0))
    q = p
    while q * y < 1:
        x = q * y
        f = isub(acos_bounds(x), iscale(sqrt_bounds(1 - x * x), x))
        collision = iadd(collision, iscale(imul(PI, f), Fraction(4 * euler_phi(q), q * q)))
        q += p

    return isub(isub(iscale(lam, s), (Fraction((s + 1) * tau(p, s, xi)),) * 2), collision)


def check_symbolic_seams(s: int) -> None:
    """Exact value/derivative checks at every floor and positive-part seam."""

    a = Fraction(2 * s + 1, 2)
    for k in range(1, s):
        x = Fraction(s, k)
        left_value = a * harmonic(k) - k * x
        right_value = a * harmonic(k - 1) - (k - 1) * x + Fraction(1, 2 * k)
        require(left_value == right_value, f"floor-boundary value seam failed at k={k}")
        require(-k == k * x - (s + k), f"floor-boundary derivative seam failed at k={k}")

    require(a - s == Fraction(1, 2), "K=0 value seam failed")
    for k in range(s):
        x = Fraction(s + 1, k + 1)
        derivative_active = (k + 1) * x - (s + k + 1)
        require(derivative_active == -k, f"positive-part derivative seam failed at k={k}")


def decimal(q: Fraction, digits: int = 18) -> str:
    with localcontext() as ctx:
        ctx.prec = digits + 15
        x = Decimal(q.numerator) / Decimal(q.denominator)
        return f"{x:.{digits}f}"


def main() -> None:
    print("INDEPENDENT THM-4089 GLOBAL-MARGIN REFEREE")
    print("exact Fraction intervals; no import from candidate/source certificate")
    print()
    print("p  s  xi*       d(L) lower          d(R) upper          tangent global upper")

    for (p, s), (left, right) in CASES.items():
        check_symbolic_seams(s)
        c = Fraction(p + 1, 12)
        xi = xi_star(p, s)
        require(1 < xi < Fraction(s + 1, s), "xi* left the first floor cell")
        require(c * xi < 1, "xi* crossed the Hasse saturation wall")
        require(s // xi == s - 1, "wrong floor at xi*")

        dl = derivative(p, s, left)
        dr = derivative(p, s, right)
        require(dl[0] > 0, "left derivative interval is not positive")
        require(dr[1] < 0, "right derivative interval is not negative")

        ml = margin(p, s, xi, left)
        correction = iscale(dl, right - left)
        upper = iadd(ml, correction)[1]
        require(upper < 0, "independent tangent upper bound is not negative")

        print(
            f"{p:1d} {s:2d}  {str(xi):9s}  {decimal(dl[0]):>20s}  "
            f"{decimal(dr[1]):>20s}  {decimal(upper):>22s}"
        )

    print()
    print("SYMBOLIC SEAMS: PASS")
    print("DERIVATIVE BRACKETS: PASS")
    print("INDEPENDENT GLOBAL NEGATIVITY: PASS")


if __name__ == "__main__":
    main()
