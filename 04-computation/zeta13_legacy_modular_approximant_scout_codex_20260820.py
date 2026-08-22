"""Exact coefficient scout for Calegari's legacy zeta_p(3) construction.

Status: FINITE-EXACT / method-specific, not an irrationality proof.

For the genus-zero primes p in {2,3,5,7,13}, Calegari expands

    E_2^*(E'_{-2} + eta) = sum_n (a_n + eta b_n) f^n,

where f=(Delta(p*tau)/Delta(tau))^(1/(p-1)).  Signs and the harmless
factor 2 depend on convention.  This script constructs the coefficient
rows exactly over Q, checks the printed p=2 coefficients from the paper,
and audits the first p=13 rational approximants for denominator/gcd savings.

The calculation deliberately separates three quantities which are easy to
conflate: p-adic analytic gain, Archimedean coefficient growth, and rational
denominator growth.  A finite prefix can suggest a gcd sidecar, but cannot
replace the asymptotic estimates or nonvanishing argument.
"""

from __future__ import annotations

from fractions import Fraction
import math


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def mul(a: list[Fraction], b: list[Fraction], degree: int) -> list[Fraction]:
    out = [Fraction(0) for _ in range(degree + 1)]
    for i, ai in enumerate(a):
        if not ai:
            continue
        stop = min(len(b) - 1, degree - i)
        for j in range(stop + 1):
            if b[j]:
                out[i + j] += ai * b[j]
    return out


def sigma_power_without_p(n: int, power: int, p: int) -> Fraction:
    ans = Fraction(0)
    for d in range(1, math.isqrt(n) + 1):
        if n % d:
            continue
        e = n // d
        if d % p:
            ans += Fraction(d**power) if power >= 0 else Fraction(1, d ** (-power))
        if e != d and e % p:
            ans += Fraction(e**power) if power >= 0 else Fraction(1, e ** (-power))
    return ans


def uniformizer(p: int, degree: int) -> list[Fraction]:
    """q * product_n (1+q^n+...+q^((p-1)n))^(24/(p-1))."""

    require(24 % (p - 1) == 0, "eta quotient exponent is not integral")
    exponent = 24 // (p - 1)
    product = [Fraction(1)] + [Fraction(0)] * degree
    for n in range(1, degree + 1):
        factor = [Fraction(0)] * (degree + 1)
        for j in range(p):
            if j * n <= degree:
                factor[j * n] = Fraction(1)
        for _ in range(exponent):
            product = mul(product, factor, degree)
    return [Fraction(0)] + product[:degree]


def decompose_in_f(target: list[Fraction], f: list[Fraction]) -> list[Fraction]:
    """Return c with target(q)=sum c_n f(q)^n modulo q^(N+1)."""

    degree = len(target) - 1
    require(f[0] == 0 and f[1] == 1, "f must be a normalized uniformizer")
    powers = [[Fraction(1)] + [Fraction(0)] * degree]
    for _ in range(degree):
        powers.append(mul(powers[-1], f, degree))
    residual = target[:]
    coeffs: list[Fraction] = []
    for n in range(degree + 1):
        cn = residual[n]
        coeffs.append(cn)
        if cn:
            for j in range(n, degree + 1):
                residual[j] -= cn * powers[n][j]
    require(all(x == 0 for x in residual), "triangular f-decomposition failed")
    return coeffs


def coefficient_rows(p: int, degree: int) -> tuple[list[Fraction], list[Fraction]]:
    f = uniformizer(p, degree)
    e2 = [Fraction(p - 1, 24)] + [
        sigma_power_without_p(n, 1, p) for n in range(1, degree + 1)
    ]
    e_minus_prime = [Fraction(0)] + [
        sigma_power_without_p(n, -3, p) for n in range(1, degree + 1)
    ]
    product = mul(e2, e_minus_prime, degree)
    b = decompose_in_f(e2, f)
    a = decompose_in_f(product, f)
    return a, b


def valuation(n: int, p: int) -> int:
    require(n != 0, "valuation of zero requested")
    n = abs(n)
    ans = 0
    while n % p == 0:
        n //= p
        ans += 1
    return ans


def legacy_exponents(p: int) -> tuple[float, float, float]:
    gain = 12 / (p - 1)
    height = 6 / (p - 1) + 3 / math.log(p)
    return gain, height, gain / height


def main() -> None:
    a2, b2 = coefficient_rows(2, 5)
    printed_b2 = [1, 24, -552, 19392, -810024, 37210944]
    printed_a2 = [0, 1, 1, Fraction(-8072, 27), Fraction(160841, 9), Fraction(-1088512616, 1125)]
    # The PDF line-break renders the last denominator poorly; the first five
    # terms are an unambiguous independent normalization check.
    require([24 * x for x in b2] == printed_b2, "p=2 E2* control failed")
    require([24 * x for x in a2[:5]] == printed_a2[:5], "p=2 product control failed")
    print("p=2 printed coefficient controls PASS")

    p = 13
    degree = 80
    raw_a, raw_b = coefficient_rows(p, degree)
    # E2* has constant (p-1)/24=1/2 at p=13.  Clear that harmless
    # normalization so the displayed finite denominators are fully reduced.
    a = [2 * x for x in raw_a]
    b = [2 * x for x in raw_b]
    gain, height, theta = legacy_exponents(p)
    print(f"p={p} analytic gain exponent A={gain:.12f}")
    print(f"p={p} legacy height exponent B={height:.12f}")
    print(f"p={p} legacy irrationality exponent theta=A/B={theta:.12f}")
    print(f"gap B-A={height-gain:.12f}; required height reduction={(height-gain)/height:.4%}")
    print()
    print("finite exact rational rows (eta approximant is a_n/b_n up to sign)")
    print(" n  v13(b_n)  log(den(a_n))/n  log(height(a_n/b_n))/(n log 13)")
    selected = {5, 10, 15, 20, 30, 40, 50, 60, 70, 80}
    for n in range(1, degree + 1):
        require(b[n].denominator == 1, f"b_{n} unexpectedly nonintegral")
        if n not in selected:
            continue
        approximant = a[n] / b[n]
        denominator_rate = math.log(a[n].denominator) / n
        height_integer = max(abs(approximant.numerator), approximant.denominator)
        height_rate = math.log(height_integer) / (n * math.log(p))
        print(
            f"{n:2d} {valuation(b[n].numerator, p):9d} "
            f"{denominator_rate:18.12f} {height_rate:34.12f}"
        )

    # A denominator-only rescue would need an asymptotic lcm exponent below
    # (1/2)log(13), instead of the universal exponent 3 used in the proof.
    denominator_gate = 0.5 * math.log(p)
    print()
    print(f"denominator-only gate d < 0.5*log(13) = {denominator_gate:.12f}")
    print("universal proof uses d=3; the finite rows above are a scout, not a limit")


if __name__ == "__main__":
    main()
