"""Exact gcd audit for the Apéry-like zeta_2(5) recurrence.

Status: FINITE-EXACT hostile audit, not an asymptotic theorem.

Lai--Sprang--Zudilin obtain integer linear forms after multiplying their two
recurrence solutions by d_n^6/Phi_n.  A common divisor g_n of the resulting
integer coefficients can be divided out.  If one could prove

    liminf log(g_n)/n >= sigma,

then Bel's lemma would improve the published height exponent by sigma.  This
script reverse-engineers the saving required by the proposed 13.806686 bound
and checks the exact finite gcds through n=1000.  A finite value of g_n is not
an admissible replacement for the liminf.  Moreover, a 2-primary divisor must
be charged against the same 2-adic smallness that drives the approximation;
only its prime-to-2 part can be treated as a free Archimedean height saving.
"""

from __future__ import annotations

from fractions import Fraction
from math import gcd, isqrt, log

from sympy import factorint, ilcm, primerange


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def phi_product(n: int) -> int:
    """Phi_n = product_{max(sqrt(2n),3)<p<=n} p."""

    ans = 1
    for prime0 in primerange(max(isqrt(2 * n) + 1, 4), n + 1):
        prime = int(prime0)
        if prime * prime > 2 * n:
            ans *= prime
    return ans


def measure_proxy(saving: float) -> float:
    alpha = 16 * log(2)
    beta = 8 * log(2) + 5 - saving
    return alpha / (alpha - beta)


def measure_proxy_with_primary_loss(height_saving: float, primary_loss: float) -> float:
    """Finite proxy after charging a 2-primary divisor to 2-adic decay."""

    alpha = 16 * log(2) - primary_loss
    beta = 8 * log(2) + 5 - height_saving
    return alpha / (alpha - beta)


def main() -> None:
    proposed_mu = 13.806686
    alpha = 16 * log(2)
    published_beta = 8 * log(2) + 5
    required_saving = published_beta - alpha * (1 - 1 / proposed_mu)
    require(abs(required_saving - 0.2580822869982988) < 1e-14, "candidate inversion drift")
    print(f"candidate mu={proposed_mu:.6f} requires liminf log(g_n)/n >= {required_saving:.12f}")
    print()

    top = 1000
    selected = {100, 200, 300, 320, 350, 400, 500, 600, 700, 800, 900, 1000}
    r0_previous, r0 = Fraction(0), Fraction(-1024)
    r3_previous, r3 = Fraction(768), Fraction(73728)
    lcm = 1
    closest: tuple[float, int, float, float] | None = None
    row320_gcd = 0
    print(" n   log(g_n)/n   finite-mu proxy")
    for n in range(1, top + 1):
        lcm = int(ilcm(lcm, n))
        scale = Fraction(lcm**6, phi_product(n))
        integer0 = scale * r0
        integer3 = scale * r3
        require(integer0.denominator == 1, f"rho_(n,0) integrality failed at n={n}")
        require(integer3.denominator == 1, f"rho_(n,3) integrality failed at n={n}")
        common = gcd(abs(integer0.numerator), abs(integer3.numerator))
        saving = log(common) / n
        proxy = measure_proxy(saving)
        error = abs(saving - required_saving)
        if closest is None or error < closest[0]:
            closest = (error, n, saving, proxy)
        if n == 320:
            row320_gcd = common
        if n in selected:
            print(f"{n:4d} {saving:14.12f} {proxy:17.9f}")

        if n == top:
            continue
        coefficient = 32 * (2 * n + 1) * (
            8 * n**4 + 16 * n**3 + 20 * n**2 + 12 * n + 3
        )
        next0 = (
            coefficient * r0 - 2**16 * n**5 * r0_previous
        ) / (n + 1) ** 5
        next3 = (
            coefficient * r3 - 2**16 * n**5 * r3_previous
        ) / (n + 1) ** 5
        determinant = r0 * next3 - next0 * r3
        expected = Fraction(3 * 2 ** (16 * n + 18), (n + 1) ** 5)
        require(determinant == expected, f"Casoratian identity failed at n={n}")
        r0_previous, r0 = r0, next0
        r3_previous, r3 = r3, next3

    require(closest is not None, "empty scan")
    _, close_n, close_saving, close_proxy = closest
    print()
    print(
        f"closest finite match through {top}: n={close_n}, saving={close_saving:.12f}, "
        f"proxy={close_proxy:.9f}"
    )
    row320_factors = factorint(row320_gcd)
    require(row320_factors.get(2) == 63, "n=320 two-primary exponent drift")
    full_rate = log(row320_gcd) / 320
    two_primary_rate = 63 * log(2) / 320
    odd_rate = full_rate - two_primary_rate
    naive_proxy = measure_proxy(full_rate)
    charged_proxy = measure_proxy_with_primary_loss(full_rate, two_primary_rate)
    odd_only_proxy = measure_proxy(odd_rate)
    require(abs(charged_proxy - 16.42673965953) < 1e-10, "charged proxy drift")
    require(abs(odd_only_proxy - 16.63138363788) < 1e-10, "odd-only proxy drift")
    print("n=320 exact gcd factorization", row320_factors)
    print(f"n=320 two-primary loss rate = {two_primary_rate:.12f}")
    print(f"n=320 odd-gcd saving rate = {odd_rate:.12f}")
    print(f"naive full-gcd finite proxy = {naive_proxy:.9f}")
    print(f"charge 2^63 against 2-adic decay -> finite proxy = {charged_proxy:.9f}")
    print(f"divide only by odd gcd -> finite proxy = {odd_only_proxy:.9f}")
    print(
        "HOSTILE: the rate falls to",
        f"{log(gcd(abs(integer0.numerator), abs(integer3.numerator))) / top:.12f}",
        "at n=1000; no positive liminf is established",
    )
    print("all integrality and exact Casoratian gates PASS through n=1000")


if __name__ == "__main__":
    main()
