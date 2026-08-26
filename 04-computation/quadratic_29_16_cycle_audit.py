#!/usr/bin/env python3
"""Exact audit for f(x)=x^2-29/16 and two finite-ring controls.

The proof that the displayed rational graph is complete is mathematical:

* an odd prime in the denominator has its negative valuation doubled;
* at 2, every preperiodic rational must have valuation exactly -2, and an
  odd numerator stays at denominator 4 because a^2-29 == 4 (mod 8);
* if |x| > (2+sqrt(33))/4, then f(x) > |x| and the real orbit escapes.

Consequently a rational preperiodic point must be a/4 with
a in {-7,-5,-3,-1,1,3,5,7}.  This script audits the resulting graph, the
arithmetic-progression/Pythagorean parameter point, and the mod-63 and
mod-43^2 positive/hostile controls without floating point arithmetic.
"""

from collections import Counter
from fractions import Fraction
from math import gcd


C = Fraction(-29, 16)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def poly_trim(poly: list[int]) -> list[int]:
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def poly_add(*polys: list[int]) -> list[int]:
    result = [0] * max(map(len, polys))
    for poly in polys:
        for index, coefficient in enumerate(poly):
            result[index] += coefficient
    return poly_trim(result)


def poly_scale(poly: list[int], scalar: int) -> list[int]:
    return poly_trim([scalar * coefficient for coefficient in poly])


def poly_mul(left: list[int], right: list[int]) -> list[int]:
    result = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] += a * b
    return poly_trim(result)


def f(x: Fraction) -> Fraction:
    return x * x + C


def scaled_T(u: int, modulus: int | None = None) -> int:
    """T(u)=(u^2-29)/4, exactly on odd u or modulo an odd modulus."""
    if modulus is None:
        numerator = u * u - 29
        require(numerator % 4 == 0, "integer scaled_T called on an invalid parity")
        return numerator // 4
    return ((u * u - 29) * pow(4, -1, modulus)) % modulus


def canonical_cycle(cycle: list[int]) -> tuple[int, ...]:
    rotations = [tuple(cycle[i:] + cycle[:i]) for i in range(len(cycle))]
    return min(rotations)


def finite_cycles(modulus: int) -> list[tuple[int, ...]]:
    """All cycles of scaled_T on Z/modulus Z, once up to rotation."""
    cycles: set[tuple[int, ...]] = set()
    for start in range(modulus):
        path: list[int] = []
        first_at: dict[int, int] = {}
        u = start
        while u not in first_at:
            first_at[u] = len(path)
            path.append(u)
            u = scaled_T(u, modulus)
        cycles.add(canonical_cycle(path[first_at[u] :]))
    return sorted(cycles, key=lambda cycle: (len(cycle), cycle))


def rational_graph() -> None:
    candidates = [Fraction(a, 4) for a in (-7, -5, -3, -1, 1, 3, 5, 7)]
    image = {x: f(x) for x in candidates}
    require(set(image.values()) <= set(candidates), "candidate graph is not closed")

    cycle = [Fraction(-7, 4), Fraction(5, 4), Fraction(-1, 4)]
    require(
        [f(x) for x in cycle] == cycle[1:] + cycle[:1],
        "claimed rational 3-cycle fails",
    )

    print("COMPLETE RATIONAL PREPERIODIC GRAPH (proof filters stated above)")
    for x in candidates:
        print(f"  {x!s:>4} -> {image[x]}")
    print("  unique rational cycle: -7/4 -> 5/4 -> -1/4 -> -7/4")
    print("  rational exact period 6: NONE")


def moduli_and_triangle() -> None:
    """Audit the standard rational 3-cycle parameter at its AP values."""

    # Low-to-high coefficient arrays for the full rational-function chart.
    # These checks certify the identities for indeterminate t, not merely at
    # the three AP specializations used below.
    D_poly = [0, 2, 2]
    A0_poly = [1, 1, 2, 1]
    A1_poly = [-1, -1, 0, 1]
    A2_poly = [-1, -3, -2, -1]
    N_poly = [1, 4, 9, 8, 4, 2, 1]
    phi3_poly = [1, 1, 1]
    zero = [0]
    require(
        poly_add(poly_mul(A0_poly, A0_poly), poly_scale(N_poly, -1), poly_scale(poly_mul(A1_poly, D_poly), -1)) == zero,
        "p0^2+c=p1 chart identity fails",
    )
    require(
        poly_add(poly_mul(A1_poly, A1_poly), poly_scale(N_poly, -1), poly_scale(poly_mul(A2_poly, D_poly), -1)) == zero,
        "p1^2+c=p2 chart identity fails",
    )
    require(
        poly_add(poly_mul(A2_poly, A2_poly), poly_scale(N_poly, -1), poly_scale(poly_mul(A0_poly, D_poly), -1)) == zero,
        "p2^2+c=p0 chart identity fails",
    )
    require(
        poly_add(poly_scale(A1_poly, 2), poly_scale(A0_poly, -1), poly_scale(A2_poly, -1))
        == poly_scale(poly_mul([-1, 1], phi3_poly), 2),
        "p1-middle AP factor fails",
    )
    require(
        poly_add(poly_scale(A0_poly, 2), poly_scale(A1_poly, -1), poly_scale(A2_poly, -1))
        == poly_scale(poly_mul([2, 1], phi3_poly), 2),
        "p0-middle AP factor fails",
    )
    require(
        poly_add(poly_scale(A2_poly, 2), poly_scale(A0_poly, -1), poly_scale(A1_poly, -1))
        == poly_scale(poly_mul([1, 2], phi3_poly), -2),
        "p2-middle AP factor fails",
    )

    def parameter_cycle(t: Fraction) -> tuple[list[Fraction], Fraction]:
        denominator = 2 * t * (t + 1)
        points = [
            (t**3 + 2 * t**2 + t + 1) / denominator,
            (t**3 - t - 1) / denominator,
            -(t**3 + 2 * t**2 + 3 * t + 1) / denominator,
        ]
        parameter = -(
            t**6
            + 2 * t**5
            + 4 * t**4
            + 8 * t**3
            + 9 * t**2
            + 4 * t
            + 1
        ) / (4 * t**2 * (t + 1) ** 2)
        return points, parameter

    expected = {Fraction(-7, 4), Fraction(-1, 4), Fraction(5, 4)}
    print("\nAP LOCUS IN THE STANDARD RATIONAL 3-CYCLE PARAMETER")
    for t in (Fraction(1), Fraction(-2), Fraction(-1, 2)):
        points, parameter = parameter_cycle(t)
        local_f = lambda x: x * x + parameter
        require(
            [local_f(x) for x in points] == points[1:] + points[:1],
            f"marked 3-cycle formula fails at t={t}",
        )
        require(set(points) == expected, f"AP point set fails at t={t}")
        require(parameter == C, f"parameter c fails at t={t}")
        print(f"  t={t}: points={tuple(map(str, points))}, c={parameter}")

    # At the AP representative t=1, Euclid's pair (t+1,t)=(2,1)
    # gives the primitive 3-4-5 triangle.  The cycle denominator is its
    # even leg, while the c denominator is the square of that leg.
    t = 1
    m, n = t + 1, t
    odd_leg = m * m - n * n
    even_leg = 2 * m * n
    hypotenuse = m * m + n * n
    require(
        (odd_leg, even_leg, hypotenuse) == (3, 4, 5),
        "Euclid specialization does not give 3-4-5",
    )
    require(16 == even_leg**2, "16=4^2 identity fails")
    require(29 == hypotenuse**2 + m**2, "29=5^2+2^2 identity fails")
    require(7 == even_leg**2 - odd_leg**2, "7=4^2-3^2 identity fails")
    print("  t=1 -> Euclid pair (2,1) -> triangle (3,4,5)")
    print("  exact identities: 16=4^2, 29=5^2+2^2, 7=4^2-3^2")

    # Cyclic relabeling of a marked 3-cycle acts on its parameter by
    # R(t)=-1/(t+1).  Its SL_2 lift A has A^3=-I and A^6=I: projectivizing
    # loses the central sign and turns a genuine linear period 6 into period 3.
    def R(value: Fraction) -> Fraction:
        return -1 / (value + 1)

    orbit = [Fraction(1)]
    for _ in range(2):
        orbit.append(R(orbit[-1]))
    require(
        orbit == [Fraction(1), Fraction(-1, 2), Fraction(-2)],
        "AP parameter relabeling orbit fails",
    )
    require(R(orbit[-1]) == orbit[0], "projective relabeling is not period 3")

    def matmul(
        left: tuple[tuple[int, int], tuple[int, int]],
        right: tuple[tuple[int, int], tuple[int, int]],
    ) -> tuple[tuple[int, int], tuple[int, int]]:
        return (
            (
                left[0][0] * right[0][0] + left[0][1] * right[1][0],
                left[0][0] * right[0][1] + left[0][1] * right[1][1],
            ),
            (
                left[1][0] * right[0][0] + left[1][1] * right[1][0],
                left[1][0] * right[0][1] + left[1][1] * right[1][1],
            ),
        )

    identity = ((1, 0), (0, 1))
    A = ((0, -1), (1, 1))
    powers = [identity]
    for _ in range(6):
        powers.append(matmul(powers[-1], A))
    require(powers[3] == ((-1, 0), (0, -1)), "A^3 != -I")
    require(powers[6] == identity, "A^6 != I")
    print(f"  relabeling R(t)=-1/(t+1) has AP orbit {tuple(map(str, orbit))}")
    print("  SL2 lift [[0,-1],[1,1]]: A^3=-I, A^6=I, charpoly=lambda^2-lambda+1")

    # There is also a unique projective interpolant on the point coordinate
    # itself.  In x it is (4x-13)/(16x+12); after the AP normalization
    # r=4x+1 it is N(r)=2(r-6)/(r+2).  The determinant-one lifts have cube -I.
    direct_cycle = [Fraction(-7, 4), Fraction(5, 4), Fraction(-1, 4)]

    def projective_x(x: Fraction) -> Fraction:
        return (4 * x - 13) / (16 * x + 12)

    require(
        [projective_x(x) for x in direct_cycle]
        == direct_cycle[1:] + direct_cycle[:1],
        "direct projective interpolant fails",
    )
    direct_lift = ((Fraction(1, 4), Fraction(-13, 16)), (1, Fraction(3, 4)))

    def frac_matmul(left, right):
        return (
            (
                left[0][0] * right[0][0] + left[0][1] * right[1][0],
                left[0][0] * right[0][1] + left[0][1] * right[1][1],
            ),
            (
                left[1][0] * right[0][0] + left[1][1] * right[1][0],
                left[1][0] * right[0][1] + left[1][1] * right[1][1],
            ),
        )

    direct_square = frac_matmul(direct_lift, direct_lift)
    direct_cube = frac_matmul(direct_square, direct_lift)
    require(direct_cube == ((-1, 0), (0, -1)), "direct lift cube is not -I")
    print("  point interpolant M(x)=(4x-13)/(16x+12); det-one lift cube=-I")


def horizontal_zero_fibre_control() -> None:
    """Typed incidence with, and separation from, THM-4134's target section."""
    a = Fraction(-48)
    x_cycle = [Fraction(-7, 4), Fraction(5, 4), Fraction(-1, 4)]
    r_cycle = [4 * x + 1 for x in x_cycle]
    require(r_cycle == [-6, 6, 0], "AP normalization r=4x+1 changed")

    # The q=0 nodal fibre is V^2=(U-a)(U+a/2)^2.  Its normalization is
    # nu(r)=(r^2+a, r(r^2+3a/2)), with weights two and three in r.
    def nu(r: Fraction) -> tuple[Fraction, Fraction]:
        return r * r + a, r * (r * r + 3 * a / 2)

    def zero_fibre_residual(point: tuple[Fraction, Fraction]) -> Fraction:
        U, V = point
        return V * V - U**3 + 3 * a * a * U / 4 + a**3 / 4

    fibre_points = [nu(r) for r in r_cycle]
    require(
        fibre_points == [(-12, 216), (-12, -216), (-48, 0)],
        "a=-48 zero-fibre points changed",
    )
    require(
        all(zero_fibre_residual(point) == 0 for point in fibre_points),
        "normalized points miss the a=-48,q=0 fibre",
    )

    # Conjugating the projective cycle by r=4x+1 gives N(r).  It does not
    # descend to an automorphism of the nodal cubic: its fixed points satisfy
    # r^2=-12, whereas the two conductor preimages satisfy r^2=72.
    def N(r: Fraction) -> Fraction:
        return 2 * (r - 6) / (r + 2)

    require([N(r) for r in r_cycle] == r_cycle[1:] + r_cycle[:1], "N orbit fails")
    require(-12 != -3 * a / 2, "projective fixed pair unexpectedly equals node pair")

    # THM-4134's polynomial section is 3P on q=a^3/2+rho^2.  On q=0,
    # rho^2=-a^3/2 and its U-coordinate is 56/3, so it is not any nu(r_i).
    rho_squared = -a**3 / 2
    section_U = a / 2 + 16 * rho_squared / (9 * a * a)
    section_V_over_rho = -1 - 64 * rho_squared / (27 * a**3)
    require(section_U == Fraction(56, 3), "q=0 horizontal-section U changed")
    require(section_V_over_rho == Fraction(5, 27), "q=0 section V/rho changed")
    require(
        all(section_U != point[0] for point in fibre_points),
        "cycle-normalization point unexpectedly lies on the 3P section",
    )

    print("\nHORIZONTAL TARGET-PENCIL CONTROL")
    print("  r=4x+1 sends the cycle to (-6,6,0)")
    print(f"  a=-48,q=0 nodal normalization points: {fibre_points}")
    print("  projective N(r)=2(r-6)/(r+2) permutes parameters but does not preserve the node pair")
    print("  THM-4134 3P section at q=0: U=56/3, V=(5/27)rho; DISJOINT from these points")
    print("  verdict: zero-fibre incidence PROVED; horizontal-carrier identification REFUTED")


def finite_ring_controls() -> None:
    print("\nFINITE-RING CONTROLS FOR T(u)=(u^2-29)/4")
    for modulus in (7, 9, 63):
        cycles = finite_cycles(modulus)
        print(f"  mod {modulus}: cycles={cycles}")
    cycles_63 = finite_cycles(63)
    require(
        cycles_63 == [(5, 62, 56), (14, 26, 20), (35, 47, 41)],
        "mod-63 cycle census changed",
    )
    require(
        all(len(cycle) == 3 for cycle in cycles_63),
        "mod-63 hostile control acquired a non-3-cycle",
    )
    print("  hostile control: mod 63 has no exact 6-cycle")

    # The rational 3-cycle multiplier is
    # (-7/2)*(5/2)*(-1/2) = 35/8, hence -1 modulo the unique good prime 43.
    p = 43
    multiplier = Fraction(35, 8)
    require(
        (35 * pow(8, -1, p)) % p == p - 1,
        "3-cycle multiplier is not -1 mod 43",
    )
    modulus = p * p
    cycles = finite_cycles(modulus)
    counts = Counter(map(len, cycles))
    require(counts[3] == 1, "mod-43^2 census lost the rational 3-cycle")
    require(counts[6] == 21, "mod-43^2 census does not have 21 six-cycles")

    # The entire tube above {-7,5,-1} mod 43 contains 3*43 points.  It is
    # exactly the old 3-cycle plus the 21 six-cycles (3 + 21*6 = 129).
    base = {(-7) % p, 5 % p, (-1) % p}
    tube_cycles = [c for c in cycles if all(u % p in base for u in c)]
    require(
        Counter(map(len, tube_cycles)) == Counter({6: 21, 3: 1}),
        "mod-43^2 residue tube cycle census fails",
    )
    require(
        sum(map(len, tube_cycles)) == 3 * p,
        "mod-43^2 residue tube does not account for all 3*43 points",
    )

    representative = [36]
    for _ in range(5):
        representative.append(scaled_T(representative[-1], modulus))
    require(
        scaled_T(representative[-1], modulus) == representative[0],
        "mod-1849 representative does not close",
    )
    require(len(set(representative)) == 6, "mod-1849 representative is not exact")
    print(f"  3-cycle multiplier = {multiplier} == -1 (mod 43)")
    print(f"  mod 43^2 tube: {dict(sorted(Counter(map(len, tube_cycles)).items()))}")
    print(f"  representative exact 6-cycle mod 1849: {representative}")
    print("  63 local F=T^3 two-cycles = 3*(43-1)/2; T groups them into 21 six-cycles")


def multiplicative_order(base: int, modulus: int) -> int:
    require(gcd(base, modulus) == 1, "multiplicative order requires a unit")
    value = 1
    for order in range(1, modulus + 1):
        value = value * base % modulus
        if value == 1:
            return order
    raise RuntimeError("multiplicative-order search exceeded Euler bound")


def doubling_cycles_63() -> None:
    """Exact period-six roots for S(z)=z^2 through exponent doubling."""
    mersenne = 2**6 - 1
    require(mersenne == 63 == 3**2 * 7, "Mersenne factorization fails")
    old_prime_support = {prime for prime in (3, 7) if mersenne % prime == 0}
    require(old_prime_support == {3, 7}, "prime support of 63 changed")
    require((2**2 - 1) % 3 == 0, "3 is not inherited from exponent 2")
    require((2**3 - 1) % 7 == 0, "7 is not inherited from exponent 3")

    orders = {modulus: multiplicative_order(2, modulus) for modulus in (9, 21, 63)}
    require(orders == {9: 6, 21: 6, 63: 6}, "conductor orders are not all six")

    def doubling_cycles(modulus: int) -> list[tuple[int, ...]]:
        cycles: set[tuple[int, ...]] = set()
        for exponent in range(modulus):
            path: list[int] = []
            at: dict[int, int] = {}
            value = exponent
            while value not in at:
                at[value] = len(path)
                path.append(value)
                value = 2 * value % modulus
            cycles.add(canonical_cycle(path[at[value] :]))
        return sorted(cycles, key=lambda cycle: (len(cycle), cycle))

    cycles = doubling_cycles(63)
    census = Counter(map(len, cycles))
    require(
        census == Counter({6: 9, 3: 2, 2: 1, 1: 1}),
        "doubling-cycle census modulo 63 fails",
    )
    exact_six = [cycle for cycle in cycles if len(cycle) == 6]
    require(sum(map(len, exact_six)) == 54, "exact-six degree is not 54")
    require(
        6 + 12 + 36 == 54,
        "cyclotomic degrees phi(9)+phi(21)+phi(63) do not sum to 54",
    )

    print("\nSQUARING MAP S(z)=z^2 ON 63RD ROOTS / EXPONENT DOUBLING")
    print("  2^6-1=63=3^2*7; 3 already divides 2^2-1 and 7 already divides 2^3-1")
    print("  primitive prime divisors at exponent 6: NONE")
    print("  exact-period-6 factor: ((z^63-1)(z-1))/((z^7-1)(z^3-1))")
    print("                       = cyclotomic Phi_9 Phi_21 Phi_63")
    print(f"  conductor orders: {orders}; degrees: 6+12+36=54")
    print(f"  exponent-doubling cycle census mod 63: {dict(sorted(census.items()))}")
    for cycle in exact_six:
        print(f"    six-cycle {cycle}")
    print("  mechanism split: conductor 9 uses 3-adic depth; conductor 21 uses CRT synchronization;")
    print("                   conductor 63 combines the two. No single mechanism is universal.")


if __name__ == "__main__":
    rational_graph()
    moduli_and_triangle()
    horizontal_zero_fibre_control()
    finite_ring_controls()
    doubling_cycles_63()
