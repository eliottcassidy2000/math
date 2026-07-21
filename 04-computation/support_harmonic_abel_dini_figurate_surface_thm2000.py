#!/usr/bin/env python3
"""Exact/numerical referee for THM-2000.

The user's convention is literal: an integer sequence contributes a *subset*
of the harmonic numbers.  Repetitions therefore disappear.  This referee
checks the support/multiplicity collision tax, Abel--Stieltjes summation,
multiplicative-block bounds, the Bertrand counterexamples to a false linear
growth dichotomy, and the two-axis reciprocal surface of the master figurate
array N(s,d,m).  It also inventories several sequence families already native
to the repository.  Analytic limit passages remain paper providers.
"""

from __future__ import annotations

import ast
from collections import Counter
from fractions import Fraction as F
from math import comb, factorial, log
from pathlib import Path

import mpmath as mp


mp.mp.dps = 50


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(f"THM-2000 referee failed: {message}")


def optimization_safety_probe() -> int:
    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    count = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
    require(count == 0, "Python assert nodes are optimization-sensitive")
    return count


def finite_signatures(values: list[int]) -> tuple[F, F, F]:
    counts = Counter(value for value in values if value > 0)
    support = sum((F(1, value) for value in counts), F(0))
    multiple = sum((F(count, value) for value, count in counts.items()), F(0))
    collision = sum((F(count - 1, value) for value, count in counts.items()), F(0))
    require(multiple == support + collision, "collision-tax decomposition")
    return support, multiple, collision


def collision_tax_audit() -> list[tuple[str, F, F, F]]:
    factorial = [1]
    for n in range(1, 13):
        factorial.append(factorial[-1] * n)

    fibonacci = [1, 1]
    for _ in range(28):
        fibonacci.append(fibonacci[-1] + fibonacci[-2])

    catalan = [comb(2 * n, n) // (n + 1) for n in range(30)]
    labeled = [2 ** comb(n, 2) for n in range(1, 14)]
    switching = [2 ** comb(n - 1, 2) for n in range(1, 15)]
    tournaments = [1, 1, 2, 4, 12, 56, 456, 6880, 191536, 9733056]
    even_graphs = [1, 1, 2, 3, 7, 16, 54, 243, 2038]
    score_sequences = [1, 1, 1, 2, 4, 9, 22, 59, 167, 490]

    families = [
        ("factorial_0_through_12", factorial),
        ("fibonacci_30_terms", fibonacci),
        ("catalan_30_terms", catalan),
        ("labeled_tournament_values", labeled),
        ("switching_values", switching),
        ("A000568_available_prefix", tournaments),
        ("A002854_available_prefix", even_graphs),
        ("A000571_available_prefix", score_sequences),
        ("THM2010_specA_n3_to_n6", [2, 3, 9, 28]),
        ("THM2010_H_value_count_n3_to_n6", [2, 3, 7, 19]),
        ("THM2010_cycle_vector_count_n3_to_n6", [2, 3, 9, 32]),
        ("THM2010_R_value_count_n3_to_n6", [2, 2, 6, 8]),
        ("THM2010_discriminant_value_count_n3_to_n6", [1, 2, 2, 5]),
        ("THM2010_arborescence_invariant_count_n3_to_n6", [2, 4, 12, 55]),
        ("THM2010_max_arborescence_n3_to_n6", [3, 10, 55, 333]),
        ("THM2010_max_R_n3_to_n6", [3, 3, 15, 15]),
        ("THM2010_metagraph_edges_n3_to_n6", [1, 5, 30, 290]),
        ("THM2010_WL_colors_n3_to_n6", [1, 2, 10, 34]),
    ]
    rows = [(name, *finite_signatures(values)) for name, values in families]

    labeled_support = set(labeled)
    switching_support = set(switching)
    require(labeled_support == switching_support,
            "labeled and switching value supports are the same shifted range")
    row_map = {name: (support, multiple, tax) for name, support, multiple, tax in rows}
    require(row_map["switching_values"][1] - row_map["labeled_tournament_values"][1] == 1,
            "the old plus-one is exactly a duplicate-one tax")
    for name in ("factorial_0_through_12", "fibonacci_30_terms", "catalan_30_terms"):
        require(row_map[name][2] == 1, f"{name} has exactly one duplicate-one tax")
    require(row_map["THM2010_R_value_count_n3_to_n6"][2] == F(1, 2),
            "THM-2010 R-value prefix duplicate-two tax")
    require(row_map["THM2010_discriminant_value_count_n3_to_n6"][2] == F(1, 2),
            "THM-2010 discriminant prefix duplicate-two tax")
    require(row_map["THM2010_max_R_n3_to_n6"][2] == F(2, 5),
            "THM-2010 max-R prefix has duplicate-three and duplicate-fifteen taxes")
    return rows


def abel_rhs(values: list[int], x: int) -> F:
    """A(x)/x + integral_1^x A(t)t^-2 dt for a finite support."""
    support = sorted({value for value in values if 1 <= value <= x})
    count = 1 if support and support[0] == 1 else 0
    previous = 1
    integral = F(0)
    for value in support:
        if value == 1:
            continue
        integral += count * (F(1, previous) - F(1, value))
        previous = value
        count += 1
    integral += count * (F(1, previous) - F(1, x))
    return F(len(support), x) + integral


def abel_stieltjes_audit() -> int:
    families = [
        [n for n in range(1, 101)],
        [n * (n + 1) // 2 for n in range(1, 40)],
        [n * n for n in range(1, 40)],
        [1, 2, 8, 64, 1024, 32768],
        [1, 1, 2, 3, 5, 8, 13, 21, 34],
    ]
    rows = 0
    for values in families:
        for x in range(1, max(values) + 7):
            support = {value for value in values if 1 <= value <= x}
            lhs = sum((F(1, value) for value in support), F(0))
            require(lhs == abel_rhs(values, x), "exact Abel--Stieltjes identity")
            rows += 1
    return rows


def block_bounds(values: list[int]) -> tuple[int, F]:
    support = sorted({value for value in values if value > 0})
    maximum = support[-1]
    rows = 0
    total = F(0)
    k = 0
    while 2**k <= maximum:
        block = [value for value in support if 2**k <= value < 2 ** (k + 1)]
        actual = sum((F(1, value) for value in block), F(0))
        lower = F(len(block), 2 ** (k + 1))
        upper = F(len(block), 2**k)
        require(lower <= actual <= upper, "multiplicative-block harmonic sandwich")
        total += F(len(block), 2**k)
        rows += 1
        k += 1
    return rows, total


def sieve_primes(limit: int) -> list[int]:
    flags = bytearray(b"\x01") * (limit + 1)
    flags[:2] = b"\x00\x00"
    for p in range(2, int(limit**0.5) + 1):
        if flags[p]:
            flags[p * p : limit + 1 : p] = b"\x00" * (((limit - p * p) // p) + 1)
    return [n for n in range(2, limit + 1) if flags[n]]


def logarithmic_block_audit() -> list[tuple[str, int, F]]:
    limit = 2**18
    naturals = list(range(1, limit + 1))
    odds = list(range(1, limit + 1, 2))
    primes = sieve_primes(limit)
    triangular = []
    n = 1
    while n * (n + 1) // 2 <= limit:
        triangular.append(n * (n + 1) // 2)
        n += 1
    fib = [1, 2]
    while fib[-1] + fib[-2] <= limit:
        fib.append(fib[-1] + fib[-2])
    bertrand_one = [int(n * log(n + 1) + 0.999999999999) for n in range(1, 30000)]
    bertrand_two = [int(n * log(n + 1) ** 2 + 0.999999999999) for n in range(1, 30000)]
    families = [
        ("naturals", naturals),
        ("odds", odds),
        ("primes", primes),
        ("triangular", triangular),
        ("fibonacci_support", fib),
        ("ceil_n_log_n", bertrand_one),
        ("ceil_n_log2_n", bertrand_two),
    ]
    rows = []
    for name, values in families:
        count, occupancy_sum = block_bounds(values)
        rows.append((name, count, occupancy_sum))
    require(bertrand_one[-1] / len(bertrand_one) > 10,
            "n log n example is genuinely superlinear")
    return rows


def master_figurate(s: int, d: int, m: int) -> int:
    return (s - 2) * comb(m + d - 2, d) + comb(m + d - 2, d - 1)


def polygonal(s: int, m: int) -> int:
    return ((s - 2) * m * m - (s - 4) * m) // 2


def polygonal_mass(s: int) -> mp.mpf:
    if s == 4:
        return mp.pi**2 / 6
    return (mp.mpf(2) / (s - 4)) * (
        mp.digamma(1) - mp.digamma(mp.mpf(2) / (s - 2))
    )


def master_surface_audit() -> tuple[int, list[tuple[int, mp.mpf]]]:
    rows = 0
    for s in range(3, 11):
        for d in range(2, 9):
            previous = None
            for m in range(1, 80):
                value = master_figurate(s, d, m)
                require(value > 0, "master figurate entries are positive")
                require(master_figurate(3, d, m) == comb(m + d - 1, d),
                        "simplex slice")
                require(master_figurate(s, 2, m) == polygonal(s, m),
                        "polygonal slice")
                if d > 2:
                    require(
                        value == sum(master_figurate(s, d - 1, i) for i in range(1, m + 1)),
                        "dimension-axis coning recurrence",
                    )
                if previous is not None:
                    require(value > previous, "each figurate column is strictly increasing")
                previous = value
                rows += 1

    # Exact termwise partial-fraction identity behind the digamma axis.
    for s in range(3, 13):
        for m in range(1, 200):
            if s == 4:
                require(F(1, polygonal(s, m)) == F(1, m * m), "square axis")
            else:
                b = s - 4
                rhs = F(2, b) * (F(1, m - F(b, s - 2)) - F(1, m))
                require(F(1, polygonal(s, m)) == rhs,
                        "polygonal partial-fraction identity")

    # Monotonicity on both axes: the common m=1 term stays one; every tail
    # term strictly decreases.
    for s in range(3, 10):
        for d in range(2, 8):
            for m in range(2, 80):
                require(master_figurate(s + 1, d, m) > master_figurate(s, d, m),
                        "shape-axis denominator increases")
                require(master_figurate(s, d + 1, m) > master_figurate(s, d, m),
                        "dimension-axis denominator increases")

    polygon_rows = [(s, polygonal_mass(s)) for s in range(3, 11)]
    require(abs(polygon_rows[0][1] - 2) < mp.mpf("1e-45"), "triangular mass is two")
    require(abs(polygon_rows[1][1] - mp.pi**2 / 6) < mp.mpf("1e-45"),
            "square mass is zeta two")
    require(abs(polygon_rows[2][1] - (3 * mp.log(3) - mp.pi / mp.sqrt(3)))
            < mp.mpf("1e-45"), "pentagonal mass closed form")
    require(abs(polygon_rows[3][1] - 2 * mp.log(2)) < mp.mpf("1e-45"),
            "hexagonal mass closed form")
    require(all(polygon_rows[i][1] > polygon_rows[i + 1][1]
                for i in range(len(polygon_rows) - 1)),
            "polygonal masses strictly descend")
    return rows, polygon_rows


def master_mass_partial(s: int, d: int, terms: int = 20000) -> mp.mpf:
    return mp.fsum(mp.mpf(1) / master_figurate(s, d, m) for m in range(1, terms + 1))


def harmonic(n: int, power: int = 1) -> F:
    return sum((F(1, k**power) for k in range(1, n + 1)), F(0))


def mp_fraction(value: F) -> mp.mpf:
    return mp.mpf(value.numerator) / value.denominator


def surface_partial_fractions(
    s: int, d: int
) -> tuple[str, dict[int, F], F, F | None]:
    """Exact coefficients for the full figurate reciprocal surface.

    With a=s-2 and beta=d/a-1,
      1/N = (d!/a)/[m(m+1)...(m+d-2)(m+beta)].
    The final return is (kind, simple coefficients indexed by integer poles,
    beta-pole/double-pole coefficient, beta).  In the collision case beta is
    an integer pole and is returned as None because the last coefficient is
    attached to (m+r)^-2.
    """
    a = s - 2
    beta = F(d, a) - 1
    scale = F(factorial(d), a)
    poles = range(d - 1)
    if beta.denominator != 1 or not 0 <= beta.numerator <= d - 2:
        simple = {
            j: scale * (-1) ** j
            / (factorial(j) * factorial(d - 2 - j) * (beta - j))
            for j in poles
        }
        beta_coefficient = scale
        for j in poles:
            beta_coefficient /= F(j) - beta
        require(sum(simple.values(), F(0)) + beta_coefficient == 0,
                "simple-pole coefficients cancel at infinity")
        return "digamma", simple, beta_coefficient, beta

    r = beta.numerator
    double = scale * (-1) ** r / (factorial(r) * factorial(d - 2 - r))
    simple = {}
    for j in poles:
        if j == r:
            simple[j] = double * (harmonic(r) - harmonic(d - 2 - r))
        else:
            simple[j] = (
                scale * (-1) ** j
                / (factorial(j) * factorial(d - 2 - j) * (r - j))
            )
    require(sum(simple.values(), F(0)) == 0,
            "double-pole simple coefficients cancel at infinity")
    return "pi2_collision", simple, double, None


def surface_closed_mass(s: int, d: int) -> tuple[str, mp.mpf]:
    kind, simple, last, beta = surface_partial_fractions(s, d)
    rational_part = -sum((coefficient * harmonic(j)
                          for j, coefficient in simple.items()), F(0))
    if kind == "digamma":
        require(beta is not None, "digamma pole is present")
        value = mp_fraction(rational_part) - mp_fraction(last) * (
            mp.digamma(mp_fraction(beta + 1)) - mp.digamma(1)
        )
    else:
        r = int(F(d, s - 2) - 1)
        value = (
            mp_fraction(rational_part - last * harmonic(r, 2))
            + mp_fraction(last) * mp.pi**2 / 6
        )
    return kind, value


def full_surface_closed_form_audit() -> tuple[int, list[tuple[int, int, str, mp.mpf, str]]]:
    rows = 0
    for s in range(3, 13):
        for d in range(2, 11):
            kind, simple, last, beta = surface_partial_fractions(s, d)
            for m in range(1, 80):
                rhs = sum((coefficient / F(m + j)
                           for j, coefficient in simple.items()), F(0))
                if kind == "digamma":
                    require(beta is not None, "simple beta pole")
                    rhs += last / (F(m) + beta)
                else:
                    r = F(d, s - 2).numerator - 1
                    rhs += last / F((m + r) ** 2)
                require(rhs == F(1, master_figurate(s, d, m)),
                        "full-surface partial fraction identity")
                rows += 1

            _, closed_value = surface_closed_mass(s, d)
            direct_value = mp.nsum(
                lambda m: 1 / (
                    (s - 2) * mp.binomial(m + d - 2, d)
                    + mp.binomial(m + d - 2, d - 1)
                ),
                [1, mp.inf],
            )
            require(abs(closed_value - direct_value) < mp.mpf("1e-40"),
                    "closed surface evaluator versus independent infinite sum")

    examples = [
        (4, 2, "pi^2/6", mp.pi**2 / 6),
        (4, 3, "18-24log(2)", 18 - 24 * mp.log(2)),
        (4, 4, "21-2pi^2", 21 - 2 * mp.pi**2),
        (4, 6, "15pi^2-1175/8", 15 * mp.pi**2 - mp.mpf(1175) / 8),
        (5, 3, "pi^2/3-2", mp.pi**2 / 3 - 2),
        (5, 4, "282/5-(162/5)log3-(18sqrt3/5)pi",
         mp.mpf(282) / 5 - mp.mpf(162) / 5 * mp.log(3)
         - mp.mpf(18) / 5 * mp.sqrt(3) * mp.pi),
        (5, 6, "1205/18-(20/3)pi^2",
         mp.mpf(1205) / 18 - mp.mpf(20) / 3 * mp.pi**2),
        (6, 3, "(72/5)log2-(12/5)pi-6/5",
         mp.mpf(72) / 5 * mp.log(2) - mp.mpf(12) / 5 * mp.pi - mp.mpf(6) / 5),
        (6, 4, "pi^2/2-15/4", mp.pi**2 / 2 - mp.mpf(15) / 4),
        (7, 5, "2pi^2/3-49/9", mp.mpf(2) / 3 * mp.pi**2 - mp.mpf(49) / 9),
        (8, 3, "(8/3)log2-2/3", mp.mpf(8) / 3 * mp.log(2) - mp.mpf(2) / 3),
        (8, 6, "5pi^2/6-1025/144",
         mp.mpf(5) / 6 * mp.pi**2 - mp.mpf(1025) / 144),
    ]
    output = []
    for s, d, expression, expected in examples:
        kind, value = surface_closed_mass(s, d)
        require(abs(value - expected) < mp.mpf("1e-45"),
                f"full-surface example s={s},d={d}")
        output.append((s, d, kind, value, expression))
    return rows, output


def power_sum(n: int, p: int) -> int:
    return sum(k**p for k in range(1, n + 1))


def powersum_axis_audit() -> list[tuple[int, mp.mpf, str]]:
    terms = 40000
    rows = []
    for p in range(1, 7):
        cumulative = 0
        reciprocal_terms = []
        for n in range(1, terms + 1):
            cumulative += n**p
            reciprocal_terms.append(mp.mpf(1) / cumulative)
        value = mp.fsum(reciprocal_terms)
        note = "numeric_partial"
        if p == 1:
            closed = mp.mpf(2)
            note = "2"
        elif p == 2:
            closed = 18 - 24 * mp.log(2)
            note = "18-24log(2)"
        elif p == 3:
            closed = 4 * mp.pi**2 / 3 - 12
            note = "4pi^2/3-12"
        elif p == 4:
            r = (mp.sqrt(21) - 3) / 6
            d_one_minus_r = mp.digamma(1 - r) + mp.euler
            d_two_plus_r = mp.digamma(2 + r) + mp.euler
            closed = (
                (-270 + 480 * mp.log(2)) / 7
                - mp.mpf(90) / 7 * (d_one_minus_r + d_two_plus_r)
            )
            note = "algebraic_argument_digamma_exact"
        elif p == 5:
            r = (mp.sqrt(3) - 1) / 2
            closed = (
                60 - 4 * mp.pi**2
                - 8 * mp.sqrt(3) * mp.pi / mp.tan(mp.pi * r)
            )
            note = "60-4pi^2-8sqrt(3)pi*cot(pi(sqrt3-1)/2)"
        else:
            closed = None
        if closed is not None:
            require(abs(value - closed) < mp.mpf("0.0001"),
                    f"power-sum p={p} closed form")
        rows.append((p, value, note))

    for n in range(1, 200):
        # 1 / sum k^2 = 6/[n(n+1)(2n+1)].
        lhs2 = F(1, power_sum(n, 2))
        rhs2 = F(6, n) + F(6, n + 1) - F(24, 2 * n + 1)
        require(lhs2 == rhs2, "square-pyramidal partial fraction")
        # sum k^3 = T_n^2.
        require(power_sum(n, 3) == (n * (n + 1) // 2) ** 2,
                "cubic Faulhaber triangular-square identity")
        require(
            F(1, power_sum(n, 4))
            == F(270 * (2 * n + 1), 7 * (3 * n**2 + 3 * n - 1))
            + F(480, 7 * (2 * n + 1)) - F(30, n) - F(30, n + 1),
            "fourth-power reciprocal rational partial fractions",
        )
        # Fifth-power Faulhaber denominator and rational partial fractions.
        require(
            power_sum(n, 5)
            == n**2 * (n + 1) ** 2 * (2 * n**2 + 2 * n - 1) // 12,
            "fifth-power Faulhaber factorization",
        )
        require(
            F(1, power_sum(n, 5))
            == F(48, 2 * n**2 + 2 * n - 1) - F(12, n**2) - F(12, (n + 1) ** 2),
            "fifth-power reciprocal partial fractions",
        )
    return rows


def triangular_balance_tower_audit() -> tuple[mp.mpf, str]:
    partial = F(0)
    for n in range(1, 501):
        center = n * (n + 1)
        left = sum(range(center - n, center + 1))
        right = sum(range(center + 1, center + n + 1))
        side_value = n * (n + 1) * (2 * n + 1) // 2
        require(left == right == side_value == 3 * power_sum(n, 2),
                "first triangular balance tower has equal side values")
        partial += F(1, side_value)
        finite_closed = (
            8 * harmonic(n) - 8 * harmonic(2 * n + 1)
            + F(2, n + 1) + 6
        )
        require(partial == finite_closed,
                "triangular balance reciprocal finite closed form")
    closed = 6 - 8 * mp.log(2)
    require(abs(mp_fraction(partial) - closed) < mp.mpf("0.002"),
            "triangular balance mass approaches 6-8log2")
    return closed, "6-8log(2)"


def maximum_cyclic_triangle_audit() -> tuple[mp.mpf, mp.mpf, mp.mpf]:
    values = []
    for n in range(3, 1001):
        if n % 2:
            maximum = n * (n * n - 1) // 24
        else:
            maximum = n * (n * n - 4) // 24
        values.append(maximum)
        if n > 3:
            k = (n - 1) // 2
            require(values[-1] - values[-2] == k * (k + 1) // 2,
                    "maximum cyclic-triangle increments repeat triangular numbers")

    for k in range(1, 500):
        odd_n = 2 * k + 1
        odd_maximum = odd_n * (odd_n * odd_n - 1) // 24
        require(odd_maximum == power_sum(k, 2),
                "odd-order maximum cyclic triangles are square-pyramidal")
        even_n = 2 * k + 2
        even_maximum = even_n * (even_n * even_n - 4) // 24
        require(even_maximum - odd_maximum == k * (k + 1) // 2,
                "odd-to-even maximum cyclic increment")

    odd_mass = 18 - 24 * mp.log(2)
    even_mass = mp.mpf(3) / 4
    all_mass = odd_mass + even_mass
    direct_odd = mp.nsum(lambda k: 6 / (k * (k + 1) * (2 * k + 1)), [1, mp.inf])
    direct_even = mp.nsum(lambda k: 3 / (k * (k * k - 1)), [2, mp.inf])
    require(abs(direct_odd - odd_mass) < mp.mpf("1e-45"),
            "odd maximum-cyclic reciprocal mass")
    require(abs(direct_even - even_mass) < mp.mpf("1e-45"),
            "even maximum-cyclic reciprocal mass")
    return odd_mass, even_mass, all_mass


def kakeya_achievement_audit() -> tuple[list[tuple[int, int, F]], list[tuple[int, str]]]:
    """Exact tail-domination audit for simplex and geometric reciprocal atoms."""
    simplex_rows = []
    for k in range(2, 10):
        total = F(k, k - 1)
        first_overlap = 2 * k - 2
        early = [F(1, comb(n, k)) for n in range(k, first_overlap)]
        interval_length = F(k, (k - 1) * comb(2 * k - 3, k - 1))

        for n in range(k, first_overlap + 30):
            atom = F(1, comb(n, k))
            remainder = F(k, (k - 1) * comb(n, k - 1))
            require(atom / remainder == F(k - 1, n - k + 1),
                    "simplex exact atom-to-tail ratio")
            require((atom > remainder) == (n < first_overlap),
                    "simplex Kakeya gap/overlap threshold")

        starts = [
            sum((early[j] for j in range(len(early)) if mask >> j & 1), F(0))
            for mask in range(1 << len(early))
        ]
        intervals = sorted((start, start + interval_length) for start in starts)
        require(len(intervals) == 2 ** (k - 2),
                "simplex achievement interval count")
        require(all(intervals[j][1] < intervals[j + 1][0]
                    for j in range(len(intervals) - 1)),
                "simplex achievement intervals are pairwise disjoint")
        require(intervals[0][0] == 0 and intervals[-1][1] == total,
                "simplex achievement-set endpoints")
        simplex_rows.append((k, len(intervals), interval_length))

    geometric_rows = []
    for base in range(2, 10):
        atom = F(1, base**7)
        remainder = F(1, (base - 1) * base**7)
        require(atom / remainder == base - 1,
                "geometric achievement atom-to-tail ratio")
        topology = "interval_[0,2]" if base == 2 else "Cantor_dimension_log_base(2)"
        geometric_rows.append((base, topology))
    return simplex_rows, geometric_rows


def g_diagonal(n: int) -> int:
    value = (
        F(n**4, 96) - F(n**3, 16) + F(n * n, 3) + F(n, 32) + F(53, 64)
        + (-1) ** n * (F(11, 64) - F(n, 32))
    )
    require(value.denominator == 1, "G quasi-polynomial is integral")
    return value.numerator


def repo_sequence_atlas() -> list[tuple[str, mp.mpf, str]]:
    terms = 30000
    moser = [1 + comb(n, 2) + comb(n, 4) for n in range(1, terms + 1)]
    g_values = [g_diagonal(n) for n in range(terms)]
    fib = [1, 1]
    for _ in range(80):
        fib.append(fib[-1] + fib[-2])

    quarter_limit = 500
    quarter_support = {
        (n * n) // 4 for n in range(2, 2 * quarter_limit + 2)
    }
    square_oblong_support = {
        k * k for k in range(1, quarter_limit + 1)
    } | {
        k * (k + 1) for k in range(1, quarter_limit + 1)
    }
    require(quarter_support == square_oblong_support,
            "quarter-square support is the disjoint square/oblong union")
    require(not ({k * k for k in range(1, quarter_limit + 1)}
                 & {k * (k + 1) for k in range(1, quarter_limit + 1)}),
            "positive square and oblong supports are disjoint")

    phi = list(range(200001))
    for p in range(2, len(phi)):
        if phi[p] == p:
            for multiple in range(p, len(phi), p):
                phi[multiple] -= phi[multiple] // p
    farey_totals = []
    totient_sum = 0
    for n in range(1, len(phi)):
        totient_sum += 1 if n == 1 else phi[n]
        farey_totals.append(2 * totient_sum - 1)

    families = [
        ("Moser_A000127", moser, "degree_4;partial_30000"),
        ("polygonal_diagonal_G", g_values, "degree_4_quasipolynomial;partial_30000"),
        ("quarter_square_A002620", [], "exact_support=zeta(2)+1"),
        ("Fibonacci_support", fib, "duplicate_1_removed"),
        ("Farey_endpoint_total", farey_totals, "2sum_phi-1;partial_200000"),
        ("A000568_support_prefix", [1, 1, 2, 4, 12, 56, 456, 6880, 191536],
         "available_prefix;duplicate_1_removed"),
        ("self_line_orbits_prefix", [2, 3, 22, 101, 852, 5582],
         "THM853_n5_to_n10_prefix"),
    ]
    rows = []
    for name, values, note in families:
        if name == "quarter_square_A002620":
            mass = mp.pi**2 / 6 + 1
        else:
            support = sorted({value for value in values if value > 0})
            mass = mp.fsum(mp.mpf(1) / value for value in support)
        rows.append((name, mass, note))
    return rows


def run_support_filtration_audit() -> tuple[mp.mpf, mp.mpf]:
    def row_sum(j: int, r: int) -> int:
        return sum(comb(r + 1, 2 * i) for i in range(j + 1) if 2 * i <= r + 1)

    for r in range(80):
        values = [row_sum(j, r) for j in range(1, 42)]
        require(all(values[j] <= values[j + 1] for j in range(len(values) - 1)),
                "run-row counts increase with allowed runs")
        require(values[-1] == 2**r, "full even-binomial row is a power of two")
        require(row_sum(1, r) == 1 + comb(r + 1, 2), "first run-row formula")

    row_one = mp.nsum(lambda r: 1 / (1 + r * (r + 1) / 2), [0, mp.inf])
    row_one_closed = 2 * mp.pi / mp.sqrt(7) * mp.tanh(mp.pi * mp.sqrt(7) / 2)
    require(abs(row_one - row_one_closed) < mp.mpf("1e-45"),
            "first run-row hyperbolic mass")

    for k in range(1, 500):
        require((2 * k) ** 2 // 4 + 1 == k * k + 1,
                "even lifted quarter-square branch")
        require((2 * k + 1) ** 2 // 4 + 1 == k * (k + 1) + 1,
                "odd lifted quarter-square branch")
    diagonal_one = (
        1
        + mp.nsum(lambda k: 1 / (k**2 + 1), [1, mp.inf])
        + mp.nsum(lambda k: 1 / (k * (k + 1) + 1), [1, mp.inf])
    )
    diagonal_one_closed = (
        mp.pi / 2 * mp.coth(mp.pi)
        + mp.pi / mp.sqrt(3) * mp.tanh(mp.pi * mp.sqrt(3) / 2)
        - mp.mpf(1) / 2
    )
    require(abs(diagonal_one - diagonal_one_closed) < mp.mpf("1e-45"),
            "first run-diagonal hyperbolic mass")
    return row_one_closed, diagonal_one_closed


def gauss_triangular_theta_audit() -> tuple[mp.mpf, mp.mpf, mp.mpf]:
    q = mp.mpf(1) / 2
    triangular_sum = mp.fsum(q ** (n * (n + 1) // 2) for n in range(80))
    q_product = mp.mpf(1)
    even_product = mp.mpf(1)
    for n in range(1, 500):
        q_product *= 1 - q**n
        even_product *= 1 - q ** (2 * n)
    gauss_product = even_product**2 / q_product
    theta_form = mp.jtheta(2, 0, mp.sqrt(q)) / (2 * q ** (mp.mpf(1) / 8))
    require(abs(triangular_sum - gauss_product) < mp.mpf("1e-45"),
            "Gauss triangular-number product")
    require(abs(triangular_sum - theta_form) < mp.mpf("1e-45"),
            "theta-two form")
    return triangular_sum, gauss_product, theta_form


def ladder_ratio_sum_product_audit() -> tuple[F, mp.mpf, mp.mpf]:
    # Ratios 2k-1 have harmonic mass divergence.  Their partial sums are
    # k^2-1, whose reciprocal mass telescopes to 3/4.  Their partial products
    # are odd double factorials and have an erf closed form.
    for n in range(2, 500):
        partial = sum((F(1, k * k - 1) for k in range(2, n + 1)), F(0))
        closed = F(3, 4) - F(1, 2 * n) - F(1, 2 * (n + 1))
        require(partial == closed, "oblong reciprocal telescoping")

    product_mass = mp.mpf(0)
    odd_double_factorial = 1
    for k in range(1, 100):
        odd_double_factorial *= 2 * k - 1
        if k >= 2:
            product_mass += mp.mpf(1) / odd_double_factorial
    product_closed = (
        mp.e ** (mp.mpf(1) / 2)
        * mp.sqrt(mp.pi / 2)
        * mp.erf(1 / mp.sqrt(2))
        - 1
    )
    require(abs(product_mass - product_closed) < mp.mpf("1e-45"),
            "odd-double-factorial reciprocal erf form")
    odd_partial = mp.fsum(mp.mpf(1) / (2 * k - 1) for k in range(1, 200000))
    return F(3, 4), product_mass, odd_partial


def tournament_analytic_reversal_audit() -> int:
    rows = 0
    for n in range(1, 80):
        # Consecutive labeled-EGF terms have absolute ratio
        # 2^n |z|/(n+1), so every z != 0 eventually makes terms grow.
        require(F(2**n, n + 1) > 1 if n >= 2 else True,
                "labeled tournament EGF ratio eventually exceeds one at z=1")

        # Orbit count V_n is at least labeled_count/n!.  This lower bound is
        # already superpolynomial, so V_n/n^s cannot tend to zero for fixed s.
        lower = F(2 ** comb(n, 2), factorial(n))
        if n >= 20:
            require(lower > n**10, "Burnside orbit lower bound beats a fixed polynomial")
        rows += 1
    return rows


def tournament_and_carrier_audit() -> None:
    print("TOURNAMENT_AND_CARRIER_AUDIT")
    print("faithful_vertices=occupied_integer_values_or_logarithmic_blocks;not_sequence_indices")
    print("predicate_preserved=support_membership,reciprocal_mass,block_occupancy,convergence")
    print("multiplicity_sidecar=collision_tax;index_shift_and_repetition_are_quotiented_out")
    print("candidate_vertices=indices,values,blocks,figurate_axes,Dirichlet_exponents,proof_obligations")
    print("pair_observable=finite_support_mass_difference;gauge=exact_rational_then_family_name")
    print("induced_tournament=transitive;cycles=0;scc_sizes=all_1;hamiltonian_paths=1")
    print("why_unfaithful=finite_mass_order_does_not_determine_tail_convergence_or_closed_form")
    print("tie_path=mass_sorted_family_list;useful_as_atlas_only_not_as_proof_quotient")


def main() -> None:
    assert_nodes = optimization_safety_probe()
    collision_rows = collision_tax_audit()
    abel_rows = abel_stieltjes_audit()
    block_rows = logarithmic_block_audit()
    surface_rows, polygon_rows = master_surface_audit()
    closed_surface_rows, closed_surface_examples = full_surface_closed_form_audit()
    power_rows = powersum_axis_audit()
    balance_mass, balance_closed = triangular_balance_tower_audit()
    max_c3_odd, max_c3_even, max_c3_all = maximum_cyclic_triangle_audit()
    simplex_achievement_rows, geometric_achievement_rows = kakeya_achievement_audit()
    atlas_rows = repo_sequence_atlas()
    run_row_mass, run_diagonal_mass = run_support_filtration_audit()
    theta_sum, theta_product, theta_form = gauss_triangular_theta_audit()
    oblong_mass, double_factorial_mass, odd_partial = ladder_ratio_sum_product_audit()
    reversal_rows = tournament_analytic_reversal_audit()

    print("THM-2000 SUPPORT-HARMONIC ABEL-DINI FIGURATE-SURFACE REFEREE")
    print(f"optimization_sensitive_assert_nodes={assert_nodes}")
    print("SUPPORT_VERSUS_MULTIPLICITY")
    for name, support, multiple, tax in collision_rows:
        print(f"{name}:support={support};multiple={multiple};collision_tax={tax}")
    print("support_law=sigma_multi=sigma_set+sum_m(max(multiplicity(m)-1,0)/m)")
    print("labeled_tournament_support=switching_support;old_plus_one_is_duplicate_one_tax")

    print(f"ABEL_STIELTJES_EXACT_ROWS={abel_rows}")
    print("abel_law=sum_(m<=x,m_in_A)1/m=A(x)/x+integral_1^x A(t)/t^2 dt")
    print("block_law=for_Nk=#(A intersect [2^k,2^(k+1))), Nk/2^(k+1)<=block_mass<=Nk/2^k")
    for name, count, occupancy in block_rows:
        print(f"block_family={name};blocks={count};finite_occupancy_sum={occupancy}")
    print("bertrand_repair=n_log_n_is_superlinear_but_reciprocal_diverges;n_log2_n_converges")

    print("KAKEYA_ACHIEVEMENT_SETS")
    for k, components, length in simplex_achievement_rows:
        print(f"simplex_k={k};components={components};common_interval_length={length};total_mass={F(k,k-1)}")
    for base, topology in geometric_achievement_rows:
        print(f"powers_base={base};atom_tail_ratio={base-1};topology={topology}")
    print("equal_mass_topology=simplex_k_has_2^(k-2)_intervals;powers_k_is_Cantor_for_k>=3;both_k2_fill_[0,2]")

    print(f"MASTER_FIGURATE_EXACT_ROWS={surface_rows}")
    print("N(s,d,m)=(s-2)C(m+d-2,d)+C(m+d-2,d-1)")
    print("mass_surface=d(d-1) double_integral (1-t)^(d-2)u^(d-1)/(1-tu^(s-2)) dtdu")
    print("simplex_axis=M(3,d)=d/(d-1);polygonal_axis=digamma_formula;both_axes_strictly_descend_to_1")
    for s, value in polygon_rows:
        print(f"polygonal_s={s};mass={mp.nstr(value, 30)}")
    for d in range(2, 7):
        partial = master_mass_partial(4, d)
        print(f"master_surface_s=4,d={d};partial20000={mp.nstr(partial, 20)}")

    print(f"FULL_SURFACE_PARTIAL_FRACTION_ROWS={closed_surface_rows}")
    print("off_resonance=rational_plus_rational_times_digamma(d/(s-2));Gauss_reduces_to_pi_and_log_algebraic")
    print("resonance=(s-2)>=2_and_(s-2)|d;one_double_pole;mass=rational_plus_nonzero_rational*pi^2")
    for s, d, kind, value, expression in closed_surface_examples:
        print(f"closed_surface_s={s},d={d};kind={kind};mass={mp.nstr(value, 30)};closed={expression}")

    print("POWER_SUM_FAULHABER_AXIS")
    for p, value, note in power_rows:
        print(f"power_p={p};partial40000={mp.nstr(value, 25)};closed={note}")
    print(f"triangular_balance_common_side_mass={mp.nstr(balance_mass, 30)};closed={balance_closed}")
    print("balance_semantics=support_counts_the_equal_left_right_value_once;two_labeled_sides_would_double")
    print(f"max_cyclic_triangles_odd_mass={mp.nstr(max_c3_odd, 30)};closed=18-24log(2)")
    print(f"max_cyclic_triangles_even_mass={mp.nstr(max_c3_even, 30)};closed=3/4")
    print(f"max_cyclic_triangles_all_mass={mp.nstr(max_c3_all, 30)};closed=75/4-24log(2)")
    print("max_cyclic_structure=odd_orders_are_square_pyramidal;successive_increments_repeat_triangular_numbers")

    print("REPO_SUPPORT_MASS_ATLAS")
    for name, value, note in atlas_rows:
        print(f"{name}:support_mass={mp.nstr(value, 25)};note={note}")
    atlas = {name: value for name, value, _ in atlas_rows}
    require(run_row_mass > atlas["Moser_A000127"] > 2,
            "row run-filtration masses descend toward powers of two")
    require(run_diagonal_mass > atlas["polygonal_diagonal_G"] > atlas["Fibonacci_support"],
            "diagonal run-filtration masses descend toward Fibonacci support")
    print(f"run_row_R1_exact_mass={mp.nstr(run_row_mass, 30)};closed=2pi/sqrt7*tanh(pi*sqrt7/2)")
    print(f"run_diagonal_G1_exact_mass={mp.nstr(run_diagonal_mass, 30)};closed=pi/2*coth(pi)+pi/sqrt3*tanh(pi*sqrt3/2)-1/2")
    print("run_filtration=row_counts_increase_to_powers2_and_diagonal_counts_increase_to_Fibonacci;reciprocal_masses_decrease")

    print("GAUSS_TRIANGULAR_THETA")
    print(f"sum_q_triangular={mp.nstr(theta_sum, 30)}")
    print(f"q_product_form={mp.nstr(theta_product, 30)};theta2_form={mp.nstr(theta_form, 30)}")
    print("labeled_and_switching_support_mass_at_q_half_are_equal_to_this_exact_Gauss_product")

    print("RATIO_SUM_PRODUCT_TRIAD")
    print(f"odd_ratio_support_partial200000={mp.nstr(odd_partial, 25)};status=diverges")
    print(f"oblong_k2_minus1_support_mass={oblong_mass}")
    print(f"odd_double_factorial_k_ge2_mass={mp.nstr(double_factorial_mass, 30)}")

    print(f"TOURNAMENT_ANALYTIC_REVERSAL_ROWS={reversal_rows}")
    print("index_Dirichlet=sum V_n/n^s diverges_for_every_fixed_s")
    print("labeled_EGF=sum 2^C(n,2)z^n/n! has_radius_zero")
    print("reciprocal_value_Dirichlet=sum V_n^(-s) has_abscissa_zero")
    print("reciprocal_OGF=sum z^n/V_n is_entire")
    print("H_spectrum_divergence=OPEN;THM-1370-h-spectrum-file_proves_only_omissions_7_21_and_finite_coverage")
    print("THM1127_affine_phase_law=N(K+194040)=N(K)+121304 => support_contains_AP => reciprocal_diverges")

    tournament_and_carrier_audit()
    print("SCOPE=no_irrationality_or_transcendence_claim_for_unnamed_numeric_constants")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
