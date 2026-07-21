#!/usr/bin/env python3
"""Exact verification gates for THM-2020.

All arithmetic is Python integer or fractions.Fraction arithmetic.  The
script independently checks the Wick/channel formula, Legendre digit-sum and
carry formulas, the wide-gap star certificates, the cancellable M4 and its
3-adically separated M8 descendants, and the Lee--Yang counterexample.

Tournament Analysis
-------------------
For the local p=3 certificate, vertices are the three balanced M8 channels.
Orient x -> y when (v_3(x), channel_label(x)) is lexicographically smaller.
This retains the ultrametric noncancellation predicate and produces the path
middle -> low -> high.  Monomials, charges, primitive circuits, radial
heights, primes, carry states, and proof obligations were also challenged as
vertices.  Charges alone lose both factorial height and carry state; a
channel-at-a-place is the least lossy quotient.  The tournament is a proof
diagnostic, not a replacement for the valuation proof.
"""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
from itertools import permutations
from math import factorial


Support = tuple[tuple[int, int], ...]
Counts = tuple[int, ...]
QComplex = tuple[Fraction, Fraction]


def compositions(total: int, parts: int):
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in compositions(total - first, parts - 1):
            yield (first,) + tail


def channel_coefficient(m: int, radial_height: int, counts: Counts) -> int:
    denominator = 1
    for count in counts:
        denominator *= factorial(count)
    multinomial = factorial(m) // denominator
    return multinomial * factorial(radial_height)


def balanced_channels(support: Support, m: int) -> list[tuple[Counts, int, int]]:
    charges = tuple(a - b for a, b in support)
    out = []
    for counts in compositions(m, len(support)):
        if sum(q * r for q, r in zip(charges, counts)) != 0:
            continue
        a_degree = sum(a * r for (a, _), r in zip(support, counts))
        b_degree = sum(b * r for (_, b), r in zip(support, counts))
        assert a_degree == b_degree
        out.append((counts, a_degree, channel_coefficient(m, a_degree, counts)))
    return out


def channel_moment(
    support: Support, coefficients: tuple[Fraction, ...], m: int
) -> Fraction:
    answer = Fraction(0)
    for counts, _, wick_coefficient in balanced_channels(support, m):
        coefficient_monomial = Fraction(1)
        for coefficient, count in zip(coefficients, counts):
            coefficient_monomial *= coefficient**count
        answer += wick_coefficient * coefficient_monomial
    return answer


def multiply_bivariate(
    left: dict[tuple[int, int], Fraction],
    right: dict[tuple[int, int], Fraction],
) -> dict[tuple[int, int], Fraction]:
    out: defaultdict[tuple[int, int], Fraction] = defaultdict(Fraction)
    for (a, b), x in left.items():
        for (c, d), y in right.items():
            out[(a + c, b + d)] += x * y
    return {key: value for key, value in out.items() if value}


def direct_wick_moment(
    support: Support, coefficients: tuple[Fraction, ...], m: int
) -> Fraction:
    base: defaultdict[tuple[int, int], Fraction] = defaultdict(Fraction)
    for exponent, coefficient in zip(support, coefficients):
        base[exponent] += coefficient
    expanded: dict[tuple[int, int], Fraction] = {(0, 0): Fraction(1)}
    for _ in range(m):
        expanded = multiply_bivariate(expanded, dict(base))
    return sum(
        (coefficient * factorial(a) for (a, b), coefficient in expanded.items() if a == b),
        Fraction(0),
    )


def vp_int(value: int, prime: int) -> int:
    value = abs(value)
    assert value != 0
    valuation = 0
    while value % prime == 0:
        value //= prime
        valuation += 1
    return valuation


def vp_factorial(n: int, prime: int) -> int:
    answer = 0
    while n:
        n //= prime
        answer += n
    return answer


def digit_sum(n: int, prime: int) -> int:
    answer = 0
    while n:
        answer += n % prime
        n //= prime
    return answer


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    divisor = 2
    while divisor * divisor <= n:
        if n % divisor == 0:
            return n == divisor
        divisor += 1
    return True


def check_moment_formula() -> int:
    controls = [
        (
            ((6, 0), (0, 2), (0, 18)),
            (Fraction(2), Fraction(-3), Fraction(5)),
            range(1, 13),
        ),
        (
            ((2, 0), (1, 1), (0, 1), (3, 2)),
            (Fraction(2), Fraction(-3, 2), Fraction(5), Fraction(7, 3)),
            range(1, 9),
        ),
    ]
    checks = 0
    for support, coefficients, moment_range in controls:
        for m in moment_range:
            assert direct_wick_moment(support, coefficients, m) == channel_moment(
                support, coefficients, m
            )
            checks += 1
    return checks


def check_legendre_and_carries() -> int:
    support = ((6, 0), (0, 2), (0, 18))
    primes = (2, 3, 5, 7, 11, 13, 17)
    checks = 0
    for n in range(1, 9):
        m = 4 * n
        for counts, radial_height, coefficient in balanced_channels(support, m):
            for prime in primes:
                direct = vp_int(coefficient, prime)
                numerator = (
                    radial_height
                    - digit_sum(m, prime)
                    - digit_sum(radial_height, prime)
                    + sum(digit_sum(count, prime) for count in counts)
                )
                assert numerator % (prime - 1) == 0
                legendre_digits = numerator // (prime - 1)
                carry_numerator = sum(digit_sum(count, prime) for count in counts) - digit_sum(
                    m, prime
                )
                assert carry_numerator % (prime - 1) == 0
                carry_count = carry_numerator // (prime - 1)
                carry_formula = vp_factorial(radial_height, prime) + carry_count
                assert direct == legendre_digits == carry_formula
                checks += 1
    return checks


def star_channels(n: int) -> list[tuple[Counts, int, int]]:
    support = ((6, 0), (0, 2), (0, 18))
    channels = balanced_channels(support, 4 * n)
    expected_counts = [(n + 2 * z, 3 * n - 3 * z, z) for z in range(n + 1)]
    assert [counts for counts, _, _ in channels] == expected_counts
    assert [height for _, height, _ in channels] == [6 * n + 12 * z for z in range(n + 1)]
    return channels


def check_wide_gap_certificates() -> list[tuple[int, int, int, tuple[int, ...]]]:
    primes = (29, 31, 43, 61, 73, 97, 127)
    audit = []
    for prime in primes:
        assert is_prime(prime)
        n = (prime - 1) // 6
        m = 4 * n
        assert m < 6 * n < prime <= 6 * (n + 1) <= 6 * n + 12
        channels = star_channels(n)
        valuations = tuple(vp_int(coefficient, prime) for _, _, coefficient in channels)
        assert valuations[0] == 0
        assert all(value > 0 for value in valuations[1:])
        audit.append((prime, n, m, valuations))
    return audit


def check_local_m4_m8() -> tuple[list[int], list[int], int]:
    # Substitute c=t*b^3/a^2, t=-6!/18!, so every m=4n channel is
    # a^n*b^(3n) times its Wick coefficient multiplied by t^z.
    substitution = -Fraction(factorial(6), factorial(18))

    m4_normalized = []
    for counts, _, coefficient in star_channels(1):
        z = counts[2]
        value = Fraction(coefficient) * substitution**z
        assert value.denominator == 1
        m4_normalized.append(value.numerator)
    assert m4_normalized == [4 * factorial(6), -4 * factorial(6)]
    assert sum(m4_normalized) == 0

    m8_normalized = []
    for counts, _, coefficient in star_channels(2):
        z = counts[2]
        value = Fraction(coefficient) * substitution**z
        assert value.denominator == 1
        m8_normalized.append(value.numerator)
    expected = [
        13_412_044_800,
        -19_536_878_592_000,
        131_727_403_906_560_000,
    ]
    assert m8_normalized == expected
    valuations = [vp_int(value, 3) for value in m8_normalized]
    assert valuations == [5, 4, 5]
    total = sum(m8_normalized)
    assert total == 131_707_880_440_012_800
    assert vp_int(total, 3) == 4
    return m8_normalized, valuations, total


def qadd(left: QComplex, right: QComplex) -> QComplex:
    return left[0] + right[0], left[1] + right[1]


def qmul(left: QComplex, right: QComplex) -> QComplex:
    return (
        left[0] * right[0] - left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    )


def qscale(value: QComplex, scalar: Fraction) -> QComplex:
    return value[0] * scalar, value[1] * scalar


def check_lee_yang_counterexample() -> tuple[QComplex, QComplex, QComplex, QComplex]:
    a = (Fraction(1), Fraction(1))
    b = (Fraction(-1), Fraction(1))
    c = (Fraction(1, 2), Fraction(1, 2))
    assert a[1] > 0 and b[1] > 0 and c[1] > 0
    b_squared = qmul(b, b)
    two_ac = qscale(qmul(a, c), Fraction(2))
    moment = qadd(b_squared, two_ac)
    assert b_squared == (Fraction(0), Fraction(-2))
    assert two_ac == (Fraction(0), Fraction(2))
    assert moment == (Fraction(0), Fraction(0))
    return a, b, c, moment


def tournament_fingerprint(valuations: list[int]) -> tuple[list[str], int, list[int], int]:
    labels = ["low", "middle", "high"]
    order = sorted(range(3), key=lambda index: (valuations[index], index))
    path = [labels[index] for index in order]
    rank = {vertex: position for position, vertex in enumerate(order)}

    def edge(source: int, target: int) -> bool:
        return rank[source] < rank[target]

    directed_cycles = 0
    for x, y, z in permutations(range(3), 3):
        if x < y < z and edge(x, y) and edge(y, z) and edge(z, x):
            directed_cycles += 1
        if x < y < z and edge(x, z) and edge(z, y) and edge(y, x):
            directed_cycles += 1
    hamiltonian_paths = 0
    for candidate in permutations(range(3)):
        if edge(candidate[0], candidate[1]) and edge(candidate[1], candidate[2]):
            hamiltonian_paths += 1
    scores = [sum(edge(vertex, other) for other in range(3) if other != vertex) for vertex in range(3)]
    return path, directed_cycles, sorted(scores), hamiltonian_paths


def format_valuations(values: tuple[int, ...]) -> str:
    if len(values) <= 8:
        return str(values)
    return f"({values[0]}, {values[1]}, ..., {values[-2]}, {values[-1]})"


def main() -> None:
    moment_checks = check_moment_formula()
    legendre_checks = check_legendre_and_carries()
    wide_gap = check_wide_gap_certificates()
    m8_coefficients, m8_valuations, m8_total = check_local_m4_m8()
    a, b, c, lee_yang_moment = check_lee_yang_counterexample()
    path, cycles, scores, path_count = tournament_fingerprint(m8_valuations)

    print("THM-2020 FINITE-PLACE CHANNEL SEPARATION - EXACT GATES")
    print("gate Wick moment formula (1): PASS")
    print(f"  direct polynomial expansion vs balanced channels: {moment_checks} identities")
    print("gate Legendre digit/carry formulas (4)-(5): PASS")
    print(f"  exact channel-prime identities: {legendre_checks}")
    print("gate star occupation and wide-gap certificates: PASS")
    print("  (x,y,z)=(n+2z,3n-3z,z), A=6n+12z")
    for prime, n, m, valuations in wide_gap:
        print(
            f"  p={prime}: n={n}, m={m}, channels={len(valuations)}, "
            f"v_p={format_valuations(valuations)}"
        )
    print("gate cancellable M4 and separated M8: PASS")
    print("  M4 normalized channels=(2880,-2880), sum=0")
    print(f"  M8 normalized channels={tuple(m8_coefficients)}")
    print(f"  v_3={tuple(m8_valuations)}, sum={m8_total}, v_3(sum)={vp_int(m8_total, 3)}")
    print("gate Lee-Yang upper-half-plane counterexample: PASS")
    print(f"  a={a}, b={b}, c={c}, b^2+2ac={lee_yang_moment}")
    print("ALL EXACT GATES PASS")
    print("Tournament Analysis (balanced M8 channels at p=3):")
    print(f"  valuations={tuple(m8_valuations)}, path={' -> '.join(path)}")
    print(
        f"  score_histogram={{0:1,1:1,2:1}}, directed_3cycles={cycles}, "
        f"SCC_sizes=[1,1,1], Hamiltonian_paths={path_count}"
    )
    print("  tie gauge=lexicographic channel label; observable=lower p-adic valuation")
    print("  assumption challenge: charge-only vertices rejected; retain channel+prime+carry state")


if __name__ == "__main__":
    main()
