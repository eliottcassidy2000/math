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


C = Fraction(-29, 16)


def f(x: Fraction) -> Fraction:
    return x * x + C


def scaled_T(u: int, modulus: int | None = None) -> int:
    """T(u)=(u^2-29)/4, exactly on odd u or modulo an odd modulus."""
    if modulus is None:
        numerator = u * u - 29
        assert numerator % 4 == 0
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
    assert set(image.values()) <= set(candidates)

    cycle = [Fraction(-7, 4), Fraction(5, 4), Fraction(-1, 4)]
    assert [f(x) for x in cycle] == cycle[1:] + cycle[:1]

    print("COMPLETE RATIONAL PREPERIODIC GRAPH (proof filters stated above)")
    for x in candidates:
        print(f"  {x!s:>4} -> {image[x]}")
    print("  unique rational cycle: -7/4 -> 5/4 -> -1/4 -> -7/4")
    print("  rational exact period 6: NONE")


def moduli_and_triangle() -> None:
    """Audit the standard rational 3-cycle parameter at its AP values."""

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
        assert [local_f(x) for x in points] == points[1:] + points[:1]
        assert set(points) == expected
        assert parameter == C
        print(f"  t={t}: points={tuple(map(str, points))}, c={parameter}")

    # At the AP representative t=1, Euclid's pair (t+1,t)=(2,1)
    # gives the primitive 3-4-5 triangle.  The cycle denominator is its
    # even leg, while the c denominator is the square of that leg.
    t = 1
    m, n = t + 1, t
    odd_leg = m * m - n * n
    even_leg = 2 * m * n
    hypotenuse = m * m + n * n
    assert (odd_leg, even_leg, hypotenuse) == (3, 4, 5)
    assert 16 == even_leg**2
    assert 29 == hypotenuse**2 + m**2
    assert 7 == even_leg**2 - odd_leg**2
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
    assert orbit == [Fraction(1), Fraction(-1, 2), Fraction(-2)]
    assert R(orbit[-1]) == orbit[0]

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
    assert powers[3] == ((-1, 0), (0, -1))
    assert powers[6] == identity
    print(f"  relabeling R(t)=-1/(t+1) has AP orbit {tuple(map(str, orbit))}")
    print("  SL2 lift [[0,-1],[1,1]]: A^3=-I, A^6=I, charpoly=lambda^2-lambda+1")


def finite_ring_controls() -> None:
    print("\nFINITE-RING CONTROLS FOR T(u)=(u^2-29)/4")
    for modulus in (7, 9, 63):
        cycles = finite_cycles(modulus)
        print(f"  mod {modulus}: cycles={cycles}")
    cycles_63 = finite_cycles(63)
    assert cycles_63 == [(5, 62, 56), (14, 26, 20), (35, 47, 41)]
    assert all(len(cycle) == 3 for cycle in cycles_63)
    print("  hostile control: mod 63 has no exact 6-cycle")

    # The rational 3-cycle multiplier is
    # (-7/2)*(5/2)*(-1/2) = 35/8, hence -1 modulo the unique good prime 43.
    p = 43
    multiplier = Fraction(35, 8)
    assert (35 * pow(8, -1, p)) % p == p - 1
    modulus = p * p
    cycles = finite_cycles(modulus)
    counts = Counter(map(len, cycles))
    assert counts[3] == 1
    assert counts[6] == 21

    # The entire tube above {-7,5,-1} mod 43 contains 3*43 points.  It is
    # exactly the old 3-cycle plus the 21 six-cycles (3 + 21*6 = 129).
    base = {(-7) % p, 5 % p, (-1) % p}
    tube_cycles = [c for c in cycles if all(u % p in base for u in c)]
    assert Counter(map(len, tube_cycles)) == Counter({6: 21, 3: 1})
    assert sum(map(len, tube_cycles)) == 3 * p

    representative = [36]
    for _ in range(5):
        representative.append(scaled_T(representative[-1], modulus))
    assert scaled_T(representative[-1], modulus) == representative[0]
    assert len(set(representative)) == 6
    print(f"  3-cycle multiplier = {multiplier} == -1 (mod 43)")
    print(f"  mod 43^2 tube: {dict(sorted(Counter(map(len, tube_cycles)).items()))}")
    print(f"  representative exact 6-cycle mod 1849: {representative}")
    print("  63 local F=T^3 two-cycles = 3*(43-1)/2; T groups them into 21 six-cycles")


if __name__ == "__main__":
    rational_graph()
    moduli_and_triangle()
    finite_ring_controls()
