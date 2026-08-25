#!/usr/bin/env python3
"""Fraction-exact hostile audit for THM-4088.

The computation checks three explicitly monotone rational approximation
sequences with rational, quadratic-irrational, and Liouville limits.  Their
labelled order tournaments agree identically.  It also checks the determinant
sidecar, the Pell recurrence used in the proof, the LRC Lipschitz inequality on
an exact hostile bank, and exact p-adic quality levels for rational controls.

Transcendence of the Liouville limit and the density/extension statements are
proved in the theorem; they are deliberately not inferred numerically here.
"""

from __future__ import annotations

from fractions import Fraction
from math import factorial


GATES = 0


def require(condition: bool, label: str) -> None:
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(f"FAILED: {label}")


def order_matrix(values: tuple[Fraction, ...]) -> tuple[tuple[int, ...], ...]:
    return tuple(
        tuple(0 if i == j else (1 if values[i] < values[j] else -1)
              for j in range(len(values)))
        for i in range(len(values))
    )


def rational_limit_sequence(length: int) -> tuple[Fraction, ...]:
    return tuple(Fraction(1, 2) - Fraction(1, n + 2) for n in range(1, length + 1))


def sqrt2_lower_sequence(length: int) -> tuple[Fraction, ...]:
    # a_k+b_k sqrt(2)=(1+sqrt(2))^(2k+1); multiplication by 3+2sqrt(2).
    a, b = 1, 1
    values: list[Fraction] = []
    for _ in range(length):
        require(a * a - 2 * b * b == -1, "negative Pell invariant")
        values.append(Fraction(a, b))
        next_a, next_b = 3 * a + 4 * b, 2 * a + 3 * b
        require(next_a * b - a * next_b == 2, "Pell consecutive determinant")
        a, b = next_a, next_b
    return tuple(values)


def liouville_truncations(length: int) -> tuple[Fraction, ...]:
    total = Fraction(0)
    values: list[Fraction] = []
    for n in range(1, length + 1):
        total += Fraction(1, 10 ** factorial(n))
        values.append(total)
        require(10 ** factorial(n) % total.denominator == 0,
                "Liouville denominator divides 10^(n!)")
    return tuple(values)


def determinant(left: Fraction, right: Fraction) -> int:
    """Numerator of right-left before reduction: a_j b_i-a_i b_j."""
    return right.numerator * left.denominator - left.numerator * right.denominator


def circle_norm(value: Fraction) -> Fraction:
    residue = value % 1
    return min(residue, 1 - residue)


def loneliness(speeds: tuple[int, ...], time: Fraction) -> Fraction:
    return min(circle_norm(speed * time) for speed in speeds)


def vp_fraction(value: Fraction, prime: int) -> int:
    require(value != 0, "finite p-adic valuation")
    numerator = abs(value.numerator)
    denominator = value.denominator
    answer = 0
    while numerator % prime == 0:
        numerator //= prime
        answer += 1
    while denominator % prime == 0:
        denominator //= prime
        answer -= 1
    return answer


def main() -> None:
    # n=7 already reaches the genuinely lacunary 10^(-5040) term while
    # keeping the exact hostile replay intentionally cheap.
    length = 7
    rational = rational_limit_sequence(length)
    quadratic = sqrt2_lower_sequence(length)
    transcendental = liouville_truncations(length)
    families = {
        "rational limit 1/2": rational,
        "quadratic limit sqrt(2)": quadratic,
        "Liouville limit": transcendental,
    }

    expected = tuple(
        tuple(0 if i == j else (1 if i < j else -1) for j in range(length))
        for i in range(length)
    )
    determinant_samples: dict[str, tuple[int, int, int]] = {}
    for name, values in families.items():
        require(all(values[i] < values[i + 1] for i in range(length - 1)),
                f"strict monotonicity: {name}")
        require(order_matrix(values) == expected,
                f"common index-order tournament prefix: {name}")
        for i in range(length):
            for j in range(i + 1, length):
                delta = determinant(values[i], values[j])
                require(delta > 0, f"positive determinant sign: {name},{i},{j}")
                require(values[j] - values[i]
                        == Fraction(delta, values[i].denominator * values[j].denominator),
                        f"determinant/height reconstruction: {name},{i},{j}")
        determinant_samples[name] = (
            determinant(values[0], values[1]),
            determinant(values[0], values[-1]).bit_length(),
            max(value.denominator for value in values).bit_length(),
        )

    require(Fraction(1, 2) - rational[-1] == Fraction(1, length + 2),
            "rational target error")
    for n in range(1, length):
        require(transcendental[n] - transcendental[n - 1]
                == Fraction(1, 10 ** factorial(n + 1)),
                "Liouville exact increment")

    # Exact Lipschitz audit: |min ||s t||-min ||s u||| <= max(s)|t-u|.
    speed_banks = ((1, 2, 5), (1, 3, 7, 13), (2, 9, 17, 31))
    times = tuple(Fraction(k, 97) for k in range(98))
    lrc_checks = 0
    for speeds in speed_banks:
        slope = max(speeds)
        for left in times:
            for right in times:
                require(abs(loneliness(speeds, left) - loneliness(speeds, right))
                        <= slope * abs(left - right),
                        f"LRC Lipschitz bank {speeds}")
                lrc_checks += 1

    # Rational p-adic controls realize any prescribed strict quality ranking.
    padic_checks = 0
    target = Fraction(3, 7)
    for prime in (2, 3, 5, 7, 11):
        approximants = tuple(target + prime ** height for height in range(1, 10))
        qualities = tuple(vp_fraction(target - value, prime) for value in approximants)
        require(qualities == tuple(range(1, 10)), f"p-adic exact quality ladder p={prime}")
        padic_checks += len(qualities)

    print("THM-4088 ORDER-TOURNAMENT ARITHMETIC-TYPE BLINDNESS AUDIT")
    print(f"prefix length: {length}")
    print("common labelled orientation: i -> j iff i < j")
    for name, sample in determinant_samples.items():
        print(
            f"{name}: first_det={sample[0]}, "
            f"long_det_bits={sample[1]}, max_den_bits={sample[2]}"
        )
    print(f"exact LRC Lipschitz pair checks: {lrc_checks}")
    print(f"exact p-adic quality checks: {padic_checks}")
    print("transcendence and density: PROOF-ONLY (not inferred by computation)")
    print(f"GATES: {GATES}")
    print("VERIFIED-EXACT: True")


if __name__ == "__main__":
    main()
