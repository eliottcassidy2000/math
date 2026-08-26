#!/usr/bin/env python3
"""Exact common-grid audit for the sharp THM-4233 resonance-zero tail.

For u=5k+2 and v=7k+3, the proof in THM-4233 derives a seven-residue
closed formula for the centered primitive oscillation.  This audit does not
stand in for that infinite macrostate proof.  It independently enumerates
every endpoint cell for k=1,...,1000 and three far controls, checks the exact
formula (including the extremizer addresses), and verifies the polynomial
sign certificates used at the 747/748 transfer-gate transition.

Reproduction:
  python3 04-computation/lrc14_resonance_zero_exact_tail_thm4233.py
"""

from __future__ import annotations

from fractions import Fraction
from math import comb, gcd


THRESHOLD = Fraction(1650, 28710227)
DENSITY_CORRECTION = (16, 0, -16, -18, -6, 6, 18)
R_COEFFICIENTS = (
    (67, 2997, 25284),
    (990, 10011, 25284),
    (2924, 17368, 25284),
    (6045, 24823, 25284),
    (10136, 32033, 25284),
    (15231, 39243, 25284),
    (21330, 46453, 25284),
)
THRESHOLD_SHIFT = (
    (17056167595726, 78191950749129, 1456258558482, 6794229750),
    (28901911461759, 78614469743901, 1459170371232, 6794229750),
    (39754122175768, 79027973077312, 1462082183982, 6794229750),
    (50360504188106, 79439494754977, 1464993996732, 6794229750),
    (61781100508659, 79858882384757, 1467905809482, 6794229750),
    (73262415238404, 80279101961037, 1470817622232, 6794229750),
    (5571348948256, 77773077303103, 1453346745732, 6794229750),
)
THRESHOLD_PREDECESSOR = (
    -59686318824671,
    -48260182140660,
    -37818562947312,
    -27620790799889,
    -16616670296366,
    -5552663330151,
    -70755175838865,
)
FNV_OFFSET = 0xCBF29CE484222325
FNV_PRIME = 0x100000001B3
MASK64 = (1 << 64) - 1


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def poly_add(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    size = max(len(left), len(right))
    answer = [0] * size
    for index in range(size):
        answer[index] = (
            (left[index] if index < len(left) else 0)
            + (right[index] if index < len(right) else 0)
        )
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return tuple(answer)


def poly_scale(poly: tuple[int, ...], scalar: int) -> tuple[int, ...]:
    return tuple(scalar * coefficient for coefficient in poly)


def poly_mul(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    answer = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            answer[i + j] += a * b
    return tuple(answer)


def poly_shift(poly: tuple[int, ...], offset: int) -> tuple[int, ...]:
    """Return coefficients of poly(n+offset), low degree first."""
    answer = (0,)
    for degree, coefficient in enumerate(poly):
        term = [0] * (degree + 1)
        for power in range(degree + 1):
            term[power] = coefficient * comb(degree, power) * offset ** (degree - power)
        answer = poly_add(answer, tuple(term))
    return answer


def poly_value(poly: tuple[int, ...], value: int) -> int:
    answer = 0
    for coefficient in reversed(poly):
        answer = answer * value + coefficient
    return answer


def fnv_word(state: int, value: int) -> int:
    value &= MASK64
    for shift in range(0, 64, 8):
        state ^= (value >> shift) & 0xFF
        state = (state * FNV_PRIME) & MASK64
    return state


def expected(k: int) -> tuple[int, int, int, int, int, Fraction]:
    q, residue = divmod(k, 7)
    u = 5 * k + 2
    v = 7 * k + 3
    grid = 14 * u * v
    r_value = poly_value(R_COEFFICIENTS[residue], q)
    amplitude = 2 * u * r_value
    if residue <= 1:
        minimum_tick = u * (98 * q - 13) % grid
    else:
        minimum_tick = u * (98 * q + 85) % grid
    maximum_tick = (grid - minimum_tick) % grid
    safe_ticks = (72 * u * v + DENSITY_CORRECTION[residue]) // 7
    omega = Fraction(r_value, 49 * u * v * v)
    return grid, safe_ticks, amplitude, minimum_tick, maximum_tick, omega


def midpoint_safe(speed: int, grid: int, left: int, right: int) -> bool:
    residue = speed * (left + right) % (2 * grid)
    return grid <= 7 * residue <= 13 * grid


def exact_pair(k: int) -> tuple[int, int, int, int, int, Fraction]:
    u = 5 * k + 2
    v = 7 * k + 3
    require(gcd(u, v) == 1 and 5 * v - 7 * u == 1, "unimodular family changed")
    grid = 14 * u * v
    points = [0, grid]
    points.extend((v * (14 * index + sign)) % grid
                  for index in range(u) for sign in (-1, 1))
    points.extend((u * (14 * index + sign)) % grid
                  for index in range(v) for sign in (-1, 1))
    points = sorted(set(points))

    cells: list[tuple[int, int, bool]] = []
    safe_ticks = 0
    for left, right in zip(points, points[1:]):
        safe = (midpoint_safe(u, grid, left, right)
                and midpoint_safe(v, grid, left, right))
        cells.append((left, right, safe))
        if safe:
            safe_ticks += right - left

    primitive = 0
    minimum = 0
    maximum = 0
    minimum_tick = 0
    maximum_tick = 0
    for left, right, safe in cells:
        primitive += (right - left) * (grid * int(safe) - safe_ticks)
        if primitive < minimum:
            minimum = primitive
            minimum_tick = right
        if primitive > maximum:
            maximum = primitive
            maximum_tick = right
    require(primitive == 0, "centered primitive failed to close")
    return (grid, safe_ticks, maximum, minimum_tick, maximum_tick,
            Fraction(maximum - minimum, grid * grid))


def threshold_polynomial(residue: int) -> tuple[int, ...]:
    u = (5 * residue + 2, 35)
    v = (7 * residue + 3, 49)
    positive = poly_scale(poly_mul(u, poly_mul(v, v)), 1650 * 49)
    negative = poly_scale(R_COEFFICIENTS[residue], -28710227)
    return poly_add(positive, negative)


def next_residue_data(residue: int) -> tuple[tuple[int, ...], tuple[int, ...], tuple[int, ...]]:
    if residue < 6:
        next_r = residue + 1
        return (
            R_COEFFICIENTS[next_r],
            (5 * next_r + 2, 35),
            (7 * next_r + 3, 49),
        )
    return (
        poly_shift(R_COEFFICIENTS[0], 1),
        (37, 35),
        (52, 49),
    )


def monotonicity_numerator(residue: int) -> tuple[int, ...]:
    u = (5 * residue + 2, 35)
    v = (7 * residue + 3, 49)
    next_r, next_u, next_v = next_residue_data(residue)
    left = poly_mul(R_COEFFICIENTS[residue],
                    poly_mul(next_u, poly_mul(next_v, next_v)))
    right = poly_mul(next_r, poly_mul(u, poly_mul(v, v)))
    return poly_add(left, poly_scale(right, -1))


def main() -> None:
    rows = list(range(1, 1001)) + [5000, 10000, 50000]
    ledger = FNV_OFFSET
    retained: dict[int, tuple[int, int, int, int, int, Fraction]] = {}
    for k in rows:
        observed = exact_pair(k)
        predicted = expected(k)
        require(observed == predicted, f"closed formula mismatch at k={k}")
        for value in observed[:5]:
            ledger = fnv_word(ledger, int(value))
        ledger = fnv_word(ledger, observed[5].numerator)
        ledger = fnv_word(ledger, observed[5].denominator)
        if k <= 6 or k in (747, 748, 1000, 5000, 10000, 50000):
            retained[k] = observed

    for residue in range(7):
        q0 = 106 if residue == 6 else 107
        shifted = poly_shift(threshold_polynomial(residue), q0)
        require(shifted == THRESHOLD_SHIFT[residue],
                f"threshold polynomial changed in residue {residue}")
        require(all(coefficient > 0 for coefficient in shifted),
                f"threshold positivity lost in residue {residue}")
        require(poly_value(threshold_polynomial(residue), q0 - 1)
                == THRESHOLD_PREDECESSOR[residue],
                f"threshold predecessor changed in residue {residue}")
        monotone = monotonicity_numerator(residue)
        require(all(coefficient > 0 for coefficient in monotone),
                f"consecutive monotonicity lost after residue {residue}")

    omega_747 = retained[747][5]
    omega_748 = retained[748][5]
    require(omega_747 == Fraction(4575651, 79563540224), "k=747 omega changed")
    require(omega_747 - THRESHOLD
            == Fraction(12591073311, 326326757250667264),
            "k=747 threshold margin changed")
    require(omega_748 == Fraction(144518186, 2516324606159), "k=748 omega changed")
    require(THRESHOLD - omega_748
            == Fraction(336393488, 8724097409553253),
            "k=748 threshold margin changed")

    print("THM-4233 resonance-zero exact common-grid audit")
    print(f"formula_rows={len(rows)} contiguous=1..1000 far=5000,10000,50000")
    for k in range(1, 7):
        grid, _, amplitude, minimum_tick, maximum_tick, omega = retained[k]
        print(f"k={k} grid={grid} extrema_ticks={maximum_tick}/{minimum_tick} "
              f"amplitude={amplitude} omega={omega}")
    print(f"k=747 omega={omega_747} omega_minus_threshold={omega_747 - THRESHOLD}")
    print(f"k=748 omega={omega_748} threshold_minus_omega={THRESHOLD - omega_748}")
    print("threshold_shift_coefficients=PASS all_positive residues=0..6")
    print("consecutive_monotonicity=PASS all_coefficients_positive k>=1")
    print(f"semantic_ledger_fnv1a64_le_words={ledger:016x}")
    print("result=PASS exact formula all k tested; exact gate iff k>=748 by analytic macrostate proof")


if __name__ == "__main__":
    main()
