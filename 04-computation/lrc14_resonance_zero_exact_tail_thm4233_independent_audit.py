#!/usr/bin/env python3
"""Independent ordered-interval audit of the THM-4233 exact tail.

Unlike the primary audit, this program never builds or classifies the union
of endpoint cells.  It constructs the two ordered families of safe teeth,
intersects them by a two-pointer scan, and evaluates the primitive only at
the resulting component starts and ends.

Reproduction:
  python3 04-computation/lrc14_resonance_zero_exact_tail_thm4233_independent_audit.py
"""

from __future__ import annotations

from fractions import Fraction
from math import gcd


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
FNV_OFFSET = 0xCBF29CE484222325
FNV_PRIME = 0x100000001B3
MASK64 = (1 << 64) - 1


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def fnv_word(state: int, value: int) -> int:
    value &= MASK64
    for shift in range(0, 64, 8):
        state ^= (value >> shift) & 0xFF
        state = (state * FNV_PRIME) & MASK64
    return state


def quadratic(coefficients: tuple[int, int, int], q: int) -> int:
    return coefficients[0] + q * (coefficients[1] + q * coefficients[2])


def safe_teeth(speed: int, other: int) -> list[tuple[int, int]]:
    return [
        (other * (14 * index + 1), other * (14 * index + 13))
        for index in range(speed)
    ]


def intersect_ordered(
    left: list[tuple[int, int]], right: list[tuple[int, int]]
) -> list[tuple[int, int]]:
    answer: list[tuple[int, int]] = []
    i = 0
    j = 0
    while i < len(left) and j < len(right):
        lo = max(left[i][0], right[j][0])
        hi = min(left[i][1], right[j][1])
        if lo < hi:
            if answer and answer[-1][1] == lo:
                answer[-1] = (answer[-1][0], hi)
            else:
                answer.append((lo, hi))
        if left[i][1] < right[j][1]:
            i += 1
        elif right[j][1] < left[i][1]:
            j += 1
        else:
            i += 1
            j += 1
    return answer


def audit(k: int) -> tuple[int, int, int, int, int, Fraction, int]:
    u = 5 * k + 2
    v = 7 * k + 3
    require(gcd(u, v) == 1 and 5 * v - 7 * u == 1, "family identity changed")
    grid = 14 * u * v
    components = intersect_ordered(safe_teeth(u, v), safe_teeth(v, u))
    safe_ticks = sum(hi - lo for lo, hi in components)

    covered = 0
    minimum = 0
    maximum = 0
    minimum_tick = 0
    maximum_tick = 0
    for lo, hi in components:
        at_start = grid * covered - safe_ticks * lo
        if at_start < minimum:
            minimum = at_start
            minimum_tick = lo
        covered += hi - lo
        at_end = grid * covered - safe_ticks * hi
        if at_end > maximum:
            maximum = at_end
            maximum_tick = hi
    require(covered == safe_ticks, "component integration changed safe mass")

    q, residue = divmod(k, 7)
    r_value = quadratic(R_COEFFICIENTS[residue], q)
    amplitude = 2 * u * r_value
    predicted_minimum_tick = (
        u * (98 * q - 13) if residue <= 1 else u * (98 * q + 85)
    ) % grid
    predicted_maximum_tick = (grid - predicted_minimum_tick) % grid
    predicted_safe_ticks = (72 * u * v + DENSITY_CORRECTION[residue]) // 7
    predicted_omega = Fraction(r_value, 49 * u * v * v)
    require(safe_ticks == predicted_safe_ticks, f"density mismatch at k={k}")
    require((minimum, maximum) == (-amplitude, amplitude),
            f"primitive amplitude mismatch at k={k}")
    require((minimum_tick, maximum_tick)
            == (predicted_minimum_tick, predicted_maximum_tick),
            f"extremizer mismatch at k={k}")
    omega = Fraction(maximum - minimum, grid * grid)
    require(omega == predicted_omega, f"omega mismatch at k={k}")
    return (grid, safe_ticks, amplitude, minimum_tick, maximum_tick, omega,
            len(components))


def main() -> None:
    rows = list(range(1, 851)) + [1000, 5000, 10000, 50000]
    ledger = FNV_OFFSET
    retained: dict[int, tuple[int, int, int, int, int, Fraction, int]] = {}
    for k in rows:
        result = audit(k)
        for value in result[:5]:
            ledger = fnv_word(ledger, int(value))
        ledger = fnv_word(ledger, result[5].numerator)
        ledger = fnv_word(ledger, result[5].denominator)
        ledger = fnv_word(ledger, result[6])
        if k in (1, 2, 747, 748, 850, 1000, 5000, 10000, 50000):
            retained[k] = result

    omega_747 = retained[747][5]
    omega_748 = retained[748][5]
    require(omega_747 > THRESHOLD and omega_748 < THRESHOLD,
            "747/748 gate orientation changed")

    print("THM-4233 independent ordered-safe-teeth tail audit")
    print(f"formula_rows={len(rows)} contiguous=1..850 far=1000,5000,10000,50000")
    for k in (1, 2, 850, 1000, 5000, 10000, 50000):
        grid, safe_ticks, amplitude, lo, hi, omega, components = retained[k]
        print(f"k={k} grid={grid} safe_ticks={safe_ticks} components={components} "
              f"extrema_ticks={hi}/{lo} amplitude={amplitude} omega={omega}")
    print(f"k=747 omega_minus_threshold={omega_747 - THRESHOLD}")
    print(f"k=748 threshold_minus_omega={THRESHOLD - omega_748}")
    print(f"semantic_ledger_fnv1a64_le_words={ledger:016x}")
    print("result=PASS independent interval intersection; formula and crossing reproduced")


if __name__ == "__main__":
    main()
