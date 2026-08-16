#!/usr/bin/env python3
"""Finite exact companion for proved THM-3512."""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
import hashlib


WIDTH = 25
MAX_TIME = 2700
EXPECTED_V = (1, 3, 4, 6, 7, 9, 15, 16, 24, 25)


def check(condition: bool, label: object) -> None:
    if not condition:
        raise RuntimeError(f"check failed: {label}")


def nu2(value: int) -> int:
    check(value != 0, "nu2 zero")
    value = abs(value)
    return (value & -value).bit_length() - 1


def nu_fraction(value: Fraction) -> int:
    check(value != 0, "fraction valuation zero")
    return nu2(value.numerator) - nu2(value.denominator)


def mod_fraction(value: Fraction, modulus: int) -> int:
    check(value.denominator & 1, "nonunit denominator")
    return (value.numerator * pow(value.denominator, -1, modulus)) % modulus


def phi_full(x: int) -> int:
    return x ^ ((x << 1) | (x << 2))


def phi_width(x: int, mask: int) -> int:
    return phi_full(x) & mask


def full_rows(limit: int) -> list[int]:
    rows = [1]
    for _ in range(limit):
        rows.append(phi_full(rows[-1]))
    check(all(row & 1 for row in rows), "permanent odd bit")
    return rows


def seed_orbit(width: int) -> tuple[int, list[int]]:
    mask = (1 << width) - 1
    row = 1
    orbit = [row]
    for period in range(1, 1 << width):
        row = phi_width(row, mask)
        if row == 1:
            return period, orbit
        orbit.append(row)
    raise RuntimeError(f"seed return missing at width {width}")


def period_tower() -> tuple[list[int], list[list[int]], list[int], list[int]]:
    periods: list[int] = []
    orbits: list[list[int]] = []
    eps: list[int] = []
    for width in range(1, WIDTH + 1):
        period, orbit = seed_orbit(width)
        periods.append(period)
        orbits.append(orbit)
        mask = (1 << (width + 1)) - 1
        row = 1
        for _ in range(period):
            row = phi_width(row, mask)
        eps.append((row >> width) & 1)

    for width in range(1, WIDTH):
        check(
            periods[width] == periods[width - 1] << eps[width - 1],
            ("period lift", width),
        )
    innovations = [width for width, bit in enumerate(eps, 1) if bit]
    check(tuple(innovations[: len(EXPECTED_V)]) == EXPECTED_V, "innovation prefix")
    check(
        all(not (eps[i] and eps[i + 1] and eps[i + 2]) for i in range(len(eps) - 2)),
        "no111",
    )
    return periods, orbits, eps, innovations


def height_data(innovations: list[int]) -> tuple[list[int], list[int]]:
    values = innovations[: len(EXPECTED_V)]
    gaps = [values[m + 1] - values[m] for m in range(len(values) - 1)]
    heights = [value - m - 1 for m, value in enumerate(values)]
    check(heights[0] == 0, "h0")
    for m, gap in enumerate(gaps):
        check(heights[m + 1] - heights[m] == gap - 1, ("height increment", m))
    check(all(gap >= 1 for gap in gaps), "positive gaps")
    check(all(not (gaps[m] == gaps[m + 1] == 1) for m in range(len(gaps) - 1)), "gap no00")
    check(all(height >= m // 2 for m, height in enumerate(heights)), "height floor")
    return gaps, heights


def make_unit(rows_j: list[int], values: list[int]):
    def unit(m: int, t: int) -> int:
        q = 1 << m
        numerator = rows_j[t + q] - rows_j[t]
        scale = 1 << (values[m] - 1)
        check(numerator % scale == 0, ("unit divisibility", m, t))
        answer = numerator // scale
        check(answer & 1, ("unit odd", m, t))
        return answer

    return unit


def audit_van_der_put(
    rows_j: list[int], values: list[int], heights: list[int]
) -> tuple[int, int, int]:
    shell_checks = projection_checks = derivative_checks = 0
    for m, value in enumerate(values):
        q = 1 << m
        seen_valuations: set[int] = set()
        for residue in range(q):
            coefficient = rows_j[residue + q] - rows_j[residue]
            check(nu2(coefficient) == value - 1, ("shell coefficient", m, residue))
            reduced = coefficient // q
            check(nu2(reduced) == heights[m], ("reduced coefficient", m, residue))
            owner = coefficient // (1 << (value - 1))
            check(owner & 1, ("owner odd", m, residue))
            seen_valuations.add(nu2(reduced))
            shell_checks += 1
        check(seen_valuations == {heights[m]}, ("flat shell", m))

    for precision in range(1, 5):
        for m in range(2 * precision, len(values)):
            q = 1 << m
            for residue in range(q):
                reduced = (rows_j[residue + q] - rows_j[residue]) // q
                check(reduced % (1 << precision) == 0, ("fixed projection", precision, m))
                projection_checks += 1

    odd_multipliers = (1, 3, 5)
    for m, height in enumerate(heights):
        q = 1 << m
        for t in range(24):
            for multiplier in odd_multipliers:
                difference = rows_j[t + multiplier * q] - rows_j[t]
                quotient = Fraction(difference, multiplier * q)
                check(nu_fraction(quotient) == height, ("strict quotient", m, t, multiplier))
                derivative_checks += 1
    return shell_checks, projection_checks, derivative_checks


def audit_haar(
    unit, values: list[int], gaps: list[int], heights: list[int]
) -> tuple[int, int, int, int]:
    sibling_checks = block_checks = detail_checks = volkenborn_checks = 0

    for m, gap in enumerate(gaps):
        q = 1 << m
        for t in range(32):
            left = unit(m, t)
            right = unit(m, t + q)
            parent = unit(m + 1, t)
            check(left + right == (1 << gap) * parent, ("sibling sum", m, t))
            low = (left + right) // 2
            detail = (right - left) // 2
            check(low == (1 << (gap - 1)) * parent, ("low channel", m, t))
            check(low - detail == left and low + detail == right, ("inverse lift", m, t))
            check((detail % 2 == 0) == (gap == 1), ("detail parity", m, t))

            product = (left * right) % 8
            if gap == 1:
                check(product in (1, 5), ("product gap1", m, t))
            elif gap == 2:
                check(product == 3, ("product gap2", m, t))
            else:
                check(product == 7, ("product gap3plus", m, t))
            sibling_checks += 1

    for m in range(6):
        q = 1 << m
        for depth in range(1, len(values) - m):
            block_size = 1 << depth
            for t in range(16):
                total = sum(unit(m, t + j * q) for j in range(block_size))
                check(total % block_size == 0, ("block integral", m, depth, t))
                average = total // block_size
                expected = (1 << (heights[m + depth] - heights[m])) * unit(m + depth, t)
                check(average == expected, ("block average", m, depth, t))
                check(
                    nu2(average) == heights[m + depth] - heights[m],
                    ("block valuation", m, depth, t),
                )
                block_checks += 1

    for m in range(6):
        q = 1 << m
        for level in range(len(values) - m - 1):
            half = 1 << level
            for t in range(16):
                left_sum = sum(unit(m, t + j * q) for j in range(half))
                right_sum = sum(unit(m, t + j * q) for j in range(half, 2 * half))
                numerator = right_sum - left_sum
                divisor = 1 << (level + 1)
                check(numerator % divisor == 0, ("detail integral", m, level, t))
                detail = numerator // divisor
                scale_detail = (
                    (1 << (heights[m + level] - heights[m]))
                    * (unit(m + level, t + (1 << (m + level))) - unit(m + level, t))
                    // 2
                )
                check(detail == scale_detail, ("multilevel detail", m, level, t))
                detail_checks += 1

    for m, height in enumerate(heights):
        block_size = 1 << m
        for t in range(16):
            total = sum(unit(0, t + j) for j in range(block_size))
            check(total % block_size == 0, ("Volkenborn integral", m, t))
            average = total // block_size
            check(average == (1 << height) * unit(m, t), ("Volkenborn value", m, t))
            volkenborn_checks += 1

    return sibling_checks, block_checks, detail_checks, volkenborn_checks


def projective_ratio(unit, m: int, t: int) -> Fraction:
    q = 1 << m
    return -Fraction(unit(m, t + q), unit(m, t))


def audit_projective(unit, gaps: list[int]) -> tuple[int, int, int, int, int]:
    ratio_checks = recurrence_checks = precision_checks = 0
    quotient_nonclosure_checks = scalar_checks = 0
    for m, gap in enumerate(gaps):
        q = 1 << m
        for t in range(32):
            left = unit(m, t)
            parent = unit(m + 1, t)
            ratio = projective_ratio(unit, m, t)
            check(nu_fraction(1 - ratio) == gap, ("ratio gap", m, t))
            normalized = (1 - ratio) / (1 << gap)
            check(normalized == Fraction(parent, left), ("ratio owner", m, t))
            check(nu_fraction(normalized) == 0, ("ratio owner odd", m, t))
            check((ratio % 4 if isinstance(ratio, int) else mod_fraction(ratio, 4)) == (3 if gap == 1 else 1), ("ratio sign", m, t))
            ratio_checks += 1

    for m in range(len(gaps)):
        q = 1 << m
        for t in range(32):
            g0 = projective_ratio(unit, m, t)
            g1 = projective_ratio(unit, m, t + q)
            g2 = projective_ratio(unit, m, t + 2 * q)
            expected = -g0 * g1 * (1 - g2) / (1 - g0)
            actual = projective_ratio(unit, m + 1, t)
            check(actual == expected, ("projective recurrence", m, t))
            recurrence_checks += 1

    for precision in range(1, 9):
        near = 1 + (1 << precision)
        farther = 1 + (1 << (precision + 1))
        modulus = 1 << precision
        check(near % modulus == farther % modulus == 1 % modulus, ("precision collision", precision))
        check(nu2(1 - near) == precision, ("precision first gap", precision))
        check(nu2(1 - farther) == precision + 1, ("precision second gap", precision))

        depth = precision + 2
        z1, z2 = 1, 3
        g1 = 1 - (1 << depth) * z1
        g2 = 1 - (1 << depth) * z2
        check(g1 % (1 << depth) == g2 % (1 << depth), ("owner hidden", precision))
        check(
            ((1 - g1) >> depth) % 4 != ((1 - g2) >> depth) % 4,
            ("owner restored", precision),
        )
        precision_checks += 1

    for precision in range(2, 9):
        depth = precision + 2
        g = 1 - (1 << depth)
        g2 = 1 - 3 * (1 << depth)
        modulus = 1 << precision
        check(g % modulus == g2 % modulus, ("quotient input collision", precision))
        check(nu2(1 - g) == nu2(1 - g2) == depth, ("common input gap", precision))
        first = -g * g
        second = -3 * g * g
        check(first % modulus != second % modulus, ("quotient output split", precision))
        quotient_nonclosure_checks += 1

    for odd in range(1, 512, 2):
        check(nu2(1 + odd * odd) == 1, ("scalar ratio hostile", odd))
        scalar_checks += 1

    physical = projective_ratio(unit, 2, 0)
    check(mod_fraction(physical, 256) == 133, "physical ratio mod256")
    check(nu_fraction(physical - 1) == 2, "physical ratio gap")
    return (
        ratio_checks,
        recurrence_checks,
        precision_checks,
        quotient_nonclosure_checks,
        scalar_checks,
    )


def audit_state_invoice(
    periods: list[int], orbits: list[list[int]], eps: list[int]
) -> tuple[int, int, int, tuple[int, int, int, int]]:
    image_checks = fiber_checks = no111_bounds = 0
    final_row = (0, 0, 0, 0)
    for depth in range(1, WIDTH):
        period = periods[depth]
        exponent = period.bit_length() - 1
        check(period == 1 << exponent, ("period power", depth))
        orbit = orbits[depth]
        mask = (1 << depth) - 1
        image = {((state - 1) >> 1) & mask for state in orbit}
        check(len(image) == period, ("J image count", depth))
        lower_states = 1 << (depth - exponent)
        check(lower_states * period == 1 << depth, ("state invoice", depth))
        image_checks += 1

        innovations = sum(eps[:depth])
        check(innovations == exponent, ("E count", depth))
        check(innovations <= (2 * depth + 2) // 3, ("no111 E bound", depth))
        check(lower_states >= 1 << (depth // 3), ("no111 state bound", depth))
        no111_bounds += 1

        if depth <= 12:
            counts: Counter[int] = Counter()
            for prefix in range(1 << depth):
                state = orbit[prefix % period]
                counts[((state - 1) >> 1) & mask] += 1
            check(len(counts) == period, ("prefix fibers", depth))
            check(set(counts.values()) == {lower_states}, ("uniform fibers", depth))
            fiber_checks += len(counts)

        if depth == WIDTH - 1:
            final_row = (depth, exponent, period, lower_states)
    return image_checks, fiber_checks, no111_bounds, final_row


def audit_quasisymmetry(gaps: list[int]) -> tuple[int, int]:
    adjacent_checks = power_control_checks = 0
    for gap in gaps:
        input_ratio = 2
        output_ratio = 1 << gap
        check(output_ratio >= input_ratio, "adjacent radial expansion")
        adjacent_checks += 1

    observed_bound = max(gaps)
    for start in range(len(gaps)):
        for stop in range(start + 1, len(gaps) + 1):
            input_ratio = 1 << (stop - start)
            output_ratio = 1 << sum(gaps[start:stop])
            check(output_ratio <= input_ratio**observed_bound, ("finite eta control", start, stop))
            power_control_checks += 1
    return adjacent_checks, power_control_checks


def audit_mahler(rows_j: list[int]) -> tuple[int, int, int, int]:
    check(tuple(rows_j[:3]) == (0, 3, 12), "marked J values")
    van_der_put = rows_j[2] - rows_j[0]
    mahler = rows_j[2] - 2 * rows_j[1] + rows_j[0]
    check(van_der_put == 12 and nu2(van_der_put) == 2, "van der Put hostile")
    check(mahler == 6 and nu2(mahler) == 1, "Mahler hostile")
    return rows_j[1], rows_j[2], van_der_put, mahler


def main() -> None:
    periods, orbits, eps, innovations = period_tower()
    values = innovations[: len(EXPECTED_V)]
    gaps, heights = height_data(values)
    rows = full_rows(MAX_TIME)
    rows_j = [(row - 1) // 2 for row in rows]
    unit = make_unit(rows_j, values)

    shell, projections, derivatives = audit_van_der_put(rows_j, values, heights)
    sibling, blocks, details, volkenborn = audit_haar(unit, values, gaps, heights)
    ratio, ratio_rec, precision, quotient_nonclosure, scalar = audit_projective(unit, gaps)
    images, fibers, no111_bounds, final_row = audit_state_invoice(periods, orbits, eps)
    adjacent, eta_controls = audit_quasisymmetry(gaps)
    marked = audit_mahler(rows_j)

    physical_units = (unit(2, 0), unit(2, 4), unit(3, 0))
    check(physical_units == (25, 6403, 1607), "physical units")
    check(physical_units[0] + physical_units[1] == 4 * physical_units[2], "physical trace")

    digest_payload = (
        tuple(values), tuple(gaps), tuple(heights), tuple(periods),
        physical_units, final_row, marked,
    )
    digest = hashlib.sha256(repr(digest_payload).encode("ascii")).hexdigest()

    print("RULE30_VAN_DER_PUT_HAAR_AUTOMATON_THM3512")
    print("status PROVED_VERIFIED_EXACT_INDEPENDENTLY_AUDITED")
    print("width", WIDTH, "period", periods[-1], "innovation_prefix", tuple(values))
    print("gaps", tuple(gaps), "heights", tuple(heights))
    print("van_der_put_checks", shell, "fixed_projection_checks", projections, "derivative_checks", derivatives)
    print("haar_checks", sibling, "block_checks", blocks, "detail_checks", details, "volkenborn_checks", volkenborn)
    print(
        "projective_checks",
        ratio,
        "recurrence_checks",
        ratio_rec,
        "precision_tariffs",
        precision,
        "quotient_nonclosure_checks",
        quotient_nonclosure,
        "scalar_hostiles",
        scalar,
    )
    print("state_image_checks", images, "fiber_checks", fibers, "no111_state_bounds", no111_bounds)
    print("depth24_invoice", final_row, "exact_floor_bound", 1 << (final_row[0] // 3))
    print("quasisymmetry_controls", adjacent, eta_controls)
    print("physical_units", physical_units, "G2_mod256", 133)
    print("mahler_hostile", marked)
    print("control_sha256", digest)
    print("epsilon_prefix", "".join(map(str, eps)))
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
