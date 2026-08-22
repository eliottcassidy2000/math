#!/usr/bin/env python3
"""Exact companion for THM-3501.

All arithmetic is over F_2.  A polynomial in q is stored as a Python integer:
bit v is the coefficient of q^v.  The verifier uses explicit ``check`` gates,
so ``python`` and ``python -O`` execute the same mathematical checks.
"""

from __future__ import annotations

from functools import lru_cache


AUDIT_DEPTHS = tuple(range(5, 17))
MAX_TIME = 320


def check(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def parity(poly: int) -> int:
    return poly.bit_count() & 1


def degree(poly: int) -> int:
    return poly.bit_length() - 1


def coeff(poly: int, exponent: int) -> int:
    return (poly >> exponent) & 1


def qpoly(exponents: tuple[int, ...] | list[int] | set[int]) -> int:
    out = 0
    for exponent in exponents:
        check(exponent >= 0, f"negative q exponent {exponent}")
        out ^= 1 << exponent
    return out


def format_q(poly: int) -> str:
    if poly == 0:
        return "0"
    terms: list[str] = []
    for exponent in range(poly.bit_length()):
        if not coeff(poly, exponent):
            continue
        if exponent == 0:
            terms.append("1")
        elif exponent == 1:
            terms.append("q")
        else:
            terms.append(f"q^{exponent}")
    return "+".join(terms)


def gf2_multiply(left: int, right: int) -> int:
    out = 0
    while right:
        if right & 1:
            out ^= left
        left <<= 1
        right >>= 1
    return out


def multiply_q_plus_one(poly: int) -> int:
    return poly ^ (poly << 1)


def divide_q_plus_one(poly: int) -> int:
    """Exact division by q+1, guarded by evaluation at q=1."""
    check(parity(poly) == 0, f"nondivisible polynomial {format_q(poly)}")
    remainder = poly
    quotient = 0
    while degree(remainder) >= 1:
        top = degree(remainder)
        quotient ^= 1 << (top - 1)
        remainder ^= (1 << top) | (1 << (top - 1))
    check(remainder == 0, f"q+1 division left remainder {remainder}")
    return quotient


@lru_cache(maxsize=None)
def trinomial_power(power: int) -> int:
    """Return (1+x+x^2)^power over F_2 as a bit polynomial."""
    out = 1
    base = 0b111
    exponent = power
    while exponent:
        if exponent & 1:
            out = gf2_multiply(out, base)
        base = gf2_multiply(base, base)
        exponent >>= 1
    return out


@lru_cache(maxsize=None)
def green(distance: int, time: int) -> int:
    """H_distance(time)=[x^(time+distance)](1+x+x^2)^time mod 2."""
    if distance < 0 or time < 0 or distance > time:
        return 0
    return (trinomial_power(time) >> (time + distance)) & 1


def green_recurrence_table(limit: int) -> list[list[int]]:
    table = [[0] * (limit + 2) for _ in range(limit + 1)]
    table[0][0] = 1
    for time in range(limit):
        for distance in range(time + 2):
            if distance == 0:
                # The two reflected distance-one contributions coincide and
                # cancel over F_2.
                value = table[time][0]
            else:
                value = (
                    table[time][distance - 1]
                    ^ table[time][distance]
                    ^ table[time][distance + 1]
                )
            table[time + 1][distance] = value
    return table


def build_centered_rows(limit: int) -> list[set[int]]:
    rows: list[set[int]] = [{0}]
    for time in range(limit):
        old = rows[-1]
        new: set[int] = set()
        for position in range(-time - 1, time + 2):
            left = int(position - 1 in old)
            center = int(position in old)
            right = int(position + 1 in old)
            value = left ^ center ^ right ^ (center & right)
            if value:
                new.add(position)
        rows.append(new)
    return rows


ROWS = build_centered_rows(MAX_TIME)


def cell(time: int, position: int) -> int:
    if time < 0 or time >= len(ROWS):
        return 0
    return int(position in ROWS[time])


def bond(time: int, position: int) -> int:
    return cell(time, position) & cell(time, position + 1)


def build_packed_rows(limit: int) -> list[int]:
    packed = [1]
    for _ in range(limit):
        state = packed[-1]
        packed.append(state ^ ((state << 1) | (state << 2)))
    return packed


PACKED = build_packed_rows(MAX_TIME)


def edge_bit(offset: int, time: int) -> int:
    return cell(time, time - offset)


def seed_period(depth: int) -> int:
    target = tuple([1] + [0] * (depth - 1))
    for period in range(1, MAX_TIME - depth):
        state = tuple(edge_bit(offset, period) for offset in range(depth))
        if state == target:
            return period
    raise RuntimeError(f"period search failed at depth {depth}")


def raw_potential(depth: int, phase: int) -> int:
    """The universal-cover Green potential E_depth(phase;q)."""
    target_time = depth + phase
    check(target_time >= 0, f"raw potential below physical cut: k={depth}, h={phase}")
    out = 0

    initial_distance = abs(phase)
    if green(initial_distance, target_time):
        initial_slack = target_time - initial_distance
        check(initial_slack >= 0, "live initial atom has negative slack")
        out ^= 1 << initial_slack

    for source_time in range(target_time):
        horizon = target_time - 1 - source_time
        for position in range(phase - horizon, phase + horizon + 1):
            distance = abs(phase - position)
            if bond(source_time, position) and green(distance, horizon):
                slack = horizon - distance
                check(slack >= 0, "live collision atom has negative slack")
                out ^= 1 << slack
    return out


def terminal_value(depth: int, phase: int) -> int:
    return cell(depth + phase, phase)


def terminal_current(depth: int, phase: int) -> int:
    time = depth + phase
    return cell(time, phase + 1) | cell(time, phase + 2)


def folded_source(source_depth: int, distance: int) -> int:
    source_time = source_depth + distance
    if source_time < 2:
        return 0
    if distance == 0:
        return bond(source_time, 0)
    return bond(source_time, distance) ^ bond(source_time, -distance)


def strip_sum(depth: int, first_source_depth: int, last_source_depth: int) -> int:
    out = 0
    for source_depth in range(first_source_depth, last_source_depth + 1):
        for distance in range(depth + 1):
            slack = depth - source_depth - 1 - 2 * distance
            if slack < 0:
                continue
            if folded_source(source_depth, distance) and green(distance, distance + slack):
                out ^= 1 << slack
    return out


def inward_polynomial(depth: int) -> int:
    backbone = green(1, depth - 2)
    return backbone ^ strip_sum(depth, 3, depth - 1)


def time_two_full_polynomial(depth: int) -> int:
    backbone = green(1, depth - 2)
    return backbone ^ strip_sum(depth, 0, depth - 1)


def boundary_strip_polynomial(depth: int) -> int:
    return strip_sum(depth, 0, 2)


def contract_to_degree(poly: int, cutoff: int) -> int:
    """Send every q^m with m>cutoff to the pointed grade 1=q^0."""
    low_mask = (1 << (cutoff + 1)) - 1
    low = poly & low_mask
    high = poly >> (cutoff + 1)
    if parity(high):
        low ^= 1
    return low


def contraction_homotopy(poly: int, cutoff: int) -> int:
    """H_N with id+C_N=(q+1)H_N."""
    out = 0
    for exponent in range(cutoff + 1, poly.bit_length()):
        if coeff(poly, exponent):
            out ^= (1 << exponent) - 1
    return out


def hasse_moment(profiles: list[int], order: int) -> int:
    """Phase Hasse moment, still as a q-polynomial."""
    out = 0
    for phase, profile in enumerate(profiles):
        if (order & ~phase) == 0:
            out ^= profile
    return out


def binomial_transform(bits: list[int]) -> list[int]:
    """Change X coefficients to Y=X+1 coefficients (an involution)."""
    out = [0] * len(bits)
    for exponent, bit in enumerate(bits):
        if not bit:
            continue
        submask = exponent
        while True:
            out[submask] ^= 1
            if submask == 0:
                break
            submask = (submask - 1) & exponent
    return out


def invert_unit_series(unit: list[int]) -> list[int]:
    """Inverse of a unit modulo Y^len(unit)."""
    check(unit[0] == 1, "attempted to invert a nonunit Y-series")
    inverse = [0] * len(unit)
    inverse[0] = 1
    for exponent in range(1, len(unit)):
        value = 0
        for index in range(1, exponent + 1):
            value ^= unit[index] & inverse[exponent - index]
        inverse[exponent] = value
    return inverse


def arc_sum(current: list[int], length: int) -> list[int]:
    period = len(current)
    out: list[int] = []
    for phase in range(period):
        value = 0
        for offset in range(length):
            value ^= current[(phase + offset) % period]
        out.append(value)
    return out


def strict_suffix_mask(mask: int, length: int) -> int:
    """Strict suffix parity on low-to-high Mealy-factor coordinates."""
    value = mask
    shift = 1
    while shift < length:
        value ^= value >> shift
        shift <<= 1
    return value >> 1


def mealy_section_masks(length: int, stop: int) -> list[int]:
    ones = (1 << length) - 1
    masks = [ones, ones]
    while len(masks) <= stop:
        masks.append(
            strict_suffix_mask(masks[-2], length)
            | strict_suffix_mask(masks[-1], length)
        )
    return masks


def audit_foundations() -> None:
    recurrence = green_recurrence_table(MAX_TIME // 2)
    for time in range(len(recurrence)):
        for distance in range(time + 1):
            check(
                recurrence[time][distance] == green(distance, time),
                f"Green engines disagree at (d,n)=({distance},{time})",
            )

    for time in range(MAX_TIME + 1):
        state = PACKED[time]
        for offset in range(2 * time + 1):
            check(
                ((state >> offset) & 1) == edge_bit(offset, time),
                f"packed/centered disagreement at (r,t)=({offset},{time})",
            )

    # Independent Mealy-section compiler versus the literal edge current,
    # including the upper/lower cut as a bit-holonomy coboundary.
    for length in range(1, 33):
        masks = mealy_section_masks(length, AUDIT_DEPTHS[-1] - 1)
        for depth in range(2, AUDIT_DEPTHS[-1] + 1):
            profile = masks[depth - 1]
            for time in range(length):
                section_bit = (profile >> (length - 1 - time)) & 1
                edge_current = edge_bit(depth - 1, time) | edge_bit(depth - 2, time)
                check(
                    section_bit == edge_current,
                    f"section/current reversal fails at k={depth}, n={length}, t={time}",
                )

        if length & 1:
            continue
        half = length // 2
        half_mask = (1 << half) - 1
        for depth in range(2, AUDIT_DEPTHS[-1] + 1):
            profile = masks[depth - 1]
            defect = (profile & half_mask) ^ (profile >> half)
            defect_parity = 0
            for time in range(half):
                reversed_defect = (defect >> (half - 1 - time)) & 1
                current_coboundary = (
                    (edge_bit(depth - 1, time + half) | edge_bit(depth - 2, time + half))
                    ^ (edge_bit(depth - 1, time) | edge_bit(depth - 2, time))
                )
                bit_holonomy_coboundary = (
                    edge_bit(depth, time + half + 1)
                    ^ edge_bit(depth, time + 1)
                    ^ edge_bit(depth, time + half)
                    ^ edge_bit(depth, time)
                )
                check(
                    reversed_defect == current_coboundary == bit_holonomy_coboundary,
                    f"section/holonomy coboundary fails at k={depth}, n={length}, t={time}",
                )
                defect_parity ^= reversed_defect
            check(
                defect_parity == edge_bit(depth, length),
                f"section defect zeroth moment fails at k={depth}, n={length}",
            )


def audit_depth(depth: int) -> dict[str, object]:
    period = seed_period(depth)
    check(depth < period, f"declared hard-regime depth wraps: k={depth}, p={period}")
    cut = -depth
    values = [raw_potential(depth, cut + index) for index in range(period + depth)]

    for index, value in enumerate(values):
        phase = cut + index
        check(degree(value) <= depth, f"degree bound fails at k={depth}, h={phase}")
        check(
            parity(value) == terminal_value(depth, phase),
            f"Duhamel scalarization fails at k={depth}, h={phase}",
        )

    raw_center = values[depth]
    check(coeff(raw_center, depth) == 1, f"raw monicity fails at k={depth}")

    current = [values[index + 1] ^ values[index] for index in range(period)]
    for index, value in enumerate(current):
        phase = cut + index
        check(
            parity(value) == terminal_current(depth, phase),
            f"current orientation fails at k={depth}, h={phase}",
        )

    cyclic = arc_sum(current, depth)
    universal = [values[index + depth] ^ values[index] for index in range(period)]
    check(cyclic[0] == raw_center, f"marked no-wrap telescope fails at k={depth}")

    for index, value in enumerate(cyclic):
        phase = cut + index
        expected = 0
        for offset in range(depth):
            expected ^= terminal_current(depth, phase + offset)
        check(parity(value) == expected, f"phase profile fails at k={depth}, r={index}")

    epsilon = edge_bit(depth, period)
    holonomy = [values[period + offset] ^ values[offset] for offset in range(depth)]
    for offset, value in enumerate(holonomy):
        check(
            parity(value)
            == (edge_bit(depth, period + offset) ^ edge_bit(depth, offset)),
            f"slack/bit holonomy bridge fails at k={depth}, a={offset}",
        )
        check(parity(value) == epsilon, f"slack holonomy scalar fails at k={depth}, a={offset}")

    seam = [universal[index] ^ cyclic[index] for index in range(period)]
    for index, value in enumerate(seam):
        if index < period - depth:
            expected = 0
        else:
            offset = index - (period - depth)
            expected = holonomy[offset] ^ holonomy[0]
        check(value == expected, f"seam formula fails at k={depth}, r={index}")

    valuation = (depth & -depth).bit_length() - 1
    ell = (1 << valuation) - 1
    for order in range(ell):
        check(hasse_moment(cyclic, order) == 0, f"cyclic image fails at k={depth}, j={order}")
        omega_moment = 0
        for offset, value in enumerate(holonomy):
            if (order & ~offset) == 0:
                omega_moment ^= value
        seam_moment = hasse_moment(seam, order)
        check(
            seam_moment == omega_moment,
            f"Hasse seam tomography fails at k={depth}, j={order}",
        )
        check(
            hasse_moment(universal, order) == seam_moment,
            f"universal image obstruction fails at k={depth}, j={order}",
        )

    inward = inward_polynomial(depth)
    full = time_two_full_polynomial(depth)
    boundary = boundary_strip_polynomial(depth)
    check(full == (inward ^ boundary), f"inward/boundary split fails at k={depth}")

    backbone = green(1, depth - 2)
    restart = qpoly((depth, depth - 2))
    if backbone:
        restart ^= qpoly((depth - 3, 0))
    check(raw_center ^ full == restart, f"restart homotopy fails at k={depth}")
    normalization = raw_center ^ inward
    check(normalization == (restart ^ boundary), f"normalization K fails at k={depth}")
    check(parity(normalization) == 0, f"normalization does not vanish at q=1, k={depth}")
    boundary_hat = divide_q_plus_one(boundary)
    homotopy_l = qpoly((depth - 2, depth - 1)) ^ boundary_hat
    if backbone:
        homotopy_l ^= (1 << (depth - 3)) - 1
    check(
        multiply_q_plus_one(homotopy_l) == normalization,
        f"explicit normalization homotopy fails at k={depth}",
    )

    cutoff = depth - 4
    check(degree(inward) == cutoff, f"inward degree fails at k={depth}")
    check(coeff(inward, cutoff) == 1, f"inward monicity fails at k={depth}")
    check(degree(boundary) <= cutoff - 1, f"boundary degree fails at k={depth}")

    contracted_current = [contract_to_degree(value, cutoff) for value in current]
    for index, value in enumerate(current):
        check(
            value ^ contracted_current[index]
            == multiply_q_plus_one(contraction_homotopy(value, cutoff)),
            f"explicit contraction homotopy fails at k={depth}, r={index}",
        )
    contracted = arc_sum(contracted_current, depth)
    check(contracted[0] == full, f"contracted marked axis fails at k={depth}")
    check(
        contracted == [contract_to_degree(value, cutoff) for value in cyclic],
        f"contraction does not commute with arc at k={depth}",
    )

    top_carrier = [coeff(value, cutoff) for value in contracted]
    top_current = [coeff(value, cutoff) for value in contracted_current]
    check(
        arc_sum(top_current, depth) == top_carrier,
        f"top carrier/current preimage fails at k={depth}",
    )
    check(top_carrier[0] == 1, f"top carrier misses marked phase at k={depth}")
    top_as_profiles = [bit for bit in top_carrier]
    for order in range(ell):
        check(
            hasse_moment(top_as_profiles, order) == 0,
            f"top carrier leaves arc image at k={depth}, j={order}",
        )

    sharp = [
        value ^ (boundary if top_carrier[index] else 0)
        for index, value in enumerate(contracted)
    ]
    sharp_current = [
        value ^ (boundary if top_current[index] else 0)
        for index, value in enumerate(contracted_current)
    ]
    check(arc_sum(sharp_current, depth) == sharp, f"sharp current repair fails at k={depth}")
    check(sharp[0] == inward, f"sharp marked inward axis fails at k={depth}")
    for index, value in enumerate(sharp):
        check(degree(value) <= cutoff, f"sharp degree fails at k={depth}, r={index}")
        check(
            parity(value) == parity(cyclic[index]),
            f"sharp phase scalarization fails at k={depth}, r={index}",
        )
    for order in range(ell):
        check(hasse_moment(sharp, order) == 0, f"sharp arc image fails at k={depth}, j={order}")

    frobenius_carrier = [0] * period
    frobenius_carrier[0] = 1
    if ell:
        frobenius_carrier[1 << valuation] = 1
        check(sum(frobenius_carrier) == 2, f"Frobenius carrier support fails at k={depth}")
    else:
        check(sum(frobenius_carrier) == 1, f"odd carrier support fails at k={depth}")
    for order in range(ell):
        check(
            hasse_moment(frobenius_carrier, order) == 0,
            f"Frobenius carrier leaves arc image at k={depth}, j={order}",
        )

    # Verify the literal preimage from A_k(X)=Y^ell U_k(Y).  The actual
    # forward-arc multiplier is X^{-(k-1)}(1+...+X^{k-1}).
    arc_multiplier_x = [0] * period
    for offset in range(depth):
        arc_multiplier_x[(offset - (depth - 1)) % period] ^= 1
    arc_multiplier_y = binomial_transform(arc_multiplier_x)
    check(
        all(bit == 0 for bit in arc_multiplier_y[:ell]),
        f"arc multiplier valuation too small at k={depth}",
    )
    unit_y = arc_multiplier_y[ell:] + [0] * ell
    unit_y = unit_y[:period]
    inverse_unit_y = invert_unit_series(unit_y)
    if ell:
        frobenius_preimage_y = [0] + inverse_unit_y[:-1]
    else:
        frobenius_preimage_y = inverse_unit_y
    frobenius_preimage = binomial_transform(frobenius_preimage_y)
    check(
        arc_sum(frobenius_preimage, depth) == frobenius_carrier,
        f"explicit Frobenius current preimage fails at k={depth}",
    )
    sparse_sharp = [
        value ^ (boundary if frobenius_carrier[index] else 0)
        for index, value in enumerate(contracted)
    ]
    sparse_current = [
        value ^ (boundary if frobenius_preimage[index] else 0)
        for index, value in enumerate(contracted_current)
    ]
    check(
        arc_sum(sparse_current, depth) == sparse_sharp,
        f"sparse Frobenius current repair fails at k={depth}",
    )
    check(sparse_sharp[0] == inward, f"sparse marked inward axis fails at k={depth}")
    for index, value in enumerate(sparse_sharp):
        check(degree(value) <= cutoff, f"sparse degree fails at k={depth}, r={index}")
        check(
            parity(value) == parity(cyclic[index]),
            f"sparse phase scalarization fails at k={depth}, r={index}",
        )
    for order in range(ell):
        check(
            hasse_moment(sparse_sharp, order) == 0,
            f"sparse arc image fails at k={depth}, j={order}",
        )

    return {
        "depth": depth,
        "period": period,
        "ell": ell,
        "cutoff": cutoff,
        "raw": raw_center,
        "full": full,
        "inward": inward,
        "boundary": boundary,
        "restart": restart,
        "normalization": normalization,
        "holonomy": holonomy,
        "seam": seam,
        "cyclic": cyclic,
        "sharp": sharp,
        "sparse_sharp": sparse_sharp,
        "top_carrier_weight": sum(top_carrier),
        "frobenius_carrier_weight": sum(frobenius_carrier),
    }


def main() -> None:
    audit_foundations()
    records = [audit_depth(depth) for depth in AUDIT_DEPTHS]

    by_depth = {int(record["depth"]): record for record in records}
    five = by_depth[5]
    six = by_depth[6]

    check(int(five["period"]) == 8, "depth-five period changed")
    check(int(five["ell"]) == 0, "depth-five invertible control lost")
    check(int(five["cutoff"]) == 1, "depth-five cutoff changed")
    check(int(five["raw"]) == qpoly((1, 3, 5)), "depth-five raw polynomial changed")
    check(int(five["full"]) == qpoly((1,)), "depth-five time-two polynomial changed")
    check(int(five["inward"]) == qpoly((1,)), "depth-five inward polynomial changed")
    check(int(five["boundary"]) == 0, "depth-five boundary current changed")
    check(int(five["restart"]) == qpoly((3, 5)), "depth-five restart changed")
    check(int(five["normalization"]) == qpoly((3, 5)), "depth-five normalization changed")
    check(
        sum(1 for value in five["seam"] if int(value)) == 2,
        "depth-five seam support changed",
    )

    check(int(six["period"]) == 8, "depth-six period changed")
    check(int(six["ell"]) == 1, "depth-six image valuation changed")
    check(int(six["cutoff"]) == 2, "depth-six cutoff changed")
    check(int(six["raw"]) == qpoly((1, 2, 4, 6)), "depth-six raw polynomial changed")
    check(int(six["full"]) == qpoly((1, 2)), "depth-six time-two polynomial changed")
    check(int(six["inward"]) == qpoly((0, 2)), "depth-six inward polynomial changed")
    check(int(six["boundary"]) == qpoly((0, 1)), "depth-six boundary current changed")
    check(int(six["restart"]) == qpoly((4, 6)), "depth-six restart changed")
    check(
        int(six["normalization"]) == qpoly((0, 1, 4, 6)),
        "depth-six normalization changed",
    )
    six_seam_m0 = hasse_moment(list(six["seam"]), 0)
    check(six_seam_m0 == qpoly((2, 4)), "depth-six seam hostile changed")
    check(list(six["sharp"]) == [
        qpoly((0, 2)),
        0,
        qpoly((1,)),
        qpoly((0, 1)),
        qpoly((0,)),
        qpoly((1,)),
        qpoly((2,)),
        qpoly((0, 1)),
    ], "depth-six sharp table changed")
    check(list(six["sparse_sharp"]) == [
        qpoly((0, 2)),
        0,
        qpoly((0,)),
        qpoly((0, 1)),
        qpoly((0,)),
        qpoly((1,)),
        qpoly((0, 1, 2)),
        qpoly((0, 1)),
    ], "depth-six sparse table changed")
    six_cut_scalar = [parity(int(value)) for value in six["sparse_sharp"]]
    check(
        six_cut_scalar == [0, 0, 1, 0, 1, 1, 1, 0],
        "depth-six cut scalar profile changed",
    )
    six_standard_scalar = [0] * int(six["period"])
    for cut_index, value in enumerate(six_cut_scalar):
        six_standard_scalar[(-6 + cut_index) % int(six["period"])] = value
    check(
        six_standard_scalar == [1, 0, 0, 0, 1, 0, 1, 1],
        "depth-six standard phase profile changed",
    )

    print("THM-3501 exact companion: PASS")
    print(f"depths={AUDIT_DEPTHS[0]}..{AUDIT_DEPTHS[-1]} hard_unwrapped={len(records)}")
    print("independent_engines=centered_rule30/packed_Phi, trinomial_bits/Green_recurrence")
    print(
        "foundation_universe=rule30_times=0..{} green_times=0..{} "
        "section_lengths=1..32 section_depths=2..{}".format(
            MAX_TIME,
            MAX_TIME // 2,
            AUDIT_DEPTHS[-1],
        )
    )
    for record in records:
        print(
            f"depth={record['depth']} p={record['period']} ell={record['ell']} "
            f"N={record['cutoff']} top_carrier_weight={record['top_carrier_weight']} "
            f"frobenius_carrier_weight={record['frobenius_carrier_weight']} "
            f"raw={format_q(int(record['raw']))} "
            f"inward={format_q(int(record['inward']))} "
            f"G={format_q(int(record['boundary']))}"
        )
    print(
        "k=5 positive: raw={} inward={} nonzero_seam={} ell=0".format(
            format_q(int(five["raw"])),
            format_q(int(five["inward"])),
            sum(1 for value in five["seam"] if int(value)),
        )
    )
    print(
        "k=6 hostile: raw={} full={} inward={} Gamma={} G={} K={} M0(W)={}".format(
            format_q(int(six["raw"])),
            format_q(int(six["full"])),
            format_q(int(six["inward"])),
            format_q(int(six["restart"])),
            format_q(int(six["boundary"])),
            format_q(int(six["normalization"])),
            format_q(six_seam_m0),
        )
    )
    print("k=6 canonical cut table=" + ",".join(format_q(int(value)) for value in six["sharp"]))
    print(
        "k=6 sparse-Frobenius cut table="
        + ",".join(format_q(int(value)) for value in six["sparse_sharp"])
    )
    print("all gates use explicit check(); optimized mode is fully active")


if __name__ == "__main__":
    main()
