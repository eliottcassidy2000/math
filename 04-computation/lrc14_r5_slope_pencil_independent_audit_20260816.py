#!/usr/bin/env python3
"""Independent exact audit of the THM-2594 r=5 slope pencil.

This program does not import either the canonical THM-2594 constructor or the
later root-shear/Fourier helper.  It reconstructs the retained common-base
table

    T_ell(u,theta) = sum_q N(u,q,ell,theta)

from the defining circle inequalities by a separate interval-set algebra and
a separate transfer-profile engine.  It then performs direct finite-field
character sums in F_547, checks all thirteen affine frequency slopes, and
compares the Fourier basis with the anchored-rectangle (mixed-Haar) basis.

The final schema audit is deliberately negative in type: an abstract binary
checkerboard embeds into any chosen 2x2 rectangle of the retained table, but
the currently visible U_full endpoint records have no base/root/sheet/address
map into the THM-2594 ancestry atoms.  No physical current is constructed.
"""

from __future__ import annotations

import hashlib
from bisect import bisect_right
from fractions import Fraction
from math import gcd


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def lcm(a: int, b: int) -> int:
    return a // gcd(a, b) * b


W = (1, 14, 27, 40, 53, 66, 13, 13**3, 2 * 13**5)
GUARD = 0
OWNER = 6
TARGET_A = 7
TARGET_B = 8
PACKET = 13**2
DEPTH = 13**5
GRID_LCM = 1
for weight in W:
    GRID_LCM = lcm(GRID_LCM, weight)
GRID = 182 * GRID_LCM
require(GRID == 297836897838480, "production grid changed")
require(GRID % (182 * max(W)) == 0, "grid does not resolve the largest comb")

I5 = Fraction(48602521488933856, 337437093630814766589)
RHO_B = Fraction(35505957232, 16132831966251)
EXPECTED_DENOMINATOR = PACKET * DEPTH * DEPTH * GRID

PENCIL_SHA256 = "e92e3f1b072db16ada1daa28925803ebd9e11658deb3532680911ed637dee85d"
OWNER_MARGINAL_SHA256 = "31475920623317779b3c6de6334f258309256fb71734cc6ecd7ea6a6476b3e68"
ABSOLUTE_LINE_SHA256 = "808df0ceb8616773cff1b5c12de2333d1495c07d119c57b9b49e68aa9bb2627f"
COMMON_BASE_TABLE_SHA256 = "18463c7af393c6090c7419c76d48b7ed40f9ac8b17d8c1be247062493e334ce8"
RETAINED_TABLE_SHA256 = "1ba72585187c014e482f89959ff4f19e2185e6e4fcd982d497923a304e0d37d8"
RECTANGLE_RAW_SHA256 = "801b278fbe7f638aa621df602672da45953da7d5b961d389797530b5517ca5f4"
RECTANGLE_FOURIER_SHA256 = "5f284f6693c1f3441f51d9a86fddf6bd5d5fcaa1f76b0a9bd8add437b5e3752f"


Interval = tuple[int, int]
Profile = tuple[list[int], list[int]]


def normalize(intervals: list[Interval]) -> list[Interval]:
    """Sort and merge an ordinary (non-circular) interval union."""

    answer: list[list[int]] = []
    for left, right in sorted(intervals):
        require(0 <= left <= right <= GRID, "interval outside the circle chart")
        if left == right:
            continue
        if answer and left <= answer[-1][1]:
            if right > answer[-1][1]:
                answer[-1][1] = right
        else:
            answer.append([left, right])
    return [(left, right) for left, right in answer]


def intersect(first: list[Interval], second: list[Interval]) -> list[Interval]:
    answer: list[Interval] = []
    i = j = 0
    while i < len(first) and j < len(second):
        left = max(first[i][0], second[j][0])
        right = min(first[i][1], second[j][1])
        if left < right:
            answer.append((left, right))
        if first[i][1] <= second[j][1]:
            i += 1
        else:
            j += 1
    return answer


def difference(first: list[Interval], removed: list[Interval]) -> list[Interval]:
    """Return first minus removed by an independently written two-pointer sweep."""

    answer: list[Interval] = []
    j = 0
    for left, right in first:
        while j < len(removed) and removed[j][1] <= left:
            j += 1
        cursor = left
        k = j
        while k < len(removed) and removed[k][0] < right:
            cut_left, cut_right = removed[k]
            if cursor < cut_left:
                answer.append((cursor, min(cut_left, right)))
            if cut_right > cursor:
                cursor = cut_right
            if cursor >= right:
                break
            k += 1
        if cursor < right:
            answer.append((cursor, right))
    return answer


def circular_comb(weight: int, period: int, low: int, high: int) -> list[Interval]:
    """Windows ((low+period*n)/(period*w),(high+period*n)/(period*w)).

    Endpoints are represented on the common integer grid.  Open/closed choices
    occur only on null boundaries and do not affect any exact integral below.
    """

    require(0 < high - low < period, "comb window must be proper")
    unit = GRID // (period * weight)
    step = period * unit
    width = (high - low) * unit
    base = (low % period) * unit
    intervals: list[Interval] = []
    for turn in range(weight):
        left = base + turn * step
        require(left < GRID, "comb representative escaped one turn")
        right = left + width
        if right <= GRID:
            intervals.append((left, right))
        else:
            intervals.append((left, GRID))
            intervals.append((0, right - GRID))
    return normalize(intervals)


NARROW_COMBS = [circular_comb(weight, 182, -13, 13) for weight in W]
WIDE_GUARD = circular_comb(W[GUARD], 91, -13, 13)


def pattern(required: tuple[int, ...], excluded: tuple[int, ...]) -> list[Interval]:
    allowed = [(0, GRID)]
    for index in required:
        allowed = intersect(allowed, NARROW_COMBS[index])
    allowed = difference(allowed, WIDE_GUARD)
    for index in excluded:
        allowed = difference(allowed, NARROW_COMBS[index])
    return allowed


# E_1, Q_{1,{b}}, and the source-safe-deleted terminal word T_b.
E_SET = pattern((OWNER,), (1, 2, 3, 4, 5, TARGET_A, TARGET_B))
QB_SET = pattern((TARGET_B,), (1, 2, 3, 4, 5, OWNER, TARGET_A))
TB_SET = pattern((TARGET_B,), (1, 2, 3, 4, 5, TARGET_A))


def measure(intervals: list[Interval]) -> int:
    return sum(right - left for left, right in intervals)


require(len(E_SET) == 57072, "independent E interval census changed")
require(Fraction(measure(E_SET), GRID) == Fraction(1882176, 28589561), "E measure changed")
require(len(QB_SET) == 131762, "independent Q_b interval census changed")
require(Fraction(measure(QB_SET), GRID) == Fraction(4839079319, 190921088358), "Q_b measure changed")
require(difference(TB_SET, NARROW_COMBS[OWNER]) == QB_SET, "T_b minus owner danger is not Q_b")


def add_event(events: dict[int, int], position: int, value: int) -> None:
    if value:
        events[position] = events.get(position, 0) + value


def cyclic_increment(events: dict[int, int], left: int, width: int, value: int) -> None:
    if width == 0 or value == 0:
        return
    right = left + width
    if right <= GRID:
        add_event(events, left, value)
        add_event(events, right, -value)
    else:
        add_event(events, left, value)
        add_event(events, GRID, -value)
        add_event(events, 0, value)
        add_event(events, right - GRID, -value)


def profile_mass(profile: Profile) -> int:
    starts, values = profile
    total = 0
    for index, left in enumerate(starts):
        right = starts[index + 1] if index + 1 < len(starts) else GRID
        total += values[index] * (right - left)
    return total


def pieces(profile: Profile) -> list[tuple[int, int, int]]:
    starts, values = profile
    return [
        (left, starts[index + 1] if index + 1 < len(starts) else GRID, values[index])
        for index, left in enumerate(starts)
        if values[index]
    ]


def transfer(weighted_pieces: list[tuple[int, int, int]], multiplier: int) -> Profile:
    """Numerator profile of sum_j h((y+j)/multiplier)."""

    baseline = 0
    events: dict[int, int] = {}
    input_mass = 0
    for left, right, value in weighted_pieces:
        if value == 0:
            continue
        require(0 <= left < right <= GRID, "bad transfer piece")
        input_mass += value * (right - left)
        wraps, remainder = divmod(multiplier * (right - left), GRID)
        baseline += wraps * value
        if remainder:
            cyclic_increment(events, (multiplier * left) % GRID, remainder, value)

    current = baseline + events.get(0, 0)
    starts = [0]
    values = [current]
    for position in sorted(position for position in events if 0 < position < GRID):
        current += events[position]
        if current != values[-1]:
            starts.append(position)
            values.append(current)
    result = (starts, values)
    require(min(values) >= 0, "negative transfer multiplicity")
    require(profile_mass(result) == multiplier * input_mass, "transfer mass failure")
    return result


def restrict_profile(profile: Profile, mask: list[Interval]) -> list[tuple[int, int, int]]:
    """Weighted nonzero pieces of profile times an indicator mask."""

    source = pieces(profile)
    answer: list[tuple[int, int, int]] = []
    i = j = 0
    while i < len(source) and j < len(mask):
        left = max(source[i][0], mask[j][0])
        right = min(source[i][1], mask[j][1])
        if left < right:
            value = source[i][2]
            if answer and answer[-1][1] == left and answer[-1][2] == value:
                answer[-1] = (answer[-1][0], right, value)
            else:
                answer.append((left, right, value))
        if source[i][1] <= mask[j][1]:
            i += 1
        else:
            j += 1
    return answer


def chart(profile: Profile, root: int) -> Profile:
    """Profile y -> h((y+root)/13) on the base chart."""

    starts, values = profile
    window = GRID // 13
    low = root * window
    high = low + window
    index = bisect_right(starts, low) - 1
    out_starts: list[int] = []
    out_values: list[int] = []
    while index < len(starts) and starts[index] < high:
        left = max(low, starts[index])
        right = min(high, starts[index + 1] if index + 1 < len(starts) else GRID)
        if left < right:
            mapped = 13 * left - root * GRID
            value = values[index]
            if not out_starts or out_values[-1] != value:
                out_starts.append(mapped)
                out_values.append(value)
        index += 1
    require(out_starts and out_starts[0] == 0, "chart does not cover the base origin")
    require(13 * high - root * GRID == GRID, "chart endpoint mismatch")
    return out_starts, out_values


def multiply(first: Profile, second: Profile) -> tuple[list[int], list[int], list[int], int]:
    """Pointwise product with a cumulative integral table."""

    s1, v1 = first
    s2, v2 = second
    i = j = 0
    starts = [0]
    values = [v1[0] * v2[0]]
    while True:
        next1 = s1[i + 1] if i + 1 < len(s1) else GRID
        next2 = s2[j + 1] if j + 1 < len(s2) else GRID
        boundary = min(next1, next2)
        if boundary >= GRID:
            break
        if boundary == next1:
            i += 1
        if boundary == next2:
            j += 1
        value = v1[i] * v2[j]
        if value != values[-1]:
            starts.append(boundary)
            values.append(value)
    cumulative = [0] * len(starts)
    for index in range(1, len(starts)):
        cumulative[index] = cumulative[index - 1] + values[index - 1] * (starts[index] - starts[index - 1])
    total = cumulative[-1] + values[-1] * (GRID - starts[-1])
    return starts, values, cumulative, total


def antiderivative(product: tuple[list[int], list[int], list[int], int], point: int) -> int:
    starts, values, cumulative, _ = product
    index = bisect_right(starts, point) - 1
    return cumulative[index] + values[index] * (point - starts[index])


def integrate(product: tuple[list[int], list[int], list[int], int], intervals: list[Interval]) -> int:
    return sum(antiderivative(product, right) - antiderivative(product, left) for left, right in intervals)


def prepare_mask(intervals: list[Interval]) -> tuple[list[int], list[int], list[int]]:
    starts: list[int] = []
    ends: list[int] = []
    mass_before: list[int] = []
    total = 0
    for left, right in intervals:
        starts.append(left)
        ends.append(right)
        mass_before.append(total)
        total += right - left
    return starts, ends, mass_before


def mask_antiderivative(mask: tuple[list[int], list[int], list[int]], point: int) -> int:
    starts, ends, mass_before = mask
    index = bisect_right(starts, point) - 1
    if index < 0:
        return 0
    return mass_before[index] + max(0, min(point, ends[index]) - starts[index])


def integrate_profile_over_mask(
    product: tuple[list[int], list[int], list[int], int],
    mask: tuple[list[int], list[int], list[int]],
) -> int:
    """Integrate by the short product profile, independently of mask intervals."""

    starts, values, _, _ = product
    total = 0
    for index, left in enumerate(starts):
        right = starts[index + 1] if index + 1 < len(starts) else GRID
        if values[index]:
            total += values[index] * (
                mask_antiderivative(mask, right) - mask_antiderivative(mask, left)
            )
    return total


def build_masks() -> list[list[list[Interval]]]:
    cells: list[list[Interval]] = []
    half = GRID // 14
    seventh = GRID // 7
    for ell in range(7):
        left = (ell * seventh - half) % GRID
        right = left + 2 * half
        cells.append(normalize([(left, min(right, GRID)), (0, max(0, right - GRID))]))

    twenty_eighth = GRID // 28
    deep: list[list[Interval]] = [[] for _ in range(13)]
    deep[0] = [(0, 13 * twenty_eighth)]
    deep[1] = [(twenty_eighth, 27 * twenty_eighth)]
    deep[2] = [(15 * twenty_eighth, GRID)]

    masks = [[[] for _ in range(13)] for _ in range(7)]
    for ell in range(7):
        word_danger = circular_comb(W[OWNER], 182, 26 * ell - 13, 26 * ell + 13)
        word_safe = difference(TB_SET, word_danger)
        cell_word = intersect(cells[ell], word_safe)
        for theta in range(3):
            masks[ell][theta] = intersect(cell_word, deep[theta])
    return masks


def common_base_table() -> tuple[list[list[list[list[int]]]], list[list[list[int]]], dict[str, object]]:
    """Build N(u,q,ell,theta), then retain q by two independent routes."""

    e_pieces = [(left, right, 1) for left, right in E_SET]
    packet_e = transfer(e_pieces, PACKET)
    f_pieces = restrict_profile(packet_e, QB_SET)
    require(Fraction(sum(value * (right - left) for left, right, value in f_pieces), PACKET * GRID) == RHO_B, "f integral changed")

    u_profile = transfer(f_pieces, DEPTH)
    v_profile = transfer(e_pieces, DEPTH)
    source_sum = transfer(pieces(v_profile), 13)
    masks = build_masks()
    prepared_masks = [[prepare_mask(masks[ell][theta]) for theta in range(13)] for ell in range(7)]

    owner_charts = [chart(u_profile, owner_root) for owner_root in range(13)]
    source_charts = [chart(v_profile, source_root) for source_root in range(13)]
    joint = [[[[0] * 13 for _ in range(7)] for _ in range(13)] for _ in range(13)]
    service_numerator = 0
    for owner_root in range(13):
        for source_root in range(13):
            product = multiply(owner_charts[owner_root], source_charts[source_root])
            service_numerator += product[3]
            for ell in range(7):
                for theta in range(3):
                    joint[owner_root][source_root][ell][theta] = integrate_profile_over_mask(
                        product, prepared_masks[ell][theta]
                    )

    table = [
        [
            [sum(joint[owner_root][source_root][ell][theta] for source_root in range(13)) for theta in range(13)]
            for owner_root in range(13)
        ]
        for ell in range(7)
    ]

    # A second route folds the source charts before multiplication.  It must
    # reproduce every retained entry, not merely the aggregate response.
    folded_table = [[[0] * 13 for _ in range(13)] for _ in range(7)]
    for owner_root in range(13):
        product = multiply(owner_charts[owner_root], source_sum)
        for ell in range(7):
            for theta in range(3):
                folded_table[ell][owner_root][theta] = integrate_profile_over_mask(
                    product, prepared_masks[ell][theta]
                )
    require(folded_table == table, "pairwise and source-folded retained tables differ")

    require(Fraction(service_numerator, EXPECTED_DENOMINATOR) == 169 * I5, "service anchor changed")
    require(all(table[0][u][t] == 0 for u in range(13) for t in range(13)), "ell=0 row should vanish")
    return joint, table, {
        "E_intervals": len(E_SET),
        "QB_intervals": len(QB_SET),
        "TB_intervals": len(TB_SET),
        "packet_profile_pieces": len(packet_e[0]),
        "U_profile_pieces": len(u_profile[0]),
        "V_profile_pieces": len(v_profile[0]),
        "source_sum_pieces": len(source_sum[0]),
        "service_numerator": service_numerator,
        "joint_entries": 13 * 13 * 7 * 13,
        "pairwise_vs_folded_retained": "PASS",
    }


PRIME = 547
ROOT_91 = 64
ROOT_7 = pow(ROOT_91, 13, PRIME)
ROOT_13 = pow(ROOT_91, 7, PRIME)
require(ROOT_7 == 81 and ROOT_13 == 475, "compatible split roots changed")
require(pow(ROOT_91, 91, PRIME) == 1 and pow(ROOT_91, 7, PRIME) != 1 and pow(ROOT_91, 13, PRIME) != 1, "64 is not primitive of order 91")
require(pow(ROOT_7, 7, PRIME) == 1 and ROOT_7 != 1, "bad order-seven root")
require(pow(ROOT_13, 13, PRIME) == 1 and ROOT_13 != 1, "bad order-thirteen root")

POW7 = [[pow(ROOT_7, (-mode * value) % 7, PRIME) for value in range(7)] for mode in range(7)]
POW13 = [[pow(ROOT_13, (-mode * value) % 13, PRIME) for value in range(13)] for mode in range(13)]


def transform(table: list[list[list[int]]], word_mode: int, owner_mode: int, root_mode: int) -> int:
    total = 0
    for ell in range(7):
        word_weight = POW7[word_mode][ell]
        for owner in range(13):
            owner_weight = POW13[owner_mode][owner]
            for root in range(13):
                total += table[ell][owner][root] * word_weight * owner_weight * POW13[root_mode][root]
    return total % PRIME


def slope_audit(table: list[list[list[int]]]) -> dict[str, object]:
    full = [
        (beta, owner_mode, root_mode, transform(table, beta, owner_mode, root_mode))
        for beta in range(1, 7)
        for owner_mode in range(13)
        for root_mode in range(1, 13)
    ]
    require(all(value for _, _, _, value in full), "primitive slope-pencil zero")
    full_digest = hashlib.sha256(repr(full).encode("ascii")).hexdigest()
    require(full_digest == PENCIL_SHA256, "candidate slope-pencil ledger mismatch")

    by_slope = tuple(
        sum(owner_mode == slope * root_mode % 13 and value != 0 for _, owner_mode, root_mode, value in full)
        for slope in range(13)
    )
    require(by_slope == (72,) * 13, "affine slope census changed")

    owner_marginal = [
        (beta, root_mode, transform(table, beta, 0, root_mode))
        for beta in range(1, 7)
        for root_mode in range(1, 13)
    ]
    absolute_line = [
        (beta, root_mode, transform(table, beta, 2 * root_mode % 13, root_mode))
        for beta in range(1, 7)
        for root_mode in range(1, 13)
    ]
    owner_digest = hashlib.sha256(repr(owner_marginal).encode("ascii")).hexdigest()
    absolute_digest = hashlib.sha256(repr(absolute_line).encode("ascii")).hexdigest()
    require(owner_digest == OWNER_MARGINAL_SHA256, "owner-marginal ledger mismatch")
    require(absolute_digest == ABSOLUTE_LINE_SHA256, "absolute-root ledger mismatch")

    # Directly verify the shear A_hat(r,s)=S_hat(r+2s,s), without constructing
    # the candidate helper's absolute table or calling any helper transform.
    for beta in range(1, 7):
        for owner_mode in range(13):
            for root_mode in range(1, 13):
                absolute_direct = 0
                for ell in range(7):
                    for owner in range(13):
                        for theta in range(13):
                            absolute_root = (theta + 2 * owner) % 13
                            absolute_direct += (
                                table[ell][owner][theta]
                                * POW7[beta][ell]
                                * POW13[owner_mode][owner]
                                * POW13[root_mode][absolute_root]
                            )
                require(absolute_direct % PRIME == transform(table, beta, (owner_mode + 2 * root_mode) % 13, root_mode), "direct absolute/slaved shear failed")

    return {
        "full_digest": full_digest,
        "owner_digest": owner_digest,
        "absolute_digest": absolute_digest,
        "by_slope": by_slope,
        "nonaxial_fourier": sum(
            transform(table, beta, owner_mode, root_mode) != 0
            for beta in range(1, 7)
            for owner_mode in range(1, 13)
            for root_mode in range(1, 13)
        ),
    }


def rectangle_audit(table: list[list[list[int]]]) -> dict[str, object]:
    rectangles = [
        [
            [
                table[ell][owner][root]
                - table[ell][owner][0]
                - table[ell][0][root]
                + table[ell][0][0]
                for root in range(1, 13)
            ]
            for owner in range(1, 13)
        ]
        for ell in range(7)
    ]
    raw_counts = tuple(
        sum(rectangles[ell][owner - 1][root - 1] != 0 for owner in range(1, 13) for root in range(1, 13))
        for ell in range(1, 7)
    )
    require(raw_counts == (132, 132, 133, 4, 4, 4), "raw anchored-rectangle census changed")
    raw_records = [
        (ell, owner, root, rectangles[ell][owner - 1][root - 1])
        for ell in range(1, 7)
        for owner in range(1, 13)
        for root in range(1, 13)
    ]
    raw_digest = hashlib.sha256(repr(raw_records).encode("ascii")).hexdigest()
    require(raw_digest == RECTANGLE_RAW_SHA256, "raw anchored-rectangle ledger mismatch")
    first_raw = next(record for record in raw_records if record[3] != 0)
    require(first_raw == (1, 1, 2, 258805184656089173356800), "first raw rectangle changed")

    fourier_records = []
    for beta in range(1, 7):
        for owner in range(1, 13):
            for root in range(1, 13):
                value = sum(rectangles[ell][owner - 1][root - 1] * POW7[beta][ell] for ell in range(7)) % PRIME
                fourier_records.append((beta, owner, root, value))
    fourier_digest = hashlib.sha256(repr(fourier_records).encode("ascii")).hexdigest()
    require(fourier_digest == RECTANGLE_FOURIER_SHA256, "anchored-rectangle Fourier ledger mismatch")
    counts = tuple(sum(beta0 == beta and value != 0 for beta0, _, _, value in fourier_records) for beta in range(1, 7))
    require(counts == (134,) * 6, "anchored-rectangle Fourier census changed")
    first_fourier = next(record for record in fourier_records if record[3] != 0)
    require(first_fourier == (1, 1, 2, 180), "first Fourier rectangle changed")
    expected_zeros = ((1, 1), (2, 1), (3, 1), (4, 1), (5, 1), (7, 1), (8, 1), (9, 1), (10, 1), (11, 1))
    zeros_by_beta = tuple(
        tuple((owner, root) for beta0, owner, root, value in fourier_records if beta0 == beta and value == 0)
        for beta in range(1, 7)
    )
    require(zeros_by_beta == (expected_zeros,) * 6, "anchored-rectangle zero locus changed")

    # Anchored rectangles and nontrivial owner/root characters are two bases
    # of the same 12x12 mixed quotient.  Verify the change of basis directly:
    # for r,s != 0, T_hat(r,s) is the DFT of Delta(u,t), u,t=1..12.
    for beta in range(1, 7):
        rectangle_beta = [
            [
                sum(rectangles[ell][owner - 1][root - 1] * POW7[beta][ell] for ell in range(7)) % PRIME
                for root in range(1, 13)
            ]
            for owner in range(1, 13)
        ]
        for owner_mode in range(1, 13):
            for root_mode in range(1, 13):
                changed = sum(
                    rectangle_beta[owner - 1][root - 1]
                    * POW13[owner_mode][owner]
                    * POW13[root_mode][root]
                    for owner in range(1, 13)
                    for root in range(1, 13)
                ) % PRIME
                require(changed == transform(table, beta, owner_mode, root_mode), "anchored/Fourier basis change failed")

    def rank_mod(matrix: list[list[int]]) -> int:
        work = [[value % PRIME for value in row] for row in matrix]
        rank = 0
        for column in range(len(work[0])):
            pivot = next((row for row in range(rank, len(work)) if work[row][column]), None)
            if pivot is None:
                continue
            work[rank], work[pivot] = work[pivot], work[rank]
            inverse = pow(work[rank][column], PRIME - 2, PRIME)
            work[rank] = [value * inverse % PRIME for value in work[rank]]
            for row in range(len(work)):
                if row != rank and work[row][column]:
                    factor = work[row][column]
                    work[row] = [(a - factor * b) % PRIME for a, b in zip(work[row], work[rank])]
            rank += 1
        return rank

    change_matrix = [[POW13[mode][coordinate] for coordinate in range(1, 13)] for mode in range(1, 13)]
    change_rank = rank_mod(change_matrix)
    require(change_rank == 12, "owner/root Fourier change matrix is singular")

    return {
        "raw_counts": raw_counts,
        "raw_digest": raw_digest,
        "first_raw": first_raw,
        "fourier_counts": counts,
        "fourier_digest": fourier_digest,
        "first_fourier": first_fourier,
        "zeros": expected_zeros,
        "one_axis_change_rank": change_rank,
        "tensor_change_rank": change_rank * change_rank,
    }


def schema_audit() -> dict[str, object]:
    """Separate formal checkerboard capacity from a lawful U_full atom map."""

    ufull_visible = (
        "circle_interval",
        "E_factor_lineage",
        "Q_factor_lineage",
        "E_component_index",
        "Q_component_index",
        "periodic_wrap_turn",
        "boundary_provenance",
    )
    ancestry_required = (
        "base",
        "root",
        "owner_sheet",
        "word_sheet",
        "source_sheet",
        "left_horizon",
        "right_horizon",
        "address",
    )
    missing = tuple(field for field in ancestry_required if field not in ufull_visible)
    require(missing == ancestry_required, "a typed ancestry field unexpectedly appeared")

    # Every selected target rectangle admits this formal, ell-independent
    # label embedding.  It preserves the checkerboard functional exactly.
    sample = ((11, 17), (23, 31))
    pushed = {(0, 0): sample[0][0], (0, 1): sample[0][1], (1, 0): sample[1][0], (1, 1): sample[1][1]}
    source_haar = sample[0][0] - sample[0][1] - sample[1][0] + sample[1][1]
    target_rectangle = pushed[(1, 1)] - pushed[(1, 0)] - pushed[(0, 1)] + pushed[(0, 0)]
    require(source_haar == target_rectangle, "formal checkerboard embedding changed sign")

    checker = ((1, -1), (-1, 1))
    row_margins = tuple(sum(row) for row in checker)
    column_margins = tuple(sum(checker[i][j] for i in range(2)) for j in range(2))
    require(row_margins == (0, 0) and column_margins == (0, 0), "checkerboard hostile gained a marginal")
    return {
        "source": "actual pre-merge U_full factor-labelled endpoint cells",
        "target": "retained THM-2594 tables T_ell(u,theta)",
        "formal_map": "(i,j)->(u_i,theta_j), with u_0=theta_0=0, fixed before ell",
        "preserved": "abstract 2x2 row/column margins and checkerboard pairing",
        "lost": "circle-cell measure, E/Q lineage, q-root, common base, sheets, horizons, F_13^3 address",
        "missing": missing,
        "hostile": checker,
        "lawful_map": "NOT_CONSTRUCTED_BY_CURRENT_API",
    }


def main() -> None:
    joint, table, construction = common_base_table()
    joint_digest = hashlib.sha256(repr(joint).encode("ascii")).hexdigest()
    table_digest = hashlib.sha256(repr(table).encode("ascii")).hexdigest()
    require(joint_digest == COMMON_BASE_TABLE_SHA256, "full common-base table differs from the canonical constructor")
    require(table_digest == RETAINED_TABLE_SHA256, "retained table differs from the canonical constructor")
    slope = slope_audit(table)
    rectangles = rectangle_audit(table)
    schema = schema_audit()

    print("LRC14 R=5 PROJECTIVE SLOPE-PENCIL INDEPENDENT AUDIT")
    print("status=ACCEPT_FINITE_EXACT_RETAINED_TABLE_SIDECAR; LRC(14) OPEN")
    print(f"constructor=standalone interval Boolean algebra + transfer profiles; imports_primary_or_helper=False")
    print(f"grid={GRID}; denominator={EXPECTED_DENOMINATOR}")
    print(f"construction={construction}")
    print(f"common_base_table_sha256={joint_digest}")
    print(f"retained_table_sha256={table_digest}")
    print(f"split_embedding=(prime={PRIME},root91={ROOT_91},root7={ROOT_7},root13={ROOT_13})")
    print(f"slope_census={slope['by_slope']}; total=936/936")
    print(f"full_slope_pencil_sha256={slope['full_digest']}")
    print(f"slope0_owner_marginal=72/72; sha256={slope['owner_digest']}; meaning=THM-2594/THM-2512 mixed block")
    print(f"slope2_fixed_absolute_root=72/72; sha256={slope['absolute_digest']}; meaning=owner-marginal after t=theta+2u shear")
    print(f"nonaxial_owner_root_Fourier={slope['nonaxial_fourier']}/864")
    print(f"anchored_rectangles_raw_by_ell_1_to_6={rectangles['raw_counts']}; first={rectangles['first_raw']}")
    print(f"anchored_rectangles_raw_sha256={rectangles['raw_digest']}")
    print(f"anchored_rectangles_C7_by_beta={rectangles['fourier_counts']}; total=804/864")
    print(f"anchored_rectangles_C7_first={rectangles['first_fourier']}; sha256={rectangles['fourier_digest']}")
    print(f"anchored_rectangles_common_zero_locus={rectangles['zeros']}")
    print(f"anchored_to_Fourier_change_of_basis=(one_axis_rank={rectangles['one_axis_change_rank']},tensor_rank={rectangles['tensor_change_rank']})")
    print("basis_verdict=all 864/864 nonaxial Fourier characters fire, while 804/864 anchored mixed-Haar coordinates fire; support is basis-dependent")
    print(f"connection_contract={schema}")
    print("typed_verdict=the table has formal mixed-Haar capacity, but no lawful ell-independent U_full cell-to-ancestry/address map is currently supplied")
    print("nonconsequence=no physical current, K4 bridge, row exclusion, or LRC(14) consequence")
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
