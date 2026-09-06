#!/usr/bin/env python3
"""Exact scratch audit for the dyadic five-ray continuation.

This file is deliberately self-contained and uses Fraction arithmetic only.
It checks:

* the equal- and mixed-radius danger-intersection formulae and their sharp
  ten-label order-statistic consequences;
* the improved q=2 one-even Haar-entry constants and a hybrid q=4 floor;
* the complete open-component atlas for the five odd-3-unit cross-combs;
* exact q=4/q=2 structured safe-set reconstruction;
* finite scale banks, endpoint coverage, component addresses, and direct
  physical witnesses for residual-band hostile bodies.

No assertion statement is used, so every gate remains active under python -O.
"""

from fractions import Fraction as Q
from math import gcd


DELTA = Q(1, 14)
FIVE = ((1, 11), (1, 23), (5, 11), (1, 37), (1, 25))


def need(condition: bool, label: object) -> None:
    if not condition:
        raise RuntimeError(label)


def word(value: Q) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def interval_word(interval: tuple[Q, Q]) -> str:
    return f"({word(interval[0])},{word(interval[1])})"


def circle_distance(value: Q) -> Q:
    value %= 1
    return min(value, 1 - value)


def bad(speed: int, phase: Q, radius: Q = DELTA) -> bool:
    return circle_distance(speed * phase) < radius


def safe(speeds: tuple[int, ...], phase: Q) -> bool:
    return all(not bad(speed, phase) for speed in speeds)


def b2(value: Q) -> Q:
    value %= 1
    return value * value - value + Q(1, 6)


def danger_overlap_formula(a: int, b: int) -> Q:
    """Haar measure of D_a intersect D_b for radius 1/14."""
    g = gcd(a, b)
    p, q = a // g, b // g
    return Q(1, 49) + (
        b2(Q(q - p, 14)) - b2(Q(q + p, 14))
    ) / (p * q)


def danger_overlap_literal(a: int, b: int) -> Q:
    walls = {Q(0), Q(1)}
    for speed in (a, b):
        for k in range(speed):
            walls.add((Q(k) + DELTA) / speed)
            walls.add((Q(k + 1) - DELTA) / speed)
    ordered = sorted(walls)
    return sum(
        right - left
        for left, right in zip(ordered, ordered[1:])
        if bad(a, (left + right) / 2) and bad(b, (left + right) / 2)
    )


def mixed_overlap_formula(c: int, r: int) -> Q:
    """Haar measure of {||cy||<1/14} intersect {||ry||<1/7}."""
    g = gcd(c, r)
    p, q = c // g, r // g
    return Q(2, 49) + (
        b2(Q(q - 2 * p, 14)) - b2(Q(q + 2 * p, 14))
    ) / (p * q)


def mixed_overlap_literal(c: int, r: int) -> Q:
    walls = {Q(0), Q(1)}
    for speed, radius in ((c, DELTA), (r, 2 * DELTA)):
        for k in range(speed):
            walls.add((Q(k) + radius) / speed)
            walls.add((Q(k + 1) - radius) / speed)
    ordered = sorted(walls)
    return sum(
        right - left
        for left, right in zip(ordered, ordered[1:])
        if bad(c, (left + right) / 2)
        and bad(r, (left + right) / 2, 2 * DELTA)
    )


def safe_components(speeds: tuple[int, ...]) -> tuple[tuple[Q, Q], ...]:
    """All closed components, including isolated equality-safe points."""
    walls = {Q(0), Q(1)}
    for speed in speeds:
        for k in range(speed):
            walls.add((Q(k) + DELTA) / speed)
            walls.add((Q(k + 1) - DELTA) / speed)
    ordered = sorted(walls)
    safe_cells = [
        safe(speeds, (left + right) / 2)
        for left, right in zip(ordered, ordered[1:])
    ]
    intervals: list[tuple[Q, Q]] = []
    for index, yes in enumerate(safe_cells):
        if not yes:
            continue
        left, right = ordered[index], ordered[index + 1]
        need(safe(speeds, left) and safe(speeds, right), ("closed safe cell", speeds, left, right))
        if intervals and intervals[-1][1] == left:
            intervals[-1] = (intervals[-1][0], right)
        else:
            intervals.append((left, right))
    for index, point in enumerate(ordered):
        if not safe(speeds, point):
            continue
        left_safe = index > 0 and safe_cells[index - 1]
        right_safe = index < len(safe_cells) and safe_cells[index]
        if not left_safe and not right_safe:
            intervals.append((point, point))
    return tuple(sorted(intervals))


def normalize_closed(intervals: list[tuple[Q, Q]]) -> tuple[tuple[Q, Q], ...]:
    result: list[tuple[Q, Q]] = []
    for left, right in sorted(intervals):
        need(left <= right, ("closed interval orientation", left, right))
        if result and left <= result[-1][1]:
            result[-1] = (result[-1][0], max(result[-1][1], right))
        else:
            result.append((left, right))
    return tuple(result)


def intersect_closed(
    first: tuple[tuple[Q, Q], ...],
    second: tuple[tuple[Q, Q], ...],
) -> tuple[tuple[Q, Q], ...]:
    pieces = []
    i = j = 0
    while i < len(first) and j < len(second):
        left = max(first[i][0], second[j][0])
        right = min(first[i][1], second[j][1])
        if left <= right:
            pieces.append((left, right))
        if first[i][1] < second[j][1]:
            i += 1
        else:
            j += 1
    return normalize_closed(pieces)


def q2_reconstruct(C: tuple[int, ...], r: int) -> tuple[tuple[Q, Q], ...]:
    return intersect_closed(safe_components(C), safe_components((r,)))


def q4_reconstruct(C: tuple[int, ...], r: int) -> tuple[tuple[Q, Q], ...]:
    lifted = []
    for left, right in safe_components(C):
        lifted.append((left / 2, right / 2))
        lifted.append(((left + 1) / 2, (right + 1) / 2))
    return intersect_closed(normalize_closed(lifted), safe_components((r,)))


def both_sheets_bad(p: int, q: int, y: Q) -> bool:
    return all(
        bad(p, (y + j) / 2) or bad(q, (y + j) / 2)
        for j in (0, 1)
    )


def cross_components(p: int, q: int) -> tuple[tuple[Q, Q], ...]:
    """Maximal open components of the primitive quotient cross-comb."""
    walls = {Q(0), Q(1)}
    for shift in (Q(0), Q(1, 2)):
        for speed in (p, q):
            for integer in range(speed):
                for sign in (-1, 1):
                    wall = (2 * ((Q(integer) + sign * DELTA) / speed - shift)) % 1
                    walls.add(wall)
    ordered = sorted(walls)
    cells = [
        (left, right)
        for left, right in zip(ordered, ordered[1:])
        if both_sheets_bad(p, q, (left + right) / 2)
    ]
    result: list[tuple[Q, Q]] = []
    for left, right in cells:
        if result and result[-1][1] == left and both_sheets_bad(p, q, left):
            result[-1] = (result[-1][0], right)
        else:
            result.append((left, right))
    need(all(not both_sheets_bad(p, q, endpoint) for interval in result for endpoint in interval),
         ("cross components are open", p, q))
    return tuple(result)


def cross_measure(p: int, q: int) -> Q:
    return sum((right - left for left, right in cross_components(p, q)), Q())


def component_address(
    component: tuple[Q, Q],
    primitive_cells: tuple[tuple[Q, Q], ...],
    scale: int,
) -> tuple[int, int] | None:
    left, right = component
    for sheet in range(scale):
        for cell_index, (a, b) in enumerate(primitive_cells):
            if (sheet + a) / scale < left and right < (sheet + b) / scale:
                return sheet, cell_index
    return None


def scale_bank(beta: Q, longest: Q) -> tuple[int, ...]:
    need(longest > 0, "positive component supplied by LRC up to eleven")
    cap = beta / longest
    return tuple(
        t for t in range(1, cap.numerator // cap.denominator + 1, 2)
        if t % 3 and Q(t) < cap
    )


def pulled_boundaries(cells: tuple[tuple[Q, Q], ...], scale: int) -> tuple[Q, ...]:
    return tuple(sorted({
        (sheet + endpoint) / scale
        for sheet in range(scale)
        for interval in cells
        for endpoint in interval
    }))


def row_witness(
    H: tuple[int, ...], p: int, q: int, scale: int,
    components: tuple[tuple[Q, Q], ...],
    cells: tuple[tuple[Q, Q], ...],
) -> tuple[Q, int, Q, Q]:
    candidates = {endpoint for interval in components for endpoint in interval}
    candidates.update((left + right) / 2 for left, right in components)
    for boundary in pulled_boundaries(cells, scale):
        if any(left <= boundary <= right for left, right in components):
            candidates.add(boundary)
    row = tuple(2 * h for h in H) + (scale * p, scale * q)
    for y in sorted(candidates):
        if not safe(H, y):
            continue
        for j in (0, 1):
            x = (y + j) / 2
            if not bad(scale * p, x) and not bad(scale * q, x):
                clearance = min(circle_distance(v * x) for v in row)
                need(clearance >= DELTA, ("direct row witness", H, p, q, scale, x, clearance))
                return y, j, x, clearance
    raise RuntimeError(("missing witness despite failed containment", H, p, q, scale))


def divisor_complete(H: tuple[int, ...]) -> bool:
    return all(any(h % modulus == 0 for h in H) for modulus in range(2, 15))


def audit_overlap_order_statistic() -> int:
    literal = 0
    for p in range(1, 61):
        for q in range(p + 1, 81):
            if gcd(p, q) != 1:
                continue
            need(danger_overlap_formula(p, q) == danger_overlap_literal(p, q),
                 ("overlap formula", p, q))
            literal += 1

    # Osc(B_2)=1/4.  The exact tenth order statistic is controlled by which
    # of 2 and 3 may divide the reduced denominator q | r.  Each finite list
    # below is global: beyond `cutoff`, the oscillation bound is already
    # strictly larger than lambda.
    tier_specs = (
        (
            "gcd1", 1, lambda q: gcd(q, 6) == 1, Q(1, 63), Q(8, 63), 55,
            {Q(1, 13), Q(1, 11), Q(2, 11), Q(3, 11), Q(10), Q(11), Q(12), Q(13)},
            {Q(9, 5), Q(9), Q(27)},
            715,
            (Q(1, 13), Q(13), Q(12), Q(1, 11), Q(2, 11),
             Q(3, 11), Q(11), Q(10), Q(9, 5), Q(9)),
        ),
        (
            "gcd2", 2, lambda q: q % 3 != 0, Q(1, 70), Q(9, 70), 40,
            {Q(1, 13), Q(1, 11), Q(2, 11), Q(3, 11), Q(11, 2), Q(11), Q(12), Q(13)},
            {Q(1, 10), Q(3, 10), Q(10)},
            1430,
            (Q(1, 13), Q(1, 11), Q(2, 11), Q(3, 11), Q(11, 2),
             Q(11), Q(12), Q(13), Q(1, 10), Q(3, 10)),
        ),
        (
            "gcd3", 3, lambda q: q % 2 != 0, Q(1, 70), Q(9, 70), 40,
            {Q(1, 13), Q(1, 11), Q(2, 11), Q(3, 11), Q(11, 3), Q(11), Q(12), Q(13)},
            {Q(10, 3), Q(10)},
            429,
            (Q(1, 13), Q(1, 11), Q(2, 11), Q(3, 11), Q(11, 3),
             Q(11), Q(12), Q(13), Q(10, 3), Q(10)),
        ),
        (
            "gcd6", 6, lambda q: True, Q(1, 77), Q(10, 77), 33,
            {Q(1, 13), Q(1, 12), Q(12), Q(13)},
            {Q(1, 11), Q(2, 11), Q(3, 11), Q(11, 3), Q(11, 2), Q(11)},
            1716,
            (Q(1, 13), Q(1, 12), Q(12), Q(13), Q(1, 11),
             Q(2, 11), Q(3, 11), Q(11, 3), Q(11, 2), Q(11)),
        ),
    )
    tier_controls: dict[str, tuple[int, tuple[int, ...]]] = {}
    for name, gcd_class, denominator_ok, lam, loss, cutoff, expected_below, expected_equal, r, ratios in tier_specs:
        need(Q(1, 49) - Q(1, 4 * (cutoff + 1)) > lam,
             ("global order-stat cutoff", name))
        below: set[Q] = set()
        equal: set[Q] = set()
        for p in range(1, cutoff + 1):
            for q in range(1, cutoff + 1):
                if p == q or p * q > cutoff or gcd(p, q) != 1 or not denominator_ok(q):
                    continue
                value = danger_overlap_formula(p, q)
                if value < lam:
                    below.add(Q(p, q))
                elif value == lam:
                    equal.add(Q(p, q))
        need(below == expected_below, ("subcritical orientations", name, below))
        need(equal == expected_equal, ("boundary orientations", name, equal))
        C = tuple(int(r * ratio) for ratio in ratios)
        need(gcd(r, 6) == gcd_class and gcd(*C) == 1 and r not in C,
             ("sharp tier arithmetic", name, r, C))
        need(len(set(C)) == 10
             and all(Q(c, r) == ratio for c, ratio in zip(C, ratios)),
             ("sharp ten-label realization", name, C))
        overlaps = tuple(danger_overlap_formula(c, r) for c in C)
        need(max(overlaps) == lam, ("sharp tenth-order overlap", name, overlaps))
        need(Q(1, 7) - lam == loss, ("q2 loss", name))
        tier_controls[name] = (r, C)
        print(
            f"overlap_tier={name};lambda={word(lam)};loss={word(loss)};"
            "below=" + ",".join(word(x) for x in sorted(below)) +
            ";equal=" + ",".join(word(x) for x in sorted(equal)) +
            f";sharp_r={r};sharp_C=" + ",".join(map(str, C))
        )

    # New q=2 one-even constants.
    pair_mass = {
        (1, 11): Q(4, 77),
        (1, 23): Q(8, 161),
        (5, 11): Q(18, 385),
        (1, 37): Q(12, 259),
        (1, 25): Q(8, 175),
    }
    q2_gates = {ratio: value + Q(8, 63) for ratio, value in pair_mass.items()}
    expected_gates = {
        (1, 11): Q(124, 693),
        (1, 23): Q(256, 1449),
        (5, 11): Q(86, 495),
        (1, 37): Q(404, 2331),
        (1, 25): Q(272, 1575),
    }
    need(q2_gates == expected_gates, ("q2 ratio gates", q2_gates))
    need(Q(4, 91) + Q(8, 63) == Q(20, 117), "q2 five-ray entry")

    expected_entries = {
        "gcd1": (Q(4, 21), Q(124, 693), Q(20, 117)),
        "gcd2": (Q(121, 630), Q(139, 770), Q(157, 910)),
        "gcd3": (Q(121, 630), Q(139, 770), Q(157, 910)),
        "gcd6": (Q(134, 693), Q(2, 11), Q(174, 1001)),
    }
    for name, _, _, _, loss, _, _, _, _, _ in tier_specs:
        actual = (Q(4, 63) + loss, Q(4, 77) + loss, Q(4, 91) + loss)
        need(actual == expected_entries[name], ("q2 entry hierarchy", name, actual))
        print(
            f"q2_entry={name};general_pair={word(actual[0])};"
            f"three_unit_pair={word(actual[1])};five_ray={word(actual[2])}"
        )

    print("overlap_formula_literal_pairs=" + str(literal))
    print("q2_floor_by_gcd=1->mu(G_C)-8/63;2,3->mu(G_C)-9/70;6->mu(G_C)-10/77")
    print("q2_ratio_gates=" + ";".join(f"{p}:{q}->{word(q2_gates[(p,q)])}" for p, q in FIVE))

    # q=4 has a mixed-radius fibre loss.  Again the Bernoulli oscillation is
    # 1/4; now pq>=49 forces overlap at least 1/28.  The oriented denominator
    # is a divisor of the odd 3-unit r and hence is coprime to 6.
    mixed_literal = 0
    for p in range(1, 41):
        for q in range(1, 61):
            if p == q or gcd(p, q) != 1:
                continue
            need(mixed_overlap_formula(p, q) == mixed_overlap_literal(p, q),
                 ("mixed overlap formula", p, q))
            mixed_literal += 1
    need(Q(2, 49) - Q(1, 4 * 49) == Q(1, 28), "mixed product-49 cutoff")
    mixed_below: set[Q] = set()
    mixed_equal: set[Q] = set()
    for p in range(1, 50):
        for q in range(1, 50):
            if p == q or p * q > 49 or gcd(p, q) != 1 or gcd(q, 6) != 1:
                continue
            value = mixed_overlap_formula(p, q)
            if value < Q(1, 28):
                mixed_below.add(Q(p, q))
            elif value == Q(1, 28):
                mixed_equal.add(Q(p, q))
    expected_mixed_below = {
        Q(6), Q(1, 11), Q(5), Q(2, 11), Q(1, 13), Q(13), Q(1, 25)
    }
    expected_mixed_equal = {Q(4, 5), Q(4), Q(12), Q(20)}
    need(mixed_below == expected_mixed_below,
         ("seven mixed subcritical orientations", mixed_below))
    need(mixed_equal == expected_mixed_equal,
         ("four mixed boundary orientations", mixed_equal))

    mixed_r = 3575
    mixed_ratios = (
        Q(6), Q(1, 11), Q(5), Q(2, 11), Q(1, 13), Q(13), Q(1, 25),
        Q(4, 5), Q(4), Q(12),
    )
    mixed_C = tuple(int(mixed_r * ratio) for ratio in mixed_ratios)
    need(len(set(mixed_C)) == 10
         and all(Q(c, mixed_r) == ratio for c, ratio in zip(mixed_C, mixed_ratios)),
         ("sharp mixed ten-label realization", mixed_C))
    mixed_overlaps = tuple(mixed_overlap_formula(c, mixed_r) for c in mixed_C)
    need(max(mixed_overlaps) == Q(1, 28),
         ("sharp mixed tenth-order overlap", mixed_overlaps))
    print("mixed_oriented_overlap_below_1/28=" + ",".join(word(x) for x in sorted(mixed_below)))
    print("mixed_oriented_overlap_equal_1/28=" + ",".join(word(x) for x in sorted(mixed_equal)))
    print("sharp_mixed_order_stat_r=3575;C=" + ",".join(map(str, mixed_C)) + ";max_overlap=1/28")
    print("q4_hybrid_floor=max(mu(G_C)/2,mu(G_C)-1/8)")
    return literal + mixed_literal


def audit_cross_atlas() -> dict[tuple[int, int], tuple[tuple[Q, Q], ...]]:
    expected = {
        (1, 11): ((Q(6, 77), Q(8, 77)), (Q(69, 77), Q(71, 77))),
        (1, 23): ((Q(6, 161), Q(8, 161)), (Q(20, 161), Q(22, 161)),
                  (Q(139, 161), Q(141, 161)), (Q(153, 161), Q(155, 161))),
        (5, 11): ((Q(6, 35), Q(15, 77)), (Q(62, 77), Q(29, 35))),
        (1, 37): ((Q(6, 259), Q(8, 259)), (Q(20, 259), Q(22, 259)),
                  (Q(34, 259), Q(36, 259)), (Q(223, 259), Q(225, 259)),
                  (Q(237, 259), Q(239, 259)), (Q(251, 259), Q(253, 259))),
        (1, 25): ((Q(6, 175), Q(8, 175)), (Q(4, 35), Q(22, 175)),
                  (Q(153, 175), Q(31, 35)), (Q(167, 175), Q(169, 175))),
    }
    expected_mass = {
        (1, 11): Q(4, 77), (1, 23): Q(8, 161), (5, 11): Q(18, 385),
        (1, 37): Q(12, 259), (1, 25): Q(8, 175),
    }
    atlas = {}
    for ratio in FIVE:
        cells = cross_components(*ratio)
        need(cells == expected[ratio], ("five-ray cells", ratio, cells))
        need(cross_measure(*ratio) == expected_mass[ratio], ("five-ray mass", ratio))
        atlas[ratio] = cells
        beta = max(right - left for left, right in cells)
        print(
            f"ray={ratio[0]}:{ratio[1]};mass={word(expected_mass[ratio])};"
            f"beta={word(beta)};cells=" + ",".join(interval_word(cell) for cell in cells)
        )
    return atlas


def audit_body(
    name: str,
    H: tuple[int, ...],
    q4_C: tuple[int, ...] | None,
    q2_r: int | None,
    expected_mass: Q,
    expected_components: int,
    expected_longest: Q,
    expected_hits: dict[tuple[int, int], int],
    expected_safe_boundaries: dict[tuple[int, int], int],
    atlas: dict[tuple[int, int], tuple[tuple[Q, Q], ...]],
) -> int:
    need(len(H) == 11 and gcd(*H) == 1 and divisor_complete(H), ("residual body sieve", name, H))
    components = safe_components(H)
    mass = sum((right - left for left, right in components), Q())
    longest = max(right - left for left, right in components)
    need((mass, len(components), longest) == (expected_mass, expected_components, expected_longest),
         ("body geometry", name, mass, len(components), longest))
    need(Q(4, 91) <= mass < Q(4, 77), ("five-ray residual band", name, mass))

    if q4_C is not None:
        odd = tuple(h for h in H if h % 2)
        need(len(q4_C) == 10 and gcd(*q4_C) == 1 and len(odd) == 1, ("exact q4 provenance", name))
        r = odd[0]
        need(tuple(sorted(tuple(2 * c for c in q4_C) + (r,))) == H, ("q4 body identity", name))
        need(q4_reconstruct(q4_C, r) == components, ("q4 fibre reconstruction", name))
    if q2_r is not None:
        C = tuple(h for h in H if h != q2_r)
        need(len(C) == 10 and gcd(*C) == 1 and q2_r % 2 == 1 and q2_r % 3,
             ("exact q2 provenance", name, C, q2_r))
        need(q2_reconstruct(C, q2_r) == components, ("q2 intersection reconstruction", name))

    checks = 0
    body_line = (
        f"body={name};H=" + ",".join(map(str, H)) +
        f";mass={word(mass)};components={len(components)};longest={word(longest)}"
    )
    print(body_line)
    for ratio in FIVE:
        cells = atlas[ratio]
        beta = max(right - left for left, right in cells)
        bank = scale_bank(beta, longest)
        need(bank == (1,), ("hostile body singleton bank", name, ratio, bank))
        t = 1
        addresses = tuple(component_address(component, cells, t) for component in components)
        hit_count = sum(address is not None for address in addresses)
        need(hit_count == expected_hits[ratio], ("component hit count", name, ratio, hit_count))
        need(hit_count < len(components), ("not a row counterexample", name, ratio))
        safe_boundary_count = sum(safe(H, point) for point in pulled_boundaries(cells, t))
        need(safe_boundary_count == expected_safe_boundaries[ratio],
             ("endpoint sieve count", name, ratio, safe_boundary_count))
        y, sheet, x, clearance = row_witness(H, *ratio, t, components, cells)
        print(
            f"bank={name}:{ratio[0]}:{ratio[1]};t=1;endpoint_safe={safe_boundary_count};"
            f"trapped_components={hit_count};escaping_components={len(components)-hit_count};"
            f"witness_y={word(y)};sheet={sheet};x={word(x)};clearance={word(clearance)}"
        )
        checks += 1
    return checks


def main() -> None:
    literal = audit_overlap_order_statistic()
    atlas = audit_cross_atlas()

    # q=4 bodies: H=2C union {r}, with gcd(C)=1.  The three bodies together
    # expose a trapped component on every one of the five primitive rays.
    q4_cases = (
        (
            "q4_A",
            (2, 12, 14, 16, 18, 20, 22, 25, 26, 34, 38),
            (1, 6, 7, 8, 9, 10, 11, 13, 17, 19),
            Q(1151237, 25865840), 44, Q(17, 3192),
            {(1, 11): 2, (1, 23): 2, (5, 11): 4, (1, 37): 0, (1, 25): 2},
            {ratio: 0 for ratio in FIVE},
        ),
        (
            "q4_B",
            (2, 6, 8, 14, 20, 22, 23, 26, 32, 34, 36),
            (1, 3, 4, 7, 10, 11, 13, 16, 17, 18),
            Q(174683549, 3945221280), 46, Q(1, 180),
            {(1, 11): 4, (1, 23): 0, (5, 11): 4, (1, 37): 2, (1, 25): 2},
            {(1, 11): 2, (1, 23): 0, (5, 11): 0, (1, 37): 0, (1, 25): 0},
        ),
    )

    # q=2 one-even bodies: H=C union {r}, gcd(C)=1.  Two bodies suffice to
    # expose trapped components on every five-ray type.
    q2_cases = (
        (
            "q2_A",
            (2, 6, 10, 14, 16, 17, 18, 22, 23, 26, 60),
            17,
            Q(12707741, 246576330), 36, Q(23, 3808),
            {(1, 11): 2, (1, 23): 2, (5, 11): 0, (1, 37): 2, (1, 25): 2},
            {ratio: 0 for ratio in FIVE},
        ),
        (
            "q2_B",
            (2, 3, 10, 12, 13, 14, 16, 17, 18, 19, 22),
            13,
            Q(10204829, 203693490), 20, Q(1, 154),
            {(1, 11): 0, (1, 23): 4, (5, 11): 2, (1, 37): 0, (1, 25): 2},
            {(1, 11): 0, (1, 23): 2, (5, 11): 0, (1, 37): 2, (1, 25): 0},
        ),
    )

    bank_checks = 0
    for name, H, C, mass, count, longest, hits, boundaries in q4_cases:
        bank_checks += audit_body(name, H, C, None, mass, count, longest, hits, boundaries, atlas)
    for name, H, r, mass, count, longest, hits, boundaries in q2_cases:
        bank_checks += audit_body(name, H, None, r, mass, count, longest, hits, boundaries, atlas)

    # Every five-ray type really appears as a component trap in each provenance
    # class, but none of the displayed finite banks is a full containment.
    for cases in (q4_cases, q2_cases):
        union = {ratio for case in cases for ratio, count in case[-2].items() if count}
        need(union == set(FIVE), ("hostile ray coverage", union))

    # For each ray there is a trapped component in a body whose every pulled
    # endpoint is already blocked, in each structured provenance class.
    for cases in (q4_cases, q2_cases):
        endpoint_hostile_union = {
            ratio for case in cases for ratio, count in case[-2].items()
            if count and case[-1][ratio] == 0
        }
        need(endpoint_hostile_union == set(FIVE),
             ("endpoint-only hostile coverage", endpoint_hostile_union))

    print(f"status=PASS;literal_overlap_checks={literal};finite_bank_checks={bank_checks};LRC14=OPEN")


if __name__ == "__main__":
    main()
