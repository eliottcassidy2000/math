#!/usr/bin/env python3
"""Clean-room exact referee for the proposed THM-4450 scratch report."""

from __future__ import annotations

from fractions import Fraction as Q
from functools import cache
from math import gcd


DELTA = Q(1, 14)
Interval = tuple[Q, Q]
CHECKS = 0


def need(condition: bool, detail: object) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(detail)


def frac(x: Q) -> Q:
    return x % 1


def gap(x: Q) -> Q:
    x %= 1
    return min(x, 1 - x)


def bernoulli2(x: Q) -> Q:
    x = frac(x)
    return x * x - x + Q(1, 6)


def merge_ae(intervals: list[Interval]) -> list[Interval]:
    out: list[list[Q]] = []
    for left, right in sorted(intervals):
        if left >= right:
            continue
        if not out or left > out[-1][1]:
            out.append([left, right])
        elif right > out[-1][1]:
            out[-1][1] = right
    return [(left, right) for left, right in out]


def intersect_ae(a: list[Interval], b: list[Interval]) -> list[Interval]:
    i = j = 0
    out: list[Interval] = []
    while i < len(a) and j < len(b):
        left = max(a[i][0], b[j][0])
        right = min(a[i][1], b[j][1])
        if left < right:
            out.append((left, right))
        if a[i][1] < b[j][1]:
            i += 1
        elif b[j][1] < a[i][1]:
            j += 1
        else:
            i += 1
            j += 1
    return merge_ae(out)


def measure(intervals: list[Interval]) -> Q:
    return sum((right - left for left, right in intervals), Q(0))


@cache
def comb(n: int, radius: Q = DELTA) -> tuple[Interval, ...]:
    arcs: list[Interval] = []
    physical_radius = radius / n
    for k in range(n):
        left = Q(k, n) - physical_radius
        right = Q(k, n) + physical_radius
        if left < 0:
            arcs.extend(((Q(0), right), (left + 1, Q(1))))
        elif right > 1:
            arcs.extend(((left, Q(1)), (Q(0), right - 1)))
        else:
            arcs.append((left, right))
    return tuple(merge_ae(arcs))


def literal_overlap(c: int, r: int, mixed: bool = False) -> Q:
    return measure(intersect_ae(list(comb(c)), list(comb(r, Q(1, 7) if mixed else DELTA))))


def overlap_formula(p: int, q: int, mixed: bool = False) -> Q:
    need(p > 0 and q > 0 and gcd(p, q) == 1, (p, q))
    if mixed:
        return Q(2, 49) + (
            bernoulli2(Q(q - 2 * p, 14)) - bernoulli2(Q(q + 2 * p, 14))
        ) / (p * q)
    return Q(1, 49) + (
        bernoulli2(Q(q - p, 14)) - bernoulli2(Q(q + p, 14))
    ) / (p * q)


def cross_pair_formula(p: int, q: int) -> Q:
    """Two-sheet cross-comb mass, distinct from equal-radius overlap."""
    alpha = (p + q) // 2
    beta = (q - p) // 2
    d = (alpha + 3) % 7 - 3
    e = (beta + 3) % 7 - 3
    return Q(2, 49) * (1 + Q(e * e - d * d, p * q))


def danger(n: int, x: Q, radius: Q = DELTA) -> bool:
    return gap(Q(n) * x) < radius


def safe(body: tuple[int, ...], x: Q) -> bool:
    return all(not danger(n, x) for n in body)


def comb_walls(n: int, radius: Q = DELTA) -> set[Q]:
    walls = {Q(0), Q(1)}
    for k in range(n):
        for sign in (-1, 1):
            walls.add(frac(Q(k, n) + sign * radius / n))
    return walls


def safe_components(body: tuple[int, ...]) -> list[tuple[Q, Q]]:
    """Complete closed components, including isolated safe points."""
    walls = sorted(set().union(*(comb_walls(n) for n in body)))
    cell_safe = [safe(body, (left + right) / 2) for left, right in zip(walls, walls[1:])]
    comps: list[tuple[Q, Q]] = []
    start: Q | None = None
    for i, live in enumerate(cell_safe):
        if live and start is None:
            start = walls[i]
            need(safe(body, start), (body, start, "closed left endpoint"))
        if start is not None and (not live or i == len(cell_safe) - 1):
            right = walls[i] if not live else walls[i + 1]
            need(safe(body, right), (body, right, "closed right endpoint"))
            comps.append((start, right))
            start = None
    for i in range(1, len(walls) - 1):
        if safe(body, walls[i]) and not cell_safe[i - 1] and not cell_safe[i]:
            comps.append((walls[i], walls[i]))
    comps.sort()
    need(not safe(body, Q(0)), (body, "unexpected circular join"))
    return comps


def pair_failure(p: int, q: int, y: Q) -> bool:
    """Quotient-y failure: both physical lifts y/2,(y+1)/2 are killed."""
    return all(
        danger(p, (y + j) / 2) or danger(q, (y + j) / 2)
        for j in (0, 1)
    )


def shifted_walls(n: int) -> set[Q]:
    return {frac(x - Q(1, 2)) for x in comb_walls(n)} | {Q(0), Q(1)}


def quotient_pair_walls(n: int) -> set[Q]:
    """Images under y=2x of all physical danger walls for tail n."""
    return {frac(2 * x) for x in comb_walls(n)} | {Q(0), Q(1)}


def pair_components(p: int, q: int) -> list[Interval]:
    """Actual strict open components, retaining a deleted touching wall."""
    need(p % 2 == q % 2 == 1, (p, q, "odd pair"))
    walls = sorted(quotient_pair_walls(p) | quotient_pair_walls(q))
    cells = [pair_failure(p, q, (left + right) / 2) for left, right in zip(walls, walls[1:])]
    comps: list[Interval] = []
    start: Q | None = None
    for i, live in enumerate(cells):
        if live:
            if start is None:
                start = walls[i]
            elif i > 0 and cells[i - 1] and not pair_failure(p, q, walls[i]):
                comps.append((start, walls[i]))
                start = walls[i]
        elif start is not None:
            comps.append((start, walls[i]))
            start = None
    if start is not None:
        comps.append((start, Q(1)))
    need(not pair_failure(p, q, Q(0)), (p, q, "unexpected circular join"))
    return comps


def pulled_cells(cells: list[Interval], t: int) -> list[Interval]:
    return sorted(((Q(n) + left) / t, (Q(n) + right) / t) for n in range(t) for left, right in cells)


def decoder(body: tuple[int, ...], p: int, q: int, t: int) -> tuple[bool, int, int]:
    components = safe_components(body)
    cells = pulled_cells(pair_components(p, q), t)
    trapped = sum(any(u < left and right < v for u, v in cells) for left, right in components)
    endpoint_safe = sum(safe(body, x) for cell in cells for x in cell)
    return trapped == len(components), trapped, endpoint_safe


def literal_containment(body: tuple[int, ...], p: int, q: int, t: int) -> bool:
    walls = set().union(*(comb_walls(n) for n in body))
    walls |= quotient_pair_walls(t * p) | quotient_pair_walls(t * q)
    ordered = sorted(walls)
    for x in ordered:
        if safe(body, x) and not pair_failure(t * p, t * q, x):
            return False
    for left, right in zip(ordered, ordered[1:]):
        x = (left + right) / 2
        if safe(body, x) and not pair_failure(t * p, t * q, x):
            return False
    return True


def safe_mass(body: tuple[int, ...]) -> Q:
    return sum((right - left for left, right in safe_components(body)), Q(0))


def check_overlap_hierarchy() -> None:
    classes = [
        (
            "g1",
            Q(1, 63),
            55,
            lambda q: gcd(q, 6) == 1,
            {Q(1, 13), Q(1, 11), Q(2, 11), Q(3, 11), Q(10), Q(11), Q(12), Q(13)},
            {Q(9, 5), Q(9), Q(27)},
        ),
        (
            "g2",
            Q(1, 70),
            40,
            lambda q: q % 3 != 0,
            {Q(1, 13), Q(1, 11), Q(2, 11), Q(3, 11), Q(11, 2), Q(11), Q(12), Q(13)},
            {Q(1, 10), Q(3, 10), Q(10)},
        ),
        (
            "g3",
            Q(1, 70),
            40,
            lambda q: q % 2 == 1,
            {Q(1, 13), Q(1, 11), Q(2, 11), Q(3, 11), Q(11, 3), Q(11), Q(12), Q(13)},
            {Q(10, 3), Q(10)},
        ),
        (
            "g6",
            Q(1, 77),
            33,
            lambda q: True,
            {Q(1, 13), Q(1, 12), Q(12), Q(13)},
            {Q(1, 11), Q(2, 11), Q(3, 11), Q(11, 3), Q(11, 2), Q(11)},
        ),
    ]
    print("EQUAL_RADIUS_ATLASES")
    literal = 0
    for label, lam, cutoff, q_allowed, expected_below, expected_equal in classes:
        below: set[Q] = set()
        equal: set[Q] = set()
        for q in range(1, cutoff + 1):
            if not q_allowed(q):
                continue
            for p in range(1, cutoff // q + 1):
                if gcd(p, q) != 1:
                    continue
                value = overlap_formula(p, q)
                need(literal_overlap(p, q) == value, (label, p, q, value))
                literal += 1
                if value < lam:
                    below.add(Q(p, q))
                elif value == lam:
                    equal.add(Q(p, q))
        need(below == expected_below, (label, "below", below ^ expected_below))
        need(equal == expected_equal, (label, "equal", equal ^ expected_equal))
        # Strict lower bound above the product atlas.
        need(Q(1, 49) - Q(1, 4 * (cutoff + 1)) > lam, (label, "cutoff"))
        print(f"{label}: lambda={lam} cutoff={cutoff} below={sorted(below)} equal={sorted(equal)}")

    controls = [
        ("g1", 715, (55, 9295, 8580, 65, 130, 195, 7865, 7150, 1287, 6435), Q(1, 63)),
        ("g2", 1430, (110, 130, 260, 390, 7865, 15730, 17160, 18590, 143, 429), Q(1, 70)),
        ("g3", 429, (33, 39, 78, 117, 1573, 4719, 5148, 5577, 1430, 4290), Q(1, 70)),
        ("g6", 1716, (132, 143, 20592, 22308, 156, 312, 468, 6292, 9438, 18876), Q(1, 77)),
    ]
    print("TENTH_ORDER_CONTROLS")
    for label, r, body, lam in controls:
        values = [overlap_formula(c // gcd(c, r), r // gcd(c, r)) for c in body]
        need(len(set(body)) == 10 and r not in body and gcd(*body) == 1, (label, "control type"))
        need(max(values) == lam, (label, max(values), lam))
        print(f"{label}: r={r} max={max(values)} equality_count={sum(v == lam for v in values)}")
    print(f"equal_radius_literal_checks={literal}")


def check_q2_arithmetic() -> None:
    lambdas = {1: Q(1, 63), 2: Q(1, 70), 3: Q(1, 70), 6: Q(1, 77)}
    losses = {g: Q(1, 7) - lam for g, lam in lambdas.items()}
    need(losses == {1: Q(8, 63), 2: Q(9, 70), 3: Q(9, 70), 6: Q(10, 77)}, losses)
    expected = {
        Q(4, 63): {1: Q(4, 21), 2: Q(121, 630), 3: Q(121, 630), 6: Q(134, 693)},
        Q(4, 77): {1: Q(124, 693), 2: Q(139, 770), 3: Q(139, 770), 6: Q(2, 11)},
        Q(4, 91): {1: Q(20, 117), 2: Q(157, 910), 3: Q(157, 910), 6: Q(174, 1001)},
    }
    for cap, table in expected.items():
        need({g: loss + cap for g, loss in losses.items()} == table, (cap, "entry table"))
    rays = [(1, 11), (1, 23), (5, 11), (1, 37), (1, 25)]
    expected_specific = [Q(124, 693), Q(256, 1449), Q(86, 495), Q(404, 2331), Q(272, 1575)]
    actual = [Q(8, 63) + cross_pair_formula(p, q) for p, q in rays]
    need(actual == expected_specific, (actual, expected_specific))
    print(f"Q2_LOSSES={losses}")
    print(f"Q2_RATIO_SPECIFIC={actual}")


def check_mixed_hierarchy() -> None:
    lam = Q(1, 28)
    cutoff = 49
    below: set[Q] = set()
    equal: set[Q] = set()
    literal = 0
    for q in range(1, cutoff + 1):
        if gcd(q, 6) != 1:
            continue
        for p in range(1, cutoff // q + 1):
            if gcd(p, q) != 1:
                continue
            value = overlap_formula(p, q, mixed=True)
            need(literal_overlap(p, q, mixed=True) == value, (p, q, "mixed formula"))
            literal += 1
            if value < lam:
                below.add(Q(p, q))
            elif value == lam:
                equal.add(Q(p, q))
    expected_below = {Q(1, 25), Q(1, 13), Q(1, 11), Q(2, 11), Q(5), Q(6), Q(13)}
    expected_equal = {Q(4, 5), Q(4), Q(12), Q(20)}
    need(below == expected_below, below ^ expected_below)
    need(equal == expected_equal, equal ^ expected_equal)
    need(Q(2, 49) - Q(1, 4 * (cutoff + 1)) > lam, "mixed cutoff")

    r = 3575
    body = (21450, 325, 17875, 650, 275, 46475, 143, 2860, 14300, 42900)
    values = [overlap_formula(c // gcd(c, r), r // gcd(c, r), mixed=True) for c in body]
    need(len(set(body)) == 10 and r not in body and gcd(*body) == 1, "mixed control type")
    need(max(values) == lam, (max(values), lam))
    need(Q(2, 7) - lam == Q(1, 4), "mixed uncovered cap")
    print(f"MIXED_ATLAS below={sorted(below)} equal={sorted(equal)} literal_checks={literal}")
    print(f"MIXED_CONTROL max={max(values)} equality_count={sum(v == lam for v in values)}")


def check_decoder_and_structured_bodies() -> None:
    rays = [(1, 11), (1, 23), (5, 11), (1, 37), (1, 25)]
    expected_cells = {
        (1, 11): [(Q(6, 77), Q(8, 77)), (Q(69, 77), Q(71, 77))],
        (1, 23): [(Q(6, 161), Q(8, 161)), (Q(20, 161), Q(22, 161)), (Q(139, 161), Q(141, 161)), (Q(153, 161), Q(155, 161))],
        (5, 11): [(Q(6, 35), Q(15, 77)), (Q(62, 77), Q(29, 35))],
        (1, 37): [(Q(6, 259), Q(8, 259)), (Q(20, 259), Q(22, 259)), (Q(34, 259), Q(36, 259)), (Q(223, 259), Q(225, 259)), (Q(237, 259), Q(239, 259)), (Q(251, 259), Q(253, 259))],
        (1, 25): [(Q(6, 175), Q(8, 175)), (Q(4, 35), Q(22, 175)), (Q(153, 175), Q(31, 35)), (Q(167, 175), Q(169, 175))],
    }
    print("PAIR_CELL_DECODER")
    for ray in rays:
        cells = pair_components(*ray)
        need(cells == expected_cells[ray], (ray, cells, expected_cells[ray]))
        for t in (1, 5, 7):
            need(pair_components(t * ray[0], t * ray[1]) == pulled_cells(cells, t), (ray, t, "pullback"))
        print(f"ray={ray} cells={cells} beta={max(b-a for a,b in cells)}")

    bodies = [
        (
            "q4_A",
            (2, 12, 14, 16, 18, 20, 22, 25, 26, 34, 38),
            "q4",
            (1, 6, 7, 8, 9, 10, 11, 13, 17, 19),
            25,
            Q(1151237, 25865840),
            44,
            Q(17, 3192),
            [2, 2, 4, 0, 2],
            [0, 0, 0, 0, 0],
        ),
        (
            "q4_B",
            (2, 6, 8, 14, 20, 22, 23, 26, 32, 34, 36),
            "q4",
            (1, 3, 4, 7, 10, 11, 13, 16, 17, 18),
            23,
            Q(174683549, 3945221280),
            46,
            Q(1, 180),
            [4, 0, 4, 2, 2],
            [2, 0, 0, 0, 0],
        ),
        (
            "q2_A",
            (2, 6, 10, 14, 16, 17, 18, 22, 23, 26, 60),
            "q2",
            (2, 6, 10, 14, 16, 18, 22, 23, 26, 60),
            17,
            Q(12707741, 246576330),
            36,
            Q(23, 3808),
            [2, 2, 0, 2, 2],
            [0, 0, 0, 0, 0],
        ),
        (
            "q2_B",
            (2, 3, 10, 12, 13, 14, 16, 17, 18, 19, 22),
            "q2",
            (2, 3, 10, 12, 14, 16, 17, 18, 19, 22),
            13,
            Q(10204829, 203693490),
            20,
            Q(1, 154),
            [0, 4, 2, 0, 2],
            [0, 2, 0, 2, 0],
        ),
    ]
    print("STRUCTURED_BODY_DECODER")
    for label, body, provenance, base, r, expected_mass, expected_count, expected_longest, expected_trapped, expected_endpoints in bodies:
        comps = safe_components(body)
        mass = safe_mass(body)
        longest = max(right - left for left, right in comps)
        need((mass, len(comps), longest) == (expected_mass, expected_count, expected_longest), (label, mass, len(comps), longest))
        if provenance == "q2":
            need(set(body) == set(base) | {r}, (label, "q2 provenance"))
            overlap = measure(intersect_ae(list(comb(r)), merge_ae([arc for n in base for arc in comb(n)])))
            need(mass == safe_mass(base) - Q(1, 7) + overlap, (label, "q2 identity"))
        else:
            need(set(body) == {2 * n for n in base} | {r}, (label, "q4 provenance"))
            base_comps = safe_components(base)
            wide = list(comb(r, Q(1, 7)))
            gc_inter_er = measure(intersect_ae(base_comps, wide))
            need(mass == safe_mass(base) - gc_inter_er / 2, (label, "q4 fibre identity"))
        trapped = []
        endpoints = []
        for p, q in rays:
            decoded, n_trapped, n_endpoint = decoder(body, p, q, 1)
            literal = literal_containment(body, p, q, 1)
            need(decoded == literal, (label, p, q, decoded, literal))
            need(not decoded, (label, p, q, "hostile became containment"))
            trapped.append(n_trapped)
            endpoints.append(n_endpoint)
        need(trapped == expected_trapped, (label, "trapped", trapped, expected_trapped))
        need(endpoints == expected_endpoints, (label, "endpoints", endpoints, expected_endpoints))
        print(f"{label}: mass={mass} components={len(comps)} longest={longest} trapped={trapped} endpoints={endpoints}")


def main() -> None:
    print("THM4450_CLEANROOM_REFEREE")
    print("STATUS=INDEPENDENT_EXACT_AUDIT")
    check_overlap_hierarchy()
    check_q2_arithmetic()
    check_mixed_hierarchy()
    check_decoder_and_structured_bodies()
    print(f"checks={CHECKS};status=PASS;LRC14=OPEN")


if __name__ == "__main__":
    main()
