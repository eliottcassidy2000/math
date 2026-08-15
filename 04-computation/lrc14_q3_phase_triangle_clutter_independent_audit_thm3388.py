#!/usr/bin/env python3
"""Independent exact audit of THM-3388.

This companion deliberately does not import the primary source.
It compares its affine-gap criterion with a direct rational endpoint/mid-cell
event sweep, rebuilds the literal clutter and all 2,793 q=3 body rows, and
checks the ternary harmonic formulas from their generating functions.
"""

from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, permutations
from math import comb, factorial, gcd
from pathlib import Path


VERTICES = (1, 2, 4, 5, 7, 8, 10, 11, 13, 14)
CORES = (1, 2, 3, 4)


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def danger(speed, time):
    residue = (speed * time) % 1
    return min(residue, 1 - residue) < F(1, 14)


def transverse_danger(speed, time, sheet):
    return danger(speed, time + F(sheet, 3))


def core_danger(clock, time):
    return danger(3 * clock, time)


def raw_boundaries(transverse, core=()):
    events = {F(0)}
    for speed in transverse:
        for sheet in range(3):
            for tooth in range(speed):
                centre = F(tooth, speed) - F(sheet, 3)
                radius = F(1, 14 * speed)
                events.add((centre - radius) % 1)
                events.add((centre + radius) % 1)
    for clock in core:
        speed = 3 * clock
        for tooth in range(speed):
            centre = F(tooth, speed)
            radius = F(1, 14 * speed)
            events.add((centre - radius) % 1)
            events.add((centre + radius) % 1)
    return sorted(events)


def boundary_events(transverse, core=()):
    ordered = raw_boundaries(transverse, core)
    samples = list(ordered)
    for index, left in enumerate(ordered):
        right = ordered[(index + 1) % len(ordered)]
        if index + 1 == len(ordered):
            right += 1
        samples.append(((left + right) / 2) % 1)
    return samples


def open_cover_components(transverse):
    ordered = raw_boundaries(transverse)
    positive_cells = []
    for index, left in enumerate(ordered):
        right = ordered[(index + 1) % len(ordered)]
        if index + 1 == len(ordered):
            right += 1
        if transverse_cover(transverse, ((left + right) / 2) % 1):
            positive_cells.append((left, right))
    components = []
    for left, right in positive_cells:
        if components and components[-1][1] == left:
            components[-1] = (components[-1][0], right)
        else:
            components.append((left, right))
    require(not components or components[-1][1] <= 1, "unexpected wrap component")
    return tuple(components)


def transverse_cover(transverse, time):
    return all(
        any(transverse_danger(speed, time, sheet) for speed in transverse)
        for sheet in range(3)
    )


def has_core_safe_transverse_cover(transverse, core=()):
    return any(
        transverse_cover(transverse, time)
        and not any(core_danger(clock, time) for clock in core)
        for time in boundary_events(transverse, core)
    )


def gap_values(left, right):
    modulus = 3 * gcd(left, right)
    residue = (left * right) % modulus
    bound = (3 * (left + right) - 1) // 14
    return tuple(
        value
        for value in range(-bound, bound + 1)
        if (value - residue) % modulus == 0
    )


def affine_witness(order):
    left, middle, right = order
    third_values = set(gap_values(right, left))
    for p in gap_values(left, middle):
        for q in gap_values(middle, right):
            numerator = -(right * p + left * q)
            if numerator % middle:
                continue
            r = numerator // middle
            if r in third_values:
                return p, q, r
    return None


def is_gap(left, right, value):
    return (
        (value - left * right) % (3 * gcd(left, right)) == 0
        and 14 * abs(value) < 3 * (left + right)
    )


def brute_glue(order, gaps):
    """Search the centre congruences directly, without generalized CRT code."""
    u, v, w = order
    p, q, r = gaps
    A = (p - u * v) // 3
    B = (q - v * w) // 3
    C = (r + 2 * w * u) // 3
    for b in range(v):
        if (A + u * b) % v or (w * b - B) % v:
            continue
        a = (A + u * b) // v
        c = (w * b - B) // v
        if u * c - w * a == C:
            return a, b, c
    return None


def prime_divisors(value):
    primes = []
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            primes.append(divisor)
            while value % divisor == 0:
                value //= divisor
        divisor += 1
    if value > 1:
        primes.append(value)
    return tuple(primes)


def valuation(value, prime):
    answer = 0
    while value % prime == 0:
        value //= prime
        answer += 1
    return answer


def local_gluing_checks(order, gaps):
    """Verify the least-valuation p-local construction modulo a deep power."""
    u, v, w = order
    p_gap, q_gap, r_gap = gaps
    A = (p_gap - u * v) // 3
    B = (q_gap - v * w) // 3
    C = (r_gap + 2 * w * u) // 3
    checks = 0
    for prime in (*prime_divisors(u * v * w), 31):
        height = min(valuation(u, prime), valuation(v, prime), valuation(w, prime))
        factor = prime**height
        uu, vv, ww = u // factor, v // factor, w // factor
        require(A % factor == B % factor == C % factor == 0, ("local divisibility", order, gaps, prime))
        AA, BB, CC = A // factor, B // factor, C // factor
        modulus = prime**6
        if uu % prime:
            a = 0
            b = (-AA * pow(uu, -1, modulus)) % modulus
            c = (CC * pow(uu, -1, modulus)) % modulus
        elif vv % prime:
            b = 0
            a = (AA * pow(vv, -1, modulus)) % modulus
            c = (-BB * pow(vv, -1, modulus)) % modulus
        else:
            require(ww % prime != 0, ("no local unit", order, prime))
            c = 0
            b = (BB * pow(ww, -1, modulus)) % modulus
            a = (-CC * pow(ww, -1, modulus)) % modulus
        require(
            (
                (vv * a - uu * b - AA) % modulus,
                (ww * b - vv * c - BB) % modulus,
                (uu * c - ww * a - CC) % modulus,
            )
            == (0, 0, 0),
            ("p-local glue", order, gaps, prime),
        )
        checks += 1
    return checks


def independent(subset, edges):
    chosen = set(subset)
    return not any(edge <= chosen for edge in edges)


def multinomial(depth, i, j):
    k = depth - i - j
    return factorial(depth) // (factorial(i) * factorial(j) * factorial(k))


def main():
    semantic = sha256()

    # The affine pair test and the integral cycle-gluing lemma are checked on
    # a wider low-speed universe before the literal [14] atlas is touched.
    pair_checks = 0
    pair_pool = tuple(value for value in range(1, 200) if value % 3)
    for left, right in combinations(pair_pool, 2):
        require(
            bool(gap_values(left, right))
            == (3 * (left + right) > 14 * gcd(left, right)),
            ("pair threshold", left, right),
        )
        pair_checks += 1

    gluing_tuples = 0
    local_checks = 0
    gluing_pool = tuple(value for value in range(1, 31) if value % 3)
    for order in combinations(gluing_pool, 3):
        u, v, w = order
        for p in gap_values(u, v):
            for q in gap_values(v, w):
                for r in gap_values(w, u):
                    closed = w * p + u * q + v * r == 0
                    glued = brute_glue(order, (p, q, r)) is not None
                    require(closed == glued, ("cycle gluing", order, (p, q, r)))
                    if closed:
                        local_checks += local_gluing_checks(order, (p, q, r))
                    gluing_tuples += 1

    # Direct geometry versus affine closure on the full literal triple pool,
    # plus every ordering gauge and direct centre-congruence reconstruction.
    literal_edges = []
    ordered_checks = 0
    glued_checks = 0
    for triple in combinations(VERTICES, 3):
        geometric = has_core_safe_transverse_cover(triple)
        decisions = []
        for order in permutations(triple):
            witness = affine_witness(order)
            decisions.append(witness is not None)
            if witness is not None:
                require(brute_glue(order, witness) is not None, ("literal glue", order, witness))
                glued_checks += 1
            ordered_checks += 1
        require(set(decisions) == {geometric}, ("literal geometry", triple, geometric, decisions))
        if geometric:
            literal_edges.append(triple)

    edges = tuple(frozenset(edge) for edge in literal_edges)
    pair_triangles = tuple(
        triple
        for triple in combinations(VERTICES, 3)
        if all(gap_values(left, right) for left, right in combinations(triple, 2))
    )
    profile = tuple(
        sum(independent(subset, edges) for subset in combinations(VERTICES, size))
        for size in range(len(VERTICES) + 1)
    )
    five_sets = tuple(
        subset
        for subset in combinations(VERTICES, 5)
        if independent(subset, edges)
    )

    require(len(literal_edges) == 48, ("edge count", len(literal_edges)))
    require(len(pair_triangles) == 82, ("pair triangle count", len(pair_triangles)))
    require(profile == (1, 10, 45, 72, 38, 6, 0, 0, 0, 0, 0), ("profile", profile))
    require(len(five_sets) == 6, ("five sets", five_sets))
    require((1, 4, 7) in pair_triangles and (1, 4, 7) not in literal_edges, "hostile 147")
    require(affine_witness((1, 4, 5)) == (1, -1, -1), "positive 145")

    single_core_containments = tuple(
        (edge, clock)
        for edge in literal_edges
        for clock in CORES
        if not has_core_safe_transverse_cover(edge, (clock,))
    )
    require(
        single_core_containments == (((8, 11, 13), 2),),
        ("single-core containments", single_core_containments),
    )
    rescue_cells = open_cover_components((8, 11, 13))
    expected_cell_starts = (583, 648, 1815, 1880, 3047, 3112)
    require(
        rescue_cells
        == tuple((F(start, 3696), F(start + 1, 3696)) for start in expected_cell_starts),
        ("rescue cells", rescue_cells),
    )
    containment_margins = []
    for left, right in rescue_cells:
        candidates = []
        for tooth in range(6):
            centre = F(tooth, 6)
            low, high = centre - F(1, 84), centre + F(1, 84)
            if low < left and right < high:
                candidates.append(min(left - low, high - right))
        require(len(candidates) == 1, ("core-2 containment", left, right, candidates))
        containment_margins.extend(candidates)
    min_core2_margin = min(containment_margins)
    require(min_core2_margin == F(1, 336), ("core-2 margin", min_core2_margin))

    # Direct rational geometry for every literal body row; the global count is
    # also recomputed independently from hypergraph independence.
    candidates = exact_rows = global_rows = 0
    rescues = []
    for core_size in range(1, 5):
        transverse_size = 6 - core_size
        for core in combinations(CORES, core_size):
            for transverse in combinations(VERTICES, transverse_size):
                candidates += 1
                globally_safe = independent(transverse, edges)
                exact = not has_core_safe_transverse_cover(transverse, core)
                global_rows += globally_safe
                exact_rows += exact
                if globally_safe:
                    require(exact, ("global safe but pointwise leaking", core, transverse))
                elif exact:
                    rescues.append((core, transverse))

    require(candidates == 2793, ("candidates", candidates))
    require(global_rows == 585, ("global rows", global_rows))
    require(exact_rows == 588, ("exact rows", exact_rows))
    require(
        tuple(rescues)
        == (
            ((1, 2, 3), (8, 11, 13)),
            ((1, 2, 4), (8, 11, 13)),
            ((2, 3, 4), (8, 11, 13)),
        ),
        ("rescues", rescues),
    )
    require(
        global_rows == 4 * profile[5] + 6 * profile[4] + 4 * profile[3] + profile[2],
        "independence decomposition",
    )

    # Independent generating-function reconstruction of the ternary shells.
    root_mass = sum((F(1, value) for value in (1, 4, 5)), F(0))
    support = []
    word = []
    shell_checks = 0
    for depth in range(13):
        support_mass = F(0)
        word_mass = F(0)
        scales = set()
        word_count = 0
        for i in range(depth + 1):
            for j in range(depth - i + 1):
                k = depth - i - j
                scale = 7**i * 11**j * 13**k
                require(scale not in scales, ("scale collision", depth, i, j))
                scales.add(scale)
                multiplicity = multinomial(depth, i, j)
                word_count += multiplicity
                support_mass += root_mass / scale
                word_mass += multiplicity * root_mass / scale
                sign = 1 if scale % 3 == 1 else -1
                scaled_gaps = (sign * scale, -sign * scale, -sign * scale)
                require(
                    all(
                        is_gap(left, right, gap)
                        for (left, right), gap in zip(
                            ((scale, 4 * scale), (4 * scale, 5 * scale), (5 * scale, scale)),
                            scaled_gaps,
                        )
                    ),
                    ("dilated gaps", depth, i, j),
                )
                require(
                    5 * scale * scaled_gaps[0]
                    + scale * scaled_gaps[1]
                    + 4 * scale * scaled_gaps[2]
                    == 0,
                    ("dilated circulation", depth, i, j),
                )
                shell_checks += 1
        require(len(scales) == comb(depth + 2, 2), ("scale count", depth))
        require(word_count == 3**depth, ("word count", depth))
        support.append(support_mass)
        word.append(word_mass)
    for depth in range(3, len(support)):
        require(
            1001 * support[depth]
            == 311 * support[depth - 1]
            - 31 * support[depth - 2]
            + support[depth - 3],
            ("support recurrence", depth),
        )
    for depth in range(len(word) - 1):
        require(1001 * word[depth + 1] == 311 * word[depth], ("word recurrence", depth))
    orbit_mass = root_mass * F(7, 6) * F(11, 10) * F(13, 12)
    require(root_mass == F(29, 20), ("root mass", root_mass))
    require(orbit_mass == F(29029, 14400), ("orbit mass", orbit_mass))

    payload = (
        tuple(literal_edges), tuple(pair_triangles), profile, five_sets,
        single_core_containments, rescue_cells, min_core2_margin,
        candidates, global_rows, exact_rows, tuple(rescues),
        tuple(support), tuple(word), orbit_mass,
    )
    semantic.update(repr(payload).encode("ascii"))

    print("THM-3388 INDEPENDENT Q3 PHASE-TRIANGLE AUDIT")
    source_hash = sha256(Path(__file__).read_bytes().replace(bytes((13, 10)), bytes((10,)))).hexdigest()
    print(f"source_sha256_lf={source_hash}")
    print(f"pair_checks={pair_checks};cycle_gluing_tuples={gluing_tuples};p_local_checks={local_checks}")
    print(f"literal_ordered_checks={ordered_checks};glued_witnesses={glued_checks}")
    print(f"hostile_147_gap_sets={(gap_values(1,4), gap_values(4,7), gap_values(7,1))};closed={affine_witness((1,4,7))}")
    print(f"positive_145_gap_sets={(gap_values(1,4), gap_values(4,5), gap_values(5,1))};closed={affine_witness((1,4,5))}")
    print(f"cover_edges={len(literal_edges)};pair_triangles={len(pair_triangles)};false_positives={len(pair_triangles)-len(literal_edges)}")
    print(f"independence_profile={profile};maximal_five_sets={five_sets}")
    print(f"single_core_containments={single_core_containments}")
    print(
        "core2_rescue_cells_den3696="
        f"{tuple((int(3696 * left), int(3696 * right)) for left, right in rescue_cells)};"
        f"min_margin={min_core2_margin}"
    )
    print(f"body_candidates={candidates};global_rows={global_rows};exact_rows={exact_rows};core_rescues={tuple(rescues)}")
    print(f"ternary_shells=0..12;checks={shell_checks};root_mass={root_mass};orbit_mass={orbit_mass}")
    print(f"audit_semantic_sha256={semantic.hexdigest()}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
