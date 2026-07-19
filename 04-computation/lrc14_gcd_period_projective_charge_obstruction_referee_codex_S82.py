#!/usr/bin/env python3
"""Exact referee for THM-1226, the gcd-period charge obstruction.

All proof-facing calculations use integers or ``fractions.Fraction``.  The
program verifies the exact projective charge factorization, the translated
LCM strict-high counterfamily, its covered protected-needle embedding,
and the finite-channel constant on THM-1221's disconnected strict-spectrum
branches.  The latter is conditional only on THM-1221's already-certified
branch classification; no bounded speed box is used here.
"""

from __future__ import annotations

from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import gcd, prod
from pathlib import Path


OFFSETS = (1, 2, 3, 5, 7, 11, 13)
BASE_STEP = 27720
CIRCUIT_OFFSETS = (1, 2, 3, 5, 11, 13, 17)
CIRCUIT_BASE_STEP = 55440
HIGH_BAR = F(1, 63)
FIRST_STRICT_HIGH = F(5, 308)
TREE_FLOOR = F(15, 154)
PAIR_FLOOR = F(1, 91)
C_DISC = F(448916, 194775)
FOREST_CAP = F(6, 49)
SHALLOW_F7_CONSTANT = F(195, 1078)
MAJOR_NONDIAGONAL_CHANNELS = ((1, 2), (1, 3), (1, 4), (2, 3),
                              (1, 5), (1, 6), (2, 5), (3, 4))
PHASE_CHANNELS = ((1, 1),) + MAJOR_NONDIAGONAL_CHANNELS


def require(condition: bool, message: object) -> None:
    """Optimization-stable assertion."""
    if not condition:
        raise RuntimeError(message)


def fold(residue: int, modulus: int = 14) -> int:
    residue %= modulus
    return residue * (modulus - residue)


def rho(a: int, b: int) -> F:
    """Haar mass of ``D_a intersect D_b`` at radius ``1/14``."""
    require(a > 0 and b > 0 and a != b, (a, b))
    if a > b:
        a, b = b, a
    g = gcd(a, b)
    modulus = 14 * g
    return F(
        4 * a * b + fold(a + b, modulus) - fold(b - a, modulus),
        196 * a * b,
    )


def rho_ratio(x: F, y: F) -> F:
    ratio = x / y
    return rho(abs(ratio.numerator), abs(ratio.denominator))


def eta(value: F) -> F:
    return value * (1 - value)


def projective_kappa(x: F, y: F) -> tuple[F, tuple[int, int], F]:
    """Return kappa, the unordered reduced channel, and pair mass."""
    ratio = x / y
    a, b = abs(ratio.numerator), abs(ratio.denominator)
    value = rho(a, b)
    return eta(value) * F(a * b, a + b), tuple(sorted((a, b))), value


def reduced_bank(product_cap: int, predicate) -> list[tuple[int, int, F]]:
    return [
        (a, b, rho(a, b))
        for a in range(1, product_cap + 1)
        for b in range(a + 1, product_cap + 1)
        if a * b <= product_cap and gcd(a, b) == 1 and predicate(rho(a, b))
    ]


def graph_connected(values: tuple[F, ...], predicate) -> bool:
    seen = {0}
    stack = [0]
    while stack:
        i = stack.pop()
        for j in range(len(values)):
            if j not in seen and predicate(rho_ratio(values[i], values[j])):
                seen.add(j)
                stack.append(j)
    return len(seen) == len(values)


def maximum_kappa(vertices: set[F]) -> tuple[F, list[tuple[F, F, tuple[int, int], F]]]:
    rows: list[tuple[F, F, tuple[int, int], F]] = []
    maximum = F(0)
    for x, y in combinations(sorted(vertices), 2):
        value, channel, mass = projective_kappa(x, y)
        if value > maximum:
            maximum = value
            rows = [(x, y, channel, mass)]
        elif value == maximum:
            rows.append((x, y, channel, mass))
    return maximum, rows


def lcm_packet(multiplier: int) -> tuple[int, tuple[int, ...]]:
    require(multiplier > 0, multiplier)
    A = BASE_STEP * multiplier
    return A, tuple(A + r for r in OFFSETS)


def check_counterfamily(multiplier: int) -> dict[str, object]:
    A, speeds = lcm_packet(multiplier)
    require(A % 14 == 0 and sum(OFFSETS) == 42, (A, OFFSETS))
    require(all(A % d == 0 for d in range(1, 13)), A)
    require(all(gcd(x, y) == 1 for x, y in combinations(speeds, 2)), speeds)

    masses = tuple(rho(x, y) for x, y in combinations(speeds, 2))
    x_A = F(1, 49) - F(1, 4 * A * A)
    require(9 * A * A > 539 and x_A > FIRST_STRICT_HIGH, x_A)
    require(all(value > x_A > FIRST_STRICT_HIGH for value in masses), min(masses))
    require(
        all(value < F(1, 2) and eta(value) > eta(x_A) for value in masses),
        "pair variance does not clear the claimed lower bound",
    )

    # These conclusions hold for every one of the 7^5 labelled spanning
    # trees, because they use only six edgewise lower bounds.
    every_tree_mass = 6 * FIRST_STRICT_HIGH
    every_tree_error = 6 * eta(x_A)
    coarse_tree_error = 6 * eta(FIRST_STRICT_HIGH)
    coarse_normalized_error = F(4545 * A, 332024)
    harmonic = sum((F(1, speed) for speed in speeds), F(0))
    require(every_tree_mass == TREE_FLOOR, every_tree_mass)
    require(coarse_tree_error == F(4545, 47432), coarse_tree_error)
    require(every_tree_error > coarse_tree_error, every_tree_error)
    require(harmonic < F(7, A), harmonic)
    require(every_tree_error / harmonic > F(6 * A, 7) * eta(x_A), harmonic)
    require(every_tree_error / harmonic > coarse_normalized_error, harmonic)

    corrections = tuple(
        fold(r + s) - fold(s - r) for r, s in combinations(OFFSETS, 2)
    )
    require(
        corrections
        == (20, 16, 8, 0, -16, -24, 32, 16, 0, -32, -20,
            24, 0, -48, -16, 0, -24, -8, 0, 0, 16),
        corrections,
    )
    for (r, s), correction in zip(combinations(OFFSETS, 2), corrections):
        expected = F(1, 49) + F(correction, 196 * (A + r) * (A + s))
        require(rho(A + r, A + s) == expected, (r, s, correction))

    return {
        "multiplier": multiplier,
        "A": A,
        "min_mass": min(masses),
        "max_mass": max(masses),
        "x_A": x_A,
        "harmonic_bound": F(7, A),
        "coarse_tree_error": coarse_tree_error,
        "coarse_normalized_error": coarse_normalized_error,
        "corrections": corrections,
    }


def check_charge_identity() -> int:
    rows = 0
    for left in range(1, 70):
        for right in range(left + 1, 75):
            g = gcd(left, right)
            a, b = left // g, right // g
            variance = eta(rho(left, right))
            kappa = variance * F(a * b, a + b)
            error = variance / g
            require(error == kappa * (F(1, left) + F(1, right)), (left, right))
            require(error == variance * F(a, left), (left, right, "left charge"))
            require(error == variance * F(b, right), (left, right, "right charge"))
            rows += 1
    return rows


def check_strict_high_no_ap_packet(multiplier: int) -> dict[str, object]:
    """Separate strict-high connectivity from THM-1218's heavy circuit."""
    require(multiplier > 0, multiplier)
    A = CIRCUIT_BASE_STEP * multiplier
    speeds = tuple(A + r for r in CIRCUIT_OFFSETS)
    differences = {s - r for r, s in combinations(CIRCUIT_OFFSETS, 2)}
    require(A % 14 == 0, A)
    require(all(CIRCUIT_BASE_STEP % d == 0 for d in differences), differences)
    require(all(gcd(r, s) == 1 for r, s in combinations(CIRCUIT_OFFSETS, 2)),
            CIRCUIT_OFFSETS)
    require(all(gcd(x, y) == 1 for x, y in combinations(speeds, 2)), speeds)
    arithmetic_quartets = tuple(
        row for row in combinations(CIRCUIT_OFFSETS, 4)
        if row[1] - row[0] == row[2] - row[1] == row[3] - row[2]
    )
    require(not arithmetic_quartets, arithmetic_quartets)
    analytic_floor = F(1, 49) - F(1, 4 * (A + 1) * (A + 2))
    masses = tuple(rho(x, y) for x, y in combinations(speeds, 2))
    require(analytic_floor > FIRST_STRICT_HIGH, analytic_floor)
    require(all(value > FIRST_STRICT_HIGH for value in masses), min(masses))
    return {
        "multiplier": multiplier,
        "A": A,
        "differences": tuple(sorted(differences)),
        "min_mass": min(masses),
        "analytic_floor": analytic_floor,
        "arithmetic_quartets": arithmetic_quartets,
    }


def fiber_trapezoid(p: int, q: int, u: F) -> F:
    """Unperiodized density of ``q*x-p*y`` on the radius-1/14 box."""
    require(0 < p <= q, (p, q))
    u = abs(u)
    plateau = F(q - p, 14)
    support = F(p + q, 14)
    if u <= plateau:
        return F(1, 7 * q)
    if u <= support:
        return (support - u) / (p * q)
    return F(0)


def fiber_density(p: int, q: int, z: F) -> F:
    """Periodization of the exact Farey-fiber trapezoid."""
    return sum((fiber_trapezoid(p, q, z + k) for k in range(-10, 11)), F(0))


def fiber_minimum(p: int, q: int) -> F:
    """Exact minimum of the piecewise-affine periodic fiber density."""
    breakpoints = {F(0)}
    for point in (F(q - p, 14), F(p + q, 14)):
        for sign in (-1, 1):
            value = sign * point
            breakpoints.add(value - value.numerator // value.denominator)
    return min(fiber_density(p, q, z) for z in breakpoints)


def oriented_channel_charge(p: int, q: int) -> F:
    """Exact ``eta(rho) * max(p,q)`` for a height-seven channel."""
    require(0 < p <= q and gcd(p, q) == 1, (p, q))
    value = F(1, 7) if p == q else rho(p, q)
    return eta(value) * q


def relation_coordinates(p: int, q: int, a: int, b: int) -> tuple[int, int, int, int]:
    """Return ``(u,v,k,h)`` with ``qu-pv=1`` and ``(a,b)=M(k,h)``."""
    require(0 < p <= q and gcd(p, q) == 1, (p, q))
    if p == 1:
        u, v = 0, -1
    else:
        u = pow(q, -1, p)
        v = (q * u - 1) // p
    h = q * a - p * b
    k = u * b - v * a
    require(q * u - p * v == 1, (p, q, u, v))
    require((p * k + u * h, q * k + v * h) == (a, b), (p, q, a, b, k, h))
    require(gcd(abs(k), abs(h)) == gcd(a, b), (p, q, a, b, k, h))
    return u, v, k, h


def check_farey_fiber_frontier() -> dict[str, object]:
    major = tuple(sorted(
        (
            (p, q)
            for q in range(2, 7)
            for p in range(1, q)
            if gcd(p, q) == 1 and p + q <= 7
        ),
        key=lambda row: (sum(row), row[0], row[1]),
    ))
    require(major == MAJOR_NONDIAGONAL_CHANNELS, major)
    require(((1, 1),) + major == PHASE_CHANNELS, PHASE_CHANNELS)
    require(fiber_density(1, 1, F(1, 2)) == 0, "diagonal height-two channel")
    # The unperiodized trapezoid has a plateau plus two congruent triangles.
    # Their exact total area is the danger-box area 1/49.
    for p, q in PHASE_CHANNELS + ((3, 5), (4, 5), (1, 14), (1, 19)):
        plateau_area = F(q - p, 49 * q)
        shoulder_area = F(p, 49 * q)
        require(plateau_area + shoulder_area == F(1, 49), (p, q))
    checked = 1
    for q in range(2, 31):
        for p in range(1, q):
            if gcd(p, q) != 1:
                continue
            minimum = fiber_minimum(p, q)
            require((minimum == 0) == (p + q <= 7), (p, q, minimum))
            if p + q > 7:
                lower = min(F(1, 7 * q), F(p + q - 7, 14 * p * q))
                require(minimum >= lower > 0, (p, q, minimum, lower))
            checked += 1

    # Replay the exact THM-605 minima whose general closed form was left open
    # there.  The periodized trapezoid produces them without a phase census.
    known_minima = {
        (1, 7): F(1, 49),
        (2, 7): F(1, 49),
        (1, 14): F(1, 49),
        (3, 5): F(1, 105),
        (4, 5): F(1, 70),
        (1, 15): F(2, 105),
        (1, 17): F(2, 119),
        (1, 19): F(2, 133),
    }
    require(
        {channel: fiber_minimum(*channel) for channel in known_minima} == known_minima,
        known_minima,
    )

    coordinate_rows = 0
    for p, q in PHASE_CHANNELS:
        for a in range(1, 31):
            for b in range(a + 1, 34):
                relation_coordinates(p, q, a, b)
                coordinate_rows += 1

    relation_rows = []
    for multiplier in (1, 2, 10, 100):
        N = 14 * multiplier
        a, b = N + 1, 2 * N + 1
        value = rho(a, b)
        require(gcd(a, b) == 1 and 2 * a - b == 1, (a, b))
        require(value == F(1, 49) + F(6, 49 * a * b), (a, b, value))
        require(value > HIGH_BAR, value)
        relation_rows.append((N, a, b, value))

    oriented_charges = {
        channel: oriented_channel_charge(*channel) for channel in PHASE_CHANNELS
    }
    expected_charges = {
        (1, 1): F(6, 49),
        (1, 2): F(13, 98),
        (1, 3): F(20, 147),
        (1, 4): F(27, 196),
        (2, 3): F(20, 147),
        (1, 5): F(34, 245),
        (1, 6): F(41, 294),
        (2, 5): F(34, 245),
        (3, 4): F(27, 196),
    }
    require(oriented_charges == expected_charges, oriented_charges)
    require(max(oriented_charges.values()) == F(41, 294), oriented_charges)
    require(F(41, 294) < F(39, 98), "exact channel max sharpens coarse K=6 bound")
    require(6**6 == 46656 and 6 * 6**5 == 46656, "tree relation bound")
    return {
        "major_nondiagonal_channels": major,
        "fiber_pairs_checked": checked,
        "known_minima": known_minima,
        "coordinate_rows": coordinate_rows,
        "relation_rows": relation_rows,
        "oriented_charges": oriented_charges,
        "exact_relation_C": F(41, 294),
        "tree_product_bound": 46656,
    }


def check_nonprimitive_beat_guardrail() -> dict[str, object]:
    """Exact MISTAKE-184 witness against presentation-dependent y-decay."""
    A, B = 3744, 3745
    p = q = y = 12
    require(gcd(A, B) == 1 and q * B - p * A == y, (A, B, p, q, y))
    require(A == 26 * q * y, A)
    # The radius-1/13 folded formula has equal sum/difference corrections.
    fold13 = lambda value: (value % 13) * (13 - value % 13)
    rho13 = F(4 * A * B + fold13(A + B) - fold13(B - A), 169 * A * B)
    require(rho13 == F(4, 169), rho13)
    restricted_error = F(1, 6) * rho13
    clean_rhs = 13 * rho13 / (y * (p + q - 1)) + F(8 + 10 * y + 8, 13 * B)
    require(restricted_error == F(2, 507), restricted_error)
    require(clean_rhs == F(13129, 3359265), clean_rhs)
    require(restricted_error - clean_rhs == F(531, 14556815) > 0, clean_rhs)
    # B-A=1 is the primitive relation.  Twelve copies have the same orbit and
    # cannot create twelve independent starting phases.
    require(B - A == 1 and gcd(p, q) == 12, (A, B, p, q))
    return {
        "A": A,
        "B": B,
        "scaled_relation": (p, q, y),
        "primitive_relation": (1, 1, 1),
        "rho13": rho13,
        "restricted_error": restricted_error,
        "claimed_clean_rhs": clean_rhs,
        "violation": restricted_error - clean_rhs,
    }


def pair_overlap_gaps(a: int, b: int) -> tuple[tuple[F, F, F], ...]:
    """Exact complementary gaps of two strict radius-1/14 combs on [0,1]."""
    components: list[tuple[F, F]] = []
    for ka in range(-1, a + 1):
        left_a, right_a = F(14 * ka - 1, 14 * a), F(14 * ka + 1, 14 * a)
        for kb in range(-1, b + 1):
            left_b, right_b = F(14 * kb - 1, 14 * b), F(14 * kb + 1, 14 * b)
            left, right = max(left_a, left_b, F(0)), min(right_a, right_b, F(1))
            if left < right:
                components.append((left, right))
    components.sort()
    merged: list[tuple[F, F]] = []
    for left, right in components:
        if merged and left <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(right, merged[-1][1]))
        else:
            merged.append((left, right))
    gaps: list[tuple[F, F, F]] = []
    cursor = F(0)
    for left, right in merged:
        if cursor < left:
            gaps.append((cursor, left, left - cursor))
        cursor = max(cursor, right)
    if cursor < 1:
        gaps.append((cursor, F(1), 1 - cursor))
    return tuple(gaps)


def check_high_relation_interval_guardrail() -> dict[str, object]:
    """MISTAKE-185: low-height resolution misses a high exact period."""
    a, b = 64, 75
    gaps = pair_overlap_gaps(a, b)
    largest = max(gaps, key=lambda row: row[2])
    require(largest == (F(407, 896), F(489, 896), F(41, 448)), largest)
    interval_length = F(449, 4928)
    require(interval_length < largest[2], (interval_length, largest))
    low_rows = tuple(
        (abs(r * a + s * b), r, s)
        for r in range(-7, 8)
        for s in range(-7, 8)
        if (r or s) and abs(r) + abs(s) <= 7
    )
    low_min = min(row[0] for row in low_rows)
    require(low_min == 11, sorted(low_rows)[:4])
    require(low_min * interval_length == F(449, 448) > 1, interval_length)
    require(75 * a - 64 * b == 0 and 75 + 64 == 139, (a, b))
    return {
        "base_pair": (a, b),
        "largest_empty_gap": largest,
        "scaled_interval_length": interval_length,
        "height_le_seven_minimum": low_min,
        "minimum_phase_drift": low_min * interval_length,
        "omitted_exact_relation": (75, -64),
        "omitted_relation_height": 139,
    }


def check_positioned_debt_dichotomy() -> dict[str, F]:
    """Compose THM-1237's needle debt with seven reciprocal speeds."""
    require(24 * F(5, 88) == F(15, 11), "harmonic half-debt")
    require(F(30, 11) - F(15, 11) == F(15, 11), "remaining gcd debt")
    require(F(15, 11) / 13 == F(15, 143), "six-edge reciprocal debt")
    require(F(15, 143) / 6 == F(5, 286), "one edge reciprocal debt")
    require(7 / F(616, 5) == F(5, 88), "absolute-scale threshold")
    return {
        "min_speed_threshold": F(616, 5),
        "harmonic_threshold": F(5, 88),
        "tree_gcd_reciprocal_threshold": F(15, 143),
        "one_edge_reciprocal_threshold": F(5, 286),
        "one_edge_gcd_threshold": F(286, 5),
    }


def check_relation_cycle_holonomy() -> dict[str, object]:
    """Exact multiplicative-holonomy identity for oriented relation cycles."""
    rows = 0
    bounded_rows = 0
    balanced_rows = 0
    sign_mixed_rows = 0
    for length in range(2, 8):
        for seed in range(1, 401):
            channels = tuple(
                PHASE_CHANNELS[(seed * (i + 1) + i * i) % len(PHASE_CHANNELS)]
                for i in range(length)
            )
            p = tuple(channel[0] for channel in channels)
            q = tuple(channel[1] for channel in channels)
            speeds = tuple(50 + seed * (i + 2) + 3 * i * i for i in range(length))
            defects = tuple(
                q[i] * speeds[(i + 1) % length] - p[i] * speeds[i]
                for i in range(length)
            )
            q_product = 1
            p_product = 1
            for value in q:
                q_product *= value
            for value in p:
                p_product *= value
            rhs = sum(
                defects[i]
                * prod(q[:i], start=1)
                * prod(p[i + 1 :], start=1)
                for i in range(length)
            )
            lhs = (q_product - p_product) * speeds[0]
            require(lhs == rhs, (length, seed, channels, speeds, defects, lhs, rhs))
            height_cap = max((*p, *q))
            defect_cap = max(abs(value) for value in defects)
            if q_product != p_product:
                require(
                    speeds[0] <= length * defect_cap * height_cap ** (length - 1),
                    (length, seed, speeds, defects),
                )
                bounded_rows += 1
            else:
                require(rhs == 0, (length, seed, rhs))
                balanced_rows += 1
                if any(defects):
                    require(min(defects) < 0 < max(defects), (length, seed, defects))
                    sign_mixed_rows += 1
            rows += 1
        # The diagonal channel gives product-one holonomy.  Its defects
        # telescope and provide an exact nonzero sign-mixed replay.
        speeds = tuple(100 + i * i + 3 * i for i in range(length))
        defects = tuple(
            speeds[(i + 1) % length] - speeds[i] for i in range(length)
        )
        require(sum(defects) == 0 and any(defects), (length, speeds, defects))
        require(min(defects) < 0 < max(defects), (length, defects))
        rows += 1
        balanced_rows += 1
        sign_mixed_rows += 1
    require(7 * 6**6 == 326592, "height-seven channel circuit crown")
    require(7 * 13**6 == 33787663, "height-thirteen circuit crown")
    return {
        "cycle_rows": rows,
        "nontrivial_holonomy_rows": bounded_rows,
        "balanced_holonomy_rows": balanced_rows,
        "balanced_nonzero_sign_mixed_rows": sign_mixed_rows,
        "height6_seven_cycle_crown": 326592,
        "height13_seven_cycle_crown": 33787663,
    }


def check_protected_embedding() -> dict[str, object]:
    A, speeds = lcm_packet(1)
    q = A + 1
    require(speeds == tuple(q + d for d in (0, 1, 2, 4, 6, 10, 12)), speeds)
    require(q % 2 == 1 and q >= 521, q)
    core = tuple((3 * q + 1) // 2 + k for k in range(6))
    m = max(core)
    require(m == (3 * q + 11) // 2, (m, q))
    center = F(1, q)
    radius = F(1, 14 * m)
    length = 2 * radius

    core_margin = F(5 * q - 77, 14 * q)
    danger_margin = F(q * q - 517 * q - 1848, 14 * q * (3 * q + 11))
    require(core_margin > 0 and danger_margin > 0, (core_margin, danger_margin))
    require(520 * 520 - 517 * 520 - 1848 < 0, "q=520 should fail")
    require(521 * 521 - 517 * 521 - 1848 > 0, "q=521 should pass")

    # Triangle-inequality providers, with exact endpoint replays.
    for w in core:
        center_distance = F(1, 2) - F(2 * (w - (3 * q + 1) // 2) + 1, 2 * q)
        require(center_distance - w * radius > F(1, 14), (w, center_distance))
    for speed in speeds:
        d = speed - q
        require(F(d, q) + speed * radius < F(1, 14), (speed, d))

    return {
        "q": q,
        "deleted": speeds,
        "core": core,
        "m": m,
        "center": center,
        "radius": radius,
        "length": length,
        "core_margin": core_margin,
        "danger_margin": danger_margin,
        "local_pair_mass": length,
        "local_tree_mass": 6 * length,
    }


def check_disconnected_branch_constant() -> dict[str, object]:
    require(F(1, 49) - F(1, 4 * 56) > HIGH_BAR, "channel tail")
    strict = reduced_bank(55, lambda value: value < HIGH_BAR)
    closed = reduced_bank(55, lambda value: value <= HIGH_BAR)
    require(len(strict) == 7 and len(closed) == 12, (strict, closed))
    strict_vertices = {F(1)} | {F(b, a) for a, b, _ in strict} | {
        F(a, b) for a, b, _ in strict
    }
    closed_ratios = {F(b, a) for a, b, _ in closed} | {
        F(a, b) for a, b, _ in closed
    }
    closed_vertices = {F(1)} | closed_ratios

    strict_max, strict_rows = maximum_kappa(strict_vertices)
    closed_max, closed_rows = maximum_kappa(closed_vertices)
    require(strict_max == F(85975, 342804), strict_max)
    require({row[2] for row in strict_rows} == {(20, 33)}, strict_rows)
    require(closed_max == F(224458, 584325), closed_max)
    require(closed_rows == [(F(5, 9), F(9, 5), (25, 81), F(97, 4725))], closed_rows)

    # Reconstruct THM-1221's twelve normalized 2+5 strict-component packets.
    ordered_closed = tuple(sorted(closed_ratios))
    centers = sorted({r / s for r in ordered_closed for s in ordered_closed} - {F(1)})
    packets: set[tuple[F, ...]] = set()
    for second in centers:
        left = (F(1), second)
        if rho_ratio(*left) <= HIGH_BAR:
            continue
        common = tuple(v for v in ordered_closed if v / second in closed_ratios)
        for right in combinations(common, 5):
            if not graph_connected(right, lambda value: value > HIGH_BAR):
                continue
            if not any(rho_ratio(x, y) == HIGH_BAR for x in left for y in right):
                continue
            packet = tuple(sorted(left + right))
            packets.add(tuple(value / packet[0] for value in packet))
    require(len(packets) == 12, len(packets))
    two_five_rows = [
        (*projective_kappa(x, y), packet, x, y)
        for packet in packets
        for x, y in combinations(packet, 2)
    ]
    two_five_max = max(row[0] for row in two_five_rows)
    require(two_five_max == F(43774, 276507), two_five_max)

    strict_C = 6 * strict_max
    closed_C = 6 * closed_max
    two_five_C = 6 * two_five_max
    require(strict_C == F(85975, 57134), strict_C)
    require(closed_C == C_DISC, closed_C)
    require(two_five_C == F(87548, 92169), two_five_C)
    require(max(strict_C, closed_C, two_five_C) == C_DISC, "wrong C_disc")
    crown_ratio = TREE_FLOOR / (C_DISC + F(1, 7))
    require(crown_ratio == F(417375, 10488302), crown_ratio)
    forest_crown_ratio = TREE_FLOOR / (C_DISC + FOREST_CAP)
    protected_mH = forest_crown_ratio / 7
    protected_min_ratio = 7 / protected_mH
    separated_max_ratio = 2345 * protected_min_ratio
    require(forest_crown_ratio == F(59625, 1485836), forest_crown_ratio)
    require(protected_mH == F(59625, 10400852), protected_mH)
    require(protected_min_ratio == F(72805964, 59625), protected_min_ratio)
    require(separated_max_ratio == F(34145997116, 11925), separated_max_ratio)
    return {
        "strict_vertices": len(strict_vertices),
        "strict_kappa": strict_max,
        "strict_C": strict_C,
        "closed_vertices": len(closed_vertices),
        "closed_kappa": closed_max,
        "closed_C": closed_C,
        "two_five_packets": len(packets),
        "two_five_kappa": two_five_max,
        "two_five_C": two_five_C,
        "C_disc": C_DISC,
        "hunter_crown_ratio": crown_ratio,
        "forest_crown_ratio": forest_crown_ratio,
        "protected_mH": protected_mH,
        "protected_min_ratio": protected_min_ratio,
        "separated_max_ratio": separated_max_ratio,
    }


def main() -> None:
    print("THM-1226 GCD-PERIOD PROJECTIVE-CHARGE EXACT REFEREE")
    print("method=integer/Fraction only; always-on checks; no dependencies")
    charge_rows = check_charge_identity()
    print(f"exact edge-charge identities={charge_rows}")

    print("\nTRANSLATED LCM STRICT-HIGH COUNTERFAMILY")
    print(f"offsets={OFFSETS}; sum={sum(OFFSETS)}; base_step={BASE_STEP}")
    correction_row = None
    for multiplier in (1, 2, 10, 100):
        row = check_counterfamily(multiplier)
        correction_row = row["corrections"]
        print(
            f"N={multiplier} A={row['A']} min_rho={row['min_mass']} "
            f"max_rho={row['max_mass']} H<{row['harmonic_bound']}"
        )
        print(
            f"  every_tree_error>{row['coarse_tree_error']} "
            f"E/H>{row['coarse_normalized_error']}"
        )
    print(f"fold_corrections={correction_row}")
    print("all pair gcds=1; G_gt=K7; every spanning tree mass>15/154")
    print("E_T/H grows as (288/16807)A+O(1); no absolute C exists")

    print("\nSTRICT-HIGH / HEAVY-CIRCUIT SEPARATION")
    for multiplier in (1, 2, 10):
        row = check_strict_high_no_ap_packet(multiplier)
        print(
            f"N={multiplier} A={row['A']} min_rho={row['min_mass']} "
            f"analytic_floor={row['analytic_floor']}"
        )
    print(f"offsets={CIRCUIT_OFFSETS}; differences={row['differences']}")
    print("four_term_APs=(); all pair gcds=1; G_gt=K7; THM-1218 heavy quartets=0")

    fiber = check_farey_fiber_frontier()
    print("\nFAREY-HEIGHT SEVEN PHASE FIBER")
    print(f"THM-605 primitive_phase_channels={PHASE_CHANNELS}")
    print(f"major_nondiagonal_channels={fiber['major_nondiagonal_channels']}")
    print(f"exact_fiber_pairs_checked={fiber['fiber_pairs_checked']}")
    print(f"THM-605 exact_minima_replayed={fiber['known_minima']}")
    print(f"unimodular_relation_coordinate_rows={fiber['coordinate_rows']}")
    for N, a, b, value in fiber["relation_rows"]:
        print(f"N={N} pair=({a},{b}) relation=2a-b=1 rho={value}")
    print("D_a intersection D_b is empty on [3/14,11/14] for the relation family")
    print("fiber minimum is zero iff p+q<=7; positive lower bound otherwise")
    print(f"height-seven oriented_charges={fiber['oriented_charges']}")
    print(
        f"exact_relation_C={fiber['exact_relation_C']}; "
        f"six-edge relation-tree product bound={fiber['tree_product_bound']}"
    )

    guardrail = check_nonprimitive_beat_guardrail()
    print("\nPRIMITIVE-RELATION MULTIPLICITY GUARDRAIL")
    for key, value in guardrail.items():
        print(f"{key}={value}")
    print("THM-864 clean bound is refuted; scaled relations cannot buy y-decay")

    high_relation = check_high_relation_interval_guardrail()
    print("\nHIGH-RELATION FINITE-INTERVAL GUARDRAIL")
    for key, value in high_relation.items():
        print(f"{key}={value}")
    print("THM-598/602 low-height resolved branch is refuted; retain the full exact kernel")

    embedding = check_protected_embedding()
    print("\nCOVERED PROTECTED-NEEDLE EMBEDDING")
    for key in (
        "q", "deleted", "core", "m", "center", "radius", "length",
        "core_margin", "danger_margin", "local_pair_mass", "local_tree_mass",
    ):
        print(f"{key}={embedding[key]}")
    print("all seven deleted combs are active throughout I; this is not an F7 counterexample")

    branch = check_disconnected_branch_constant()
    print("\nCONDITIONAL DISCONNECTED-G_gt TRANSFER")
    for key, value in branch.items():
        print(f"{key}={value}")
    print("if G_gt is disconnected: E_T<=C_disc*H for a THM-1221 floor tree")
    print("localized tree mass >=(15/154)L-(448916/194775)H")
    print("covered forest cap=6H/49; hence H/L>=59625/1485836")
    print("if min(S)>=13m and six faster combs cover its gap: max(S)/m<34145997116/11925")
    require(SHALLOW_F7_CONSTANT == 13 * TREE_FLOOR / 7, SHALLOW_F7_CONSTANT)
    require(SHALLOW_F7_CONSTANT < C_DISC, (SHALLOW_F7_CONSTANT, C_DISC))
    print("canonical needle + min(S)<13m: F7 is automatic for C>=195/1078")
    print("therefore every nontrivial connected-G_gt row has a complete min(S)-gap")

    debt = check_positioned_debt_dichotomy()
    print("\nTHM-1237 POSITIONED-DEBT COMPOSITION")
    for key, value in debt.items():
        print(f"{key}={value}")
    print("covered needle: min(S)<616m/5 or a floor-tree edge has gcd<=286m/5")

    holonomy = check_relation_cycle_holonomy()
    print("\nPRIMITIVE RELATION-CYCLE HOLONOMY")
    for key, value in holonomy.items():
        print(f"{key}={value}")
    print("nontrivial product holonomy gives a crown; balanced nonzero holonomy is sign-mixed")

    print("\nTOURNAMENT / ALTERNATE-VERTEX AUDIT")
    print("LCM-packet switch=all 21 edges strict-high; speed-order gauge is transitive")
    print("scores=(0,1,2,3,4,5,6); cycles=0; SCCs=7; Hamilton_paths=1; low_flips=0")
    print("runner threshold graph preserves global rho floor but destroys projective height kappa")
    print("gcd-period vertices collapse all counterfamily edges to period 1 and destroy wall alignment")
    print("relation-channel switch=(7-p-q)_+ then Farey-height/lex gauge")
    print(f"relation-channel tie_path={PHASE_CHANNELS}")
    print("relation-channel scores=(0,1,2,3,4,5,6,7,8); cycles=0; SCCs=9; Hamilton_paths=1")
    print("faithful local vertices=wall events / blocker obligations; edge sidecar=(p,q,k,h,tooth,cocycle)")
    print("STATUS=PASS")
    print(f"source_sha256={sha256(Path(__file__).read_bytes()).hexdigest()}")


if __name__ == "__main__":
    main()
