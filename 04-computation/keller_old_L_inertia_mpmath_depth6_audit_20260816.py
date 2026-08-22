#!/usr/bin/env python3
"""High-precision ancestry-local continuation audit at old L, depth six.

This is a bounded numerical sidecar, not an exact inertia theorem.  It tracks
all 3^6 inverse leaves of the fixed sporadic Keller map around the canonical
transverse loop

    (a,b,c) = (2/27 + rho*exp(2*pi*i*s), 1, 1).

Unlike the older double-precision scout, continuation is performed locally in
each labelled ternary ancestry block.  Every cubic root is Newton-refined at
high precision, supplied with an a-posteriori Rouche disk, checked against the
other two roots, and reconstructed as a full inverse point whose forward image
is tested.  Endpoint matching is also ancestry-local and must project through
all six levels.  Repeated precision, radius, and step controls are required.

The script deliberately uses explicit ``require`` failures and contains no
executable ``assert`` statements, so ``python -O`` preserves every gate.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from itertools import permutations
import json
from math import lcm

from mpmath import mp


Point = tuple[mp.mpc, mp.mpc, mp.mpc]
Permutation = tuple[int, ...]

DEPTH = 6
EXPECTED_DEPTH6_HISTOGRAM = ((1, 243), (12, 9), (24, 9), (54, 1), (108, 1))
EXPECTED_DEPTH6_EXPONENT = 466
EXPECTED_REFLECTION_GAUGE_A_HISTOGRAM = ((1, 1), (2, 121))
EXPECTED_REFLECTION_GAUGE_B_HISTOGRAM = ((1, 7), (2, 118))
EXPECTED_C_HISTOGRAM = ((6, 9), (12, 9), (27, 1), (54, 1))
EXPECTED_C_ORBIT_COUNT = 20
EXPECTED_C_ORDER = 108
EXPECTED_SEMANTIC_SHA256 = (
    "dc22c0fd73fa57c92fdb62d41acce8edc9127fcf36e93d24bc5bf5e49006ceb4"
)


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def l_value(point: Point) -> mp.mpc:
    a, b, c = point
    return 27 * a * a * c * c - 18 * a * b * c + 16 * a + b**3 * c - b * b


def forward(point: Point) -> Point:
    x, y, z = point
    unit = 1 + x * y
    return (
        unit**3 * z + y * y * unit * (4 + 3 * x * y),
        y + 3 * x * unit * unit * z + 3 * x * y * y * (4 + 3 * x * y),
        2 * x - 3 * x * x * y - x**3 * z,
    )


def relative_error(left: Point, right: Point) -> mp.mpf:
    return max(
        abs(x - y) / max(mp.mpf(1), abs(x), abs(y))
        for x, y in zip(left, right)
    )


def chordal(left: mp.mpc, right: mp.mpc) -> mp.mpf:
    return abs(left - right) / (
        mp.sqrt(1 + abs(left) ** 2) * mp.sqrt(1 + abs(right) ** 2)
    )


def point_distance(left: Point, right: Point) -> mp.mpf:
    return mp.sqrt(sum(chordal(x, y) ** 2 for x, y in zip(left, right)))


def core_coefficients(target: Point) -> tuple[mp.mpc, mp.mpc, mp.mpc]:
    _, b, c = target
    return l_value(target), 4 - 3 * b * c, c


@dataclass(frozen=True)
class RootCertificate:
    point: Point
    root_residual: mp.mpf
    forward_residual: mp.mpf
    rouche_ratio: mp.mpf
    projective_radius: mp.mpf
    denominator_ratio: mp.mpf
    newton_iterations: int
    coordinate_size: mp.mpf


def rouche_disk(
    z: mp.mpc,
    f: mp.mpc,
    fp: mp.mpc,
    fpp: mp.mpc,
    cubic_coefficient: mp.mpc,
) -> tuple[mp.mpf, mp.mpf]:
    """Certify one simple root in a disk by a linear-term Rouche bound."""

    require(abs(fp) > 0, ("zero derivative", z, f))
    beta = abs(f) / abs(fp)
    radius = 4 * beta + 32 * mp.eps * (1 + abs(z))
    linear_boundary = abs(fp) * radius
    remainder = (
        abs(f)
        + abs(fpp) * radius**2 / 2
        + abs(cubic_coefficient) * radius**3
    )
    require(linear_boundary > remainder, (
        "Rouche disk failed", z, radius, linear_boundary, remainder
    ))
    return radius, remainder / linear_boundary


def refine_x_root(
    target: Point, seed_x: mp.mpc
) -> tuple[mp.mpc, int, mp.mpf, mp.mpf, mp.mpf]:
    """Continue one labelled cubic root, using x or u=1/x projectively."""

    lead, linear, c = core_coefficients(target)
    tolerance = mp.power(10, -(mp.dps - 28))
    reciprocal = abs(seed_x) > 2
    if reciprocal:
        value = 1 / seed_x
        for iteration in range(1, 21):
            f = lead + linear * value**2 - 2 * c * value**3
            fp = 2 * linear * value - 6 * c * value**2
            require(abs(fp) > 0, ("reciprocal Newton derivative", target, value))
            correction = f / fp
            value -= correction
            if abs(correction) <= tolerance * max(mp.mpf(1), abs(value)):
                break
        else:
            raise RuntimeError(("reciprocal Newton did not converge", target, seed_x))
        f = lead + linear * value**2 - 2 * c * value**3
        fp = 2 * linear * value - 6 * c * value**2
        fpp = 2 * linear - 12 * c * value
        radius, ratio = rouche_disk(value, f, fp, fpp, -2 * c)
        require(abs(value) > radius, ("reciprocal disk contains zero", value, radius))
        root = 1 / value
        projective_radius = radius
    else:
        value = seed_x
        for iteration in range(1, 21):
            f = lead * value**3 + linear * value - 2 * c
            fp = 3 * lead * value**2 + linear
            require(abs(fp) > 0, ("affine Newton derivative", target, value))
            correction = f / fp
            value -= correction
            if abs(correction) <= tolerance * max(mp.mpf(1), abs(value)):
                break
        else:
            raise RuntimeError(("affine Newton did not converge", target, seed_x))
        f = lead * value**3 + linear * value - 2 * c
        fp = 3 * lead * value**2 + linear
        fpp = 6 * lead * value
        radius, ratio = rouche_disk(value, f, fp, fpp, lead)
        root = value
        projective_radius = radius
    scale = max(
        mp.mpf(1), abs(lead * root**3), abs(linear * root), abs(2 * c)
    )
    residual = abs(lead * root**3 + linear * root - 2 * c) / scale
    return root, iteration, residual, ratio, projective_radius


def reconstruct_point(target: Point, root: mp.mpc) -> tuple[Point, mp.mpf]:
    a, b, c = target
    terms = (12 * a * root * root, -b * b * root * root, b * root, mp.mpc(2))
    denominator = sum(terms)
    denominator_scale = max(mp.mpf(1), *(abs(term) for term in terms))
    denominator_ratio = abs(denominator) / denominator_scale
    require(denominator_ratio > mp.power(10, -(mp.dps // 3)), (
        "inverse denominator cancellation", target, root, denominator_ratio
    ))
    y = b - 3 * a * root * (9 * a * c * root - b * root + 2) / denominator
    require(abs(root) > mp.power(10, -(mp.dps // 3)), ("inverse x near zero", target, root))
    z = (2 * root - 3 * root * root * y - c) / root**3
    return (root, y, z), denominator_ratio


def certify_from_seed(target: Point, seed_x: mp.mpc) -> RootCertificate:
    root, iterations, root_residual, rouche_ratio, projective_radius = (
        refine_x_root(target, seed_x)
    )
    point, denominator_ratio = reconstruct_point(target, root)
    forward_residual = relative_error(forward(point), target)
    return RootCertificate(
        point=point,
        root_residual=root_residual,
        forward_residual=forward_residual,
        rouche_ratio=rouche_ratio,
        projective_radius=projective_radius,
        denominator_ratio=denominator_ratio,
        newton_iterations=iterations,
        coordinate_size=max(abs(value) for value in point),
    )


def scaled_initial_roots(target: Point) -> tuple[mp.mpc, mp.mpc, mp.mpc]:
    lead, linear, c = core_coefficients(target)
    require(abs(lead) > 0, ("zero leading coefficient on loop", target))
    candidates = [mp.mpf(1)]
    if abs(linear) > 0:
        candidates.append(mp.sqrt(abs(linear / lead)))
    if abs(c) > 0:
        candidates.append(abs(2 * c / lead) ** (mp.mpf(1) / 3))
    scale_x = max(candidates)
    coefficients = [lead * scale_x**3, mp.mpc(0), linear * scale_x, -2 * c]
    coefficient_scale = max(abs(value) for value in coefficients)
    normalized = [value / coefficient_scale for value in coefficients]
    roots = mp.polyroots(normalized, maxsteps=300, extraprec=80, error=False)
    require(len(roots) == 3, ("initial cubic root count", target, roots))
    refined = [certify_from_seed(target, scale_x * root) for root in roots]
    points = [item.point for item in refined]
    order = sorted(
        range(3),
        key=lambda index: (
            mp.arg(points[index][0]),
            abs(points[index][0]),
            mp.re(points[index][0]),
            mp.im(points[index][0]),
        ),
    )
    return tuple(points[index][0] for index in order)


@dataclass
class Metrics:
    max_root_residual: mp.mpf = mp.mpf(0)
    max_forward_residual: mp.mpf = mp.mpf(0)
    max_rouche_ratio: mp.mpf = mp.mpf(0)
    max_projective_radius: mp.mpf = mp.mpf(0)
    min_denominator_ratio: mp.mpf = mp.inf
    max_coordinate_size: mp.mpf = mp.mpf(0)
    max_newton_iterations: int = 0
    max_step_cost: mp.mpf = mp.mpf(0)
    min_step_margin: mp.mpf = mp.inf
    min_root_separation_ratio: mp.mpf = mp.inf
    endpoint_cost: mp.mpf = mp.mpf(0)
    endpoint_margin: mp.mpf = mp.inf

    def absorb(self, certificate: RootCertificate) -> None:
        self.max_root_residual = max(self.max_root_residual, certificate.root_residual)
        self.max_forward_residual = max(
            self.max_forward_residual, certificate.forward_residual
        )
        self.max_rouche_ratio = max(self.max_rouche_ratio, certificate.rouche_ratio)
        self.max_projective_radius = max(
            self.max_projective_radius, certificate.projective_radius
        )
        self.min_denominator_ratio = min(
            self.min_denominator_ratio, certificate.denominator_ratio
        )
        self.max_coordinate_size = max(
            self.max_coordinate_size, certificate.coordinate_size
        )
        self.max_newton_iterations = max(
            self.max_newton_iterations, certificate.newton_iterations
        )


def assignment_costs(left: list[Point], right: list[Point]) -> list[tuple[mp.mpf, tuple[int, ...]]]:
    require(len(left) == len(right) == 3, ("assignment block size", len(left), len(right)))
    costs = [
        [point_distance(left[row], right[column]) for column in range(3)]
        for row in range(3)
    ]
    rows = []
    for permutation in permutations(range(3)):
        cost = sum(costs[row][permutation[row]] for row in range(3))
        rows.append((cost, tuple(permutation)))
    return sorted(rows, key=lambda row: row[0])


def update_step_matching(previous: list[Point], current: list[Point], metrics: Metrics) -> None:
    rows = assignment_costs(previous, current)
    require(rows[0][1] == (0, 1, 2), ("continuation branch jump", rows[:2]))
    best = rows[0][0]
    second = rows[1][0]
    margin = second / max(best, mp.eps)
    require(margin > 2, ("continuation assignment margin", best, second, margin))
    metrics.max_step_cost = max(
        metrics.max_step_cost,
        max(point_distance(previous[index], current[index]) for index in range(3)),
    )
    metrics.min_step_margin = min(metrics.min_step_margin, margin)


def certify_distinct_roots(certificates: list[RootCertificate], metrics: Metrics) -> None:
    require(len(certificates) == 3, ("root certificate count", len(certificates)))
    for left in range(3):
        for right in range(left + 1, 3):
            separation = chordal(
                certificates[left].point[0], certificates[right].point[0]
            )
            radii = (
                certificates[left].projective_radius
                + certificates[right].projective_radius
            )
            ratio = separation / max(radii, mp.eps)
            require(ratio > 100, ("root disks not separated", left, right, ratio))
            metrics.min_root_separation_ratio = min(
                metrics.min_root_separation_ratio, ratio
            )


def initial_block(target: Point, metrics: Metrics) -> list[Point]:
    roots = scaled_initial_roots(target)
    certificates = [certify_from_seed(target, root) for root in roots]
    certify_distinct_roots(certificates, metrics)
    for certificate in certificates:
        metrics.absorb(certificate)
    return [certificate.point for certificate in certificates]


def continued_block(target: Point, previous: list[Point], metrics: Metrics) -> list[Point]:
    certificates = [certify_from_seed(target, point[0]) for point in previous]
    certify_distinct_roots(certificates, metrics)
    current = [certificate.point for certificate in certificates]
    update_step_matching(previous, current, metrics)
    for certificate in certificates:
        metrics.absorb(certificate)
    return current


def initial_levels(target: Point, depth: int, metrics: Metrics) -> list[list[Point]]:
    levels: list[list[Point]] = []
    parents = [target]
    for _ in range(depth):
        children: list[Point] = []
        for parent in parents:
            children.extend(initial_block(parent, metrics))
        levels.append(children)
        parents = children
    return levels


def continue_levels(
    target: Point, previous: list[list[Point]], depth: int, metrics: Metrics
) -> list[list[Point]]:
    levels: list[list[Point]] = []
    first = continued_block(target, previous[0], metrics)
    levels.append(first)
    for level_index in range(1, depth):
        children: list[Point] = []
        for parent_index, parent in enumerate(levels[level_index - 1]):
            old_block = previous[level_index][
                3 * parent_index : 3 * parent_index + 3
            ]
            children.extend(continued_block(parent, old_block, metrics))
        levels.append(children)
    return levels


def endpoint_block(
    final: list[Point], initial: list[Point], metrics: Metrics
) -> tuple[int, int, int]:
    rows = assignment_costs(final, initial)
    best, permutation = rows[0]
    second = rows[1][0]
    margin = second / max(best, mp.eps)
    require(margin > 100, ("endpoint assignment margin", best, second, margin))
    metrics.endpoint_cost = max(metrics.endpoint_cost, best)
    metrics.endpoint_margin = min(metrics.endpoint_margin, margin)
    return permutation


def endpoint_tree_permutations(
    initial: list[list[Point]], final: list[list[Point]], metrics: Metrics
) -> tuple[Permutation, ...]:
    root = endpoint_block(final[0], initial[0], metrics)
    answer: list[Permutation] = [root]
    for level_index in range(1, len(initial)):
        parent = answer[-1]
        child = [-1] * len(initial[level_index])
        for source_parent, target_parent in enumerate(parent):
            final_block = final[level_index][
                3 * source_parent : 3 * source_parent + 3
            ]
            initial_block = initial[level_index][
                3 * target_parent : 3 * target_parent + 3
            ]
            local = endpoint_block(final_block, initial_block, metrics)
            for source_child, target_child in enumerate(local):
                child[3 * source_parent + source_child] = (
                    3 * target_parent + target_child
                )
        require(-1 not in child, ("incomplete endpoint level", level_index + 1))
        require(sorted(child) == list(range(len(child))), (
            "endpoint is not a permutation", level_index + 1
        ))
        for source_child, target_child in enumerate(child):
            require(target_child // 3 == parent[source_child // 3], (
                "endpoint block projection", level_index + 1, source_child, target_child
            ))
        answer.append(tuple(child))
    return tuple(answer)


def permutation_cycles(permutation: Permutation) -> tuple[tuple[int, ...], ...]:
    seen = [False] * len(permutation)
    cycles: list[tuple[int, ...]] = []
    for start in range(len(permutation)):
        if seen[start]:
            continue
        cycle: list[int] = []
        cursor = start
        while not seen[cursor]:
            seen[cursor] = True
            cycle.append(cursor)
            cursor = permutation[cursor]
        cycles.append(tuple(cycle))
    return tuple(cycles)


def cycle_histogram(permutation: Permutation) -> tuple[tuple[int, int], ...]:
    lengths = [len(cycle) for cycle in permutation_cycles(permutation)]
    return tuple((length, lengths.count(length)) for length in sorted(set(lengths)))


def permutation_order(permutation: Permutation) -> int:
    value = 1
    for cycle in permutation_cycles(permutation):
        value = lcm(value, len(cycle))
    return value


def compose(left: Permutation, right: Permutation) -> Permutation:
    require(len(left) == len(right), ("permutation composition", len(left), len(right)))
    return tuple(left[right[index]] for index in range(len(right)))


def inverse_permutation(permutation: Permutation) -> Permutation:
    inverse = [-1] * len(permutation)
    for source, target in enumerate(permutation):
        inverse[target] = source
    require(-1 not in inverse, "incomplete inverse permutation")
    return tuple(inverse)


def reflection_factorization(rotation: Permutation) -> tuple[Permutation, Permutation]:
    """Return deterministic involutions B,A with rotation = B o A.

    Equal-length rotation cycles are paired whenever possible.  On the one
    unpaired cycle of each occurring length, the reflection axis is chosen so
    that an even cycle contributes zero A-fixed points and two B-fixed points;
    an odd cycle contributes one fixed point to each.  For the depth-six C
    profile this gives A=(2^121,1) and B=(2^118,1^7).
    """

    cycles_by_length: dict[int, list[tuple[int, ...]]] = {}
    for cycle in permutation_cycles(rotation):
        cycles_by_length.setdefault(len(cycle), []).append(cycle)
    reflection_a = [-1] * len(rotation)
    for length in sorted(cycles_by_length):
        cycles = cycles_by_length[length]
        paired_count = len(cycles) - (len(cycles) % 2)
        for pair_start in range(0, paired_count, 2):
            left = cycles[pair_start]
            right = cycles[pair_start + 1]
            for index in range(length):
                reflection_a[left[index]] = right[-index % length]
                reflection_a[right[index]] = left[-index % length]
        if paired_count != len(cycles):
            cycle = cycles[-1]
            offset = 1 if length % 2 == 0 else 0
            for index in range(length):
                reflection_a[cycle[index]] = cycle[(offset - index) % length]
    require(-1 not in reflection_a, "incomplete reflection factorization")
    a = tuple(reflection_a)
    b = compose(rotation, a)
    identity = tuple(range(len(rotation)))
    require(compose(a, a) == identity, "normalized A is not involutive")
    require(compose(b, b) == identity, "normalized B is not involutive")
    require(compose(b, a) == rotation, "normalized reflections do not multiply to C")
    return a, b


def root_action_and_sections(
    permutation: Permutation, depth: int
) -> tuple[Permutation, tuple[Permutation, Permutation, Permutation]]:
    block_size = 3 ** (depth - 1)
    require(len(permutation) == 3 * block_size, ("root section degree", depth))
    action: list[int] = []
    sections: list[Permutation] = []
    for source_root in range(3):
        images = permutation[
            source_root * block_size : (source_root + 1) * block_size
        ]
        target_roots = {image // block_size for image in images}
        require(len(target_roots) == 1, ("root block image", source_root, target_roots))
        target_root = next(iter(target_roots))
        section = tuple(image % block_size for image in images)
        require(sorted(section) == list(range(block_size)), (
            "root section permutation", source_root
        ))
        action.append(target_root)
        sections.append(section)
    require(sorted(action) == [0, 1, 2], ("root action", action))
    return tuple(action), (sections[0], sections[1], sections[2])


def target_point(radius: mp.mpf, step: int, steps: int) -> Point:
    if step == steps:
        parameter = radius
    else:
        angle = 2 * mp.pi * step / steps
        parameter = radius * mp.mpc(mp.cos(angle), mp.sin(angle))
    return (mp.mpf(2) / 27 + parameter, mp.mpc(1), mp.mpc(1))


@dataclass(frozen=True)
class RunSummary:
    radius: str
    steps: int
    dps: int
    permutations: tuple[Permutation, ...]
    histograms: tuple[tuple[tuple[int, int], ...], ...]
    exponents: tuple[int, ...]
    root_action: Permutation
    a_histogram: tuple[tuple[int, int], ...]
    b_histogram: tuple[tuple[int, int], ...]
    normalized_a_histogram: tuple[tuple[int, int], ...]
    normalized_b_histogram: tuple[tuple[int, int], ...]
    c_histogram: tuple[tuple[int, int], ...]
    c_orbit_count: int
    c_order: int
    metrics: Metrics


def run_control(radius_text: str, steps: int, dps: int, depth: int = DEPTH) -> RunSummary:
    mp.dps = dps
    radius = mp.mpf(radius_text)
    metrics = Metrics()
    initial_target = target_point(radius, 0, steps)
    initial = initial_levels(initial_target, depth, metrics)
    tracked = [list(level) for level in initial]
    for step in range(1, steps + 1):
        tracked = continue_levels(
            target_point(radius, step, steps), tracked, depth, metrics
        )
    endpoint = endpoint_tree_permutations(initial, tracked, metrics)
    histograms = tuple(cycle_histogram(permutation) for permutation in endpoint)
    exponents = tuple(
        len(permutation) - len(permutation_cycles(permutation))
        for permutation in endpoint
    )
    root_action, sections = root_action_and_sections(endpoint[-1], depth)
    fixed_roots = [index for index, image in enumerate(root_action) if index == image]
    swapped_roots = [index for index, image in enumerate(root_action) if index != image]
    require(len(fixed_roots) == 1 and len(swapped_roots) == 2, (
        "old-L root action is not a transposition", root_action
    ))
    require(root_action[swapped_roots[0]] == swapped_roots[1], root_action)
    require(root_action[swapped_roots[1]] == swapped_roots[0], root_action)
    fixed_section = sections[fixed_roots[0]]
    source_a, source_b = swapped_roots
    section_a = sections[source_a]
    section_b = sections[source_b]
    identity = tuple(range(len(fixed_section)))
    require(fixed_section == identity, "depth-six fixed root section")
    # A section from one exchanged block lands in the other block.  Its raw
    # cycle type is therefore not invariant under independent ancestry gauges
    # on those two blocks.  The return map C = section_b o section_a is the
    # restriction of the square of the global monodromy and is conjugacy
    # invariant.  The return map based in the other block must have the same
    # cycle data and order.
    rotation = compose(section_b, section_a)
    reverse_rotation = compose(section_a, section_b)
    require(cycle_histogram(rotation) == cycle_histogram(reverse_rotation), (
        "two-block return histograms disagree",
        cycle_histogram(rotation), cycle_histogram(reverse_rotation)
    ))
    require(permutation_order(rotation) == permutation_order(reverse_rotation), (
        "two-block return orders disagree",
        permutation_order(rotation), permutation_order(reverse_rotation)
    ))
    normalized_a, normalized_b = reflection_factorization(rotation)
    # Explicitly realize this factorization as a change of coordinates on the
    # second exchanged ancestry block.  The first block keeps its coordinates.
    second_block_transport = compose(section_a, normalized_a)
    transported_a = compose(
        inverse_permutation(second_block_transport), section_a
    )
    transported_b = compose(section_b, second_block_transport)
    require(transported_a == normalized_a, "A gauge transport failed")
    require(transported_b == normalized_b, "B gauge transport failed")
    summary = RunSummary(
        radius=radius_text,
        steps=steps,
        dps=dps,
        permutations=endpoint,
        histograms=histograms,
        exponents=exponents,
        root_action=root_action,
        a_histogram=cycle_histogram(section_a),
        b_histogram=cycle_histogram(section_b),
        normalized_a_histogram=cycle_histogram(normalized_a),
        normalized_b_histogram=cycle_histogram(normalized_b),
        c_histogram=cycle_histogram(rotation),
        c_orbit_count=len(permutation_cycles(rotation)),
        c_order=permutation_order(rotation),
        metrics=metrics,
    )
    residual_gate = mp.power(10, -(dps // 3))
    require(metrics.max_forward_residual < residual_gate, (
        "forward residual gate", radius_text, steps, dps,
        metrics.max_forward_residual, residual_gate
    ))
    require(metrics.max_root_residual < residual_gate, (
        "root residual gate", radius_text, metrics.max_root_residual, residual_gate
    ))
    require(metrics.max_rouche_ratio < mp.mpf("0.5"), (
        "Rouche ratio gate", radius_text, metrics.max_rouche_ratio
    ))
    if depth == DEPTH:
        require(histograms[-1] == EXPECTED_DEPTH6_HISTOGRAM, (
            "depth-six cycle hostile", radius_text, histograms[-1]
        ))
        require(exponents[-1] == EXPECTED_DEPTH6_EXPONENT, (
            "depth-six exponent hostile", radius_text, exponents[-1]
        ))
        require(summary.c_histogram == EXPECTED_C_HISTOGRAM, (
            "C section hostile", radius_text, summary.c_histogram
        ))
        require(
            summary.normalized_a_histogram
            == EXPECTED_REFLECTION_GAUGE_A_HISTOGRAM,
            ("normalized A section hostile", radius_text,
             summary.normalized_a_histogram),
        )
        require(
            summary.normalized_b_histogram
            == EXPECTED_REFLECTION_GAUGE_B_HISTOGRAM,
            ("normalized B section hostile", radius_text,
             summary.normalized_b_histogram),
        )
        require(summary.c_orbit_count == EXPECTED_C_ORBIT_COUNT, (
            "C orbit count hostile", radius_text, summary.c_orbit_count
        ))
        require(summary.c_order == EXPECTED_C_ORDER, (
            "C order hostile", radius_text, summary.c_order
        ))
    return summary


def sci(value: mp.mpf, digits: int = 8) -> str:
    return mp.nstr(value, digits, min_fixed=0, max_fixed=0)


def run_record(run: RunSummary) -> dict[str, object]:
    metrics = run.metrics
    return {
        "radius": run.radius,
        "steps": run.steps,
        "dps": run.dps,
        "histograms": run.histograms,
        "exponents": run.exponents,
        "root_action": run.root_action,
        "A": run.a_histogram,
        "B": run.b_histogram,
        "reflection_gauge_A": run.normalized_a_histogram,
        "reflection_gauge_B": run.normalized_b_histogram,
        "C": run.c_histogram,
        "C_orbit_count": run.c_orbit_count,
        "C_order": run.c_order,
        "max_root_residual": sci(metrics.max_root_residual),
        "max_forward_residual": sci(metrics.max_forward_residual),
        "max_rouche_ratio": sci(metrics.max_rouche_ratio),
        "max_projective_radius": sci(metrics.max_projective_radius),
        "min_denominator_ratio": sci(metrics.min_denominator_ratio),
        "max_coordinate_size": sci(metrics.max_coordinate_size),
        "max_newton_iterations": metrics.max_newton_iterations,
        "max_step_cost": sci(metrics.max_step_cost),
        "min_step_margin": sci(metrics.min_step_margin),
        "min_root_separation_ratio": sci(metrics.min_root_separation_ratio),
        "endpoint_cost": sci(metrics.endpoint_cost),
        "endpoint_margin": sci(metrics.endpoint_margin),
    }


def semantic_digest(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()


def main() -> None:
    controls = (
        run_control("1e-3", 80, 90),
        run_control("1e-3", 160, 120),
        run_control("3e-4", 160, 120),
    )
    require(controls[0].permutations == controls[1].permutations, (
        "same-radius step/precision permutation mismatch",
        controls[0].histograms, controls[1].histograms
    ))
    structural = tuple(
        (
            run.histograms,
            run.exponents,
            run.normalized_a_histogram,
            run.normalized_b_histogram,
            run.c_histogram,
            run.c_orbit_count,
            run.c_order,
        )
        for run in controls
    )
    require(structural[0] == structural[1] == structural[2], (
        "radius/step structural mismatch", structural
    ))
    require(controls[0].c_order == EXPECTED_C_ORDER, (
        "C order hostile", controls[0].c_order, EXPECTED_C_ORDER
    ))
    record = {
        "depth": DEPTH,
        "controls": [run_record(run) for run in controls],
        "depth6_histogram": EXPECTED_DEPTH6_HISTOGRAM,
        "depth6_exponent": EXPECTED_DEPTH6_EXPONENT,
        "reflection_gauge_A": EXPECTED_REFLECTION_GAUGE_A_HISTOGRAM,
        "reflection_gauge_B": EXPECTED_REFLECTION_GAUGE_B_HISTOGRAM,
        "C": EXPECTED_C_HISTOGRAM,
        "C_orbit_count": EXPECTED_C_ORBIT_COUNT,
        "C_order": EXPECTED_C_ORDER,
        "status": "VERIFIED_NUMERICAL_BOUNDED_DEPTH6",
    }
    digest = semantic_digest(record)
    print("== fixed Keller old-L mpmath depth-six ancestry audit ==")
    for run in controls:
        row = run_record(run)
        print(
            f"control=radius:{run.radius},steps:{run.steps},dps:{run.dps};"
            f"max_step_cost={row['max_step_cost']};"
            f"min_step_margin={row['min_step_margin']};"
            f"max_forward_residual={row['max_forward_residual']};"
            f"max_root_residual={row['max_root_residual']}"
        )
        print(
            f"rouche=max_ratio:{row['max_rouche_ratio']},"
            f"max_projective_radius:{row['max_projective_radius']},"
            f"min_root_separation_ratio:{row['min_root_separation_ratio']};"
            f"endpoint=max_cost:{row['endpoint_cost']},"
            f"min_margin:{row['endpoint_margin']}"
        )
        print(
            f"conditioning=min_denominator_ratio:{row['min_denominator_ratio']},"
            f"max_coordinate_size:{row['max_coordinate_size']},"
            f"max_newton_iterations:{row['max_newton_iterations']}"
        )
        print(f"level_histograms={run.histograms};exponents={run.exponents}")
        print(
            f"root_action={run.root_action};raw_gauge_sections="
            f"A:{run.a_histogram},B:{run.b_histogram};"
            f"reflection_gauge=A:{run.normalized_a_histogram},"
            f"B:{run.normalized_b_histogram};"
            f"C=return_after_two_blocks:{run.c_histogram};"
            f"C_orbit_count={run.c_orbit_count};C_order={run.c_order}"
        )
    print(
        f"depth6_cycles={EXPECTED_DEPTH6_HISTOGRAM};"
        f"tame_permutation_exponent={EXPECTED_DEPTH6_EXPONENT}"
    )
    print(f"semantic_sha256={digest}")
    print(
        "section_gauge_boundary=raw A and B cycle types are not intrinsic;"
        "the displayed involutive A,B are a checked reflection-gauge normalization;"
        "the two-block return C is the intrinsic section datum"
    )
    print(
        "status=VERIFIED-NUMERICAL BOUNDED DEPTH6 SIDE-CAR ONLY;"
        "not an exact inertia theorem,not an all-level recurrence,"
        "not a coordinate-index statement"
    )
    if EXPECTED_SEMANTIC_SHA256:
        require(digest == EXPECTED_SEMANTIC_SHA256, ("semantic hash", digest))
    print("all high-precision continuation gates passed")


if __name__ == "__main__":
    main()
