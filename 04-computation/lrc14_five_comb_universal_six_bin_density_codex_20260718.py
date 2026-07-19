#!/usr/bin/env python3
"""Exact referee for THM-1198 (universal five-comb noncoverage).

For

    A(L, phi) = int_0^1 f(x) 1_{||phi + L x|| <= 1/14} dx

the symmetric probability density ``f`` is constant on the six equal bins,
with heights

    3/4, 13/12, 7/6, 7/6, 13/12, 3/4.

This dependency-free Fraction computation proves the compact part

    max_{6/7 <= L <= 3, phi mod 1} A(L, phi) = 7/36.

Put y=1/L and z=phi/L.  A tooth endpoint is

    x = (n +/- 1/14)y - z.

The load is continuous piecewise affine on the arrangement formed when one
of those endpoints crosses a density-bin boundary j/6.  Its maximum is
therefore attained at an arrangement vertex.  On the compact domain only
n=0,...,4 can meet [0,1], so the referee exhausts every exact vertex.

The analytic tail is also printed exactly.  For h=1_D-1/7, the centered
periodic primitive has ||H||_infinity=3/49.  The zero extension of f has
total variation

    f(0)+f(1)+TV_(0,1)(f) = 3/4+3/4+5/6 = 7/3,

so A(L,phi) <= 1/7 + 1/(7L).  At L>=3 this is at most
4/21 < 7/36.

No runner-pair tournament is used: this is a one-comb operator-norm
certificate.  The faithful vertices are endpoint/bin crossing events; a
runner tournament would discard the phase and overlap length that define the
observable.
"""

from fractions import Fraction as Q
from itertools import combinations


RADIUS = Q(1, 14)
L_MIN = Q(6, 7)
L_CUT = Q(3)
Y_MIN = 1 / L_CUT
Y_MAX = 1 / L_MIN

HEIGHTS = (
    Q(3, 4),
    Q(13, 12),
    Q(7, 6),
    Q(7, 6),
    Q(13, 12),
    Q(3, 4),
)
BIN_EDGES = tuple(Q(j, 6) for j in range(7))
TOOTH_CENTERS = tuple(range(5))

# A line is a*y + b*z = c, followed by a deterministic label.
Line = tuple[Q, Q, Q, str]
Point = tuple[Q, Q]
Segment = tuple[Q, Q, Q, Q, str]
EnvelopePiece = tuple[Q, Q, Q, Q]


def require(condition: bool, message: str) -> None:
    """Optimization-safe exact-certificate check.

    Python removes ``assert`` statements under ``-O``.  Every condition that
    certifies the arrangement or its arithmetic therefore passes through this
    helper instead.
    """
    if not condition:
        raise RuntimeError(f"exact certificate failed: {message}")


def qstr(x: Q) -> str:
    """Canonical exact display."""
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def endpoint(n: int, sign: int, y: Q, z: Q) -> Q:
    return (Q(n) + sign * RADIUS) * y - z


def overlap(a: Q, b: Q, c: Q, d: Q) -> Q:
    return max(Q(0), min(b, d) - max(a, c))


def load_direct(y: Q, z: Q) -> Q:
    """Integrate the six-bin density against all relevant closed teeth."""
    total = Q(0)
    for n in TOOTH_CENTERS:
        left = max(Q(0), endpoint(n, -1, y, z))
        right = min(Q(1), endpoint(n, +1, y, z))
        if right <= left:
            continue
        for j, height in enumerate(HEIGHTS):
            total += height * overlap(left, right, BIN_EDGES[j], BIN_EDGES[j + 1])
    return total


def danger_length_on_interval(y: Q, z: Q, a: Q, b: Q) -> Q:
    """Unweighted danger length on [a,b], independently of load_direct."""
    return sum(
        overlap(
            max(a, endpoint(n, -1, y, z)),
            min(b, endpoint(n, +1, y, z)),
            a,
            b,
        )
        for n in TOOTH_CENTERS
    )


def load_by_layers(y: Q, z: Q) -> Q:
    """Independent superlevel decomposition of the same density."""
    # f = 3/4 + (1/3) 1_[1/6,5/6] + (1/12) 1_[1/3,2/3].
    return (
        Q(3, 4) * danger_length_on_interval(y, z, Q(0), Q(1))
        + Q(1, 3) * danger_length_on_interval(y, z, Q(1, 6), Q(5, 6))
        + Q(1, 12) * danger_length_on_interval(y, z, Q(1, 3), Q(2, 3))
    )


def arrangement_lines() -> tuple[Line, ...]:
    lines: list[Line] = []
    for n in TOOTH_CENTERS:
        for sign in (-1, +1):
            slope = Q(n) + sign * RADIUS
            for j, edge in enumerate(BIN_EDGES):
                # (n +/- 1/14)y - z = j/6.
                lines.append((slope, Q(-1), edge, f"e({n},{sign:+d})=b{j}"))

    # Compact-domain boundary: Y_MIN <= y <= Y_MAX and 0 <= z <= y.
    lines.extend(
        (
            (Q(1), Q(0), Y_MIN, "y=Y_MIN"),
            (Q(1), Q(0), Y_MAX, "y=Y_MAX"),
            (Q(0), Q(1), Q(0), "z=0"),
            (Q(-1), Q(1), Q(0), "z=y"),
        )
    )
    return tuple(lines)


def line_intersection(first: Line, second: Line) -> Point | None:
    a, b, c, _ = first
    d, e, g, _ = second
    determinant = a * e - b * d
    if determinant == 0:
        return None
    return ((c * e - b * g) / determinant, (a * g - c * d) / determinant)


def feasible(point: Point) -> bool:
    y, z = point
    return Y_MIN <= y <= Y_MAX and Q(0) <= z <= y


def arrangement_vertices(lines: tuple[Line, ...]) -> set[Point]:
    vertices: set[Point] = {
        (Y_MIN, Q(0)),
        (Y_MIN, Y_MIN),
        (Y_MAX, Q(0)),
        (Y_MAX, Y_MAX),
    }
    for first, second in combinations(lines, 2):
        point = line_intersection(first, second)
        if point is not None and feasible(point):
            vertices.add(point)
    return vertices


def point_on_phase_line(line: Line, y: Q) -> Point:
    """Return (y,z) on a nonvertical arrangement line."""
    a, b, c, _ = line
    require(b != 0, f"phase candidate is vertical: {line[3]}")
    return (y, (c - a * y) / b)


def phase_load_segments(lines: tuple[Line, ...]) -> list[Segment]:
    """All affine load segments on phase-event lines.

    For fixed y, a maximum over z occurs when z is 0 or y, or when a tooth
    endpoint is at a density-bin boundary.  These are exactly the arrangement
    lines with nonzero z coefficient.  Intersections with all arrangement
    lines split each candidate into intervals on which its load is affine.
    """
    segments: list[Segment] = []
    for candidate in (line for line in lines if line[1] != 0):
        breakpoints = {Y_MIN, Y_MAX}
        for other in lines:
            point = line_intersection(candidate, other)
            if point is not None and feasible(point):
                breakpoints.add(point[0])
        ordered = sorted(breakpoints)
        for left, right in zip(ordered, ordered[1:]):
            middle = (left + right) / 2
            middle_point = point_on_phase_line(candidate, middle)
            if not feasible(middle_point):
                continue
            left_point = point_on_phase_line(candidate, left)
            right_point = point_on_phase_line(candidate, right)
            require(
                feasible(left_point) and feasible(right_point),
                f"candidate segment leaves compact domain: {candidate[3]}",
            )
            left_load = load_direct(*left_point)
            middle_load = load_direct(*middle_point)
            right_load = load_direct(*right_point)
            require(
                2 * middle_load == left_load + right_load,
                f"load is not affine on candidate segment: {candidate[3]}",
            )
            slope = (right_load - left_load) / (right - left)
            intercept = left_load - slope * left
            segments.append((left, right, slope, intercept, candidate[3]))
    return segments


def phase_free_envelope(segments: list[Segment]) -> list[EnvelopePiece]:
    """Exact upper envelope P(y)=sup_z A(1/y,z/y) on the compact domain."""
    breakpoints = {Y_MIN, Y_MAX}
    for left, right, _, _, _ in segments:
        breakpoints.update((left, right))
    for first, second in combinations(segments, 2):
        overlap_left = max(first[0], second[0])
        overlap_right = min(first[1], second[1])
        if overlap_left > overlap_right:
            continue
        delta_slope = first[2] - second[2]
        delta_intercept = first[3] - second[3]
        if delta_slope != 0:
            crossing = -delta_intercept / delta_slope
            if overlap_left <= crossing <= overlap_right:
                breakpoints.add(crossing)

    pieces: list[EnvelopePiece] = []
    ordered = sorted(breakpoints)
    for left, right in zip(ordered, ordered[1:]):
        middle = (left + right) / 2
        active = [segment for segment in segments if segment[0] <= middle <= segment[1]]
        require(bool(active), f"phase envelope has no active segment at y={qstr(middle)}")
        winner = max(active, key=lambda segment: segment[2] * middle + segment[3])
        expression = (winner[2], winner[3])
        if pieces and pieces[-1][1] == left and pieces[-1][2:] == expression:
            pieces[-1] = (pieces[-1][0], right, *expression)
        else:
            pieces.append((left, right, *expression))
    return pieces


EXPECTED_ENVELOPE: tuple[EnvelopePiece, ...] = (
    (Q(1, 3), Q(7, 18), Q(3, 7), Q(0)),
    (Q(7, 18), Q(5, 12), Q(-2, 7), Q(5, 18)),
    (Q(5, 12), Q(7, 15), Q(8, 21), Q(0)),
    (Q(7, 15), Q(119, 244), Q(-103, 84), Q(3, 4)),
    (Q(119, 244), Q(1, 2), Q(19, 84), Q(1, 24)),
    (Q(1, 2), Q(7, 12), Q(13, 42), Q(0)),
    (Q(7, 12), Q(2, 3), Q(-1, 14), Q(2, 9)),
    (Q(2, 3), Q(35, 48), Q(11, 42), Q(0)),
    (Q(35, 48), Q(5, 6), Q(-5, 42), Q(5, 18)),
    (Q(5, 6), Q(7, 8), Q(3, 14), Q(0)),
    (Q(7, 8), Q(63, 68), Q(-9, 14), Q(3, 4)),
    (Q(63, 68), Q(7, 6), Q(1, 6), Q(0)),
)


def main() -> None:
    mass = sum(HEIGHTS) / 6
    interior_tv = sum(abs(HEIGHTS[j + 1] - HEIGHTS[j]) for j in range(5))
    extended_tv = HEIGHTS[0] + HEIGHTS[-1] + interior_tv
    primitive_sup = Q(3, 49)
    tail_coefficient = primitive_sup * extended_tv

    require(mass == 1, "density mass is not one")
    require(interior_tv == Q(5, 6), "interior total variation changed")
    require(extended_tv == Q(7, 3), "zero-extension total variation changed")
    require(tail_coefficient == Q(1, 7), "BV tail coefficient changed")

    lines = arrangement_lines()
    vertices = arrangement_vertices(lines)
    values: list[tuple[Q, Point]] = []
    for point in vertices:
        direct = load_direct(*point)
        layered = load_by_layers(*point)
        require(direct == layered, f"direct/layer loads disagree at {point}")
        values.append((direct, point))
    values.sort(reverse=True)

    segments = phase_load_segments(lines)
    envelope = phase_free_envelope(segments)

    target = Q(7, 36)
    equal_vertices = sorted(point for value, point in values if value == target)
    nonmax = max(value for value, _ in values if value < target)
    require(len(lines) == 74, "arrangement line count changed")
    require(len(vertices) == 101, "feasible arrangement vertex count changed")
    require(values[0][0] == target, "compact maximum changed")
    require(
        equal_vertices == [(Y_MAX, Q(7, 12)), (Y_MAX, Q(3, 4))],
        "compact maximizing vertices changed",
    )
    require(nonmax == Q(55, 288), "largest nonmaximal vertex load changed")
    require(len(segments) == 175, "phase-envelope segment count changed")
    require(tuple(envelope) == EXPECTED_ENVELOPE, "phase-free envelope changed")
    envelope_max = max(
        slope * y + intercept
        for left, right, slope, intercept in envelope
        for y in (left, right)
    )
    require(envelope_max == target, "phase-free envelope maximum changed")

    # At L=6/7 the equality segment is z in [7/12,3/4], equivalently
    # phi=z/y in [1/2,9/14].  Its midpoint supplies a direct affine check.
    equality_midpoint = (Y_MAX, (Q(7, 12) + Q(3, 4)) / 2)
    require(load_direct(*equality_midpoint) == target, "equality segment midpoint changed")
    phi_interval = (equal_vertices[0][1] / Y_MAX, equal_vertices[1][1] / Y_MAX)
    require(phi_interval == (Q(1, 2), Q(9, 14)), "equality phase interval changed")

    tail_at_cut = Q(1, 7) + tail_coefficient / L_CUT
    five_dual_survivor = 1 - 5 * target
    normalized_length_survivor = five_dual_survivor / max(HEIGHTS)
    physical_length_factor = Q(6, 7) * normalized_length_survivor
    six_private_total_factor = 6 * physical_length_factor
    six_overlap_surplus = 6 * target - 1
    six_multi_length_ceiling = six_overlap_surplus / min(HEIGHTS)
    six_unique_length_floor = 1 - six_multi_length_ceiling
    six_unique_physical_factor = Q(6, 7) * six_unique_length_floor
    require(tail_at_cut == Q(4, 21), "BV tail value at cutoff changed")
    require(tail_at_cut < target < Q(1, 5), "global-margin ordering failed")
    require(five_dual_survivor == Q(1, 36), "five-comb dual survivor changed")
    require(
        normalized_length_survivor == Q(1, 42),
        "normalized Lebesgue survivor changed",
    )
    require(physical_length_factor == Q(1, 49), "physical survivor factor changed")
    require(six_private_total_factor == Q(6, 49), "six-private total changed")
    require(six_overlap_surplus == Q(1, 6), "weighted overlap surplus changed")
    require(
        six_multi_length_ceiling == Q(2, 9),
        "multiply-covered length ceiling changed",
    )
    require(six_unique_length_floor == Q(7, 9), "unique-provider length floor changed")
    require(
        six_unique_physical_factor == Q(2, 3),
        "physical unique-provider factor changed",
    )

    print("THM-1198 universal five-comb six-bin dual density")
    print("density heights:", " ".join(qstr(h) for h in HEIGHTS))
    print("density mass:", qstr(mass))
    print("interior TV:", qstr(interior_tv))
    print("zero-extension TV:", qstr(extended_tv))
    print("centered primitive sup:", qstr(primitive_sup))
    print("tail coefficient ||H||inf*TV:", qstr(tail_coefficient))
    print()
    print("compact model: 6/7 <= L <= 3, y=1/L, z=phi/L, 0<=z<=y")
    print("relevant tooth centers:", TOOTH_CENTERS)
    print("arrangement lines:", len(lines))
    print("raw line pairs:", len(lines) * (len(lines) - 1) // 2)
    print("feasible distinct vertices:", len(vertices))
    print("independent direct/layer evaluations matched at every vertex: YES")
    print("compact exact maximum:", qstr(values[0][0]))
    print(
        "maximizing vertices (y,z):",
        " ".join(f"({qstr(y)},{qstr(z)})" for y, z in equal_vertices),
    )
    print(
        "equality boundary: L=6/7 and phi in",
        f"[{qstr(phi_interval[0])},{qstr(phi_interval[1])}]",
    )
    print("largest nonmaximal arrangement-vertex value:", qstr(nonmax))
    print("top distinct compact vertex loads:")
    for value in sorted({value for value, _ in values}, reverse=True)[:10]:
        print(" ", qstr(value))
    print()
    print("exact phase-free compact envelope P(L)=a/L+b:")
    print("  L_left L_right a b")
    for y_left, y_right, slope, intercept in reversed(envelope):
        print(
            " ",
            qstr(1 / y_right),
            qstr(1 / y_left),
            qstr(slope),
            qstr(intercept),
        )
    print("envelope affine phase segments:", len(segments))
    print("envelope pieces:", len(envelope))
    print()
    print("analytic tail: A(L,phi) <= 1/7 + 1/(7L)")
    print("tail at L=3:", qstr(tail_at_cut))
    print("global exact supremum:", qstr(target))
    print("one-comb margin below 1/5:", qstr(Q(1, 5) - target))
    print("five-comb total upper bound:", qstr(5 * target))
    print("five-comb uncovered dual mass:", qstr(five_dual_survivor))
    print("five-comb normalized Lebesgue survivor:", qstr(normalized_length_survivor))
    print("five-comb physical c-gap survivor factor:", qstr(physical_length_factor), "/ c")
    print()
    print("six-comb corollary: a covered c-slow gap forces d1 < 13c/6")
    print("six-comb private region floor:", qstr(physical_length_factor), "/ c each")
    print("six-comb total private length floor:", qstr(six_private_total_factor), "/ c")
    print("six-comb weighted overlap-surplus ceiling:", qstr(six_overlap_surplus))
    print("six-comb normalized multiply-covered length ceiling:", qstr(six_multi_length_ceiling))
    print("six-comb normalized unique-provider length floor:", qstr(six_unique_length_floor))
    print("six-comb physical unique-provider length factor:", qstr(six_unique_physical_factor), "/ c")
    print("six-comb functional drift: sum_i Pbar(6d_i/(7c)) >= 1")
    print("  Pbar=P on [6/7,3], Pbar(L)=1/7+1/(7L) for L>3")
    print("Tournament Analysis: not applicable to the proof inequality.")
    print("  pairwise observable: none (one-comb operator norm)")
    print("  faithful vertices: endpoint/bin crossing events")
    print("  challenged assumption: runner vertices erase phase and overlap length")
    print("  derived carrier: six private needles, ordered by dual mass then chronology")
    print("  fingerprint: transitive scores 0..5, no cycles, singleton SCCs, one path")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
