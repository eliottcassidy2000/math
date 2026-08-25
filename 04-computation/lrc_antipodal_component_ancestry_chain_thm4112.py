#!/usr/bin/env python3
"""Fraction-exact primary referee for THM-4112.

The theorem is analytic; this finite referee freezes its exact algebra and
the boundary/hostile controls used to audit the statement.  Danger teeth are
literal open intervals.  Consequently two teeth merge only under strict
overlap: endpoint contact remains a safe separator.  Every gate uses
``require`` and remains active under ``python -O``.

Its finite universes are the explicit AP7 four-outlier bank and the named
chain/hostile rows below.
"""

from __future__ import annotations

from fractions import Fraction as F
from hashlib import sha256
from heapq import merge
from itertools import combinations
from math import gcd


DELTA = F(1, 14)

AP7_CORE = tuple(range(1, 8))
AP7_J = (F(4, 35), F(13, 98))
AP7_L = F(9, 490)

D0_CORE = (3, 4, 5, 6, 8, 10, 12)
D0_J = (F(1, 42), F(13, 168))
D0_L = F(3, 56)

AP8_CORE = tuple(range(1, 9))
AP8_J = (F(11, 49), F(13, 56))
AP8_L = F(3, 392)

# Filled after the semantic payload was first frozen.  It deliberately hashes
# theorem-bearing exact data, not presentation prose or source bytes.
EXPECTED_SEMANTIC_SHA256 = "9834738b13362e78c30f577bd4dc1219221c9c1d8413a18102dc1b56d65cf093"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def omega(speed: int) -> int:
    """Number of distinct antipodal danger phases for this parity."""
    return 1 if speed % 2 == 0 else 2


def circle_norm(value: F) -> F:
    residue = value % 1
    return min(residue, 1 - residue)


def literal_safe(speed: int, theta: F) -> bool:
    return all(
        circle_norm(speed * (theta + F(half, 2))) >= DELTA
        for half in (0, 1)
    )


def tooth_width(speed: int) -> F:
    return F(1, 7 * speed)


def parent_gap(speed: int) -> F:
    """Gap between consecutive open teeth in the parity-reduced comb."""
    return F(1, omega(speed) * speed) - tooth_width(speed)


def pair_envelope(lower: int, upper: int) -> F:
    require(lower < upper, "pair speeds are not strictly ordered")
    return min(
        F(2, 7 * lower),
        tooth_width(lower) + F(2, 7 * upper),
    )


def triple_envelope(lower: int, middle: int, upper: int) -> F:
    require(lower < middle < upper, "triple speeds are not strictly ordered")
    return tooth_width(lower) + 2 * pair_envelope(middle, upper)


def four_envelope(root: int, residual: tuple[int, int, int]) -> tuple[F, F]:
    lower, middle, upper = sorted(residual)
    child = triple_envelope(lower, middle, upper)
    return child, tooth_width(root) + 2 * child


def ap7_root_choice(outliers: tuple[int, int, int, int]) -> tuple[int, str]:
    """Root split for a<b<c<d with a>=84 and at least one of a,b even."""
    lower, second, _, _ = outliers
    require(all(left < right for left, right in zip(outliers, outliers[1:])), "AP7 row is unordered")
    require(lower >= 84, "AP7 root split used below its threshold")
    require(lower % 2 == 0 or second % 2 == 0, "AP7 root split has no even a or b")
    if lower % 2 == 0:
        return lower, "a-even"
    if second < 2 * lower:
        return second, "b-even-near"
    return lower, "b-even-separated"


def chain_coefficients(depth: int) -> tuple[int, ...]:
    """Symbolic coefficients produced by h_i=A_i+2h_{i+1}."""
    require(depth >= 1, "chain depth must be positive")
    coefficients = (1,)
    for _ in range(depth - 1):
        coefficients = (1,) + tuple(2 * value for value in coefficients)
    return coefficients


def chain_levels(speeds: tuple[int, ...]) -> tuple[F, ...]:
    require(bool(speeds), "empty component chain")
    require(
        all(left < right for left, right in zip(speeds, speeds[1:])),
        "component chain is not strictly ordered",
    )
    reversed_levels = [tooth_width(speeds[-1])]
    for speed in reversed(speeds[:-1]):
        reversed_levels.append(tooth_width(speed) + 2 * reversed_levels[-1])
    return tuple(reversed(reversed_levels))


def expanded_chain_levels(speeds: tuple[int, ...]) -> tuple[F, ...]:
    return tuple(
        sum(
            (F(2 ** (index - start), 7 * speeds[index]) for index in range(start, len(speeds))),
            F(0),
        )
        for start in range(len(speeds))
    )


def geometric_chain_constants(scale: F) -> tuple[F, F]:
    """Coefficients in h_{i+1}<=tail/v_i and h_i<=head/v_i."""
    require(scale > 2, "geometric scale must exceed two")
    tail = F(1, 7) / (scale - 2)
    head = scale / (7 * (scale - 2))
    require(
        F(1, 7) * (F(1, 1) / scale) / (1 - F(2, 1) / scale) == tail,
        "infinite geometric tail identity changed",
    )
    require(tooth_width(1) + 2 * tail == head, "geometric head identity changed")
    return tail, head


def lifted_teeth(speed: int, left: F = F(-1), right: F = F(2)) -> tuple[tuple[F, F], ...]:
    """All literal open teeth meeting a lifted interval, in endpoint order."""
    count = omega(speed) * speed
    radius = F(1, 14 * speed)
    lower_index = count * (left - radius)
    upper_index = count * (right + radius)
    first = lower_index.numerator // lower_index.denominator - 2
    last = upper_index.numerator // upper_index.denominator + 3
    return tuple(
        (F(index, count) - radius, F(index, count) + radius)
        for index in range(first, last + 1)
        if F(index, count) + radius > left
        and F(index, count) - radius < right
    )


def circular_components(speeds: tuple[int, ...]) -> tuple[tuple[F, F, F], ...]:
    """One full lifted representative of every literal circular component."""
    require(bool(speeds), "empty component family")
    merged_components: list[tuple[F, F]] = []
    active_left: F | None = None
    active_right: F | None = None
    for left, right in merge(*(lifted_teeth(speed) for speed in speeds)):
        if active_left is not None and left < active_right:
            active_right = max(active_right, right)
        else:
            if active_left is not None:
                merged_components.append((active_left, active_right))
            active_left, active_right = left, right
    if active_left is not None:
        merged_components.append((active_left, active_right))

    # Lifted component starts repeat with period one.  Selecting starts in
    # [0,1) counts a seam component once without truncating its length.
    return tuple(
        (left, right, right - left)
        for left, right in merged_components
        if 0 <= left < 1
    )


def overlapping_teeth(speed: int, component: tuple[F, F, F]) -> tuple[tuple[F, F], ...]:
    left, right, _ = component
    return tuple(
        tooth
        for tooth in lifted_teeth(speed)
        if tooth[0] < right and tooth[1] > left
    )


def even_presentation_denominator(value: F) -> int:
    return value.denominator if value.denominator % 2 == 0 else 2 * value.denominator


def arrangement_witness(
    outliers: tuple[int, ...],
    core: tuple[int, ...],
    window: tuple[F, F],
) -> tuple[F, int]:
    walls = {window[0], window[1]}
    for speed in outliers:
        for left, right in lifted_teeth(speed, *window):
            if window[0] <= left <= window[1]:
                walls.add(left)
            if window[0] <= right <= window[1]:
                walls.add(right)
    body = core + outliers
    for theta in sorted(walls):
        if all(literal_safe(speed, theta) for speed in body):
            return theta, even_presentation_denominator(theta)
    raise RuntimeError(f"no arrangement witness for {body}")


def ceil_fraction(value: F) -> int:
    return -((-value.numerator) // value.denominator)


def norm_range_on_interval(speed: int, window: tuple[F, F]) -> tuple[F, F]:
    """Exact range of ||speed*theta|| from all affine extrema."""
    left, right = (speed * endpoint for endpoint in window)
    doubled_first = ceil_fraction(2 * left)
    doubled_last = (2 * right).numerator // (2 * right).denominator
    candidates = {left, right}
    candidates.update(F(index, 2) for index in range(doubled_first, doubled_last + 1))
    values = tuple(circle_norm(candidate) for candidate in candidates)
    return min(values), max(values)


def is_arithmetic_progression(speeds: tuple[int, ...]) -> bool:
    if len(speeds) <= 2:
        return True
    step = speeds[1] - speeds[0]
    return step > 0 and all(
        speeds[index + 1] - speeds[index] == step
        for index in range(1, len(speeds) - 1)
    )


def common_gcd(speeds: tuple[int, ...]) -> int:
    value = 0
    for speed in speeds:
        value = gcd(value, speed)
    return value


def canonical(value: object) -> str:
    """Version-stable encoding of the exact semantic payload."""
    if isinstance(value, F):
        return f"Q({value.numerator},{value.denominator})"
    if isinstance(value, tuple):
        return "T(" + ",".join(canonical(item) for item in value) + ")"
    if isinstance(value, str):
        return "S(" + value.replace("\\", "\\\\").replace(")", "\\)") + ")"
    if isinstance(value, bool):
        return "B(1)" if value else "B(0)"
    if isinstance(value, int):
        return f"Z({value})"
    raise TypeError(f"unsupported semantic value: {type(value)!r}")


def main() -> None:
    # Symbolic and Fraction-exact recursion algebra.  The final one-comb level
    # h_r=A_r is an equality; only a reinserted parent gains the strict
    # component inequality supplied by strict-open overlap.
    coefficient_rows = tuple(chain_coefficients(depth) for depth in range(1, 13))
    require(
        coefficient_rows == tuple(tuple(2**index for index in range(depth)) for depth in range(1, 13)),
        "symbolic chain coefficients changed",
    )
    algebra_numeric_rows = 0
    for depth in range(2, 7):
        for speeds in combinations(range(2, 13), depth):
            require(chain_levels(speeds) == expanded_chain_levels(speeds), "expanded h recursion failed")
            algebra_numeric_rows += 1
    require(algebra_numeric_rows == 1474, "numeric chain algebra census changed")
    base_component_count = 0
    for speed in range(1, 101):
        expected_gap = F(6, 7 * speed) if speed % 2 == 0 else F(5, 14 * speed)
        require(parent_gap(speed) == expected_gap, "parity parent-gap formula changed")
        base_components = circular_components((speed,))
        require(len(base_components) == omega(speed) * speed, "one-comb component count changed")
        require(
            all(component[2] == tooth_width(speed) for component in base_components),
            "one-comb base span is not equality",
        )
        base_component_count += len(base_components)
    require(base_component_count == 7550, "one-comb equality census changed")

    # Arbitrary finite depth under geometric scale separation.  The displayed
    # tail is the infinite geometric majorant, so no finite-depth scan is used
    # to justify the corollary.
    scale_12_5 = F(12, 5)
    tail_12_5, head_12_5 = geometric_chain_constants(scale_12_5)
    require(
        (tail_12_5, head_12_5, head_12_5 / D0_L, head_12_5 / AP8_L)
        == (F(5, 14), F(6, 7), F(16), F(112)),
        "scale-12/5 arbitrary-depth constants changed",
    )
    require(tail_12_5 <= F(5, 14), "scale-12/5 does not fit the odd parent gap")
    scale_3 = F(3)
    tail_3, head_3 = geometric_chain_constants(scale_3)
    require(
        (tail_3, head_3, head_3 / AP8_L) == (F(1, 7), F(3, 7), F(56)),
        "scale-three arbitrary-depth constants changed",
    )
    require(tail_3 <= F(5, 14), "scale-three does not fit the odd parent gap")

    # AP7 core and four-outlier boundary.  The row at 82 is only a criterion
    # near miss; it is not asserted to be an actual covering family.
    require(AP7_J[1] - AP7_J[0] == AP7_L, "AP7 interval length changed")
    require(
        all(not lifted_teeth(speed, *AP7_J) for speed in AP7_CORE),
        "AP7 interval lost literal safety",
    )
    child84, bound84 = four_envelope(84, (85, 86, 87))
    child82, bound82 = four_envelope(82, (83, 84, 85))
    require(
        (pair_envelope(86, 87), child84, parent_gap(84), bound84, AP7_L - bound84)
        == (F(1, 301), F(213, 25585), F(1, 98), F(39439, 2149140), F(1, 61404)),
        "AP7 boundary-84 arithmetic changed",
    )
    require(child84 <= parent_gap(84) and bound84 <= AP7_L, "AP7 boundary-84 gate failed")
    require(child82 <= parent_gap(82) and bound82 > AP7_L, "AP7 boundary-82 criterion changed")
    require(bound82 - AP7_L == F(151, 357315), "AP7 boundary-82 deficit changed")
    witness84 = arrangement_witness((84, 85, 86, 87), AP7_CORE, AP7_J)
    require(witness84 == (F(74, 609), 1218), "AP7 boundary witness changed")

    # Exact finite positive bank around the threshold, using the theorem's
    # parity/scale root split.  If a is odd and b is even, root b only while
    # b<2a; once b>=2a the safe root is a.  Every literal component and every
    # arrangement witness is checked.
    ap7_positive_rows = 0
    ap7_root_case_counts = {"a-even": 0, "b-even-near": 0, "b-even-separated": 0}
    ap7_maximum_ratio = F(0)
    ap7_maximum_row: tuple[object, ...] | None = None
    for outliers in combinations(range(84, 101), 4):
        lower, second, _, _ = outliers
        if lower % 2 and second % 2:
            continue
        root, root_case = ap7_root_choice(outliers)
        residual = tuple(speed for speed in outliers if speed != root)
        child, bound = four_envelope(root, residual)
        require(child <= parent_gap(root), f"AP7 parent gap failed at {outliers}")
        require(bound <= AP7_L, f"AP7 core-length gate failed at {outliers}")
        components = circular_components(outliers)
        longest = max(components, key=lambda component: component[2])
        require(longest[2] < bound, f"AP7 literal component bound failed at {outliers}")
        arrangement_witness(outliers, AP7_CORE, AP7_J)
        ratio = longest[2] / bound
        if ratio > ap7_maximum_ratio:
            ap7_maximum_ratio = ratio
            ap7_maximum_row = (outliers, longest, bound)
        ap7_root_case_counts[root_case] += 1
        ap7_positive_rows += 1
    expected_ap7_maximum = (
        (88, 95, 97, 99),
        (F(251, 1330), F(137, 693), F(1181, 131670)),
        F(93167, 5676440),
    )
    require(ap7_positive_rows == 1932, "AP7 positive-row census changed")
    require(
        tuple(ap7_root_case_counts[key] for key in ("a-even", "b-even-near", "b-even-separated"))
        == (1344, 588, 0),
        "AP7 root-case census changed",
    )
    require(ap7_maximum_ratio == F(458228, 838503), "AP7 maximum ratio changed")
    require(ap7_maximum_row == expected_ap7_maximum, "AP7 maximum row changed")

    # Exact representatives of both odd-a/even-b branches.  The far row also
    # freezes the sharp warning that rooting at b is invalid when b is widely
    # separated, even though the correct root a succeeds.
    near_split_row = (85, 86, 87, 88)
    near_root, near_case = ap7_root_choice(near_split_row)
    near_child, near_bound = four_envelope(near_root, tuple(v for v in near_split_row if v != near_root))
    require(
        (near_root, near_case, near_child, parent_gap(near_root), near_bound, AP7_L - near_bound)
        == (86, "b-even-near", F(61, 7395), F(3, 301), F(80839, 4451790), F(650, 3116253)),
        "near odd-a/even-b root split changed",
    )
    require(
        arrangement_witness(near_split_row, AP7_CORE, AP7_J) == (F(141, 1190), 1190),
        "near root-split witness changed",
    )
    separated_split_row = (85, 170, 171, 172)
    separated_root, separated_case = ap7_root_choice(separated_split_row)
    separated_child, separated_bound = four_envelope(
        separated_root,
        tuple(v for v in separated_split_row if v != separated_root),
    )
    require(
        (
            separated_root,
            separated_case,
            separated_child,
            parent_gap(separated_root),
            separated_bound,
            AP7_L - separated_bound,
        )
        == (85, "b-even-separated", F(851, 203490), F(1, 238), F(146, 14535), F(2371, 284886)),
        "separated odd-a/even-b root split changed",
    )
    require(
        arrangement_witness(separated_split_row, AP7_CORE, AP7_J) == (F(137, 1197), 2394),
        "separated root-split witness changed",
    )
    far_split_row = (85, 1000, 1001, 1002)
    far_root, far_case = ap7_root_choice(far_split_row)
    far_child, far_bound = four_envelope(far_root, tuple(v for v in far_split_row if v != far_root))
    wrong_far_child, wrong_far_bound = four_envelope(1000, (85, 1001, 1002))
    require(
        (far_root, far_case, far_child, far_bound, AP7_L - far_bound)
        == (85, "b-even-separated", F(5001, 7007000), F(185117, 59559500), F(908833, 59559500)),
        "far separated correct-root control changed",
    )
    require(
        (wrong_far_child, parent_gap(1000), wrong_far_child - parent_gap(1000), wrong_far_bound)
        == (F(1341, 595595), F(3, 3500), F(83049, 59559500), F(553417, 119119000)),
        "far separated wrong-root hostile changed",
    )
    require(
        arrangement_witness(far_split_row, AP7_CORE, AP7_J) == (F(4, 35), 70),
        "far root-split witness changed",
    )

    # Two scalar-ceiling crossings: an adaptive weight-seven root with no
    # residual odd gap 4, 8, or 12, and a fully odd weight-eight row.
    adaptive_row = (85, 86, 91, 101)
    adaptive_child, adaptive_bound = four_envelope(86, (85, 91, 101))
    adaptive_witness = arrangement_witness(adaptive_row, AP7_CORE, AP7_J)
    require(
        (
            pair_envelope(91, 101),
            adaptive_child,
            parent_gap(86),
            adaptive_bound,
            AP7_L - adaptive_bound,
            adaptive_witness,
        )
        == (
            F(2, 637),
            F(431, 54145),
            F(3, 301),
            F(81867, 4656470),
            F(366, 465647),
            (F(81, 707), 1414),
        ),
        "adaptive AP7 weight-seven control changed",
    )
    require(
        tuple(sorted((91 - 85, 101 - 91, 101 - 85))) == (6, 10, 16),
        "adaptive residual gap signature changed",
    )
    odd_row = (47, 95, 97, 99)
    odd_child, odd_bound = four_envelope(47, (95, 97, 99))
    odd_witness = arrangement_witness(odd_row, AP7_CORE, AP7_J)
    require(
        (
            pair_envelope(97, 99),
            odd_child,
            parent_gap(47),
            odd_bound,
            AP7_L - odd_bound,
            odd_witness,
        )
        == (
            F(2, 679),
            F(477, 64505),
            F(5, 658),
            F(54053, 3031735),
            F(22847, 42444290),
            (F(4, 35), 70),
        ),
        "all-odd AP7 weight-eight control changed",
    )

    # Parent-gap hostile: the residual bridges adjacent U_85 gaps, the final
    # union meets three root teeth, and the naive formal envelope is false.
    gap_hostile = (85, 87, 89, 91)
    hostile_child, hostile_bound = four_envelope(85, (87, 89, 91))
    require(
        (
            pair_envelope(89, 91),
            hostile_child,
            parent_gap(85),
            hostile_bound,
        )
        == (F(2, 623), F(437, 54201), F(1, 238), F(11719, 658155)),
        "parent-gap hostile arithmetic changed",
    )
    require(hostile_child > parent_gap(85), "parent-gap hostile now passes its missing hypothesis")
    residual_components = circular_components((87, 89, 91))
    residual_bridge_endpoints = (
        (F(38, 637), F(13, 203)),
        (F(83, 1274), F(85, 1218)),
    )
    residual_bridges = tuple(
        component
        for component in residual_components
        if component[:2] in residual_bridge_endpoints
    )
    require(
        residual_bridges
        == (
            (F(38, 637), F(13, 203), F(81, 18473)),
            (F(83, 1274), F(85, 1218), F(257, 55419)),
        ),
        "parent-gap hostile residual bridges changed",
    )
    require(
        all(len(overlapping_teeth(85, component)) >= 2 for component in residual_bridges),
        "parent-gap hostile residual lost a two-tooth bridge",
    )
    hostile_component = max(circular_components(gap_hostile), key=lambda component: component[2])
    require(
        hostile_component == (F(69, 1274), F(46, 637), F(23, 1274)),
        "parent-gap hostile component changed",
    )
    require(hostile_component[2] > hostile_bound, "naive parent recursion no longer fails")
    require(
        hostile_component[2] - hostile_bound == F(207559, 838489470),
        "parent-gap hostile excess changed",
    )
    require(len(overlapping_teeth(85, hostile_component)) == 3, "three-root-tooth ancestry changed")

    # Strong method hostile: every possible choice of top root has a residual
    # component meeting two consecutive root teeth.  Yet AP7 survives at its
    # left endpoint, so this is an ancestry-method obstruction, not a cover.
    all_root_hostile = (43, 45, 47, 49)
    expected_all_root = (
        (
            43,
            (F(41, 686), F(43, 630), F(131, 15435)),
            ((F(17, 301), F(18, 301)), (F(41, 602), F(1, 14))),
        ),
        (
            45,
            (F(17, 301), F(43, 658), F(251, 28294)),
            ((F(17, 315), F(2, 35)), (F(41, 630), F(43, 630))),
        ),
        (
            47,
            (F(17, 315), F(43, 686), F(269, 30870)),
            ((F(17, 329), F(18, 329)), (F(41, 658), F(43, 658))),
        ),
        (
            49,
            (F(17, 329), F(18, 301), F(115, 14147)),
            ((F(17, 343), F(18, 343)), (F(41, 686), F(43, 686))),
        ),
    )
    all_root_rows: list[tuple[int, tuple[F, F, F], tuple[tuple[F, F], ...]]] = []
    for root in all_root_hostile:
        residual = tuple(speed for speed in all_root_hostile if speed != root)
        bridges = tuple(
            (component, overlapping_teeth(root, component))
            for component in circular_components(residual)
            if len(overlapping_teeth(root, component)) >= 2
        )
        require(bool(bridges), f"all-root hostile lost bridge at root {root}")
        component, hits = bridges[0]
        all_root_rows.append((root, component, hits[:2]))
    require(tuple(all_root_rows) == expected_all_root, "all-root hostile exact bridges changed")
    require(
        all(literal_safe(speed, AP7_J[0]) for speed in AP7_CORE + all_root_hostile),
        "all-root method hostile lost its AP7 survivor",
    )

    # Non-AP safe core D0.  Exact affine extrema prove the whole interval,
    # rather than merely sampling endpoints or scanning a clock.
    require(D0_J[1] - D0_J[0] == D0_L, "D0 interval length changed")
    require(not is_arithmetic_progression(D0_CORE), "D0 became an arithmetic progression")
    require(common_gcd(D0_CORE) == 1, "D0 lost primitivity")
    expected_d0_ranges = (
        (3, F(1, 14), F(13, 56)),
        (4, F(2, 21), F(13, 42)),
        (5, F(5, 42), F(65, 168)),
        (6, F(1, 7), F(13, 28)),
        (8, F(4, 21), F(1, 2)),
        (10, F(19, 84), F(1, 2)),
        (12, F(1, 14), F(1, 2)),
    )
    d0_ranges = tuple((speed,) + norm_range_on_interval(speed, D0_J) for speed in D0_CORE)
    require(d0_ranges == expected_d0_ranges, "D0 exact norm ranges changed")
    for speed, minimum, maximum in d0_ranges:
        require(minimum >= DELTA, f"D0 lower range safety failed at speed {speed}")
        if speed % 2:
            require(maximum <= F(3, 7), f"D0 odd upper range safety failed at speed {speed}")
    require(
        all(not lifted_teeth(speed, *D0_J) for speed in D0_CORE),
        "D0 interval and literal-tooth proofs disagree",
    )

    # Four-outlier controls over D0, including the adjacent criterion near
    # miss and an all-odd weight-eight crossing.
    d0_child28, d0_bound28 = four_envelope(28, (29, 30, 31))
    d0_child26, d0_bound26 = four_envelope(26, (27, 28, 29))
    require(d0_child28 <= parent_gap(28) and d0_bound28 <= D0_L, "D0 boundary-28 gate failed")
    require(D0_L - d0_bound28 == F(89, 170520), "D0 boundary-28 slack changed")
    require(d0_child26 <= parent_gap(26) and d0_bound26 > D0_L, "D0 boundary-26 criterion changed")
    require(d0_bound26 - D0_L == F(457, 137592), "D0 boundary-26 deficit changed")
    require(
        arrangement_witness((28, 29, 30, 31), D0_CORE, D0_J) == (F(1, 42), 42),
        "D0 boundary-28 witness changed",
    )
    d0_odd_row = (15, 31, 33, 35)
    d0_odd_child, d0_odd_bound = four_envelope(15, (31, 33, 35))
    require(d0_odd_child <= parent_gap(15) and d0_odd_bound < D0_L, "D0 all-odd gate failed")
    require(D0_L - d0_odd_bound == F(19, 95480), "D0 all-odd slack changed")
    require(
        arrangement_witness(d0_odd_row, D0_CORE, D0_J) == (F(1, 42), 42),
        "D0 all-odd witness changed",
    )

    # D0 plus six doubling outliers: exact equality h_1=L0 is admissible
    # because at least one strict-open parent insertion makes the final
    # component shorter than h_1.
    d0_six = (16, 32, 64, 128, 256, 512)
    d0_six_levels = chain_levels(d0_six)
    expected_d0_six_levels = (F(3, 56), F(5, 224), F(1, 112), F(3, 896), F(1, 896), F(1, 3584))
    require(d0_six_levels == expected_d0_six_levels, "D0 six-chain levels changed")
    require(d0_six_levels[0] == D0_L, "D0 six-chain equality boundary changed")
    require(
        all(d0_six_levels[index + 1] <= parent_gap(d0_six[index]) for index in range(5)),
        "D0 six-chain parent gap failed",
    )
    d0_six_components = circular_components(d0_six)
    d0_six_maximum = max(component[2] for component in d0_six_components)
    require(
        (len(d0_six_components), d0_six_maximum) == (416, F(1, 112)),
        "D0 six-chain literal component control changed",
    )
    d0_six_witness = arrangement_witness(d0_six, D0_CORE, D0_J)
    require(d0_six_witness == (F(43, 1792), 1792), "D0 six-chain witness changed")

    # AP8 plus five doubling outliers: the strict positive-slack boundary row.
    ap8_five = (94, 188, 376, 752, 1504)
    ap8_five_levels = chain_levels(ap8_five)
    expected_ap8_five_levels = (F(5, 658), F(1, 329), F(3, 2632), F(1, 2632), F(1, 10528))
    require(AP8_J[1] - AP8_J[0] == AP8_L, "AP8 interval length changed")
    require(
        all(not lifted_teeth(speed, *AP8_J) for speed in AP8_CORE),
        "AP8 interval lost literal safety",
    )
    require(ap8_five_levels == expected_ap8_five_levels, "AP8 five-chain levels changed")
    require(AP8_L - ap8_five_levels[0] == F(1, 18424), "AP8 five-chain slack changed")
    require(
        all(ap8_five_levels[index + 1] <= parent_gap(ap8_five[index]) for index in range(4)),
        "AP8 five-chain parent gap failed",
    )
    ap8_five_components = circular_components(ap8_five)
    ap8_five_maximum = max(component[2] for component in ap8_five_components)
    require(
        (len(ap8_five_components), ap8_five_maximum) == (1316, F(1, 658)),
        "AP8 five-chain literal component control changed",
    )
    ap8_five_witness = arrangement_witness(ap8_five, AP8_CORE, AP8_J)
    require(ap8_five_witness == (F(11, 49), 98), "AP8 five-chain witness changed")

    semantic_payload = (
        ("chain-coefficients", coefficient_rows),
        ("chain-algebra", algebra_numeric_rows, base_component_count),
        (
            "geometric-arbitrary-depth",
            scale_12_5,
            tail_12_5,
            head_12_5,
            head_12_5 / D0_L,
            head_12_5 / AP8_L,
            scale_3,
            tail_3,
            head_3,
            head_3 / AP8_L,
        ),
        ("ap7-boundary", child84, bound84, AP7_L - bound84, child82, bound82 - AP7_L, witness84),
        (
            "ap7-positive",
            ap7_positive_rows,
            tuple(ap7_root_case_counts[key] for key in ("a-even", "b-even-near", "b-even-separated")),
            ap7_maximum_ratio,
            expected_ap7_maximum,
        ),
        (
            "ap7-root-split",
            near_split_row,
            near_root,
            near_child,
            near_bound,
            separated_split_row,
            separated_root,
            separated_child,
            separated_bound,
            far_split_row,
            far_root,
            far_child,
            far_bound,
            wrong_far_child,
            wrong_far_bound,
        ),
        ("ap7-adaptive", adaptive_row, adaptive_child, adaptive_bound, adaptive_witness),
        ("ap7-all-odd", odd_row, odd_child, odd_bound, odd_witness),
        ("gap-hostile", gap_hostile, hostile_child, hostile_bound, hostile_component, residual_bridges),
        ("all-root-hostile", all_root_hostile, expected_all_root, AP7_J[0]),
        ("d0-ranges", D0_CORE, D0_J, D0_L, d0_ranges),
        ("d0-four", d0_child28, d0_bound28, d0_child26, d0_bound26, d0_odd_child, d0_odd_bound),
        ("d0-six", d0_six, d0_six_levels, len(d0_six_components), d0_six_maximum, d0_six_witness),
        ("ap8-five", ap8_five, ap8_five_levels, len(ap8_five_components), ap8_five_maximum, ap8_five_witness),
    )
    semantic_sha256 = sha256(canonical(semantic_payload).encode("ascii")).hexdigest()
    require(semantic_sha256 == EXPECTED_SEMANTIC_SHA256, "semantic payload hash changed")

    print("THM-4112 ANTIPODAL COMPONENT-ANCESTRY CHAIN EXACT REFEREE")
    print("status=PASS")
    print("open_interval_semantics=strict overlap only; endpoint touching is a safe separator")
    print(
        "chain_formula=h_i=sum_{j=i}^r 2^(j-i)/(7*v_j);"
        "base_components_have_span_equal_to_h_r;reinserted_levels_are_strict"
    )
    print(
        f"symbolic_coefficient_depths=1..12;numeric_algebra_rows={algebra_numeric_rows};"
        f"one_comb_equality_components={base_component_count}"
    )
    print(
        "geometric_arbitrary_depth=lambda>=12/5;"
        "h_tail<=1/[7*v_i*(lambda-2)]<=5/(14*v_i);"
        "h_i<=lambda/[7*v_i*(lambda-2)]"
    )
    print(
        "geometric_thresholds=lambda12/5:D0_v1>=16,AP8_v1>=112;"
        "lambda3:AP8_v1>=56"
    )
    print(
        "AP7_four_boundary=(84,85,86,87);P=1/301;Q=213/25585;g84=1/98;"
        "h=39439/2149140;slack=1/61404;theta=74/609;clock=1218"
    )
    print("AP7_criterion_near_miss=(82,83,84,85);h_minus_L=151/357315;not_a_cover_claim")
    print(
        "AP7_root_split=a_even:root_a;odd_a_even_b_and_b<2a:root_b;"
        "odd_a_even_b_and_b>=2a:root_a"
    )
    print(
        "AP7_positive_bank=84<=a<b<c<d<=100,(a_or_b_even);rows=1932;"
        "root_cases=a_even:1344,b_even_near:588,b_even_separated:0;"
        "max_ratio=458228/838503;row=(88,95,97,99);component=[251/1330,137/693];"
        "length=1181/131670;bound=93167/5676440"
    )
    print(
        "AP7_root_split_controls=near(85,86,87,88)->root86,h=80839/4451790;"
        "boundary(85,170,171,172)->root85,h=146/14535;"
        "far(85,1000,1001,1002)->root85,h=185117/59559500"
    )
    print(
        "AP7_wrong_root_hostile=(85,1000,1001,1002),root1000;"
        "Q=1341/595595>g1000=3/3500;excess=83049/59559500"
    )
    print(
        "AP7_adaptive_weight7=(85,86,91,101),root86;P=2/637;Q=431/54145;"
        "g=3/301;h=81867/4656470;slack=366/465647;theta=81/707;clock=1414"
    )
    print(
        "AP7_all_odd_weight8=(47,95,97,99),root47;P=2/679;Q=477/64505;"
        "g=5/658;h=54053/3031735;slack=22847/42444290;theta=4/35;clock=70"
    )
    print(
        "parent_gap_hostile=(85,87,89,91);Q=437/54201>g85=1/238;"
        "component=[69/1274,46/637];length=23/1274;formal=11719/658155;"
        "excess=207559/838489470;root_teeth=3"
    )
    print("all_root_hostile=(43,45,47,49);every_root_has_a_two_consecutive_tooth_bridge;AP7_theta=4/35")
    print(
        "D0=(3,4,5,6,8,10,12);J=[1/42,13/168];L=3/56;gcd=1;nonAP;"
        "norm_ranges=3:[1/14,13/56],4:[2/21,13/42],5:[5/42,65/168],"
        "6:[1/7,13/28],8:[4/21,1/2],10:[19/84,1/2],12:[1/14,1/2]"
    )
    print(
        "D0_four_controls=(28,29,30,31):slack=89/170520;"
        "criterion_near_miss(26,27,28,29):deficit=457/137592;"
        "all_odd(15,31,33,35):slack=19/95480"
    )
    print(
        "D0_six_doubling=(16,32,64,128,256,512);"
        "levels=(3/56,5/224,1/112,3/896,1/896,1/3584);h1=L;"
        "components=416;max=1/112;theta=43/1792;clock=1792"
    )
    print(
        "AP8_five_doubling=(94,188,376,752,1504);"
        "levels=(5/658,1/329,3/2632,1/2632,1/10528);slack=1/18424;"
        "components=1316;max=1/658;theta=11/49;clock=98"
    )
    print("scope=conditional component ancestry plus structured positive suppliers;LRC(14) remains OPEN")
    print(f"semantic_sha256={semantic_sha256}")


if __name__ == "__main__":
    main()
