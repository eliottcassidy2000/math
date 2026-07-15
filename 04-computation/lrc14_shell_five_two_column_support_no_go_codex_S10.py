#!/usr/bin/env python3
"""Exact referee for the shell-five danger-support classification.

The relaxed shell-five lift model has

    d mod 52 in {15,37,41},  B=(13d+5)/2,
    free raw classes R={1,2,4,5,8,11,12},
    forced speeds B-3,B-2,B.

For a unit endpoint column p/q, a free raw class is in its danger support if
some lift u<=B in that class has ||pu/q||<1/11.  This script certifies:

* on q=13d every column has at least three danger classes;
* on q=5d every column not already killed by a forced speed has at least two
  danger classes, except

      d=41 mod 52, p=+/- (45d+1)/26 mod 5d,

  whose danger support is exactly {4}; and
* any two fixed unit endpoint columns can therefore be killed simultaneously
  by some choice in the relaxed lift model.

The q=5d all-size step is a finite-state proof, not an extrapolation.  A
150-speed exact skeleton confines a hypothetical singleton support to within
42/[143(B-38)] of one of six rational centres.  The endpoint grid then leaves
432 periodic affine rows.  Exact same-danger-cell inequalities certify 429;
the remaining three residue rows are precisely the explicit exceptional
family above.

Only integers and fractions.Fraction are used.  The final two-column no-go is
scoped to the relaxed signed-complement shell model.  The constructed lift set
need not satisfy the additional divisor-completeness or parity-support
conditions of THM-772 and THM-803.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, combinations_with_replacement
from math import gcd


FREE = (1, 2, 4, 5, 8, 11, 12)
OPEN_CLASSES = (15, 37, 41)
BASE_D = (171, 67, 119, 37, 89, 141, 41, 93, 145)
BASE_CENTRES = (F(2, 13), F(4, 13), F(9, 26))
OMIT_FOUR_CENTRES = (
    F(4, 39),
    F(2, 13),
    F(3, 13),
    F(4, 13),
    F(9, 26),
    F(17, 39),
)
DANGER_RADIUS = F(1, 11)


EXPECTED_POSITIVE_COMPONENTS = {
    1: (
        (F(116, 759), F(241, 1562)),
        (F(496, 1617), F(483, 1562)),
        (F(270, 781), F(80, 231)),
    ),
    2: (
        (F(27, 176), F(241, 1562)),
        (F(496, 1617), F(483, 1562)),
        (F(270, 781), F(80, 231)),
    ),
    4: (
        (F(133, 1298), F(131, 1276)),
        (F(27, 176), F(241, 1562)),
        (F(375, 1628), F(117, 506)),
        (F(485, 1584), F(483, 1562)),
        (F(270, 781), F(61, 176)),
        (F(628, 1441), F(681, 1562)),
    ),
    5: (
        (F(27, 176), F(241, 1562)),
        (F(496, 1617), F(483, 1562)),
        (F(270, 781), F(80, 231)),
    ),
    8: (
        (F(27, 176), F(241, 1562)),
        (F(496, 1617), F(483, 1562)),
        (F(270, 781), F(80, 231)),
    ),
    11: (
        (F(27, 176), F(241, 1562)),
        (F(496, 1617), F(483, 1562)),
        (F(270, 781), F(80, 231)),
    ),
    12: (
        (F(27, 176), F(63, 407)),
        (F(496, 1617), F(494, 1595)),
        (F(19, 55), F(80, 231)),
    ),
}


def shell_B(d: int) -> int:
    assert d > 0 and d % 2 == 1
    return (13 * d + 5) // 2


def balanced(value: int, modulus: int) -> int:
    residue = value % modulus
    if 2 * residue > modulus:
        residue -= modulus
    return residue


def danger(p: int, u: int, q: int) -> bool:
    return 11 * abs(balanced(p * u, q)) < q


def floor_fraction(value: F) -> int:
    return value.numerator // value.denominator


def ceil_fraction(value: F) -> int:
    return -floor_fraction(-value)


def reflect(interval: tuple[F, F]) -> tuple[F, F]:
    return (1 - interval[1], 1 - interval[0])


def interval_intersection(
    left: tuple[tuple[F, F], ...], right: tuple[tuple[F, F], ...]
) -> tuple[tuple[F, F], ...]:
    answer: list[tuple[F, F]] = []
    i = j = 0
    while i < len(left) and j < len(right):
        lo = max(left[i][0], right[j][0])
        hi = min(left[i][1], right[j][1])
        if lo <= hi:
            answer.append((lo, hi))
        if left[i][1] < right[j][1]:
            i += 1
        else:
            j += 1
    return tuple(answer)


def deep_components(speeds: tuple[int, ...]) -> tuple[tuple[F, F], ...]:
    """Closed components of {t in [0,1]: ||ut||>=1/11 for every u}."""

    answer = ((F(0), F(1)),)
    for speed in speeds:
        safe = tuple(
            (F(11 * k + 1, 11 * speed), F(11 * k + 10, 11 * speed))
            for k in range(speed)
        )
        answer = interval_intersection(answer, safe)
    return answer


def carrier_for_side(
    omitted: int,
    centre: F,
    component: tuple[F, F],
    side: int,
) -> tuple[int, int, int, F]:
    """Choose an unwrapped safe-arc carrier for x<0 or x>0.

    For x<0 we need a positive centre phase, and for x>0 a negative one.
    The endpoint check is essential: it says the first retained lift stays in
    the intended unwrapped safe arc over the whole skeleton component.
    """

    assert side in (-1, 1)
    modulus = centre.denominator
    endpoint = component[0] if side < 0 else component[1]
    candidates = []
    for state in range(1, modulus + 1):
        if state % 13 not in FREE or state % 13 == omitted:
            continue
        z = balanced(centre.numerator * state, modulus)
        phase = F(z, modulus) + state * (endpoint - centre)
        if side < 0:
            on_arc = z > 0 and DANGER_RADIUS <= phase <= F(1, 2)
        else:
            on_arc = z < 0 and -F(1, 2) <= phase <= -DANGER_RADIUS
        margin_numerator = 11 * abs(z) - modulus
        if on_arc and margin_numerator > 0:
            candidates.append((margin_numerator, state, z, phase))
    assert candidates
    return min(candidates)


def skeleton_audit() -> tuple[object, ...]:
    rows = []
    all_carriers = []
    for omitted in FREE:
        skeleton = tuple(
            u for u in range(1, 151) if u % 13 in FREE and u % 13 != omitted
        )
        components = deep_components(skeleton)
        positive = EXPECTED_POSITIVE_COMPONENTS[omitted]
        expected = positive + tuple(reflect(item) for item in reversed(positive))
        assert components == expected

        centres = OMIT_FOUR_CENTRES if omitted == 4 else BASE_CENTRES
        assert len(positive) == len(centres)
        carrier_rows = []
        for centre, component in zip(centres, positive):
            lo, hi = component
            modulus = centre.denominator
            assert lo < centre < hi < F(1, 2)
            step_bound = modulus * max(centre - lo, hi - centre)
            assert step_bound < F(2, 11)

            left = carrier_for_side(omitted, centre, component, -1)
            right = carrier_for_side(omitted, centre, component, 1)
            for carrier in (left, right):
                margin_numerator, state, _, _ = carrier
                assert margin_numerator <= 42
                assert 1 <= state <= modulus
                # The last available lift in this state is at worst m-1
                # below B; every relevant centre denominator is <=39.
                assert modulus - 1 <= 38
            carrier_rows.append((centre, component, step_bound, left, right))
            all_carriers.extend((left, right))
        rows.append((omitted, len(skeleton), tuple(carrier_rows)))

    # If all retained lifts are deep, no step of size m|x| can cross the
    # open danger gap of width 2/11.  The last carrier u is >=B-38, whence
    # |x| <= M/(11*m*u) <=42/[143(B-38)].  On q=5d this is below one grid
    # spacing for every d>=37:
    #
    #   5d*42/[143(B-38)] < 1  <=>  1439d>10153.
    assert 1439 * 37 > 10153
    assert 1439 > 0
    assert max(item[0] for item in all_carriers) == 42
    return (
        tuple(rows),
        len(all_carriers) // 2,
        max(item[0] for item in all_carriers),
        38,
    )


def support_row(d: int, p: int) -> tuple[tuple[int, ...], tuple[int, ...]]:
    q = 5 * d
    B = shell_B(d)
    support = tuple(
        residue
        for residue in FREE
        if any(danger(p, u, q) for u in range(residue, B + 1, 13))
    )
    forced = tuple(u for u in (B - 3, B - 2, B) if danger(p, u, q))
    return support, forced


def d15_audit() -> tuple[object, ...]:
    d = 15
    q = 5 * d
    rows = []
    for p in range(1, q // 2 + 1):
        if gcd(p, q) != 1:
            continue
        support, forced = support_row(d, p)
        if not forced:
            rows.append((p, support))
    histogram = Counter(len(support) for _, support in rows)
    assert len(rows) == 14
    assert histogram == Counter({7: 7, 6: 5, 5: 1, 2: 1})
    assert min(len(support) for _, support in rows) == 2
    assert next(row for row in rows if len(row[1]) == 2) == (29, (2, 5))
    return (tuple(rows), tuple(sorted(histogram.items())))


def same_danger_cell(y0: F, y_infinity: F) -> tuple[int, F] | None:
    """Return the common integer cell and its strict endpoint margin."""

    integer = floor_fraction(y0 + F(1, 2))
    margins = (
        DANGER_RADIUS - abs(y0 - integer),
        DANGER_RADIUS - abs(y_infinity - integer),
    )
    if min(margins) > 0:
        return integer, min(margins)
    return None


def affine_witness(
    d0: int, omitted: int, centre: F, p0: int
) -> tuple[object, ...] | None:
    """Find one of the two frozen affine witness types.

    Along d=d0+156h, q=q0+780h and B=B0+1014h.  If
    u=u0+alpha*h, c*alpha is integral, and both endpoint phases lie in the
    same open danger cell, monotonicity of (u0+alpha*h)/(q0+780h) proves
    danger for every h>=0.
    """

    q0 = 5 * d0
    B0 = shell_B(d0)
    delta = F(p0) - q0 * centre

    def check(
        kind: str, label: int, u0: int, alpha: int
    ) -> tuple[object, ...] | None:
        assert alpha % 13 == 0
        assert 0 <= alpha <= 1014
        assert (centre * alpha).denominator == 1
        y0 = F(p0 * u0, q0)
        y_infinity = centre * u0 + delta * F(alpha, 780)
        common = same_danger_cell(y0, y_infinity)
        if common is None:
            return None
        integer, margin = common

        # Direct replays are guards only; the same-cell monotonicity is the
        # source of the infinite quantifier.
        for h in (0, 1, 17, 1000):
            d = d0 + 156 * h
            q = 5 * d
            p = p0 + int(780 * centre) * h
            u = u0 + alpha * h
            assert 1 <= u <= shell_B(d)
            assert danger(p, u, q)
        return (kind, label, u0, alpha, integer, margin)

    # Forced speeds have fixed offsets below B.
    for offset in (3, 2, 0):
        witness = check("forced", offset, B0 - offset, 1014)
        if witness is not None:
            return witness

    # Almost every remaining row is killed by a free lift with a fixed
    # offset below B (alpha=1014).  Three rows use the slower raw-1 family
    # alpha=117.  No optimizer or unbounded search is involved.
    for alpha in (1014, 117):
        if (centre * alpha).denominator != 1:
            continue
        for u0 in range(1, B0 + 1):
            if u0 % 13 not in FREE or u0 % 13 == omitted:
                continue
            witness = check("free", u0 % 13, u0, alpha)
            if witness is not None:
                return witness
    return None


def affine_candidate_audit() -> tuple[object, ...]:
    certified = []
    unresolved = []
    for d0 in BASE_D:
        assert d0 >= 37
        centres_for_class = []
        for omitted in FREE:
            centres = OMIT_FOUR_CENTRES if omitted == 4 else BASE_CENTRES
            for centre in centres:
                grid_centre = 5 * d0 * centre
                assert grid_centre.denominator > 1
                for endpoint in (floor_fraction(grid_centre), ceil_fraction(grid_centre)):
                    witness = affine_witness(d0, omitted, centre, endpoint)
                    row = (d0, omitted, centre, endpoint, endpoint - grid_centre)
                    if witness is None:
                        unresolved.append(row)
                    else:
                        certified.append(row + (witness,))
                    centres_for_class.append((omitted, centre, endpoint))
        assert len(centres_for_class) == 48

    expected_unresolved = (
        (41, 4, F(9, 26), 71, F(1, 26)),
        (93, 4, F(9, 26), 161, F(1, 26)),
        (145, 4, F(9, 26), 251, F(1, 26)),
    )
    assert len(certified) == 429
    assert tuple(unresolved) == expected_unresolved
    type_histogram = Counter((row[-1][0], row[-1][3]) for row in certified)
    assert type_histogram == Counter(
        {("forced", 1014): 254, ("free", 1014): 172, ("free", 117): 3}
    )
    minimum_margin = min(row[-1][-1] for row in certified)
    assert minimum_margin == F(1, 6545)
    return (
        len(certified) + len(unresolved),
        len(certified),
        tuple(unresolved),
        tuple(sorted(type_histogram.items())),
        minimum_margin,
        sha256(repr(certified).encode()).hexdigest(),
    )


def exception_audit(limit: int = 100) -> tuple[object, ...]:
    """Replay and algebraically guard the unique singleton family."""

    phase_table = tuple(
        (
            residue,
            tuple(
                balanced(9 * state, 26)
                for state in (residue, residue + 13)
            ),
        )
        for residue in FREE
    )
    assert phase_table == (
        (1, (9, -4)),
        (2, (-8, 5)),
        (4, (10, -3)),
        (5, (-7, 6)),
        (8, (-6, 7)),
        (11, (-5, 8)),
        (12, (4, -9)),
    )

    # Symbolic inequalities used in the proof, checked at the least d;
    # every displayed left-minus-right numerator then has positive slope.
    d = 41
    assert 37 * d - 55 > 0
    assert 27 * d - 5 > 0  # 9/26+B/(130d)<1/2
    assert 1339 * d - 143 > 0
    assert 7 * d - 1 > 0  # 11/26+(B-2)/(130d)<1/2
    assert 87 * d - 5 > 0  # 3/26+B/(130d)<1/2
    assert 1898 * d - 8866 > 0
    assert 17 * d + 31 > 0  # (B-18)/(130d)<3/26

    rows = []
    for k in range(limit):
        d = 41 + 52 * k
        q = 5 * d
        B = shell_B(d)
        numerator = (45 * d + 1) // 26
        assert 26 * numerator - 45 * d == 1
        assert numerator % 5 == 1
        assert gcd(numerator, q) == 1
        assert B % 26 == 9

        support, forced = support_row(d, numerator)
        reflected_support, reflected_forced = support_row(d, q - numerator)
        assert support == reflected_support == (4,)
        assert forced == reflected_forced == ()

        raw_four_witness = B - 18
        assert raw_four_witness % 26 == 17
        assert danger(numerator, raw_four_witness, q)
        rows.append((d, B, numerator, raw_four_witness))
    return (limit, phase_table, rows[0], rows[-1], sha256(repr(rows).encode()).hexdigest())


def q13_three_pair_witnesses(d: int, p: int) -> tuple[tuple[int, int, int], ...]:
    q = 13 * d
    assert gcd(p, q) == 1
    inverse = pow(p, -1, q)
    rows = []
    supports = []
    for signed_pair_representative in (1, 2, 5):
        j = balanced(p * signed_pair_representative, 13)
        v = balanced(inverse * j, q)
        u = abs(v)
        assert v % 13 == signed_pair_representative
        assert 1 <= abs(j) <= 6
        assert 1 <= u <= q // 2 < shell_B(d)
        assert u % 13 in (
            signed_pair_representative,
            -signed_pair_representative % 13,
        )
        assert danger(p, u, q)
        rows.append((signed_pair_representative, j, u))
        supports.append(u % 13)
    assert len(set(supports)) == 3
    assert all(residue in FREE for residue in supports)
    return tuple(rows)


def q13_audit(class_depth: int = 6) -> tuple[object, ...]:
    hasher = sha256()
    unit_columns = 0
    for residue_class in OPEN_CLASSES:
        for k in range(class_depth):
            d = residue_class + 52 * k
            q = 13 * d
            assert 6 * 11 < q
            for p in range(1, q // 2 + 1):
                if gcd(p, q) != 1:
                    continue
                witnesses = q13_three_pair_witnesses(d, p)
                hasher.update(repr((d, p, witnesses)).encode())
                unit_columns += 1
    return (class_depth, unit_columns, hasher.hexdigest())


def danger_lifts(d: int, p: int, q: int, residue: int) -> tuple[int, ...]:
    B = shell_B(d)
    return tuple(u for u in range(residue, B + 1, 13) if danger(p, u, q))


def endpoint_column(d: int, q: int, p: int) -> tuple[object, ...]:
    B = shell_B(d)
    if q == 5 * d:
        support, forced_speeds = support_row(d, p)
    else:
        assert q == 13 * d
        support = tuple(
            residue
            for residue in FREE
            if danger_lifts(d, p, q, residue)
        )
        forced_speeds = tuple(
            u for u in (B - 3, B - 2, B) if danger(p, u, q)
        )
        assert len(support) >= 3
    return (q, p, support, forced_speeds)


def cover_two_columns(
    d: int, first: tuple[object, ...], second: tuple[object, ...]
) -> tuple[int, ...]:
    """Construct ten relaxed shell lifts which kill both columns."""

    B = shell_B(d)
    selected: dict[int, int] = {}
    pending = [column for column in (first, second) if not column[3]]

    if len(pending) == 1:
        q, p, support, _ = pending[0]
        residue = support[0]
        selected[residue] = danger_lifts(d, p, q, residue)[0]
    elif len(pending) == 2:
        q1, p1, support1, _ = pending[0]
        q2, p2, support2, _ = pending[1]
        if len(support1) == len(support2) == 1:
            assert support1 == support2 == (4,)
            common = tuple(
                sorted(
                    set(danger_lifts(d, p1, q1, 4))
                    & set(danger_lifts(d, p2, q2, 4))
                )
            )
            assert common and B - 18 in common
            selected[4] = B - 18
        else:
            # One support has at least two classes.  Put it first and choose
            # an SDR for the two column obligations.
            if len(support1) < 2:
                q1, p1, support1, q2, p2, support2 = (
                    q2,
                    p2,
                    support2,
                    q1,
                    p1,
                    support1,
                )
            residue2 = support2[0]
            residue1 = next(residue for residue in support1 if residue != residue2)
            selected[residue1] = danger_lifts(d, p1, q1, residue1)[0]
            selected[residue2] = danger_lifts(d, p2, q2, residue2)[0]

    free_lifts = tuple(selected.get(residue, residue) for residue in FREE)
    U = tuple(sorted(free_lifts + (B - 3, B - 2, B)))
    assert len(U) == 10 and len(set(U)) == 10
    assert {u % 13 for u in free_lifts} == set(FREE)
    assert all(1 <= u <= B for u in U)
    for q, p, _, _ in (first, second):
        assert any(danger(p, u, q) for u in U)
    return U


def finite_two_column_replay() -> tuple[object, ...]:
    rows = []
    hasher = sha256()
    for d in OPEN_CLASSES:
        columns = []
        for q in (5 * d, 13 * d):
            for p in range(1, q // 2 + 1):
                if gcd(p, q) == 1:
                    columns.append(endpoint_column(d, q, p))
        pairs = 0
        for first, second in combinations_with_replacement(columns, 2):
            U = cover_two_columns(d, first, second)
            hasher.update(repr((d, first[:2], second[:2], U)).encode())
            pairs += 1
        rows.append((d, len(columns), pairs))
    return (tuple(rows), hasher.hexdigest())


def tournament_from_weights(
    weights: dict[int, int], tie_path: tuple[int, ...]
) -> tuple[object, ...]:
    position = {vertex: index for index, vertex in enumerate(tie_path)}
    edge = {}
    for left, right in combinations(tie_path, 2):
        if weights[left] != weights[right]:
            winner = left if weights[left] > weights[right] else right
        else:
            winner = left if position[left] < position[right] else right
        edge[(left, right)] = winner

    def beats(left: int, right: int) -> bool:
        if left == right:
            return False
        key = (left, right) if left < right else (right, left)
        return edge[key] == left

    scores = {vertex: sum(beats(vertex, other) for other in tie_path if other != vertex) for vertex in tie_path}
    cycles = sum(
        beats(a, b) and beats(b, c) and beats(c, a)
        or beats(a, c) and beats(c, b) and beats(b, a)
        for a, b, c in combinations(tie_path, 3)
    )

    reach = {vertex: {other for other in tie_path if beats(vertex, other)} | {vertex} for vertex in tie_path}
    for pivot in tie_path:
        for vertex in tie_path:
            if pivot in reach[vertex]:
                reach[vertex] |= reach[pivot]
    sccs = []
    unseen = set(tie_path)
    while unseen:
        vertex = min(unseen, key=position.get)
        component = tuple(other for other in tie_path if other in reach[vertex] and vertex in reach[other])
        sccs.append(component)
        unseen -= set(component)

    # Exact Hamiltonian-path count by subset DP.
    n = len(tie_path)
    dp = [[0] * n for _ in range(1 << n)]
    for index in range(n):
        dp[1 << index][index] = 1
    for mask in range(1 << n):
        for last in range(n):
            if not dp[mask][last]:
                continue
            for nxt in range(n):
                if mask & (1 << nxt) or not beats(tie_path[last], tie_path[nxt]):
                    continue
                dp[mask | (1 << nxt)][nxt] += dp[mask][last]
    paths = sum(dp[-1])
    return (
        tuple(weights[vertex] for vertex in tie_path),
        tuple(sorted(scores.values())),
        cycles,
        tuple(sccs),
        paths,
        edge,
    )


def tournament_telemetry() -> tuple[object, ...]:
    """A lossy raw-class quotient of the exact danger incidence object."""

    d = 41
    weights = {}
    for q in (5 * d, 13 * d):
        columns = [
            endpoint_column(d, q, p)
            for p in range(1, q // 2 + 1)
            if gcd(p, q) == 1
        ]
        weights[q] = {
            residue: sum(not column[3] and residue in column[2] for column in columns)
            for residue in FREE
        }
    q5 = tournament_from_weights(weights[5 * d], FREE)
    q13 = tournament_from_weights(weights[13 * d], FREE)
    flips = sum(q5[-1][edge] != q13[-1][edge] for edge in q5[-1])
    return (
        d,
        tuple(weights[5 * d].items()),
        tuple(weights[13 * d].items()),
        q5[:-1],
        q13[:-1],
        flips,
    )


def main() -> None:
    skeleton = skeleton_audit()
    base_case = d15_audit()
    affine = affine_candidate_audit()
    exception = exception_audit()
    q13 = q13_audit()
    pairs = finite_two_column_replay()
    tournament = tournament_telemetry()
    certificate = sha256(
        repr((skeleton, base_case, affine, exception, q13, pairs, tournament)).encode()
    ).hexdigest()

    print("LRC14 SHELL-FIVE TWO-COLUMN SUPPORT NO-GO")
    print("arithmetic=integer+fractions.Fraction optimizer=none SAT=none floating_point=none")
    print("model=relaxed_shell_admissible_lifts_only")
    print("scope_guard=does_not_impose_THM-772_divisor_or_THM-803_parity_conditions")
    print()
    print("LEAVE_ONE_CLASS_SKELETON")
    print("speeds=1..150_in_retained_free_raw_classes")
    print(f"omission_rows={len(skeleton[0])} positive_components={skeleton[1]}")
    print("centres_non4=(2/13,4/13,9/26)")
    print("centres_omit4=(4/39,2/13,3/13,4/13,9/26,17/39)")
    print(f"max_wall_margin_M={skeleton[2]} max_top_lift_deficit={skeleton[3]}")
    print("radius_bound=42/[143(B-38)] q5_radius_below_one_for_d>=37")
    print(f"skeleton_sha256={sha256(repr(skeleton[0]).encode()).hexdigest()}")
    print()
    print("D15_BASE_CASE")
    print(f"forced_uncovered_rows={len(base_case[0])} support_histogram={base_case[1]}")
    print("minimum_row=(p=29,support={2,5})")
    print(f"rows_sha256={sha256(repr(base_case[0]).encode()).hexdigest()}")
    print()
    print("AFFINE_ALL_SIZE_CERTIFICATE")
    print("period=d0+156h q=q0+780h B=B0+1014h")
    print(f"base_d={BASE_D} candidate_rows={affine[0]} certified_rows={affine[1]}")
    print(f"witness_type_histogram={affine[3]} minimum_strict_margin={affine[4]}")
    print(f"unresolved_rows={affine[2]}")
    print("unresolved_union=d_mod52_41,p=(45d+1)/26,omitted_raw_class=4")
    print(f"certified_rows_sha256={affine[5]}")
    print()
    print("EXCEPTION_FAMILY")
    print("classification=d_mod52_41,p=+/-(45d+1)/26,support={4},forced_uncovered")
    print(f"phase_table={exception[1]}")
    print(f"replay_limit={exception[0]} first={exception[2]} last={exception[3]}")
    print(f"rows_sha256={exception[4]}")
    print()
    print("Q13_THREE_SIGNED_PAIRS")
    print("pairs={1,12},{2,11},{5,8};inverse_lift_gives_one_danger_class_in_each")
    print(f"class_depth={q13[0]} replayed_positive_half_unit_columns={q13[1]} rows_sha256={q13[2]}")
    print()
    print("TWO_COLUMN_REPLAY")
    print(f"rows={pairs[0]} rows_sha256={pairs[1]}")
    print("columns=positive_half_representatives reflection_supplies_full_unit_grids")
    print("all_size_carrier=two-support_SDR_plus_common_raw4_lift_for_double_exception")
    print()
    print("TOURNAMENT_ANALYSIS")
    print("vertices=free_raw_classes_not_runners")
    print("pairwise_observable=number_of_nonforced_endpoint_columns_supported_at_d41")
    print("switch=q5d_coverage_gauge_to_q13d_coverage_gauge tie_Hamiltonian_path=(1,2,4,5,8,11,12)")
    print(f"q5_weights={tournament[1]} q13_weights={tournament[2]}")
    print(f"q5_fingerprint={tournament[3]} q13_fingerprint={tournament[4]} edge_flips={tournament[5]}")
    print("preserves=marginal_raw_class_coverage_order")
    print("destroys=support_intersections_and_common_concrete_lifts")
    print("theorem_carrier=bipartite_column_to_raw_class_incidence_with_lift_labels")
    print()
    print("VERDICT")
    print("q5_nonforced_support_size>=2_except_explicit_reflection_pair_with_support={4}")
    print("q13_support_size>=3")
    print("any_two_fixed_unit_endpoint_columns_are_jointly_killable_in_the_relaxed_shell_model")
    print(f"certificate_sha256={certificate}")
    print("PASS")


if __name__ == "__main__":
    main()
