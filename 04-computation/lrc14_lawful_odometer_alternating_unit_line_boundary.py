#!/usr/bin/env python3
"""Exact THM-2640 unit-line audit at the lawful odometer two-cycle.

The alternating THM-2584 rails through

    x_+=1/2+(13^5+1)/(14*13^6),
    x_-=1/2-(13^5+1)/(14*13^6)

have metadata ``(1,3,0)`` and ``(1,4,12)``.  Both cycle points have future
coordinate ``{13^6*x}=1/2``, hence ``h=6`` and (with the canonical half-open
choice) ``kappa=1``; their predecessor carries are 7 and 5.

This companion rebuilds just those two THM-2640 rows.  The same-edge unit
profiles have exact clock-reversal scalar laws, but those coefficient laws
do not transport the underlying physical packet near the cycle.  On the
positive side of the future-half wall, the coefficient-covariant edge 0
halves begin one deep subcell beyond the cycle; local geometry instead uses
edge 1.  On the negative side, the analogous edge-1 lines miss and local
geometry uses edge 0.  The locally geometric alternatives are still units,
but lose the scalar clock-reversal law.

More decisively, an explicit open neighbourhood of ``x_+`` is outside the
union of *all seven* THM-2640 present factors at ``h=6``.  The alternating
THM-2584 H=2 cylinder still meets that fixed ``h=6`` union, whereas the whole
H=3 cylinder lies strictly in the gap.  Dynamically retyping
``h(x)=floor(13*{13^6*x})`` leaves smaller escape pieces through H=5, but an
independent exact digit-cell intersection shows that H=6 is the sharp first
horizon with no current-present support for *any* ``h``.  (The reflected point
``x_-`` is present, so this is deliberately a plus-current obstruction, not a
symmetric claim.)  The targeted unit coefficients are integrated witnesses
elsewhere on their rows, not pointwise full packets transported by the
repelling rail cycle.
"""

from bisect import bisect_right
from fractions import Fraction
from math import floor

import lrc14_predecessor_carry_private_root_atlas_thm2640 as m
import lrc14_slope7_twelve_chart_component_witness_thm2672 as witness


P = 13
R = P**6
S = P**5
K = S + 1
ROW_BOUNDS = (7, 9)
METADATA = ((1, 3, 0), (1, 4, 12))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def normalize(values, root):
    scale = pow(root, -1, P)
    return tuple((value // 26) * scale % P for value in values)


def reverse_clock(values):
    return tuple(values[(-ell) % 7] for ell in range(7))


def private_interval(root, edge):
    if edge == 0:
        return Fraction(14 * root - 13, 182), Fraction(14 * root, 182)
    return Fraction(14 * root, 182), Fraction(14 * root + 13, 182)


def in_open(value, interval):
    return interval[0] < value < interval[1]


def merge_intervals(intervals):
    out = []
    for left, right in sorted(intervals):
        if out and left <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], right))
        else:
            out.append((left, right))
    return tuple(out)


def complement_component(intervals, point, modulus):
    require(not any(left <= point < right for left, right in intervals),
            "point lies in the present-factor union")
    before = max((right for left, right in intervals if right <= point),
                 default=0)
    after = min((left for left, right in intervals if left >= point),
                default=modulus)
    require(before < point < after,
            "point entered the present-factor union")
    return before, after


def intersect_open_with_grid_union(interval, intervals, modulus):
    """Positive-length pieces of a rational interval in an integer-grid union."""
    left, right = interval
    pieces = []
    for grid_left, grid_right in intervals:
        overlap_left = max(left, Fraction(grid_left, modulus))
        overlap_right = min(right, Fraction(grid_right, modulus))
        if overlap_left < overlap_right:
            pieces.append((overlap_left, overlap_right))
    return tuple(pieces)


def intersect_sorted_unions(first, second):
    out = []
    i = j = 0
    while i < len(first) and j < len(second):
        left = max(first[i][0], second[j][0])
        right = min(first[i][1], second[j][1])
        if left < right:
            out.append((left, right))
        if first[i][1] <= second[j][1]:
            i += 1
        else:
            j += 1
    return merge_intervals(out)


def alternating_translation_coefficients(horizon):
    """Coefficients c_i in x_i={P^i*x+c_i/R}, starting at x_+."""
    coefficients = [0]
    lifts = (-K, K)
    for step in range(horizon - 1):
        coefficients.append(P * coefficients[-1] + lifts[step % 2])
    return tuple(coefficients)


def alternating_rail_cylinder(plus, outer_margin, inner_margin, horizon):
    scale = P ** (horizon - 1)
    if horizon % 2:
        return (
            plus - outer_margin / scale,
            plus + inner_margin / scale,
        )
    return (
        plus - inner_margin / scale,
        plus + outer_margin / scale,
    )


def dynamic_present_intersection(cylinder, horizon, present_by_h, base_grid):
    """Intersect after exact local typing by ``floor(P*{R*x})``.

    On a common integer grid, digit ``h`` occupies exactly the cells
    ``[(P*j+h)u,(P*j+h+1)u)``, where ``u=grid/(P*R)``.  Only cells meeting the
    cylinder are scanned; each is intersected with the present union for its
    own digit.  In particular, this never freezes the centre's boundary value
    ``h=6`` across a neighbourhood.
    """
    exponent = max(1, horizon - 1)
    grid = P**exponent * base_grid
    require(grid % (P * R) == 0,
            "dynamic future-digit cells left the common integer grid")
    endpoints = tuple(endpoint * grid for endpoint in cylinder)
    require(all(endpoint.denominator == 1 for endpoint in endpoints),
            "rail cylinder left the asserted common integer grid")
    cylinder_left, cylinder_right = (
        endpoint.numerator for endpoint in endpoints
    )
    unit = grid // (P * R)
    scale_present = grid // base_grid
    scaled_present = tuple(tuple(
        (left * scale_present, right * scale_present)
        for left, right in present_by_h[h]
    ) for h in range(P))
    present_ends = tuple(tuple(right for _, right in intervals)
                         for intervals in scaled_present)

    first_digit = max(0, cylinder_left // unit)
    last_digit = min(P * R - 1, (cylinder_right - 1) // unit)
    tagged_pieces = []
    for digit_index in range(first_digit, last_digit + 1):
        h = digit_index % P
        digit_left = max(cylinder_left, digit_index * unit)
        digit_right = min(cylinder_right, (digit_index + 1) * unit)
        intervals = scaled_present[h]
        index = bisect_right(present_ends[h], digit_left)
        while index < len(intervals) and intervals[index][0] < digit_right:
            left = max(digit_left, intervals[index][0])
            right = min(digit_right, intervals[index][1])
            if left < right:
                tagged_pieces.append((left, right, h))
            index += 1

    merged = merge_intervals((left, right)
                             for left, right, _ in tagged_pieces)
    mass_by_h = tuple(
        Fraction(sum(right - left for left, right, digit
                     in tagged_pieces if digit == h), grid)
        for h in range(P)
    )
    return (
        tuple((Fraction(left, grid), Fraction(right, grid))
              for left, right in merged),
        sum(mass_by_h, Fraction(0)),
        mass_by_h,
    )


def main():
    result = m.shard(ROW_BOUNDS)
    content, metadata, rows = result[1], result[5], result[6]
    require(content == 26 and metadata == METADATA and len(rows) == 2,
            "targeted THM-2640 row reconstruction changed")

    # (name,row,carry,kappa,edge,root,expected normalized profile)
    covariant_specs = (
        ("plus_positive", 0, 7, 1, 0, 3, (0, 0, 0, 0, 8, 8, 7)),
        ("minus_positive", 1, 5, 1, 0, 12, (0, 5, 2, 2, 0, 0, 0)),
        ("plus_negative", 0, 7, 0, 1, 1, (0, 0, 0, 0, 11, 11, 8)),
        ("minus_negative", 1, 5, 0, 1, 10, (0, 6, 5, 5, 0, 0, 0)),
    )
    covariant_profiles = {}
    for name, row, carry, kappa, edge, root, expected in covariant_specs:
        computed_root = (
            2 * carry + (12 + kappa) // P + (edge == 0)
        ) % P
        require(computed_root == root,
                "targeted private root typing changed")
        values = rows[row][0][edge][carry][kappa][6]
        profile = normalize(values, root)
        require(profile == expected and m.is_unit(values, root, 26),
                "targeted normalized unit profile changed")
        covariant_profiles[name] = profile
    require(
        covariant_profiles["minus_positive"]
        == tuple(10 * value % P for value in reverse_clock(
            covariant_profiles["plus_positive"]
        )), "positive-side clock-reversal scalar law changed",
    )
    require(
        covariant_profiles["minus_negative"]
        == tuple(4 * value % P for value in reverse_clock(
            covariant_profiles["plus_negative"]
        )), "negative-side clock-reversal scalar law changed",
    )

    # The edge which actually contains a sufficiently near one-sided cycle
    # point is opposite to the coefficient-covariant edge above.
    local_specs = (
        ("plus_positive", 0, 7, 1, 1, 2, (0, 0, 0, 0, 9, 5, 4)),
        ("minus_positive", 1, 5, 1, 1, 11, (0, 3, 9, 9, 0, 0, 0)),
        ("plus_negative", 0, 7, 0, 0, 2, (0, 0, 0, 0, 4, 4, 10)),
        ("minus_negative", 1, 5, 0, 0, 11, (0, 9, 8, 4, 0, 0, 0)),
    )
    local_profiles = {}
    for name, row, carry, kappa, edge, root, expected in local_specs:
        values = rows[row][0][edge][carry][kappa][6]
        profile = normalize(values, root)
        require(profile == expected and m.is_unit(values, root, 26),
                "locally geometric unit profile changed")
        local_profiles[name] = profile
    for suffix in ("positive", "negative"):
        plus_profile = local_profiles[f"plus_{suffix}"]
        minus_profile = local_profiles[f"minus_{suffix}"]
        require(not any(
            minus_profile
            == tuple(scalar * value % P
                     for value in reverse_clock(plus_profile))
            for scalar in range(1, P)
        ), "local geometric profiles unexpectedly gained scalar covariance")

    amplitude = Fraction(K, 14 * R)
    plus = Fraction(1, 2) + amplitude
    minus = Fraction(1, 2) - amplitude
    require((R * plus) % 1 == (R * minus) % 1 == Fraction(1, 2),
            "cycle future coordinate left the half wall")
    require((floor(P * ((S * plus) % 1)),
             floor(P * ((S * minus) % 1))) == (7, 5),
            "cycle predecessor carries changed")

    # A radius valid for the three-event scout keeps the same sign under all
    # expanding iterates.  Check the claimed and local private halves at both
    # cycle states on each one-sided branch.
    horizon = 3
    radius = Fraction(1, 3 * R * P**horizon)
    geometric_rows = []
    for sign_name, sign, kappa, covariant_edge, local_edge in (
        ("positive", 1, 1, 0, 1),
        ("negative", -1, 0, 1, 0),
    ):
        for state_name, centre, row, carry in (
            ("plus", plus, 0, 7), ("minus", minus, 1, 5)
        ):
            point = centre + sign * radius / 2
            y = (R * point) % 1
            digit = floor(2 * P * y)
            h_seen, kappa_seen = divmod(digit, 2)
            require((h_seen, kappa_seen) == (6, kappa),
                    "one-sided point left the targeted future half-digit")
            deep = (2 * S * point) % 1
            covariant_root = (
                2 * carry + (12 + kappa) // P
                + (covariant_edge == 0)
            ) % P
            local_root = (
                2 * carry + (12 + kappa) // P
                + (local_edge == 0)
            ) % P
            covariant_ok = in_open(
                deep, private_interval(covariant_root, covariant_edge)
            )
            local_ok = in_open(deep, private_interval(local_root, local_edge))
            require(not covariant_ok and local_ok,
                    "one-sided private-edge geometry changed")
            require(m.is_unit(
                rows[row][0][local_edge][carry][kappa][6],
                local_root, 26,
            ), "local private edge lost its unit")
            geometric_rows.append((
                state_name, sign_name, carry, kappa,
                covariant_edge, covariant_root, covariant_ok,
                local_edge, local_root, local_ok,
            ))

    module, _, _, _, rails, present, _ = m.core.build_carrier_data()
    require(tuple(rails[index][:3] for index in range(*ROW_BOUNDS))
            == METADATA,
            "targeted row indices stopped matching their THM-2584 rails")
    all_present_h6 = merge_intervals(
        interval
        for ell5 in range(m.Q7)
        for interval in present[ell5, (-6) % P]
    )
    plus_grid = plus * m.T
    minus_grid = minus * m.T
    require(plus_grid.denominator == minus_grid.denominator == 1,
            "cycle points left the exact present-factor grid")
    plus_gap_grid = complement_component(
        all_present_h6, plus_grid.numerator, m.T
    )
    plus_gap = tuple(Fraction(endpoint, m.T)
                     for endpoint in plus_gap_grid)
    require(plus_gap == (
        Fraction(202105, 399854), Fraction(202131, 399854)
    ), "plus-cycle present-free gap changed")
    minus_present_clocks = tuple(
        ell5 for ell5 in range(m.Q7)
        if witness.in_interval_union(
            minus * m.T, present[ell5, (-6) % P]
        )
    )
    require(minus_present_clocks == (3,),
            "reflected cycle point's present-clock witness changed")
    largest_standard_radius = Fraction(1, 3 * R * P)
    require(plus_gap[0] < plus - largest_standard_radius
            < plus + largest_standard_radius < plus_gap[1],
            "standard finite-horizon neighbourhood escaped the present gap")
    require(not any(
        witness.in_interval_union(plus * m.T,
                                  present[ell5, (-6) % P])
        for ell5 in range(m.Q7)
    ), "cycle point unexpectedly gained a present clock")

    # Use the exact alternating THM-2584 rail margins, independently of the
    # companion rail-horizon program.  H=2 is a hostile positive control;
    # H=3 is the first cylinder wholly swallowed by the plus present-free gap.
    outer_margin = Fraction(14_281, 7 * R)
    inner_margin = Fraction(2_040, R)
    horizon_intersections = []
    for horizon in (1, 2, 3):
        cylinder = alternating_rail_cylinder(
            plus, outer_margin, inner_margin, horizon
        )
        pieces = intersect_open_with_grid_union(
            cylinder, all_present_h6, m.T
        )
        horizon_intersections.append((
            horizon, cylinder, len(pieces),
            sum((right - left for left, right in pieces), Fraction(0)),
        ))
    require(horizon_intersections == [
        (1, (Fraction(1195, 2366), Fraction(171, 338)),
         569, Fraction(3411, 5198102)),
        (2, (Fraction(63433983, 125497034),
             Fraction(444095003, 878479238)),
         12, Fraction(11337, 878479238)),
        (3, (Fraction(5772835171, 11420230094),
             Fraction(824698899, 1631461442)),
         0, Fraction(0)),
    ], "exact current-present horizon census changed")
    h3 = horizon_intersections[2][1]
    require(plus_gap[0] < h3[0] < h3[1] < plus_gap[1]
            and h3[1] - h3[0] == Fraction(1, 7 * P**4),
            "three-event rail cylinder escaped the strict present-free gap")

    # Now restore the pointwise future digit instead of freezing h=6 at the
    # boundary centre.  This is the physically typed current-present test.
    all_present_by_h = tuple(merge_intervals(
        interval
        for ell5 in range(m.Q7)
        for interval in present[ell5, (-h) % P]
    ) for h in range(P))
    dynamic_rows = []
    dynamic_piece_rows = []
    for horizon in range(1, 7):
        cylinder = alternating_rail_cylinder(
            plus, outer_margin, inner_margin, horizon
        )
        pieces, mass, mass_by_h = dynamic_present_intersection(
            cylinder, horizon, all_present_by_h, m.T
        )
        dynamic_piece_rows.append(pieces)
        dynamic_rows.append((
            horizon, len(pieces), mass,
            tuple((h, value) for h, value in enumerate(mass_by_h) if value),
        ))
    require(tuple(row[2] for row in dynamic_rows) == (
        Fraction(7_021_099, 12_298_709_332),
        Fraction(76_227, 1_756_958_476),
        Fraction(5_807, 1_756_958_476),
        Fraction(405, 1_756_958_476),
        Fraction(1_883, 275_716_983_698),
        Fraction(0),
    ), "dynamically typed present-support masses changed")
    require(dynamic_rows[4][1] == 1
            and dynamic_rows[4][3] == (
                (7, Fraction(1_883, 275_716_983_698)),
            ) and dynamic_rows[5][1:] == (0, Fraction(0), ()),
            "sharp H5-to-H6 dynamic present boundary changed")
    h6_cylinder = alternating_rail_cylinder(
        plus, outer_margin, inner_margin, 6
    )
    h7_cylinder = alternating_rail_cylinder(
        plus, outer_margin, inner_margin, 7
    )
    require(P * inner_margin > outer_margin
            and P * outer_margin > inner_margin
            and h6_cylinder[0] < h7_cylinder[0]
            < h7_cylinder[1] < h6_cylinder[1],
            "post-H6 rail cylinders stopped nesting")
    h5_survivor = dynamic_piece_rows[4][0]
    require(h5_survivor[1] - h5_survivor[0]
            == Fraction(1_883, 275_716_983_698),
            "last dynamic present survivor length changed")
    h1_dynamic_pieces = dynamic_piece_rows[0]
    require(not any(left <= plus < right
                    for left, right in h1_dynamic_pieces),
            "plus cycle point entered dynamically typed present support")
    dynamic_gap = (
        max(right for left, right in h1_dynamic_pieces if right <= plus),
        min(left for left, right in h1_dynamic_pieces if left >= plus),
    )
    require(dynamic_gap == (
        Fraction(31_719_030, 62_748_517),
        Fraction(31_719_032, 62_748_517),
    ) and (plus - dynamic_gap[0], dynamic_gap[1] - plus) == (
        Fraction(3, 26 * R), Fraction(1, 26 * R)
    ), "two-future-digit dynamic present-free gap changed")
    require(h5_survivor[0] == dynamic_gap[1]
            and dynamic_gap[0] < h6_cylinder[0]
            < h6_cylinder[1] < dynamic_gap[1],
            "H5/H6 boundary stopped being controlled by the h=7 endpoint")

    # Stronger physical control: require dynamically typed present support at
    # every state of the alternating affine orbit, not only at state zero.
    # The state-zero support above binds at H=6, but the positive H<=5 rows
    # verify that earlier event requirements do not manufacture an earlier
    # empty intersection.
    joint_rows = []
    joint_piece_rows = []
    for horizon in range(1, 7):
        cylinder = alternating_rail_cylinder(
            plus, outer_margin, inner_margin, horizon
        )
        joint = dynamic_piece_rows[horizon - 1]
        coefficients = alternating_translation_coefficients(horizon)
        for step in range(1, horizon):
            coefficient = coefficients[step]
            raw_centre = (P**step * plus
                          + Fraction(coefficient, R))
            branch = floor(raw_centre)
            state_centre = raw_centre - branch
            require(state_centre == (plus if step % 2 == 0 else minus),
                    "alternating affine state centre changed")
            state_interval = (
                P**step * cylinder[0] + Fraction(coefficient, R) - branch,
                P**step * cylinder[1] + Fraction(coefficient, R) - branch,
            )
            require(Fraction(0) < state_interval[0]
                    < state_interval[1] < Fraction(1),
                    "local affine state interval crossed a circle branch")
            state_support, _, _ = dynamic_present_intersection(
                state_interval, horizon - step, all_present_by_h, m.T
            )
            pulled_support = tuple((
                (left + branch - Fraction(coefficient, R)) / P**step,
                (right + branch - Fraction(coefficient, R)) / P**step,
            ) for left, right in state_support)
            joint = intersect_sorted_unions(joint, pulled_support)
            if not joint:
                break
        joint_rows.append((
            horizon, len(joint),
            sum((right - left for left, right in joint), Fraction(0)),
        ))
        joint_piece_rows.append(joint)
    require(alternating_translation_coefficients(6) == (
        0, -371_294, -4_455_528, -58_293_158,
        -757_439_760, -9_847_088_174,
    ), "alternating affine translation coefficients changed")
    require(tuple(joint_rows) == (
        (1, 6380, Fraction(7_021_099, 12_298_709_332)),
        (2, 4298, Fraction(687_145, 22_840_460_188)),
        (3, 3207, Fraction(780_376, 519_620_469_277)),
        (4, 1645, Fraction(67_780, 965_009_442_943)),
        (5, 68, Fraction(1_259, 7_168_641_576_148)),
        (6, 0, Fraction(0)),
    ), "full-event dynamically typed present census changed")
    h5_joint_witness = (
        Fraction(905_927_272_952, 1_792_160_394_037),
        Fraction(905_927_272_953, 1_792_160_394_037),
    )
    require(h5_joint_witness in joint_piece_rows[4]
            and h5_joint_witness[1] - h5_joint_witness[0]
            == Fraction(1, P**11),
            "explicit open five-event joint witness changed")
    h5_midpoint = sum(h5_joint_witness, Fraction(0)) / 2
    h5_coefficients = alternating_translation_coefficients(5)
    h5_digits = tuple(floor(P * ((R * (
        P**step * h5_midpoint + Fraction(coefficient, R)
    )) % 1)) for step, coefficient in enumerate(h5_coefficients))
    require(h5_digits == (7, 0, 0, 0, 0),
            "five-event joint witness digit word changed")

    print("LRC14 lawful odometer alternating THM2640 unit-line boundary")
    print(f"target_rows={ROW_BOUNDS} metadata={metadata} global_content={content}")
    print(f"cycle=(plus:{plus},minus:{minus}) carries=(7,5) future=(h:6,wall:1/2)")
    print(f"coefficient_covariant_profiles={tuple(covariant_profiles.items())}")
    print("coefficient_covariance=positive side Yminus(ell)=10*Yplus(-ell); negative side scalar=4; all four are units")
    print(f"one_sided_private_geometry={tuple(geometric_rows)}")
    print(f"local_geometric_unit_profiles={tuple(local_profiles.items())}")
    print("local_profile_covariance=NONE: no nonzero scalar times reversed plus profile equals minus")
    print(f"all_clock_present_free_gap_plus={plus_gap}")
    print(f"minus_cycle_present_clocks={minus_present_clocks} symmetric_gap_claim=False")
    print(f"largest_standard_radius_Hge1={largest_standard_radius} contained_in_plus_gap=True")
    print(f"plus_current_present_horizon_census={tuple(horizon_intersections)}")
    print("fixed_h6_boundary=H2 has positive mass; H3 has zero mass and lies strictly inside the plus gap")
    print(f"dynamic_present_horizon_census={tuple(dynamic_rows)}")
    print(f"last_dynamic_survivor=(H:5,h:7,interval:{h5_survivor})")
    print(f"first_empty_dynamic_cylinder=(H:6,interval:{h6_cylinder})")
    print(f"dynamic_present_free_gap={dynamic_gap} relative_margins=(3/(26R),1/(26R)) missing_future_digits=(5,6)")
    print(f"all_event_dynamic_present_census={tuple(joint_rows)}")
    print(f"open_H5_all_event_witness={h5_joint_witness} length=1/13^11 digit_word={h5_digits}")
    print("all_event_boundary=earlier states split the H5 survivor into 68 components but do not kill it; state zero binds at H6")
    print("sharp_dynamic_boundary=H5 survives only in one h=7 component; H6 is empty; nesting forces all H>=6 empty")
    print("verdict=coefficient-line covariance is exact, but it is not physical packet transport")
    print("first_failures=claimed covariant edge misses the near-cycle private half; repaired local edge is unit but loses scalar covariance; fixed h=6 dies at H3 and dynamically typed current-present support dies at H6")
    print("scope=selected alternating rail cycle and all current-present event factors; delayed words,private geometry over surviving components,configuration switching,semantic transition,row exclusion,and LRC14 remain untested")


if __name__ == "__main__":
    main()
