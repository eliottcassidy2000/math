#!/usr/bin/env python3
"""Exact locally packet-typed transverse-17 two-cycle for the LRC14 carrier.

The two affine maps have the old base-thirteen linear part but translations
on the enlarged ``17*13^6`` grid.  Their delayed phases are 12/17 and 5/17,
which exchange the strict THM-2693 raw-word points 4/17 and 13/17.  Suitable
13-primary fibre components place the lifted circle points inside the actual
THM-2584/2640 rail, present, sector, carry, root, and primitive-unit packets.

The script certifies a positive packet cylinder at every prescribed finite
horizon.  It does not type the new transverse phase as a THM-2657 root lift
and supplies no semantic endpoint, row exclusion, or LRC(14) conclusion.
"""

from bisect import bisect_right
from fractions import Fraction
from math import gcd
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc14_odometer_alternating_lift_labelled_tail_scout_20260728 as old


P = 13
R = P**6
L = 17 * R

# Integer numerators on the enlarged C_(17R) stalk.
POINT_NUMERATORS = (39_123_022, 41_305_508)
LIFTS = (25_040_740, 76_541_689)
TRANSVERSE_PHASES = tuple(m % 17 for m in LIFTS)
FIBRE_LIFTS = tuple((m - (m % 17)) // 17 for m in LIFTS)

# (shallow, owner, deep, carry, h, kappa, right root, D-right root)
EXPECTED_STATES = (
    (1, 4, 12, 3, 3, 0, 6, 6),
    (4, 1, 0, 9, 9, 1, 6, 6),
)
EXPECTED_ACTIVE_SOURCES = (
    (1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12),
    (1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12),
)
EXPECTED_UNIT_SOURCES = (
    (1, 2, 3, 4, 5, 7, 8, 9, 10, 12),
    (1, 3, 4, 5, 6, 8, 9, 10, 11, 12),
)
EXPECTED_RIGHT_VECTORS = (
    (0, 8, 9, 12, 0, 0, 0),
    (3, 0, 0, 0, 12, 2, 0),
)
EXPECTED_RIGHT_DETERMINANTS = (4, 3)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def frac(value):
    return value - value.numerator // value.denominator


def strict_in_unweighted(point, intervals):
    return any(Fraction(left) < point < Fraction(right)
               for left, right in intervals)


def strict_in_weighted(point, pieces):
    return any(weight and Fraction(left) < point < Fraction(right)
               for left, right, weight in pieces)


def clock_label(x, speed):
    return int((7 * frac(speed * x) + Fraction(1, 2)) // 1) % 7


def clock_component(x, speed):
    """Open x-component on which the nearest-seven clock is constant."""
    image = speed * x
    branch = image.numerator // image.denominator
    walls = []
    for n in range(branch - 2, branch + 3):
        for j in range(7):
            walls.append(Fraction(14 * n + 2 * j + 1, 14 * speed))
    return (max(w for w in walls if w < x),
            min(w for w in walls if w > x))


def component_unweighted(point, intervals):
    require(strict_in_unweighted(point, intervals),
            "requested unweighted component is not strict")
    return old.containing_unweighted(point, intervals)


def component_weighted(point, pieces):
    require(strict_in_weighted(point, pieces),
            "requested weighted component is not strict")
    left, right, weight = old.containing_weighted(point, pieces)
    require(weight, "selected weighted component has zero weight")
    return left, right


def add_absolute_constraint(rows, name, x, left, right):
    require(left < x < right, f"{name} is not strict at the cycle point")
    rows.append((left - x, right - x, name))


def add_relative_constraint(rows, name, left, right):
    require(left < 0 < right, f"{name} relative component is not strict")
    rows.append((left, right, name))


def state_digits(x):
    rx = R * x
    predecessor = rx.numerator // rx.denominator
    y = frac(rx)
    h = (P * y).numerator // (P * y).denominator
    kappa = (2 * P * y).numerator // (2 * P * y).denominator - 2 * h
    return predecessor, predecessor % P, h, kappa, y


def state_constraints(module, pair_prefixes, present, rail, point,
                      state, sector=0):
    """All strict geometric factors needed to freeze one selected state."""
    shallow, owner, _deep, carry, h, kappa, right_root, d_right_root = state
    point_grid = point * old.T
    predecessor, actual_carry, actual_h, actual_kappa, y = state_digits(point)
    require((actual_carry, actual_h, actual_kappa) == (carry, h, kappa),
            "selected physical digit state changed")
    y_grid = y * old.T
    rows = []

    left, right = component_weighted(point_grid, rail[3])
    add_absolute_constraint(rows, "rail", point,
                            Fraction(left, old.T), Fraction(right, old.T))

    present_word = present[shallow, (-h) % P]
    left, right = component_unweighted(point_grid, present_word)
    add_absolute_constraint(rows, "present", point,
                            Fraction(left, old.T), Fraction(right, old.T))

    add_absolute_constraint(rows, "carry", point,
                            Fraction(predecessor, R),
                            Fraction(predecessor + 1, R))
    half_digit = 2 * h + kappa
    add_absolute_constraint(
        rows, "half_digit", point,
        Fraction(26 * predecessor + half_digit, 26 * R),
        Fraction(26 * predecessor + half_digit + 1, 26 * R),
    )

    delayed_word = old.prefix_intervals(
        pair_prefixes[sector][shallow][h][kappa]
    )
    left, right = component_unweighted(y_grid, delayed_word)
    add_absolute_constraint(
        rows, "delayed_Q", point,
        Fraction(predecessor, R) + Fraction(left, old.T * R),
        Fraction(predecessor, R) + Fraction(right, old.T * R),
    )

    left, right = old.comb_component(
        point_grid, module.C3, 14 * right_root, 14 * right_root + 13
    )
    add_absolute_constraint(rows, "current_right_root", point,
                            Fraction(left, old.T), Fraction(right, old.T))

    d_point = frac(P * point)
    left, right = old.comb_component(
        d_point * old.T, module.C3,
        14 * d_right_root, 14 * d_right_root + 13,
    )
    add_relative_constraint(
        rows, "D_right_root",
        (Fraction(left, old.T) - d_point) / P,
        (Fraction(right, old.T) - d_point) / P,
    )

    for name, speed in (("shallow_clock", P), ("owner_clock", P**2)):
        left, right = clock_component(point, speed)
        add_absolute_constraint(rows, name, point, left, right)

    return tuple(rows)


def prepare_unweighted(intervals):
    rows = tuple((17 * left, 17 * right) for left, right in intervals)
    return tuple(left for left, _ in rows), rows


def prepare_weighted(pieces):
    rows = tuple((17 * left, 17 * right, weight)
                 for left, right, weight in pieces)
    return tuple(left for left, _, _ in rows), rows


def inside_prepared_unweighted(point, prepared):
    starts, rows = prepared
    index = bisect_right(starts, point) - 1
    return index >= 0 and rows[index][0] < point < rows[index][1]


def inside_prepared_weighted(point, prepared):
    starts, rows = prepared
    index = bisect_right(starts, point) - 1
    return (index >= 0 and rows[index][2]
            and rows[index][0] < point < rows[index][1])


def clock_from_stalk_numerator(numerator, speed):
    residue = speed * numerator % L
    return ((14 * residue + L) // (2 * L)) % 7


def targeted_fibre_audit(module, pair_prefixes, rails,
                         present, present_starts):
    """Exhaust the literal-root subuniverse over 4/17 and 13/17.

    Literal right-root gluing forces carry 3 over 4/17 and carry 9 over
    13/17, because both unshifted D-right roots are 6.  Within those fibres
    retain a nonconstant intrinsic clock edge, sector zero, the matching
    present packet, at least one physical rail, and at least one primitive
    right-edge unit.  Pair the two sides by reverse clock labels and a common
    unit source.  Any such pair has a unique pair of affine interpolants.
    """
    require(old.T % R == 0, "carrier grid stopped resolving the R fibres")
    unit = old.T // R
    rail_prepared = tuple(prepare_weighted(rail[3]) for rail in rails)
    rails_by_owner = {owner: [] for owner in range(7)}
    for index, rail in enumerate(rails):
        rails_by_owner[rail[1]].append(index)
    present_prepared = {
        (shallow, h): prepare_unweighted(present[shallow, (-h) % P])
        for shallow in range(7) for h in range(P)
    }

    def geometry(a, carry, h, kappa):
        y_grid = Fraction(a * old.T, 17)
        sector_ok = {
            shallow: strict_in_unweighted(
                y_grid,
                old.prefix_intervals(
                    pair_prefixes[0][shallow][h][kappa]
                ),
            )
            for shallow in range(7)
        }
        rows = []
        combinations = set()
        for predecessor in range(carry, R, P):
            numerator = 17 * predecessor + a
            point_scaled_17 = numerator * unit
            shallow = clock_from_stalk_numerator(numerator, P)
            owner = clock_from_stalk_numerator(numerator, P**2)
            if shallow == owner or not sector_ok[shallow]:
                continue
            if not inside_prepared_unweighted(
                    point_scaled_17, present_prepared[shallow, h]):
                continue
            active = tuple(
                index for index in rails_by_owner[owner]
                if inside_prepared_weighted(
                    point_scaled_17, rail_prepared[index]
                )
            )
            if not active:
                continue
            rows.append((numerator, shallow, owner, active))
            combinations.update((index, shallow) for index in active)
        return tuple(rows), tuple(sorted(combinations))

    geometry_rows = []
    unit_rows = []
    parameters = ((4, 3, 3, 0), (13, 9, 9, 1))
    for a, carry, h, kappa in parameters:
        rows, combinations = geometry(a, carry, h, kappa)
        geometry_rows.append(rows)
        unit_table = {}
        for rail_index, shallow in combinations:
            root, _vector, determinant = old.unit_vector(
                module, pair_prefixes, present, present_starts,
                rails[rail_index], 0, 1, carry, h, kappa,
            )
            require(root == 6,
                    "literal-root fibre acquired the wrong right root")
            unit_table[rail_index, shallow] = determinant
        typed = []
        for numerator, shallow, owner, active in rows:
            unit_indices = tuple(
                index for index in active
                if unit_table[index, shallow]
            )
            sources = tuple(sorted({rails[index][0]
                                    for index in unit_indices}))
            if sources:
                typed.append((numerator, shallow, owner,
                              sources, unit_indices))
        unit_rows.append(tuple(typed))

    require(tuple(len(rows) for rows in geometry_rows) == (6260, 6252),
            "literal-root pre-unit fibre census changed")
    require(tuple(len(rows) for rows in unit_rows) == (6260, 5966),
            "literal-root primitive-unit fibre census changed")

    following_by_edge = {}
    for row in unit_rows[1]:
        following_by_edge.setdefault((row[1], row[2]), []).append(row)
    paired_initial_states = 0
    pair_count = 0
    first_pair = None
    for current in unit_rows[0]:
        candidates = following_by_edge.get((current[2], current[1]), ())
        has_pair = False
        for following in candidates:
            common = tuple(sorted(
                set(current[3]).intersection(following[3])
            ))
            if not common:
                continue
            has_pair = True
            pair_count += 1
            if first_pair is None:
                first_pair = (current, following, common)
        paired_initial_states += has_pair
    require((paired_initial_states, pair_count) == (3441, 902_300),
            "literal-root reverse-clock pair census changed")
    require(first_pair is not None
            and (first_pair[0][0], first_pair[1][0]) == POINT_NUMERATORS
            and first_pair[2]
            == (1, 3, 4, 5, 8, 9, 10, 12),
            "canonical literal-root pair stopped being the first witness")
    return {
        "geometry_counts": tuple(len(rows) for rows in geometry_rows),
        "unit_counts": tuple(len(rows) for rows in unit_rows),
        "paired_initial_states": paired_initial_states,
        "pair_count": pair_count,
        "first_pair_sources": first_pair[2],
    }


def main():
    require(L == 82_055_753,
            "enlarged transverse odometer denominator changed")
    require(TRANSVERSE_PHASES == (12, 5),
            "transverse phase residues changed")
    require(FIBRE_LIFTS == (1_472_984, 4_502_452)
            and tuple(k % P for k in FIBRE_LIFTS) == (6, 6),
            "13-primary fibre lifts changed")

    points = tuple(Fraction(n, L) for n in POINT_NUMERATORS)
    for index in range(2):
        following = frac(P * points[index] + Fraction(LIFTS[index], L))
        require(following == points[1 - index],
                "transverse affine maps stopped forming the exact two-cycle")

    delayed = tuple(frac(R * point) for point in points)
    require(delayed == (Fraction(4, 17), Fraction(13, 17)),
            "delayed chamber two-cycle changed")
    require(frac(P * delayed[0] + Fraction(12, 17)) == delayed[1]
            and frac(P * delayed[1] + Fraction(5, 17)) == delayed[0],
            "12/17,5/17 delayed arrows changed")

    # Exact enlarged-odometer skew product.  Write m=17k+a and Rx=n+y:
    # y'={13y+a/17}, n'=13n+k+floor(13y+a/17).
    skew_rows = []
    for index, (point, m, a, k) in enumerate(zip(
            points, LIFTS, TRANSVERSE_PHASES, FIBRE_LIFTS)):
        rx = R * point
        n = rx.numerator // rx.denominator
        y = frac(rx)
        phase_value = P * y + Fraction(a, 17)
        digit = phase_value.numerator // phase_value.denominator
        next_n = P * n + k + digit
        next_y = frac(phase_value)
        following = points[1 - index]
        require(R * following == (next_n % R) + next_y,
                "enlarged-odometer skew product failed")
        skew_rows.append((a, k, k % P, digit, next_n % P, next_y))
    require(tuple(row[4] for row in skew_rows) == (9, 3),
            "transverse carry handoff changed")

    # The correct compressed state is nu=N mod 221.  It simultaneously
    # retains the delayed C17 colour, predecessor carry, and exact high-probe
    # microphase.  The symmetric stalk is forced by
    #
    #     13 nu + m = -nu,       2m=7                 (mod 221).
    #
    # Both 2 and 14 are units, so m=7/2 and nu=-1/4 are unique.
    stalk_modulus = 13 * 17
    nu = tuple(n % stalk_modulus for n in POINT_NUMERATORS)
    lift_stalk = tuple(m % stalk_modulus for m in LIFTS)
    require(nu == (55, 166) and nu[1] == (-nu[0]) % stalk_modulus,
            "reflection-symmetric C221 states changed")
    require(lift_stalk == (114, 107)
            and lift_stalk[1] == (-lift_stalk[0]) % stalk_modulus,
            "opposite C221 lifts changed")
    forced_m = 7 * pow(2, -1, stalk_modulus) % stalk_modulus
    forced_nu = -forced_m * pow(P + 1, -1, stalk_modulus) % stalk_modulus
    require((forced_m, forced_nu) == (114, 55)
            and forced_nu == (-pow(4, -1, stalk_modulus)) % stalk_modulus,
            "slope-seven/quarter-torsion uniqueness changed")
    high_states = tuple(2 * value % stalk_modulus for value in nu)
    d_high_states = tuple(P * value % stalk_modulus for value in high_states)
    high_shifts = tuple(2 * value % stalk_modulus for value in lift_stalk)
    require(high_states == (110, 111)
            and d_high_states == (104, 117)
            and high_shifts == (7, 214),
            "refined high-probe microphase ledger changed")
    require(tuple((d_high_states[i] + high_shifts[i]) % stalk_modulus
                  for i in range(2)) == tuple(reversed(high_states)),
            "refined high-probe phases stopped transporting literally")

    # All four affine interpolants exist at packet-support level.  Intrinsic
    # clock gluing later deletes the two loops and retains only the cross maps.
    interpolants = tuple(
        tuple((POINT_NUMERATORS[j] - P * POINT_NUMERATORS[i]) % L
              for j in range(2))
        for i in range(2)
    )
    require(tuple(tuple(m % 17 for m in row) for row in interpolants)
            == ((3, 12), (5, 14)),
            "complete two-state transverse phase graph changed")

    # The enlarged cyclic stalk has kernel C_(13^5) over C_221, but the
    # shared factor 13 prevents a cyclic splitting of this exact sequence.
    require(L == stalk_modulus * P**5,
            "C221 stalk/kernel factorization changed")
    subgroup_generator = L // stalk_modulus
    require(gcd(subgroup_generator, stalk_modulus) == P,
            "order-221 subgroup unexpectedly became quotient-surjective")

    (module, _prefixes, _whole_prefixes, _digit_masses, rails,
     present, present_starts) = old.cross.build_carrier_data()
    pair_prefixes = old.private.build_pair_prefixes(module)

    selected_rails = []
    unit_rows = []
    all_constraints = []
    for event, (point, expected) in enumerate(zip(points, EXPECTED_STATES)):
        shallow, owner, deep, carry, h, kappa, right_root, d_right_root = expected
        require((clock_label(point, P), clock_label(point, P**2))
                == (shallow, owner) and shallow != owner,
                f"event {event} intrinsic nonconstant clock edge changed")
        require(owner == EXPECTED_STATES[1 - event][0],
                f"event {event} owner no longer glues to following shallow")

        point_grid = point * old.T
        predecessor, actual_carry, actual_h, actual_kappa, y = state_digits(point)
        require((actual_carry, actual_h, actual_kappa)
                == (carry, h, kappa),
                f"event {event} carry/future labels changed")

        require(strict_in_unweighted(
            point_grid, present[shallow, (-h) % P]
        ), f"event {event} left its present packet")

        sectors = tuple(
            sector for sector in range(2)
            if strict_in_unweighted(
                y * old.T,
                old.prefix_intervals(
                    pair_prefixes[sector][shallow][h][kappa]
                ),
            )
        )
        require(sectors == (0,),
                f"event {event} delayed sector typing changed")

        require(old.root_memberships(module, point_grid)
                == ((0, right_root + 1), (1, right_root)),
                f"event {event} current private roots changed")
        d_point_grid = frac(P * point) * old.T
        require((1, d_right_root)
                in old.root_memberships(module, d_point_grid),
                f"event {event} D-right root changed")

        active = tuple(
            (index, rail) for index, rail in enumerate(rails)
            if rail[1] == owner and strict_in_weighted(point_grid, rail[3])
        )
        active_sources = tuple(rail[0] for _, rail in active)
        require(active_sources == EXPECTED_ACTIVE_SOURCES[event],
                f"event {event} active rail sources changed")
        rail_index, selected = next(
            row for row in active
            if row[1][0] == 1 and row[1][2] == deep
        )
        selected_rails.append((rail_index, selected))

        per_edge_sources = {}
        per_edge_selected = {}
        for edge in (0, 1):
            sources = []
            for index, rail in active:
                root, vector, determinant = old.unit_vector(
                    module, pair_prefixes, present, present_starts,
                    rail, 0, edge, carry, h, kappa,
                )
                require((edge, root)
                        in old.root_memberships(module, point_grid),
                        f"event {event} unit root is not physically present")
                if determinant:
                    sources.append(rail[0])
                if index == rail_index:
                    per_edge_selected[edge] = (root, vector, determinant)
            per_edge_sources[edge] = tuple(sources)
        require(per_edge_sources[0] == per_edge_sources[1]
                == EXPECTED_UNIT_SOURCES[event],
                f"event {event} primitive-unit source census changed")
        require(per_edge_selected[1]
                == (right_root, EXPECTED_RIGHT_VECTORS[event],
                    EXPECTED_RIGHT_DETERMINANTS[event]),
                f"event {event} selected right-unit row changed")
        require(per_edge_selected[0][2],
                f"event {event} selected left edge stopped being a unit")
        unit_rows.append((active_sources, per_edge_sources[1],
                          per_edge_selected))

        constraints = state_constraints(
            module, pair_prefixes, present, selected, point, expected
        )
        all_constraints.extend(
            (left, right, f"event{event}:{name}")
            for left, right, name in constraints
        )

    common_unit_sources = tuple(sorted(
        set(EXPECTED_UNIT_SOURCES[0]).intersection(EXPECTED_UNIT_SOURCES[1])
    ))
    require(common_unit_sources == (1, 3, 4, 5, 8, 9, 10, 12),
            "same-source unit intersection changed")

    fibre_audit = targeted_fibre_audit(
        module, pair_prefixes, rails, present, present_starts
    )

    eta = min(min(-left, right) for left, right, _ in all_constraints)
    require(eta == Fraction(11, 853_068_347_561_612),
            "uniform strict packet radius changed")
    limiting = tuple(
        name for left, right, name in all_constraints
        if min(-left, right) == eta
    )
    require(limiting == ("event0:delayed_Q", "event1:delayed_Q"),
            "uniform radius stopped being delayed-word limited")

    # Exact all-H cylinder.  If x=x0+d and |d|<eta/13^(H-1), then the j-th
    # event is x_(j mod 2)+13^j d.  Every retained factor is strict on the
    # radius-eta neighbourhood of its parity point.
    for horizon in range(1, 13):
        radius = eta / P ** (horizon - 1)
        for displacement in (-radius / 2, Fraction(0), radius / 2):
            point = points[0] + displacement
            for event in range(horizon):
                parity = event % 2
                expected_point = points[parity] + P**event * displacement
                require(point == expected_point,
                        "affine perturbation law changed")
                for left, right, _name in (
                        row[:3] for row in all_constraints
                        if row[2].startswith(f"event{parity}:")):
                    require(left < point - points[parity] < right,
                            "all-H probe left a selected strict factor")
                point = frac(
                    P * point + Fraction(LIFTS[parity], L)
                )

    print("LRC14 transverse mod-17 locally packet-typed two-cycle probe")
    print("status=FINITE-EXACT DESIGN TARGET; positive packet cylinders for every finite H")
    print(f"parameters=R:{R},L=17R:{L}")
    print(f"cycle_numerators={POINT_NUMERATORS}; cycle_points={points}")
    print(f"lifts={LIFTS}; transverse_phases={TRANSVERSE_PHASES}; fibre_lifts={FIBRE_LIFTS}; fibre_mod13={tuple(k % P for k in FIBRE_LIFTS)}")
    print(f"delayed_cycle={delayed}; delayed_arrows=(+13,+12/17),(+13,+5/17)")
    print("two_step_C17_holonomies=phase_word(12,5):r->-r+8 fixes r=4 uniquely; phase_word(5,12):r->-r+9 fixes r=13 uniquely")
    print("skew_product=m=17k+a; y'={13y+a/17}; n'=13n+k+floor(13y+a/17)")
    print(f"skew_rows=(a,k,k_mod13,phase_floor,next_carry,next_y)={tuple(skew_rows)}")
    print(f"C221_states={nu}; C221_lifts={lift_stalk}; high_states_2nu={high_states}; unshifted_D_high={d_high_states}; transverse_high_shifts={high_shifts}")
    print("forced_C221_law=13nu+m=-nu and 2m=7 has unique solution m=7/2=114, nu=-1/4=55 mod221")
    print(f"all_packet_interpolants={interpolants}; phase_matrix=((3,12),(5,14)); clocks_delete_loops=True")
    print("stalk_extension=0->C_(13^5)->C_(17*13^6)->C_221->0; nonsplit_on_shared_C13")
    print("event_columns=(shallow,owner,deep,carry,h,kappa,right_root,D_right_root,rail_index,active_sources,unit_sources,right_vector,right_det,left_det)")
    for state, (rail_index, _rail), row in zip(
            EXPECTED_STATES, selected_rails, unit_rows):
        active_sources, unit_sources, selected = row
        print(state + (rail_index, active_sources, unit_sources,
                       selected[1][1], selected[1][2], selected[0][2]))
    print(f"same_source_unit_intersection={common_unit_sources}; count={len(common_unit_sources)}")
    print(f"targeted_literal_root_fibre_audit={fibre_audit}")
    print(f"uniform_strict_radius_eta={eta}; limiting_factors={limiting}")
    print("all_H_start_cylinder=(x0-eta/13^(H-1),x0+eta/13^(H-1))")
    print(f"all_H_cylinder_width=({2 * eta})/13^(H-1)")
    print("root_interface=current_right_roots=(6,6); D_right_roots=(6,6); literal coarse root gluing on both edges")
    print("refined_root_transport=(104+7=111,117-7=110) mod221; this is literal C221 microphase transport, not the inherited coarse F13 additive-root law")
    print("typing_boundary=transverse phases a/17 are not THM2657 integer/R translations; the C221 state plus C_(13^5) lift, not C17 colour alone, determines the physical packet")
    print("scope=rail,present,sector,carry,half-digit,current/D-root,primitive-unit,clock-glued finite-H support; no semantic endpoint, scalar row, row exclusion, or LRC14")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
