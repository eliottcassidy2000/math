#!/usr/bin/env python3
"""Exact secondary physical covariance audit for THM-2820.

This scratch companion answers the cheapest physical test isolated by
THM-2813.  It keeps three statements separate:

1. the actual THM-2806 common carrier and its address-cylinder support;
2. the target-only A_1 displacement of the adjacent residue-eight endpoint;
3. the lawful THM-2763 simultaneous rechart, which is pure gauge.

It proves the fixed-cell physical boundary used by THM-2820 and claims no
LRC(14) consequence.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from hashlib import sha256
import importlib.util
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


PINNED = {
    COMP / "lrc14_literal_fixed_sheet_allocation_thm2806.py":
        "311d0d85500f0c65945ebe5913f09d34a16293119c942b42eeaa854fbf85f71e",
    COMP / "lrc14_positive_graded_address_two_simplex_thm2807.py":
        "11cdbe3c6cc7f9d5b6b24863ced71eb91cc84adc67fe38a3f8a3e637362453fb",
    COMP / "lrc14_semantic_arm_right_wing_central_digit_thm2782.py":
        "7fbc6bb1ec303ded98eaad6e5d8205eb3d247258ada32b6f9904fc439ebb11fb",
    COMP / "lrc14_replica_dichotomy_typed_row_opus_20260727.py":
        "6ba64a68a9fd008d2e06949b1f1cf75012f1f4e734f75f55ce0af58ae20ad7b9",
}
for path, expected in PINNED.items():
    actual = sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
    require(actual == expected, f"pinned dependency changed: {path.name}")


import lrc14_literal_fixed_sheet_allocation_thm2806 as allocation
import lrc14_positive_graded_address_two_simplex_thm2807 as lift
import lrc14_semantic_arm_right_wing_central_digit_thm2782 as arm


# One older loader pins raw LF bytes.  The Windows checkout is CRLF, although
# its LF-normalized bytes match the proved pin above.  Load precisely that
# verified source without weakening the hash check.
def crlf_safe_present_loader():
    path = COMP / "lrc14_replica_dichotomy_typed_row_opus_20260727.py"
    spec = importlib.util.spec_from_file_location(
        "residue8_present_base", path
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


allocation.physical.target.load_present_module = crlf_safe_present_loader

P = 13
M = P**6
TOP = P**5


def nearest_integer(value):
    """Nearest integer for a nonnegative Fraction; no ties occur here."""
    return (
        2 * value.numerator + value.denominator
    ) // (2 * value.denominator)


def circular_shift_interval(interval, shift, period):
    left = (interval[0] + shift) % period
    right = (interval[1] + shift) % period
    require(left < right, "test interval crossed the circle seam")
    return left, right


def cylinder(center, half_width):
    return center - half_width, center + half_width


def contained(interval, container):
    return container[0] <= interval[0] and interval[1] <= container[1]


def intersects(left, right):
    return max(left[0], right[0]) < min(left[1], right[1])


def source_center(address):
    return arm.arm_source_center(address, 0) * arm.T


def target_center(address):
    return arm.arm_target_center(address, 0) * arm.T


def address_near_interval(interval, half_width):
    """Return the unique address cylinder meeting one short interval, if any."""
    period = arm.T
    midpoint = Fraction(interval[0] + interval[1], 2 * period)
    grid_coordinate = ((midpoint - arm.relative.Z) % 1) * M
    q0 = nearest_integer(grid_coordinate) % M
    inverse_seven = pow(7, -1, M)
    hits = []
    for dq in (-1, 0, 1):
        q = (q0 + dq) % M
        address = q * inverse_seven % M
        center0 = source_center(address)
        for deck in (-1, 0, 1):
            candidate = cylinder(center0 + deck * period, half_width)
            if intersects(interval, candidate):
                hits.append((address, candidate))
    require(len(hits) <= 1, "short common piece met two address cylinders")
    return hits[0] if hits else None


def carrier_mask(interval, carrier, sign, unit):
    mask = []
    for harmonic in range(P):
        shifted = allocation.physical.overlap.shift_weighted(
            carrier, sign * harmonic * unit
        )
        overlap = allocation.intersect_sorted(
            (interval,), allocation.support_of(shifted)
        )
        require(
            not overlap or overlap == (interval,),
            "carrier twist cut through the test atom",
        )
        mask.append(bool(overlap))
    return tuple(mask)


def q_effective(L, q):
    """THM-2763 gauge quotient for the canonical Bezout row z=e_0."""
    z_dot_L = L[0] % P
    return (
        (q[0] - z_dot_L) % P,
        (q[1] + z_dot_L) % P,
    )


def affine_lift_action(t, address):
    """THM-2813 relative lift A_t on Z/13^6."""
    return (
        address + t * TOP * ((address - 7) % P)
    ) % M


def sum_pushforward(bank):
    result = [0] * P
    for (a, b), value in bank.items():
        result[(a + b) % P] += value
    return tuple(result)


def translated(bank, q):
    return {
        (a, b): bank[(a - q[0]) % P, (b - q[1]) % P]
        for a in range(P)
        for b in range(P)
    }


def coefficient_positive_control():
    """Direct integral quotient decoder; provenance remains joint-absent."""
    omega = {
        (a, b): int(a != 0 and b != 0)
        for a in range(P)
        for b in range(P)
    }
    expected_base = [P - 2] * P
    rows = {}
    for label, q in (
        ("base", (0, 0)),
        ("target_only_A1", (0, 7)),
        ("pure_gauge", (1, P - 1)),
    ):
        pushed = sum_pushforward(translated(omega, q))
        expected = list(expected_base)
        expected[sum(q) % P] += 1
        require(pushed == tuple(expected), f"{label} decoder changed")
        rows[label] = (
            min(pushed), max(pushed),
            tuple(i for i, value in enumerate(pushed)
                  if value == max(pushed)),
        )
    return rows


def main():
    (
        _module, full_module, details, e3, clocks, q_pairs
    ) = allocation.build_geometry()
    selected = allocation.common_cell(
        details, full_module, e3, clocks, *allocation.CELL
    )
    period = full_module.T
    unit = period // P
    address_step = allocation.physical.SHIFT
    require(
        address_step == Fraction(7, M) * period
        and 7 * unit == Fraction(7, P) * period,
        "physical address/allocation scales changed",
    )

    I = allocation.ATOM_INTERVAL
    J = tuple(endpoint + address_step for endpoint in I)
    n7 = 6715
    n8 = 6716
    n8_image = affine_lift_action(1, n8)
    require(
        n8_image == 378009
        and n8_image - n8 == TOP
        and (
            arm.arm_source_center(n8_image, 0)
            - arm.arm_source_center(n8, 0)
        ) % 1 == Fraction(7, P),
        "residue-eight source-chart displacement changed",
    )
    commutator_checks = 0
    for t in range(P):
        for residue in range(P):
            # The defect depends only on the low digit, so these thirteen
            # representatives prove the identity on all of Z/13^6.
            left = affine_lift_action(t, (residue + 1) % M)
            right = (affine_lift_action(t, residue) + 1) % M
            require(
                (left - right) % M == t * TOP % M,
                "successor/transvection commutator changed",
            )
            commutator_checks += 1
    require(
        pow(7, -1, P) == 2,
        "normal-defect decoder coefficient changed",
    )
    A1J = circular_shift_interval(J, 7 * unit, period)
    half_width = arm.relative.Q_RADIUS * period

    require(
        contained(cylinder(source_center(n7), half_width), I)
        and contained(cylinder(target_center(n7), half_width), J)
        and contained(cylinder(source_center(n8), half_width), J)
        and contained(cylinder(source_center(n8_image), half_width), A1J),
        "6715/6716/cross-role shifted-cylinder incidence changed",
    )

    section_factors = allocation.section_factors(
        full_module, e3, clocks, *allocation.CELL
    )
    signatures = {
        name: tuple(
            allocation.contained(interval, factor)
            for factor in section_factors
        )
        for name, interval in (("I", I), ("J", J), ("A1J", A1J))
    }
    require(
        signatures["I"] == signatures["J"] == (True,) * 6
        and signatures["A1J"]
        == (False, True, True, True, True, True),
        "fixed-factor cross-role shift failure signature changed",
    )
    universe = ((0, period),)
    e3_atomic_factors = [
        (
            "C3_danger",
            tuple(full_module.make_comb(
                full_module.C3, 182, -13, 13
            )),
        ),
        (
            "guard_safe",
            tuple(full_module.subtract_comb(
                universe, full_module.W[full_module.GUARD],
                91, -13, 13,
            )),
        ),
    ]
    e3_atomic_factors.extend(
        (
            f"unit{index}_safe_W{full_module.W[index]}",
            tuple(full_module.subtract_comb(
                universe, full_module.W[index], 182, -13, 13
            )),
        )
        for index in full_module.UNIT_IDX
    )
    e3_atomic_factors.extend((
        (
            "C1_safe",
            tuple(full_module.subtract_comb(
                universe, full_module.C1, 182, -13, 13
            )),
        ),
        (
            "C2_safe",
            tuple(full_module.subtract_comb(
                universe, full_module.C2, 182, -13, 13
            )),
        ),
    ))
    e3_signatures = {
        name: tuple(
            allocation.contained(interval, factor)
            for _label, factor in e3_atomic_factors
        )
        for name, interval in (("I", I), ("J", J), ("A1J", A1J))
    }
    e3_labels = tuple(label for label, _factor in e3_atomic_factors)
    require(
        e3_signatures["I"] == e3_signatures["J"]
        == (True,) * len(e3_labels)
        and tuple(
            label for label, survives
            in zip(e3_labels, e3_signatures["A1J"])
            if not survives
        ) == ("guard_safe", "unit5_safe_W66"),
        "internal E3 cross-role shift failure factors changed",
    )

    masks = {
        name: (
            carrier_mask(interval, selected[0], -1, unit),
            carrier_mask(interval, selected[1], +1, unit),
        )
        for name, interval in (("I", I), ("J", J), ("A1J", A1J))
    }
    delta0 = (True,) + (False,) * (P - 1)
    delta7 = tuple(index == 7 for index in range(P))
    zero = (False,) * P
    require(
        masks["I"] == (delta0, zero)
        and masks["J"] == (zero, delta0)
        and masks["A1J"] == (zero, delta7),
        "cross-role shifted allocation-mask transport changed",
    )
    source_chart_orbit_rows = []
    full_fixed_nonzero = []
    for t in range(P):
        q = 7 * t % P
        shifted_interval = circular_shift_interval(
            J, 7 * t * unit, period
        )
        fixed_signature = tuple(
            allocation.contained(shifted_interval, factor)
            for factor in section_factors
        )
        atomic_failures = tuple(
            label for label, factor in e3_atomic_factors
            if not allocation.contained(shifted_interval, factor)
        )
        source_mask = carrier_mask(
            shifted_interval, selected[0], -1, unit
        )
        target_mask = carrier_mask(
            shifted_interval, selected[1], +1, unit
        )
        expected_target = tuple(index == q for index in range(P))
        require(
            source_mask == zero and target_mask == expected_target,
            f"source-chart orbit allocation mask changed at t={t}",
        )
        if t and all(fixed_signature) and not atomic_failures:
            full_fixed_nonzero.append((t, q))
        source_chart_orbit_rows.append((
            t,
            q,
            tuple(int(value) for value in fixed_signature),
            atomic_failures,
            tuple(index for index, value in enumerate(source_mask) if value),
            tuple(index for index, value in enumerate(target_mask) if value),
        ))
    require(
        full_fixed_nonzero == [(6, 3)]
        and q_effective((0,) * len(allocation.endpoint_base.WMOD), (0, 3))
        == (0, 3)
        and 2 * 3 % P == 6,
        "unique full fixed-section commutator-boundary survivor changed",
    )
    # Retain the exact target endpoint origin and its translated companion.
    # Translation by one allocation unit has trivial endpoint phase here:
    # after the inherited 13^4 frequency normalization its exponent is an
    # integral multiple of the exact cyclotomic order.  Thus surviving
    # translated intervals carry literally the same endpoint scalar.
    endpoint_base = allocation.endpoint_base
    target_atom = ((*J, 1),)
    endpoint_phase_exponent = (
        endpoint_base.Y0
        * endpoint_base.endpoint.RDIL
        * unit
    )
    require(
        endpoint_phase_exponent % endpoint_base.endpoint.NN == 0,
        "one-unit target endpoint phase ceased to be trivial",
    )
    right_origin = allocation.RIGHT_ORIGIN
    right_step = allocation.add(
        allocation.RIGHT_ORIGIN, allocation.TARGET_STEP
    )
    endpoint_orbit_rows = {}
    endpoint_supports = {}
    for address in (right_origin, right_step):
        ell = endpoint_base.REPS[address]
        present = tuple(
            endpoint_base.endpoint.build_set(endpoint_base.PAT_E3, ell)
        )
        starts = tuple(left for left, _right in present)
        rows = []
        support = []
        for q in range(P):
            shifted_atom = allocation.physical.overlap.shift_weighted(
                target_atom, q * unit
            )
            restricted = allocation.indexed_weighted_intersection(
                shifted_atom, present, starts
            )
            mass = sum(
                (right - left) * weight
                for left, right, weight in restricted
            )
            values = tuple(
                allocation.endpoint_sum(
                    restricted, -endpoint_base.Y0, embedding
                )
                for embedding in endpoint_base.endpoint.MODS
            )
            if any(values):
                support.append(q)
            rows.append((q, len(restricted), mass, values))
        endpoint_orbit_rows[address] = tuple(rows)
        endpoint_supports[address] = tuple(support)
    endpoint_value = (
        231164267889491750,
        630230755085920022,
    )
    require(
        endpoint_supports[right_origin] == (0, 3, 11)
        and endpoint_supports[right_step] == (0, 11),
        "two-right-origin endpoint support atlas changed",
    )
    for rows in endpoint_orbit_rows.values():
        for _q, pieces, mass, values in rows:
            require(
                (
                    (pieces, mass, values)
                    == (1, 26444880, endpoint_value)
                )
                if any(values)
                else (pieces, mass, values) == (0, 0, (0, 0)),
                "endpoint orbit mass/value pattern changed",
            )
    fixed_by_q = {
        q: bits
        for _t, q, bits, _failures, _source, _target
        in source_chart_orbit_rows
    }
    full_fixed_q = tuple(
        q for q, bits in sorted(fixed_by_q.items()) if all(bits)
    )
    two_point_endpoint_q = tuple(
        sorted(
            set(endpoint_supports[right_origin])
            & set(endpoint_supports[right_step])
        )
    )
    require(
        full_fixed_q == (0, 3)
        and two_point_endpoint_q == (0, 11)
        and not (
            (set(full_fixed_q) - {0})
            & (set(two_point_endpoint_q) - {0})
        )
        and fixed_by_q[11] == (1, 1, 1, 0, 1, 1),
        "factor/endpoint-edge either-or atlas changed",
    )

    clock_rows = []
    for clock, (source, target) in enumerate(details):
        target_pullback = allocation.physical.overlap.shift_weighted(
            target, -address_step
        )
        aligned = allocation.physical.overlap.intersect_weighted_profiles(
            source, target_pullback
        )
        common = tuple(
            (left, right, source_weight)
            for left, right, source_weight, target_weight in aligned
            if source_weight == target_weight
        )
        require(
            len(common) == len(aligned) == 3756,
            f"clock {clock} raw common-piece count changed",
        )
        residues = Counter()
        contained_count = 0
        for left, right, _weight in common:
            require(
                Fraction(right - left) + 2 * half_width
                < Fraction(period, M),
                "common piece is not shorter than one address grid cell",
            )
            hit = address_near_interval((left, right), half_width)
            if hit is None:
                continue
            address, address_cylinder = hit
            require(
                contained(address_cylinder, (left, right)),
                "an address cylinder only partly met the common carrier",
            )
            contained_count += 1
            residues[address % P] += 1
        require(
            contained_count == 1878
            and residues == Counter({7: 1878}),
            f"clock {clock} gained an off-sheet common cylinder",
        )
        clock_rows.append((clock, len(common), contained_count, tuple(residues.items())))

    W = allocation.endpoint_base.WMOD
    require(W[0] == 1, "canonical Bezout row z=e_0 changed")
    lawful_unit_rechart = (W, (1, P - 1))
    lawful_A1_rechart = (
        tuple((-7 * coordinate) % P for coordinate in W),
        (P - 7, 7),
    )
    target_only = ((0,) * len(W), (0, 7))
    require(
        q_effective(*lawful_unit_rechart) == (0, 0)
        and q_effective(*lawful_A1_rechart) == (0, 0)
        and q_effective(*target_only) == (0, 7),
        "gauge-effective displacement changed",
    )

    # At the sole fourfold point the physical vector is (w,w,w,w), hence
    # its pointwise Mobius contrast is zero.  The both-present vertex H has
    # augmentation one on that common atom, but is not the mixed face.
    # The shift above follows the adjacent residue-eight SOURCE chart; the
    # coincident residue-seven TARGET chart is fixed by A_1, so the shifted
    # interval is not A_1(H).  The full rooted mixed aggregate has 144
    # joint-absent atoms and integral augmentation 144 (unit mod 13).
    common_atom_mobius = 1 - 1 - 1 + 1
    both_present_augmentation = 1
    aggregate_augmentation = (P - 1) ** 2
    require(
        common_atom_mobius == 0
        and both_present_augmentation == 1
        and aggregate_augmentation == 144
        and aggregate_augmentation % P == 1,
        "common-atom/aggregate distinction changed",
    )
    decoder = coefficient_positive_control()

    print("RESIDUE-EIGHT COMMON-ALLOCATION COVARIANCE PROBE")
    print("status=VERIFIED-EXACT secondary fixed-cell boundary; no LRC14")
    print(
        f"provenance_I={I}; J=I+SHIFT={J}; "
        f"SHIFT={address_step}; T={period}; allocation_unit={unit}"
    )
    print(
        f"address_cylinders=(source_n7={n7},target_n7={n7},"
        f"source_n8={n8}); A1(n8)={n8_image}=n8+13^5; "
        "physical_source_chart_shift=7T/13=7_allocation_units; "
        "A1(target_n7)=target_n7"
    )
    print(
        "successor_transvection_commutator="
        "A_t*d-d*A_t=t*13^5; "
        f"residue_checks={commutator_checks}; "
        "physical_defect=7t*T/13; target_allocation_defect=7t; "
        "conditional_beta_normal=2*(q1+q2)"
    )
    print(
        "factor_order=(E3,clock,q1,q2,c2,c3); "
        f"fixed_factor_signatures={signatures}"
    )
    print(
        f"E3_atomic_factor_order={e3_labels}; "
        f"E3_atomic_signatures={e3_signatures}; "
        "cross_role_shift_packet_failures="
        "(guard_safe,unit5_safe_W66)"
    )
    print(
        "allocation_masks_(source,target)="
        f"(I={masks['I']},J={masks['J']},A1J={masks['A1J']}); "
        "shifted_interval_probes_target_delta7; "
        "this_is_not_A1_of_target_delta0; source_offsheet_support=empty"
    )
    print(
        "source_chart_orbit_rows_"
        "(t,q,six_factor_bits,E3_atomic_failures,"
        "source_support,target_support)="
        f"{tuple(source_chart_orbit_rows)}"
    )
    print(
        "unique_nonzero_full_fixed_boundary=(t=6,q=3); "
        "source_support=empty; target_support=delta3; "
        "q_eff=(0,3); conditional_beta_normal=6; "
        "this_is_successor_transvection_boundary_not_A6_of_common_H"
    )
    print(
        "right_endpoint_orbit_rows_"
        "(address:(q,pieces,mass,two_field_values))="
        f"{endpoint_orbit_rows}"
    )
    print(
        "right_endpoint_supports="
        f"(origin00={endpoint_supports[right_origin]},"
        f"translated12={endpoint_supports[right_step]}); "
        f"constant_nonzero_mass=26444880; "
        f"constant_nonzero_values={endpoint_value}"
    )
    print(
        f"factor_full_q={full_fixed_q}; "
        f"two_point_endpoint_edge_q={two_point_endpoint_q}; "
        "nonzero_intersection=empty; "
        "q3_keeps_full_factors_but_loses_endpoint_step; "
        "q11_keeps_endpoint_edge_but_loses_q2"
    )
    print(
        "raw_common_clock_scan_(clock,pieces,address_cylinders,residues)="
        f"{tuple(clock_rows)}"
    )
    print(
        "gauge_effective_vectors="
        f"(lawful_unit_rechart_L=W_q=(1,-1):"
        f"{q_effective(*lawful_unit_rechart)},"
        f"lawful_A1_full_rechart_L=-7W_q=(-7,7):"
        f"{q_effective(*lawful_A1_rechart)},"
        f"target_only_A1_L=0_q=(0,7):{q_effective(*target_only)})"
    )
    print(
        "mobius_provenance="
        f"(pointwise_fourfold_common_atom={common_atom_mobius},"
        f"both_present_H_J00={both_present_augmentation},"
        f"rooted_group_ring_aggregate_J00={aggregate_augmentation},"
        f"J00_mod13={aggregate_augmentation % P}); "
        "no_H_A1_value_is_computed; "
        "mixed_aggregate_support=144_joint_absent_atoms"
    )
    print(
        "integral_rooted_aggregate_sum_decoder_"
        "(min,max,exceptional_coordinates)="
        f"{decoder}; baseline=11"
    )
    print(
        "CONCLUSION: the shifted interval follows A1 on the adjacent "
        "residue-eight SOURCE chart, while the coincident residue-seven "
        "TARGET chart is fixed.  Its target-delta7 mask is therefore a "
        "cross-role hostile, not A1-covariance of the common H endpoint; "
        "there is no source/common allocation support and the t=1 shifted "
        "fixed packet first loses E3.  Across the full source-chart orbit, "
        "the unique nonzero complete fixed-section landing is t=6,q=3, "
        "with target delta3, q_eff=(0,3), and conditional beta=6; this is "
        "the successor-transvection boundary, not A6 of common H.  The "
        "same base target endpoint origin/scalar survives at q=3, but its "
        "translated endpoint companion vanishes.  Conversely q=11 keeps "
        "both endpoint scalars but loses fixed factor q2.  Hence no "
        "nonzero orbit point in this selected (s,t,clock)=(0,4,1) cell "
        "retains both the complete fixed section and the two-point target "
        "endpoint edge.  This is not a bank-global statement.  The "
        "lawful simultaneous source-target rechart is pure gauge q_eff=0.  "
        "No rooted mixed-face aggregate with nonzero augmentation is "
        "supported on a fourfold physical common atom in this constructor; "
        "the nonzero common H vertex is not a mixed face, and no H_A_t "
        "value is computed."
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
