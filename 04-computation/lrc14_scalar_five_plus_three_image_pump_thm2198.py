#!/usr/bin/env python3
"""Exact audit for THM-2198's scalar 5+3 image pump.

The load-bearing finite assertion is the exclusion of actual blocker
valuation profile (1,1,2).  At N=13^3, after normalizing the guard to one,
the depth-two blocker is safe on the primitive guard-safe annulus.  The two
depth-one blockers are represented by the 78 sign classes modulo 13^2, and
the five unit masks by the 1014 sign classes modulo 13^3.

For every unordered depth-one pair, including a diagonal pair, we remove its
two masks from the annulus.  The sum of the five largest individual unit-mask
intersections with the residual is still smaller than the residual.  Hence
even five optimally chosen unit masks cannot cover it.

The script also reconstructs every mask by the root-image integer-window
formula, independently of the direct torsion inequalities.  Arithmetic is
integer-exact; ``require`` remains active under ``python -O``.
"""

from collections import Counter, defaultdict
from hashlib import sha256


P = 13
M = 3
N = P ** M
Q = P ** (M - 1)


def require(condition, message):
    """Raise on a failed exact check, including under optimized Python."""
    if not condition:
        raise RuntimeError(message)


def norm_numerator(value, modulus):
    """Return modulus times the circle norm of value/modulus."""
    residue = value % modulus
    return min(residue, modulus - residue)


def sign_classes(modulus):
    """Represent the unit classes modulo sign by 1 <= a < modulus/2."""
    return tuple(
        a
        for a in range(1, (modulus + 1) // 2)
        if a % P != 0
    )


def direct_universe():
    """Primitive guard-safe annulus at depth three, guard coefficient one."""
    return tuple(
        z
        for z in range(N)
        if z % P != 0 and 7 * norm_numerator(z, N) > N
    )


def direct_mask(universe, multiplier):
    """Definition-level strict danger mask as a Python integer bitset."""
    answer = 0
    for index, z in enumerate(universe):
        if 14 * norm_numerator(multiplier * z, N) < N:
            answer |= 1 << index
    return answer


def integers_in_open_window(center_num, center_den, radius_num, radius_den):
    """Integers d with |center_num/center_den-d| < radius_num/radius_den."""
    center_floor = center_num // center_den
    answer = []
    # Every radius used below is less than two.
    for d in range(center_floor - 2, center_floor + 4):
        if (
            radius_den * abs(center_num - d * center_den)
            < radius_num * center_den
        ):
            answer.append(d)
    return tuple(answer)


def pump_guard_safe_roots(phase):
    """Depth-three numerators reconstructed from the guard root window."""
    unsafe_sheets = {
        d % P
        for d in integers_in_open_window(phase, Q, 13, 7)
    }
    # With guard one, the rooted sheet label is s=-k mod 13.
    return tuple(
        phase + k * Q
        for k in range(P)
        if (-k) % P not in unsafe_sheets
    )


def pump_unit_roots(phase, multiplier):
    """Unit-mask roots from the moving singleton/chord integer window."""
    inverse = pow(multiplier % P, -1, P)
    active_sheets = {
        (inverse * d) % P
        for d in integers_in_open_window(multiplier * phase, Q, 13, 14)
    }
    return {
        phase + k * Q
        for k in range(P)
        if (-k) % P in active_sheets
    }


def pump_unit_mask(universe_index, multiplier):
    """Root-image reconstruction of a unit mask on the annulus."""
    answer = 0
    for phase in range(Q):
        if phase % P == 0:
            continue
        active = pump_unit_roots(phase, multiplier)
        for z in pump_guard_safe_roots(phase):
            if z in active:
                answer |= 1 << universe_index[z]
    return answer


def pump_depth_one_mask(universe_index, unit_part):
    """A 13*unit_part mask is constant on each root fibre."""
    answer = 0
    for phase in range(Q):
        if phase % P == 0:
            continue
        active = 14 * norm_numerator(unit_part * phase, Q) < Q
        if active:
            for z in pump_guard_safe_roots(phase):
                answer |= 1 << universe_index[z]
    return answer


def expected_annulus_size():
    """Closed formula from primitive root counting at odd depth."""
    epsilon = -1
    numerator = 60 * N - 130 * epsilon
    require(numerator % 91 == 0, "annulus formula is not integral")
    return numerator // 91


def normalization_audit(universe):
    """Check all guard-unit changes of variables and sample coefficients."""
    normalized = set(universe)
    hostile_coefficients = (1, 6, 13, 183, 799, 169, 1007, 2196)
    hostile_points = (1, 2, 7, 13, 84, 168, 1098, 2196)
    checks = 0
    for guard in range(1, N):
        if guard % P == 0:
            continue
        inverse = pow(guard, -1, N)
        original = {
            z
            for z in range(N)
            if z % P != 0
            and 7 * norm_numerator(guard * z, N) > N
        }
        require(
            {guard * z % N for z in original} == normalized,
            "guard normalization did not carry its annulus to U",
        )
        for coefficient in hostile_coefficients:
            for z in hostile_points:
                require(
                    coefficient * z % N
                    == coefficient * inverse * (guard * z % N) % N,
                    "coefficient normalization identity failed",
                )
                checks += 1
    return checks


def run():
    universe = direct_universe()
    universe_index = {z: i for i, z in enumerate(universe)}
    full = (1 << len(universe)) - 1
    units = sign_classes(N)
    depth_one_units = sign_classes(Q)

    require(len(universe) == expected_annulus_size() == 1450,
            "bad primitive annulus size")
    require(len(units) == 1014, "bad depth-three unit sign-class count")
    require(len(depth_one_units) == 78,
            "bad depth-one blocker sign-class count")
    require(N % 7 != 0 and N % 14 != 0 and Q % 14 != 0,
            "a strict torsion endpoint appeared")
    normalization_checks = normalization_audit(universe)

    # Independent guard disintegration: every primitive depth-two phase has
    # 13 primitive lifts, and the root window selects exactly the direct
    # guard-safe ones.
    pump_universe = set()
    fibre_histogram = defaultdict(int)
    for phase in range(Q):
        if phase % P == 0:
            continue
        roots = pump_guard_safe_roots(phase)
        fibre_histogram[len(roots)] += 1
        pump_universe.update(roots)
    require(pump_universe == set(universe), "guard root pump mismatch")
    require(dict(sorted(fibre_histogram.items())) == {9: 110, 10: 46},
            "bad 9/10-sheet fibre histogram")

    direct_units = tuple(direct_mask(universe, a) for a in units)
    direct_depth_one = tuple(
        direct_mask(universe, P * a) for a in depth_one_units
    )
    require(
        all(
            mask == direct_mask(universe, N - label)
            for label, mask in zip(units, direct_units)
        ),
        "unit sign symmetry failed",
    )
    require(
        all(
            mask == direct_mask(universe, P * (Q - label))
            for label, mask in zip(depth_one_units, direct_depth_one)
        ),
        "depth-one sign symmetry failed",
    )

    # Reconstruct all 1092 masks without using the direct depth-three
    # multiplier inequality.
    pump_units = tuple(
        pump_unit_mask(universe_index, a) for a in units
    )
    pump_depth_one = tuple(
        pump_depth_one_mask(universe_index, a)
        for a in depth_one_units
    )
    require(pump_units == direct_units, "unit root-image reconstruction failed")
    require(
        pump_depth_one == direct_depth_one,
        "depth-one all-or-nothing fibre reconstruction failed",
    )

    # At actual valuation two the blocker is a nonzero thirteenth root on
    # every primitive numerator, so all six sign classes are empty.
    depth_two_sizes = tuple(
        direct_mask(universe, Q * a).bit_count() for a in range(1, 13)
    )
    require(depth_two_sizes == (0,) * 12,
            "a depth-two blocker was not safe")

    # Work with distinct mask sets.  Distinct positive integer coefficients
    # may reduce to the same sign class, but repetitions never enlarge a
    # union.  The only duplicate among the 1014 representatives is the empty
    # mask, so this quotient cannot hide a useful candidate.
    mask_to_label = {}
    duplicate_labels = defaultdict(list)
    for label, mask in zip(units, direct_units):
        duplicate_labels[mask].append(label)
        mask_to_label.setdefault(mask, label)
    duplicated_nonempty = {
        mask: tuple(labels)
        for mask, labels in duplicate_labels.items()
        if mask and len(labels) > 1
    }
    require(not duplicated_nonempty, "unexpected duplicated nonempty masks")
    empty_labels = tuple(duplicate_labels[0])
    require(empty_labels == (1, 1098), "bad duplicated-empty-mask control")
    distinct_unit_masks = tuple(mask_to_label)

    # Recover THM-2138's one-active-mask capacity direction as an inherited
    # positive control before testing the genuinely new two-active case.
    single_residual_histogram = Counter()
    single_best_margin = None
    single_best_record = None
    for label, shallow_mask in zip(depth_one_units, direct_depth_one):
        residual = full & ~shallow_mask
        residual_size = residual.bit_count()
        intersections = tuple(
            (residual & mask).bit_count() for mask in direct_units
        )
        maximum = max(intersections)
        margin = residual_size - 5 * maximum
        single_residual_histogram[residual_size] += 1
        maximizing_labels = tuple(
            unit_label
            for unit_label, size in zip(units, intersections)
            if size == maximum
        )
        record = (label, residual_size, maximum, maximizing_labels)
        if single_best_margin is None or margin < single_best_margin:
            single_best_margin = margin
            single_best_record = record
    require(
        single_residual_histogram
        == Counter({1210: 2, 1218: 1, 1222: 3, 1226: 21,
                    1228: 40, 1230: 11}),
        "single-active residual histogram changed",
    )
    require(single_best_margin == 30, "bad single-active minimum margin")
    require(
        single_best_record == (14, 1230, 240, (183,)),
        "bad single-active worst record",
    )

    best_margin = None
    best_record = None
    minimizers = 0
    pair_count = 0
    diagonal_count = 0
    crude_nonpositive = []
    margin_histogram = Counter()
    maximum_margin = None
    transcript = sha256()

    for i, left_label in enumerate(depth_one_units):
        for j in range(i, len(depth_one_units)):
            right_label = depth_one_units[j]
            pair_count += 1
            if i == j:
                diagonal_count += 1
            residual = full & ~(
                direct_depth_one[i] | direct_depth_one[j]
            )
            residual_size = residual.bit_count()
            intersections = sorted(
                (
                    ((residual & mask).bit_count(), mask_to_label[mask])
                    for mask in distinct_unit_masks
                ),
                key=lambda item: (-item[0], item[1]),
            )
            top_five = intersections[:5]
            margin = residual_size - sum(size for size, _ in top_five)
            require(margin > 0, "found a nonpositive five-mask margin")
            margin_histogram[margin] += 1
            maximum_margin = (
                margin if maximum_margin is None
                else max(maximum_margin, margin)
            )
            crude_margin = residual_size - 5 * intersections[0][0]
            if crude_margin <= 0:
                crude_nonpositive.append(
                    (
                        left_label,
                        right_label,
                        crude_margin,
                        margin,
                        intersections[0][0],
                    )
                )

            record = (
                left_label,
                right_label,
                residual_size,
                tuple(top_five),
            )
            if best_margin is None or margin < best_margin:
                best_margin = margin
                best_record = record
                minimizers = 1
            elif margin == best_margin:
                minimizers += 1

            all_intersections = sorted(
                (
                    ((residual & mask).bit_count(), label)
                    for label, mask in zip(units, direct_units)
                ),
                key=lambda item: (-item[0], item[1]),
            )
            require(
                tuple(all_intersections[:5]) == tuple(top_five),
                "the duplicate-mask quotient changed the top five",
            )
            transcript.update(
                (
                    f"{left_label},{right_label},{residual_size},"
                    f"{sum(size for size, _ in top_five)},{margin}:"
                ).encode("ascii")
            )
            for size, label in all_intersections:
                transcript.update(f"{label}={size},".encode("ascii"))
            transcript.update(b"\n")

    require(pair_count == 3081, "bad unordered-pair count")
    require(diagonal_count == 78, "diagonal pairs were not all retained")
    require(best_margin == 86, "bad global minimum margin")
    require(maximum_margin == 218, "bad global maximum margin")
    require(len(margin_histogram) == 66, "bad margin histogram support")
    require(sum(margin_histogram.values()) == 3081,
            "bad margin histogram mass")
    require(minimizers == 1, "minimum margin is not unique")
    require(
        tuple(crude_nonpositive)
        == (
            (7, 14, -22, 94, 228),
            (14, 61, 0, 102, 216),
            (23, 46, -2, 98, 224),
        ),
        "bad hostile rows for the crude repeated-maximum bound",
    )
    require(
        best_record
        == (
            14,
            46,
            1046,
            ((204, 183), (200, 799), (190, 599),
             (186, 1000), (180, 1007)),
        ),
        "bad worst-pair certificate",
    )
    table_digest = transcript.hexdigest()
    require(
        table_digest
        == "d31e5c874d8b5893ff33fc35095c18dcdf865c7f8296285c1b3441a8d8d679d9",
        "full conditional-capacity table digest changed",
    )
    worst_phase_histogram = Counter()
    worst_phases = []
    for phase in range(Q):
        if phase % P == 0:
            continue
        both_inactive = (
            14 * norm_numerator(14 * phase, Q) >= Q
            and 14 * norm_numerator(46 * phase, Q) >= Q
        )
        if both_inactive:
            worst_phases.append(phase)
            worst_phase_histogram[len(pump_guard_safe_roots(phase))] += 1
    require(
        worst_phase_histogram == Counter({9: 74, 10: 38}),
        "bad worst-residual image-fibre disintegration",
    )
    drift_rows = []
    for capacity, label in best_record[3]:
        singleton_events = 0
        outside_endpoints = 0
        reconstructed_capacity = 0
        inverse = pow(label % P, -1, P)
        for phase in worst_phases:
            integers = integers_in_open_window(
                label * phase, Q, 13, 14
            )
            singleton_events += len(integers) == 1
            active_sheets = {(inverse * d) % P for d in integers}
            unsafe_sheets = {
                d % P
                for d in integers_in_open_window(phase, Q, 13, 7)
            }
            outside_endpoints += len(active_sheets & unsafe_sheets)
            reconstructed_capacity += len(active_sheets - unsafe_sheets)
        require(
            reconstructed_capacity == capacity,
            "H-drift capacity reconstruction failed",
        )
        drift_rows.append(
            (label, capacity, singleton_events, outside_endpoints)
        )
    require(
        tuple(drift_rows)
        == (
            (183, 204, 0, 20),
            (799, 200, 0, 24),
            (599, 190, 10, 24),
            (1000, 186, 0, 38),
            (1007, 180, 12, 32),
        ),
        "bad H-drift component rows",
    )
    guard_c_frames = worst_phase_histogram[9]
    h_drift = (
        sum(singletons + outside for _, _, singletons, outside in drift_rows)
        - guard_c_frames
    )
    require(h_drift == best_margin == 86, "bad H-drift identity")

    # Positive-direction control for the capacity inequality: five named
    # masks do cover the union constructed from those same masks.
    control_labels = (3, 4, 5, 7, 11)
    control_masks = [
        direct_units[units.index(label)] for label in control_labels
    ]
    control_union = 0
    for mask in control_masks:
        control_union |= mask
    require(
        all((mask & control_union) == mask for mask in control_masks),
        "positive union control failed",
    )
    require(control_union.bit_count() == 1024,
            "positive union control size changed")

    print("THM-2198 scalar five-plus-three image-pump audit")
    print("arithmetic: exact integers only; no floats; no assert dependence")
    print()
    print("depth-three primitive guard annulus")
    print(f"  N={N} Q={Q} |U_N|={len(universe)}")
    print(
        "  primitive image phases=156 "
        f"fibre_histogram={tuple(sorted(fibre_histogram.items()))}"
    )
    print(
        "  unit_sign_classes=1014 depth_one_sign_classes=78 "
        "depth_two_nonzero_unit_parts=12"
    )
    print(
        "  direct/root-pump parity: "
        "guard=PASS units=1014/1014 depth_one=78/78"
    )
    print(f"  all_guard_normalization_checks={normalization_checks}")
    print(f"  depth_two_mask_sizes={depth_two_sizes}")
    print()
    print("first-depth capacity certificate")
    print(
        f"  unordered_depth_one_pairs={pair_count} "
        f"including_diagonal={diagonal_count}"
    )
    print(
        f"  distinct_unit_masks={len(distinct_unit_masks)} "
        f"duplicated_empty_labels={empty_labels}"
    )
    left, right, residual_size, top_five = best_record
    print(
        f"  unique_worst_pair=({left},{right}) "
        f"residual_size={residual_size}"
    )
    print(
        "  worst_residual_image_fibres="
        f"{tuple(sorted(worst_phase_histogram.items()))} "
        f"phases={sum(worst_phase_histogram.values())}"
    )
    print(f"  top_five_(intersection,label)={top_five}")
    print(
        f"  conditional_margin_range={best_margin}..{maximum_margin} "
        f"histogram_support={len(margin_histogram)}"
    )
    print(
        "  H_drift_(label,capacity,singletons,outside)="
        f"{tuple(drift_rows)}"
    )
    print(
        f"  H_drift=sum(singletons+outside)-C_frames="
        f"{h_drift} with C_frames={guard_c_frames}"
    )
    print(f"  full_table_sha256={table_digest}")
    print(
        "  consequence: every residual needs at least six distinct "
        "unit masks"
    )
    print()
    print("hostile controls")
    print(
        "  inherited single-active check: "
        f"minimum_margin={single_best_margin} "
        f"worst={single_best_record} PASS"
    )
    print(f"  crude_5max_nonpositive_rows={tuple(crude_nonpositive)}")
    print(
        "  diagonal blocker pairs retained; repeated coefficient residues "
        "quotiented only after union-idempotence"
    )
    print(
        f"  positive five-mask constructed-union control: "
        f"labels={control_labels} size={control_union.bit_count()} PASS"
    )
    print()
    print("CONCLUSION")
    print(
        "  actual blocker valuation profile (1,1,2) cannot cover the "
        "primitive guard-safe annulus"
    )
    print("  after THM-2192, every scalar 5+3 survivor has deepest depth >=3")
    print()
    print("reproduce:")
    print(
        "  python3 "
        "04-computation/lrc14_scalar_five_plus_three_image_pump_thm2198.py"
    )
    print(
        "  python3 -O "
        "04-computation/lrc14_scalar_five_plus_three_image_pump_thm2198.py"
    )


if __name__ == "__main__":
    run()
