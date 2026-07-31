#!/usr/bin/env python3
"""Exact rail-eight common/cofiber copy audit for THM-2806/THM-2771.

The selected common interval I and the two right-cofiber intervals J-, J+
have equal folded rail weight and equal delayed coefficient.  This audit
tests how far that equality lifts before a physical map fails:

* literal A, B, M=A∩B, and R=B\\A interval geometry;
* all native/pulled semantic-factor masks;
* delayed carry/root values and the literal THM-2584 ancestry label sets;
* carrier-twist masks and all 169 inherited endpoint-present masks;
* the exact split of the THM-2771 raw coefficient and intrinsic Bockstein.

No map is inferred from equal scalar content.
"""

from __future__ import annotations

from hashlib import sha256
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


PINNED = {
    "lrc14_fully_marked_root_zero_target_profile_thm2749.py":
        "d67c852c52f88feaadb2fcaa0a9a07a212f2e47018040b455855df886200595e",
    "lrc14_root_zero_clutch_mayer_vietoris_wing_shear_thm2751.py":
        "25cbed38026d61891173c687006250a69fe38aea56d67439406bd8bb60fa2552",
    "lrc14_root_zero_full_target_semantic_clutch_20260728.py":
        "208f71020efa19fa47f66d2da061ab03fa7bc87beeb077b4008c069f499736d8",
    "lrc14_full_arm_orbit_path_sheet_audit_thm2791.py":
        "1e00b6711db0d878fca70047f5f1532518084dbf6928551cd28fe51283fde543",
    "lrc14_extended_carrier_endpoint_lib.py":
        "4b3f9f195b1634e1e84a1bc8bccb878a1c8c44aec13f24d197f92547c9e36c57",
}
for name, expected in PINNED.items():
    actual = sha256(
        (COMP / name).read_bytes().replace(b"\r\n", b"\n")
    ).hexdigest()
    require(actual == expected, f"pinned dependency changed: {name}")


import lrc14_extended_carrier_endpoint_lib as endpoint_base
import lrc14_fully_marked_root_zero_target_profile_thm2749 as marked
import lrc14_full_arm_orbit_path_sheet_audit_thm2791 as ancestry
import lrc14_root_zero_clutch_mayer_vietoris_wing_shear_thm2751 as wing
import lrc14_root_zero_full_target_semantic_clutch_20260728 as physical


P = 13
T = physical.T
SHIFT = physical.SHIFT
DEPTH = P**5
CELL = (0, 4, 1)
I = (142004992589460, 142005019034340)
J_MINUS = tuple(endpoint - T // DEPTH for endpoint in I)
J_PLUS = tuple(endpoint + 96 * T // DEPTH for endpoint in I)
SHEETS = (("M", I), ("R-", J_MINUS), ("R+", J_PLUS))
W = 27581135604
C = 103478815440
G0 = 5905329039529920
K1 = 483303
FIBRE = 57068
PATH = (59162, 26, 56658)
FACTOR_NAMES = ("E3", "clock", "q1", "q2", "c2", "c3")


def intersect_sorted(left, right):
    out = []
    i = j = 0
    while i < len(left) and j < len(right):
        a = max(left[i][0], right[j][0])
        b = min(left[i][1], right[j][1])
        if a < b:
            out.append((a, b))
        if left[i][1] <= right[j][1]:
            i += 1
        else:
            j += 1
    return tuple(out)


def contained(interval, intervals):
    return intersect_sorted((interval,), intervals) == (interval,)


def support_of(weighted):
    return physical.overlap.merge_intervals(
        (left, right) for left, right, weight in weighted if weight
    )


def weighted_intersection(pieces, intervals):
    return tuple(
        physical.relative.private.old.intersect_weighted_union(
            pieces, intervals
        )
    )


def section_factors(full_module, e3, clocks):
    s, t, clock = CELL
    universe = ((0, full_module.T),)
    return (
        e3,
        clocks[clock],
        full_module.subtract_comb(
            universe, full_module.W[1], 182,
            -14 * s - 13, -14 * s + 13,
        ),
        full_module.subtract_comb(
            universe, full_module.W[2], 182,
            -14 * t - 13, -14 * t + 13,
        ),
        full_module.subtract_comb(
            universe, full_module.C2, 182,
            14 * s - 13, 14 * s + 13,
        ),
        full_module.subtract_comb(
            universe, full_module.C3, 182,
            14 * t - 13, 14 * t + 13,
        ),
    )


def build_physical_geometry():
    (
        module, _prefixes, _whole, _masses, rails, present, _starts
    ) = physical.relative.lift.m.core.build_carrier_data()
    pair_prefixes = physical.relative.private.build_pair_prefixes(module)
    _sv, _tv, rail_pairs, details = physical.overlap.overlap_vectors(
        module, pair_prefixes, rails, present, rail_index=8
    )
    full_module = physical.target.load_present_module()
    e3 = physical.target.exclusive_source(full_module, 3)
    fork = physical.target.deepest_fork(full_module)
    clocks = tuple(
        full_module.make_comb(
            full_module.C1, 182, 26 * clock - 13, 26 * clock + 13
        )
        for clock in range(7)
    )
    q_pairs = physical.q_restricted_pair_prefixes(
        full_module, pair_prefixes, fork
    )
    source_base, target_base = details[CELL[2]]
    section = physical.target.source_present_section(
        full_module, e3, CELL[2], CELL[0], CELL[1], clocks
    )
    source = weighted_intersection(source_base, section)
    target = weighted_intersection(target_base, section)
    target_pullback = physical.overlap.shift_weighted(target, -SHIFT)
    aligned = physical.overlap.intersect_weighted_profiles(
        source, target_pullback
    )
    require(
        all(a == b for _left, _right, a, b in aligned),
        "common profile lost equal source/target weights",
    )
    common = tuple(
        (left, right, a)
        for left, right, a, _b in aligned
    )
    right = physical.subtract_weighted(target_pullback, common)
    require(
        tuple(rails[8][:3]) == (1, 4, 12)
        and all(a == b for _left, _right, a, b in rail_pairs),
        "rail-eight metadata or aligned weight law changed",
    )
    return (
        module, rails, present, full_module, e3, clocks, q_pairs,
        source, target, target_pullback, common, right,
    )


def thm2771_geometry(module, rails):
    source = marked.two.exclusive_source(module, 3)
    clocks = tuple(
        module.make_comb(module.C1, 182, 26 * e - 13, 26 * e + 13)
        for e in range(7)
    )
    raw = tuple(marked.two.intersect_sorted(source, clocks[CELL[2]]))
    raw = tuple(module.subtract_comb(
        raw, module.W[1], 182, -13, 13
    ))
    raw = tuple(module.subtract_comb(
        raw, module.C2, 182, -13, 13
    ))
    raw = tuple(module.subtract_comb(
        raw, module.W[2], 182,
        -14 * CELL[1] - 13, -14 * CELL[1] + 13,
    ))
    raw = tuple(module.subtract_comb(
        raw, module.C3, 182,
        14 * CELL[1] - 13, 14 * CELL[1] + 13,
    ))
    _source_weight, _target_weight, rail_common = marked.rail_data(
        rails, marked.RAIL
    )
    a = marked.intersect(rail_common, raw)
    b = marked.intersect(
        rail_common, marked.shift_union(raw, marked.SHIFT)
    )
    common = marked.intersect(a, b)
    right = wing.difference(b, a)
    return a, b, common, right


def factor_masks(interval, factors):
    target_interval = tuple(endpoint + SHIFT for endpoint in interval)
    return (
        tuple(contained(interval, factor) for factor in factors),
        tuple(
            contained(
                interval,
                physical.overlap.shift_union(factor, -SHIFT),
            )
            for factor in factors
        ),
        tuple(contained(target_interval, factor) for factor in factors),
        tuple(
            contained(
                target_interval,
                physical.overlap.shift_union(factor, SHIFT),
            )
            for factor in factors
        ),
    )


def semantic_profile(interval, q_pair):
    source = ((*interval, 1),)
    target = physical.overlap.shift_weighted(source, SHIFT)
    source_values = physical.relative.private.delayed_carry_pair(
        source, q_pair, {}
    )
    target_values = physical.relative.private.delayed_carry_pair(
        target, q_pair, {}
    )
    source_nonzero = tuple(
        (carry, value)
        for carry, value in enumerate(source_values)
        if value != (0, 0)
    )
    target_nonzero = tuple(
        (carry, value)
        for carry, value in enumerate(target_values)
        if value != (0, 0)
    )
    midpoint = sum(interval) // 2
    arrival_root = P * midpoint // T
    deep_digit = 2 * P * midpoint // T % P
    source_root = (arrival_root - 1) % P
    address = P**6 * midpoint // T
    return (
        source_nonzero,
        target_nonzero,
        (source_root, arrival_root, deep_digit),
        address,
    )


def ancestry_data(intervals):
    e_set = tuple(
        ancestry.base.build_set(ancestry.base.PAT_E, ancestry.base.ZELL)
    )
    q_set = tuple(
        ancestry.base.build_set(
            ancestry.host.PAT_QB, ancestry.base.ZELL
        )
    )
    arguments = (
        *ancestry.scaled_intervals(q_set, DEPTH),
        *ancestry.scaled_intervals(e_set, DEPTH * P**2),
        *ancestry.scaled_intervals(e_set, DEPTH),
    )
    rows = []
    for interval in intervals:
        coordinate = sum(interval) // 2
        u_labels, v_labels = ancestry.contributor_sets(
            coordinate, *arguments
        )
        rows.append((u_labels, v_labels))
    first_u, first_v = rows[0]
    require(
        all(u == first_u and v == first_v for u, v in rows[1:]),
        "literal ancestry label sets differ",
    )
    packed_path = PATH[0] * P**2 + PATH[1]
    require(
        packed_path in first_u and PATH[2] in first_v,
        "supplied THM-2584 path is absent",
    )
    return (
        len(first_u), len(first_v),
        ancestry.path_digest(first_u, first_v),
        True,
    )


def carrier_twist_masks(interval, source, target):
    target_interval = tuple(endpoint + SHIFT for endpoint in interval)
    unit = T // P
    source_mask = tuple(
        bool(intersect_sorted(
            (interval,),
            support_of(physical.overlap.shift_weighted(
                source, -twist * unit
            )),
        ))
        for twist in range(P)
    )
    target_mask = tuple(
        bool(intersect_sorted(
            (target_interval,),
            support_of(physical.overlap.shift_weighted(
                target, twist * unit
            )),
        ))
        for twist in range(P)
    )
    return source_mask, target_mask


def endpoint_mask(interval, present_sets):
    values = []
    partial = []
    for address in endpoint_base.KEYS:
        overlap = intersect_sorted((interval,), present_sets[address])
        if overlap == (interval,):
            values.append(1)
        elif not overlap:
            values.append(0)
        else:
            values.append(2)
            partial.append((address, overlap))
    return tuple(values), tuple(partial)


def translated_masks(source, target):
    source_bank = dict(zip(endpoint_base.KEYS, source))
    target_bank = dict(zip(endpoint_base.KEYS, target))
    answers = []
    for delta in endpoint_base.KEYS:
        if all(
            target_bank[address]
            == source_bank[(
                (address[0] + delta[0]) % P,
                (address[1] + delta[1]) % P,
            )]
            for address in endpoint_base.KEYS
        ):
            answers.append(delta)
    return tuple(answers)


def endpoint_masks():
    present_sets = endpoint_base.present_cache()
    rows = {}
    for name, interval in SHEETS:
        source_mask, source_partial = endpoint_mask(
            interval, present_sets
        )
        target_interval = tuple(
            endpoint + SHIFT for endpoint in interval
        )
        target_mask, target_partial = endpoint_mask(
            target_interval, present_sets
        )
        require(
            not source_partial and not target_partial,
            "an endpoint mask cuts through a selected interval",
        )
        rows[name] = (source_mask, target_mask)
    expected = {
        ("M", "S"): (
            (88, 81),
            "681c43d56526b5b52b1f0bb3521f08cf8bf8ba9f3b51967a2283be760565e4ee",
        ),
        ("M", "T"): (
            (88, 81),
            "681c43d56526b5b52b1f0bb3521f08cf8bf8ba9f3b51967a2283be760565e4ee",
        ),
        ("R-", "S"): (
            (169, 0),
            "605d47a6802a6ba6675ce2970606011e1d53eebdd846effd6f47bd0903d7ed13",
        ),
        ("R-", "T"): (
            (88, 81),
            "681c43d56526b5b52b1f0bb3521f08cf8bf8ba9f3b51967a2283be760565e4ee",
        ),
        ("R+", "S"): (
            (79, 90),
            "e58ebf4b115e803796a8fd10516811a739a445a5db8b879f799a9db9f921c8de",
        ),
        ("R+", "T"): (
            (70, 99),
            "5243cf1be640698c0960708c6ae595b93b735a1188893f269d82fffbe18c34ae",
        ),
    }
    reports = {}
    for name in ("M", "R-", "R+"):
        for side, mask in zip(("S", "T"), rows[name]):
            counts = (mask.count(0), mask.count(1))
            digest = sha256(bytes(mask)).hexdigest()
            require(
                (counts, digest) == expected[name, side],
                f"endpoint mask changed at {(name, side)}",
            )
            reports[name, side] = (counts, digest)
    comparisons = {}
    for side_index, side in enumerate(("S", "T")):
        common = rows["M"][side_index]
        for name in ("R-", "R+"):
            cofiber = rows[name][side_index]
            differences = sum(
                left != right
                for left, right in zip(common, cofiber)
            )
            translations = translated_masks(common, cofiber)
            comparisons[side, name] = (differences, translations)
    require(
        comparisons == {
            ("S", "R-"): (81, ()),
            ("S", "R+"): (27, ()),
            ("T", "R-"): (0, ((0, 0),)),
            ("T", "R+"): (18, ()),
        },
        "endpoint-mask comparison changed",
    )
    return reports, comparisons


def coefficient_split(module, rails, present, right_support):
    delayed = marked.marked_prefixes(
        module,
        marked.private.build_pair_prefixes(module),
        marked.two.deepest_fork(module),
    )
    _source_weight, target_weight, _rail_common = marked.rail_data(
        rails, marked.RAIL
    )
    rows = []
    expected_vector = (0,) + (W * C,) * 6
    for _name, interval in SHEETS[1:]:
        vector, masses = wing.coefficient_vector(
            module, delayed, present, target_weight, (interval,), 12,
            [{} for _ in range(7)],
        )
        require(
            vector == expected_vector
            and masses == ((interval[1] - interval[0]) * W,) * 7,
            "one cofiber interval contribution changed",
        )
        rows.append((vector, masses))
    total_vector, _total_masses = wing.coefficient_vector(
        module, delayed, present, target_weight, right_support, 12,
        [{} for _ in range(7)],
    )
    require(
        total_vector == tuple(
            rows[0][0][index] + rows[1][0][index]
            for index in range(7)
        ),
        "right-cofiber coefficient is not the disjoint-copy sum",
    )
    require(
        FIBRE * K1 == W
        and FIBRE * C == G0
        and G0 * K1 == W * C
        and wing.scalar_amplitude(total_vector) == 2 * W * C
        and wing.scalar_amplitude(total_vector) // G0 == 2 * K1,
        "factor-two arithmetic changed",
    )
    copy_beta = (W * C // P) % P
    total_beta = (2 * W * C // P) % P
    require(
        W * C % P == 0
        and W * C % P**2 != 0
        and (copy_beta, total_beta) == (9, 5)
        and total_beta == 2 * copy_beta % P,
        "copywise intrinsic Bockstein changed",
    )
    return (
        expected_vector,
        total_vector,
        copy_beta,
        total_beta,
    )


def main():
    require(
        (P, T, SHIFT, DEPTH)
        == (13, 297836897838480, 431933040, 371293),
        "canonical scales changed",
    )
    require(
        J_MINUS == (142004190428100, 142004216872980)
        and J_PLUS == (142082000080020, 142082026524900),
        "depth-five translated intervals changed",
    )
    (
        module, rails, present, full_module, e3, clocks, q_pairs,
        source, target, _target_pullback, common_weighted, right_weighted,
    ) = build_physical_geometry()
    a, b, common, right = thm2771_geometry(module, rails)
    expected_right = (J_MINUS, J_PLUS)
    require(
        I in common
        and len(right) == 241
        and support_of(right_weighted) == expected_right
        and right_weighted
        == ((*J_MINUS, W), (*J_PLUS, W)),
        "effective THM-2771 right cofiber is not the two claimed copies",
    )
    require(
        I in tuple((left, right) for left, right, weight in common_weighted
                   if weight == W)
        and contained(I, a) and contained(I, b)
        and all(contained(interval, b) and not contained(interval, a)
                for interval in expected_right),
        "M/R typing changed",
    )

    factors = section_factors(full_module, e3, clocks)
    masks = {
        name: factor_masks(interval, factors)
        for name, interval in SHEETS
    }
    require(
        masks == {
            "M": ((True,) * 6,) * 4,
            "R-": (
                (False, True, True, True, True, True),
                (True,) * 6,
                (True,) * 6,
                (False, True, True, True, True, True),
            ),
            "R+": (
                (False, True, True, True, False, True),
                (True,) * 6,
                (True,) * 6,
                (False, True, True, True, False, True),
            ),
        },
        "native/pulled factor masks changed",
    )

    semantics = {
        name: semantic_profile(interval, q_pairs[CELL[2]])
        for name, interval in SHEETS
    }
    for name in semantics:
        source_nonzero, target_nonzero, roots, _address = semantics[name]
        require(
            source_nonzero == ((12, (0, C)),)
            and target_nonzero == ((6, (0, C)),)
            and roots == (5, 6, 12),
            f"semantic/carry/root data changed at {name}",
        )
    require(
        tuple(semantics[name][3] for name, _interval in SHEETS)
        == (2301363, 2301350, 2302611),
        "full addresses changed",
    )

    ancestry_report = ancestry_data(
        tuple(interval for _name, interval in SHEETS)
    )
    require(
        ancestry_report
        == (
            966606, 28534,
            "15c804c7cea9f61feab3b641eccdc035d937142b446d1cc14e059210eb1534fd",
            True,
        )
        and ancestry_report[0] * ancestry_report[1] == W,
        "ancestry report changed",
    )

    twist_masks = {
        name: carrier_twist_masks(interval, source, target)
        for name, interval in SHEETS
    }
    delta_zero = (True,) + (False,) * 12
    zero = (False,) * 13
    require(
        twist_masks == {
            "M": (delta_zero, delta_zero),
            "R-": (zero, delta_zero),
            "R+": (zero, delta_zero),
        },
        "carrier-twist masks changed",
    )

    endpoint_reports, endpoint_comparisons = endpoint_masks()
    coefficient = coefficient_split(
        module, rails, present, right
    )

    print("THM-2806/THM-2771 RAIL-EIGHT COFIBER COPY AUDIT")
    print("status=FINITE-EXACT scratch; no LRC conclusion")
    print(
        f"cell={CELL}; factor_order={FACTOR_NAMES}; "
        f"I={I}; Jminus={J_MINUS}; Jplus={J_PLUS}"
    )
    print(
        "translations=(Jminus-I=-T/13^5,Jplus-I=96T/13^5); "
        "full_address_shifts=(-13,+1248)"
    )
    print(
        f"objects=(I in M=A_intersect_B; "
        f"(Jminus,Jplus) in R=B_minus_A); "
        f"raw_R_piece_count={len(right)};"
        f"coefficient_effective_piece_count={len(right_weighted)}"
    )
    for name, _interval in SHEETS:
        print(
            f"{name}_factor_masks="
            f"(source_native={masks[name][0]},"
            f"source_pulled={masks[name][1]},"
            f"target_native={masks[name][2]},"
            f"target_adjacent={masks[name][3]})"
        )
    print(
        "first_factor_failure=(Rminus:E3,Rplus:E3_then_c2); "
        "both_pass_all_pulled_target_factors"
    )
    print(
        f"ancestry=(U={ancestry_report[0]},V={ancestry_report[1]},"
        f"product={ancestry_report[0] * ancestry_report[1]},"
        f"digest={ancestry_report[2]},identity_sets_on_all_three=yes,"
        f"supplied_path={PATH}:active)"
    )
    for name, _interval in SHEETS:
        print(
            f"{name}_semantic=(source_nonzero={semantics[name][0]},"
            f"target_nonzero={semantics[name][1]},"
            f"roots_deep={semantics[name][2]},"
            f"full_address={semantics[name][3]})"
        )
        print(
            f"{name}_carrier_twist_masks="
            f"(source={twist_masks[name][0]},"
            f"target={twist_masks[name][1]})"
        )
        for side in ("S", "T"):
            print(
                f"{name}_{side}_endpoint_mask="
                f"(zero_one_counts={endpoint_reports[name, side][0]},"
                f"digest={endpoint_reports[name, side][1]})"
            )
    print(f"endpoint_mask_comparisons={endpoint_comparisons}")
    print(
        f"copy_vector={coefficient[0]}; "
        f"right_total_vector={coefficient[1]}"
    )
    print(
        f"arithmetic=57068*k1=w:{FIBRE * K1};"
        f"57068*c=g0:{FIBRE * C};"
        f"each_copy=wc=g0*k1:{W * C};"
        f"right_total=2wc:{2 * W * C};"
        f"primitive=2k1:{2 * K1}"
    )
    print(
        f"intrinsic_Bockstein_mod13="
        f"(each_copy={coefficient[2]},sum={coefficient[3]}=9+9)"
    )
    print(
        "lawful_positive=after the inherited coefficient filters, R is the "
        "literal disjoint union of two equal coefficient/Bockstein summands"
    )
    print(
        "minimal_obstruction=no cospan from M to the copies in the typed "
        "physical category: native E3 already fails on both; Rplus also "
        "fails c2; source carrier-twist support is empty and endpoint masks "
        "are not translates"
    )
    print(
        "boundary=the additive Bockstein split survives, but it supplies "
        "neither THM-2771 target convolution nor a target-role/root-deck "
        "intertwiner; equal content is not a common-ancestry map"
    )


if __name__ == "__main__":
    main()
