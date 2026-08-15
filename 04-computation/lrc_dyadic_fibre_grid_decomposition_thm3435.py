"""Exact companion for the dyadic fibre-grid decomposition of half twists.

For Q=2^a N with N odd, write r=2^b s after truncating b at a.  If
b=a, the literal half-twist block is a full 2^a-fold pullback from N.  If
b<a, remove the b inactive levels and put M=2^(a-b)=7q+c.  On every
remaining M-fibre the odd coefficient has q points, plus one exactly over
the smaller-radius block H_N^(c/14)(s+qN).  The local points are consecutive
in the s-coordinate on Z/MZ.

The same phase-grid calculation gives a sharp Boolean branch boundary.
Two and four coefficient preimages have pairwise disjoint masks (OR=XOR),
whereas the eight-preimage family covers every sheet with multiplicity one
or two.  On each sheet its possible double hit is a neighbouring pair in a
sheet-dependent cyclic order; globally only opposite-parity labels can meet.
"""

from __future__ import annotations

import ast
from collections import Counter
from hashlib import sha256
from math import gcd
from pathlib import Path


EXPECTED_SEMANTIC_SHA256 = "bb9fce50f5ee51048eb02a79aef83151a1593b8b16ef73992243b88d8726545e"


def require(condition: bool, payload) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    data = path.read_bytes().replace(b"\r\n", b"\n")
    return sha256(data).hexdigest()


def v2(value: int) -> int:
    if value == 0:
        return 10**9
    answer = 0
    while value % 2 == 0:
        value //= 2
        answer += 1
    return answer


def odd_part_of_modulus(modulus: int) -> tuple[int, int]:
    exponent = v2(modulus)
    return exponent, modulus >> exponent


def danger_mask(modulus: int, residue: int, radius_numerator: int = 1) -> int:
    """H_modulus^(radius_numerator/14)(residue), with strict endpoints."""
    require(modulus >= 1 and 0 <= radius_numerator <= 7, (modulus, radius_numerator))
    period = 2 * modulus
    residue %= period
    answer = 0
    for sheet in range(modulus):
        phase = residue * (2 * sheet + 1) % period
        distance = min(phase, period - phase)
        if 14 * distance < 2 * modulus * radius_numerator:
            answer |= 1 << sheet
    return answer


def dyadic_data(modulus: int, residue: int):
    exponent, odd_base = odd_part_of_modulus(modulus)
    residue %= 2 * modulus
    inactive_depth = min(v2(residue), exponent)
    if inactive_depth == exponent:
        return {
            "a": exponent,
            "N": odd_base,
            "b": inactive_depth,
            "inactive": True,
            "s": residue >> exponent,
        }
    odd_residue = residue >> inactive_depth
    active_size = 1 << (exponent - inactive_depth)
    quotient, remainder = divmod(active_size, 7)
    require(odd_residue % 2 == 1 and remainder in (1, 2, 4), (modulus, residue))
    return {
        "a": exponent,
        "N": odd_base,
        "b": inactive_depth,
        "inactive": False,
        "s": odd_residue,
        "M": active_size,
        "q": quotient,
        "c": remainder,
        "extra_residue": (odd_residue + quotient * odd_base) % (2 * odd_base),
    }


def local_active_n_set(odd_base: int, active_size: int, odd_residue: int, base: int) -> tuple[int, ...]:
    """Occupied points after the coordinate change n=s*t on Z/MZ."""
    period = 2 * odd_base * active_size
    answer = []
    for coordinate in range(active_size):
        numerator = odd_residue * (2 * base + 1) + 2 * odd_base * coordinate
        phase = numerator % period
        distance = min(phase, period - phase)
        if 14 * distance < 2 * odd_base * active_size:
            answer.append(coordinate)
    return tuple(answer)


def is_cyclic_interval(values: tuple[int, ...], modulus: int) -> bool:
    value_set = set(values)
    if not values or len(values) == modulus:
        return True
    length = len(values)
    return any(
        value_set == {(start + offset) % modulus for offset in range(length)}
        for start in values
    )


def predicted_mask(modulus: int, residue: int) -> int:
    data = dyadic_data(modulus, residue)
    odd_base = data["N"]
    if data["inactive"]:
        downstairs = danger_mask(odd_base, data["s"])
        return sum(
            1 << (base + odd_base * lift)
            for base in range(odd_base)
            if downstairs >> base & 1
            for lift in range(1 << data["a"])
        )

    active_size = data["M"]
    inactive_multiplicity = 1 << data["b"]
    inverse = pow(data["s"], -1, active_size)
    answer = 0
    for base in range(odd_base):
        for coordinate in local_active_n_set(odd_base, active_size, data["s"], base):
            reduced_lift = inverse * coordinate % active_size
            for copy in range(inactive_multiplicity):
                lift = reduced_lift + active_size * copy
                answer |= 1 << (base + odd_base * lift)
    return answer


def formula_audit(limit: int = 192):
    sheet_cells = 0
    fibre_rows = 0
    active_rows = 0
    inactive_rows = 0
    inactive_endpoint_checks = 0
    remainder_histogram = Counter()
    parity_cycle = Counter()
    active_endpoint_hits = 0
    inactive_endpoint_example = None

    for modulus in range(1, limit + 1):
        exponent, odd_base = odd_part_of_modulus(modulus)
        for residue in range(2 * modulus):
            direct = danger_mask(modulus, residue)
            predicted = predicted_mask(modulus, residue)
            require(direct == predicted, ("mask", modulus, residue, direct, predicted))
            sheet_cells += modulus
            data = dyadic_data(modulus, residue)

            endpoint_sheets = tuple(
                sheet
                for sheet in range(modulus)
                if 14 * min(
                    residue * (2 * sheet + 1) % (2 * modulus),
                    2 * modulus - residue * (2 * sheet + 1) % (2 * modulus),
                )
                == 2 * modulus
            )

            if data["inactive"]:
                inactive_rows += 1
                divisor = gcd(odd_base, data["s"])
                quotient_order = odd_base // divisor
                normalized = data["s"] // divisor
                endpoint_expected = quotient_order % 7 == 0 and normalized % 2 == 1
                expected_endpoint_count = (1 << (exponent + 1)) * divisor if endpoint_expected else 0
                require(len(endpoint_sheets) == expected_endpoint_count,
                        ("inactive_endpoint_law", modulus, residue, endpoint_sheets, data,
                         divisor, quotient_order, normalized, expected_endpoint_count))
                inactive_endpoint_checks += 1
                if modulus % 2 == 0 and endpoint_sheets and inactive_endpoint_example is None:
                    inactive_endpoint_example = (modulus, residue, endpoint_sheets[0])
                continue

            active_rows += 1
            active_endpoint_hits += len(endpoint_sheets)
            extra = danger_mask(odd_base, data["extra_residue"], data["c"])
            counts = []
            for base in range(odd_base):
                local = local_active_n_set(odd_base, data["M"], data["s"], base)
                require(is_cyclic_interval(local, data["M"]), ("interval", modulus, residue, base, local))
                expected_reduced = data["q"] + ((extra >> base) & 1)
                require(len(local) == expected_reduced, ("count", modulus, residue, base, local, data))
                full_count = sum(
                    bool(direct >> (base + odd_base * lift) & 1)
                    for lift in range(1 << exponent)
                )
                require(full_count == (1 << data["b"]) * expected_reduced,
                        ("full_count", modulus, residue, base, full_count, data))
                counts.append(full_count)
                fibre_rows += 1
            expected_total = (1 << data["b"]) * (data["q"] * odd_base + extra.bit_count())
            require(direct.bit_count() == expected_total,
                    ("total", modulus, residue, direct.bit_count(), expected_total, data))
            remainder_histogram[(data["M"], data["q"], data["c"])] += odd_base
            parity_cycle[(data["M"].bit_length() - 1) % 3, data["c"], data["q"] % 2] += 1

    require(active_endpoint_hits == 0, active_endpoint_hits)
    require(inactive_endpoint_example == (14, 2, 0), inactive_endpoint_example)
    return (
        limit,
        sheet_cells,
        fibre_rows,
        active_rows,
        inactive_rows,
        inactive_endpoint_checks,
        tuple(sorted(remainder_histogram.items())),
        tuple(sorted(parity_cycle)),
        active_endpoint_hits,
        inactive_endpoint_example,
    )


def branch_family(modulus: int, seed: int, depth: int) -> tuple[int, ...]:
    require(depth >= 1 and modulus % (1 << (depth - 1)) == 0, (modulus, depth))
    step = modulus // (1 << (depth - 1))
    return tuple((seed + branch * step) % (2 * modulus) for branch in range(1 << depth))


def branch_audit(limit: int = 192):
    rows = Counter()
    first_seven_fail = None
    first_seven_cover = None
    first_depth_two = None
    first_depth_three = None

    for depth in (1, 2, 3):
        divisor = 1 << (depth - 1)
        for modulus in range(divisor, limit + 1, divisor):
            step = modulus // divisor
            for seed in range(step):
                residues = branch_family(modulus, seed, depth)
                masks = tuple(danger_mask(modulus, residue) for residue in residues)
                multiplicities = tuple(
                    sum(bool(mask >> sheet & 1) for mask in masks)
                    for sheet in range(modulus)
                )
                joined = 0
                for mask in masks:
                    joined |= mask
                rows[(depth, min(multiplicities), max(multiplicities))] += 1

                if depth <= 2:
                    require(max(multiplicities) <= 1, ("xor", modulus, seed, depth, multiplicities))
                    fused = danger_mask(modulus, (1 << depth) * seed, 1 << depth)
                    require(joined == fused, ("fusion", modulus, seed, depth, joined, fused))
                    if depth == 2 and first_depth_two is None and joined:
                        first_depth_two = (modulus, seed, residues, joined.bit_count())
                    continue

                full = (1 << modulus) - 1
                require(joined == full, ("depth3_cover", modulus, seed, joined.bit_count()))
                require(min(multiplicities) >= 1 and max(multiplicities) <= 2,
                        ("depth3_multiplicity", modulus, seed, multiplicities))
                for left in range(8):
                    for right in range(left + 1, 8):
                        # On sheet ell the branch order is multiplied by the
                        # odd unit 2ell+1.  A double hit is adjacent in that
                        # local order.  After forgetting ell this permits any
                        # odd label difference, but never an even one.
                        if (right - left) % 2 == 0:
                            require(masks[left] & masks[right] == 0,
                                    ("same_parity_overlap", modulus, seed, left, right))
                if first_depth_three is None:
                    first_depth_three = (modulus, seed, residues, Counter(multiplicities))
                transverse_nonempty = all(
                    residue % modulus != 0 and mask != 0
                    for residue, mask in zip(residues, masks)
                )
                for dropped in range(8):
                    seven_union = 0
                    for index, mask in enumerate(masks):
                        if index != dropped:
                            seven_union |= mask
                    row = (modulus, seed, dropped, seven_union.bit_count())
                    if transverse_nonempty and seven_union != full and first_seven_fail is None:
                        first_seven_fail = row
                    if transverse_nonempty and seven_union == full and first_seven_cover is None:
                        first_seven_cover = row

    require(first_depth_two is not None and first_depth_three is not None,
            (first_depth_two, first_depth_three))
    require(first_seven_fail is not None and first_seven_cover is not None,
            (first_seven_fail, first_seven_cover))
    return (
        limit,
        tuple(sorted(rows.items())),
        first_depth_two,
        first_depth_three,
        first_seven_fail,
        first_seven_cover,
    )


def hostile_audit():
    # Odd active masks are sections only at dyadic depths one and two.
    section_failure = None
    for modulus in range(2, 100, 2):
        exponent, odd_base = odd_part_of_modulus(modulus)
        for residue in range(1, 2 * modulus, 2):
            counts = tuple(
                sum(bool(danger_mask(modulus, residue) >> (base + odd_base * lift) & 1)
                    for lift in range(1 << exponent))
                for base in range(odd_base)
            )
            if max(counts, default=0) >= 2:
                section_failure = (modulus, residue, counts)
                break
        if section_failure is not None:
            break
    require(section_failure == (8, 1, (2,)), section_failure)

    # Counts alone do not locate the affine lift.
    affine_loss = None
    for modulus in range(2, 100, 2):
        profiles = {}
        exponent, odd_base = odd_part_of_modulus(modulus)
        for residue in range(1, 2 * modulus, 2):
            mask = danger_mask(modulus, residue)
            profile = tuple(
                sum(bool(mask >> (base + odd_base * lift) & 1)
                    for lift in range(1 << exponent))
                for base in range(odd_base)
            )
            if profile in profiles and profiles[profile][0] != mask:
                affine_loss = (modulus, profiles[profile][1], residue, profile,
                               profiles[profile][0], mask)
                break
            profiles[profile] = (mask, residue)
        if affine_loss is not None:
            break
    require(affine_loss is not None, "missing affine-loss hostile")

    # The active branch has no strict endpoint, while inactive pullbacks do.
    strict_inactive = (14, 2, tuple(
        sheet for sheet in range(14)
        if 14 * min(2 * (2 * sheet + 1) % 28, 28 - 2 * (2 * sheet + 1) % 28) == 28
    ))
    require(strict_inactive == (14, 2, (0, 6, 7, 13)), strict_inactive)
    return section_failure, affine_loss, strict_inactive


def main() -> None:
    source = Path(__file__)
    tree = ast.parse(source.read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "assert found")

    formula = formula_audit()
    branches = branch_audit()
    hostiles = hostile_audit()
    semantic_surface = (formula, branches, hostiles)
    semantic_digest = sha256(repr(semantic_surface).encode("ascii")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_digest == EXPECTED_SEMANTIC_SHA256,
                (semantic_digest, EXPECTED_SEMANTIC_SHA256))

    print("THM-3435 dyadic fibre-grid decomposition exact companion")
    print("status=VERIFIED_EXACT_COMPANION_FOR_PROVED_THM3435;all_modulus_literal_half_twist_identity;no_rank7_or_LRC14_conclusion")
    print("decomposition=Q=2^a*N_Nodd;r=2^b*s_with_b=min(v2(r),a);b=a_gives_full_2^a_pullback;"
          "b<a_gives_M=2^(a-b)=7q+c_and_each_N_fibre_has_2^b*(q+extra_y)_points")
    print("extra_mask=extra_y_is_H_N^(c/14)(s+qN);c_cycles_1,2,4_with_active_depth_mod3;"
          "local_points_are_a_cyclic_interval_in_the_n=s*t_coordinate")
    print("strict_boundary=active_dyadic_branches_never_hit_1/14_endpoints_by_2_adic_valuation;"
          "inactive_pullbacks_can_hit_and_retain_strict_exclusion")
    print("boolean_branch=depth1_and_depth2_preimage_masks_are_pairwise_disjoint_and_fuse_to_radius_2^k/14;"
          "depth3_eight_branch_family_covers_every_sheet_with_multiplicity_1_or_2;each_sheet_has_a_local_C8_order_but_global_overlap_is_only_across_the_4plus4_parity_bipartition")
    print(f"formula_audit=(limit,sheet_cells,fibre_rows,active_rows,inactive_rows,inactive_endpoint_checks,remainder_histogram,parity_cycle,active_endpoint_hits,inactive_endpoint_example)={formula}")
    print(f"branch_audit=(limit,multiplicity_histogram,first_depth2,first_depth3,first_seven_fail,first_seven_cover)={branches}")
    print(f"hostiles=(first_not_a_section,first_equal_count_different_affine_lift,inactive_strict_endpoint)={hostiles}")
    print("information_loss=fibre_counts_and_extra_mask_do_not_locate_the_affine_cyclic_interval;"
          "a_tournament_orientation_is_not_intrinsic;depth2_has_missing_edges_and_depth3_global_conflict_graph_is_a_subgraph_of_K4,4_with_a_sheet-dependent_local_C8_sidecar")
    print("scope=literal_half_twist_masks_only;no_owner_cap_classification;no_arbitrary_common_time;no_physical_current;LRC14_remains_open")
    print(f"semantic_sha256={semantic_digest}")
    print(f"script_sha256_lf={lf_sha256(source)}")
    print("reproducibility=standard_library_only;no_elapsed_fields;all_truth_gates_survive_python_O")
    print("commands=python -B 04-computation/lrc_dyadic_fibre_grid_decomposition_thm3435.py;"
          "python -B -O 04-computation/lrc_dyadic_fibre_grid_decomposition_thm3435.py")


if __name__ == "__main__":
    main()
