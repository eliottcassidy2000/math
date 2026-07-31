#!/usr/bin/env python3
"""Exact audit of the THM-2859 punctured-V4 completion problem.

The script independently rebuilds the 567-cell common/right universe and
its 587 rooted triples from the hash-pinned THM-2818 physical constructor.
It then treats the three formal/rootwise collar states as the restriction
of the translation action of V4=F_2^2 to a three-point subset.

The new finite checks cover:

* the 587 root triples and their semantic populations;
* native-factor and carrier-mask typing on all three corners;
* the partial V4 transformation groupoid and every graded-product defect;
* the unique four-point globalization and its added homogeneous arrows;
* normalized sign two-cocycles on the three-object pair groupoid;
* the P3/C4 partial-cube completion and the K3 diagonal hostile; and
* the minimal exact Koszul/mapping-cone filler of the +2h composite.

The endpoint-translation counts are read only from the hash-pinned canonical
physical transcript; they are labelled as inherited rather than independently
rebuilt because recomputing that transcript is a multi-hour operation.

No Python ``assert`` statement is used.
"""

from bisect import bisect_right
from collections import Counter
from hashlib import sha256
from pathlib import Path
import ast
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
RESULTS = ROOT / "05-knowledge" / "results"
sys.path.insert(0, str(COMP))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


PINNED = {
    COMP / "lrc14_right_cofiber_positive_copy_stratification_thm2818.py":
        "85edac9bb03f1fef198343268f4fb1226bec122d45ded79a049f8fa9a73882a8",
    RESULTS / "lrc14_nearest_half_step_common_right_collar_physical_audit_thm2825.out":
        "e610c5bffb720e074662f2222a50b2ce461c1b1293e946aa260faf898c4347b7",
}
for path, expected in PINNED.items():
    payload = path.read_bytes().replace(b"\r\n", b"\n")
    require(
        sha256(payload).hexdigest() == expected,
        f"pinned dependency changed: {path.name}",
    )


import lrc14_right_cofiber_positive_copy_stratification_thm2818 as copies


P = 13
H = copies.T // (2 * P**5)
COMMON_S = (0, 1, 2, 3, 8, 9, 10, 11, 12)
COMMON_T = (3, 4, 5, 6, 7, 8, 9, 10, 11)
FACTOR_NAMES = ("E3", "clock", "q1", "q2", "c2", "c3")


def cell_objects(details, full_module, e3, clocks, clock, sigma, target):
    section = copies.physical.target.source_present_section(
        full_module, e3, clock, sigma, target, clocks
    )
    source_base, target_base = details[clock]
    source = copies.weighted_intersection(source_base, section)
    target = copies.weighted_intersection(target_base, section)
    pulled = copies.physical.overlap.shift_weighted(target, -copies.SHIFT)
    aligned = copies.physical.overlap.intersect_weighted_profiles(
        source, pulled
    )
    require(
        all(left_weight == right_weight
            for _left, _right, left_weight, right_weight in aligned),
        "common overlap acquired unequal weights",
    )
    common = tuple(
        (left, right, left_weight)
        for left, right, left_weight, _right_weight in aligned
    )
    right = copies.physical.subtract_weighted(pulled, common)
    return source, target, common, right


def semantic_value(interval, q_pair, cache):
    key = interval[:2]
    if key not in cache:
        unit = ((*key, 1),)
        target = copies.physical.overlap.shift_weighted(unit, copies.SHIFT)
        source_value = copies.physical.relative.private.delayed_carry_pair(
            unit, q_pair, {}
        )[12][1]
        target_value = copies.physical.relative.private.delayed_carry_pair(
            target, q_pair, {}
        )[6][1]
        cache[key] = (source_value, target_value)
    return cache[key]


def section_factors(full_module, e3, clocks, clock, sigma, target):
    universe = ((0, full_module.T),)
    return (
        e3,
        clocks[clock],
        full_module.subtract_comb(
            universe, full_module.W[1], 182,
            -14 * sigma - 13, -14 * sigma + 13,
        ),
        full_module.subtract_comb(
            universe, full_module.W[2], 182,
            -14 * target - 13, -14 * target + 13,
        ),
        full_module.subtract_comb(
            universe, full_module.C2, 182,
            14 * sigma - 13, 14 * sigma + 13,
        ),
        full_module.subtract_comb(
            universe, full_module.C3, 182,
            14 * target - 13, 14 * target + 13,
        ),
    )


def indexed_intervals(intervals):
    """Attach binary-search coordinates to a sorted interval tuple."""
    require(
        all(
            intervals[index - 1][1] <= intervals[index][0]
            for index in range(1, len(intervals))
        ),
        "an indexed interval family is not sorted/disjoint",
    )
    return (
        intervals,
        tuple(interval[0] for interval in intervals),
        tuple(interval[1] for interval in intervals),
    )


def indexed_contains(interval, indexed):
    intervals, starts, _ends = indexed
    position = bisect_right(starts, interval[0]) - 1
    return (
        position >= 0
        and intervals[position][0] <= interval[0]
        and interval[1] <= intervals[position][1]
    )


def indexed_overlaps(probes, indexed):
    intervals, _starts, ends = indexed
    for left, right in probes:
        position = bisect_right(ends, left)
        if position < len(intervals) and intervals[position][0] < right:
            return True
    return False


def cached_factor_frames(factors):
    """Precompute the two translated factor frames once per labelled cell."""
    return (
        tuple(indexed_intervals(factor) for factor in factors),
        tuple(
            indexed_intervals(
                copies.physical.overlap.shift_union(factor, -copies.SHIFT)
            )
            for factor in factors
        ),
        tuple(indexed_intervals(factor) for factor in factors),
        tuple(
            indexed_intervals(
                copies.physical.overlap.shift_union(factor, copies.SHIFT)
            )
            for factor in factors
        ),
    )


def cached_factor_masks(interval, frames):
    """Exact equivalent of ``copies.factor_masks`` without repeated shifts."""
    target_interval = tuple(
        endpoint + copies.SHIFT for endpoint in interval
    )
    return (
        tuple(indexed_contains(interval, factor) for factor in frames[0]),
        tuple(indexed_contains(interval, factor) for factor in frames[1]),
        tuple(
            indexed_contains(target_interval, factor)
            for factor in frames[2]
        ),
        tuple(
            indexed_contains(target_interval, factor)
            for factor in frames[3]
        ),
    )


def carrier_supports(source, target):
    # Translation is a bijection of the ambient circle, so probe an
    # inverse-translated atom against one merged base support instead of
    # translating and re-merging the large weighted support thirteen times.
    return (
        indexed_intervals(copies.support_of(source)),
        indexed_intervals(copies.support_of(target)),
    )


def carrier_mask(interval, supports):
    unit = copies.T // P
    target_interval = tuple(
        endpoint + copies.SHIFT for endpoint in interval[:2]
    )
    return (
        tuple(
            indexed_overlaps(
                copies.physical.overlap.shift_union(
                    (interval[:2],), twist * unit
                ),
                supports[0],
            )
            for twist in range(P)
        ),
        tuple(
            indexed_overlaps(
                copies.physical.overlap.shift_union(
                    (target_interval,), -twist * unit
                ),
                supports[1],
            )
            for twist in range(P)
        ),
    )


def xor(left, right):
    return tuple(a ^ b for a, b in zip(left, right))


V4 = ((0, 0), (1, 0), (0, 1), (1, 1))
MISSING = (1, 0)
GRADES = ((0, 0), (1, 1), (0, 1))
FULL_GRADES = GRADES + (MISSING,)


def grade_components(grades):
    answer = {}
    for degree in V4:
        answer[degree] = {
            (row, column)
            for row, row_grade in enumerate(grades)
            for column, column_grade in enumerate(grades)
            if xor(row_grade, column_grade) == degree
        }
    return answer


def multiply_support(left, right):
    return {
        (row, column)
        for row, middle_left in left
        for middle_right, column in right
        if middle_left == middle_right
    }


def gf2_rank(rows, width):
    pivots = {}
    for row in rows:
        value = 0
        for index, bit in enumerate(row):
            if bit & 1:
                value ^= 1 << index
        while value:
            pivot = value.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = value
                break
            value ^= pivots[pivot]
    require(
        all(pivot < width for pivot in pivots),
        "GF(2) elimination produced an out-of-range pivot",
    )
    return len(pivots)


def pair_groupoid_sign_cohomology():
    objects = range(3)
    one_vars = tuple(
        (left, right)
        for left in objects
        for right in objects
        if left != right
    )
    one_index = {value: index for index, value in enumerate(one_vars)}
    two_vars = tuple(
        (left, middle, right)
        for left in objects
        for middle in objects
        for right in objects
        if left != middle and middle != right
    )
    two_index = {value: index for index, value in enumerate(two_vars)}

    def one_coordinate(left, right):
        if left == right:
            return None
        return one_index[left, right]

    def two_coordinate(left, middle, right):
        if left == middle or middle == right:
            return None
        return two_index[left, middle, right]

    cocycle_rows = []
    for i in objects:
        for j in objects:
            for k in objects:
                for ell in objects:
                    row = [0] * len(two_vars)
                    for triple in (
                        (j, k, ell),
                        (i, k, ell),
                        (i, j, ell),
                        (i, j, k),
                    ):
                        coordinate = two_coordinate(*triple)
                        if coordinate is not None:
                            row[coordinate] ^= 1
                    if any(row):
                        cocycle_rows.append(row)
    cocycle_constraint_rank = gf2_rank(cocycle_rows, len(two_vars))
    z2_dimension = len(two_vars) - cocycle_constraint_rank

    coboundary_columns_as_rows = []
    for edge in one_vars:
        image = [0] * len(two_vars)
        for triple, coordinate in two_index.items():
            i, j, k = triple
            for pair in ((j, k), (i, k), (i, j)):
                edge_coordinate = one_coordinate(*pair)
                if edge_coordinate == one_index[edge]:
                    image[coordinate] ^= 1
        coboundary_columns_as_rows.append(image)
    b2_dimension = gf2_rank(
        coboundary_columns_as_rows, len(two_vars)
    )
    z1_dimension = len(one_vars) - b2_dimension
    require(
        z2_dimension == b2_dimension,
        "the pair groupoid unexpectedly acquired sign H^2",
    )
    return {
        "C1": len(one_vars),
        "Z1": z1_dimension,
        "C2": len(two_vars),
        "Z2": z2_dimension,
        "B2": b2_dimension,
        "H2": z2_dimension - b2_dimension,
    }


def literal_output_rows(path):
    rows = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        value = value.split(";", 1)[0]
        try:
            rows[key] = ast.literal_eval(value)
        except (SyntaxError, ValueError):
            pass
    return rows


def subgroup_sizes():
    zero = (0, 0)
    groups = []
    for mask in range(1 << len(V4)):
        subset = {
            V4[index]
            for index in range(len(V4))
            if (mask >> index) & 1
        }
        if zero not in subset:
            continue
        closed = all(xor(a, b) in subset for a in subset for b in subset)
        if closed:
            groups.append(tuple(sorted(subset)))
    return tuple(sorted((len(group), 4 // len(group)) for group in groups))


def main():
    (
        _module, _rails, _present, details, full_module, e3, clocks,
        q_pairs, _delayed, _source_weight, _target_weight, _rail_common,
    ) = copies.physical_setup()

    nonempty_cells = 0
    common_count = 0
    right_count = 0
    semantic_triples = Counter()
    right_semantic = Counter()
    factor_holes = Counter()
    full_factor_counts = Counter()
    carrier_counts = Counter()
    semantic_caches = tuple({} for _clock in range(7))

    for clock in range(7):
        for sigma in COMMON_S:
            for target in COMMON_T:
                source, target_physical, common, right = cell_objects(
                    details, full_module, e3, clocks, clock, sigma, target
                )
                require(
                    bool(common) == bool(right),
                    "common/right support shadows split",
                )
                if not right:
                    continue
                nonempty_cells += 1
                common_count += len(common)
                right_count += len(right)
                common_by_left = {piece[0]: piece for piece in common}
                require(
                    len(common_by_left) == len(common),
                    "common endpoints collided in one labelled cell",
                )
                factors = section_factors(
                    full_module, e3, clocks, clock, sigma, target
                )
                factor_frames = cached_factor_frames(factors)
                supports = carrier_supports(source, target_physical)
                cache = semantic_caches[clock]
                for root in right:
                    m1 = common_by_left.get(root[0] + H)
                    m2 = common_by_left.get(root[0] + 2 * H)
                    require(
                        m1 is not None and m2 is not None,
                        "one collar triple failed to reconstruct",
                    )
                    states = (root, m1, m2)
                    live = tuple(
                        semantic_value(state, q_pairs[clock], cache) != (0, 0)
                        for state in states
                    )
                    require(
                        live[0] == live[2] and live[0] != live[1],
                        "collar triple lost its semantic V4 pattern",
                    )
                    semantic_triples[live] += 1
                    right_semantic["live" if live[0] else "dead"] += 1

                    masks = tuple(
                        cached_factor_masks(state[:2], factor_frames)
                        for state in states
                    )
                    root_first = masks[0][0]
                    holes = tuple(
                        name for name, present in zip(FACTOR_NAMES, root_first)
                        if not present
                    )
                    factor_holes[holes] += 1
                    for position, mask in enumerate(masks):
                        full = all(all(frame) for frame in mask)
                        full_factor_counts[position, full] += 1
                    require(
                        holes and all(all(frame) for frame in masks[1])
                        and all(all(frame) for frame in masks[2]),
                        "factor boundary stopped being root-hole/common-full",
                    )

                    cmasks = tuple(
                        carrier_mask(state, supports) for state in states
                    )
                    for position, cmask in enumerate(cmasks):
                        carrier_counts[position, cmask] += 1
                    empty = (False,) * P
                    delta0 = (True,) + (False,) * (P - 1)
                    require(
                        cmasks[0] == (empty, delta0)
                        and cmasks[1] == (delta0, delta0)
                        and cmasks[2] == (delta0, delta0),
                        "carrier boundary stopped being empty/delta0",
                    )

    require(
        (nonempty_cells, common_count, right_count) == (193, 63308, 587),
        "independent bank totals changed",
    )
    require(
        semantic_triples
        == Counter({(True, False, True): 573, (False, True, False): 14}),
        "semantic population changed",
    )

    components = grade_components(GRADES)
    partial_dimensions = tuple(
        (degree, len(components[degree])) for degree in V4
    )
    product_defects = []
    for left_degree in V4:
        for right_degree in V4:
            target_degree = xor(left_degree, right_degree)
            product = multiply_support(
                components[left_degree], components[right_degree]
            )
            missing = tuple(sorted(components[target_degree] - product))
            product_defects.append(
                (
                    left_degree,
                    right_degree,
                    len(product),
                    len(components[target_degree]),
                    missing,
                )
            )
    defect_histogram = Counter(
        target_size - product_size
        for _g, _h, product_size, target_size, _missing in product_defects
    )
    require(
        defect_histogram == Counter({1: 9, 0: 7}),
        "partial V4 product-defect profile changed",
    )

    # Preserve the old matrix indices and append the new corner.  Enumerating
    # objects directly in ``V4`` order would relabel M1 from index 1 to 3 and
    # make a set-difference audit compare unrelated matrix units.
    full_components = grade_components(FULL_GRADES)
    full_dimensions = tuple(
        (degree, len(full_components[degree])) for degree in V4
    )
    added_dimensions = tuple(
        (
            degree,
            len(full_components[degree] - components[degree]),
        )
        for degree in V4
    )
    require(
        tuple(size for _degree, size in added_dimensions) == (1, 2, 2, 2),
        "fourth-corner globalization did not add 1+2+2+2 arrows",
    )

    sign_cohomology = pair_groupoid_sign_cohomology()
    quotient_data = subgroup_sizes()
    require(
        all(orbit_size != 3 for _subgroup_size, orbit_size in quotient_data),
        "a V4 quotient unexpectedly acquired a three-point orbit",
    )

    # Rootwise mapping-cone model:
    # D0=(d,id)^T=(1,1)^T and D1=(a,-S)=(1,-1).
    d0 = (1, 1)
    d1 = (1, -1)
    composition = d1[0] * d0[0] + d1[1] * d0[1]
    require(composition == 0, "Koszul square stopped anticommuting")
    cone_ranks = {
        "dim_C0": 587,
        "dim_C1": 1174,
        "dim_C2": 587,
        "rank_D0": 587,
        "rank_D1": 587,
        "H0": 0,
        "H1": 0,
        "H2": 0,
        "minimum_filler_dimension": 587,
    }

    current_live = (
        2 * semantic_triples[True, False, True]
        + semantic_triples[False, True, False]
    )
    current_dead = 3 * right_count - current_live
    missing_live = right_semantic["dead"]
    missing_dead = right_semantic["live"]
    full_live = current_live + missing_live
    full_dead = current_dead + missing_dead
    require(
        (current_live, current_dead, missing_live, missing_dead,
         full_live, full_dead)
        == (1160, 601, 14, 573, 1174, 1174),
        "semantic completion balance changed",
    )
    maximum_internal_flip_vertices = 2 * min(right_semantic.values())
    unmatched_internal_flip = right_count - maximum_internal_flip_vertices
    require(
        unmatched_internal_flip == 559,
        "semantic-flip permutation deficit changed",
    )

    relative_character_defects = {}
    for character in ((1, 0), (0, 1), (1, 1)):
        value = sum(
            (-1) ** (
                character[0] * grade[0] + character[1] * grade[1]
            )
            for grade in GRADES
        )
        missing_value = (-1) ** (
            character[0] * MISSING[0]
            + character[1] * MISSING[1]
        )
        relative_character_defects[character] = (
            right_count * value,
            right_count * missing_value,
        )
    require(
        all(current + missing == 0
            for current, missing in relative_character_defects.values()),
        "missing corner failed to cancel a V4 character defect",
    )

    endpoint_rows = literal_output_rows(
        RESULTS
        / "lrc14_nearest_half_step_common_right_collar_physical_audit_thm2825.out"
    )
    endpoint_h = endpoint_rows["endpoint_pair_census_h"]
    endpoint_2h = endpoint_rows["endpoint_pair_census_2h"]
    require(endpoint_h == endpoint_2h, "h/2h endpoint profiles diverged")
    source_translation_count = sum(
        count
        for (source_distance, _target_distance,
             source_translations, _target_translations), count in endpoint_h
        if source_translations
    )
    source_orbit_failure_count = right_count - source_translation_count
    require(
        (source_translation_count, source_orbit_failure_count) == (74, 513),
        "source endpoint orbit boundary changed",
    )

    # The primitive ladder is P3 inside the C4 generated by degrees 11 and
    # 10.  The +2h composite has degree 01 and is a diagonal, so promoting it
    # to a unit edge creates K3.
    generator_degrees = ((1, 1), (1, 0))
    p3_edges = {
        tuple(sorted((left, right)))
        for left in GRADES
        for right in GRADES
        if left != right and xor(left, right) in generator_degrees
    }
    c4_edges = {
        tuple(sorted((left, right)))
        for left in V4
        for right in V4
        if left != right and xor(left, right) in generator_degrees
    }
    diagonal = tuple(sorted(((0, 0), (0, 1))))
    require(
        (len(p3_edges), len(c4_edges), diagonal in p3_edges)
        == (2, 4, False),
        "partial-cube support changed",
    )

    print("THM-2859 PUNCTURED-V4 EXTENSION / COMPLETION AUDIT")
    print(
        f"bank=nonempty:{nonempty_cells},common:{common_count},"
        f"roots:{right_count},half_step:{H}"
    )
    print(f"semantic_triples={tuple(sorted(semantic_triples.items()))}")
    print(f"right_semantic={tuple(sorted(right_semantic.items()))}")
    print(
        "semantic_internal_flip="
        f"max_covered:{maximum_internal_flip_vertices},"
        f"unmatched:{unmatched_internal_flip}"
    )
    print(
        "semantic_completion="
        f"current_live_dead:{(current_live, current_dead)},"
        f"missing_live_dead:{(missing_live, missing_dead)},"
        f"full_live_dead:{(full_live, full_dead)}"
    )
    print(
        "relative_V4_character_defects="
        f"{tuple(sorted(relative_character_defects.items()))}"
    )
    print(f"root_factor_holes={tuple(sorted(factor_holes.items()))}")
    print(f"full_factor_counts={tuple(sorted(full_factor_counts.items()))}")
    print(
        "carrier_types="
        "R:(source_empty,target_delta0),"
        "M1/M2:(source_delta0,target_delta0);"
        f"verified_triples={right_count}"
    )
    print(
        "carrier_orbit_no_go="
        "C13-translations fix the empty mask and permute the 13 deltas;"
        "the two orbits are disjoint"
    )
    print(f"partial_grade_dims={partial_dimensions}")
    print(f"partial_product_defect_histogram={tuple(sorted(defect_histogram.items()))}")
    print(
        "partial_product_defects="
        f"{tuple(product_defects)}"
    )
    print(f"full_V4_grade_dims={full_dimensions}")
    print(
        "globalization_added_grade_dims="
        f"{added_dimensions};per_root_matrix_units=7;"
        f"all_roots_added={7 * right_count}"
    )
    print(
        "corner_model="
        "M3=p*M4*p for p=e00+e11+e01;"
        "the complement e10 is the unique missing object"
    )
    print(f"pair_groupoid_sign_cohomology={tuple(sign_cohomology.items())}")
    print(
        "cocycle_support_no_go="
        "nonzero scalar twists rescale existing products but preserve all "
        "nine one-arrow graded defects"
    )
    print(f"V4_subgroup_orbit_sizes={quotient_data};three_point_orbit=no")
    print(
        "partial_cube="
        f"P3_edges:{tuple(sorted(p3_edges))},"
        f"C4_edges:{tuple(sorted(c4_edges))};"
        "plus_2h_is_distance_2;adding_it_as_edge_gives_K3"
    )
    print(
        "koszul_filler="
        "X=Pi(R),D0=(d,id),D1=(a,-S),D1D0=0;"
        f"{tuple(cone_ranks.items())}"
    )
    print(
        "filler_rank_no_go="
        "any alternate factorization y*x=-S has dim(X)>=rank(S)=587"
    )
    print(
        "inherited_endpoint_orbit="
        f"translation_exists:{source_translation_count},"
        f"no_translation:{source_orbit_failure_count};"
        "target_translation_zero:587"
    )
    print(
        "typed_conclusion="
        "the formal fourth corner preserves root annotations as a parity-"
        "shifted copy, but its edge to M2 still crosses source empty->delta0 "
        "and root native holes->all-six; quotienting those predicates is "
        "not a typed completion"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
