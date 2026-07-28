#!/usr/bin/env python3
"""Exact replay of the refuted THM-2751 clock-blind wing calculation.

This companion reconstructs the legacy clock-blind source and target sheets
and subtracts THM-2749's marked common section.  MISTAKE-313 proves that these
objects are not lawfully comparable: the legacy source helper omits the
source-one clock factor.  The program remains as an exact hostile replay of
the arithmetic that produced the false global interpretation.

The final C3 calculation is deliberately typed as a boundary: the displayed
internal chart strata are arm-blind.  A hypothetical external selector with
the same three numerical gains would be charged, but no such selector is
constructed here.

All arithmetic is integral/rational.  Dependency pins use LF-normalized
bytes, and there are no truth-bearing Python assertions.
"""

from hashlib import sha256
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))

DEPENDENCIES = {
    "lrc14_fully_marked_root_zero_target_profile_thm2749.py":
        "d67c852c52f88feaadb2fcaa0a9a07a212f2e47018040b455855df886200595e",
    "lrc14_root_zero_overlap_clutch_20260728.py":
        "e10fa7c9a5a238461ef422ea314dc334f7e65ec1787cf65d4e4bea12b96aefb8",
    "lrc14_semantic_root_zero_clutch_refinement_probe_20260728.py":
        "b3b623a4c016b1303ac3a74c72f9ae0bbd69cdb97a553f613143f0005a3dd286",
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_bytes(path):
    return path.read_bytes().replace(bytes((13, 10)), bytes((10,)))


for dependency, expected in DEPENDENCIES.items():
    actual = sha256(lf_bytes(COMP / dependency)).hexdigest()
    require(actual == expected,
            f"audited dependency changed: {dependency}: {actual} != {expected}")


import lrc14_fully_marked_root_zero_target_profile_thm2749 as marked
import lrc14_semantic_root_zero_clutch_refinement_probe_20260728 as natural


clutch = natural.clutch
private = marked.private
two = marked.two

P = 13
T = marked.T
SHIFT = marked.SHIFT
CONTENT = marked.CONTENT
RAIL = marked.RAIL
TARGET_LABEL = 3

COMMON_AMPLITUDE = 339633525654239542165440
NATURAL_SOURCE_AMPLITUDE = 1812281403506324508838080
NATURAL_TARGET_AMPLITUDE = 1826551436254490256030720
LEFT_WING_AMPLITUDE = 1472647877852084966672640
RIGHT_WING_AMPLITUDE = 1486917910600250713865280
SHEAR_AMPLITUDE = 14270032748165747192640


def canonical_weighted(pieces):
    """Canonicalize a disjoint weighted step profile.

    The upstream exact interval routines return disjoint pieces.  We check
    that invariant explicitly and merge only touching pieces of equal weight.
    """
    result = []
    for left, right, weight in sorted(pieces):
        require(left < right, "empty weighted interval")
        if weight == 0:
            continue
        if result:
            old_left, old_right, old_weight = result[-1]
            require(left >= old_right, "weighted profile pieces overlap")
            if left == old_right and weight == old_weight:
                result[-1] = (old_left, right, weight)
                continue
        result.append((left, right, weight))
    return tuple(result)


def support(pieces):
    return marked.merge(
        (left, right) for left, right, weight in pieces if weight
    )


def weighted_mass(pieces):
    return sum((right - left) * weight for left, right, weight in pieces)


def restrict_weighted(pieces, intervals):
    return canonical_weighted(marked.restrict_weighted(
        canonical_weighted(pieces), marked.merge(intervals)
    ))


def subtract_weighted(pieces, intervals):
    """Remove an interval union from a weighted profile on the T-circle."""
    return restrict_weighted(pieces, marked.complement(marked.merge(intervals)))


def disjoint_union_weighted(left, right):
    left = canonical_weighted(left)
    right = canonical_weighted(right)
    require(not marked.intersect(support(left), support(right)),
            "purported weighted disjoint union overlaps")
    return canonical_weighted(left + right)


def coefficient(pieces, prefix, carry):
    return private.delayed_carry_pair(pieces, prefix, {})[carry][1]


def c3_profile(values):
    """Normalized C3 Fourier profile over F13 with omega=3."""
    omega = 3
    inv3 = pow(3, -1, P)
    return tuple(
        inv3 * sum(value * pow(omega, -frequency * arm, P)
                   for arm, value in enumerate(values)) % P
        for frequency in range(3)
    )


def c3_charged(values):
    inv3 = pow(3, -1, P)
    mean = inv3 * sum(values) % P
    return tuple((value - mean) % P for value in values)


def matrix_product(left, right):
    return tuple(tuple(
        sum(left[i][k] * right[k][j] for k in range(len(right))) % P
        for j in range(len(right[0])))
        for i in range(len(left))
    )


def common_rows(module, delayed_present, source_weight, semantic_common, ell):
    """Reconstruct THM-2749's source and direct-target common pieces."""
    source_safe = marked.complement(delayed_present[ell, 7])
    target_safe = marked.complement(
        marked.shift_union(delayed_present[ell, 7], SHIFT)
    )
    row_common = marked.intersect(semantic_common, source_safe)
    row_common = marked.intersect(row_common, target_safe)

    source = marked.restrict_weighted(source_weight, row_common)
    source = private.old.intersect_weighted_comb(
        source, module.C3, 182, 169, 181
    )

    # marked.shift_union is a pullback convention, so the argument -SHIFT
    # pushes the source-coordinate common set into the forward target chart.
    direct_common = marked.shift_union(row_common, -SHIFT)
    target = marked.restrict_weighted(source_weight, direct_common)
    target = private.old.intersect_weighted_comb(
        target, module.C3, 182, 1, 13
    )
    return canonical_weighted(source), canonical_weighted(target)


def valuation(number, prime):
    result = 0
    while number and number % prime == 0:
        number //= prime
        result += 1
    return result


def main():
    require((P, T, SHIFT, CONTENT, RAIL, TARGET_LABEL)
            == (13, 297836897838480, 431933040, 26, 8, 3),
            "canonical clutch constants changed")

    module, _prefixes, _a, _b, rails, present, _starts = (
        marked.lift.m.core.build_carrier_data()
    )
    require(module.C3 == 742586 and len(rails) == 162,
            "canonical carrier changed")

    natural_prefixes = natural.build_q3_pair_prefixes(module)
    marked_prefix_bank = marked.marked_prefixes(
        module,
        private.build_pair_prefixes(module),
        two.deepest_fork(module),
    )
    require(all(natural_prefixes[ell][6] == marked_prefix_bank[ell]
                for ell in range(7)),
            "natural and two-sided marked prefixes differ")

    target_pullback = clutch.shift_weighted(rails[RAIL][3], -SHIFT)
    rail_pairs = clutch.intersect_weighted_profiles(
        rails[RAIL][3], target_pullback
    )
    require(rail_pairs and all(row[2] == row[3] for row in rail_pairs),
            "rail-eight source/target weights stopped agreeing")

    source_weight, _target_weight, rail_common = marked.rail_data(
        rails, RAIL
    )
    semantic_source = two.exclusive_source(module, TARGET_LABEL)
    clock_comb = tuple(
        module.make_comb(module.C1, 182, 26 * ell - 13, 26 * ell + 13)
        for ell in range(7)
    )
    static_common = marked.static_semantic_common(
        module, semantic_source, clock_comb, rail_common
    )
    semantic_common = marked.target_label_common(
        module, static_common, TARGET_LABEL
    )

    source_natural_vector = []
    target_natural_vector = []
    source_common_vector = []
    target_common_vector = []
    left_wing_vector = []
    right_wing_vector = []

    source_natural_masses = []
    target_natural_masses = []
    source_common_masses = []
    target_common_masses = []
    left_wing_masses = []
    right_wing_masses = []

    piece_counts = []
    for ell in range(7):
        source_natural, target_natural = clutch.restrict_to_relative_overlap(
            module, present, rail_pairs, ell
        )
        source_natural = canonical_weighted(natural.restrict_e3_and_sheet(
            source_natural, module, *natural.SHEET
        ))
        target_natural = canonical_weighted(natural.restrict_e3_and_sheet(
            target_natural, module, *natural.SHEET
        ))

        source_common, target_common = common_rows(
            module, present, source_weight, semantic_common, ell
        )

        # Verify the common carriers as actual weighted subprofiles of the
        # natural carriers, before constructing either wing.
        require(restrict_weighted(source_natural, support(source_common))
                == source_common,
                f"source common profile is not a natural subprofile at {ell}")
        require(restrict_weighted(target_natural, support(target_common))
                == target_common,
                f"target common profile is not a natural subprofile at {ell}")

        left_wing = subtract_weighted(
            source_natural, support(source_common)
        )
        right_wing = subtract_weighted(
            target_natural, support(target_common)
        )

        require(disjoint_union_weighted(source_common, left_wing)
                == source_natural,
                f"A=M disjoint-union L failed at clock {ell}")
        require(disjoint_union_weighted(target_common, right_wing)
                == target_natural,
                f"B=M disjoint-union R failed at clock {ell}")

        prefix = marked_prefix_bank[ell]
        source_natural_value = coefficient(source_natural, prefix, 12)
        target_natural_value = coefficient(target_natural, prefix, 6)
        source_common_value = coefficient(source_common, prefix, 12)
        target_common_value = coefficient(target_common, prefix, 6)
        left_wing_value = coefficient(left_wing, prefix, 12)
        right_wing_value = coefficient(right_wing, prefix, 6)

        require(source_natural_value
                == source_common_value + left_wing_value,
                f"source functional additivity failed at {ell}")
        require(target_natural_value
                == target_common_value + right_wing_value,
                f"target functional additivity failed at {ell}")

        source_natural_vector.append(source_natural_value)
        target_natural_vector.append(target_natural_value)
        source_common_vector.append(source_common_value)
        target_common_vector.append(target_common_value)
        left_wing_vector.append(left_wing_value)
        right_wing_vector.append(right_wing_value)

        source_natural_masses.append(weighted_mass(source_natural))
        target_natural_masses.append(weighted_mass(target_natural))
        source_common_masses.append(weighted_mass(source_common))
        target_common_masses.append(weighted_mass(target_common))
        left_wing_masses.append(weighted_mass(left_wing))
        right_wing_masses.append(weighted_mass(right_wing))
        piece_counts.append((
            len(source_natural), len(source_common), len(left_wing),
            len(target_natural), len(target_common), len(right_wing),
        ))

    source_natural_vector = tuple(source_natural_vector)
    target_natural_vector = tuple(target_natural_vector)
    source_common_vector = tuple(source_common_vector)
    target_common_vector = tuple(target_common_vector)
    left_wing_vector = tuple(left_wing_vector)
    right_wing_vector = tuple(right_wing_vector)

    expected = lambda amplitude: (0,) + (amplitude,) * 6
    require(source_natural_vector == expected(NATURAL_SOURCE_AMPLITUDE),
            "natural source vector changed")
    require(target_natural_vector == expected(NATURAL_TARGET_AMPLITUDE),
            "natural target vector changed")
    require(source_common_vector == target_common_vector
            == expected(COMMON_AMPLITUDE),
            "two-sided common vector changed")
    require(left_wing_vector == expected(LEFT_WING_AMPLITUDE),
            "left-wing vector changed")
    require(right_wing_vector == expected(RIGHT_WING_AMPLITUDE),
            "right-wing vector changed")

    require(tuple(source_natural_masses)
            == (929934280541992017600,) * 7,
            "natural source weighted mass changed")
    require(tuple(target_natural_masses)
            == (929934304688494607040,) * 7,
            "natural target weighted mass changed")
    require(all(natural_mass == common_mass + wing_mass
                for natural_mass, common_mass, wing_mass in zip(
                    source_natural_masses,
                    source_common_masses,
                    left_wing_masses,
                )), "source weighted masses are not additive")
    require(all(natural_mass == common_mass + wing_mass
                for natural_mass, common_mass, wing_mass in zip(
                    target_natural_masses,
                    target_common_masses,
                    right_wing_masses,
                )), "target weighted masses are not additive")

    source_common_unit = marked.unit_reduction(source_common_vector, 12)
    target_common_unit = marked.unit_reduction(target_common_vector, 1)
    left_wing_unit = marked.unit_reduction(left_wing_vector, 12)
    right_wing_unit = marked.unit_reduction(right_wing_vector, 1)
    source_natural_unit = marked.unit_reduction(source_natural_vector, 12)
    target_natural_unit = marked.unit_reduction(target_natural_vector, 1)
    units = (
        source_common_unit, target_common_unit,
        left_wing_unit, right_wing_unit,
        source_natural_unit, target_natural_unit,
    )
    expected_units = (
        ((9, 0, 0, 0, 0, 0), 1),
        ((4, 0, 0, 0, 0, 0), 1),
        ((9, 0, 0, 0, 0, 0), 1),
        ((5, 0, 0, 0, 0, 0), 12),
        ((5, 0, 0, 0, 0, 0), 12),
        ((9, 0, 0, 0, 0, 0), 1),
    )
    require(units == expected_units, "root-normalized unit table changed")

    common_gain = 4 * pow(9, -1, P) % P
    wing_gain = 5 * pow(9, -1, P) % P
    natural_gain = 9 * pow(5, -1, P) % P
    require((common_gain, wing_gain, natural_gain) == (12, 2, 7),
            "common/wing/natural normalized ratios changed")
    require((common_gain - wing_gain) % P == 10,
            "internal grading contrast changed")

    shear_from_natural = (
        target_natural_vector[1] - source_natural_vector[1]
    )
    shear_from_wings = right_wing_vector[1] - left_wing_vector[1]
    require(shear_from_natural == shear_from_wings == SHEAR_AMPLITUDE,
            "natural shear is not the wing boundary current")
    require(valuation(SHEAR_AMPLITUDE, P) == 1
            and SHEAR_AMPLITUDE // CONTENT % P == 12,
            "wing shear valuation/content residue changed")

    # External-arm boundary.  These three gains would charge both nontrivial
    # C3 modes if a physical arm selector supplied them.  The actual internal
    # M/L/R operator is arm-blind, so Q(I_arm tensor T)Pi vanishes instead.
    hypothetical_arm_gains = (common_gain, wing_gain, wing_gain)
    hypothetical_fourier = c3_profile(hypothetical_arm_gains)
    mean = hypothetical_fourier[0]
    q_profile = tuple((value - mean) % P
                      for value in hypothetical_arm_gains)
    require((mean, q_profile, hypothetical_fourier)
            == (1, (11, 1, 1), (1, 12, 12)),
            "hypothetical external C3 selector profile changed")

    inv3 = pow(3, -1, P)
    identity = tuple(tuple(int(i == j) for j in range(3))
                     for i in range(3))
    projection = tuple(tuple(inv3 for _j in range(3))
                       for _i in range(3))
    charged = tuple(tuple((identity[i][j] - projection[i][j]) % P
                          for j in range(3))
                    for i in range(3))
    require(matrix_product(matrix_product(charged, identity), projection)
            == ((0, 0, 0),) * 3,
            "arm-blind invariant-to-charged block became nonzero")

    # Conditional 2-by-3 ANOVA diagnostic.  This is an exhaustive algebraic
    # test, not an assertion that the physical packet supplies these rows.
    anova_checks = 0
    for source_total in range(P):
        for m0 in range(P):
            for m1 in range(P):
                for m2 in range(P):
                    common_by_arm = (m0, m1, m2)
                    wing_by_arm = tuple(
                        (source_total - value) % P
                        for value in common_by_arm
                    )
                    target_by_arm = tuple(
                        (common_gain * common_by_arm[index]
                         + wing_gain * wing_by_arm[index]) % P
                        for index in range(3)
                    )
                    rectangles = tuple(
                        (common_by_arm[i] + wing_by_arm[j]
                         - common_by_arm[j] - wing_by_arm[i]) % P
                        for i, j in ((0, 1), (0, 2), (1, 2))
                    )
                    rectangle_zero = all(value == 0 for value in rectangles)
                    common_constant = len(set(common_by_arm)) == 1
                    target_charged = c3_charged(target_by_arm)
                    require(rectangle_zero == common_constant
                            == (target_charged == (0, 0, 0)),
                            "conditional rectangle criterion failed")
                    require(target_charged == tuple(
                        10 * value % P
                        for value in c3_charged(common_by_arm)
                    ), "conditional Q3 target formula failed")
                    anova_checks += 1
    require(anova_checks == P**4,
            "conditional ANOVA universe changed")

    source_masses = (
        source_natural_masses[0],
        source_common_masses[0],
        left_wing_masses[0],
    )
    target_masses = (
        target_natural_masses[0],
        target_common_masses[0],
        right_wing_masses[0],
    )
    require(len(set(source_common_masses)) == 1
            and len(set(target_common_masses)) == 1
            and len(set(left_wing_masses)) == 1
            and len(set(right_wing_masses)) == 1,
            "weighted masses stopped being clock-constant")
    require(len(set(piece_counts)) == 1,
            "weighted-piece census stopped being clock-constant")

    print("LRC14 ROOT-ZERO CLOCK-BLIND WING REPLAY")
    print("status=REFUTED_GLOBAL_CLAIM;VERIFIED_EXACT_CLOCK_BLIND_CALCULATION")
    print(f"dependency_pins=LF-normalized count={len(DEPENDENCIES)}")
    print("prefix_identity=natural_Q_(3,{1,2}) equals THM2749_marked_prefix for all_7_clocks")
    print(f"piece_counts=(A,M,L;B,M,R)={piece_counts[0]}")
    print(f"weighted_masses=(A,M,L)={source_masses}")
    print(f"weighted_masses=(B,M,R)={target_masses}")
    print("weighted_piece_identities=A=M_disjoint_union_L; B=M_disjoint_union_R for all_7_clocks")
    print(f"natural_vectors=source:{source_natural_vector} target:{target_natural_vector}")
    print(f"common_vectors=source:{source_common_vector} target:{target_common_vector}")
    print(f"wing_vectors=left:{left_wing_vector} right:{right_wing_vector}")
    print("functional_identities=S(A)=S(M)+S(L); T(B)=T(M)+T(R)")
    print("root_profiles=common_source12:9/det1 common_target1:4/det1")
    print("root_profiles=left12:9/det1 right1:5/det12 natural_source12:5/det12 natural_target1:9/det1")
    print(f"normalized_ratios=common_clutch_gain:{common_gain} formal_wing_ratio:{wing_gain} folded_natural_ratio:{natural_gain}")
    print("internal_grading=source(M,L):(9,9) target(M,R):(4,5) formal_diagonal_comparison:(12,2); ratio7_after_fold_only")
    print(f"wing_shear={SHEAR_AMPLITUDE} v13=1 (shear/26)_mod13=12")
    print("shear_identity=T(B)-S(A)=T(R)-S(L)")
    print(f"hypothetical_external_arm_gains={hypothetical_arm_gains} mean={mean} Qg={q_profile} normalized_C3={hypothetical_fourier}")
    print("actual_external_arm_typing=I_arm_tensor_internal_T; Q(I_tensor_T)Pi=0")
    print("conditional_external_formula=if_Q3(S_i)=0_on_one_common_carrier_then_Q3(T_i)=10*Q3(M_i); carrier_not_constructed")
    print(f"conditional_ANOVA_exhaustive={anova_checks}_2x3_tables rectangles_zero_iff_Q3T_zero")
    print("FIRST_FAILURE: legacy natural source omits the source-one clock factor of the marked common section")
    print("SCOPE: exact hostile replay on mismatched carriers; no physical wing decomposition, external C3 selector, endpoint current, row exclusion, or LRC14")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
