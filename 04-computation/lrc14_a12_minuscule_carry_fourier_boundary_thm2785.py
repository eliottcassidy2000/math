#!/usr/bin/env python3
"""Exact companion for THM-2785.

The script rebuilds two THM-2672 carry-rebase facet banks:

* the canonical rail-zero bank used in THM-2672; and
* the rail-eight, source-root-12 configuration matching THM-2749's local
  source/target typing.

It verifies the type-A12 minuscule-weight orbit, the collapse of that orbit
to one class in P(A12)/Q(A12), the genuine ordinary spatial Fourier
character at the physical deep frequency 2*13^5, and the exact failure of
the rail-eight selected component to meet THM-2749's semantic section.

The Fourier calculation is ordinary interval integration.  It is not a
THM-2334 endpoint coefficient, a virtual nerve homology class, or a physical
carrier-allocation current.
"""

from fractions import Fraction
import hashlib
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


DEPENDENCIES = {
    "lrc14_slope7_rebase_facet_torsor_thm2672.py":
        "722c86b1e74df79cbc5d78da47be7c03317437c25d2783748e5789a9143c1347",
    "lrc14_fully_marked_root_zero_clutch_thm2749.py":
        "93b46b2701db8f72d00fa2ae131f9d9afd3200f32998959af3bb2e1fa2f56841",
}

for dependency, expected_hash in DEPENDENCIES.items():
    payload = (COMP / dependency).read_bytes().replace(b"\r\n", b"\n")
    actual_hash = hashlib.sha256(payload).hexdigest()
    require(actual_hash == expected_hash,
            f"audited dependency changed: {dependency}")


import lrc14_predecessor_carry_private_root_atlas_thm2640 as atlas
import lrc14_slope7_rebase_facet_torsor_thm2672 as rebase
import lrc14_fully_marked_root_zero_clutch_thm2749 as marked


P = 13
R = P**6
S = P**5
DEEP = 2 * S
MISSING_CARRY = 6


def build_rebase_facets(module, rails, present, prefixes, *,
                        rail, sector, edge, kappa, h,
                        expected_metadata, expected_length):
    """Rebuild all thirteen exact selected components for one configuration."""
    result = atlas.shard((rail, rail + 1))
    content, metadata, rows = result[1], result[5], result[6]
    require(content == 26, "primitive content changed")
    require(metadata == (expected_metadata,), "selected rail metadata changed")

    roots = tuple(
        (2 * carry + (2 * h + kappa) // P + (edge == 0)) % P
        for carry in range(P)
    )
    flags = tuple(
        atlas.is_unit(
            rows[0][sector][edge][carry][kappa][h],
            roots[carry],
            content,
        )
        for carry in range(P)
    )
    unit_carries = tuple(c for c, flag in enumerate(flags) if flag)
    require(
        unit_carries == (0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12),
        "selected twelve-carry row changed",
    )

    pieces = rails[rail][3]
    rail_support = [
        (left, right) for left, right, weight in pieces if weight
    ]
    delayed_cache = [
        [[{} for _ in range(P)] for _ in range(atlas.Q7)]
        for _ in range(2)
    ]
    shifted_rail = {}
    shifted_present = {}
    facets = []

    for source_carry in range(P):
        missing_label = (
            2 * (MISSING_CARRY - source_carry)
        ) % P
        labels = tuple(
            label for label in range(P) if label != missing_label
        )
        require(
            tuple(sorted(
                (source_carry + 7 * label) % P for label in labels
            )) == unit_carries,
            "active labels stopped landing on the unit-carry set",
        )
        require(
            all(
                (roots[source_carry] + label) % P != 0
                for label in labels
            )
            and (roots[source_carry] + missing_label) % P == 0,
            "private-root/missing-label law changed",
        )

        anchor = labels[0]
        anchor_shift = 7 * anchor * atlas.T // R
        common_rail = rebase.shift_weighted(
            pieces, anchor_shift, atlas.T
        )
        for label in labels[1:]:
            key = (rail, label)
            if key not in shifted_rail:
                shift = 7 * label * atlas.T // R
                shifted_rail[key] = rebase.shift_union(
                    rail_support, shift, atlas.T
                )
            common_rail = atlas.old.intersect_weighted_union(
                common_rail, shifted_rail[key]
            )
        require(common_rail, "translated rail supports became disjoint")

        total = 0
        first_component = None
        for clock in range(atlas.Q7):
            source_present = present[clock, (-h) % P]
            common = atlas.old.intersect_weighted_union(
                common_rail,
                rebase.shift_union(
                    source_present, anchor_shift, atlas.T
                ),
            )
            for label in labels[1:]:
                if not common:
                    break
                key = (clock, h, label)
                if key not in shifted_present:
                    shift = 7 * label * atlas.T // R
                    shifted_present[key] = rebase.shift_union(
                        source_present, shift, atlas.T
                    )
                common = atlas.old.intersect_weighted_union(
                    common, shifted_present[key]
                )
            if not common:
                continue

            half = rebase.intersect_root_half(
                common, module.C3, edge, roots[source_carry]
            )
            component = rebase.first_delayed_component(
                half,
                prefixes[sector][clock][h][kappa],
                source_carry,
            )
            if component is not None:
                candidate = (clock, component)
                if (
                    first_component is None
                    or component < first_component[1]
                ):
                    first_component = candidate
            values = atlas.delayed_carry_pair(
                half,
                prefixes[sector][clock][h],
                delayed_cache[sector][clock][h],
            )
            total += values[source_carry][kappa]

        require(
            total > 0 and first_component is not None,
            "one carry-rebased facet lost its positive component",
        )
        clock, (left_raw, right_raw) = first_component
        left = Fraction(left_raw, R * atlas.T)
        right = Fraction(right_raw, R * atlas.T)
        facets.append(
            (source_carry, missing_label, total, clock, left, right)
        )

    require(
        tuple(sorted(row[1] for row in facets)) == tuple(range(P)),
        "source-carry/missing-label map stopped being bijective",
    )
    require(
        {row[5] - row[4] for row in facets} == {expected_length},
        "selected component length changed",
    )
    require(
        all(row[3] == 1 for row in facets),
        "selected first components left delayed clock one",
    )

    # The section is affine after choosing carry zero as origin.  The wrap
    # is deliberately placed between c=0 and c=1.
    left0, right0 = facets[0][4], facets[0][5]
    for carry in range(1, P):
        shift = Fraction(carry - P, R)
        require(
            facets[carry][4] == left0 + shift
            and facets[carry][5] == right0 + shift,
            "canonical component section lost its affine translation law",
        )
    return tuple(facets), roots


def scaled_minuscule_weight(missing_label):
    """Return 13*((1/13)1-e_m), an integral model of -lambda_m."""
    vector = [1] * P
    vector[missing_label] -= P
    return tuple(vector)


def root_vector(positive, negative):
    vector = [0] * P
    vector[positive] = 1
    vector[negative] = -1
    return tuple(vector)


def interval_overlap(interval, intervals, denominator):
    """Physical length of an interval's overlap with scaled intervals."""
    low = interval[0] * denominator
    high = interval[1] * denominator
    total = Fraction(0)
    for left, right in intervals:
        total += max(
            Fraction(0),
            min(high, Fraction(right)) - max(low, Fraction(left)),
        )
    return total / denominator


def weighted_interval_overlap(interval, pieces, denominator):
    """Physical length covered by positive pieces of a weighted partition."""
    low = interval[0] * denominator
    high = interval[1] * denominator
    total = Fraction(0)
    for left, right, weight in pieces:
        if not weight:
            continue
        total += max(
            Fraction(0),
            min(high, Fraction(right)) - max(low, Fraction(left)),
        )
    return total / denominator


def constant_cover_weight(interval, pieces, denominator):
    """Return the common positive weight covering an entire interval."""
    low = interval[0] * denominator
    high = interval[1] * denominator
    hits = []
    covered = Fraction(0)
    for left, right, weight in pieces:
        overlap = max(
            Fraction(0),
            min(high, Fraction(right)) - max(low, Fraction(left)),
        )
        if overlap:
            hits.append(weight)
            covered += overlap
    require(covered == high - low, "weighted profile does not cover interval")
    require(hits and len(set(hits)) == 1 and hits[0] > 0,
            "interval does not have one positive profile weight")
    return hits[0]


def verify_a12(facets):
    """Verify the minuscule orbit and its constant P/Q class."""
    missing = tuple(row[1] for row in facets)
    weights = tuple(scaled_minuscule_weight(label) for label in missing)
    require(len(set(weights)) == P, "minuscule orbit lost a weight")
    require(all(sum(weight) == 0 for weight in weights),
            "projected facet weight left the A12 standard plane")

    # For x in P(A12), all entries of 13*x are congruent modulo 13.
    # The convention gamma=[lambda_0] assigns class
    # class(x)=-13*x_i mod13.  Every -lambda_m therefore has class -gamma.
    quotient_classes = []
    for weight in weights:
        residues = {coordinate % P for coordinate in weight}
        require(len(residues) == 1, "weight left P(A12)")
        quotient_classes.append((-weight[0]) % P)
    require(tuple(quotient_classes) == (P - 1,) * P,
            "minuscule orbit stopped collapsing to -gamma")

    for left in range(P):
        for right in range(P):
            dot_scaled = sum(
                weights[left][slot] * weights[right][slot]
                for slot in range(P)
            )
            expected_dot = 12 * P if left == right else -P
            require(
                dot_scaled == expected_dot,
                "minuscule orbit stopped being a regular A12 simplex",
            )
            if left == right:
                continue
            difference = tuple(
                (weights[left][slot] - weights[right][slot]) // P
                for slot in range(P)
            )
            require(
                all(
                    weights[left][slot] - weights[right][slot]
                    == P * difference[slot]
                    for slot in range(P)
                )
                and sorted(difference) == [-1] + [0] * 11 + [1],
                "two minuscule weights stopped differing by an A12 root",
            )

    # One slope-seven step c->c+7 sends m->m-1.  Its weight increment is
    # e_m-e_(m-1), a simple cyclic A12 root.
    for carry in range(P):
        successor = (carry + 7) % P
        label = missing[carry]
        next_label = missing[successor]
        require(next_label == (label - 1) % P,
                "slope-seven missing-label step changed")
        difference = tuple(
            (weights[successor][slot] - weights[carry][slot]) // P
            for slot in range(P)
        )
        require(
            difference == root_vector(label, next_label),
            "slope-seven minuscule increment stopped being the expected root",
        )

    # For the actual weights w_c=scaled_weight/13, the cyclic Gram sums are
    # A_0=12 and A_d=-1 for d!=0.  Hence every nontrivial unnormalized DFT
    # has squared norm 12-sum_(d!=0)zeta^(kd)=13.
    gram_sums = []
    for displacement in range(P):
        gram_scaled = sum(
            sum(
                weights[carry][slot]
                * weights[(carry + displacement) % P][slot]
                for slot in range(P)
            )
            for carry in range(P)
        )
        require(
            gram_scaled % (P * P) == 0,
            "cyclic Gram sum left the weight-lattice denominator",
        )
        gram_sums.append(gram_scaled // (P * P))
    require(tuple(gram_sums) == (12,) + (-1,) * 12,
            "regular-simplex cyclic Gram profile changed")

    # Exact cyclotomic reduction separately confirms that no nontrivial
    # coefficient-space DFT vector vanishes.
    for frequency in range(1, P):
        reduced_coordinates = []
        for slot in range(P):
            coefficients = [0] * (P - 1)
            for carry, weight in enumerate(weights):
                exponent = (-frequency * carry) % P
                if exponent == P - 1:
                    for basis in range(P - 1):
                        coefficients[basis] -= weight[slot]
                else:
                    coefficients[exponent] += weight[slot]
            reduced_coordinates.append(tuple(coefficients))
        require(any(any(row) for row in reduced_coordinates),
                "one nontrivial minuscule-orbit DFT vanished")
    return missing, tuple(quotient_classes), tuple(gram_sums)


def verify_fourier_section(facets, component_length, frequency_multiple):
    """Verify the exact C13 character after annihilating the odometer kernel."""
    require(frequency_multiple != 0, "zero frequency is not a charged mode")
    frequency = frequency_multiple * S
    require(R == P * S and frequency % S == 0,
            "Fourier descent scale changed")

    # The exact component integral is nonzero because
    # 0 < |N|*length < 1, so sin(pi*N*length) cannot vanish.
    phase_length = abs(frequency) * component_length
    require(0 < phase_length < 1,
            "selected component nonvanishing test left its exact range")

    exponents = tuple(
        (-frequency_multiple * carry) % P for carry in range(P)
    )
    require(
        all(
            exponents[(carry + 7) % P] - exponents[carry]
            == (-7 * frequency_multiple) % P
            or (
                exponents[(carry + 7) % P] - exponents[carry]
            ) % P == (-7 * frequency_multiple) % P
            for carry in range(P)
        ),
        "slope step lost its constant Fourier multiplier",
    )

    # The chosen component section is not cyclic.  Its slope lift differs
    # from the selected successor by either 0 or the THM-2657 kernel 1/S.
    # All S-divisible frequencies annihilate both possibilities.
    kernel_steps = []
    for carry in range(P):
        successor = (carry + 7) % P
        source_left = facets[carry][4]
        source_right = facets[carry][5]
        kernel_left = (
            source_left + Fraction(7, R) - facets[successor][4]
        )
        kernel_right = (
            source_right + Fraction(7, R) - facets[successor][5]
        )
        require(kernel_left == kernel_right,
                "component lift changed interval length")
        require(kernel_left in (Fraction(0), Fraction(1, S)),
                "component lift left the odometer kernel")
        if kernel_left:
            require((frequency * kernel_left).denominator == 1,
                    "Fourier mode did not annihilate the odometer correction")
        kernel_steps.append(kernel_left)
    require(
        sum(kernel_steps, Fraction(0)) == Fraction(7, S),
        "thirteen slope lifts lost the THM-2657 cocycle-seven total",
    )

    # The iff criterion is algebraic: the internal ratio is exp(-2pi i N/R);
    # it is a thirteenth root exactly when S divides the integer N.
    require(R // P == S, "character-descent divisibility changed")
    test_frequencies = (
        -2 * S - 3, -2 * S, -S - 1, -S, -1,
        0, 1, S, S + 1, 2 * S, 2 * S + 3,
    )
    for test_frequency in test_frequencies:
        algebraic_descent = (
            (P * test_frequency) % R == 0
        )
        require(
            algebraic_descent == (test_frequency % S == 0),
            "finite character-descent arithmetic control failed",
        )
    return frequency, phase_length, exponents, tuple(kernel_steps)


def main():
    require(
        (atlas.P, atlas.R, atlas.S, atlas.T, marked.SHIFT)
        == (
            13,
            4_826_809,
            371_293,
            297_836_897_838_480,
            431_933_040,
        ),
        "canonical scale data changed",
    )

    module, _, _, _, rails, present, _ = (
        atlas.core.build_carrier_data()
    )
    require(module.C3 == DEEP, "physical deep speed changed")
    prefixes = atlas.build_pair_prefixes(module)

    canonical_length = Fraction(3, 7 * P**11)
    canonical, canonical_roots = build_rebase_facets(
        module,
        rails,
        present,
        prefixes,
        rail=0,
        sector=0,
        edge=0,
        kappa=0,
        h=6,
        expected_metadata=(1, 0, 12),
        expected_length=canonical_length,
    )
    missing, quotient_classes, gram_sums = verify_a12(canonical)
    canonical_modes = tuple(
        verify_fourier_section(canonical, canonical_length, multiple)
        for multiple in range(1, P)
    )
    canonical_fourier = canonical_modes[1]
    require(
        canonical_roots
        == (1, 3, 5, 7, 9, 11, 0, 2, 4, 6, 8, 10, 12),
        "canonical private-root row changed",
    )
    require(
        canonical_fourier[2]
        == tuple((-2 * carry) % P for carry in range(P)),
        "physical deep-frequency character changed",
    )
    require(
        (-7 * 2) % P == P - 1,
        "slope step stopped multiplying the deep character by zeta^-1",
    )

    rail8_length = Fraction(1, 28 * P**11)
    rail8, rail8_roots = build_rebase_facets(
        module,
        rails,
        present,
        prefixes,
        rail=8,
        sector=0,
        edge=1,
        kappa=1,
        h=6,
        expected_metadata=(1, 4, 12),
        expected_length=rail8_length,
    )
    rail8_missing, rail8_classes, rail8_gram_sums = verify_a12(rail8)
    rail8_modes = tuple(
        verify_fourier_section(rail8, rail8_length, multiple)
        for multiple in range(1, P)
    )
    rail8_fourier = rail8_modes[1]
    require(
        rail8_roots == canonical_roots
        and rail8_missing == missing
        and rail8_classes == quotient_classes,
        "rail-eight facet orbit left the canonical A12 typing",
    )
    require(
        rail8_gram_sums == gram_sums
        and tuple(
            mode[2][1] for mode in canonical_modes
        ) == tuple((-multiple) % P for multiple in range(1, P))
        and tuple(
            mode[2][1] for mode in rail8_modes
        ) == tuple((-multiple) % P for multiple in range(1, P)),
        "all-character physical carry bank changed",
    )

    # THM-2749 typing uses source carry 12/root 12 and target carry 6/root 0
    # in this configuration.  The selected carry-12 component is a genuine
    # positive rail component, but it misses E3 and hence the full section
    # S_(1,0,4)=E3 intersect F_(1,0,4) exactly.
    source_row = rail8[12]
    target_row = rail8[6]
    source_interval = (source_row[4], source_row[5])
    target_interval = (target_row[4], target_row[5])
    translated_source_interval = (
        source_interval[0] + Fraction(7, R),
        source_interval[1] + Fraction(7, R),
    )
    require(
        rail8_roots[12] == 12
        and rail8_roots[6] == 0
        and source_row[2] > 0
        and target_row[2] > 0,
        "rail-eight source/target carry typing changed",
    )

    source_e3 = marked.semantic.exclusive_source(module, 3)
    terminal_fork = marked.semantic.deepest_fork(module)
    semantic_prefixes = marked.build_semantic_prefixes(
        module, terminal_fork
    )
    sections = marked.semantic_sections(module, source_e3, 0, 4)
    source_vector, target_vector, details = marked.fully_marked_vectors(
        module,
        rails,
        present,
        semantic_prefixes,
        sections,
        8,
    )
    require(
        source_vector == target_vector and source_vector[1] > 0,
        "THM-2749 rail-eight marked positive control changed",
    )

    target_rail_pullback = marked.clutch.shift_weighted(
        rails[8][3], -marked.SHIFT
    )
    source_rail_weight = constant_cover_weight(
        source_interval, rails[8][3], marked.T
    )
    target_rail_weight = constant_cover_weight(
        source_interval, target_rail_pullback, marked.T
    )
    require(
        source_rail_weight == target_rail_weight == 27_581_135_604,
        "rail-eight equal-weight overlap changed",
    )
    require(
        interval_overlap(source_interval, source_e3, marked.T) == 0
        and interval_overlap(source_interval, sections[1], marked.T) == 0
        and weighted_interval_overlap(
            source_interval, details[1][0], marked.T
        ) == 0,
        "selected source-root-12 component acquired a semantic atom",
    )
    require(
        weighted_interval_overlap(
            target_interval, details[1][1], marked.T
        ) == 0,
        "selected carry-six component acquired a marked target atom",
    )
    require(
        weighted_interval_overlap(
            translated_source_interval, details[1][1], marked.T
        ) == 0,
        "actual translate of source component acquired a marked target atom",
    )

    # The selected section itself witnesses the nonsplit extension:
    # I_12+7/R is I_6 shifted by the kernel 1/S.  The deep Fourier mode
    # annihilates that discrepancy, but the semantic section does not.
    require(
        source_interval[0] + Fraction(7, R) - target_interval[0]
        == Fraction(1, S)
        and source_interval[1] + Fraction(7, R) - target_interval[1]
        == Fraction(1, S),
        "rail-eight source/target selected components lost the kernel shift",
    )
    require(
        DEEP * Fraction(1, S) == 2,
        "deep Fourier character stopped annihilating the kernel shift",
    )

    print("THM2785 A12 MINUSCULE CARRY FOURIER BOUNDARY EXACT REFEREE")
    print(f"scales=(P={P},R={R},S={S},C3={DEEP})")
    print("canonical_config=(rail0,sector0,edge0,kappa0,h6) "
          "metadata=(1,0,12)")
    print(f"canonical_missing_labels={missing}")
    print(f"canonical_component_length={canonical_length}")
    print(f"a12_weight_orbit_size={P} pair_differences=A12_roots")
    print(f"a12_regular_simplex_cyclic_gram={gram_sums}")
    print("a12_nontrivial_vector_DFT_norm_squared=13 "
          "for_each_of_12_characters")
    print(f"a12_P_over_Q_classes={quotient_classes} "
          "meaning=all_minus_gamma")
    print("direct_carry_to_P_over_Q=constant; "
          "no facet-induced quotient identification")
    print("fourier_descent_criterion=integer_N descends_to_C13_character "
          "iff 13^5 divides N")
    print("physical_component_spectrum=N=a*13^5 realizes character "
          "-a*gamma for every a=1,...,12")
    print(f"canonical_C3_character_exponents={canonical_fourier[2]}")
    print("canonical_C3_character_class=-2gamma "
          "magnitude=sin(6*pi/(7*13^6))/(2*pi*13^5)")
    print("slope_step=c_to_c_plus_7 character_multiplier=zeta^-1")
    print(f"canonical_odometer_kernel_steps="
          f"{tuple(str(value) for value in canonical_fourier[3])}")
    print("canonical_thirteen_step_kernel_sum=7/13^5 "
          "annihilated_by_C3")
    print("rail8_config=(rail8,sector0,edge1,kappa1,h6) "
          "metadata=(1,4,12)")
    print(f"rail8_component_length={rail8_length}")
    print("rail8_C3_character_class=-2gamma "
          "magnitude=sin(pi/(14*13^6))/(2*pi*13^5)")
    print(f"rail8_source_I12={source_interval} root=12 "
          f"positive_total={source_row[2]}")
    print(f"rail8_target_I6={target_interval} root=0 "
          f"positive_total={target_row[2]}")
    print(f"rail8_equal_weight={source_rail_weight}")
    print("rail8_I12_intersection_E3=0 "
          "intersection_S_(1,0,4)=0 "
          "intersection_full_marked_source=0")
    print("rail8_I6_and_translated_I12_intersection_"
          "full_marked_target=0")
    print("rail8_selected_transport=I12+7/R=I6+1/S; "
          "Fourier_C3_forgets_1/S but semantic_carrier_does_not")
    print("comparison=THM2771 target/Bockstein coordinates and THM2772 "
          "allocation K4 are not supplied")
    print("scope=physical component-indicator Fourier current only; "
          "no virtual homology, THM2334 endpoint current, common atom, "
          "row decrement, or LRC14")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
