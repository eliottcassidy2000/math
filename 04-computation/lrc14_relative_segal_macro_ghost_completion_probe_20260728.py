#!/usr/bin/env python3
"""Exact relative-Segal audit of the transverse-17 packet two-cycle.

The transverse C221 scout supplies two strict packet points over delayed
phases 4/17 and 13/17.  This probe forgets the transverse translation and
asks what remains in the ordinary C_(13^6) affine category.

There is an exact nonconstant two-cycle under state-dependent ``D^2`` macro
arrows.  Every factorization of either macro into two ordinary degree-one
arrows has a forced delayed midpoint outside the inherited word.  Requiring
literal clock gluing forces the two midpoint labels to be the diagonal
self-lines (4,4) and (1,1).  Exact lifts of those diagonal ghosts can retain
rail, present, and root data, but none retains the delayed sector.  Moreover
the two macro endpoints have no THM-2305 exclusive-source leg.  Thus the
construction is a minimal relative Segal completion of packet support, not a
semantic endpoint/current cospan.

No scalar row, row exclusion, or LRC(14) conclusion is asserted.
"""

from fractions import Fraction
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc14_mod17_transverse_phase_typed_cycle_probe_20260728 as transverse
import lrc14_odometer_alternating_lift_labelled_tail_scout_20260728 as old
import lrc14_phase_cycle_semantic_gate_probe_20260728 as semantic


P = 13
R = P**6
L = 17 * R
POINT_NUMERATORS = (39_123_022, 41_305_508)
EXPECTED_MACRO_TRANSLATIONS = (4_472_399, 1_954_775)
ENDPOINT_CLOCKS = ((1, 4), (4, 1))
GHOST_CLOCKS = ((4, 4), (1, 1))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def frac(value):
    return value - value.numerator // value.denominator


def clock_pair_from_numerator(numerator):
    return (
        transverse.clock_from_stalk_numerator(numerator, P),
        transverse.clock_from_stalk_numerator(numerator, P**2),
    )


def macro_translation(source_numerator, target_numerator):
    difference = target_numerator - P**2 * source_numerator
    require(difference % 17 == 0,
            "D^2 endpoint interpolation stopped descending to C_R")
    return (difference // 17) % R


def prepare_carrier():
    (module, _prefixes, _whole_prefixes, _digit_masses, rails,
     present, present_starts) = old.cross.build_carrier_data()
    pair_prefixes = old.private.build_pair_prefixes(module)
    unit = old.T // R
    require(old.T % R == 0, "carrier grid stopped resolving R fibres")
    prepared_present = {
        (shallow, h): transverse.prepare_unweighted(
            present[shallow, (-h) % P]
        )
        for shallow in range(7) for h in range(P)
    }
    prepared_rails = tuple(
        transverse.prepare_weighted(rail[3]) for rail in rails
    )
    rails_by_owner = {
        owner: tuple(index for index, rail in enumerate(rails)
                     if rail[1] == owner)
        for owner in range(7)
    }
    return (module, rails, present, present_starts, pair_prefixes, unit,
            prepared_present, prepared_rails, rails_by_owner)


def audit_macros():
    points = tuple(Fraction(numerator, L)
                   for numerator in POINT_NUMERATORS)
    translations = tuple(
        macro_translation(POINT_NUMERATORS[index],
                          POINT_NUMERATORS[1 - index])
        for index in range(2)
    )
    require(translations == EXPECTED_MACRO_TRANSLATIONS,
            "ordinary D^2 macro translations changed")

    for index, (point, translation) in enumerate(zip(points, translations)):
        target = frac(P**2 * point + Fraction(translation, R))
        require(target == points[1 - index],
                "ordinary D^2 macro stopped exchanging the packet endpoints")

    delayed = tuple(frac(R * point) for point in points)
    require(delayed == (Fraction(4, 17), Fraction(13, 17)),
            "macro endpoint phases changed")
    orbit = [delayed[0]]
    for _ in range(4):
        orbit.append(semantic.B0(orbit[-1]))
    require(tuple(orbit) == (
        Fraction(4, 17), Fraction(1, 17), Fraction(13, 17),
        Fraction(16, 17), Fraction(4, 17),
    ), "ordinary delayed four-cycle changed")
    require(semantic.B0(semantic.B0(delayed[0])) == delayed[1]
            and semantic.B0(semantic.B0(delayed[1])) == delayed[0],
            "B0^2 stopped exchanging the macro endpoints")
    require(tuple(clock_pair_from_numerator(numerator)
                  for numerator in POINT_NUMERATORS) == ENDPOINT_CLOCKS,
            "endpoint reciprocal clock labels changed")

    endpoint_rows = []
    for point, phase in zip(points, delayed):
        bits = semantic.literal_bits(phase)
        sources = tuple(owner + 1 for owner in range(3)
                        if semantic.exclusive_source(point, owner))
        terminal = semantic.iterate(semantic.B0, point, 6)
        word = semantic.terminal_word(terminal, 2)
        require(bits[0] and all(bits[1]) and not bits[2]
                and bits[3] and not bits[4],
                "macro endpoint semantic debt changed")
        require(not sources and word == (1, 2),
                "macro endpoint acquired a THM-2305 source leg")
        endpoint_rows.append((phase, bits, sources, terminal, word))

    return points, translations, tuple(delayed), tuple(orbit), tuple(endpoint_rows)


def audit_all_factorizations(translations, carrier):
    (module, rails, _present, _present_starts, pair_prefixes, unit,
     prepared_present, prepared_rails, rails_by_owner) = carrier

    rows = []
    selected_rows = []
    for index, (source_numerator, macro_translation_value,
                expected_clock) in enumerate(zip(
                    POINT_NUMERATORS, translations, GHOST_CLOCKS)):
        target_numerator = POINT_NUMERATORS[1 - index]
        expected_phase = (Fraction(1, 17), Fraction(16, 17))[index]
        counts = {
            "factorizations": 0,
            "clock_glued": 0,
            "present": 0,
            "rail": 0,
            "present_and_rail": 0,
            "present_rail_and_both_roots": 0,
        }
        first_full_local = None

        # Every factorization has the form
        #
        #   (A_b D) o (A_a D) = A_(b+13a) D^2,
        #
        # and for each a in C_R there is a unique b.  This loop exhausts the
        # full finite universe rather than sampling a preferred digit section.
        for inner_translation in range(R):
            outer_translation = (
                macro_translation_value - P * inner_translation
            ) % R
            midpoint_numerator = (
                P * source_numerator + 17 * inner_translation
            ) % L
            composite_numerator = (
                P * midpoint_numerator + 17 * outer_translation
            ) % L
            require(composite_numerator == target_numerator,
                    "ordinary two-edge factorization identity failed")
            require(midpoint_numerator % 17
                    == (1, 16)[index],
                    "an integral factorization escaped the forced ghost phase")
            counts["factorizations"] += 1

            if clock_pair_from_numerator(midpoint_numerator) != expected_clock:
                continue
            counts["clock_glued"] += 1
            shallow, owner = expected_clock
            h = P * (midpoint_numerator % 17) // 17
            scaled_point = midpoint_numerator * unit
            present_ok = transverse.inside_prepared_unweighted(
                scaled_point, prepared_present[shallow, h]
            )
            active_rails = tuple(
                rail_index for rail_index in rails_by_owner[owner]
                if transverse.inside_prepared_weighted(
                    scaled_point, prepared_rails[rail_index]
                )
            )
            rail_ok = bool(active_rails)
            counts["present"] += present_ok
            counts["rail"] += rail_ok
            counts["present_and_rail"] += present_ok and rail_ok

            if not (present_ok and rail_ok):
                continue
            midpoint = Fraction(midpoint_numerator, L)
            roots = old.root_memberships(module, midpoint * old.T)
            d_roots = old.root_memberships(
                module, frac(P * midpoint) * old.T
            )
            both_roots = bool(roots) and bool(d_roots)
            counts["present_rail_and_both_roots"] += both_roots
            if both_roots and first_full_local is None:
                first_full_local = (
                    inner_translation, outer_translation,
                    midpoint_numerator, active_rails, roots, d_roots,
                )

        require(counts == (
            {
                "factorizations": R,
                "clock_glued": 106_080,
                "present": 27_249,
                "rail": 8_160,
                "present_and_rail": 2_925,
                "present_rail_and_both_roots": 2_925,
            } if index == 0 else {
                "factorizations": R,
                "clock_glued": 106_080,
                "present": 15_988,
                "rail": 8_160,
                "present_and_rail": 2_479,
                "present_rail_and_both_roots": 2_479,
            }
        ), "ghost lift census changed")
        require(first_full_local is not None,
                "no rail/present/root ghost lift survived")

        (inner_translation, outer_translation, midpoint_numerator,
         active_rails, roots, d_roots) = first_full_local
        midpoint = Fraction(midpoint_numerator, L)
        predecessor, carry, h, kappa, phase = transverse.state_digits(midpoint)
        shallow, owner = clock_pair_from_numerator(midpoint_numerator)
        delayed_sectors = tuple(
            sector for sector in range(2)
            if transverse.strict_in_unweighted(
                phase * old.T,
                old.prefix_intervals(
                    pair_prefixes[sector][shallow][h][kappa]
                ),
            )
        )
        literal_bits = semantic.literal_bits(phase)
        exclusive_sources = tuple(
            source + 1 for source in range(3)
            if semantic.exclusive_source(midpoint, source)
        )

        expected_selected = (
            (1_485_215, 4_471_840, 41_513_423,
             2_441_966, 7, 0, 1, (4, 4), 11,
             ((0, 2), (1, 1)), ((0, 2), (1, 1)), (), ()),
            (4_459_459, 1_903_516, 38_392_136,
             2_258_360, 0, 12, 0, (1, 1), 11,
             ((0, 2), (1, 1)), ((0, 12), (1, 11)), (), (2,)),
        )[index]
        selected = (
            inner_translation, outer_translation, midpoint_numerator,
            predecessor, carry, h, kappa, (shallow, owner),
            len(active_rails), roots, d_roots, delayed_sectors,
            exclusive_sources,
        )
        require(selected == expected_selected,
                "selected maximally local ghost witness changed")
        require(phase == expected_phase,
                "selected ghost phase changed")
        require(literal_bits == (
            False, (True, True, True, True, True), True, False, False
        ), "ghost semantic bit pattern changed")
        require(not delayed_sectors,
                "a ghost midpoint entered the inherited delayed sector")

        rows.append((index, expected_phase, expected_clock, counts))
        selected_rows.append((
            index, inner_translation, outer_translation,
            midpoint_numerator, phase, (shallow, owner),
            (predecessor, carry, h, kappa), len(active_rails),
            roots, d_roots, delayed_sectors, literal_bits,
            exclusive_sources,
        ))

    return tuple(rows), tuple(selected_rows)


def main():
    require((P, R, L) == (13, 4_826_809, 82_055_753),
            "ambient lattice changed")
    points, translations, delayed, orbit, endpoint_rows = audit_macros()
    carrier = prepare_carrier()
    factor_rows, selected_rows = audit_all_factorizations(
        translations, carrier
    )

    # The reciprocal endpoint labels force diagonal midpoint labels:
    # (a,b)->(b,a) can be subdivided only through (b,b), and conversely.
    require((ENDPOINT_CLOCKS[0][1], ENDPOINT_CLOCKS[1][0]) == (4, 4)
            and (ENDPOINT_CLOCKS[1][1], ENDPOINT_CLOCKS[0][0]) == (1, 1),
            "reciprocal-clock diagonal factorization changed")

    print("LRC14 relative Segal macro and ghost-midpoint exact probe")
    print("status=VERIFIED-EXACT packet-support completion; no semantic cospan")
    print(f"parameters=P:{P},R:{R},L:{L}")
    print(f"packet_points={points}; numerators={POINT_NUMERATORS}")
    print(f"ordinary_D2_macro_translations={translations}; root_steps={tuple(2*k % P for k in translations)}")
    print(f"ordinary_D2_packet_cycle=x0--A_{translations[0]}D2-->x1--A_{translations[1]}D2-->x0")
    print(f"delayed_B0_orbit={orbit}")
    print("orbit_grading=grade0 endpoints {4/17,13/17}; grade1 ghosts {1/17,16/17}; B0 toggles grade and B0^2 is reflection")
    print(f"clock_grading=off_diagonal_endpoints:{ENDPOINT_CLOCKS}; forced_diagonal_ghosts:{GHOST_CLOCKS}")
    print(f"endpoint_semantic_rows=(phase,bits,exclusive_sources,D6endpoint,deep_word)={endpoint_rows}")
    print(f"all_ordinary_two_D_factorizations={factor_rows}")
    print(f"selected_clock_glued_ghost_lifts={selected_rows}")
    print("minimality=zero ghosts impossible because forced phases are outside the delayed word; one ghost cannot cover distinct phases 1/17 and16/17; the two selected lifts suffice for the clock-glued four-cycle")
    print("fixed_grading_cycle=(1,4)->(4,4)->(4,1)->(1,1)->(1,4)")
    print("physical_boundary=some diagonal ghosts retain rail,present,current-root,and D-root, but every ghost has no inherited delayed sector and lies outside the nonconstant clock class")
    print("semantic_boundary=macro endpoints fail literal Qa only at c1-safe and lie in no exclusive E_j; adjoining ghosts is a relative Segal factorization completion, not a THM2305 endpoint/current cospan")
    print("next_target=new mapping-cone current absorbing c1/guard debt, or re-root at an actual exclusive source; endpoint equality alone cannot supply either")
    print("SCOPE: no terminal current, scalar row, row exclusion, or LRC14 conclusion")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
