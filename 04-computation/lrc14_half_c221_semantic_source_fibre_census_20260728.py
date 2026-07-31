#!/usr/bin/env python3
"""Full exact semantic-source census of the half and C221 packet fibres.

The previously displayed physical two-cycle centres are not exclusive
THM-2305 sources.  This scout audits the *complete declared candidate
universes* of the THM-2698 half-odometer companion and the transverse-C221
literal-root companion.  It adds only the canonical scalar predicates

    E_j at x,                 Q_(j,sigma) at D^(2j)(x),

where D(x)={13x}.  Affine half/C221 arrows are kept separate from this
ordinary semantic chronology.  Every retained reciprocal pair still has the
original rail, present, clock, carry, private-root, and primitive-unit data.
"""

from __future__ import annotations

from bisect import bisect_right
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path
from typing import NamedTuple
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc14_central_half_odometer_full_local_cycle_thm2698 as half
import lrc14_mod17_transverse_phase_typed_cycle_probe_20260728 as c221


P = 13
R = P**6
L = 17 * R
GUARD = 1
UNITS = (14, 27, 40, 53, 66)
BLOCKERS = (P, P**3, 2 * P**5)
CLOCKS = (2, 4, 6)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def frac(value):
    return value - value.numerator // value.denominator


def centered(value):
    residue = frac(value)
    return residue if residue < Fraction(1, 2) else residue - 1


def distance(value):
    return abs(centered(value))


def danger(speed, value, denominator=14):
    phase = centered(speed * value)
    radius = Fraction(1, denominator)
    return -radius <= phase < radius


def common_safe(value):
    return (
        not danger(GUARD, value, 7)
        and all(not danger(speed, value) for speed in UNITS)
    )


def exclusive_source_label(value):
    if not common_safe(value):
        return 0
    hits = tuple(danger(speed, value) for speed in BLOCKERS)
    if sum(hits) != 1:
        return 0
    return hits.index(True) + 1


def terminal_word(value, source_label):
    """THM-2305 terminal word, with blocker labels numbered 1,2,3."""
    require(source_label in (1, 2, 3), "unknown source label")
    source = source_label - 1
    if not common_safe(value) or danger(BLOCKERS[source], value):
        return None
    word = tuple(
        index + 1
        for index, speed in enumerate(BLOCKERS)
        if index != source and danger(speed, value)
    )
    return word or None


def ordinary_dilate(value, steps):
    return frac(P**steps * value)


def semantic_record(value):
    source = exclusive_source_label(value)
    if not source:
        return 0, None
    endpoint = ordinary_dilate(value, CLOCKS[source - 1])
    return source, terminal_word(endpoint, source)


def danger_numerator(speed, numerator, denominator, danger_denominator=14):
    """Exact half-open danger predicate without constructing Fractions."""
    residue = speed * numerator % denominator
    if 2 * residue < denominator:
        return danger_denominator * residue < denominator
    return danger_denominator * (denominator - residue) <= denominator


def common_safe_numerator(numerator, denominator):
    return (
        not danger_numerator(GUARD, numerator, denominator, 7)
        and all(not danger_numerator(speed, numerator, denominator)
                for speed in UNITS)
    )


def exclusive_source_label_numerator(numerator, denominator):
    if not common_safe_numerator(numerator, denominator):
        return 0
    hits = tuple(danger_numerator(speed, numerator, denominator)
                 for speed in BLOCKERS)
    return hits.index(True) + 1 if sum(hits) == 1 else 0


def terminal_word_numerator(numerator, denominator, source):
    source_index = source - 1
    if (not common_safe_numerator(numerator, denominator)
            or danger_numerator(
                BLOCKERS[source_index], numerator, denominator)):
        return None
    word = tuple(
        index + 1 for index, speed in enumerate(BLOCKERS)
        if index != source_index
        and danger_numerator(speed, numerator, denominator)
    )
    return word or None


def semantic_record_numerator(numerator, denominator):
    source = exclusive_source_label_numerator(numerator, denominator)
    if not source:
        return 0, None
    endpoint_numerator = (
        P**CLOCKS[source - 1] * numerator % denominator
    )
    return (
        source,
        terminal_word_numerator(endpoint_numerator, denominator, source),
    )


def semantic_cospan_radius(value, source, word):
    """Symmetric physical-x radius preserving strict source and endpoint."""
    require(source in (1, 2, 3), "cospan radius needs an exclusive source")
    require(word, "cospan radius needs a nonempty terminal word")
    source_index = source - 1
    rows = []

    # Strict source E_j.
    rows.append((distance(value) - Fraction(1, 7)) / GUARD)
    rows.extend(
        (distance(speed * value) - Fraction(1, 14)) / speed
        for speed in UNITS
    )
    for index, speed in enumerate(BLOCKERS):
        slack = (
            Fraction(1, 14) - distance(speed * value)
            if index == source_index
            else distance(speed * value) - Fraction(1, 14)
        )
        rows.append(slack / speed)

    # Strict Q_(j,word) after ordinary D^k, pulled back to x.
    steps = CLOCKS[source_index]
    endpoint = ordinary_dilate(value, steps)
    expansion = P**steps
    rows.append(
        (distance(endpoint) - Fraction(1, 7)) / (GUARD * expansion)
    )
    rows.extend(
        (distance(speed * endpoint) - Fraction(1, 14))
        / (speed * expansion)
        for speed in UNITS
    )
    for index, speed in enumerate(BLOCKERS):
        hit = index + 1 in word
        slack = (
            Fraction(1, 14) - distance(speed * endpoint)
            if hit
            else distance(speed * endpoint) - Fraction(1, 14)
        )
        rows.append(slack / (speed * expansion))

    require(all(row > 0 for row in rows),
            "semantic source/fork witness is not strict")
    return min(rows)


class HalfSignature(NamedTuple):
    carry: int
    edge: int
    root: int
    d_roots_left: tuple
    d_roots_right: tuple
    source: int
    word: tuple | None


def generate_half_candidates(module, rails, present, present_starts,
                             pair_prefixes):
    source_rows = half.private.shard((0, 14))
    require(source_rows[1] == 26
            and all(meta[0] == 1 for meta in source_rows[5]),
            "source-one THM-2640 bank changed")
    rows = source_rows[6]
    rail_starts = [[left for left, _, _ in rail[3]]
                   for rail in rails[:14]]
    candidates = []
    semantic_counts = Counter()
    filter_counts = Counter()
    reroot_counts = Counter()
    raw_e3_witness = None

    for phase in half.FIXED_PHASES:
        require(terminal_word(phase, 3) == (1, 2),
                "half delayed phase left the deepest-owner fork")
        h = half.floor_fraction(P * phase)
        kappa = half.floor_fraction(2 * P * phase) - 2 * h
        first_n = half.floor_fraction(Fraction(6 * R, P) - phase) + 1
        last_n = half.floor_fraction(Fraction(7 * R, P) - phase)
        for n in range(first_n, last_n + 1):
            value = Fraction(n, R) + Fraction(phase, R)
            raw_denominator = R * phase.denominator
            raw_numerator = n * phase.denominator + phase.numerator
            source, word = semantic_record_numerator(
                raw_numerator, raw_denominator
            )
            filter_counts["raw", phase, source] += 1
            if source == 3 and raw_e3_witness is None:
                raw_e3_witness = (value, phase, word)
            shallow = half.shallow(value)
            if shallow == 0:
                continue
            filter_counts["shallow_nonzero", phase, source] += 1
            coordinate = value * half.T
            rail_indices_without_present = tuple(
                index for index in range(14)
                if half.strict_weighted_member(
                    coordinate, rails[index][3], rail_starts[index]
                )
            )
            if rail_indices_without_present:
                filter_counts["rail_without_present", phase, source] += 1
            if source == 3:
                present_labels = tuple(
                    label for label in range(P)
                    if half.strict_interval_member(
                        coordinate, present[shallow, label],
                        present_starts[shallow, label]
                    )
                )
                if present_labels:
                    reroot_counts["E3_any_present_label", phase] += 1
                for label in present_labels:
                    reroot_counts["E3_present_label", phase, label] += 1
                carry_probe = n % P
                root_probe = any(
                    (
                        root := (
                            2 * carry_probe + (2 * h + kappa) // P
                            + int(edge == 0)
                        ) % P
                    ) != 0
                    and half.half_membership(
                        module, value, edge, root
                    )
                    for _rail_index in rail_indices_without_present
                    for edge in (0, 1)
                )
                if root_probe:
                    reroot_counts[
                        "E3_rail_and_private_root_without_present", phase
                    ] += 1
            if not half.strict_interval_member(
                    coordinate, present[shallow, (-h) % P],
                    present_starts[shallow, (-h) % P]):
                continue
            filter_counts["present", phase, source] += 1
            rail_indices = rail_indices_without_present
            if not rail_indices:
                continue
            filter_counts["rail", phase, source] += 1
            carry = n % P
            root_survived = False
            unit_survived = False
            for rail_index in rail_indices:
                require(rails[rail_index][0] == 1
                        and rails[rail_index][1] == half.owner(value),
                        "half rail typing changed")
                for edge in (0, 1):
                    root = (
                        2 * carry + (2 * h + kappa) // P
                        + int(edge == 0)
                    ) % P
                    if (root == 0
                            or not half.half_membership(
                                module, value, edge, root)):
                        continue
                    root_survived = True
                    values = rows[rail_index][0][edge][carry][kappa][h]
                    if not half.private.is_unit(values, root, 26):
                        continue
                    unit_survived = True
                    d_value = frac(P * value)
                    d_roots = tuple(
                        half.half_roots(module, d_value, target_edge)
                        for target_edge in (0, 1)
                    )
                    candidate = {
                        "x": value,
                        "N": n,
                        "phase": phase,
                        "rail_index": rail_index,
                        "shallow": shallow,
                        "owner": half.owner(value),
                        "carry": carry,
                        "h": h,
                        "kappa": kappa,
                        "edge": edge,
                        "root": root,
                        "d_roots": d_roots,
                        "source": source,
                        "word": word,
                    }
                    candidates.append(candidate)
                    semantic_counts[phase, source, word] += 1
            if root_survived:
                filter_counts["private_root", phase, source] += 1
            if unit_survived:
                filter_counts["primitive_unit", phase, source] += 1

    require(len(candidates) == 332_668,
            "half full-packet candidate universe changed")
    require(tuple(sum(count for (phase, _source, _word), count
                      in semantic_counts.items() if phase == fixed)
                  for fixed in half.FIXED_PHASES) == (156_798, 175_870),
            "half candidate phase census changed")
    require(all(word == (1, 2)
                for (_phase, source, word) in semantic_counts
                if source == 3),
            "half E3 row lost its prescribed deepest fork")
    require(tuple(filter_counts["raw", phase, 3]
                  for phase in half.FIXED_PHASES) == (12_848, 12_848)
            and tuple(filter_counts["shallow_nonzero", phase, 3]
                      for phase in half.FIXED_PHASES)
            == (12_848, 12_848)
            and tuple(filter_counts["present", phase, 3]
                      for phase in half.FIXED_PHASES) == (0, 0),
            "half E3 present-factor stopping row changed")
    require(not any(key[0] == "E3_any_present_label"
                    for key in reroot_counts),
            "half E3 unexpectedly entered an alternative present label")
    require(tuple(reroot_counts[
                        "E3_rail_and_private_root_without_present", phase
                    ] for phase in half.FIXED_PHASES)
            == (12_848, 12_848),
            "half E3 rail/root survivor changed")
    require(raw_e3_witness is not None
            and raw_e3_witness[2] == (1, 2),
            "half raw fibre lost the E3-to-fork intersection")
    return (
        candidates, semantic_counts, rail_starts,
        filter_counts, reroot_counts, raw_e3_witness,
    )


def half_signature(candidate):
    return HalfSignature(
        candidate["carry"], candidate["edge"], candidate["root"],
        candidate["d_roots"][0], candidate["d_roots"][1],
        candidate["source"], candidate["word"],
    )


def half_signature_compatible(first, second, phase):
    delayed_carry = half.floor_fraction(P * phase + Fraction(1, 2))
    delta_to_second = (2 * (second.carry - delayed_carry) + 1) % P
    delta_to_first = (2 * (first.carry - delayed_carry) + 1) % P
    if not delta_to_second or not delta_to_first:
        return False
    first_d_roots = (
        first.d_roots_left if second.edge == 0 else first.d_roots_right
    )
    second_d_roots = (
        second.d_roots_left if first.edge == 0 else second.d_roots_right
    )
    return (
        any((root + delta_to_second) % P == second.root
            for root in first_d_roots)
        and any((root + delta_to_first) % P == first.root
                for root in second_d_roots)
    )


def half_signature_conditions(first, second, phase):
    delayed_carry = half.floor_fraction(P * phase + Fraction(1, 2))
    delta_to_second = (2 * (second.carry - delayed_carry) + 1) % P
    delta_to_first = (2 * (first.carry - delayed_carry) + 1) % P
    delta_ok = bool(delta_to_second and delta_to_first)
    first_d_roots = (
        first.d_roots_left if second.edge == 0 else first.d_roots_right
    )
    second_d_roots = (
        second.d_roots_left if first.edge == 0 else second.d_roots_right
    )
    first_root_ok = bool(delta_to_second) and any(
        (root + delta_to_second) % P == second.root
        for root in first_d_roots
    )
    second_root_ok = bool(delta_to_first) and any(
        (root + delta_to_first) % P == first.root
        for root in second_d_roots
    )
    return delta_ok, first_root_ok, second_root_ok


def census_half_pairs(candidates):
    groups = defaultdict(Counter)
    examples = {}
    displayed = {}
    for candidate in candidates:
        key = (candidate["phase"], candidate["shallow"], candidate["owner"])
        signature = half_signature(candidate)
        groups[key][signature] += 1
        examples.setdefault((key, signature), candidate)
        displayed[candidate["x"]] = candidate

    # Full reciprocal-pair census, counted once by shallow<owner.
    cycle_count = 0
    source_pair_counts = Counter()
    participating = set()
    witness = None
    both_e3_witness = None
    stage_counts = Counter()
    e3_stage_counts = Counter()
    e3_stage_participation = defaultdict(set)
    for phase in half.FIXED_PHASES:
        for shallow in range(7):
            for owner in range(shallow + 1, 7):
                left_key = (phase, shallow, owner)
                right_key = (phase, owner, shallow)
                for left_signature, left_count in groups[left_key].items():
                    for right_signature, right_count in groups[right_key].items():
                        multiplicity = left_count * right_count
                        contains_e3 = (
                            left_signature.source == 3
                            or right_signature.source == 3
                        )
                        delta_ok, left_root_ok, right_root_ok = (
                            half_signature_conditions(
                                left_signature, right_signature, phase
                            )
                        )
                        stages = {
                            "reverse_clock": True,
                            "nonzero_deltas": delta_ok,
                            "one_root_leg": (
                                delta_ok and (left_root_ok or right_root_ok)
                            ),
                            "both_root_legs": (
                                delta_ok and left_root_ok and right_root_ok
                            ),
                        }
                        for stage, passed in stages.items():
                            if not passed:
                                continue
                            stage_counts[stage] += multiplicity
                            if contains_e3:
                                e3_stage_counts[stage] += multiplicity
                                if left_signature.source == 3:
                                    e3_stage_participation[stage].add(
                                        (left_key, left_signature)
                                    )
                                if right_signature.source == 3:
                                    e3_stage_participation[stage].add(
                                        (right_key, right_signature)
                                    )
                        if not stages["both_root_legs"]:
                            continue
                        cycle_count += multiplicity
                        source_pair_counts[
                            tuple(sorted((left_signature.source,
                                          right_signature.source)))
                        ] += multiplicity
                        participating.add((left_key, left_signature))
                        participating.add((right_key, right_signature))
                        pair = (
                            examples[left_key, left_signature],
                            examples[right_key, right_signature],
                        )
                        if (left_signature.source == 3
                                or right_signature.source == 3):
                            witness = witness or pair
                        if (left_signature.source == 3
                                and right_signature.source == 3):
                            both_e3_witness = both_e3_witness or pair

    participating_counts = Counter()
    for key, signatures in groups.items():
        for signature, count in signatures.items():
            if (key, signature) in participating:
                participating_counts[signature.source, signature.word] += count

    original_points = (
        Fraction(55_232_507, 24 * R),
        Fraction(58_313_459, 24 * R),
    )
    original_pair = tuple(displayed[value] for value in original_points)
    require(half_signature_compatible(
                half_signature(original_pair[0]),
                half_signature(original_pair[1]), Fraction(11, 24)),
            "THM-2698 displayed pair failed the aggregated pair predicate")
    e3_participating_candidates = {
        stage: sum(
            groups[key][signature]
            for key, signature in signatures
        )
        for stage, signatures in e3_stage_participation.items()
    }
    return {
        "cycle_count": cycle_count,
        "source_pair_counts": source_pair_counts,
        "participating_counts": participating_counts,
        "witness": witness,
        "both_e3_witness": both_e3_witness,
        "original_sources": tuple(row["source"] for row in original_pair),
        "stage_counts": stage_counts,
        "e3_stage_counts": e3_stage_counts,
        "e3_stage_participating_candidates": e3_participating_candidates,
    }


def generate_c221_candidates(module, rails, present, present_starts,
                             pair_prefixes):
    require(c221.old.T % R == 0,
            "carrier grid stopped resolving the C221 fibres")
    unit = c221.old.T // R
    rail_prepared = tuple(c221.prepare_weighted(rail[3]) for rail in rails)
    rails_by_owner = {owner: [] for owner in range(7)}
    for index, rail in enumerate(rails):
        rails_by_owner[rail[1]].append(index)
    present_prepared = {
        (shallow, h): c221.prepare_unweighted(
            present[shallow, (-h) % P]
        )
        for shallow in range(7) for h in range(P)
    }
    filter_counts = Counter()
    reroot_counts = Counter()
    raw_e3_witness = None

    # The literal-root generator below fixes one predecessor carry on each
    # side.  Audit the entire D^6 fibre first so that this re-rooting choice
    # is visible as a separate filter rather than baked into "raw".
    for a in (4, 13):
        phase = Fraction(a, 17)
        for predecessor in range(R):
            numerator = 17 * predecessor + a
            source, word = semantic_record_numerator(numerator, L)
            filter_counts["full_fibre", phase, source] += 1
            if source == 3 and raw_e3_witness is None:
                raw_e3_witness = (Fraction(numerator, L), phase, word)
            if source == 3:
                reroot_counts[
                    "E3_predecessor_carry", phase, predecessor % P
                ] += 1

    def geometry(a, carry, h, kappa):
        phase = Fraction(a, 17)
        y_grid = Fraction(a * c221.old.T, 17)
        sector_ok = {
            shallow: c221.strict_in_unweighted(
                y_grid,
                c221.old.prefix_intervals(
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
            source, word = semantic_record_numerator(numerator, L)
            filter_counts["raw_forced_carry", phase, source] += 1
            shallow = c221.clock_from_stalk_numerator(numerator, P)
            owner = c221.clock_from_stalk_numerator(numerator, P**2)
            if shallow == owner or not sector_ok[shallow]:
                continue
            filter_counts["nonconstant_clock_sector", phase, source] += 1
            if not c221.inside_prepared_unweighted(
                    point_scaled_17,
                    present_prepared[shallow, h]):
                continue
            filter_counts["present", phase, source] += 1
            active = tuple(
                index for index in rails_by_owner[owner]
                if c221.inside_prepared_weighted(
                    point_scaled_17, rail_prepared[index]
                )
            )
            if not active:
                continue
            filter_counts["rail", phase, source] += 1
            rows.append((numerator, shallow, owner, active))
            combinations.update((index, shallow) for index in active)
        return rows, combinations

    parameters = ((4, 3, 3, 0), (13, 9, 9, 1))
    sides = []
    geometry_counts = []
    semantic_counts = Counter()
    for a, carry, h, kappa in parameters:
        phase = Fraction(a, 17)
        require(terminal_word(phase, 3) == (1, 2),
                "C221 delayed phase left the deepest-owner fork")
        rows, combinations = geometry(a, carry, h, kappa)
        geometry_counts.append(len(rows))
        unit_table = {}
        for rail_index, shallow in combinations:
            root, _vector, determinant = c221.old.unit_vector(
                module, pair_prefixes, present, present_starts,
                rails[rail_index], 0, 1, carry, h, kappa,
            )
            require(root == 6,
                    "C221 literal-root fibre acquired the wrong root")
            unit_table[rail_index, shallow] = determinant
        typed = []
        for numerator, shallow, owner, active in rows:
            unit_indices = tuple(
                index for index in active
                if unit_table[index, shallow]
            )
            unit_sources = tuple(sorted(
                {rails[index][0] for index in unit_indices}
            ))
            if not unit_sources:
                continue
            value = Fraction(numerator, L)
            source, word = semantic_record_numerator(numerator, L)
            candidate = {
                "x": value,
                "numerator": numerator,
                "phase": phase,
                "shallow": shallow,
                "owner": owner,
                "carry": carry,
                "h": h,
                "kappa": kappa,
                "root": 6,
                "d_root": 6,
                "active_indices": active,
                "active_sources": tuple(sorted(
                    {rails[index][0] for index in active}
                )),
                "unit_indices": unit_indices,
                "unit_sources": unit_sources,
                "source": source,
                "word": word,
            }
            typed.append(candidate)
            semantic_counts[phase, source, word] += 1
            filter_counts["primitive_unit", phase, source] += 1
        sides.append(typed)

    require(tuple(geometry_counts) == (6_260, 6_252),
            "C221 pre-unit candidate universe changed")
    require(tuple(map(len, sides)) == (6_260, 5_966),
            "C221 primitive-unit candidate universe changed")
    require(all(word == (1, 2)
                for (_phase, source, word) in semantic_counts
                if source == 3),
            "C221 E3 row lost its prescribed deepest fork")
    require(raw_e3_witness is not None
            and raw_e3_witness[2] == (1, 2),
            "C221 full fibre lost the E3-to-fork intersection")
    return (
        sides, semantic_counts, filter_counts,
        reroot_counts, raw_e3_witness,
    )


def census_c221_pairs(sides):
    following_by_edge = defaultdict(list)
    for index, candidate in enumerate(sides[1]):
        following_by_edge[
            candidate["shallow"], candidate["owner"]
        ].append((index, candidate))

    pair_count = 0
    paired_initial_states = 0
    source_pair_counts = Counter()
    participating_left = set()
    participating_right = set()
    witness = None
    both_e3_witness = None
    stage_counts = Counter()
    e3_stage_counts = Counter()
    e3_stage_participation = defaultdict(set)
    for left_index, current in enumerate(sides[0]):
        has_pair = False
        for right_index, following in following_by_edge[
                current["owner"], current["shallow"]]:
            contains_e3 = (
                current["source"] == 3 or following["source"] == 3
            )
            stages = {
                "reverse_clock": True,
                "common_active_rail_source": bool(
                    set(current["active_sources"])
                    .intersection(following["active_sources"])
                ),
                "common_primitive_unit_source": bool(
                    set(current["unit_sources"])
                    .intersection(following["unit_sources"])
                ),
            }
            for stage, passed in stages.items():
                if not passed:
                    continue
                stage_counts[stage] += 1
                if contains_e3:
                    e3_stage_counts[stage] += 1
                    if current["source"] == 3:
                        e3_stage_participation[stage].add((0, left_index))
                    if following["source"] == 3:
                        e3_stage_participation[stage].add((1, right_index))
            common = tuple(sorted(
                set(current["unit_sources"])
                .intersection(following["unit_sources"])
            ))
            if not common:
                continue
            has_pair = True
            pair_count += 1
            participating_left.add(left_index)
            participating_right.add(right_index)
            source_pair_counts[
                tuple(sorted((current["source"], following["source"])))
            ] += 1
            pair = (current, following, common)
            if current["source"] == 3 or following["source"] == 3:
                witness = witness or pair
            if current["source"] == following["source"] == 3:
                both_e3_witness = both_e3_witness or pair
        paired_initial_states += int(has_pair)

    require((paired_initial_states, pair_count) == (3_441, 902_300),
            "C221 reciprocal pair universe changed")
    e3_stage_participating_candidates = {}
    for stage, indices in e3_stage_participation.items():
        e3_stage_participating_candidates[stage] = len(indices)
    participating_counts = Counter()
    for index, candidate in enumerate(sides[0]):
        if index in participating_left:
            participating_counts[
                candidate["source"], candidate["word"]
            ] += 1
    for index, candidate in enumerate(sides[1]):
        if index in participating_right:
            participating_counts[
                candidate["source"], candidate["word"]
            ] += 1
    original = {
        candidate["numerator"]: candidate
        for side in sides for candidate in side
        if candidate["numerator"] in c221.POINT_NUMERATORS
    }
    require(tuple(original[numerator]["source"]
                  for numerator in c221.POINT_NUMERATORS) == (0, 0),
            "displayed C221 source failure changed")
    return {
        "pair_count": pair_count,
        "paired_initial_states": paired_initial_states,
        "source_pair_counts": source_pair_counts,
        "participating_counts": participating_counts,
        "witness": witness,
        "both_e3_witness": both_e3_witness,
        "stage_counts": stage_counts,
        "e3_stage_counts": e3_stage_counts,
        "e3_stage_participating_candidates": (
            e3_stage_participating_candidates
        ),
    }


def generate_c221_e3_reroot_candidates(module, rails, present,
                                        present_starts, pair_prefixes):
    """Reopen the two E3 carry classes and their extremal roots 1/12."""
    unit = c221.old.T // R
    rail_prepared = tuple(c221.prepare_weighted(rail[3]) for rail in rails)
    rails_by_owner = {owner: [] for owner in range(7)}
    for index, rail in enumerate(rails):
        rails_by_owner[rail[1]].append(index)
    present_prepared = {
        (shallow, h): c221.prepare_unweighted(
            present[shallow, (-h) % P]
        )
        for shallow in range(7) for h in range(P)
    }
    parameters = (
        (4, 3, 0, (0, 6)),
        (13, 9, 1, (6, 12)),
    )
    sides = []
    stage_counts = Counter()
    root_counts = Counter()

    for a, h, kappa, carries in parameters:
        phase = Fraction(a, 17)
        y_grid = Fraction(a * c221.old.T, 17)
        sector_ok = {
            shallow: c221.strict_in_unweighted(
                y_grid,
                c221.old.prefix_intervals(
                    pair_prefixes[0][shallow][h][kappa]
                ),
            )
            for shallow in range(7)
        }
        geometry = []
        combinations = set()
        for carry in carries:
            for predecessor in range(carry, R, P):
                numerator = 17 * predecessor + a
                source, word = semantic_record_numerator(numerator, L)
                if source != 3:
                    continue
                require(word == (1, 2),
                        "C221 re-root E3 row lost its fork")
                stage_counts["E3_carry", phase, carry] += 1
                point_scaled_17 = numerator * unit
                shallow = c221.clock_from_stalk_numerator(numerator, P)
                owner = c221.clock_from_stalk_numerator(numerator, P**2)
                if shallow == owner or not sector_ok[shallow]:
                    continue
                stage_counts["clock_sector", phase, carry] += 1
                if not c221.inside_prepared_unweighted(
                        point_scaled_17,
                        present_prepared[shallow, h]):
                    continue
                stage_counts["present", phase, carry] += 1
                active = tuple(
                    index for index in rails_by_owner[owner]
                    if c221.inside_prepared_weighted(
                        point_scaled_17, rail_prepared[index]
                    )
                )
                if not active:
                    continue
                stage_counts["rail", phase, carry] += 1
                geometry.append(
                    (numerator, shallow, owner, carry, active)
                )
                combinations.update(
                    (index, shallow, carry, edge)
                    for index in active for edge in (0, 1)
                    if (
                        2 * carry + (2 * h + kappa) // P
                        + int(edge == 0)
                    ) % P
                )

        unit_table = {}
        for rail_index, shallow, carry, edge in combinations:
            root, _vector, determinant = c221.old.unit_vector(
                module, pair_prefixes, present, present_starts,
                rails[rail_index], 0, edge, carry, h, kappa,
            )
            unit_table[rail_index, shallow, carry, edge] = (
                root, determinant
            )

        typed = []
        for numerator, shallow, owner, carry, active in geometry:
            value = Fraction(numerator, L)
            d_value = frac(P * value)
            d_memberships = tuple(
                (edge, root)
                for edge in (0, 1)
                for root in half.half_roots(module, d_value, edge)
            )
            if not d_memberships:
                continue
            stage_counts["D_private_root", phase, carry] += 1
            for edge in (0, 1):
                root = (
                    2 * carry + (2 * h + kappa) // P
                    + int(edge == 0)
                ) % P
                if (root == 0
                        or not half.half_membership(
                            module, value, edge, root)):
                    continue
                stage_counts[
                    "current_private_root", phase, carry, edge, root
                ] += 1
                unit_indices = tuple(
                    index for index in active
                    if unit_table.get(
                        (index, shallow, carry, edge), (root, 0)
                    )[1]
                )
                unit_sources = tuple(sorted(
                    {rails[index][0] for index in unit_indices}
                ))
                if not unit_sources:
                    continue
                stage_counts[
                    "primitive_unit", phase, carry, edge, root
                ] += 1
                root_counts[phase, carry, edge, root] += 1
                typed.append({
                    "x": value,
                    "numerator": numerator,
                    "phase": phase,
                    "shallow": shallow,
                    "owner": owner,
                    "carry": carry,
                    "h": h,
                    "kappa": kappa,
                    "edge": edge,
                    "root": root,
                    "d_memberships": d_memberships,
                    "active_indices": active,
                    "unit_indices": unit_indices,
                    "unit_sources": unit_sources,
                    "source": 3,
                    "word": (1, 2),
                })
        sides.append(typed)

    require(all(candidate["root"] in (1, 12)
                for side in sides for candidate in side),
            "C221 E3 re-root acquired a non-extremal root")
    expected_stage = {
        ("E3_carry", Fraction(4, 17), 0): 65_653,
        ("E3_carry", Fraction(4, 17), 6): 66_099,
        ("E3_carry", Fraction(13, 17), 6): 66_099,
        ("E3_carry", Fraction(13, 17), 12): 65_653,
        ("clock_sector", Fraction(4, 17), 0): 53_227,
        ("clock_sector", Fraction(4, 17), 6): 53_590,
        ("clock_sector", Fraction(13, 17), 6): 53_590,
        ("clock_sector", Fraction(13, 17), 12): 53_227,
    }
    require(dict(stage_counts) == expected_stage,
            "C221 variable-root E3 filter ledger changed")
    require(not any(key[0] == "present" for key in stage_counts)
            and not any(sides) and not root_counts,
            "C221 variable-root E3 unexpectedly passed present typing")
    return sides, stage_counts, root_counts


def census_c221_e3_reroot_pairs(sides):
    following_by_edge = defaultdict(list)
    for index, candidate in enumerate(sides[1]):
        following_by_edge[
            candidate["shallow"], candidate["owner"]
        ].append((index, candidate))
    stage_counts = Counter()
    pair_count = 0
    participating_left = set()
    participating_right = set()
    witness = None
    for left_index, current in enumerate(sides[0]):
        for right_index, following in following_by_edge[
                current["owner"], current["shallow"]]:
            stage_counts["reverse_clock"] += 1
            common = tuple(sorted(
                set(current["unit_sources"])
                .intersection(following["unit_sources"])
            ))
            if not common:
                continue
            stage_counts["common_primitive_unit_source"] += 1
            # The enlarged C221 state is the lawful root sidecar.  Its exact
            # affine transport follows from N' = 13N+m modulo 17R and is
            # checked explicitly for every retained pair.
            lift0 = (
                following["numerator"] - P * current["numerator"]
            ) % L
            lift1 = (
                current["numerator"] - P * following["numerator"]
            ) % L
            for first, second, lift in (
                    (current, following, lift0),
                    (following, current, lift1)):
                require(
                    (P * (2 * first["numerator"] % 221)
                     + 2 * lift) % 221
                    == 2 * second["numerator"] % 221,
                    "C221 re-root microphase transport failed",
                )
            pair_count += 1
            participating_left.add(left_index)
            participating_right.add(right_index)
            witness = witness or (current, following, common, lift0, lift1)
    return {
        "pair_count": pair_count,
        "stage_counts": stage_counts,
        "participating_counts": (
            len(participating_left), len(participating_right)
        ),
        "witness": witness,
    }


def generalized_c221_state_radius(module, pair_prefixes, present, rail,
                                   candidate):
    """Strict x-radius for one variable-edge C221 packet state."""
    value = candidate["x"]
    point_grid = value * c221.old.T
    predecessor, carry, h, kappa, y = c221.state_digits(value)
    require((carry, h, kappa) == (
                candidate["carry"], candidate["h"], candidate["kappa"]),
            "C221 re-root digit state changed")
    rows = []

    left, right = c221.component_weighted(point_grid, rail[3])
    rows.append((Fraction(left, c221.old.T) - value,
                 Fraction(right, c221.old.T) - value))
    word = present[candidate["shallow"], (-h) % P]
    left, right = c221.component_unweighted(point_grid, word)
    rows.append((Fraction(left, c221.old.T) - value,
                 Fraction(right, c221.old.T) - value))
    rows.append((Fraction(predecessor, R) - value,
                 Fraction(predecessor + 1, R) - value))
    half_digit = 2 * h + kappa
    rows.append((
        Fraction(26 * predecessor + half_digit, 26 * R) - value,
        Fraction(26 * predecessor + half_digit + 1, 26 * R) - value,
    ))
    delayed = c221.old.prefix_intervals(
        pair_prefixes[0][candidate["shallow"]][h][kappa]
    )
    left, right = c221.component_unweighted(y * c221.old.T, delayed)
    rows.append((
        Fraction(predecessor, R) + Fraction(left, c221.old.T * R) - value,
        Fraction(predecessor, R) + Fraction(right, c221.old.T * R) - value,
    ))
    low = (14 * candidate["root"] - 13
           if candidate["edge"] == 0 else 14 * candidate["root"])
    left, right = c221.old.comb_component(
        point_grid, module.C3, low, low + 13
    )
    rows.append((Fraction(left, c221.old.T) - value,
                 Fraction(right, c221.old.T) - value))
    d_edge, d_root = candidate["d_memberships"][0]
    d_value = frac(P * value)
    d_grid = d_value * c221.old.T
    d_low = 14 * d_root - 13 if d_edge == 0 else 14 * d_root
    left, right = c221.old.comb_component(
        d_grid, module.C3, d_low, d_low + 13
    )
    rows.append(((Fraction(left, c221.old.T) - d_value) / P,
                 (Fraction(right, c221.old.T) - d_value) / P))
    for speed in (P, P**2):
        left, right = c221.clock_component(value, speed)
        rows.append((left - value, right - value))
    require(all(left < 0 < right for left, right in rows),
            "C221 re-root packet state is not strict")
    return min(min(-left, right) for left, right in rows)


def c221_e3_reroot_witness_ledger(module, rails, present,
                                  pair_prefixes, pair):
    first, second, common, lift0, lift1 = pair
    radii = []
    selected_rails = []
    for candidate in (first, second):
        source = common[0]
        rail_index = next(
            index for index in candidate["unit_indices"]
            if rails[index][0] == source
        )
        selected_rails.append(rail_index)
        radii.append(generalized_c221_state_radius(
            module, pair_prefixes, present, rails[rail_index], candidate
        ))
        require(semantic_record(candidate["x"]) == (3, (1, 2)),
                "C221 re-root witness lost its semantic cospan")
    semantic_radii = tuple(
        semantic_cospan_radius(candidate["x"], 3, (1, 2))
        for candidate in (first, second)
    )
    cycle_radius = min(
        semantic_radii[0] / P**2,
        semantic_radii[1] / P,
        radii[0] / P**2,
        radii[1] / P,
    )
    require(cycle_radius > 0,
            "C221 E3 re-root cycle lost positive width")
    return {
        "first": first,
        "second": second,
        "common_unit_sources": common,
        "selected_rails": tuple(selected_rails),
        "lifts": (lift0, lift1),
        "packet_radii": tuple(radii),
        "semantic_radii": semantic_radii,
        "one_cycle_radius": cycle_radius,
    }


def half_witness_ledger(module, rails, rail_starts, present,
                        present_starts, pair_prefixes, pair):
    first, second = pair
    if first["source"] != 3:
        first, second = second, first
    require(first["source"] == 3 and first["word"] == (1, 2),
            "half witness is not an E3-to-fork event")
    phase = first["phase"]
    delayed_carry = half.floor_fraction(P * phase + Fraction(1, 2))
    k0 = (second["N"] - P * first["N"] - delayed_carry) % R
    k1 = (first["N"] - P * second["N"] - delayed_carry) % R
    require(half.half_handoff(first["x"], k0) == second["x"]
            and half.half_handoff(second["x"], k1) == first["x"],
            "half E3 witness affine cycle failed")
    require(ordinary_dilate(first["x"], 6) == phase,
            "half semantic chronology stopped at the wrong endpoint")
    projected_next = frac(P * phase + Fraction(1, 2))
    require(ordinary_dilate(second["x"], 6) == projected_next == phase,
            "half affine edge stopped commuting with delayed projection")
    semantic_radius = semantic_cospan_radius(first["x"], 3, (1, 2))
    first_slack, _ = half.packet_slack(
        module, rails, rail_starts, present, present_starts,
        pair_prefixes, first,
    )
    second_slack, _ = half.packet_slack(
        module, rails, rail_starts, present, present_starts,
        pair_prefixes, second,
    )
    cycle_radius = min(
        semantic_radius, first_slack / P**2, second_slack / P
    )
    require(cycle_radius > 0,
            "half semantic-attached cycle lost positive width")
    return {
        "first": first,
        "second": second,
        "lifts": (k0, k1),
        "semantic_radius": semantic_radius,
        "one_cycle_radius": cycle_radius,
    }


def c221_witness_ledger(module, rails, present, present_starts,
                        pair_prefixes, pair):
    left, right, common_sources = pair
    start, following = (left, right)
    if start["source"] != 3:
        start, following = following, start
    require(start["source"] == 3 and start["word"] == (1, 2),
            "C221 witness is not an E3-to-fork event")
    lift0 = (following["numerator"] - P * start["numerator"]) % L
    lift1 = (start["numerator"] - P * following["numerator"]) % L
    require(frac(P * start["x"] + Fraction(lift0, L)) == following["x"]
            and frac(P * following["x"] + Fraction(lift1, L)) == start["x"],
            "C221 E3 witness affine cycle failed")
    semantic_endpoint = ordinary_dilate(start["x"], 6)
    require(semantic_endpoint == start["phase"]
            and terminal_word(semantic_endpoint, 3) == (1, 2),
            "C221 semantic chronology stopped at the wrong endpoint")
    projected_next = frac(P * start["phase"] + Fraction(lift0 % 17, 17))
    require(projected_next == following["phase"],
            "C221 affine edge stopped commuting with delayed projection")

    selected_source = common_sources[0]
    constraints = []
    for event, candidate in enumerate((start, following)):
        rail_index = next(
            index for index in candidate["unit_indices"]
            if rails[index][0] == selected_source
        )
        rail = rails[rail_index]
        state = (
            candidate["shallow"], candidate["owner"], rail[2],
            candidate["carry"], candidate["h"], candidate["kappa"],
            candidate["root"], candidate["d_root"],
        )
        rows = c221.state_constraints(
            module, pair_prefixes, present, rail, candidate["x"], state
        )
        constraints.append(rows)
    eta = tuple(
        min(min(-left_bound, right_bound)
            for left_bound, right_bound, _name in rows)
        for rows in constraints
    )
    semantic_radius = semantic_cospan_radius(start["x"], 3, (1, 2))
    cycle_radius = min(semantic_radius, eta[0] / P**2, eta[1] / P)
    require(cycle_radius > 0,
            "C221 semantic-attached cycle lost positive width")
    return {
        "start": start,
        "following": following,
        "common_unit_sources": common_sources,
        "selected_unit_source": selected_source,
        "lifts": (lift0, lift1),
        "semantic_radius": semantic_radius,
        "packet_radii": eta,
        "one_cycle_radius": cycle_radius,
    }


def sorted_counter(counter):
    return tuple(sorted(counter.items(), key=lambda row: repr(row[0])))


def main():
    require((half.P, half.R, c221.P, c221.R, c221.L)
            == (P, R, P, R, L), "inherited carrier scales changed")
    (module, _prefixes, _whole_prefixes, _digit_masses, rails,
     present, present_starts) = half.core.build_carrier_data()
    pair_prefixes = half.private.build_pair_prefixes(module)

    (half_candidates, half_semantics, rail_starts, half_filters,
     half_reroots, half_raw_e3) = generate_half_candidates(
         module, rails, present, present_starts, pair_prefixes
     )
    half_pairs = census_half_pairs(half_candidates)
    half_pair = half_pairs["both_e3_witness"] or half_pairs["witness"]
    half_ledger = None
    if half_pair is not None:
        half_ledger = half_witness_ledger(
            module, rails, rail_starts, present, present_starts,
            pair_prefixes, half_pair,
        )

    (c221_sides, c221_semantics, c221_filters,
     c221_reroots, c221_raw_e3) = generate_c221_candidates(
         module, rails, present, present_starts, pair_prefixes
     )
    c221_pairs = census_c221_pairs(c221_sides)
    c221_pair = c221_pairs["both_e3_witness"] or c221_pairs["witness"]
    c221_ledger = None
    if c221_pair is not None:
        c221_ledger = c221_witness_ledger(
            module, rails, present, present_starts,
            pair_prefixes, c221_pair
        )
    c221_e3_sides, c221_e3_filters, c221_e3_roots = (
        generate_c221_e3_reroot_candidates(
            module, rails, present, present_starts, pair_prefixes
        )
    )
    c221_e3_pairs = census_c221_e3_reroot_pairs(c221_e3_sides)
    c221_e3_ledger = None
    if c221_e3_pairs["witness"] is not None:
        c221_e3_ledger = c221_e3_reroot_witness_ledger(
            module, rails, present, pair_prefixes,
            c221_e3_pairs["witness"],
        )

    require(half_pairs["original_sources"] == (0, 0),
            "displayed half source failure changed")
    require(c221_e3_pairs["pair_count"] == 0
            and c221_e3_pairs["witness"] is None,
            "C221 E3 variable-root stopping boundary changed")
    half_e3_cospan = half_raw_e3
    c221_e3_cospan = c221_raw_e3
    half_e3_cospan_radius = semantic_cospan_radius(
        half_e3_cospan[0], 3, (1, 2)
    )
    c221_e3_cospan_radius = semantic_cospan_radius(
        c221_e3_cospan[0], 3, (1, 2)
    )

    print("LRC14 HALF/C221 FULL PACKET-FIBRE SEMANTIC SOURCE CENSUS")
    print("status=FINITE-EXACT STOPPING/RE-ROOT CERTIFICATE; canonical typed row")
    print(f"row={(GUARD,) + UNITS + BLOCKERS}; clocks={CLOCKS}")
    print("semantic_chronology=E_j at x -> ordinary D^(2j)(x) in "
          "Q_(j,word); affine cycle edges are separate arrows")
    print(f"half_candidate_semantics={sorted_counter(half_semantics)}")
    print(f"half_unique_filter_counts={sorted_counter(half_filters)}")
    print(f"half_E3_reroot_counts={sorted_counter(half_reroots)}")
    print("half_E3_present_union=EMPTY across all 13 present labels; "
          "all E3 rows retain rail+private-root support when present is "
          "removed")
    print(f"half_reciprocal_cycle_count={half_pairs['cycle_count']}")
    print("half_cycle_source_pairs="
          f"{sorted_counter(half_pairs['source_pair_counts'])}")
    print("half_participating_semantics="
          f"{sorted_counter(half_pairs['participating_counts'])}")
    print(f"half_displayed_sources={half_pairs['original_sources']}")
    print(f"half_pair_stage_counts={sorted_counter(half_pairs['stage_counts'])}")
    print("half_E3_pair_stage_counts="
          f"{sorted_counter(half_pairs['e3_stage_counts'])}")
    print("half_E3_stage_participating_candidates="
          f"{tuple(sorted(half_pairs['e3_stage_participating_candidates'].items()))}")
    print("half_E3_witness=(first,second,lifts,semantic_radius,"
          f"one_cycle_radius)={half_ledger}")
    print("half_strict_E3_to_fork_cospan=(candidate,radius)="
          f"{(half_e3_cospan, half_e3_cospan_radius)}")
    print(f"C221_candidate_semantics={sorted_counter(c221_semantics)}")
    print(f"C221_unique_filter_counts={sorted_counter(c221_filters)}")
    print(f"C221_E3_reroot_counts={sorted_counter(c221_reroots)}")
    print(f"C221_reciprocal_pair_count={c221_pairs['pair_count']}; "
          f"paired_initial_states={c221_pairs['paired_initial_states']}")
    print("C221_cycle_source_pairs="
          f"{sorted_counter(c221_pairs['source_pair_counts'])}")
    print("C221_participating_semantics="
          f"{sorted_counter(c221_pairs['participating_counts'])}")
    print(f"C221_pair_stage_counts={sorted_counter(c221_pairs['stage_counts'])}")
    print("C221_E3_pair_stage_counts="
          f"{sorted_counter(c221_pairs['e3_stage_counts'])}")
    print("C221_E3_stage_participating_candidates="
          f"{tuple(sorted(c221_pairs['e3_stage_participating_candidates'].items()))}")
    print("C221_E3_witness=(start,following,common/lifts/radii)="
          f"{c221_ledger}")
    print("C221_variable_root_E3_filter_counts="
          f"{sorted_counter(c221_e3_filters)}")
    print("C221_variable_root_E3_candidate_counts="
          f"{sorted_counter(c221_e3_roots)}")
    print("C221_extremal_root_law=phase4/17:carries(0,6)->"
          "(left1,right12);phase13/17:carries(6,12)->"
          "(left1,right12);the complementary edge has root0")
    print("C221_variable_root_E3_pair_census="
          f"{c221_e3_pairs['pair_count'], sorted_counter(c221_e3_pairs['stage_counts']), c221_e3_pairs['participating_counts']}")
    print("C221_variable_root_E3_witness="
          f"{c221_e3_ledger}")
    print("C221_strict_E3_to_fork_cospan=(candidate,radius)="
          f"{(c221_e3_cospan, c221_e3_cospan_radius)}")
    print("bridge_verdict=both raw D^6 fibres contain strict E3-to-fork "
          "semantic cospans; inherited half/C221 root-6 candidate universes "
          "contain no E3, and the complete extremal-root 1/12 C221 re-root "
          "also dies exactly at the present factor")
    print("THM2305_vs_sidecars=only E3, ordinary D^6, and the terminal fork "
          "are required by the semantic cospan; rail/present/clock/carry/"
          "root/unit conditions are auxiliary local-cycle compatibility "
          "sidecars and may be re-rooted without changing THM2305")
    print("attachment_scope=the ordinary E3-to-D^6-fork cospan is strict "
          "and positive, but no retained affine-cycle state shares its E3 "
          "source; half first loses E3 at present typing, inherited C221 "
          "first loses it at forced carry/root6, and all-carry C221 again "
          "loses it at present typing")
    print("SCOPE: canonical typed row and declared packet fibres only; no "
          "uniform row reduction, endpoint current, row exclusion, or LRC14")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
