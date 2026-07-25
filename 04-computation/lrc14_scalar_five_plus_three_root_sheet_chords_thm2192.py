#!/usr/bin/env python3
"""Exact finite audit for THM-2192's scalar 5+3 chord invoice.

The proof reduces a generic guard-safe thirteen-root fibre to one of three
anchored chord carriers on F_13 with four consecutive forbidden vertices:

M: one singleton plus a perfect matching of the other eight safe vertices;
X: one safe-to-forbidden chord plus a matching of the other eight safe ones;
D: five safe-safe chords with one double-owned safe vertex.

This script enumerates all three carriers, records their unoriented chord
length multisets, and checks the eight excluded profiles.  It also enumerates
the guard-danger carrier: five chords perfectly match the ten vertices outside
an anchored three-consecutive block, producing 216 of the 252 length profiles
and the exact 36-profile fat-guard fork.  It separately checks the
one-root/two-root count and the H*q^{-1} chord-step law for every nonzero
residue pair.  Finally it audits the odd-image density inequality, the
11/13 divisibility trap, every quotient triple through 64, and all exact
fragmentation constants used to empty the fat-guard fork in the actual
scalar lane.  It also separates the 36 missing profiles into 17 detected
by three endpoint-potential inequalities and 19 genuine holes inside the
convex hull of the matching support.  It uses no floating point and no
Python ``assert``.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations, combinations_with_replacement


P = 13
FORBIDDEN = frozenset({0, 1, 2, 3})
SAFE = frozenset(range(P)) - FORBIDDEN
GUARD_DANGER_FORBIDDEN = frozenset({0, 1, 2})
GUARD_DANGER_SAFE = frozenset(range(P)) - GUARD_DANGER_FORBIDDEN


def require(condition: bool, message: str) -> None:
    """Runtime check which remains active under ``python -O``."""
    if not condition:
        raise RuntimeError(message)


def circle_norm(x: Fraction) -> Fraction:
    """Exact distance of a rational number from the nearest integer."""
    residue = x % 1
    return min(residue, 1 - residue)


def matchings(vertices: tuple[int, ...]):
    """Generate every perfect matching of an even ordered vertex tuple."""
    if not vertices:
        yield ()
        return
    first = vertices[0]
    for index in range(1, len(vertices)):
        second = vertices[index]
        remainder = vertices[1:index] + vertices[index + 1 :]
        for tail in matchings(remainder):
            yield ((first, second),) + tail


def chord_length(edge: tuple[int, int]) -> int:
    """Unoriented cyclic length in {1,...,6} on F_13."""
    left, right = edge
    difference = (right - left) % P
    return min(difference, P - difference)


def profile(edges: tuple[tuple[int, int], ...], extra: int | None = None):
    """Sorted multiset of chord lengths, optionally with a singleton label."""
    lengths = [chord_length(edge) for edge in edges]
    if extra is not None:
        lengths.append(extra)
    return tuple(sorted(lengths))


def enumerate_carriers():
    """Return raw state counts and profile sets for M, X, and D."""
    safe = tuple(sorted(SAFE))
    monomer_profiles = set()
    cross_profiles = set()
    double_profiles = set()
    monomer_states = 0
    cross_states = 0
    double_states = 0

    # M: choose the singleton sheet and its otherwise invisible residue
    # length, then match the other eight safe sheets.
    for singleton in safe:
        remainder = tuple(x for x in safe if x != singleton)
        for matching in matchings(remainder):
            for singleton_length in range(1, 7):
                monomer_states += 1
                monomer_profiles.add(profile(matching, singleton_length))

    # X: choose the unique safe-to-forbidden chord, then match the eight
    # remaining safe sheets.
    for safe_endpoint in safe:
        remainder = tuple(x for x in safe if x != safe_endpoint)
        for forbidden_endpoint in sorted(FORBIDDEN):
            for matching in matchings(remainder):
                cross_states += 1
                edges = ((safe_endpoint, forbidden_endpoint),) + matching
                cross_profiles.add(profile(edges))

    # D: choose the unique double-owned vertex and its two distinct
    # neighbours, then match the remaining six safe sheets.
    for repeated in safe:
        others = tuple(x for x in safe if x != repeated)
        for first, second in combinations(others, 2):
            remainder = tuple(
                x for x in others if x != first and x != second
            )
            for matching in matchings(remainder):
                double_states += 1
                edges = (
                    (repeated, first),
                    (repeated, second),
                ) + matching
                double_profiles.add(profile(edges))

    return (
        (monomer_states, cross_states, double_states),
        (monomer_profiles, cross_profiles, double_profiles),
    )


def enumerate_guard_danger_matchings():
    """Enumerate perfect matchings off one anchored three-vertex block."""
    profiles = set()
    states = 0
    for matching in matchings(tuple(sorted(GUARD_DANGER_SAFE))):
        states += 1
        profiles.add(profile(matching))
    return states, profiles


def profile_counts(length_profile: tuple[int, ...]) -> tuple[int, ...]:
    """Return the six length multiplicities of a chord profile."""
    return tuple(length_profile.count(length) for length in range(1, 7))


def profile_word(word: str) -> tuple[int, ...]:
    """Decode the compact five-digit profile notation used in certificates."""
    require(len(word) == 5, "profile certificate word has wrong length")
    result = tuple(int(digit) for digit in word)
    require(
        tuple(sorted(result)) == result and set(result) <= set(range(1, 7)),
        "profile certificate word is malformed",
    )
    return result


def matching_shadow_audit(
    feasible_profiles: set[tuple[int, ...]],
    missing_profiles: tuple[tuple[int, ...], ...],
):
    """Audit the linear-shadow/nonlinear-hole split of the matching support."""
    vertices = tuple(sorted(GUARD_DANGER_SAFE))
    potentials = (
        # Chord length itself, certified by distance from the anchor 8.
        ((1, 2, 3, 4, 5, 6), (5, 4, 3, 2, 1, 0, 1, 2, 3, 4), 25),
        # The indicator of length two.
        ((0, 1, 0, 0, 0, 0), (0, 0, 1, 1, 0, 0, 1, 1, 0, 0), 4),
        # The mixed length-two/three/four obstruction.
        ((0, 1, 3, 2, 0, 0), (1, 0, 0, 2, 3, 3, 2, 0, 0, 1), 12),
    )

    edge_checks = 0
    for weights, endpoint_potential, total in potentials:
        require(sum(endpoint_potential) == total, "wrong potential total")
        potential = dict(zip(vertices, endpoint_potential, strict=True))
        for left, right in combinations(vertices, 2):
            weight = weights[chord_length((left, right)) - 1]
            require(
                weight <= potential[left] + potential[right],
                "endpoint potential does not dominate its edge weight",
            )
            edge_checks += 1

    def linearly_visible(length_profile: tuple[int, ...]) -> bool:
        counts = profile_counts(length_profile)
        return (
            sum(length_profile) > 25
            or counts[1] > 4
            or counts[1] + 3 * counts[2] + 2 * counts[3] > 12
        )

    visible = tuple(
        profile for profile in missing_profiles if linearly_visible(profile)
    )
    holes = tuple(
        profile for profile in missing_profiles if not linearly_visible(profile)
    )
    expected_holes = tuple(
        map(
            profile_word,
            (
                "11112",
                "11123",
                "11222",
                "11233",
                "12223",
                "12333",
                "13444",
                "14445",
                "15555",
                "22233",
                "22244",
                "23555",
                "24444",
                "24555",
                "25556",
                "34444",
                "35555",
                "44445",
                "45555",
            ),
        )
    )
    require(len(visible) == 17, "linear-visible missing count changed")
    require(holes == expected_holes, "nonlinear matching holes changed")

    # Every hole is an exact convex combination of at most three feasible
    # profiles.  Thus no affine-linear inequality in the six length counts
    # which is valid on the feasible support can separate any one of them.
    convex_certificates = {
        "11112": ((Fraction(3, 4), "11111"), (Fraction(1, 4), "12222")),
        "11123": (
            (Fraction(1, 2), "11111"),
            (Fraction(1, 4), "12222"),
            (Fraction(1, 4), "13333"),
        ),
        "11222": ((Fraction(1, 4), "11111"), (Fraction(3, 4), "12222")),
        "11233": (
            (Fraction(1, 4), "11111"),
            (Fraction(1, 4), "12222"),
            (Fraction(1, 2), "13333"),
        ),
        "12223": ((Fraction(3, 4), "12222"), (Fraction(1, 4), "13333")),
        "12333": ((Fraction(1, 2), "11333"), (Fraction(1, 2), "22333")),
        "13444": (
            (Fraction(3, 20), "11111"),
            (Fraction(1, 4), "13333"),
            (Fraction(3, 5), "44444"),
        ),
        "14445": (
            (Fraction(1, 5), "11111"),
            (Fraction(3, 5), "44444"),
            (Fraction(1, 5), "55555"),
        ),
        "15555": ((Fraction(1, 5), "11111"), (Fraction(4, 5), "55555")),
        "22233": ((Fraction(1, 2), "22223"), (Fraction(1, 2), "22333")),
        "22244": ((Fraction(3, 4), "22224"), (Fraction(1, 4), "44444")),
        "23555": (
            (Fraction(1, 4), "22225"),
            (Fraction(1, 4), "33335"),
            (Fraction(1, 2), "55555"),
        ),
        "24444": ((Fraction(1, 4), "22224"), (Fraction(3, 4), "44444")),
        "24555": (
            (Fraction(1, 4), "22225"),
            (Fraction(1, 5), "44444"),
            (Fraction(11, 20), "55555"),
        ),
        "25556": (
            (Fraction(1, 12), "22225"),
            (Fraction(1, 3), "22666"),
            (Fraction(7, 12), "55555"),
        ),
        "34444": ((Fraction(1, 2), "33444"), (Fraction(1, 2), "44444")),
        "35555": ((Fraction(1, 4), "33335"), (Fraction(3, 4), "55555")),
        "44445": ((Fraction(4, 5), "44444"), (Fraction(1, 5), "55555")),
        "45555": ((Fraction(1, 5), "44444"), (Fraction(4, 5), "55555")),
    }
    require(
        set(map(profile_word, convex_certificates)) == set(holes),
        "convex-certificate universe changed",
    )
    maximum_terms = 0
    maximum_denominator = 1
    for hole_word, terms in convex_certificates.items():
        hole = profile_word(hole_word)
        require(sum(weight for weight, _ in terms) == 1, "weights do not sum to one")
        anchors = tuple(profile_word(anchor_word) for _, anchor_word in terms)
        require(
            all(anchor in feasible_profiles for anchor in anchors),
            "convex certificate uses an infeasible anchor",
        )
        target_counts = profile_counts(hole)
        mixed_counts = tuple(
            sum(
                weight * profile_counts(anchor)[length_index]
                for (weight, _), anchor in zip(terms, anchors, strict=True)
            )
            for length_index in range(6)
        )
        require(mixed_counts == target_counts, "convex certificate is inexact")
        maximum_terms = max(maximum_terms, len(terms))
        maximum_denominator = max(
            maximum_denominator,
            *(weight.denominator for weight, _ in terms),
        )

    return (
        edge_checks,
        visible,
        holes,
        maximum_terms,
        maximum_denominator,
    )


def root_step_audit() -> int:
    """Check root counts and the H*q^{-1} doubleton step for all residues."""
    checks = 0

    # A lower dangerous phase gives one root; either safe sign gives two
    # adjacent roots in terminal coordinates.
    for lower_phase, expected_size in (
        (Fraction(0), 1),
        (Fraction(1, 4), 2),
        (Fraction(-1, 4), 2),
    ):
        for guard in range(1, P):
            guard_inverse = pow(guard, -1, P)
            for terminal in range(1, P):
                terminal_step = terminal * guard_inverse % P
                # Write the terminal root with numerator
                # lower_phase + terminal_step*label, reduced modulo 13.
                roots = tuple(
                    label
                    for label in range(P)
                    if circle_norm(
                        (lower_phase + terminal_step * label) / P
                    )
                    < Fraction(1, 14)
                )
                require(
                    len(roots) == expected_size,
                    "terminal root count changed",
                )
                if expected_size == 2:
                    difference = (roots[1] - roots[0]) % P
                    expected = guard * pow(terminal, -1, P) % P
                    require(
                        difference in {expected, (-expected) % P},
                        "H*q^{-1} chord-step law failed",
                    )
                checks += 1

    # A guard-safe lower phase has four consecutive guard-unsafe roots.
    guard_bad = tuple(
        label
        for label in range(P)
        if circle_norm((Fraction(1, 4) + label) / P)
        < Fraction(1, 7)
    )
    require(len(guard_bad) == 4, "guard root count changed")
    guard_bad_set = frozenset(guard_bad)
    consecutive = any(
        guard_bad_set
        == frozenset((start + offset) % P for offset in range(4))
        for start in range(P)
    )
    require(
        consecutive,
        "guard-unsafe roots are not consecutive",
    )

    # A guard-danger lower phase has three consecutive guard-unsafe roots.
    guard_danger_bad = tuple(
        label
        for label in range(P)
        if circle_norm(Fraction(label, P)) < Fraction(1, 7)
    )
    require(len(guard_danger_bad) == 3, "guard-danger root count changed")
    guard_danger_bad_set = frozenset(guard_danger_bad)
    guard_danger_consecutive = any(
        guard_danger_bad_set
        == frozenset((start + offset) % P for offset in range(3))
        for start in range(P)
    )
    require(
        guard_danger_consecutive,
        "guard-danger roots are not consecutive",
    )
    return checks + 2


def interval_covered_by_dangers(
    speeds: tuple[int, ...],
    left: Fraction,
    right: Fraction,
) -> bool:
    """Decide exact a.e. coverage of a positive interval by closed combs."""
    intervals = []
    for speed in speeds:
        require(speed > 0, "interval audit expects positive speeds")
        for tooth in range(speed + 1):
            lower = Fraction(14 * tooth - 1, 14 * speed)
            upper = Fraction(14 * tooth + 1, 14 * speed)
            if upper < left or lower > right:
                continue
            intervals.append((max(left, lower), min(right, upper)))

    reach = left
    for lower, upper in sorted(intervals):
        if upper < reach:
            continue
        if lower > reach:
            return False
        reach = max(reach, upper)
        if reach >= right:
            return True
    return reach >= right


def fat_guard_obstruction_audit():
    """Audit the exact constants in the odd-guard no-cover theorem."""
    odd_orders = tuple(range(3, 10_000, 2))
    for order in odd_orders:
        require(
            2 * (order // 7 + 1) < order,
            "odd translated two-band density bound failed",
        )

    divisibility_checks = 0
    for modulus in (11, 13):
        for residue in range(modulus):
            in_closed_danger = (
                circle_norm(Fraction(residue, modulus))
                <= Fraction(1, 14)
            )
            require(
                in_closed_danger == (residue == 0),
                "11/13 divisibility test failed",
            )
            divisibility_checks += 1

    whole_interval_threshold = Fraction(8, 7)
    no_unit_maximum = Fraction(1, 2) + Fraction(1, 3) + Fraction(1, 4)
    side_interval_threshold = Fraction(5, 14)
    different_holder_maximum = Fraction(1, 11) + Fraction(1, 13)
    same_holder_other_floor = side_interval_threshold - Fraction(1, 143)
    terminal_fragmentation = Fraction(1, 98) + Fraction(1, 1001)

    require(
        no_unit_maximum < whole_interval_threshold,
        "whole-interval harmonic obstruction failed",
    )
    require(
        different_holder_maximum < side_interval_threshold,
        "different-holder divisibility obstruction failed",
    )
    require(
        same_holder_other_floor > Fraction(1, 3),
        "same-holder branch no longer forces the other quotient to two",
    )
    require(
        terminal_fragmentation < Fraction(1, 14),
        "terminal one-comb fragmentation obstruction failed",
    )

    quotient_triples_checked = 0
    quotient_covers = 0
    for speeds in combinations(range(1, 65), 3):
        quotient_triples_checked += 1
        if interval_covered_by_dangers(
            speeds,
            Fraction(0),
            Fraction(1, 7),
        ):
            quotient_covers += 1
    require(
        quotient_triples_checked == 41_664,
        "quotient hostile-sweep universe changed",
    )
    require(quotient_covers == 0, "quotient hostile sweep found a cover")

    return (
        len(odd_orders),
        divisibility_checks,
        whole_interval_threshold,
        no_unit_maximum,
        side_interval_threshold,
        different_holder_maximum,
        same_holder_other_floor,
        terminal_fragmentation,
        quotient_triples_checked,
        quotient_covers,
    )


def main() -> None:
    raw_counts, profile_sets = enumerate_carriers()
    monomer_profiles, cross_profiles, double_profiles = profile_sets
    guard_danger_states, guard_danger_profiles = (
        enumerate_guard_danger_matchings()
    )
    all_profiles = set(combinations_with_replacement(range(1, 7), 5))
    allowed = monomer_profiles | cross_profiles | double_profiles
    excluded = tuple(sorted(all_profiles - allowed))
    guard_danger_excluded = tuple(
        sorted(all_profiles - guard_danger_profiles)
    )

    expected_excluded = (
        (3, 3, 3, 3, 3),
        (3, 6, 6, 6, 6),
        (4, 5, 6, 6, 6),
        (4, 6, 6, 6, 6),
        (5, 5, 5, 6, 6),
        (5, 5, 6, 6, 6),
        (5, 6, 6, 6, 6),
        (6, 6, 6, 6, 6),
    )
    expected_guard_danger_excluded = (
        (1, 1, 1, 1, 2),
        (1, 1, 1, 2, 3),
        (1, 1, 2, 2, 2),
        (1, 1, 2, 3, 3),
        (1, 2, 2, 2, 3),
        (1, 2, 3, 3, 3),
        (1, 3, 4, 4, 4),
        (1, 4, 4, 4, 5),
        (1, 5, 5, 5, 5),
        (2, 2, 2, 2, 2),
        (2, 2, 2, 3, 3),
        (2, 2, 2, 4, 4),
        (2, 3, 3, 3, 3),
        (2, 3, 5, 5, 5),
        (2, 4, 4, 4, 4),
        (2, 4, 5, 5, 5),
        (2, 5, 5, 5, 6),
        (2, 6, 6, 6, 6),
        (3, 3, 3, 3, 3),
        (3, 3, 3, 3, 4),
        (3, 3, 3, 4, 4),
        (3, 4, 4, 4, 4),
        (3, 5, 5, 5, 5),
        (3, 5, 6, 6, 6),
        (3, 6, 6, 6, 6),
        (4, 4, 4, 4, 5),
        (4, 4, 6, 6, 6),
        (4, 5, 5, 5, 5),
        (4, 5, 5, 6, 6),
        (4, 5, 6, 6, 6),
        (4, 6, 6, 6, 6),
        (5, 5, 5, 5, 6),
        (5, 5, 5, 6, 6),
        (5, 5, 6, 6, 6),
        (5, 6, 6, 6, 6),
        (6, 6, 6, 6, 6),
    )

    require(raw_counts == (5670, 3780, 3780), "raw carrier counts changed")
    require(
        guard_danger_states == 945,
        "guard-danger matching count changed",
    )
    require(len(all_profiles) == 252, "ambient profile count changed")
    require(
        tuple(map(len, profile_sets)) == (244, 244, 231),
        "individual profile counts changed",
    )
    require(len(allowed) == 244, "allowed profile union changed")
    require(excluded == expected_excluded, "excluded profile list changed")
    require(
        len(guard_danger_profiles) == 216,
        "guard-danger profile count changed",
    )
    require(
        guard_danger_excluded == expected_guard_danger_excluded,
        "guard-danger excluded profile list changed",
    )
    require(
        set(excluded) <= set(guard_danger_excluded),
        "unconditional exclusions escaped the guard-danger list",
    )
    arithmetic_checks = root_step_audit()
    require(arithmetic_checks == 434, "root arithmetic check count changed")
    fat_guard_audit = fat_guard_obstruction_audit()
    matching_shadow = matching_shadow_audit(
        guard_danger_profiles,
        guard_danger_excluded,
    )

    print("THM-2192 scalar five-plus-three root-sheet chord audit")
    print(f"raw_carrier_states=M:{raw_counts[0]},X:{raw_counts[1]},D:{raw_counts[2]}")
    print(
        "carrier_profile_counts="
        f"M:{len(monomer_profiles)},X:{len(cross_profiles)},"
        f"D:{len(double_profiles)}"
    )
    print(f"all_profiles={len(all_profiles)}; allowed={len(allowed)}; excluded={len(excluded)}")
    print("excluded_profiles=" + repr(excluded))
    print(
        "guard_danger_matching_states="
        f"{guard_danger_states}; allowed={len(guard_danger_profiles)}; "
        f"excluded={len(guard_danger_excluded)}"
    )
    print(
        "guard_danger_excluded_profiles="
        + repr(guard_danger_excluded)
    )
    print(
        "unconditional_exclusions_subset_guard_danger="
        + repr(set(excluded) <= set(guard_danger_excluded))
    )
    print(f"root_sheet_arithmetic_checks={arithmetic_checks}")
    print(
        "fat_guard_odd_image_orders="
        f"{fat_guard_audit[0]}; divisibility_residues={fat_guard_audit[1]}"
    )
    print(
        "fat_guard_whole_threshold="
        f"{fat_guard_audit[2]}; no_unit_max={fat_guard_audit[3]}"
    )
    print(
        "fat_guard_side_threshold="
        f"{fat_guard_audit[4]}; different_holder_max={fat_guard_audit[5]}; "
        f"same_holder_other_floor={fat_guard_audit[6]}"
    )
    print(
        "fat_guard_terminal_fragmentation="
        f"{fat_guard_audit[7]}; target=1/14"
    )
    print(
        "fat_guard_quotient_triples_checked="
        f"{fat_guard_audit[8]}; covers={fat_guard_audit[9]}"
    )
    print("fat_guard_actual_missing_profiles_empty=36; residual_profiles=216")
    print(
        "matching_linear_shadow="
        f"edge_checks:{matching_shadow[0]}; visible_missing:{len(matching_shadow[1])}; "
        f"convex_support_holes:{len(matching_shadow[2])}"
    )
    print("matching_convex_support_holes=" + repr(matching_shadow[2]))
    print(
        "matching_convex_certificates="
        f"max_terms:{matching_shadow[3]}; max_denominator:{matching_shadow[4]}"
    )
    print("PASS")


if __name__ == "__main__":
    main()
