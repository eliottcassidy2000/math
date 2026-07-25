#!/usr/bin/env python3
"""Exact finite audit for THM-2192's scalar 5+3 chord invoice.

The proof reduces a generic guard-safe thirteen-root fibre to one of three
anchored chord carriers on F_13 with four consecutive forbidden vertices:

M: one singleton plus a perfect matching of the other eight safe vertices;
X: one safe-to-forbidden chord plus a matching of the other eight safe ones;
D: five safe-safe chords with one double-owned safe vertex.

This script enumerates all three carriers, records their unoriented chord
length multisets, and checks the eight excluded profiles.  It separately
checks the one-root/two-root count and the H*q^{-1} chord-step law for every
nonzero residue pair.  It uses no floating point and no Python ``assert``.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations, combinations_with_replacement


P = 13
FORBIDDEN = frozenset({0, 1, 2, 3})
SAFE = frozenset(range(P)) - FORBIDDEN


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
    return checks + 1


def main() -> None:
    raw_counts, profile_sets = enumerate_carriers()
    monomer_profiles, cross_profiles, double_profiles = profile_sets
    all_profiles = set(combinations_with_replacement(range(1, 7), 5))
    allowed = monomer_profiles | cross_profiles | double_profiles
    excluded = tuple(sorted(all_profiles - allowed))

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

    require(raw_counts == (5670, 3780, 3780), "raw carrier counts changed")
    require(len(all_profiles) == 252, "ambient profile count changed")
    require(
        tuple(map(len, profile_sets)) == (244, 244, 231),
        "individual profile counts changed",
    )
    require(len(allowed) == 244, "allowed profile union changed")
    require(excluded == expected_excluded, "excluded profile list changed")
    arithmetic_checks = root_step_audit()
    require(arithmetic_checks == 433, "root arithmetic check count changed")

    print("THM-2192 scalar five-plus-three root-sheet chord audit")
    print(f"raw_carrier_states=M:{raw_counts[0]},X:{raw_counts[1]},D:{raw_counts[2]}")
    print(
        "carrier_profile_counts="
        f"M:{len(monomer_profiles)},X:{len(cross_profiles)},"
        f"D:{len(double_profiles)}"
    )
    print(f"all_profiles={len(all_profiles)}; allowed={len(allowed)}; excluded={len(excluded)}")
    print("excluded_profiles=" + repr(excluded))
    print(f"root_sheet_arithmetic_checks={arithmetic_checks}")
    print("PASS")


if __name__ == "__main__":
    main()
