#!/usr/bin/env python3
"""Exact transports of the proved rank-seven atom-generated monoid.

The seven-base antichain is not known to be complete.  Nevertheless, every
multiple of one of its replayed atoms that avoids the proved rank-at-most-six
bases has zero-mode-cochain rank exactly seven.  This companion computes that
unconditional subfamily's density and its Fibonacci/Berggren transports.
"""

from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from math import gcd, lcm
from pathlib import Path

from lrc_zero_mode_cochain_rank6_support_probe_20260815 import (
    BERGGREN_MATRICES,
    danger_mask,
    matrix_child,
)


ROOT = Path(__file__).resolve().parents[1]
LOWER_BASES = (8, 9, 10, 11, 12, 15, 23, 25)
ATOM_BASES = (13, 14, 29, 38, 51, 68, 148)
EXPECTED_SEMANTIC_DIGEST = "64e9032fc79cf2c5014497baddbbccbc3439650903b25447505f93bc66fb1c26"
PINNED = (
    (
        "THM-3416",
        ROOT / "01-canon/theorems/THM-3416-zero-mode-cochain-global-rank-six-support.md",
        "42a9309145de51d1bb6fca0b7c1945302ff37a63a3183e1dfed838c07118e8bf",
    ),
    (
        "THM-3420",
        ROOT / "01-canon/theorems/THM-3420-prime-rank-seven-zero-and-half-twist-splitter-closures.md",
        "5f307c906c054850f2589b2e9fb7c0e902b3571747a6c61ebc2c2f2628212b2f",
    ),
    (
        "rank-seven-scout-script",
        ROOT / "04-computation/lrc_zero_mode_cochain_rank7_ancestor_scout_20260815.py",
        "bcf80194bd5697b121b08b239f096d67d7e0214bbc81abf4a9fb1f6d939fa560",
    ),
    (
        "rank-seven-scout-output",
        ROOT / "05-knowledge/results/lrc_zero_mode_cochain_rank7_ancestor_scout_20260815.out",
        "95012438ffa2e36183ef3e56a9881c3daa52ddfe118abebfd02106b1b28f65a2",
    ),
)

ATOMS = (
    (13, 1, (1, 2, 3, 5, 7, 9, 11)),
    (14, 1, (1, 3, 4, 5, 9, 11, 13)),
    (29, 0, (1, 4, 5, 6, 7, 9, 13)),
    (29, 1, (1, 5, 7, 8, 12, 13, 22)),
    (38, 1, (1, 9, 17, 20, 21, 29, 37)),
    (51, 1, (1, 11, 12, 18, 23, 34, 35)),
    (68, 1, (8, 11, 23, 24, 45, 56, 57)),
    (148, 1, (8, 33, 41, 100, 107, 115, 140)),
)


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def lf_hash(path):
    return sha256(path.read_bytes().replace(bytes((13, 10)), bytes((10,)))).hexdigest()


def atom_generated_rank_seven(value):
    return (
        any(value % base == 0 for base in ATOM_BASES)
        and not any(value % base == 0 for base in LOWER_BASES)
    )


def union_minus_union_density(positive_bases, negative_bases):
    answer = Fraction(0)
    for positive_size in range(1, len(positive_bases) + 1):
        for positive_packet in combinations(positive_bases, positive_size):
            positive_lcm = lcm(*positive_packet)
            for negative_size in range(len(negative_bases) + 1):
                for negative_packet in combinations(negative_bases, negative_size):
                    modulus = lcm(positive_lcm, *negative_packet) if negative_packet else positive_lcm
                    answer += (-1) ** (positive_size + negative_size + 1) * Fraction(1, modulus)
    return answer


def first_zero(modulus):
    left, right = 0, 1
    for index in range(1, 6 * modulus + 1):
        left, right = right, (left + right) % modulus
        if left == 0:
            return index
    raise RuntimeError((modulus, "no zero in guarded range"))


def atom_replay():
    rows = []
    for q, epsilon, residues in ATOMS:
        full = (1 << q) - 1
        joined = 0
        for residue in residues:
            joined |= danger_mask(q, residue, epsilon)
        require(joined == full, (q, epsilon, residues, joined, full))
        modulus = q if epsilon == 0 else 2 * q
        require(gcd(modulus, *residues) == 1, (q, epsilon, residues))
        require(atom_generated_rank_seven(q), (q, "not in generated exact-rank-seven set"))
        rows.append((q, epsilon, residues))
    return tuple(rows)


def fibonacci_transport(modulus):
    apparition = tuple(
        (base, first_zero(base))
        for base in LOWER_BASES + ATOM_BASES
    )
    expected = (
        (8, 6), (9, 12), (10, 15), (11, 10), (12, 12), (15, 20),
        (23, 24), (25, 25), (13, 7), (14, 24), (29, 14), (38, 18),
        (51, 36), (68, 18), (148, 114),
    )
    require(apparition == expected, apparition)

    period = 1050
    residues = []
    left, right = 0, 1
    for index in range(period):
        require(
            atom_generated_rank_seven(left)
            == (
                index % 7 == 0
                and index % 6 != 0
                and index % 10 != 0
                and index % 15 != 0
                and index % 25 != 0
            ),
            (index, left),
        )
        if atom_generated_rank_seven(left):
            residues.append(index)
        left, right = right, (left + right) % modulus
    require(len(residues) == 108, len(residues))
    return apparition, period, tuple(residues), Fraction(len(residues), period)


def berggren_transport(modulus):
    spine_period = lcm(51, 9, 11)
    require(spine_period == 1683, spine_period)
    roots_51 = tuple(index for index in range(51) if (4 * index * index + 12 * index + 11) % 51 == 0)
    require(roots_51 == (2, 19, 29, 46), roots_51)
    require(not any((4 * index * index + 12 * index + 11) % 13 == 0 for index in range(13)), "unexpected mod13 root")
    require(not any((4 * index * index + 12 * index + 11) % 29 == 0 for index in range(29)), "unexpected mod29 root")
    spine = tuple(
        index
        for index in range(spine_period)
        if atom_generated_rank_seven(4 * index * index + 12 * index + 11)
    )
    require(
        all(
            atom_generated_rank_seven(4 * index * index + 12 * index + 11)
            == (
                index % 51 in roots_51
                and index % 9 not in (1, 5)
                and index % 11 not in (0, 8)
            )
            for index in range(spine_period)
        ),
        "U-spine residue law",
    )
    require((len(spine), Fraction(len(spine), spine_period)) == (72, Fraction(8, 187)), len(spine))

    level = ((3, 4, 5),)
    levels = []
    for depth in range(11):
        accepting = sum(atom_generated_rank_seven(2 * vector[2] + 1) for vector in level)
        levels.append((depth, len(level), len(set(level)), accepting))
        level = tuple(
            matrix_child(vector, matrix, modulus)
            for vector in level
            for matrix in BERGGREN_MATRICES
        )
    return spine_period, roots_51, spine, Fraction(8, 187), tuple(levels)


def main():
    dependencies = tuple((label, lf_hash(path)) for label, path, _ in PINNED)
    for label, path, expected in PINNED:
        require(lf_hash(path) == expected, (label, lf_hash(path), expected))

    atoms = atom_replay()
    modulus = lcm(*(LOWER_BASES + ATOM_BASES))
    require(modulus == 14_362_718_970_600, modulus)
    density = union_minus_union_density(ATOM_BASES, LOWER_BASES)
    require(density == Fraction(165_741_596, 1_554_406_815), density)
    cumulative_lower_bound = Fraction(149, 345) + density
    require(cumulative_lower_bound == Fraction(837_065_119, 1_554_406_815), cumulative_lower_bound)

    fibonacci = fibonacci_transport(modulus)
    berggren = berggren_transport(modulus)
    semantic_surface = (atoms, modulus, density, cumulative_lower_bound, fibonacci, berggren)
    semantic_digest = sha256(repr(semantic_surface).encode("ascii")).hexdigest()
    if EXPECTED_SEMANTIC_DIGEST is not None:
        require(semantic_digest == EXPECTED_SEMANTIC_DIGEST, (semantic_digest, EXPECTED_SEMANTIC_DIGEST))

    print("Rank-seven atom-monoid recurrence transport")
    print(f"dependency_sha256_lf={dependencies}")
    print("status=PROVED_EXACT_INCLUSION atom_generated_degrees_have_rho_ZMC_7;FINITE_EXACT_Berggren_depth10_prefix;NO_completeness_claim;no_LRC14_decrement")
    print(f"atoms_replayed={atoms}")
    print(f"atom_monoid=(period,density_and_harmonic_coefficient,cumulative_rank_le7_lower_bound)={(modulus,density,cumulative_lower_bound)}")
    print(f"Fibonacci=(apparition,period,count,density,law)={(fibonacci[0],fibonacci[1],len(fibonacci[2]),fibonacci[3],'7|n and 6,10,15,25 do not divide n')}")
    print(f"Berggren_U_spine=(period,roots_mod51,count,density,law)={(berggren[0],berggren[1],len(berggren[2]),berggren[3],'n mod51 in roots; n mod9 not 1,5; n mod11 not 0,8')}")
    print(f"Berggren_ternary_atom_monoid_levels=(depth,nodes,states,accepting)={berggren[4]}")
    print(f"semantic_sha256={semantic_digest}")
    print("scope=verified_rank7_atoms_times_fibre_degrees_outside_THM3416_rank_le6_support;density_is_unconditional_for_this_subfamily_and_only_a_lower_bound_for_full_rank7;harmonic_coefficients_refer_to_periodic_degree_or_index_subsets;LRC14_open")


if __name__ == "__main__":
    main()
