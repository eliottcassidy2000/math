#!/usr/bin/env python3
"""Finite-exact scout for the next zero-mode-cochain grade after THM-3416.

This file does not claim an all-q rank-seven theorem.  It imports the frozen
augmented-bank and rare-coordinate solver from the audited THM-3416 companion,
restricts to moduli outside the proved rank-at-most-six support, and searches
both twists through Q=200 with cap seven.  The resulting seven-element divisor
antichain is a named hypothesis and a proof target, not an extrapolation.
"""

from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from math import gcd, lcm
from pathlib import Path

from lrc_zero_mode_cochain_rank6_support_probe_20260815 import (
    augmented_bank,
    danger_mask,
    rare_coordinate_cover,
)


ROOT = Path(__file__).resolve().parents[1]
SCAN_LIMIT = 200
LOWER_BASES = (8, 9, 10, 11, 12, 15, 23, 25)
EXPECTED_BASES = (13, 14, 29, 38, 51, 68, 148)
PINNED = (
    (
        "THM-3416",
        ROOT / "01-canon/theorems/THM-3416-zero-mode-cochain-global-rank-six-support.md",
        "42a9309145de51d1bb6fca0b7c1945302ff37a63a3183e1dfed838c07118e8bf",
    ),
    (
        "THM-3416-primary-script",
        ROOT / "04-computation/lrc_zero_mode_cochain_rank6_support_probe_20260815.py",
        "ed7672fd75b7c5ede23c0b7752e06849faa9c693d18e0639eee5b39f17d03a21",
    ),
    (
        "THM-3416-primary-output",
        ROOT / "05-knowledge/results/lrc_zero_mode_cochain_rank6_support_probe_20260815.out",
        "3eb932643e76f2ac7836d8d435bf259183e2c3778216b0458bf4487b5894102a",
    ),
)

EXPECTED_EVENT_DIGEST = "c888a5f6a066d35abffaa60eed5ea1ac20dcdabaaf33237cdfdafab9355e9849"
EXPECTED_SEMANTIC_DIGEST = "dec1ecff2442edb765582647ba529b338a99aa0e04fb7f27313fe97ac92c406a"


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def lf_hash(path):
    return sha256(path.read_bytes().replace(bytes((13, 10)), bytes((10,)))).hexdigest()


def target_free(q):
    return not any(q % base == 0 for base in LOWER_BASES)


def order_profile(q, residues):
    return tuple(sorted(q // gcd(q, residue) for residue in residues))


def scan_record(q, epsilon):
    modulus, primes, raw, unique, maximal, full = augmented_bank(q, epsilon)
    witness, nodes, branches = rare_coordinate_cover(full, maximal, 7)
    if witness:
        require(gcd(modulus, *witness) == 1, (q, epsilon, witness))
    return (
        q,
        epsilon,
        raw,
        len(unique),
        len(maximal),
        witness,
        order_profile(q, witness),
        nodes,
        branches,
        len(primes),
    )


def minimal_divisor_antichain(keys):
    return tuple(q for q in keys if not any(q % smaller == 0 for smaller in keys if smaller < q))


def atom_record(q, epsilon, residues):
    masks = tuple(danger_mask(q, residue, epsilon) for residue in residues)
    full = (1 << q) - 1
    union = 0
    xor = 0
    for mask in masks:
        union |= mask
        xor ^= mask
    require(union == full, (q, epsilon, residues, union, full))
    multiplicities = tuple(sum(mask >> sheet & 1 for mask in masks) for sheet in range(q))
    collisions = tuple(
        (sheet, tuple(index for index, mask in enumerate(masks) if mask >> sheet & 1))
        for sheet in range(q)
        if multiplicities[sheet] > 1
    )
    pair_graph = tuple(
        (left, right, (masks[left] & masks[right]).bit_count())
        for left in range(len(masks))
        for right in range(left + 1, len(masks))
        if masks[left] & masks[right]
    )
    return (
        q,
        epsilon,
        residues,
        order_profile(q, residues),
        tuple(mask.bit_count() for mask in masks),
        tuple(sorted(Counter(multiplicities).items())),
        xor.bit_count(),
        collisions,
        pair_graph,
    )


def conditional_density(new_bases, lower_bases):
    answer = Fraction(0)
    for new_size in range(1, len(new_bases) + 1):
        for new_packet in combinations(new_bases, new_size):
            new_lcm = lcm(*new_packet)
            for lower_size in range(len(lower_bases) + 1):
                for lower_packet in combinations(lower_bases, lower_size):
                    modulus = lcm(new_lcm, *lower_packet) if lower_packet else new_lcm
                    sign = 1 if (new_size + lower_size) % 2 else -1
                    answer += sign * Fraction(1, modulus)
    return answer


def main():
    dependencies = tuple((label, lf_hash(path)) for label, path, _ in PINNED)
    for label, path, expected in PINNED:
        require(lf_hash(path) == expected, (label, lf_hash(path), expected))

    records = tuple(
        scan_record(q, epsilon)
        for q in range(2, SCAN_LIMIT + 1)
        if target_free(q)
        for epsilon in (0, 1)
    )
    positives = tuple((row[0], row[1], row[5], row[6]) for row in records if row[5])
    keys = tuple(sorted({row[0] for row in records if row[5]}))
    expected_keys = (13, 14, 26, 28, 29, 38, 42, 51, 52, 58, 68, 148)
    require(keys == expected_keys, keys)
    bases = minimal_divisor_antichain(keys)
    require(bases == EXPECTED_BASES, bases)

    global_rank_seven = tuple(
        q
        for q in range(2, SCAN_LIMIT + 1)
        if target_free(q) and any(q % base == 0 for base in bases)
    )
    observed_by_divisor = tuple(
        q
        for q in range(2, SCAN_LIMIT + 1)
        if target_free(q) and any(q % key == 0 for key in keys)
    )
    require(global_rank_seven == observed_by_divisor, (global_rank_seven, observed_by_divisor))
    require(
        global_rank_seven
        == (13, 14, 26, 28, 29, 38, 39, 42, 51, 52, 58, 65, 68, 76, 78,
            87, 91, 98, 102, 114, 116, 145, 148, 169, 174, 182, 196),
        global_rank_seven,
    )

    witnesses = {(q, epsilon): residues for q, epsilon, residues, _ in positives}
    atoms = (
        atom_record(13, 1, witnesses[(13, 1)]),
        atom_record(14, 1, witnesses[(14, 1)]),
        atom_record(29, 0, witnesses[(29, 0)]),
        atom_record(29, 1, witnesses[(29, 1)]),
        atom_record(38, 1, witnesses[(38, 1)]),
        atom_record(51, 1, witnesses[(51, 1)]),
        atom_record(68, 1, witnesses[(68, 1)]),
        atom_record(148, 1, witnesses[(148, 1)]),
    )
    require(atoms[0][5:] == (((1, 13),), 13, (), ()), atoms[0])
    require(atoms[1][5:] == (((1, 14),), 14, (), ()), atoms[1])
    require(atoms[2][5] == ((1, 28), (7, 1)) and atoms[2][6] == 29, atoms[2])
    require(atoms[3][5] == ((1, 28), (3, 1)) and atoms[3][6] == 29, atoms[3])
    require(atoms[4][5] == ((1, 34), (2, 4)) and atoms[4][6] == 34, atoms[4])
    require(atoms[6][5] == ((1, 64), (3, 4)) and atoms[6][6] == 68, atoms[6])
    require(atoms[7][5] == ((1, 140), (2, 8)) and atoms[7][6] == 140, atoms[7])

    hostiles = tuple(
        q
        for q in (17, 19, 31, 37, 43, 67, 74, 95, 127, 149, 199)
        if target_free(q) and q not in keys
    )
    require(hostiles == (17, 19, 31, 37, 43, 67, 74, 95, 127, 149, 199), hostiles)

    hypothesis_density = conditional_density(bases, LOWER_BASES)
    require(hypothesis_density == Fraction(165741596, 1554406815), hypothesis_density)
    cumulative_density = Fraction(149, 345) + hypothesis_density
    require(cumulative_density == Fraction(837065119, 1554406815), cumulative_density)

    event_surface = tuple((row[0], row[1], row[5], row[6]) for row in records)
    event_digest = sha256(repr(event_surface).encode("ascii")).hexdigest()
    semantic_surface = (
        positives,
        keys,
        bases,
        global_rank_seven,
        atoms,
        hostiles,
        hypothesis_density,
        cumulative_density,
    )
    semantic_digest = sha256(repr(semantic_surface).encode("ascii")).hexdigest()
    if EXPECTED_EVENT_DIGEST is not None:
        require(event_digest == EXPECTED_EVENT_DIGEST, (event_digest, EXPECTED_EVENT_DIGEST))
    if EXPECTED_SEMANTIC_DIGEST is not None:
        require(semantic_digest == EXPECTED_SEMANTIC_DIGEST, (semantic_digest, EXPECTED_SEMANTIC_DIGEST))

    print("Zero-mode-cochain rank-seven divisor-ancestor scout")
    print(f"dependency_sha256_lf={dependencies}")
    print("status=FINITE_EXACT target_free_Q2_Q200_both_twists_cap7;HYPOTHESIS_ONLY all_q_rank7_bases_13,14,29,38,51,68,148;no_LRC14_decrement")
    print(f"scan_counts=(records,raw,unique,maximal,nodes,branches)={(len(records),sum(r[2] for r in records),sum(r[3] for r in records),sum(r[4] for r in records),sum(r[7] for r in records),sum(r[8] for r in records))}")
    print(f"primitive_exact7_target_free_Q2_Q200={positives}")
    print(f"candidate_minimal_divisor_antichain={bases}")
    print(f"global_exact7_Q2_Q200={global_rank_seven}")
    print(f"rank7_atom_collision_atlas=(q,epsilon,residues,orders,sizes,multiplicity_hist,xor_size,collision_hyperedges,pair_graph)={atoms}")
    print(f"hostile_no_primitive_cap7_controls={hostiles}")
    print(f"conditional_if_all_q_hypothesis=(rank7_harmonic_coefficient,cumulative_through7)={(hypothesis_density,cumulative_density)}")
    print(f"event_sha256={event_digest}")
    print(f"semantic_sha256={semantic_digest}")
    print("scope=primitive_augmented_search_only_outside_proved_rank_le6_support;Q200_finite_boundary;candidate_density_is_conditional;collision_hypergraph_not_tournament;LRC14_open")


if __name__ == "__main__":
    main()
