#!/usr/bin/env python3
"""Finite exact scout for the prime half-twist classes left open by THM-3420.

THM-3420 closes fixed-zero primes, rules out half-twist prime classes 3 and 5
modulo 14, and closes class 13.  This companion searches only the remaining
classes 1, 9, and 11 through p=500.  It is a hostile finite census, not an
all-prime theorem.
"""

from collections import Counter
from hashlib import sha256
from math import isqrt
from pathlib import Path

from lrc_zero_mode_cochain_rank6_support_probe_20260815 import (
    augmented_bank,
    danger_mask,
    rare_coordinate_cover,
)


ROOT = Path(__file__).resolve().parents[1]
SCAN_LIMIT = 500
LIVE_CLASSES = (1, 9, 11)
EXPECTED_EVENT_DIGEST = "6111e8a3484a1775999f58ce355d4cecaf34d7f7af40f7d76dbda9a9ec964f49"
EXPECTED_SEMANTIC_DIGEST = "ee528fc41630c7a9146fe62c5976fbf5421c3c4f427de07d2d7a6e8b9e777959"
PINNED = (
    (
        "THM-3420",
        ROOT / "01-canon/theorems/THM-3420-prime-rank-seven-zero-and-half-twist-splitter-closures.md",
        "5f307c906c054850f2589b2e9fb7c0e902b3571747a6c61ebc2c2f2628212b2f",
    ),
    (
        "THM-3420-script",
        ROOT / "04-computation/lrc_prime_rank7_splitter_closures_thm3420.py",
        "f95aaf081ccd5c92cb7474104f242623bba9246bdc19eb04ce8da81f3c4e2af6",
    ),
    (
        "THM-3420-output",
        ROOT / "05-knowledge/results/lrc_prime_rank7_splitter_closures_thm3420.out",
        "22740907984a9e64e75208c0a01ce0222fd844dfe65826685ff11e64274ef959",
    ),
    (
        "THM-3416-primary-solver",
        ROOT / "04-computation/lrc_zero_mode_cochain_rank6_support_probe_20260815.py",
        "ed7672fd75b7c5ede23c0b7752e06849faa9c693d18e0639eee5b39f17d03a21",
    ),
)


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def lf_hash(path):
    return sha256(path.read_bytes().replace(bytes((13, 10)), bytes((10,)))).hexdigest()


def is_prime(value):
    if value < 2:
        return False
    return all(value % divisor for divisor in range(2, isqrt(value) + 1))


def atom_profile(p, witness):
    masks = tuple(danger_mask(p, residue, 1) for residue in witness)
    multiplicities = tuple(
        sum(mask >> sheet & 1 for mask in masks)
        for sheet in range(p)
    )
    full = (1 << p) - 1
    joined = 0
    xor = 0
    for mask in masks:
        joined |= mask
        xor ^= mask
    require(joined == full, (p, witness, joined, full))
    even_blocks = sum(residue % 2 == 0 for residue in witness)
    fixed_sheet = (p - 1) // 2
    require(multiplicities[fixed_sheet] == even_blocks, (p, witness, multiplicities[fixed_sheet]))
    return (
        even_blocks,
        tuple(mask.bit_count() for mask in masks),
        tuple(sorted(Counter(multiplicities).items())),
        p - xor.bit_count(),
    )


def main():
    dependencies = tuple((label, lf_hash(path)) for label, path, _ in PINNED)
    for label, path, expected in PINNED:
        require(lf_hash(path) == expected, (label, lf_hash(path), expected))

    primes = tuple(
        p for p in range(29, SCAN_LIMIT + 1)
        if is_prime(p) and p % 14 in LIVE_CLASSES
    )
    records = []
    for p in primes:
        modulus, factors, raw, unique, maximal, full = augmented_bank(p, 1)
        witness, nodes, branches = rare_coordinate_cover(full, maximal, 7)
        require(modulus == 2 * p and factors == (2, p), (p, modulus, factors))
        profile = atom_profile(p, witness) if witness else ()
        records.append(
            (p, p % 14, raw, len(unique), len(maximal), witness, profile, nodes, branches)
        )
    records = tuple(records)

    positives = tuple((row[0], row[5], row[6]) for row in records if row[5])
    require(
        positives
        == ((
            29,
            (1, 5, 7, 8, 12, 13, 22),
            (3, (4, 4, 4, 5, 5, 4, 5), ((1, 28), (3, 1)), 0),
        ),),
        positives,
    )
    require(tuple(row[0] for row in records[:4]) == (29, 37, 43, 53), records[:4])
    require(tuple(row[1] for row in records[:4]) == (1, 9, 1, 11), records[:4])

    hostile_controls = tuple(
        (row[0], row[1], row[7], row[8])
        for row in records
        if row[0] in (37, 43, 53, 211, 499)
    )
    require(
        tuple(row[:2] for row in hostile_controls)
        == ((37, 9), (43, 1), (53, 11), (211, 1), (499, 9)),
        hostile_controls,
    )

    event_surface = tuple((row[0], row[1], row[5]) for row in records)
    semantic_surface = (primes, positives, hostile_controls)
    event_digest = sha256(repr(event_surface).encode("ascii")).hexdigest()
    semantic_digest = sha256(repr(semantic_surface).encode("ascii")).hexdigest()
    if EXPECTED_EVENT_DIGEST is not None:
        require(event_digest == EXPECTED_EVENT_DIGEST, (event_digest, EXPECTED_EVENT_DIGEST))
    if EXPECTED_SEMANTIC_DIGEST is not None:
        require(semantic_digest == EXPECTED_SEMANTIC_DIGEST, (semantic_digest, EXPECTED_SEMANTIC_DIGEST))

    print("Prime rank-seven remaining half-twist classes scout")
    print(f"dependency_sha256_lf={dependencies}")
    print("status=FINITE_EXACT primes_p29_p500_classes_1,9,11_mod14_half_twist_cap7;NO_all_prime_theorem;no_LRC14_decrement")
    print(f"scan_counts=(primes,raw,unique,maximal,nodes,branches)={(len(records),sum(row[2] for row in records),sum(row[3] for row in records),sum(row[4] for row in records),sum(row[7] for row in records),sum(row[8] for row in records))}")
    print(f"positives=(p,residues,(even_blocks,sizes,multiplicity_hist,xor_defect))={positives}")
    print(f"negative_primes={tuple(row[0] for row in records if not row[5])}")
    print(f"hostile_controls=(p,class,nodes,branches)={hostile_controls}")
    print(f"event_sha256={event_digest}")
    print(f"semantic_sha256={semantic_digest}")
    print("scope=THM3420_live_prime_half_classes_only;finite_boundary_500;generic_exact_augmented_solver;fixed_zero_and_other_classes_not_rescanned;LRC14_open")


if __name__ == "__main__":
    main()
