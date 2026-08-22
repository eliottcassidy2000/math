#!/usr/bin/env python3
"""Exact collision/breaker state and recurrence transport for rank-seven atoms.

The seven-base antichain is not known to be complete.  This companion only
refines the already proved atom-generated subfamily.  It retains multiplicity
polynomials, collision genus, quotient-order joint period, the half-twist
parity breaker, and the multiplicative fibre defect that support-only
automata forget.
"""

from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from math import gcd, lcm
from pathlib import Path

from lrc_zero_mode_cochain_rank6_support_probe_20260815 import (
    BERGGREN_MATRICES,
    danger_mask,
)


ROOT = Path(__file__).resolve().parents[1]
LOWER_BASES = (8, 9, 10, 11, 12, 15, 23, 25)
ATOM_BASES = (13, 14, 29, 38, 51, 68, 148)
PARITY_PERFECT_BASES = (13, 14, 29, 68)
EVEN_DEFECT_BASES = (38, 51, 148)
EXPECTED_SEMANTIC_DIGEST = "2dc0c1826a03b19b797ddd702c9c877cd77d0a8ec9e0883a5ca13bee22936d0d"

PINNED = (
    (
        "THM-3416",
        ROOT / "01-canon/theorems/THM-3416-zero-mode-cochain-global-rank-six-support.md",
        "42a9309145de51d1bb6fca0b7c1945302ff37a63a3183e1dfed838c07118e8bf",
    ),
    (
        "THM-3421",
        ROOT / "01-canon/theorems/THM-3421-prime-half-twist-rank-seven-classification.md",
        "2f577354a06628660d90f70aca34378a186b7e67a55795e6c4a5c32a255d9736",
    ),
    (
        "rank-seven-scout",
        ROOT / "04-computation/lrc_zero_mode_cochain_rank7_ancestor_scout_20260815.py",
        "bcf80194bd5697b121b08b239f096d67d7e0214bbc81abf4a9fb1f6d939fa560",
    ),
    (
        "rank-seven-scout-output",
        ROOT / "05-knowledge/results/lrc_zero_mode_cochain_rank7_ancestor_scout_20260815.out",
        "95012438ffa2e36183ef3e56a9881c3daa52ddfe118abebfd02106b1b28f65a2",
    ),
    (
        "rank-seven-transport",
        ROOT / "04-computation/lrc_rank7_atom_monoid_recurrence_transport_20260815.py",
        "5248721148ca18bd82cae9bda3f8fa86b3bf2fc35fbce88369c5d3c1c0d3c36a",
    ),
    (
        "rank-seven-transport-output",
        ROOT / "05-knowledge/results/lrc_rank7_atom_monoid_recurrence_transport_20260815.out",
        "572c001aaa647774d19b53a1309ad10f7c8044426db15d145fbcfbd7a81e5cc6",
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

EXPECTED_ATOM_STATES = (
    (13, 1, ((13, 7),), ((1, 13),), 0, 0, 0, 13, 1, 1),
    (14, 1, ((7, 1), (14, 6)), ((1, 14),), 0, 0, 0, 14, 1, 1),
    (29, 0, ((29, 7),), ((1, 28), (7, 1)), 6, 0, 3, 29, 1, 1),
    (29, 1, ((29, 7),), ((1, 28), (3, 1)), 2, 0, 1, 29, 1, 1),
    (38, 1, ((19, 1), (38, 6)), ((1, 34), (2, 4)), 4, 4, 0, 38, 1, 1),
    (51, 1, ((3, 1), (17, 2), (51, 4)), ((1, 42), (2, 4), (3, 3), (4, 2)), 16, 6, 5, 51, 1, 1),
    (68, 1, ((17, 3), (68, 4)), ((1, 64), (3, 4)), 8, 0, 4, 68, 1, 1),
    (148, 1, ((37, 3), (148, 4)), ((1, 140), (2, 8)), 8, 8, 0, 148, 1, 1),
)


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def lf_hash(path):
    return sha256(path.read_bytes().replace(bytes((13, 10)), bytes((10,)))).hexdigest()


def order_histogram(q, residues):
    orders = (q // gcd(q, residue) for residue in residues)
    return tuple(sorted(Counter(orders).items()))


def certificate(q, epsilon, residues):
    masks = tuple(danger_mask(q, residue, epsilon) for residue in residues)
    full = (1 << q) - 1
    joined = 0
    for mask in masks:
        joined |= mask
    require(joined == full, (q, epsilon, residues, "not a cover"))

    multiplicities = tuple(
        sum((mask >> sheet) & 1 for mask in masks)
        for sheet in range(q)
    )
    polynomial = tuple(sorted(Counter(multiplicities).items()))
    omega = sum(multiplicities) - q
    even_defect = sum(value % 2 == 0 for value in multiplicities)
    collision_genus = sum((value - 1) // 2 for value in multiplicities)
    require(omega == even_defect + 2 * collision_genus,
            (q, omega, even_defect, collision_genus))

    orders = tuple(q // gcd(q, residue) for residue in residues)
    joint_period = lcm(*orders)
    period_defect = q // joint_period
    odd_residue = int(any(residue % 2 for residue in residues))
    augmented_modulus = q if epsilon == 0 else 2 * q
    augmented_gcd = gcd(augmented_modulus, *residues)
    return (
        q,
        epsilon,
        order_histogram(q, residues),
        polynomial,
        omega,
        even_defect,
        collision_genus,
        joint_period,
        period_defect,
        odd_residue,
        augmented_gcd,
    )


def atom_state_audit():
    rows = tuple(certificate(*atom) for atom in ATOMS)
    surface = tuple(row[:10] for row in rows)
    require(surface == EXPECTED_ATOM_STATES, surface)
    require(all(row[10] == 1 for row in rows), rows)
    require(tuple(row[5] == 0 for row in rows) ==
            (True, True, True, True, False, False, True, False), rows)
    return rows


def fibre_audit(atom_rows):
    records = []
    for atom, base in zip(ATOMS, atom_rows):
        q, epsilon, residues = atom
        for fibre in (2, 3, 5):
            lifted = certificate(fibre * q, epsilon, tuple(fibre * r for r in residues))
            expected_polynomial = tuple((power, fibre * count) for power, count in base[3])
            require(lifted[2] == base[2], (atom, fibre, lifted[2], base[2]))
            require(lifted[3] == expected_polynomial, (atom, fibre, lifted[3]))
            require(lifted[4:7] == tuple(fibre * value for value in base[4:7]),
                    (atom, fibre, lifted[4:7], base[4:7]))
            require(lifted[7] == q and lifted[8] == fibre, (atom, fibre, lifted[7:9]))
            expected_eta = int(bool(fibre % 2 and base[9]))
            require(lifted[9] == expected_eta, (atom, fibre, lifted[9], expected_eta))
            records.append((q, epsilon, fibre, *lifted[3:10]))
    return tuple(records)


def promotion_hostiles():
    scaled_26 = certificate(26, 1, tuple(2 * r for r in ATOMS[0][2]))
    primitive_26 = certificate(26, 1, (1, 4, 6, 7, 10, 19, 25))
    require(scaled_26[3:10] == (((1, 26),), 0, 0, 0, 13, 2, 0), scaled_26)
    require(primitive_26[3:10] == (((1, 26),), 0, 0, 0, 26, 1, 1), primitive_26)

    scaled_58 = certificate(58, 1, tuple(2 * r for r in ATOMS[3][2]))
    primitive_58 = certificate(58, 1, (4, 21, 25, 33, 37, 48, 54))
    require(scaled_58[3:10] == (((1, 56), (3, 2)), 4, 0, 2, 29, 2, 0), scaled_58)
    require(primitive_58[3:10] == (((1, 56), (2, 2)), 2, 2, 0, 58, 1, 1), primitive_58)
    return scaled_26[2:11], primitive_26[2:11], scaled_58[2:11], primitive_58[2:11]


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


def density_split():
    parity_perfect = union_minus_union_density(PARITY_PERFECT_BASES, LOWER_BASES)
    even_only = union_minus_union_density(
        EVEN_DEFECT_BASES,
        LOWER_BASES + PARITY_PERFECT_BASES,
    )
    total = union_minus_union_density(ATOM_BASES, LOWER_BASES)
    require(parity_perfect == Fraction(67404, 737035), parity_perfect)
    require(even_only == Fraction(4717312, 310881363), even_only)
    require(total == Fraction(165741596, 1554406815), total)
    require(parity_perfect + even_only == total, (parity_perfect, even_only, total))
    return parity_perfect, even_only, total


def fibonacci_refined_states():
    apparition = {
        8: 6, 9: 12, 10: 15, 11: 10, 12: 12, 15: 20, 23: 24, 25: 25,
        13: 7, 14: 24, 29: 14, 38: 18, 51: 36, 68: 18, 148: 114,
    }
    full_period = lcm(*apparition.values())
    require(full_period == 239400, full_period)
    labels = Counter()
    for index in range(full_period):
        accepted = (
            index % 7 == 0
            and index % 6 != 0
            and index % 10 != 0
            and index % 15 != 0
            and index % 25 != 0
        )
        if not accepted:
            continue
        active = tuple(base for base in (13, 29) if index % apparition[base] == 0)
        fibre_parity = "even" if index % 3 == 0 else "odd"
        labels[(active, fibre_parity)] += 1

    expected = (
        ((13,), "odd", 10944, Fraction(8, 175)),
        ((13,), "even", 4560, Fraction(2, 105)),
        ((13, 29), "odd", 9120, Fraction(4, 105)),
    )
    key_order = (((13,), "odd"), ((13,), "even"), ((13, 29), "odd"))
    observed = tuple(
        (active, parity, labels[(active, parity)],
         Fraction(labels[(active, parity)], full_period))
        for active, parity in key_order
    )
    require(observed == expected, observed)
    require(sum(row[3] for row in observed) == Fraction(18, 175), observed)

    f14 = 377
    lifted_13 = certificate(f14, 1, tuple(29 * r for r in ATOMS[0][2]))
    lifted_29_zero = certificate(f14, 0, tuple(13 * r for r in ATOMS[2][2]))
    lifted_29_half = certificate(f14, 1, tuple(13 * r for r in ATOMS[3][2]))
    certificates = (
        (13, f14 // 13, *lifted_13[3:9]),
        (29, f14 // 29, *lifted_29_zero[3:9]),
        (29, f14 // 29, *lifted_29_half[3:9]),
    )
    require(certificates[0][1] == 29 and certificates[1][1] == 13, certificates)
    return full_period, observed, certificates


def determinant_three(matrix):
    return (
        matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
    )


def berggren_refined_states():
    odd_modulus = lcm(9, 11, 15, 23, 25, 13, 29, 51)
    require(odd_modulus == 364832325, odd_modulus)
    determinants = tuple(determinant_three(matrix) for matrix in BERGGREN_MATRICES)
    require(all(abs(value) == 1 for value in determinants), determinants)

    spine_period = 1683
    rows = []
    for index in range(spine_period):
        q = 4 * index * index + 12 * index + 11
        accepted = (
            index % 51 in (2, 19, 29, 46)
            and index % 9 not in (1, 5)
            and index % 11 not in (0, 8)
        )
        if not accepted:
            continue
        require(q % 51 == 0 and q % 2 == 1, (index, q))
        require(q % 13 and q % 29, (index, q))
        fibre = q // 51
        require(fibre % 2 == 1, (index, q, fibre))
        rows.append((index, q, fibre, 16 * fibre, 6 * fibre, 5 * fibre, 51, fibre, 1))
    require(len(rows) == 72, len(rows))
    require(rows[0][:3] == (2, 51, 1), rows[0])
    row_29 = next(row for row in rows if row[0] == 29)
    require(row_29 == (29, 3723, 73, 1168, 438, 365, 51, 73, 1), row_29)

    # Each determinant is a unit modulo odd_modulus.  Hence the three child
    # maps are permutations of the finite root orbit.  The Cesaro statement
    # recorded by the reflection follows from the finite irreducible Markov
    # chain obtained by averaging these permutations; no orbit enumeration is
    # claimed here.
    return odd_modulus, determinants, len(rows), rows[0], row_29


def main():
    dependencies = tuple((label, lf_hash(path)) for label, path, _ in PINNED)
    for label, path, expected in PINNED:
        require(lf_hash(path) == expected, (label, lf_hash(path), expected))

    atoms = atom_state_audit()
    fibres = fibre_audit(atoms)
    hostiles = promotion_hostiles()
    densities = density_split()
    fibonacci = fibonacci_refined_states()
    berggren = berggren_refined_states()
    semantic_surface = (atoms, fibres, hostiles, densities, fibonacci, berggren)
    semantic_digest = sha256(repr(semantic_surface).encode("ascii")).hexdigest()
    if EXPECTED_SEMANTIC_DIGEST is not None:
        require(semantic_digest == EXPECTED_SEMANTIC_DIGEST,
                (semantic_digest, EXPECTED_SEMANTIC_DIGEST))

    print("Rank-seven collision/breaker recurrence state")
    print(f"dependency_sha256_lf={dependencies}")
    print("status=PROVED_EXACT_REFINEMENT_OF_ATOM_GENERATED_SUBFAMILY;FINITE_STATE_PLUS_INTEGER_FIBRE_COCYCLE;NO_ANTICHAIN_COMPLETENESS;no_LRC14_decrement")
    print(f"atom_states=(q,epsilon,order_hist,P_hist,Omega,even_defect,G,L,Delta,eta,gcd_aug)={atoms}")
    print(f"fibre_law_controls=(q,epsilon,k,P_hist,Omega,even_defect,G,L,Delta,eta)={fibres}")
    print(f"period_promotion_hostiles=(scaled26,primitive26,scaled58,primitive58)={hostiles}")
    print(f"atom_density_split=(parity_perfect_available,even_defect_only,total)={densities}")
    print(f"Fibonacci_refined=(label_period,states,F14_certificates)={fibonacci}")
    print(f"Berggren_refined=(odd_modulus,determinants,accepted_U_spine_count,first,n29)={berggren}")
    print(f"semantic_sha256={semantic_digest}")
    print("scope=certificate_valued_refinement_only;collision_anatomy_does_not_determine_joint_period;finite_control_needs_unbounded_multiplicative_fibre_cocycle;ternary_Cesaro_uniformity_is_abstract_not_an_orbit_census;composite_rank7_and_LRC14_open")


if __name__ == "__main__":
    main()
