#!/usr/bin/env python3
"""Exact controls for the all-q atomized coset-cochain theorem.

The analytic theorem is unbounded.  This verifier declares three finite
universes: every multiplier residue and phase cell for q<=10; every distinct
atom tuple of rank <=3 with owner speeds <=10, and rank 4 with owner speeds
<=5; and every owner subset of {1,2,3,4}, with all capacity-lawful atom
selections, for q<=10.  Repeated-owner atoms are retained literally.
"""

from __future__ import annotations

import ast
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
from math import gcd
from pathlib import Path


EXPECTED_SEMANTIC = "39aade82ac505a78e51d779e500acfc1d1b916c5350632311e9274b727e55aec"


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def floor_q(value):
    return value.numerator // value.denominator


def ceil_q(value):
    return -floor_q(-value)


def mod_one(value):
    return value - floor_q(value)


def distance_to_integer(value):
    residue = mod_one(value)
    return min(residue, 1 - residue)


def strict_danger(q, speed, label, time):
    return 14 * distance_to_integer(speed * (time + F(label, q))) < 1


def image_order(q, speed):
    return q // gcd(q, speed)


def capacity_atoms(q, speed):
    m = image_order(q, speed)
    return (m + 6) // 7


def coset_mask(q, speed, label):
    m = image_order(q, speed)
    return sum(1 << sheet for sheet in range(q) if sheet % m == label % m)


def owner_mask(q, speed, time):
    return sum(
        1 << sheet
        for sheet in range(q)
        if strict_danger(q, speed, sheet, time)
    )


def circle_event_samples(q, speeds, include_endpoints=False):
    events = {F(0)}
    for speed in speeds:
        for label in range(q):
            for tooth in range(speed):
                centre = F(tooth, speed) - F(label, q)
                radius = F(1, 14 * speed)
                events.add(mod_one(centre - radius))
                events.add(mod_one(centre + radius))
    ordered = tuple(sorted(events))
    midpoints = []
    for index, left in enumerate(ordered):
        right = ordered[(index + 1) % len(ordered)]
        if index + 1 == len(ordered):
            right += 1
        midpoints.append(mod_one((left + right) / 2))
    if include_endpoints:
        return ordered + tuple(midpoints)
    return tuple(midpoints)


def phase_cell_samples(m):
    events = {F(0)}
    for residue in range(m):
        centre = -F(residue, m)
        events.add(mod_one(centre - F(1, 14)))
        events.add(mod_one(centre + F(1, 14)))
    ordered = tuple(sorted(events))
    midpoints = []
    for index, left in enumerate(ordered):
        right = ordered[(index + 1) % len(ordered)]
        if index + 1 == len(ordered):
            right += 1
        midpoints.append(mod_one((left + right) / 2))
    return ordered + tuple(midpoints)


def exact_local_multiplicity_audit():
    checks = 0
    profile = Counter()
    for q in range(2, 11):
        for residue_speed in range(q):
            g = gcd(q, residue_speed)
            m = q // g
            counts = set()
            for theta in phase_cell_samples(m):
                image_hits = tuple(
                    residue
                    for residue in range(m)
                    if 14 * distance_to_integer(theta + F(residue, m)) < 1
                )
                sheet_hits = tuple(
                    sheet
                    for sheet in range(q)
                    if 14 * distance_to_integer(theta + F(residue_speed * sheet, q)) < 1
                )
                expected_sheets = tuple(
                    sheet
                    for sheet in range(q)
                    if (residue_speed // g * sheet) % m in image_hits
                )
                require(sheet_hits == expected_sheets, ("pullback cosets", q, residue_speed, theta))
                require(len(sheet_hits) == g * len(image_hits), ("atom weight", q, residue_speed, theta))
                counts.add(len(image_hits))
                checks += 1
            cap = (m + 6) // 7
            require(counts == {cap - 1, cap}, ("exact multiplicities", q, residue_speed, counts, cap))
            profile[(m, tuple(sorted(counts)))] += 1
    return checks, tuple(sorted(profile.items()))


def ceil_div(left, right):
    return -((-left) // right)


def gap_values(q, left_speed, right_speed, left_label, right_label):
    modulus = q * gcd(left_speed, right_speed)
    residue = ((right_label - left_label) * left_speed * right_speed) % modulus
    bound = (q * (left_speed + right_speed) - 1) // 14
    first = ceil_div(-bound - residue, modulus)
    last = (bound - residue) // modulus
    return tuple(residue + modulus * index for index in range(first, last + 1))


def cochain_witness(q, atoms):
    """Return one skew complete cochain, allowing repeated owner speeds."""
    if not atoms:
        return None
    if len(atoms) == 1:
        return ()
    speeds = tuple(atom[0] for atom in atoms)
    labels = tuple(atom[1] for atom in atoms)
    star_sets = tuple(
        gap_values(q, speeds[0], speeds[index], labels[0], labels[index])
        for index in range(1, len(atoms))
    )
    if any(not values for values in star_sets):
        return None
    for star in product(*star_sets):
        edges = {(0, index): star[index - 1] for index in range(1, len(atoms))}
        good = True
        for left, right in combinations(range(1, len(atoms)), 2):
            numerator = speeds[left] * star[right - 1] - speeds[right] * star[left - 1]
            if numerator % speeds[0]:
                good = False
                break
            value = numerator // speeds[0]
            if value not in gap_values(q, speeds[left], speeds[right], labels[left], labels[right]):
                good = False
                break
            edges[(left, right)] = value
        if good:
            return tuple((left, right, edges[(left, right)]) for left, right in combinations(range(len(atoms)), 2))
    return None


def integers_strictly_between(left, right):
    return range(floor_q(left) + 1, ceil_q(right))


def direct_lifted_interval_witness(q, atoms):
    """Independent real-line tooth enumeration, with no gap variables."""
    if not atoms:
        return None
    first_speed, first_label = atoms[0]
    first_radius = F(1, 14 * first_speed)
    for first_tooth in range(first_speed):
        first_centre = F(first_tooth, first_speed) - F(first_label, q)
        candidate_teeth = []
        for speed, label in atoms[1:]:
            radius = F(1, 14 * speed)
            reach = first_radius + radius
            low = speed * (first_centre + F(label, q) - reach)
            high = speed * (first_centre + F(label, q) + reach)
            choices = tuple(integers_strictly_between(low, high))
            if not choices:
                break
            candidate_teeth.append(choices)
        else:
            for teeth in product(*candidate_teeth):
                centres = (first_centre,) + tuple(
                    F(tooth, speed) - F(label, q)
                    for tooth, (speed, label) in zip(teeth, atoms[1:])
                )
                left = max(
                    centre - F(1, 14 * atom[0])
                    for centre, atom in zip(centres, atoms)
                )
                right = min(
                    centre + F(1, 14 * atom[0])
                    for centre, atom in zip(centres, atoms)
                )
                if left < right:
                    return centres, (left, right)
    return None


def atom_universe(q, speed_bound):
    return tuple(
        (speed, label)
        for speed in range(1, speed_bound + 1)
        for label in range(image_order(q, speed))
    )


def exhaustive_atom_tuple_audit():
    checks = 0
    positives = 0
    repeated_owner = 0
    repeated_positive = 0
    rank_profile = Counter()
    for q in range(2, 11):
        atoms = atom_universe(q, 10)
        for rank in (1, 2, 3):
            for chosen in combinations(atoms, rank):
                direct = direct_lifted_interval_witness(q, chosen) is not None
                cochain = cochain_witness(q, chosen) is not None
                require(direct == cochain, ("atom tuple", q, chosen, direct, cochain))
                repeated = len({speed for speed, _ in chosen}) < rank
                checks += 1
                positives += direct
                repeated_owner += repeated
                repeated_positive += repeated and direct
                rank_profile[(q, rank, direct, repeated)] += 1

        atoms4 = atom_universe(q, 5)
        for chosen in combinations(atoms4, 4):
            direct = direct_lifted_interval_witness(q, chosen) is not None
            cochain = cochain_witness(q, chosen) is not None
            require(direct == cochain, ("rank four", q, chosen, direct, cochain))
            repeated = len({speed for speed, _ in chosen}) < 4
            checks += 1
            positives += direct
            repeated_owner += repeated
            repeated_positive += repeated and direct
            rank_profile[(q, 4, direct, repeated)] += 1
    compact_profile = tuple(
        (
            q,
            rank,
            sum(count for (qq, rr, _, _), count in rank_profile.items() if (qq, rr) == (q, rank)),
            sum(count for (qq, rr, positive, _), count in rank_profile.items() if (qq, rr) == (q, rank) and positive),
            sum(count for (qq, rr, _, repeated), count in rank_profile.items() if (qq, rr) == (q, rank) and repeated),
            sum(count for (qq, rr, positive, repeated), count in rank_profile.items() if (qq, rr) == (q, rank) and positive and repeated),
        )
        for q in range(2, 11)
        for rank in range(1, 5)
    )
    return checks, positives, repeated_owner, repeated_positive, compact_profile


def owner_atom_choices(q, speed):
    labels = tuple(range(image_order(q, speed)))
    cap = capacity_atoms(q, speed)
    choices = [()]
    for size in range(1, cap + 1):
        for chosen in combinations(labels, size):
            atoms = tuple((speed, label) for label in chosen)
            if direct_lifted_interval_witness(q, atoms) is not None:
                choices.append(chosen)
    return tuple(choices)


def event_full_cover(q, owners):
    full = (1 << q) - 1
    for time in circle_event_samples(q, owners):
        mask = 0
        for speed in owners:
            mask |= owner_mask(q, speed, time)
        if mask == full:
            return True
    return False


def atomized_full_cover(q, owners):
    if sum(gcd(q, speed) * capacity_atoms(q, speed) for speed in owners) < q:
        return False
    choices = tuple(owner_atom_choices(q, speed) for speed in owners)
    full = (1 << q) - 1
    for grouped_labels in product(*choices):
        atoms = tuple(
            (speed, label)
            for speed, labels in zip(owners, grouped_labels)
            for label in labels
        )
        if not atoms:
            continue
        mask = 0
        for speed, label in atoms:
            mask |= coset_mask(q, speed, label)
        if mask != full:
            continue
        if cochain_witness(q, atoms) is not None:
            return True
    return False


def exhaustive_owner_cover_audit():
    checks = 0
    positives = 0
    profile = Counter()
    for q in range(2, 11):
        for size in range(1, 5):
            for owners in combinations((1, 2, 3, 4), size):
                physical = event_full_cover(q, owners)
                atomic = atomized_full_cover(q, owners)
                require(physical == atomic, ("owner cover", q, owners, physical, atomic))
                checks += 1
                positives += physical
                profile[(q, size, physical)] += 1
    return checks, positives, tuple(sorted(profile.items()))


def strict_boundary_controls():
    q7_atoms = ((1, 0), (1, 1))
    q8_atoms = ((1, 0), (1, 1))
    q14_atoms = ((1, 0), (1, 1), (1, 2))
    require(gap_values(7, 1, 1, 0, 1) == (), "q7 strict touch")
    require(direct_lifted_interval_witness(7, q7_atoms) is None, "q7 direct touch")
    require(gap_values(8, 1, 1, 0, 1) == (1,), "q8 adjacent gap")
    require(direct_lifted_interval_witness(8, q8_atoms) is not None, "q8 adjacent direct")
    require(capacity_atoms(7, 1) == 1 and capacity_atoms(8, 1) == 2, "q7 q8 capacity")
    require(capacity_atoms(14, 1) == 2, "q14 capacity")
    require(direct_lifted_interval_witness(14, q14_atoms) is None, "q14 third atom")
    require(direct_lifted_interval_witness(14, q14_atoms[:2]) is not None, "q14 two atoms")
    return (
        (7, q7_atoms, False, "14|p|=q(u+v) strict equality"),
        (8, q8_atoms, True, "14|p|<q(u+v)"),
        (14, q14_atoms, False, "three atoms span exactly 1/7"),
    )


def multiatom_full_cover_controls():
    controls = (
        (
            8,
            F(55, 7056),
            (4, 7, 9, 11, 13),
            ((4, 0), (7, 0), (7, 1), (9, 0), (9, 7), (11, 5), (13, 3)),
        ),
        (
            9,
            F(41, 7056),
            (3, 7, 8, 10, 11, 13, 14),
            ((3, 0), (7, 0), (7, 5), (8, 0), (8, 1), (10, 0), (10, 8), (11, 0), (11, 4), (13, 2), (14, 7)),
        ),
        (
            10,
            F(9, 1960),
            (2, 5, 7, 9, 11, 13),
            ((2, 0), (5, 0), (7, 0), (7, 7), (9, 0), (9, 1), (11, 0), (11, 9), (13, 0), (13, 3)),
        ),
    )
    for q, time, owners, expected_atoms in controls:
        full = 0
        actual_atoms = []
        for speed in owners:
            mask = owner_mask(q, speed, time)
            full |= mask
            for label in range(image_order(q, speed)):
                atom_mask = coset_mask(q, speed, label)
                if atom_mask & mask:
                    require(atom_mask & mask == atom_mask, ("partial coset", q, speed, label))
                    actual_atoms.append((speed, label))
        require(full == (1 << q) - 1, ("multiatom cover", q, time, owners))
        require(tuple(actual_atoms) == expected_atoms, ("multiatom atoms", q, tuple(actual_atoms)))
        require(len({speed for speed, _ in actual_atoms}) < len(actual_atoms), ("no repeat", q))
        require(cochain_witness(q, expected_atoms) is not None, ("multiatom cochain", q))
        require(direct_lifted_interval_witness(q, expected_atoms) is not None, ("multiatom direct", q))
    return controls


def source_hash():
    data = Path(__file__).read_bytes().replace(bytes((13, 10)), bytes((10,)))
    return sha256(data).hexdigest()


def main():
    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "assert node")
    require(
        not any(isinstance(node, ast.Constant) and isinstance(node.value, float) for node in ast.walk(tree)),
        "floating literal",
    )

    local_checks, local_profile = exact_local_multiplicity_audit()
    atom_audit = exhaustive_atom_tuple_audit()
    owner_audit = exhaustive_owner_cover_audit()
    boundaries = strict_boundary_controls()
    multiatom_covers = multiatom_full_cover_controls()

    semantic = sha256()
    for item in (
        (local_checks, local_profile),
        atom_audit,
        owner_audit,
        boundaries,
        multiatom_covers,
    ):
        semantic.update(repr(item).encode("ascii") + bytes((10,)))
    digest = semantic.hexdigest()
    if EXPECTED_SEMANTIC != "TO_BE_FROZEN":
        require(digest == EXPECTED_SEMANTIC, ("semantic", digest))

    print("ALL-SHEET ATOMIZED COSET-COCHAIN THM-3402 PRIMARY VERIFIER")
    print(f"source_sha256_lf={source_hash()}")
    print("status=PROVED_analytic_controls_plus_declared_FINITE-EXACT_universes;THM-3402")
    print(f"local_q2_q10_multiplier_phase_checks={local_checks};local_profile={local_profile}")
    print(
        "atom_tuple_checks="
        f"{atom_audit[0]};positive={atom_audit[1]};"
        f"repeated_owner={atom_audit[2]};repeated_owner_positive={atom_audit[3]}"
    )
    print(f"atom_tuple_rank_profile={atom_audit[4]}")
    print(f"owner_cover_checks={owner_audit[0]};positive={owner_audit[1]};profile={owner_audit[2]}")
    print(f"strict_boundary_controls={boundaries}")
    print(f"multiatom_full_cover_controls={multiatom_covers}")
    print("body_relative=pointwise_B_subset_A;aligned_grid_B_minus_A_subset_Gamma;not_B_empty")
    print("scope=no_ledger_decrement,no_physical_transport,no_LRC14")
    print(f"semantic_sha256={digest}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()

