#!/usr/bin/env python3
"""Exact audit for THM-3701's radial successor-mass gate.

The proof itself is the composition of THM-2365's successor-marginal
nonconstancy implication with the thirteen shared masses of THM-3670.  This
companion verifies the two finite linear frames, every packet-count statement,
and the improved exact energy invoice on THM-3672's typed control.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from itertools import combinations, product
from math import gcd
from pathlib import Path
import sys


sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
P = 13
Q = 6
DIM = 1 + 2 * Q
CONTROL_DENOMINATOR = 50334435734703120

PINS = {
    "04-computation/lrc_successor_mass_all_pair_swap_control_thm3672.py":
        "a191f934b20494b98e878fc5504a328c47d4cb6a92ea332246f782fac80c01c8",
    "05-knowledge/results/lrc_successor_mass_all_pair_swap_control_thm3672.out":
        "0d72ff662c908d9d9934976edce67895135c17a1eb0ba409e0abc382563b52c6",
}

# THM-3672 successor numerators, ordered as S0, A_0,...,A_5, B_0,...,B_5.
CONTROL = (
    110219232915792,
    99299969997228,
    99299969997228,
    99299969997228,
    99050735597814,
    96962693739912,
    96754696164126,
    113880033524816,
    113880033524816,
    113880033524816,
    113525486798282,
    111210186545380,
    111203438664916,
)


def require(condition, payload):
    if condition is not True:
        raise RuntimeError(payload)


def dot(left, right):
    return sum(x * y for x, y in zip(left, right))


def matvec(matrix, vector):
    return tuple(dot(row, vector) for row in matrix)


def transpose(matrix):
    return tuple(tuple(row[j] for row in matrix) for j in range(len(matrix[0])))


def gram(matrix):
    columns = transpose(matrix)
    return tuple(tuple(dot(x, y) for y in columns) for x in columns)


def rank(matrix):
    work = [[Fraction(value) for value in row] for row in matrix]
    if not work:
        return 0
    pivot_row = 0
    for column in range(len(work[0])):
        pivot = next(
            (row for row in range(pivot_row, len(work)) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        scale = work[pivot_row][column]
        work[pivot_row] = [value / scale for value in work[pivot_row]]
        for row in range(len(work)):
            if row == pivot_row or not work[row][column]:
                continue
            scale = work[row][column]
            work[row] = [
                value - scale * pivot_value
                for value, pivot_value in zip(work[row], work[pivot_row])
            ]
        pivot_row += 1
    return pivot_row


def radial_rows():
    # Coordinates are S0,A0,...,A5,B0,...,B5.
    rows = []
    for leaf in range(1, DIM):
        row = [0] * DIM
        row[0] = -1
        row[leaf] = 1
        rows.append(tuple(row))
    return tuple(rows)


def forward_rows():
    rows = []
    for k in range(Q):
        for ell in range(Q):
            if k == ell:
                continue
            row = [0] * DIM
            row[0] = 1
            row[1 + k] = 1
            row[1 + Q + ell] = -2
            rows.append(tuple(row))
    return tuple(rows)


def reverse_rows():
    rows = []
    for k in range(Q):
        for ell in range(Q):
            if k == ell:
                continue
            row = [0] * DIM
            row[0] = 1
            row[1 + k] = -2
            row[1 + Q + ell] = 1
            rows.append(tuple(row))
    return tuple(rows)


def eigencheck(matrix, vector, eigenvalue):
    return matvec(matrix, vector) == tuple(eigenvalue * x for x in vector)


def squared_distance_to_scalars(vector):
    mean = Fraction(sum(vector), DIM)
    return sum((Fraction(value) - mean) ** 2 for value in vector)


def frame_energy(rows, vector):
    return sum(value * value for value in matvec(rows, vector))


def matrices_equal(left, right):
    return tuple(tuple(row) for row in left) == tuple(tuple(row) for row in right)


def matrix_add(left, right):
    return tuple(
        tuple(x + y for x, y in zip(left_row, right_row))
        for left_row, right_row in zip(left, right)
    )


def valuation(value, prime):
    answer = 0
    while value % prime == 0:
        value //= prime
        answer += 1
    return answer


def main():
    for relative, expected in PINS.items():
        actual = sha256((ROOT / relative).read_bytes()).hexdigest()
        require(actual == expected, ("parent hash drift", relative, actual))

    radial = radial_rows()
    forward = forward_rows()
    reverse = reverse_rows()
    doubled = forward + reverse
    radial_gram = gram(radial)
    doubled_gram = gram(doubled)

    require((len(radial), len(forward), len(reverse)) == (12, 30, 30), "row counts")
    require(rank(radial) == 12, "radial rank")
    require(rank(forward) == 11, "forward rank")
    require(rank(doubled) == 12, "doubled rank")

    scalar = (1,) * DIM
    radial_basis = [scalar]
    for leaf in range(1, DIM - 1):
        vector = [0] * DIM
        vector[leaf] = 1
        vector[DIM - 1] = -1
        radial_basis.append(tuple(vector))
        require(eigencheck(radial_gram, vector, 1), ("radial eigenvalue 1", leaf))
    radial_top = (-12,) + (1,) * 12
    radial_basis.append(radial_top)
    require(eigencheck(radial_gram, scalar, 0), "radial scalar kernel")
    require(eigencheck(radial_gram, radial_top, 13), "radial eigenvalue 13")
    require(rank(radial_basis) == DIM, "radial eigenbasis")

    doubled_basis = [scalar]
    for sign, eigenvalue in ((1, 29), (-1, 21)):
        for index in range(Q - 1):
            vector = [0] * DIM
            vector[1 + index] = 1
            vector[1 + Q - 1] = -1
            vector[1 + Q + index] = sign
            vector[1 + Q + Q - 1] = -sign
            doubled_basis.append(tuple(vector))
            require(
                eigencheck(doubled_gram, vector, eigenvalue),
                ("doubled standard eigenvalue", sign, index),
            )
    doubled_45 = (0,) + (-1,) * Q + (1,) * Q
    doubled_65 = (-12,) + (1,) * 12
    doubled_basis.extend((doubled_45, doubled_65))
    require(eigencheck(doubled_gram, scalar, 0), "doubled scalar kernel")
    require(eigencheck(doubled_gram, doubled_45, 45), "doubled eigenvalue 45")
    require(eigencheck(doubled_gram, doubled_65, 65), "doubled eigenvalue 65")
    require(rank(doubled_basis) == DIM, "doubled eigenbasis")

    # Exhaust all vectors in {-1,0,1}^13 with support at most three.  This is
    # independent finite control for both sharp frame inequalities.
    frame_controls = 0
    for support_size in range(4):
        for support in combinations(range(DIM), support_size):
            for signs in product((-1, 1), repeat=support_size):
                vector = [0] * DIM
                for index, sign in zip(support, signs):
                    vector[index] = sign
                distance = squared_distance_to_scalars(vector)
                radial_energy = frame_energy(radial, vector)
                doubled_energy = frame_energy(doubled, vector)
                require(distance <= radial_energy <= 13 * distance, ("radial frame", vector))
                require(21 * distance <= doubled_energy <= 65 * distance, ("doubled frame", vector))
                frame_controls += 1

    # For fixed first graft k there are five legal second grafts and four
    # omitted labels; the same count holds after fixing the second graft.
    legal_counts = []
    for fixed in range(Q):
        legal = [
            (other, omitted)
            for other in range(Q)
            for omitted in range(Q)
            if other != fixed and omitted not in (fixed, other)
        ]
        legal_counts.append(len(legal))
    require(legal_counts == [20] * Q, legal_counts)

    derangements = []
    for permutation in __import__("itertools").permutations(range(Q)):
        if all(permutation[k] != k for k in range(Q)):
            derangements.append(permutation)
            covered_a = set(range(Q))
            covered_b = set(permutation)
            require(covered_a == covered_b == set(range(Q)), permutation)
            pairwise_rows = []
            cross_rows = []
            for k, ell in enumerate(permutation):
                center_a = [0] * DIM
                center_a[0] = 1
                center_a[1 + k] = -1
                a_b = [0] * DIM
                a_b[1 + k] = 1
                a_b[1 + Q + ell] = -1
                b_center = [0] * DIM
                b_center[1 + Q + ell] = 1
                b_center[0] = -1
                pairwise_rows.extend((center_a, a_b, b_center))
                cross_rows.append(a_b)
            require(
                matrices_equal(
                    gram(pairwise_rows),
                    matrix_add(radial_gram, gram(cross_rows)),
                ),
                ("derangement star identity", permutation),
            )
    require(len(derangements) == 265, len(derangements))

    radial_control = tuple(value - CONTROL[0] for value in CONTROL[1:])
    require(all(value != 0 for value in radial_control), "control radial zero")
    leaf = max(range(12), key=lambda index: abs(radial_control[index]))
    require(leaf == 5, ("strongest leaf", leaf))
    delta = Fraction(radial_control[leaf], CONTROL_DENOMINATOR)

    # A two-site difference has mask norm kappa=2.  THM-3674's lawful
    # diagonal-zero estimate then gives these exact p=13 tariffs.
    drift_floor = delta * delta / (2 * P**3 * (P - 1))
    target_floor = delta * delta / (2 * P**4 * (P - 1))
    require(
        drift_floor
        == Fraction(
            20682636670561347,
            15240344288762454200781478400,
        ),
        drift_floor,
    )
    require(
        target_floor
        == Fraction(
            20682636670561347,
            198124475753911904610159219200,
        ),
        target_floor,
    )

    old_delta = Fraction(20786137969714, CONTROL_DENOMINATOR)
    old_drift_floor = old_delta * old_delta / (6 * P**3 * (P - 1))
    require(drift_floor > old_drift_floor, (drift_floor, old_drift_floor))
    improvement = drift_floor / old_drift_floor

    # Three prescribed values of one legal packet give a stronger joint
    # variance invoice.  The least variance of a 169-entry table containing
    # x0,x1,x2 is their pairwise squared dispersion divided by 3*169.
    joint_rows = []
    for k in range(Q):
        for ell in range(Q):
            if k == ell:
                continue
            values = (CONTROL[0], CONTROL[1 + k], CONTROL[1 + Q + ell])
            pairwise = sum(
                (values[i] - values[j]) ** 2
                for i, j in combinations(range(3), 2)
            )
            joint_rows.append((k, ell, pairwise))
    joint_k, joint_ell, joint_square_numerator = max(
        joint_rows,
        key=lambda row: row[2],
    )
    require((joint_k, joint_ell) == (5, 0), (joint_k, joint_ell))
    joint_square = Fraction(
        joint_square_numerator,
        CONTROL_DENOMINATOR * CONTROL_DENOMINATOR,
    )
    joint_drift_floor = joint_square / (3 * P**3 * (P - 1))
    joint_target_floor = joint_square / (3 * P**4 * (P - 1))
    require(
        joint_drift_floor
        == Fraction(
            360926324521774869394441,
            148212992112761067316289860457462400,
        ),
        joint_drift_floor,
    )
    require(
        joint_target_floor
        == Fraction(
            360926324521774869394441,
            1926768897465893875111768185947011200,
        ),
        joint_target_floor,
    )
    require(joint_drift_floor > drift_floor > old_drift_floor, "invoice order")

    # THM-2367 (19) gives a sharp all-packet bare-owner hostile.  Every
    # adjacent quotient in this divisibility chain is a multiple of seven,
    # so all lower d/g/u_2 factors peel at their Haar means for arbitrary
    # packet phases.  Only the deep d*g profile J(r-t) remains.
    hostile_units = (1,) + tuple(7**power for power in range(1, 6))
    hostile_c1 = P * 7 * hostile_units[-1]
    hostile_c2 = hostile_c1 * 7 * P
    hostile_c3 = hostile_c2 * 7 * P**3
    hostile_chain = hostile_units + (hostile_c1, hostile_c2, hostile_c3)
    hostile_ratios = tuple(
        right // left for left, right in zip(hostile_chain, hostile_chain[1:])
    )
    require(
        all(right % left == 0 for left, right in zip(hostile_chain, hostile_chain[1:])),
        hostile_chain,
    )
    require(all(ratio % 7 == 0 for ratio in hostile_ratios), hostile_ratios)
    require(tuple(value % P for value in hostile_units) == (1, 7, 10, 5, 9, 11), hostile_units)
    require(
        tuple(valuation(value, P) for value in (hostile_c1, hostile_c2, hostile_c3))
        == (1, 2, 5),
        (hostile_c1, hostile_c2, hostile_c3),
    )
    require(gcd(*hostile_chain) == 1, hostile_chain)
    hostile_prefactor = Fraction(5, 7) * Fraction(6, 7) ** 5 * Fraction(1, 7) * Fraction(6, 7)
    hostile_j = tuple(
        Fraction(0)
        if difference == 0
        else Fraction(1, 13)
        if difference in (1, P - 1)
        else Fraction(1, 7)
        for difference in range(P)
    )
    hostile_mass = hostile_prefactor * sum(hostile_j)
    hostile_uncovered = Fraction(5, 7) * Fraction(6, 7) ** 8
    require(hostile_mass == Fraction(33592320, 524596891), hostile_mass)
    require(hostile_uncovered == Fraction(8398080, 40353607), hostile_uncovered)

    print("radial_rows=12 rank=12 nullity=1")
    print("radial_spectrum=0^1,1^11,13^1")
    print("radial_frame=dist2<=energy<=13*dist2")
    print("forward_rows=30 rank=11 nullity=2")
    print("forward_plus_reverse_rows=60 rank=12 nullity=1")
    print("forward_plus_reverse_spectrum=0^1,21^5,29^5,45^1,65^1")
    print("doubled_frame=21*dist2<=energy<=65*dist2")
    print(f"sparse_frame_controls={frame_controls}")
    print("legal_packets_per_fixed_radial_graft=20")
    print("minimal_derangement_chart_banks=265")
    print("derangement_star_transfer=max_joint_Q>=R_star/6")
    print("control_all_12_radial_differences_nonzero=True")
    print(f"control_strongest=A5-S0={delta}")
    print(f"control_D_floor={drift_floor}")
    print(f"control_E_dt_floor={target_floor}")
    print(f"old_best_D_floor={old_drift_floor}")
    print(f"improvement_ratio={improvement}")
    print(f"control_strongest_joint_chart=({joint_k},{joint_ell})")
    print(f"control_joint_D_floor={joint_drift_floor}")
    print(f"control_joint_E_dt_floor={joint_target_floor}")
    print(f"bare_owner_hostile_units={hostile_units}")
    print(f"bare_owner_hostile_blockers={(hostile_c1, hostile_c2, hostile_c3)}")
    print(f"bare_owner_hostile_mass={hostile_mass}")
    print(f"bare_owner_hostile_uncovered_mass={hostile_uncovered}")
    print("bare_owner_hostile_scope=all_packets_not_delayed_not_cover")
    print("PASS")


if __name__ == "__main__":
    main()
