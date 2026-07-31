#!/usr/bin/env python3
"""Exact finite-field certificate for THM-2634.

The theorem is an algebraic typing/no-go result.  This companion verifies:

* invariance of c=floor(169{x}) mod 13 under the THM-2625 quotient shift;
* the two-carry correlation / character factorization on the full endpoint
  plane over an exact field containing a primitive thirteenth root;
* recovery of the old determinant bank in the neutral carry character;
* a sharp full-endpoint/full-determinant hostile with zero matched carry;
* exact pair-twist coefficient extraction; and
* the transverse affine-section criterion on every q!=0 and every v.

No floating-point arithmetic is used.  Every logical gate uses require(), so
optimized Python cannot erase the checks.

Script: 04-computation/lrc14_endpoint_pair_two_carry_cospan_thm2634.py
Output: 05-knowledge/results/lrc14_endpoint_pair_two_carry_cospan_thm2634.out
"""

from collections import Counter
from hashlib import sha256


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


P = 53
ROOT = 16
require(all(P % divisor for divisor in range(2, 8)), "modulus is not prime")
INV13 = pow(13, P - 2, P)
INV169 = INV13 * INV13 % P
require(pow(ROOT, 13, P) == 1 and ROOT != 1, "root has wrong order")
require(all(pow(ROOT, k, P) != 1 for k in range(1, 13)), "root is not primitive")

XY = [(x, y) for x in range(13) for y in range(13)]


def det(a, b):
    return (a[0] * b[1] - a[1] * b[0]) % 13


def digest(values):
    return sha256(",".join(str(value) for value in values).encode()).hexdigest()


def sector_bank(left, right):
    """Unnormalized S(q,Delta)=sum_det left(r+q)right(r)."""
    sectors = [0] * (169 * 13)
    for q_index, q in enumerate(XY):
        for r_index, r in enumerate(XY):
            l = ((r[0] + q[0]) % 13, (r[1] + q[1]) % 13)
            l_index = 13 * l[0] + l[1]
            delta = det(q, r)
            sectors[q_index * 13 + delta] = (
                sectors[q_index * 13 + delta] + left[l_index] * right[r_index]
            ) % P
    return sectors


def difference_bank(left, right):
    """C_d=(1/13)sum_(det fibre,c) L_(c+d)R_c."""
    bank = [0] * (169 * 13 * 13)
    for q_index, q in enumerate(XY):
        for r_index, r in enumerate(XY):
            l = ((r[0] + q[0]) % 13, (r[1] + q[1]) % 13)
            l_index = 13 * l[0] + l[1]
            delta = det(q, r)
            offset = (q_index * 13 + delta) * 13
            for c in range(13):
                rv = right[r_index][c]
                if not rv:
                    continue
                for d in range(13):
                    bank[offset + d] = (
                        bank[offset + d] + INV13 * left[l_index][(c + d) % 13] * rv
                    ) % P
    return bank


def character_bank(left, right):
    """Chat_k=sum_det Lhat_k(r+q)Rhat_k(r), normalized in carry."""
    left_hat = [[0] * 169 for _ in range(13)]
    right_hat = [[0] * 169 for _ in range(13)]
    for k in range(13):
        for point in range(169):
            left_hat[k][point] = INV13 * sum(
                left[point][c] * pow(ROOT, (-k * c) % 13, P) for c in range(13)
            ) % P
            right_hat[k][point] = INV13 * sum(
                right[point][c] * pow(ROOT, (k * c) % 13, P) for c in range(13)
            ) % P

    bank = [0] * (169 * 13 * 13)
    for q_index, q in enumerate(XY):
        for r_index, r in enumerate(XY):
            l = ((r[0] + q[0]) % 13, (r[1] + q[1]) % 13)
            l_index = 13 * l[0] + l[1]
            delta = det(q, r)
            offset = (q_index * 13 + delta) * 13
            for k in range(13):
                bank[offset + k] = (
                    bank[offset + k] + left_hat[k][l_index] * right_hat[k][r_index]
                ) % P
    return bank


def normalized_dft_difference(difference):
    output = [0] * len(difference)
    for sector in range(169 * 13):
        offset = sector * 13
        for k in range(13):
            output[offset + k] = INV13 * sum(
                difference[offset + d] * pow(ROOT, (-k * d) % 13, P)
                for d in range(13)
            ) % P
    return output


def inverse_difference(character):
    output = [0] * len(character)
    for sector in range(169 * 13):
        offset = sector * 13
        for d in range(13):
            output[offset + d] = sum(
                character[offset + k] * pow(ROOT, (k * d) % 13, P)
                for k in range(13)
            ) % P
    return output


def support(values):
    return sum(value != 0 for value in values)


def main():
    print("THM-2634 exact endpoint-pair two-carry cospan certificate")
    print(f"field=F_{P}; primitive13root={ROOT}; endpoint_plane=F13^2")

    # Right-continuous endpoint digit on an exact toy N=169*T_DEN grid.  A
    # representative shift x->x+s/13 adds 13s*T_DEN to the numerator.
    carry_checks = 0
    toy_denominator = 17
    toy_order = 169 * toy_denominator
    for j in range(169):
        for remainder in (0, 1, toy_denominator - 1):
            numerator = j * toy_denominator + remainder
            for s in range(13):
                shifted = (numerator + 13 * s * toy_denominator) % toy_order
                original_carry = (numerator // toy_denominator) % 13
                shifted_carry = (shifted // toy_denominator) % 13
                require(shifted_carry == original_carry, "endpoint carry gauge failure")
                carry_checks += 1
    print(f"right-continuous carry quotient-gauge checks={carry_checks}: PASS")

    # Deterministic dense control for the exact factorization.
    left = [
        [((3 * x + 5 * y + 7 * c + x * y + 2 * c * x + 1) % P) for c in range(13)]
        for x, y in XY
    ]
    right = [
        [((11 * x + 13 * y + 17 * c + 3 * x * y + c * y + 2) % P) for c in range(13)]
        for x, y in XY
    ]

    difference = difference_bank(left, right)
    character = character_bank(left, right)
    require(normalized_dft_difference(difference) == character, "carry DFT factorization")
    require(inverse_difference(character) == difference, "carry DFT inversion")

    total_left = [sum(row) % P for row in left]
    total_right = [sum(row) % P for row in right]
    old_sectors = sector_bank(total_left, total_right)
    neutral_character = [character[13 * sector] for sector in range(169 * 13)]
    require(
        neutral_character == [INV169 * value % P for value in old_sectors],
        "neutral carry character is not the old determinant bank",
    )
    print(
        "dense factorization: support(old,difference,character,matched)="
        f"({support(old_sectors)},{support(difference)},{support(character)},"
        f"{support(difference[13 * sector] for sector in range(169 * 13))})"
    )
    print(
        "dense sha256(old,difference,character)="
        f"({digest(old_sectors)},{digest(difference)},{digest(character)})"
    )

    # Sharp full-support hostile.  The old endpoint/determinant array and both
    # singleton norms are unchanged when the right carry moves from 0 to 1,
    # but the matched-carry contraction changes from full support to zero.
    left_single = [[1] + [0] * 12 for _ in XY]
    right_matched = [[1] + [0] * 12 for _ in XY]
    right_shifted = [[0, 1] + [0] * 11 for _ in XY]
    totals = [1] * 169
    old_full = sector_bank(totals, totals)
    base_support = support(old_full)
    require(base_support == 2185, "full determinant hostile support")

    matched_difference = difference_bank(left_single, right_matched)
    shifted_difference = difference_bank(left_single, right_shifted)
    require(
        all(
            matched_difference[13 * sector] == old_full[sector] * INV13 % P
            for sector in range(169 * 13)
        ),
        "matched hostile lost diagonal",
    )
    require(
        all(shifted_difference[13 * sector] == 0 for sector in range(169 * 13)),
        "shifted hostile retained diagonal",
    )
    require(
        all(
            shifted_difference[13 * sector + d] == (old_full[sector] * INV13 % P if d == 12 else 0)
            for sector in range(169 * 13) for d in range(13)
        ),
        "shifted hostile difference support",
    )
    shifted_character = character_bank(left_single, right_shifted)
    require(support(shifted_character) == 2185 * 13, "shifted hostile character support")
    print(
        "sharp hostile: old endpoint cells=28561; old determinant sectors=2185; "
        "matched supports=(2185,0); shifted character support=28405"
    )

    # Exact pair-twist coefficient extraction in the t-variable.  The second
    # cross coefficient is deliberately unrelated; orthogonality removes it.
    pair_twist_checks = 0
    for cross in range(P):
        other = (7 * cross + 3) % P
        norm = (5 * cross + 11) % P
        intensities = [
            (norm + pow(ROOT, -t % 13, P) * cross + pow(ROOT, t, P) * other) % P
            for t in range(13)
        ]
        recovered = INV13 * sum(
            intensities[t] * pow(ROOT, t, P) for t in range(13)
        ) % P
        require(recovered == cross, "pair-twist quadrature")
        pair_twist_checks += 1
    print(f"pair-twist exact quadrature checks={pair_twist_checks}: PASS")

    # Exhaust the transverse section law for every nonzero q and every v.
    section_checks = 0
    alpha_histogram = Counter()
    for q in XY:
        if q == (0, 0):
            continue
        for v in XY:
            alpha = det(q, v)
            alpha_histogram[alpha] += 1
            r0 = ((3 * q[0] + 5 * v[0] + 1) % 13, (7 * q[1] + 2 * v[1] + 4) % 13)
            beta = det(q, r0)
            for c in range(13):
                r = ((r0[0] + c * v[0]) % 13, (r0[1] + c * v[1]) % 13)
                delta = det(q, r)
                require(delta == (beta + alpha * c) % 13, "affine section law")
                if alpha:
                    recovered = (delta - beta) * pow(alpha, 11, 13) % 13
                    require(recovered == c, "transverse carry inversion")
                section_checks += 1
    require(alpha_histogram[0] == 168 * 13, "parallel-section count")
    require(all(alpha_histogram[a] == 168 * 13 for a in range(1, 13)), "transverse-section count")
    print(
        f"affine section checks={section_checks}; per-alpha sections={168 * 13}; "
        "alpha=0 parallel, alpha!=0 invertible: PASS"
    )
    print("single-c typing: ABSENT without pair twist/common section")
    print("all exact checks passed")


if __name__ == "__main__":
    main()
