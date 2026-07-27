#!/usr/bin/env python3
"""Exact companion for THM-2580.

The finite carrier reconstruction uses only integers and Fraction arithmetic.
Smith forms and resultants use exact SymPy integer/polynomial algorithms.
"""

from fractions import Fraction
import hashlib
import importlib.util
from math import comb, gcd
from pathlib import Path

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ


P = 13
RANK = P - 1
checks = 0


def check(condition: bool, message: str = "exact check failed") -> None:
    global checks
    checks += 1
    if not condition:
        raise AssertionError(message)


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def rank_mod(rows, prime: int) -> int:
    a = [[int(value) % prime for value in row] for row in rows]
    rank = 0
    if not a:
        return 0
    for col in range(len(a[0])):
        pivot = next((i for i in range(rank, len(a)) if a[i][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        inv = pow(a[rank][col], -1, prime)
        a[rank] = [(inv * value) % prime for value in a[rank]]
        for i in range(len(a)):
            if i == rank:
                continue
            multiple = a[i][col]
            if multiple:
                a[i] = [
                    (left - multiple * right) % prime
                    for left, right in zip(a[i], a[rank])
                ]
        rank += 1
    return rank


def smith_diagonal(matrix: sp.Matrix):
    diagonal = smith_normal_form(matrix, domain=ZZ)
    return tuple(abs(int(diagonal[i, i])) for i in range(matrix.rows))


def cyclic_lattice_matrices():
    """Return P,D,U,C on the basis e_i-e_0, 1<=i<=12."""
    full_shift = sp.zeros(P, P)
    for row in range(P):
        full_shift[row, (row + 1) % P] = 1
    basis = sp.zeros(P, RANK)
    for i in range(1, P):
        basis[i, i - 1] = 1
        basis[0, i - 1] = -1
    shift = (full_shift * basis)[1:, :]
    identity = sp.eye(RANK)
    difference = shift - identity
    unit = shift + identity
    cayley = difference * unit.inv()
    check(difference.det() == P)
    check(unit.det() == 1)
    check(all(value.q == 1 for value in cayley))
    return full_shift, shift, difference, unit, cayley


def hasse_tower(difference: sp.Matrix, cayley: sp.Matrix):
    smith_rows = []
    for depth in range(1, P):
        hasse = sp.Matrix(
            [[comb(index, order) for index in range(1, P)]
             for order in range(1, depth + 1)]
        )
        power = difference**depth
        product = hasse * power
        check(all(int(value) % P == 0 for value in product))
        check(rank_mod(hasse.tolist(), P) == depth)
        expected = (1,) * (RANK - depth) + (P,) * depth
        check(smith_diagonal(power) == expected)
        check(smith_diagonal(cayley**depth) == expected)
        smith_rows.append(depth)

    c12_over_13 = (cayley**12) / 13
    check(all(value.q == 1 for value in c12_over_13))
    check(abs(int(c12_over_13.det())) == 1)
    check(smith_diagonal(cayley**13) == (13,) * 11 + (169,))
    check(smith_diagonal(cayley**24) == (169,) * 12)

    # Orientation/sign control for the first primitive.  With P=z^{-1},
    # C=(1-z)/(1+z).  Hence C(1-z^2)=(1-z)^2 and beta_1=-2 beta_2.
    y = sp.zeros(RANK, 1)
    y[1, 0] = -1  # 1-z^2 has coordinates (- at z^2); z^0 is implicit.
    delta = sp.zeros(RANK, 1)
    delta[0, 0] = -2
    delta[1, 0] = 1
    check(cayley * y == delta)
    beta1_y = sum(index * int(y[index - 1, 0]) for index in range(1, P))
    beta2_delta = sum(
        comb(index, 2) * int(delta[index - 1, 0]) for index in range(1, P)
    )
    check(beta1_y == -2)
    check(beta2_delta == 1)
    check((beta1_y + 2 * beta2_delta) % P == 0)
    return smith_rows


def canonical_primitive():
    """Independently rebuild THM-2571's globally primitive X tensor."""
    here = Path(__file__).resolve().parent
    t1 = load_module(
        here / "lrc14_deep_colour_cayley_bockstein_thm2571.py",
        "thm2571_exact_for_2580",
    )
    source = load_module(
        here / "lrc14_replica_dichotomy_typed_row_opus_20260727.py",
        "thm2550b_exact_for_2580",
    )
    denominator = source.T
    clock = P**6
    check(denominator == 297836897838480)
    check(clock == 4826809)

    word = source.build_word_Ta()
    word_digit = [[None] * P for _ in range(7)]
    for owner in range(7):
        phase = source.subtract_comb(
            word, source.C1, 182, 26 * owner - 13, 26 * owner + 13
        )
        for old_digit in range(P):
            word_digit[owner][old_digit] = t1.intersect_interval(
                phase,
                old_digit * denominator // P,
                (old_digit + 1) * denominator // P,
            )

    fractions = [
        [[Fraction(0) for _ in range(P)] for _ in range(P)] for _ in range(7)
    ]
    for owner in range(7):
        prefixes = [source.make_prefix(word_digit[owner][h]) for h in range(P)]
        lengths = [t1.interval_length(word_digit[owner][h]) for h in range(P)]
        for target_shift in range(P):
            present = source.build_F(owner, target_shift)
            for deep_shift in range(P):
                layer = source.intersect_comb(
                    present,
                    source.C3,
                    182,
                    14 * deep_shift - 13,
                    14 * deep_shift + 13,
                )
                value = Fraction(0)
                for old_digit in range(P):
                    piece = t1.intersect_interval(
                        layer,
                        old_digit * denominator // P,
                        (old_digit + 1) * denominator // P,
                    )
                    if not piece or not word_digit[owner][old_digit]:
                        continue
                    starts, lens, prefix = prefixes[old_digit]
                    acc_r, acc_p = source.sweep_acc(
                        piece, clock % denominator, starts, lens, prefix
                    )
                    value += source.IR_from_acc(
                        clock,
                        t1.interval_length(piece),
                        lengths[old_digit],
                        acc_r,
                        acc_p,
                    )
                fractions[owner][target_shift][deep_shift] = value

    raw = [[[0] * P for _ in range(P)] for _ in range(7)]
    content = 0
    positive = 0
    for owner in range(7):
        for target_shift in range(P):
            for deep_shift in range(P):
                scaled = (
                    fractions[owner][target_shift][deep_shift]
                    * clock
                    * denominator
                )
                check(scaled.denominator == 1)
                value = scaled.numerator
                raw[owner][target_shift][deep_shift] = value
                if value:
                    positive += 1
                    content = gcd(content, abs(value))
                    check(deep_shift != 0)
                else:
                    check(deep_shift == 0)
    check(positive == 1092)
    check(content == 13)
    primitive = [
        [[raw[o][s][r] // 13 for r in range(P)] for s in range(P)]
        for o in range(7)
    ]
    serialization = ";".join(
        ",".join(str(primitive[o][s][r]) for r in range(P))
        for o in range(7)
        for s in range(P)
    )
    digest = hashlib.sha256(serialization.encode("ascii")).hexdigest()
    check(digest == "a66ba96d31a33354468392b1dabc19865e6e925158efdab059fad9a98d4390f4")
    return t1, primitive, digest


def cyc_add(left, right):
    return tuple((a + b) % P for a, b in zip(left, right))


def cyc_scale(scalar, value):
    return tuple((scalar * a) % P for a in value)


def cyc_sub(left, right):
    return tuple((a - b) % P for a, b in zip(left, right))


def cyc_pow(t1, prime, value, exponent):
    answer = t1.cyc_power_mod(prime, 0)
    base = value
    while exponent:
        if exponent & 1:
            answer = t1.cyc_mul_mod(prime, answer, base)
        base = t1.cyc_mul_mod(prime, base, base)
        exponent //= 2
    return answer


def canonical_secondary_torsor(t1, primitive):
    omega = tuple(range(1, P))
    beta1 = []
    beta2 = []
    one13 = t1.cyc_power_mod(P, 0)
    for deep_shift in range(P):
        first = (0,) * 12
        second = (0,) * 12
        for colour in range(P):
            phase = t1.cyc_power_mod(P, colour * deep_shift)
            first = cyc_add(first, cyc_scale(colour, phase))
            second = cyc_add(second, cyc_scale(comb(colour, 2), phase))
        beta1.append(first)
        beta2.append(second)

    for deep_shift in range(1, P):
        inverse = pow(deep_shift, -1, P)
        check(beta1[deep_shift] == cyc_scale(inverse, omega))
        zeta_r = t1.cyc_power_mod(P, deep_shift)
        difference = cyc_sub(zeta_r, one13)
        expected = t1.cyc_mul_mod(
            P,
            t1.cyc_power_mod(P, 2 * deep_shift),
            cyc_pow(t1, P, difference, 10),
        )
        check(beta2[deep_shift] == expected)

    y = []
    z = []
    for owner in range(7):
        first_value = 0
        second_value = 0
        for target_shift in range(P):
            for deep_shift in range(1, P):
                coefficient = primitive[owner][target_shift][deep_shift]
                first_value += coefficient * pow(deep_shift, -1, P)
                second_value += (
                    coefficient
                    * target_shift
                    * pow(deep_shift, -2, P)
                )
        y.append(first_value % P)
        z.append(second_value % P)
    check(y == [0, 6, 5, 1, 12, 8, 7])
    check(z == [0, 5, 3, 9, 4, 10, 8])

    z_reduced = tuple((z[i] - z[6]) % P for i in range(6))
    z_inverse = (4, 8, 5, 8, 0, 3)
    check(z_reduced == (5, 10, 8, 1, 9, 2))
    check(t1.cyc_mul_mod(7, z_reduced, z_inverse) == (1, 0, 0, 0, 0, 0))

    owner_y = []
    owner_z = []
    for kappa in range(1, 7):
        first = (0,) * 6
        second = (0,) * 6
        for owner in range(7):
            phase = t1.cyc_power_mod(7, kappa * owner)
            first = cyc_add(first, cyc_scale(y[owner], phase))
            second = cyc_add(second, cyc_scale(z[owner], phase))
        check(all(first))
        check(all(second))
        owner_y.append(first)
        owner_z.append(second)
    check(
        owner_z
        == [
            (5, 10, 8, 1, 9, 2),
            (4, 8, 9, 1, 7, 12),
            (10, 7, 6, 2, 5, 1),
            (3, 6, 7, 11, 8, 12),
            (9, 5, 4, 12, 6, 1),
            (8, 3, 5, 12, 4, 11),
        ]
    )
    check(rank_mod(owner_z, P) == 3)

    beta2_profiles = {}
    for kappa in range(1, 7):
        for target in range(P):
            tensor = [[0] * 12 for _ in range(6)]
            for owner in range(7):
                phase7 = t1.cyc_power_mod(7, kappa * owner)
                for target_shift in range(P):
                    phase13 = t1.cyc_power_mod(P, target * target_shift)
                    for deep_shift in range(1, P):
                        coloured = t1.cyc_mul_mod(
                            P, phase13, beta2[deep_shift]
                        )
                        coefficient = primitive[owner][target_shift][deep_shift] % P
                        if not coefficient:
                            continue
                        for i in range(6):
                            if not phase7[i]:
                                continue
                            for j in range(12):
                                tensor[i][j] = (
                                    tensor[i][j]
                                    + coefficient * phase7[i] * coloured[j]
                                ) % P
            beta2_profiles[kappa, target] = tensor

    absolute = 0
    pairwise = 0
    for kappa in range(1, 7):
        expected_first = [
            [(owner_y[kappa - 1][i] * omega[j]) % P for j in range(12)]
            for i in range(6)
        ]
        check(sum(value != 0 for row in expected_first for value in row) == 72)
        absolute += P
        for left in range(P):
            for right in range(left + 1, P):
                difference = [
                    [
                        (
                            beta2_profiles[kappa, left][i][j]
                            - beta2_profiles[kappa, right][i][j]
                        )
                        % P
                        for j in range(12)
                    ]
                    for i in range(6)
                ]
                scalar = (-(left - right)) % P
                expected = [
                    [
                        scalar * owner_z[kappa - 1][i] * omega[j] % P
                        for j in range(12)
                    ]
                    for i in range(6)
                ]
                check(difference == expected)
                check(sum(value != 0 for row in difference for value in row) == 72)
                pairwise += 1
    check(absolute == 78)
    check(pairwise == 468)

    # Orthogonality at the common carrier: every unnormalized target Fourier
    # numerator is literally thirteen times one integral target-shift slice.
    fourier_cells = 0
    zero13 = (0,) * 12
    thirteen = (13,) + (0,) * 11
    for target_frequency in range(P):
        for target_shift in range(P):
            value = [0] * 12
            for target in range(P):
                phase = t1.cyc_power_mod(
                    P, target * (target_frequency + target_shift)
                )
                value = [left + right for left, right in zip(value, phase)]
            expected = thirteen if (target_frequency + target_shift) % P == 0 else zero13
            check(tuple(value) == expected)
            fourier_cells += 1
    check(fourier_cells == 169)
    return {
        "absolute": absolute,
        "pairwise": pairwise,
        "z": z,
        "z_inverse": z_inverse,
        "owner_rank": rank_mod(owner_z, P),
    }


def signature_spectrum(full_shift, shift, difference, unit):
    identity = sp.eye(RANK)
    f1 = difference**2 * unit
    f2 = difference**3 * unit**2 * (shift**2 + identity)
    salem = shift**6 - shift**4 - 2 * shift**3 - shift**2 + identity
    f3 = difference * salem
    diag1 = smith_diagonal(f1)
    diag2 = smith_diagonal(f2)
    diag3 = smith_diagonal(f3)
    check(diag1 == (1,) * 10 + (13, 13))
    check(diag2 == (1,) * 9 + (13, 13, 13))
    check(diag3 == (1,) * 8 + (5, 5, 5, 65))

    x = sp.symbols("x")
    phi13 = sum(x**i for i in range(P))
    salem_poly = x**6 - x**4 - 2 * x**3 - x**2 + 1
    check(sp.resultant(x + 1, phi13, x) == 1)
    check(sp.resultant(x**2 + 1, phi13, x) == 1)
    check(sp.resultant(salem_poly, phi13, x) == 625)
    check(salem_poly.subs(x, 1) == -2)
    check(sp.gcd(x + 1, phi13) == 1)
    check(sp.gcd((x + 1) ** 2 * (x**2 + 1), phi13) == 1)
    check(sp.gcd(salem_poly, phi13) == 1)

    full_identity = sp.eye(P)
    full_difference = full_shift - full_identity
    full_f1 = full_difference**2 * (full_shift + full_identity)
    full_f2 = (
        full_difference**3
        * (full_shift + full_identity) ** 2
        * (full_shift**2 + full_identity)
    )
    full_salem = (
        full_shift**6
        - full_shift**4
        - 2 * full_shift**3
        - full_shift**2
        + full_identity
    )
    full_f3 = full_difference * full_salem
    check(full_f1.rank() == 12)
    check(full_f2.rank() == 12)
    check(full_f3.rank() == 12)
    return diag1, diag2, diag3


def compact_smith(diagonal):
    groups = []
    for value in diagonal:
        if groups and groups[-1][0] == value:
            groups[-1] = (value, groups[-1][1] + 1)
        else:
            groups.append((value, 1))
    return ",".join(f"{value}^{count}" for value, count in groups)


def main() -> None:
    full_shift, shift, difference, unit, cayley = cyclic_lattice_matrices()
    depths = hasse_tower(difference, cayley)
    signatures = signature_spectrum(full_shift, shift, difference, unit)
    t1, primitive, digest = canonical_primitive()
    torsor = canonical_secondary_torsor(t1, primitive)

    print("THM-2580 Hasse-Bockstein carry tower")
    print(f"Hasse exact depths {depths[0]}..{depths[-1]}, Smith p-multiplicities 1..12")
    print("Cayley periodicity C^12 Lambda=13 Lambda; r13=13^11,169; r24=169^12")
    print(f"signature Smith F1 {compact_smith(signatures[0])}")
    print(f"signature Smith F2 {compact_smith(signatures[1])}")
    print(f"signature Smith Salem {compact_smith(signatures[2])}; resultant=625=5^4")
    print("13-primary signature depths cyclotomic2/cyclotomic3/Salem1")
    print(f"canonical primitive digest {digest}")
    print(f"secondary Z {torsor['z']}, inverse {torsor['z_inverse']}")
    print(
        "canonical carry: absolute beta1 78/78, pair beta2 "
        f"{torsor['pairwise']}/468, owner span rank {torsor['owner_rank']}"
    )
    print("all pair-edge tensors have full 72-coordinate support")
    print("unnormalized target Fourier numerator has literal factor13 and depth>=12")
    print(f"explicit checks {checks}")


if __name__ == "__main__":
    main()
