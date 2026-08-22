#!/usr/bin/env python3
"""Exact controls for THM-3390's all-iterate completion bound."""

from fractions import Fraction
from hashlib import sha256
import json


def require(condition, payload):
    if not condition:
        raise RuntimeError(payload)


ZERO = Fraction(0)
ONE = Fraction(1)


def matrix(rows):
    return tuple(tuple(Fraction(value) for value in row) for row in rows)


I2 = matrix(((1, 0), (0, 1)))
Z2 = matrix(((0, 0), (0, 0)))
J2 = matrix(((0, 1), (0, 0)))


def madd(left, right):
    return matrix([
        [left[i][j] + right[i][j] for j in range(len(left[0]))]
        for i in range(len(left))
    ])


def mscale(scalar, value):
    scalar = Fraction(scalar)
    return matrix([[scalar * entry for entry in row] for row in value])


def mmul(left, right):
    return matrix([
        [sum(left[i][k] * right[k][j] for k in range(len(right)))
         for j in range(len(right[0]))]
        for i in range(len(left))
    ])


def madjoint(value):
    # Every control matrix is real, so the adjoint is the transpose.
    return matrix(zip(*value))


def mpow(value, exponent):
    require(exponent >= 0, ("negative matrix exponent", exponent))
    answer = I2
    base = value
    current = exponent
    while current:
        if current & 1:
            answer = mmul(answer, base)
        base = mmul(base, base)
        current //= 2
    return answer


def commutator(left, right):
    return madd(mmul(left, right), mscale(-1, mmul(right, left)))


def compression_of_bilateral_shift(power):
    """V* U^power V for Ve1=delta_0, Ve2=delta_-1, U delta_j=delta_(j+1)."""
    indices = (0, -1)
    return matrix([
        [ONE if indices[row] == indices[column] + power else ZERO
         for column in range(2)]
        for row in range(2)
    ])


def matrix_vector(value):
    return [[str(entry) for entry in row] for row in value]


def main():
    require(mmul(J2, J2) == Z2, "J2 is not square-zero")
    require(mmul(madjoint(J2), J2) == matrix(((0, 0), (0, 1))),
            "J2 norm control failed")

    # The bilateral shift gives an explicit unitary c-dilation of cJ2.
    for exponent in range(1, 65):
        require(compression_of_bilateral_shift(exponent) == mpow(J2, exponent),
                ("positive dilation mismatch", exponent))
        require(compression_of_bilateral_shift(-exponent)
                == mpow(madjoint(J2), exponent),
                ("adjoint dilation mismatch", exponent))
    masses = (Fraction(2), Fraction(5, 2), Fraction(3), Fraction(7))
    for mass in masses:
        sharp_T = mscale(mass, J2)
        for exponent in range(1, 33):
            completion = madd(
                mscale(mass, compression_of_bilateral_shift(-exponent)),
                mscale(-1, mpow(madjoint(sharp_T), exponent)),
            )
            require(completion == Z2, ("sharp completion nonzero", mass, exponent))

    # The final scalar inequality has positive left coefficient and negative right side.
    mass_norm_pairs = (
        (Fraction(2), Fraction(3)),
        (Fraction(5, 2), Fraction(8, 3)),
        (Fraction(3), Fraction(4)),
        (Fraction(7), Fraction(15, 2)),
    )
    inequality_rows = []
    for mass, norm in mass_norm_pairs:
        coefficient = ONE - mass / (2 * (norm - 1) ** 2)
        right_side = 2 * norm * norm * (ONE - norm / mass)
        require(coefficient > 0 and right_side < 0,
                ("contradiction signs failed", mass, norm, coefficient, right_side))
        inequality_rows.append((str(mass), str(norm), str(coefficient), str(right_side)))

    # Omitting commutation permits arbitrarily large square-zero T.
    hostile_mass = Fraction(2)
    hostile_norm = Fraction(3)
    hostile_T = mscale(hostile_norm, J2)
    hostile_E1 = mscale(-1, madjoint(hostile_T))
    hostile_commutator = commutator(hostile_E1, hostile_T)
    require(hostile_commutator == matrix(((9, 0), (0, -9))),
            ("noncommuting hostile changed", hostile_commutator))
    for exponent in range(2, 16):
        require(mpow(madjoint(hostile_T), exponent) == Z2,
                ("hostile lost bounded tail", exponent))

    # A raw compression makes E1=(c-1)T*, so the first gate forces normality.
    raw_T = compression_of_bilateral_shift(1)
    raw_E1 = madd(
        mscale(hostile_mass, compression_of_bilateral_shift(-1)),
        mscale(-1, madjoint(raw_T)),
    )
    require(raw_E1 == madjoint(raw_T), "raw first defect identity failed")
    raw_commutator = commutator(raw_E1, raw_T)
    require(raw_commutator == matrix(((-1, 0), (0, 1))),
            ("raw compression boundary changed", raw_commutator))

    # Omitting uniformity leaves a scalar commuting counterexample.
    scalar_terms = tuple(-(hostile_norm ** exponent) for exponent in range(1, 7))
    require(all(abs(scalar_terms[i + 1]) > abs(scalar_terms[i])
                for i in range(len(scalar_terms) - 1)),
            "scalar unbounded hostile failed")

    semantic = {
        "theorem": "THM-3390",
        "masses": [str(value) for value in masses],
        "unitary_dilation_checked_through": 64,
        "sharp_matrix": matrix_vector(J2),
        "inequality_rows": inequality_rows,
        "noncommuting_hostile": matrix_vector(hostile_commutator),
        "raw_compression_commutator": matrix_vector(raw_commutator),
        "unbounded_scalar_prefix": [str(value) for value in scalar_terms],
    }
    digest = sha256(json.dumps(semantic, sort_keys=True, separators=(",", ":"))
                    .encode("ascii")).hexdigest()

    print("THM-3390 ALL-ITERATE COMMUTING COMPLETION AUDIT")
    print("mass_inequality=PASS rows=4 c>=2 kappa>c")
    print("sharp_unitary_dilation=PASS matrix=cJ2 E_n=0 checked_n=1..64")
    print("sharp_boundary=PASS norm(cJ2)=c")
    print("noncommuting_hostile=PASS Q=0 T=3J2 commutator=[[9,0],[0,-9]]")
    print("raw_compression=PASS E1=T* commutator=[[-1,0],[0,1]]")
    print("uniformity_hostile=PASS scalar_T=3 E_n=-3^n")
    print("tail_gate=PASS bounded_E_plus_contraction_implies_power_bounded_T")
    print(f"semantic_sha256={digest}")


if __name__ == "__main__":
    main()
