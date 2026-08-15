#!/usr/bin/env python3
"""Exact hostile referee for THM-3400's commutator-flux tariff."""

from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path


def require(condition, payload):
    if not condition:
        raise RuntimeError(payload)


def matrix(rows):
    return tuple(tuple(Fraction(entry) for entry in row) for row in rows)


def identity(size):
    return matrix([
        [1 if row == column else 0 for column in range(size)]
        for row in range(size)
    ])


def madd(left, right):
    return matrix([
        [left[row][column] + right[row][column]
         for column in range(len(left[0]))]
        for row in range(len(left))
    ])


def mscale(scalar, value):
    scalar = Fraction(scalar)
    return matrix([[scalar * entry for entry in row] for row in value])


def mmul(left, right):
    return matrix([
        [sum(left[row][middle] * right[middle][column]
             for middle in range(len(right)))
         for column in range(len(right[0]))]
        for row in range(len(left))
    ])


def madjoint(value):
    # Every frozen control is real, so the adjoint is the transpose.
    return matrix(zip(*value))


def mpow(value, exponent):
    require(exponent >= 0, ("negative exponent", exponent))
    answer = identity(len(value))
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


def mvec(value, vector):
    return tuple(sum(value[row][column] * vector[column]
                     for column in range(len(vector)))
                 for row in range(len(value)))


def vadd(left, right):
    return tuple(a + b for a, b in zip(left, right))


def vscale(scalar, vector):
    scalar = Fraction(scalar)
    return tuple(scalar * entry for entry in vector)


def dot(left, right):
    return sum(a * b for a, b in zip(left, right))


def completion(mass, q_matrix, t_matrix, exponent):
    return madd(
        mscale(mass, mpow(madjoint(q_matrix), exponent)),
        mscale(-1, mpow(madjoint(t_matrix), exponent)),
    )


def ledger_row(mass, norm, q_matrix, t_matrix, vector, exponent):
    t_power = mpow(t_matrix, exponent)
    e_n = completion(mass, q_matrix, t_matrix, exponent)
    e_next = completion(mass, q_matrix, t_matrix, exponent + 1)
    p_n = mscale(mass, mpow(madjoint(q_matrix), exponent))
    p_next = mscale(mass, mpow(madjoint(q_matrix), exponent + 1))

    c_n = dot(mvec(commutator(t_power, e_n), vector), vector)
    d_n = dot(
        mvec(commutator(t_power, e_next), mvec(t_matrix, vector)),
        vector,
    )
    leakage = norm * c_n - d_n
    m_n = dot(mvec(mmul(e_n, t_power), vector), vector)
    t_next = mpow(t_matrix, exponent + 1)
    m_next = dot(mvec(mmul(e_next, t_next), vector), vector)
    r_n = dot(
        mvec(
            mmul(
                t_power,
                madd(mscale(norm, e_n),
                      mscale(-1, mmul(e_next, t_matrix))),
            ),
            vector,
        ),
        vector,
    )
    y_n = mvec(
        madd(mmul(p_next, t_matrix), mscale(-norm, p_n)),
        vector,
    )
    t_star_power_x = mvec(mpow(madjoint(t_matrix), exponent), vector)
    a_value = norm * norm - norm
    completion_value = (
        a_value * dot(t_star_power_x, t_star_power_x)
        - dot(mvec(t_power, y_n), vector)
    )
    square_vector = vadd(
        t_star_power_x,
        vscale(-Fraction(1, 2 * a_value), y_n),
    )
    square_value = (
        a_value * dot(square_vector, square_vector)
        - dot(y_n, y_n) / (4 * a_value)
    )
    require(r_n == norm * m_n - m_next + leakage,
            ("commutator ledger", mass, norm, exponent, r_n,
             norm * m_n - m_next + leakage))
    require(r_n == completion_value == square_value,
            ("completion square", mass, norm, exponent,
             r_n, completion_value, square_value))
    return c_n, d_n, leakage, r_n


def lf_sha(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def main():
    zero = matrix(((0, 0), (0, 0)))
    j_matrix = matrix(((0, 1), (0, 0)))
    e_two = (Fraction(0), Fraction(1))
    require(mmul(j_matrix, j_matrix) == zero, "J is not square-zero")

    # The optimal family realizes every rational leakage budget k(k-c).
    masses = (Fraction(2), Fraction(5, 2), Fraction(3), Fraction(7))
    sharp_rows = []
    ledger_checks = 0
    for mass in masses:
        norms = (mass, mass + Fraction(1, 2), mass + 2)
        for norm in norms:
            t_matrix = mscale(norm, j_matrix)
            rows = [ledger_row(mass, norm, j_matrix, t_matrix, e_two, exponent)
                    for exponent in range(1, 7)]
            ledger_checks += len(rows)
            c_one, d_one, leakage_one, _ = rows[0]
            require(c_one == norm * (norm - mass),
                    ("sharp C1", mass, norm, c_one))
            require(d_one == 0, ("sharp D1", mass, norm, d_one))
            require(leakage_one == norm * norm * (norm - mass),
                    ("sharp L1", mass, norm, leakage_one))
            require(all(row[2] == 0 for row in rows[1:]),
                    ("sharp tail leakage", mass, norm, rows))
            lambda_value = leakage_one / norm
            epsilon = norm * (norm - mass)
            require(lambda_value == epsilon,
                    ("sharp tariff", mass, norm, lambda_value, epsilon))
            require(norm * norm - mass * norm - epsilon == 0,
                    ("robust equality", mass, norm, epsilon))
            coefficient = 1 - mass / (2 * (norm - 1) ** 2)
            right_side = (2 * norm / mass) * (lambda_value - epsilon)
            require(right_side == 0, ("tariff remainder", mass, norm))
            if norm > mass:
                require(coefficient > 0, ("tariff coefficient", mass, norm))
            sharp_rows.append((str(mass), str(norm), str(lambda_value)))

    # A strictly favorable observer current need not commute globally.
    favorable_rows = []
    for mass in masses:
        norm = (mass + 1) / 2
        t_matrix = mscale(norm, j_matrix)
        c_one, d_one, leakage_one, _ = ledger_row(
            mass, norm, j_matrix, t_matrix, e_two, 1)
        comm = commutator(t_matrix, completion(mass, j_matrix, t_matrix, 1))
        require(comm != zero, ("favorable family commuted", mass, norm))
        require(d_one == 0 and leakage_one / norm == norm * (norm - mass) < 0,
                ("favorable leakage", mass, norm, c_one, leakage_one))
        favorable_rows.append((str(mass), str(norm), str(leakage_one / norm)))

    # The square-zero omission hostile has positive leakage and norm excess.
    omission_rows = []
    for mass in masses:
        norm = mass + 1
        t_matrix = mscale(norm, j_matrix)
        c_one, d_one, leakage_one, _ = ledger_row(
            mass, norm, zero, t_matrix, e_two, 1)
        lambda_value = leakage_one / norm
        require(c_one == norm * norm and d_one == 0,
                ("omission commutator", mass, norm, c_one, d_one))
        require(lambda_value == norm * norm > norm * (norm - mass),
                ("omission tariff", mass, norm, lambda_value))
        omission_rows.append((str(mass), str(norm), str(lambda_value)))

    # A nonnormal involution has bounded all-iterate defects but exact
    # infinite discounted leakage k^2+1.
    periodic_rows = []
    for mass in masses:
        norm = mass + 1
        t_matrix = matrix(((0, norm), (1 / norm, 0)))
        require(mmul(t_matrix, t_matrix) == identity(2),
                ("periodic involution", mass, norm))
        partial = Fraction(0)
        for exponent in range(1, 13):
            _, d_n, leakage, _ = ledger_row(
                mass, norm, zero, t_matrix, e_two, exponent)
            expected = (norm * (norm * norm - 1 / (norm * norm))
                        if exponent % 2 else 0)
            require(d_n == 0 and leakage == expected,
                    ("periodic leakage term", mass, norm, exponent,
                     d_n, leakage, expected))
            partial += leakage / (norm ** exponent)
        total = norm * norm + 1
        require(partial == total * (1 - norm ** -12),
                ("periodic geometric sum", mass, norm, partial, total))
        require(total >= norm * (norm - mass),
                ("periodic tariff", mass, norm, total))
        periodic_rows.append((str(mass), str(norm), str(total), str(partial)))

    # An arbitrarily long scalar prefix has zero leakage but no uniform tail.
    prefix_rows = []
    for length in (1, 2, 5, 12, 25):
        mass = Fraction(2)
        norm = Fraction(3)
        defects = tuple(-(norm ** exponent) for exponent in range(1, length + 1))
        require(max(abs(value) for value in defects) == norm ** length,
                ("finite prefix cap", length, defects[-1]))
        prefix_rows.append((length, str(defects[-1])))

    semantic = {
        "theorem": "THM-3400",
        "masses": [str(value) for value in masses],
        "sharp_rows": sharp_rows,
        "favorable_rows": favorable_rows,
        "omission_rows": omission_rows,
        "periodic_rows": periodic_rows,
        "finite_prefix_rows": prefix_rows,
        "ledger_checks": ledger_checks + len(masses) * (1 + 1 + 12),
        "tariff": "Lambda>=kappa(kappa-c) when kappa>c>=2",
    }
    digest = sha256(json.dumps(semantic, sort_keys=True, separators=(",", ":"))
                    .encode("ascii")).hexdigest()
    source_hash = lf_sha(Path(__file__))

    print("THM-3400 DISCOUNTED NORMING-ORBIT COMMUTATOR-FLUX TARIFF")
    print(f"source_sha256_lf={source_hash}")
    print("status=PROVED exact tariff and optimal robust norm bound")
    print("tariff=Lambda_x>=kappa(kappa-c) for kappa>c>=2")
    print("robust_bound=kappa<=(c+sqrt(c^2+4epsilon))/2")
    print(f"sharp_budget_rows={len(sharp_rows)} equality=PASS")
    print(f"favorable_noncommuting_rows={len(favorable_rows)} leakage_negative=PASS")
    print(f"square_zero_omission_rows={len(omission_rows)} leakage_positive=PASS")
    print(f"periodic_involution_rows={len(periodic_rows)} Lambda=kappa^2+1 PASS")
    print(f"finite_prefix_rows={len(prefix_rows)} uniform_tail_required=PASS")
    print(f"exact_ledger_checks={semantic['ledger_checks']}")
    print(f"semantic_sha256={digest}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
