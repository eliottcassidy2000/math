#!/usr/bin/env python3
"""Exact high-order formal scout for the exceptional constant normal.

Work directly in K=Q[alpha]/(F_4) with Fraction arithmetic.  For a prescribed
pencil Q_s=Q+s*h, recursively solve the six labelled-triple equations using
the constant implicit Jacobian from the formal-local synthesis.  The sixth
unknown at each order is the compensating x coefficient.  This determines
the first order at which an uncorrected straight pencil bends away from the
triple locus, if such an order occurs within the audited cutoff.
"""

from __future__ import annotations

from fractions import Fraction as F
import hashlib
import sys


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


ORDER = 12
CHECKS = 0
K = tuple[F, F, F, F]


def check(condition: bool, label: str) -> None:
    global CHECKS
    if not condition:
        raise AssertionError(label)
    CHECKS += 1


def kval(value=0) -> K:
    return (F(value), F(0), F(0), F(0))


K_ZERO = kval(0)
K_ONE = kval(1)
K_ALPHA = (F(0), F(1), F(0), F(0))


# alpha^4=(77822208*alpha^3+28419741*alpha^2-7849770*alpha+1276420)/72783360.
ALPHA4 = (
    F(1276420, 72783360),
    F(-7849770, 72783360),
    F(28419741, 72783360),
    F(77822208, 72783360),
)


def kadd(left: K, right: K) -> K:
    return tuple(a + b for a, b in zip(left, right))  # type: ignore[return-value]


def kneg(value: K) -> K:
    return tuple(-item for item in value)  # type: ignore[return-value]


def ksub(left: K, right: K) -> K:
    return kadd(left, kneg(right))


def kscale(value: K, scalar: F | int) -> K:
    scalar = F(scalar)
    return tuple(scalar * item for item in value)  # type: ignore[return-value]


def kmul(left: K, right: K) -> K:
    work = [F(0)] * 7
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            work[i + j] += a * b
    for degree in range(6, 3, -1):
        coefficient = work[degree]
        if coefficient:
            work[degree] = F(0)
            for j, replacement in enumerate(ALPHA4):
                work[degree - 4 + j] += coefficient * replacement
    return tuple(work[:4])  # type: ignore[return-value]


def kpow(value: K, exponent: int) -> K:
    answer = K_ONE
    for _ in range(exponent):
        answer = kmul(answer, value)
    return answer


def kformat(value: K) -> str:
    if value[1:] == (F(0), F(0), F(0)):
        return str(value[0])
    return "[" + ",".join(map(str, value)) + "]"


def szero() -> list[K]:
    return [K_ZERO for _ in range(ORDER + 1)]


def sadd(left: list[K], right: list[K]) -> list[K]:
    return [kadd(a, b) for a, b in zip(left, right)]


def sneg(value: list[K]) -> list[K]:
    return [kneg(item) for item in value]


def ssub(left: list[K], right: list[K]) -> list[K]:
    return sadd(left, sneg(right))


def sscale(value: list[K], scalar: K) -> list[K]:
    return [kmul(item, scalar) for item in value]


def smul(left: list[K], right: list[K]) -> list[K]:
    answer = szero()
    for degree in range(ORDER + 1):
        total = K_ZERO
        for i in range(degree + 1):
            total = kadd(total, kmul(left[i], right[degree - i]))
        answer[degree] = total
    return answer


def spow(value: list[K], exponent: int) -> list[K]:
    answer = szero()
    answer[0] = K_ONE
    for _ in range(exponent):
        answer = smul(answer, value)
    return answer


def sconstant(value: K) -> list[K]:
    answer = szero()
    answer[0] = value
    return answer


def seval(coefficients: list[K], value: list[K]) -> list[K]:
    answer = szero()
    for coefficient in reversed(coefficients):
        answer = sadd(smul(answer, value), sconstant(coefficient))
    return answer


def pmul(left: list[K], right: list[K]) -> list[K]:
    answer = [K_ZERO] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            answer[i + j] = kadd(answer[i + j], kmul(a, b))
    return answer


def padd(left: list[K], right: list[K]) -> list[K]:
    size = max(len(left), len(right))
    answer = [K_ZERO] * size
    for i in range(size):
        answer[i] = kadd(
            left[i] if i < len(left) else K_ZERO,
            right[i] if i < len(right) else K_ZERO,
        )
    return answer


def solve_fraction(matrix: tuple[tuple[F, ...], ...], rhs: tuple[F, ...]) -> tuple[F, ...]:
    size = len(rhs)
    work = [list(matrix[row]) + [rhs[row]] for row in range(size)]
    for column in range(size):
        pivot = next(row for row in range(column, size) if work[row][column])
        work[column], work[pivot] = work[pivot], work[column]
        scale = work[column][column]
        work[column] = [entry / scale for entry in work[column]]
        for row in range(size):
            if row == column:
                continue
            multiple = work[row][column]
            if multiple:
                work[row] = [
                    a - multiple * b
                    for a, b in zip(work[row], work[column])
                ]
    return tuple(work[row][-1] for row in range(size))


def determinant_fraction(matrix: tuple[tuple[F, ...], ...]) -> F:
    """Exact determinant by elimination, independent of the printed constant."""
    size = len(matrix)
    work = [list(row) for row in matrix]
    answer = F(1)
    for column in range(size):
        pivot = next(row for row in range(column, size) if work[row][column])
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            answer = -answer
        pivot_value = work[column][column]
        answer *= pivot_value
        for row in range(column + 1, size):
            multiple = work[row][column] / pivot_value
            for index in range(column, size):
                work[row][index] -= multiple * work[column][index]
    return answer


J = (
    (F(3), 0, 0, -1, 0, -2),
    (F(-9), 0, 0, 0, -1, 2),
    (0, F(3), 0, -1, 0, 0),
    (0, F(4), 0, 0, -1, 0),
    (0, 0, F(3), -1, 0, -2),
    (0, 0, F(9), 0, -1, -2),
)


def solve_k(rhs: tuple[K, ...]) -> tuple[K, ...]:
    coordinate_solutions = [
        solve_fraction(J, tuple(rhs[row][coordinate] for row in range(6)))
        for coordinate in range(4)
    ]
    return tuple(
        tuple(coordinate_solutions[coordinate][row] for coordinate in range(4))
        for row in range(6)
    )  # type: ignore[return-value]


def exceptional_q() -> list[K]:
    p_alpha = kadd(
        kadd(kval(F(-5717, 729)), kscale(K_ALPHA, F(-1688, 81))),
        kscale(kpow(K_ALPHA, 2), F(520, 9)),
    )
    a = kadd(kadd(kval(F(-259, 36)), p_alpha), kscale(K_ALPHA, 4))
    b = kscale(K_ALPHA, -9)
    c = kneg(p_alpha)
    auxiliary = [a, b, c]
    P = [K_ZERO, K_ZERO, K_ONE, K_ZERO, kval(-2), K_ZERO, K_ONE]
    Q1 = [
        kval(F(-3, 4)),
        K_ONE,
        kval(F(-27, 4)),
        kval(-2),
        kval(F(9, 2)),
        K_ONE,
    ]
    return padd(Q1, pmul(P, auxiliary))


Q = exceptional_q()


def compiler_at(
    xseries: list[K],
    compensator: list[K],
    forcing: list[K],
) -> tuple[list[K], list[K]]:
    qseries = seval(Q, xseries)
    s_parameter = szero()
    s_parameter[1] = K_ONE
    qseries = sadd(qseries, smul(s_parameter, seval(forcing, xseries)))
    qseries = sadd(qseries, smul(compensator, xseries))
    D = sadd(sconstant(K_ONE), smul(spow(xseries, 2), qseries))
    C = smul(xseries, smul(D, sadd(D, sconstant(kval(2)))))
    E = smul(qseries, sadd(D, sconstant(kval(3))))
    return C, E


def recurse(forcing: list[K], label: str) -> tuple[list[list[K]], list[K], list[K], list[K]]:
    endpoints = [sconstant(kval(point)) for point in (-1, 0, 1)]
    common_c = sconstant(K_ZERO)
    common_e = sconstant(kval(-3))
    compensator = sconstant(K_ZERO)

    for order in range(1, ORDER + 1):
        residual: list[K] = []
        for endpoint in endpoints:
            C, E = compiler_at(endpoint, compensator, forcing)
            residual.extend((ksub(C[order], common_c[order]), ksub(E[order], common_e[order])))
        solution = solve_k(tuple(kneg(value) for value in residual))
        for index in range(3):
            endpoints[index][order] = solution[index]
        common_c[order] = solution[3]
        common_e[order] = solution[4]
        compensator[order] = solution[5]

        for endpoint in endpoints:
            C, E = compiler_at(endpoint, compensator, forcing)
            check(C[order] == common_c[order], f"{label} C order {order}")
            check(E[order] == common_e[order], f"{label} E order {order}")

    return endpoints, common_c, common_e, compensator


def main() -> None:
    # Sanity checks on the exceptional field relation and base polynomial.
    check(determinant_fraction(J) == -288, "implicit Jacobian determinant")
    alpha4_direct = kpow(K_ALPHA, 4)
    check(alpha4_direct == ALPHA4, "quartic field reduction")
    for point in (-1, 0, 1):
        xseries = sconstant(kval(point))
        qvalue = seval(Q, xseries)[0]
        D = kadd(K_ONE, kmul(kval(point * point), qvalue))
        C = kmul(kval(point), kmul(D, kadd(D, kval(2))))
        E = kmul(qvalue, kadd(D, kval(3)))
        check(C == K_ZERO and E == kval(-3), f"base collision point {point}")

    constant = [K_ONE]
    _, _, _, constant_compensator = recurse(constant, "constant")
    check(constant_compensator[1] == K_ZERO, "constant first compensator")
    check(constant_compensator[2] == K_ZERO, "constant second compensator")
    check(constant_compensator[3] == K_ZERO, "constant third compensator")
    check(constant_compensator[4] == K_ZERO, "constant fourth compensator")
    expected_constant_fifth = (
        F(259188338368, 129140163),
        F(-46584993664, 4782969),
        F(23019960448, 531441),
        F(9180348416, 177147),
    )
    check(
        constant_compensator[5] == expected_constant_fifth,
        "constant fifth compensator regression",
    )
    # THM-4046's order-eight constant-response obstruction in the same basis.
    kappa_order_eight = (
        F(-5183766767360, 3**19),
        F(931699873280, 3**16),
        F(-460399208960, 3**14),
        F(-183606968320, 3**13),
    )
    check(
        constant_compensator[5] == kscale(kappa_order_eight, F(-9, 20)),
        "constant fifth compensator equals -9/20 times THM-4046 kappa",
    )

    quadratic = [K_ZERO, kval(-9), kval(4)]
    _, _, _, quadratic_compensator = recurse(quadratic, "quadratic")
    check(quadratic_compensator[1] == K_ZERO, "quadratic first compensator")
    check(quadratic_compensator[2] != K_ZERO, "quadratic second-order hostile")
    expected_quadratic_second = (F(-1089575, 1296), F(-11696, 9), F(2080), F(0))
    check(
        quadratic_compensator[2] == expected_quadratic_second,
        "quadratic second-order synthesis control",
    )

    first_constant_nonzero = next(
        (order for order in range(1, ORDER + 1) if constant_compensator[order] != K_ZERO),
        None,
    )
    check(first_constant_nonzero == 5, "constant exact first nonzero order")
    payload = ";".join(kformat(value) for value in constant_compensator)
    print("exceptional_quartic_constant_normal_high_order_scout")
    print(f"field=Q[alpha]/F4;formal_cutoff={ORDER};arithmetic=Fraction_tuple_mod_F4")
    print("implicit_unknowns=three_endpoints,common_C,common_E,x_compensator")
    print("implicit_J_determinant=-288;recursion=exact_order_by_order")
    for order in range(1, ORDER + 1):
        print(f"constant_x_compensator_s^{order}={kformat(constant_compensator[order])}")
    print(f"constant_first_nonzero_through_cutoff={first_constant_nonzero}")
    print(
        "constant_s5_common_denominator="
        "64*(4049817787-19653044202*alpha+87403912326*alpha^2+"
        "104569906176*alpha^3)/3^17"
    )
    print("constant_s5=-9/20*THM4046_order8_kappa")
    print("constant_straight_pencil=collision_mod_s^5;not_collision_mod_s^6")
    print(f"constant_compensator_hash={hashlib.sha256(payload.encode('ascii')).hexdigest()}")
    print(f"quadratic_x_compensator_s^1={kformat(quadratic_compensator[1])}")
    print(f"quadratic_x_compensator_s^2={kformat(quadratic_compensator[2])}")
    print("scope=formal_labelled_triple_at_exceptional_fibre;finite_order_scout_only")
    print("NO_CLAIM=all_order_constant_pencil_or_polynomial_termination_or_global_graph_or_JC2")
    print(f"checks={CHECKS};result=PASS")


if __name__ == "__main__":
    main()
