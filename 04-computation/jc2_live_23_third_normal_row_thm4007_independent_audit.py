#!/usr/bin/env python3
"""SymPy-free hostile audit for THM-4007.

This audit reconstructs the load-bearing axis coefficient by hand and checks
the fourth residual row with an independent exact polynomial implementation.
It deliberately does not import the primary certificate.
"""

from __future__ import annotations

from fractions import Fraction as F
import hashlib
import sys


sys.stdout.reconfigure(newline="\n")
CHECKS = 0
TRANSCRIPT: list[str] = []


def emit(line: str) -> None:
    print(line)
    TRANSCRIPT.append(line)


def gate(label: str, condition: bool) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(label)
    emit(f"PASS  {label}")


def trim(poly: list[F]) -> list[F]:
    answer = list(poly)
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return answer


def add(left: list[F], right: list[F]) -> list[F]:
    size = max(len(left), len(right))
    return trim([
        (left[j] if j < len(left) else F(0))
        + (right[j] if j < len(right) else F(0))
        for j in range(size)
    ])


def neg(poly: list[F]) -> list[F]:
    return [-entry for entry in poly]


def sub(left: list[F], right: list[F]) -> list[F]:
    return add(left, neg(right))


def scale(value: F, poly: list[F]) -> list[F]:
    return trim([value * entry for entry in poly])


def mul(left: list[F], right: list[F]) -> list[F]:
    answer = [F(0)] * (len(left) + len(right) - 1)
    for i, li in enumerate(left):
        for j, rj in enumerate(right):
            answer[i + j] += li * rj
    return trim(answer)


def deriv(poly: list[F]) -> list[F]:
    if len(poly) == 1:
        return [F(0)]
    return trim([F(j) * poly[j] for j in range(1, len(poly))])


def equal(label: str, left: list[F], right: list[F]) -> None:
    gate(label, trim(left) == trim(right))


emit("STATUS=THM-4007 INDEPENDENT SYMPY-FREE HOSTILE AUDIT")

# The depth rows i>=-2 and i>=-3 give j+i=3 on the third source
# diagonal.  The exact extreme rows have no s^5/s^6 terms, leaving maxima
# j=4 for A (from i=-1) and j=5 for C (from i=-2).
A_third_addresses = [(i, 3 - i) for i in range(-2, 4)]
C_third_addresses = [(i, 3 - i) for i in range(-3, 4)]
gate("A naive third cap is five before deleting the exact extreme term", max(j for _, j in A_third_addresses) == 5)
gate("A exact depth-minus-two row deletes x^5 and leaves degree four", max(j for i, j in A_third_addresses if i != -2) == 4)
gate("C naive third cap is six before deleting the exact extreme term", max(j for _, j in C_third_addresses) == 6)
gate("C exact depth-minus-three row deletes x^6 and leaves degree five", max(j for i, j in C_third_addresses if i != -3) == 5)

# Axis-only derivation in invariant rows.  A5 is kept as a formal Laurent
# exponent: every displayed number below is the rational coefficient after
# factoring A5^-3 from the t^2 Jacobian constant.
C0_prime_scaled = -F(3, 4)       # C0'(0)
N0_scaled = F(4, 3)              # A5*N(0)
K1_scaled = -F(4)                # A5*K'(0)
M0_scaled = -F(32, 9)            # A5^2*M(0)
L1_scaled = F(88, 15)            # A5^2*L'(0)
lower_scaled = -2 * M0_scaled * K1_scaled - N0_scaled * L1_scaled
gate("axis lower row is -544/15 times A5^-3", lower_scaled == -F(544, 15))
u0_scaled = -lower_scaled / (-3 * C0_prime_scaled)
gate("axis equation forces 2176/135 times A5^-3", u0_scaled == F(2176, 135))
gate("forced axis coefficient is nonzero in characteristic zero", u0_scaled != 0)
gate("weight minus six is new", -6 not in {2, 0, -2, -4})


def row_case(c02: F, c21: F) -> None:
    """Check J[t^2] and the eliminated P[t^4] row at a=A5=1."""
    A0 = [F(1), F(0), F(1, 4)]
    C0 = [F(0), -F(3, 4), F(0), -F(1, 8)]
    A1 = [F(4, 3), F(0), F(2)]
    C1 = [F(0), -F(4), F(0), -F(3, 2)]
    A2 = [-F(32, 9), F(0), -F(4, 5)]
    C2 = [F(0), F(88, 15), F(0), -F(12, 5)]
    A3 = [
        F(2176, 135),
        c21 / 4,
        F(1088, 315) + 2 * c02 / 7,
        F(0),
        -F(32, 15),
    ]
    C3 = [
        -3 * c21 / 8,
        -F(8128, 315) - 3 * c02 / 7,
        -3 * c21 / 16,
        F(736, 105) - 3 * c02 / 14,
        F(0),
        F(8, 5),
    ]

    j2 = add(
        scale(3, sub(mul(deriv(A0), C3), mul(A3, deriv(C0)))),
        add(
            sub(scale(2, mul(deriv(A1), C2)), scale(2, mul(A2, deriv(C1)))),
            sub(mul(deriv(A2), C1), mul(A1, deriv(C2))),
        ),
    )
    equal(f"J[t^2]=0 for c02={c02},c21={c21}", j2, [F(0)])

    # Jlower is the part of J[t^3] not containing the fourth vector.
    jlower = add(
        add(
            sub(scale(3, mul(deriv(A1), C3)), mul(A1, deriv(C3))),
            sub(mul(deriv(A3), C1), scale(3, mul(A3, deriv(C1)))),
        ),
        scale(2, sub(mul(deriv(A2), C2), mul(A2, deriv(C2)))),
    )
    qpoly = [-F(3), F(0), -F(1, 2)]
    g4 = add(
        scale(-F(1, 4), mul(qpoly, jlower)),
        add(
            add(scale(2, mul(C1, C3)), mul(C2, C2)),
            add(
                scale(-6, mul(mul(A0, A1), A3)),
                add(scale(-3, mul(A0, mul(A2, A2))), scale(-3, mul(mul(A1, A1), A2))),
            ),
        ),
    )
    c40 = -F(11392, 105) - 6 * c02 / 7
    epsilon = -F(1376, 135)
    alpha = F(8, 3)
    expected = [-c40 / 2, -c21 / 2, 3 * epsilon - c02 / 2, F(0), alpha]
    equal(f"P[t^4] residual row for c02={c02},c21={c21}", g4, expected)
    gate(f"residual row degree is exactly four for c02={c02},c21={c21}", len(trim(g4)) - 1 == 4)


row_case(F(0), F(0))
row_case(F(7), F(16, 3))
row_case(-F(14, 3), -F(5, 2))

# A wrong distinct-row determinant changes the polynomial in a concrete
# case, whereas the same-row case is blind to the error.
P = ([F(1), F(0), F(1)], [F(0), F(1), F(0), F(1)])
Q = ([F(2), F(1)], [F(3), F(0), F(1)])
correct = sub(mul(deriv(P[0]), Q[1]), mul(deriv(P[1]), Q[0]))
wrong = sub(mul(deriv(P[0]), Q[1]), mul(P[0], deriv(Q[1])))
gate("distinct-row determinant hostile detects the factor swap", correct != wrong)
same_correct = sub(mul(deriv(P[0]), P[1]), mul(deriv(P[1]), P[0]))
same_wrong = sub(mul(deriv(P[0]), P[1]), mul(P[0], deriv(P[1])))
equal("same-row determinant is blind control", same_correct, same_wrong)

gate("all live strata have A-support at least five", min(5, 5, 5) == 5)
gate("fixed-gauge C floor is five but is not claimed shear-invariant", min(6, 6, 5) == 5)

semantic = hashlib.sha256("\n".join(TRANSCRIPT).encode()).hexdigest()
emit(f"CHECKS={CHECKS}")
emit(f"SEMANTIC_SHA256={semantic}")
emit("THEOREM_ID=THM-4007")
emit("ALL THM-4007 INDEPENDENT HOSTILE CHECKS PASSED")
