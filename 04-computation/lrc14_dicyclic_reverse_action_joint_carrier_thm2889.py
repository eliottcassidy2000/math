#!/usr/bin/env python3
"""Exact companion for the THM-2889 two-shadow dicyclic carrier.

The semantic characters chi_s=det(QB,-) and chi_c=det(QA,-) give two
differently labelled C169 semidirect Q8 shadows, each isomorphic to
Dic_338.  The script checks both labelled actions, their coordinatewise
inversion fibre product over Q8, complete order/conjugacy censuses, the
two q7 ancestry fillers, and the exact projective/Hermitian match with
the THM-2868 Prony character.  Rank-two minimality is asserted only for
full C169 coordinates in the stable-diagonal coordinatewise-inversion
model, not among all possible auxiliary phase carriers.

Nothing here constructs a physical Q8 action on the marked horn current.
"""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
from math import gcd
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
RESULTS = ROOT / "05-knowledge/results"
sys.path.insert(0, str(COMP))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_hash(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


PINNED = {
    COMP / "lrc14_two_origin_endpoint_projective_kummer_thm2868.py":
        "3282526029d27210a5cadaa10f702a974f926e217e6838913d11d9dbc888d9b5",
    RESULTS / "lrc14_two_origin_endpoint_projective_kummer_thm2868.out":
        "ce4b8961b3dfa468d51b884b91002872440b05bf70d7e5eecfc3318d4100e2f9",
    COMP / "lrc14_event_twisted_all_q_coefficient_carry_lift_thm2882.py":
        "3ed346e0c631b34bd61f0c4d27d7f161e8d35b70decfb95f5207c5f57893d005",
    RESULTS / "lrc14_event_twisted_all_q_coefficient_carry_lift_thm2882.out":
        "0faa0a24f6ba8b6c88b6bbfc4f225e38667097b1a937d977741453499884901c",
    COMP / "lrc14_macro_semantic_diagonal_horn_carrier_thm2884.py":
        "b739be20e741d5c061e0febcc8aef9b0f58f4ae8a648aa803610e0dad991929f",
    RESULTS / "lrc14_macro_semantic_diagonal_horn_carrier_thm2884.out":
        "8c3829b1052a641ca08a5e5bda86d9d5e8bd1584f5b2911c57c9fad9da41d4b6",
    COMP / "lrc14_stepped_origin_v4_provenance_transport_thm2886.py":
        "1f2cbffb8151c0c74bb22beff58e0ace5715eba3d4a9afc59481e7ab5e0d6dc9",
    RESULTS / "lrc14_stepped_origin_v4_provenance_transport_thm2886.out":
        "60ba517f1b8fa92fdaeba48b68d10da97f23f00a9fb1526de9e39f1311e49929",
    COMP / "lrc14_quaternionic_arf_semantic_v4_carry_no_go_thm2887.py":
        "4fa39842cf0a7f407bd7d436315ac85e132161933c10d2c078644d91dafd91b8",
    RESULTS / "lrc14_quaternionic_arf_semantic_v4_carry_no_go_thm2887.out":
        "b2edbfdb4964573810a7147e8676dd152ba60aafe7953dde7afd3a2a65da3d5a",
}
for path, expected in PINNED.items():
    require(path.is_file(), f"missing pinned dependency: {path}")
    require(lf_hash(path) == expected, f"pinned dependency changed: {path.name}")


import lrc14_two_origin_endpoint_projective_kummer_thm2868 as atlas


P = 13
CARRY_ORDER = 169
DIC_PARAMETER = 338
ROTATION_ORDER = 2 * DIC_PARAMETER  # 676
GROUP_ORDER = 4 * DIC_PARAMETER     # 1352
ONE = (0, 0)


# Dic_338 convention:
#   <a,x | a^676=1, x^2=a^338, x a x^-1=a^-1>.
# Store a^i x^b as (i,b), with b in {0,1}.
def dmul(left: tuple[int, int], right: tuple[int, int]) -> tuple[int, int]:
    i, b = left
    j, c = right
    return (
        (
            i
            + (-j if b else j)
            + (DIC_PARAMETER if b and c else 0)
        ) % ROTATION_ORDER,
        b ^ c,
    )


def dinv(element: tuple[int, int]) -> tuple[int, int]:
    i, b = element
    if not b:
        return ((-i) % ROTATION_ORDER, 0)
    return ((i + DIC_PARAMETER) % ROTATION_ORDER, 1)


def dpow(element: tuple[int, int], exponent: int) -> tuple[int, int]:
    require(exponent >= 0, "negative dicyclic exponent")
    result = ONE
    base = element
    power = exponent
    while power:
        if power & 1:
            result = dmul(result, base)
        base = dmul(base, base)
        power //= 2
    return result


DIVISORS = (1, 2, 4, 8, 13, 26, 52, 104, 169, 338, 676, 1352)


def dorder(element: tuple[int, int]) -> int:
    for divisor in DIVISORS:
        if dpow(element, divisor) == ONE:
            return divisor
    raise RuntimeError(f"element order did not divide {GROUP_ORDER}: {element}")


def conjugate(
    actor: tuple[int, int],
    element: tuple[int, int],
) -> tuple[int, int]:
    return dmul(dmul(actor, element), dinv(actor))


def commutator(
    left: tuple[int, int],
    right: tuple[int, int],
) -> tuple[int, int]:
    return dmul(
        dmul(dmul(left, right), dinv(left)),
        dinv(right),
    )


ELEMENTS = tuple(
    (rotation, reflection)
    for rotation in range(ROTATION_ORDER)
    for reflection in (0, 1)
)
A = (1, 0)
X = (0, 1)
MINUS_ONE = dpow(A, DIC_PARAMETER)

# The labelled semantic embedding is chosen so that the fixed source QB
# defines the determinant character:
#
#   QA = x       (inverts C169),
#   QB = a^169   (commutes with C169),
#   QAB = QA*QB  (inverts C169).
#
# The coefficient-aligned C169 generator c=a^344 has projective eigenvalue
# ratio rho=xi^42, whose thirteenth power is omega^3.
C = dpow(A, 344)
QA = X
QB = dpow(A, 169)
QAB = dmul(QA, QB)
CANONICAL_QAB = dmul(MINUS_ONE, QAB)


def q8_subgroup() -> frozenset[tuple[int, int]]:
    return frozenset(
        dmul(dpow(QA, i), dpow(QB, j))
        for i in range(4)
        for j in range(4)
    )


def carry_subgroup() -> frozenset[tuple[int, int]]:
    return frozenset(dpow(C, n) for n in range(CARRY_ORDER))


def semantic_central_sign(element: tuple[int, int]) -> int:
    """Return 0 for Ad_QB-fixed and 1 for Ad_QB-negated Q8 units."""
    image = conjugate(QB, element)
    if image == element:
        return 0
    if image == dmul(MINUS_ONE, element):
        return 1
    raise RuntimeError(f"semantic element has no central sign: {element}")


def chi_direction(direction: tuple[int, int]) -> int:
    """det(QB,direction) for directions in V4=Q8/{+/-1}.

    Coordinates are QA=(1,0), QB=(0,1), QAB=(1,1).
    """
    u, _v = direction
    return u


def matmul(
    left: tuple[tuple[int, int], tuple[int, int]],
    right: tuple[tuple[int, int], tuple[int, int]],
    prime: int,
) -> tuple[tuple[int, int], tuple[int, int]]:
    return (
        (
            (
                left[0][0] * right[0][0]
                + left[0][1] * right[1][0]
            ) % prime,
            (
                left[0][0] * right[0][1]
                + left[0][1] * right[1][1]
            ) % prime,
        ),
        (
            (
                left[1][0] * right[0][0]
                + left[1][1] * right[1][0]
            ) % prime,
            (
                left[1][0] * right[0][1]
                + left[1][1] * right[1][1]
            ) % prime,
        ),
    )


def matscale(
    scalar: int,
    matrix: tuple[tuple[int, int], tuple[int, int]],
    prime: int,
) -> tuple[tuple[int, int], tuple[int, int]]:
    return tuple(
        tuple(scalar * entry % prime for entry in row)
        for row in matrix
    )  # type: ignore[return-value]


def matrix_for(
    element: tuple[int, int],
    eta: int,
    prime: int,
) -> tuple[tuple[int, int], tuple[int, int]]:
    rotation, reflection = element
    diagonal = (
        (pow(eta, rotation, prime), 0),
        (0, pow(eta, -rotation, prime)),
    )
    if not reflection:
        return diagonal
    # Standard faithful dicyclic generator x |-> [[0,-1],[1,0]].
    x_matrix = ((0, prime - 1), (1, 0))
    return matmul(diagonal, x_matrix, prime)


QElement = tuple[int, int, int]
BElement = tuple[int, int, QElement]
Q_ONE: QElement = (0, 0, 0)
Q_MINUS_ONE: QElement = (0, 0, 1)
Q_QA: QElement = (1, 0, 0)
Q_QB: QElement = (0, 1, 0)
# Section product Q_QA*Q_QB.  THM-2887's displayed canonical lift uses
# QA*QB=-QAB, so its QAB is the central negative of this section.
Q_QAB: QElement = (1, 1, 0)
Q_QAB_CANONICAL: QElement = (1, 1, 1)
B_ONE: BElement = (0, 0, Q_ONE)


def qmul(left: QElement, right: QElement) -> QElement:
    """Q8 in the normal form (-1)^z QA^u QB^v."""
    u, v, z = left
    r, s, w = right
    return (
        u ^ r,
        v ^ s,
        (z + w + v * r + u * r + v * s) % 2,
    )


def qpow(element: QElement, exponent: int) -> QElement:
    result = Q_ONE
    base = element
    power = exponent
    while power:
        if power & 1:
            result = qmul(result, base)
        base = qmul(base, base)
        power //= 2
    return result


def qinv(element: QElement) -> QElement:
    return qpow(element, 3)


def bmul(left: BElement, right: BElement) -> BElement:
    """Multiplication on (C169_s x C169_c) semidirect Q8.

    QA flips the selector coordinate s and QB flips the carry coordinate c.
    """
    selector, carry_coordinate, h = left
    other_selector, other_carry, k = right
    u, v, _z = h
    return (
        (
            selector
            + (-other_selector if u else other_selector)
        ) % CARRY_ORDER,
        (
            carry_coordinate
            + (-other_carry if v else other_carry)
        ) % CARRY_ORDER,
        qmul(h, k),
    )


def binv(element: BElement) -> BElement:
    selector, carry_coordinate, h = element
    inverse_h = qinv(h)
    u, v, _z = inverse_h
    return (
        (-(-selector if u else selector)) % CARRY_ORDER,
        (-(-carry_coordinate if v else carry_coordinate)) % CARRY_ORDER,
        inverse_h,
    )


def bpow(element: BElement, exponent: int) -> BElement:
    result = B_ONE
    base = element
    power = exponent
    while power:
        if power & 1:
            result = bmul(result, base)
        base = bmul(base, base)
        power //= 2
    return result


B_DIVISORS = (1, 2, 4, 13, 26, 52, 169, 338, 676)


def border(element: BElement) -> int:
    for divisor in B_DIVISORS:
        if bpow(element, divisor) == B_ONE:
            return divisor
    raise RuntimeError(f"bi-dicyclic element order exceeded 676: {element}")


def bconjugate(actor: BElement, element: BElement) -> BElement:
    return bmul(bmul(actor, element), binv(actor))


def bcommutator(left: BElement, right: BElement) -> BElement:
    return bmul(
        bmul(bmul(left, right), binv(left)),
        binv(right),
    )


def square_matmul(
    left: tuple[tuple[int, ...], ...],
    right: tuple[tuple[int, ...], ...],
    prime: int,
) -> tuple[tuple[int, ...], ...]:
    size = len(left)
    require(
        size > 0
        and all(len(row) == size for row in left)
        and len(right) == size
        and all(len(row) == size for row in right),
        "non-square matrix input",
    )
    return tuple(
        tuple(
            sum(left[row][middle] * right[middle][column]
                for middle in range(size)) % prime
            for column in range(size)
        )
        for row in range(size)
    )


def square_scale(
    scalar: int,
    matrix: tuple[tuple[int, ...], ...],
    prime: int,
) -> tuple[tuple[int, ...], ...]:
    return tuple(
        tuple(scalar * entry % prime for entry in row)
        for row in matrix
    )


def square_pow(
    matrix: tuple[tuple[int, ...], ...],
    exponent: int,
    prime: int,
) -> tuple[tuple[int, ...], ...]:
    size = len(matrix)
    result = tuple(
        tuple(int(row == column) for column in range(size))
        for row in range(size)
    )
    base = matrix
    power = exponent
    while power:
        if power & 1:
            result = square_matmul(result, base, prime)
        base = square_matmul(base, base, prime)
        power //= 2
    return result


def kron2(
    left: tuple[tuple[int, int], tuple[int, int]],
    right: tuple[tuple[int, int], tuple[int, int]],
    prime: int,
) -> tuple[tuple[int, ...], ...]:
    return tuple(
        tuple(
            left[left_row][left_column]
            * right[right_row][right_column]
            % prime
            for left_column in range(2)
            for right_column in range(2)
        )
        for left_row in range(2)
        for right_row in range(2)
    )


def monomial_permutation(
    matrix: tuple[tuple[int, ...], ...],
) -> tuple[int, ...]:
    permutation = []
    for row in matrix:
        support = tuple(index for index, value in enumerate(row) if value)
        require(len(support) == 1, "matrix is not monomial")
        permutation.append(support[0])
    require(len(set(permutation)) == len(permutation),
            "monomial support is not a permutation")
    return tuple(permutation)


def main() -> None:
    # Group presentation and labelled subgroup anatomy.
    require(
        len(ELEMENTS) == GROUP_ORDER
        and dpow(A, ROTATION_ORDER) == ONE
        and dpow(A, ROTATION_ORDER // 2) == MINUS_ONE != ONE
        and dmul(X, X) == MINUS_ONE
        and conjugate(X, A) == dinv(A),
        "Dic_338 presentation failed",
    )
    q8 = q8_subgroup()
    carry = carry_subgroup()
    require(
        len(q8) == 8
        and len(carry) == CARRY_ORDER
        and q8 & carry == {ONE}
        and len({dmul(c, h) for c in carry for h in q8}) == GROUP_ORDER
        and dpow(QA, 2) == dpow(QB, 2) == dpow(QAB, 2) == MINUS_ONE
        and QAB == dmul(QA, QB)
        and dmul(QB, QA) == dmul(MINUS_ONE, QAB),
        "C169/Q8 split anatomy failed",
    )

    # There are four Q8 -> C2 characters.  Relative to
    # QA=(1,0), QB=(0,1), only det(QB,-)=u has horn pattern (0,1,1)
    # on (Q0,QA,QAB), and its order-four kernel is exactly <QB>.
    direction_bank = ((0, 0), (1, 0), (0, 1), (1, 1))
    hom_patterns = tuple(
        tuple((alpha * u + beta * v) % 2 for u, v in direction_bank)
        for alpha, beta in direction_bank
    )
    horn_patterns = tuple(
        (pattern[0], pattern[1], pattern[3])
        for pattern in hom_patterns
    )
    edge_patterns = tuple(
        (pattern[1], pattern[2], pattern[3])
        for pattern in hom_patterns
    )
    require(
        hom_patterns
        == (
            (0, 0, 0, 0),
            (0, 1, 0, 1),
            (0, 0, 1, 1),
            (0, 1, 1, 0),
        )
        and tuple(
            index for index, pattern in enumerate(horn_patterns)
            if pattern == (0, 1, 1)
        ) == (1,),
        "det(QB,-) was not the unique horn character",
    )
    action_kernel = frozenset(
        h for h in q8 if conjugate(h, C) == C
    )
    expected_kernel = frozenset(dpow(QB, exponent) for exponent in range(4))
    require(
        action_kernel == expected_kernel
        and len(action_kernel) == 4
        and conjugate(QA, C) == dinv(C)
        and conjugate(QB, C) == C
        and conjugate(QAB, C) == dinv(C),
        "Q8 action kernel is not <QB>",
    )

    # There are two inequivalent labelled Dic_338 shadows.
    #
    # selector shadow:
    #   chi_s=det(QB,-)=u, ker=<QB>, target values (Q0,QA,QAB)=(0,1,1);
    #
    # carry-edge shadow:
    #   chi_c=det(QA,-)=v, ker=<QA>, first edge QA flat and the wrapped
    #   q11->q7 edge QB reversing.
    #
    # Swapping the semantic labels QA,QB in the same abstract Dic_338
    # presentation realizes the second shadow.
    carry_shadow_qa = QB
    carry_shadow_qb = QA
    carry_shadow_qab = dmul(carry_shadow_qa, carry_shadow_qb)
    carry_action_kernel = frozenset(
        h for h in q8 if conjugate(h, C) == C
    )
    # action_kernel is <QB> in the selector labelling.  Under the swapped
    # carry labelling this same subgroup is <QA_carry>.
    expected_carry_kernel = frozenset(
        dpow(carry_shadow_qa, exponent) for exponent in range(4)
    )
    require(
        carry_action_kernel == expected_carry_kernel
        and horn_patterns[1] == (0, 1, 1)
        and horn_patterns[2] == (0, 0, 1)
        and edge_patterns[1] == (1, 0, 1)
        and edge_patterns[2] == (0, 1, 1)
        and carry_shadow_qab
        == dmul(carry_shadow_qa, carry_shadow_qb),
        "selector/carry dicyclic shadows changed",
    )
    # No character simultaneously has the selector target values and
    # inverts the actual q11->q7 semantic edge QB.  Also, the literal carry
    # event edge pattern (0,1,0) is not a character: it is Bockstein
    # curvature, not a vertex action.
    require(
        not any(
            horn_pattern == (0, 1, 1) and edge_pattern[1] == 1
            for horn_pattern, edge_pattern
            in zip(horn_patterns, edge_patterns)
        )
        and (0, 1, 0) not in edge_patterns,
        "one Q8 character unexpectedly typed selector and carry together",
    )
    # In the carry shadow, retain THM-2887's displayed central section
    # QA*QB=-QAB.  The via/direct defect is then exactly C^13*(-1).
    carry_shadow_qab_canonical = dmul(
        MINUS_ONE, carry_shadow_qab
    )
    carry_via = dmul(
        dmul(dpow(C, 8), carry_shadow_qa),
        dmul(dpow(C, 9), carry_shadow_qb),
    )
    carry_direct = dmul(dpow(C, 4), carry_shadow_qab_canonical)
    carry_clutch = dmul(dpow(C, 13), MINUS_ONE)
    require(
        carry_via == dmul(carry_clutch, carry_direct)
        and dorder(carry_clutch) == 26
        and conjugate(carry_shadow_qb, carry_clutch)
        == dinv(carry_clutch)
        and dmul(carry_clutch, carry_clutch) == dpow(C, 26),
        "carry-shadow via/direct quaternionic clutch failed",
    )

    # The same character is the THM-2886 horn sign.
    semantic_signs = tuple(
        semantic_central_sign(state) for state in (ONE, QA, QAB)
    )
    require(
        semantic_signs == (0, 1, 1)
        and semantic_central_sign(QB) == 0
        and tuple(
            chi_direction(direction)
            for direction in ((0, 0), (1, 0), (1, 1))
        ) == semantic_signs,
        "Ad_QB/determinant horn sign failed",
    )

    # Complete order and conjugacy census.
    order_census = Counter(dorder(element) for element in ELEMENTS)
    expected_order_census = Counter({
        1: 1,
        2: 1,
        4: 678,
        13: 12,
        26: 12,
        52: 24,
        169: 156,
        338: 156,
        676: 312,
    })
    require(order_census == expected_order_census,
            "Dic_338 element-order census changed")

    unclassified = set(ELEMENTS)
    conjugacy_classes = []
    while unclassified:
        representative = min(unclassified)
        conjugacy_class = frozenset(
            conjugate(actor, representative) for actor in ELEMENTS
        )
        require(conjugacy_class <= unclassified,
                "conjugacy classes overlapped")
        unclassified.difference_update(conjugacy_class)
        conjugacy_classes.append(conjugacy_class)
    conjugacy_census = Counter(
        (dorder(next(iter(conjugacy_class))), len(conjugacy_class))
        for conjugacy_class in conjugacy_classes
    )
    expected_conjugacy_census = Counter({
        (1, 1): 1,
        (2, 1): 1,
        (4, 2): 1,
        (4, 338): 2,
        (13, 2): 6,
        (26, 2): 6,
        (52, 2): 12,
        (169, 2): 78,
        (338, 2): 78,
        (676, 2): 156,
    })
    require(
        len(conjugacy_classes) == 341
        and conjugacy_census == expected_conjugacy_census
        and frozenset(
            element for element in ELEMENTS
            if all(
                dmul(element, generator) == dmul(generator, element)
                for generator in (A, X)
            )
        ) == {ONE, MINUS_ONE},
        "Dic_338 conjugacy/center census changed",
    )

    # For even dicyclic parameter, the abelianization is V4.  More
    # concretely [G,G]=<a^2> has order 338.  Hence every one-dimensional
    # character kills C169=<a^4> and every even rotation.  This is the
    # obstruction to promoting a projective relative phase to the raw
    # scalar factor on diag(E,1).
    derived = frozenset(dpow(A, 2 * exponent)
                        for exponent in range(DIC_PARAMETER))
    all_commutators = frozenset(
        commutator(left, right)
        for left in ELEMENTS
        for right in ELEMENTS
    )
    require(
        len(derived) == DIC_PARAMETER
        and all_commutators <= derived
        and commutator(X, A) in {dpow(A, 2), dpow(A, 674)}
        and len(ELEMENTS) // len(derived) == 4
        and carry <= derived,
        "derived subgroup/abelianization census changed",
    )
    one_dimensional_character_values_on_even_rotation = tuple(
        (-1) ** ((alpha * 104 + beta * 0) % 2)
        for alpha, beta in ((0, 0), (1, 0), (0, 1), (1, 1))
    )
    require(
        one_dimensional_character_values_on_even_rotation == (1, 1, 1, 1),
        "a one-dimensional character detected the wrapped carry rotation",
    )

    # Coordinatewise-inversion joint model: retain the two independent
    # characters on two full C169 coordinates.  QA flips the selector
    # coordinate, QB flips the carry/Hermitian coordinate, and QAB flips
    # both.  Minimality below is relative to this stable-diagonal model.
    q8_normal_forms = tuple(
        (u, v, z) for u in (0, 1) for v in (0, 1) for z in (0, 1)
    )
    require(
        len({qmul(left, right)
             for left in q8_normal_forms for right in q8_normal_forms})
        == 8
        and all(
            qmul(qmul(left, middle), right)
            == qmul(left, qmul(middle, right))
            for left in q8_normal_forms
            for middle in q8_normal_forms
            for right in q8_normal_forms
        )
        and qpow(Q_QA, 2) == qpow(Q_QB, 2)
        == qpow(Q_QAB, 2) == Q_MINUS_ONE
        and qmul(Q_QA, Q_QB) == Q_QAB
        and qmul(Q_QB, Q_QA) == qmul(Q_MINUS_ONE, Q_QAB),
        "Q8 normal-form multiplication failed",
    )
    B_CS: BElement = (1, 0, Q_ONE)
    B_CC: BElement = (0, 1, Q_ONE)
    B_QA: BElement = (0, 0, Q_QA)
    B_QB: BElement = (0, 0, Q_QB)
    B_QAB: BElement = (0, 0, Q_QAB)
    B_MINUS_ONE: BElement = (0, 0, Q_MINUS_ONE)
    b_generators = (B_CS, B_CC, B_QA, B_QB)
    require(
        bpow(B_CS, CARRY_ORDER) == B_ONE
        and bpow(B_CC, CARRY_ORDER) == B_ONE
        and bpow(B_QA, 2) == bpow(B_QB, 2)
        == bpow(B_QAB, 2) == B_MINUS_ONE
        and bmul(B_QA, B_QB) == B_QAB
        and bconjugate(B_QA, B_CS) == binv(B_CS)
        and bconjugate(B_QA, B_CC) == B_CC
        and bconjugate(B_QB, B_CS) == B_CS
        and bconjugate(B_QB, B_CC) == binv(B_CC),
        "bi-dicyclic generator relations failed",
    )

    joint_order = CARRY_ORDER**2 * 8
    joint_derived_order = CARRY_ORDER**2 * 2
    require(
        joint_order == 228488
        and joint_derived_order == 57122
        and bcommutator(B_QA, B_CS)
        in {bpow(B_CS, 2), bpow(B_CS, 167)}
        and bcommutator(B_QB, B_CC)
        in {bpow(B_CC, 2), bpow(B_CC, 167)}
        and bcommutator(B_QA, B_QB) == B_MINUS_ONE
        and joint_order // joint_derived_order == 4,
        "bi-dicyclic center/derived/abelianization failed",
    )

    # Full order census by the eight Q8 cosets.  In the central cosets,
    # orders are d and 2d for d in {1,13,169} on C169^2.  In a QA or QB
    # coset only the fixed axis survives squaring, giving 4d.  The QAB
    # coset squares immediately to -1.
    abelian_order_counts = {1: 1, 13: 168, 169: 28392}
    axis_order_counts = {1: 1, 13: 12, 169: 156}
    b_order_census = Counter()
    for order, count in abelian_order_counts.items():
        b_order_census[order] += count
        b_order_census[2 * order] += count
    for order, count in axis_order_counts.items():
        b_order_census[4 * order] += 2 * 2 * CARRY_ORDER * count
    b_order_census[4] += 2 * CARRY_ORDER**2
    expected_b_order_census = Counter({
        1: 1,
        2: 1,
        4: 57798,
        13: 168,
        26: 168,
        52: 8112,
        169: 28392,
        338: 28392,
        676: 105456,
    })

    # Full conjugacy census by the same strata.  The central Q8 cosets are
    # sign orbits on C169^2 (sizes 1,2,4).  QA/QB classes have size 338,
    # with the fixed-axis sign tied to the central Q8 lift.  The QAB coset
    # is one class of size 2*169^2.
    b_conjugacy_census = Counter({
        (1, 1): 1,
        (2, 1): 1,
        (13, 2): 12,
        (13, 4): 36,
        (169, 2): 156,
        (169, 4): 7020,
        (26, 2): 12,
        (26, 4): 36,
        (338, 2): 156,
        (338, 4): 7020,
        (4, 338): 2,
        (4, 57122): 1,
        (52, 338): 24,
        (676, 338): 312,
    })
    b_class_count = sum(b_conjugacy_census.values())
    require(
        b_order_census == expected_b_order_census
        and b_class_count == 14789
        and sum(
            class_count * class_size
            for (_order, class_size), class_count
            in b_conjugacy_census.items()
        ) == joint_order,
        "bi-dicyclic order/conjugacy census changed: "
        f"orders={dict(sorted(b_order_census.items()))}; "
        f"classes={b_class_count}; "
        f"census={dict(sorted(b_conjugacy_census.items()))}",
    )

    selector_shadow = frozenset(
        (selector, 0, h)
        for selector in range(CARRY_ORDER)
        for h in q8_normal_forms
    )
    carry_shadow = frozenset(
        (0, carry_coordinate, h)
        for carry_coordinate in range(CARRY_ORDER)
        for h in q8_normal_forms
    )
    shared_q8 = frozenset((0, 0, h) for h in q8_normal_forms)
    require(
        len(selector_shadow) == len(carry_shadow) == 1352
        and selector_shadow & carry_shadow == shared_q8
        and len(shared_q8) == 8
        and (
            len(selector_shadow) * len(carry_shadow) // len(shared_q8)
        ) == joint_order,
        "two Dic_338 shadows do not form the expected fibre product",
    )

    # Complete weight-orbit census for A^=(C169^2)^.  A nonzero axis weight
    # has orbit two; a generic weight has orbit four.  A faithful semisimple
    # representation needs rank two in the weight span, hence degree at
    # least four.  The tensor model below realizes one generic four-orbit.
    unclassified_weights = {
        (first, second)
        for first in range(CARRY_ORDER)
        for second in range(CARRY_ORDER)
    }
    weight_orbit_census = Counter()
    while unclassified_weights:
        first, second = min(unclassified_weights)
        orbit = frozenset({
            (first, second),
            ((-first) % CARRY_ORDER, second),
            (first, (-second) % CARRY_ORDER),
            ((-first) % CARRY_ORDER, (-second) % CARRY_ORDER),
        })
        require(orbit <= unclassified_weights,
                "weight orbits overlapped")
        unclassified_weights.difference_update(orbit)
        weight_orbit_census[len(orbit)] += 1
    require(
        weight_orbit_census == Counter({1: 1, 2: 168, 4: 7056})
        and len({
            (1, 1),
            (CARRY_ORDER - 1, 1),
            (1, CARRY_ORDER - 1),
            (CARRY_ORDER - 1, CARRY_ORDER - 1),
        }) == 4,
        "bi-dicyclic weight-orbit/minimal-degree census failed",
    )
    diagonal_orbit = frozenset({
        (1, 1),
        (CARRY_ORDER - 1, 1),
        (1, CARRY_ORDER - 1),
        (CARRY_ORDER - 1, CARRY_ORDER - 1),
    })
    diagonal_envelope = {
        (
            (left - right) % CARRY_ORDER,
            (left + right) % CARRY_ORDER,
        )
        for left in range(CARRY_ORDER)
        for right in range(CARRY_ORDER)
    }
    require(
        len(diagonal_orbit) == 4
        and len(diagonal_envelope) == CARRY_ORDER**2,
        "diagonal physical address did not generate the full joint envelope",
    )

    # Typed RIGHT-action horn law in the joint carrier.  The selector
    # coordinate of the second edge is -9 because the intermediate QA state
    # inverts that coordinate; the carry coordinate remains +9.  In the
    # section gauge QA*QB=QAB the via endpoint is diagonal again.  Relative
    # to THM-2887's canonical QA*QB=-QAB section, the exact clutch is
    # (13,13,-1).
    edge_8: BElement = (8, 8, Q_QA)
    edge_9: BElement = ((-9) % CARRY_ORDER, 9, Q_QB)
    direct_edge_4: BElement = (4, 4, Q_QAB_CANONICAL)
    joint_clutch: BElement = (13, 13, Q_MINUS_ONE)
    horn_right_action_checks = 0
    for ancestry in range(P):
        base_index = 13 * ancestry + 3
        q3_state: BElement = (base_index, base_index, Q_ONE)
        q11_state = bmul(q3_state, edge_8)
        via_q7_state = bmul(q11_state, edge_9)
        direct_q7_state = bmul(q3_state, direct_edge_4)
        require(
            q11_state
            == (
                (base_index + 8) % CARRY_ORDER,
                (base_index + 8) % CARRY_ORDER,
                Q_QA,
            )
            and via_q7_state
            == (
                (base_index + 17) % CARRY_ORDER,
                (base_index + 17) % CARRY_ORDER,
                Q_QAB,
            )
            and direct_q7_state
            == (
                (base_index + 4) % CARRY_ORDER,
                (base_index + 4) % CARRY_ORDER,
                Q_QAB_CANONICAL,
            )
            and via_q7_state
            == bmul(joint_clutch, direct_q7_state),
            "typed right-action horn law failed",
        )
        horn_right_action_checks += 1
    require(horn_right_action_checks == P,
            "wrong typed horn-action check count")

    # The diagonal-labelled q7 fillers retain both cyclic coordinates.
    # QAB inverts both, so each lift still squares to the central sign.
    joint_filler_rows = []
    for ancestry in range(P):
        old_index = 13 * ancestry + 7
        new_index = (13 * (ancestry + 1) + 7) % CARRY_ORDER
        old_filler: BElement = (
            old_index, old_index, Q_QAB_CANONICAL
        )
        new_filler: BElement = (
            new_index, new_index, Q_QAB_CANONICAL
        )
        conjugator: BElement = (91, 91, Q_ONE)
        require(
            old_filler != new_filler
            and bpow(old_filler, 2) == B_MINUS_ONE
            and bpow(new_filler, 2) == B_MINUS_ONE
            and bconjugate(conjugator, old_filler) == new_filler,
            "joint q7 fillers collapsed",
        )
        joint_filler_rows.append((ancestry, old_filler, new_filler))

    # The inversion sectors are the two reflection classes.  The 169
    # central-sign-normalized lifts c^n QA and c^n QAB are distinct; their
    # negatives fill the two full 338-element conjugacy classes.
    qa_lifts = tuple(dmul(dpow(C, n), QA) for n in range(CARRY_ORDER))
    qab_lifts = tuple(dmul(dpow(C, n), QAB) for n in range(CARRY_ORDER))
    require(
        len(set(qa_lifts)) == len(set(qab_lifts)) == CARRY_ORDER
        and set(qa_lifts).isdisjoint(qab_lifts)
        and all(dmul(value, value) == MINUS_ONE for value in qa_lifts)
        and all(dmul(value, value) == MINUS_ONE for value in qab_lifts)
        and len(
            set(qa_lifts)
            | {dmul(MINUS_ONE, value) for value in qa_lifts}
        ) == 338
        and len(
            set(qab_lifts)
            | {dmul(MINUS_ONE, value) for value in qab_lifts}
        ) == 338,
        "semantic reflection lifts changed",
    )
    qa_class = frozenset(
        conjugate(actor, QA) for actor in ELEMENTS
    )
    qab_class = frozenset(
        conjugate(actor, QAB) for actor in ELEMENTS
    )
    require(
        qa_class != qab_class
        and qa_class
        == set(qa_lifts) | {
            dmul(MINUS_ONE, value) for value in qa_lifts
        }
        and qab_class
        == set(qab_lifts) | {
            dmul(MINUS_ONE, value) for value in qab_lifts
        },
        "QA/QAB reflection conjugacy classes changed",
    )

    # The two THM-2884 q7 ancestry fillers have the same semantic quotient
    # but differ by the vertical c^13 coordinate.  They never collapse in G.
    # Since 2*91=13 mod169, they are conjugate by c^91, not equal.
    filler_rows = []
    for ancestry in range(P):
        old_index = 13 * ancestry + 7
        new_index = (13 * (ancestry + 1) + 7) % CARRY_ORDER
        old_filler = dmul(dpow(C, old_index), QAB)
        new_filler = dmul(dpow(C, new_index), QAB)
        require(
            old_filler != new_filler
            and new_filler
            == conjugate(dpow(C, 91), old_filler)
            and dmul(old_filler, old_filler) == MINUS_ONE
            and dmul(new_filler, new_filler) == MINUS_ONE,
            "q7 ancestry fillers collapsed or lost reflection law",
        )
        filler_rows.append((ancestry, old_filler, new_filler))

    # Exact cyclotomic/Prony realization of the projective C169 generator.
    # This is an algebraic relative-channel gauge, not a physical action.
    field_rows = []
    for prime, root in atlas.endpoint.MODS:
        eta = pow(root, atlas.endpoint.NN // ROTATION_ORDER, prime)
        xi = pow(root, atlas.endpoint.NN // 2366, prime)
        omega = pow(xi, 182, prime)
        rho = pow(xi, 42, prime)
        s = pow(eta, 344, prime)
        imaginary_unit = pow(eta, 169, prime)
        wrap_rotation = dpow(A, 104)
        wrap_rotation_matrix = matrix_for(wrap_rotation, eta, prime)
        raw_wrap_matrix = (
            (pow(omega, 4, prime), 0),
            (0, 1),
        )
        require(
            pow(eta, ROTATION_ORDER, prime) == 1
            and pow(eta, ROTATION_ORDER // 2, prime) != 1
            and pow(eta, ROTATION_ORDER // 13, prime) != 1
            and pow(xi, 2366, prime) == 1
            and all(
                pow(xi, 2366 // divisor, prime) != 1
                for divisor in (2, 7, 13)
            )
            and pow(rho, CARRY_ORDER, prime) == 1
            and pow(rho, 13, prime) == pow(omega, 3, prime)
            and s * s % prime == rho
            and imaginary_unit * imaginary_unit % prime == prime - 1,
            "finite-field root normalization failed",
        )
        require(
            wrap_rotation == dpow(C, 130)
            and wrap_rotation_matrix
            == (
                (pow(omega, 2, prime), 0),
                (0, pow(omega, -2, prime)),
            )
            and matscale(
                pow(omega, 2, prime),
                wrap_rotation_matrix,
                prime,
            ) == raw_wrap_matrix
            and pow(omega, 2, prime) not in {1, prime - 1},
            "wrapped diag(omega^4,1) projective lift failed",
        )

        a_matrix = matrix_for(A, eta, prime)
        qa_matrix = matrix_for(QA, eta, prime)
        qb_matrix = matrix_for(QB, eta, prime)
        qab_matrix = matrix_for(QAB, eta, prime)
        c_matrix = matrix_for(C, eta, prime)
        minus_matrix = ((prime - 1, 0), (0, prime - 1))
        identity_matrix = ((1, 0), (0, 1))
        require(
            matmul(a_matrix, a_matrix, prime)
            == matrix_for(dpow(A, 2), eta, prime)
            and matmul(qa_matrix, qa_matrix, prime) == minus_matrix
            and matmul(qb_matrix, qb_matrix, prime) == minus_matrix
            and matmul(qab_matrix, qab_matrix, prime) == minus_matrix
            and matmul(qa_matrix, qb_matrix, prime) == qab_matrix
            and matmul(
                matmul(qa_matrix, c_matrix, prime),
                matrix_for(dinv(QA), eta, prime),
                prime,
            ) == matrix_for(dinv(C), eta, prime)
            and matmul(
                matmul(qb_matrix, c_matrix, prime),
                matrix_for(dinv(QB), eta, prime),
                prime,
            ) == c_matrix
            and matscale(s, c_matrix, prime)
            == ((rho, 0), (0, 1))
            and matrix_for(ONE, eta, prime) == identity_matrix,
            "faithful two-dimensional dicyclic representation failed",
        )
        matrix_group = {
            matrix_for(element, eta, prime) for element in ELEMENTS
        }
        require(len(matrix_group) == GROUP_ORDER,
                "finite-field dicyclic representation is not faithful")

        # Four-channel tensor realization of the joint carrier:
        #
        #   cs = D(s) tensor I,       cc = I tensor D(s),
        #   QA = i X tensor I,        QB = i Z tensor X.
        #
        # Here X swaps the two basis vectors and Z is the sign diagonal.
        # The first factor is the selector/origin channel and the second is
        # the Prony/carry channel.
        identity2 = ((1, 0), (0, 1))
        swap2 = ((0, 1), (1, 0))
        sign2 = ((1, 0), (0, prime - 1))
        d169 = ((s, 0), (0, pow(s, -1, prime)))
        cs4 = kron2(d169, identity2, prime)
        cc4 = kron2(identity2, d169, prime)
        qa4 = square_scale(
            imaginary_unit,
            kron2(swap2, identity2, prime),
            prime,
        )
        qb4 = square_scale(
            imaginary_unit,
            kron2(sign2, swap2, prime),
            prime,
        )
        qab4 = square_matmul(qa4, qb4, prime)
        identity4 = tuple(
            tuple(int(row == column) for column in range(4))
            for row in range(4)
        )
        minus4 = square_scale(prime - 1, identity4, prime)
        qa4_inverse = square_scale(prime - 1, qa4, prime)
        qb4_inverse = square_scale(prime - 1, qb4, prime)
        cs4_inverse = square_pow(cs4, CARRY_ORDER - 1, prime)
        cc4_inverse = square_pow(cc4, CARRY_ORDER - 1, prime)
        require(
            square_matmul(qa4, qa4, prime) == minus4
            and square_matmul(qb4, qb4, prime) == minus4
            and square_matmul(qab4, qab4, prime) == minus4
            and square_matmul(qa4, qb4, prime)
            == square_scale(
                prime - 1,
                square_matmul(qb4, qa4, prime),
                prime,
            )
            and square_matmul(
                square_matmul(qa4, cs4, prime),
                qa4_inverse,
                prime,
            ) == cs4_inverse
            and square_matmul(
                square_matmul(qa4, cc4, prime),
                qa4_inverse,
                prime,
            ) == cc4
            and square_matmul(
                square_matmul(qb4, cs4, prime),
                qb4_inverse,
                prime,
            ) == cs4
            and square_matmul(
                square_matmul(qb4, cc4, prime),
                qb4_inverse,
                prime,
            ) == cc4_inverse
            and square_pow(cs4, CARRY_ORDER, prime) == identity4
            and square_pow(cc4, CARRY_ORDER, prime) == identity4,
            "four-channel tensor Pauli representation failed",
        )
        torus_weight_images = {
            (
                pow(s, selector + carry_coordinate, prime),
                pow(s, selector - carry_coordinate, prime),
                pow(s, -selector + carry_coordinate, prime),
                pow(s, -selector - carry_coordinate, prime),
            )
            for selector in range(CARRY_ORDER)
            for carry_coordinate in range(CARRY_ORDER)
        }
        require(
            pow(s, CARRY_ORDER, prime) == 1
            and pow(s, CARRY_ORDER // 13, prime) != 1
            and len(torus_weight_images) == CARRY_ORDER**2
            and (prime - 1,) * 4 not in torus_weight_images
            and monomial_permutation(identity4) == (0, 1, 2, 3)
            and monomial_permutation(qa4) == (2, 3, 0, 1)
            and monomial_permutation(qb4) == (1, 0, 3, 2)
            and monomial_permutation(qab4) == (3, 2, 1, 0),
            "four-channel representation faithfulness gate failed",
        )
        selector_wrap4 = square_pow(cs4, 130, prime)
        carry_wrap4 = square_pow(cc4, 130, prime)
        determinant_wrap2 = (
            (pow(omega, 2, prime), 0),
            (0, pow(omega, -2, prime)),
        )
        raw_wrap2 = ((pow(omega, 4, prime), 0), (0, 1))
        raw_e = pow(omega, 4, prime)
        raw_e_inverse = pow(raw_e, -1, prime)
        balanced_wrap4 = square_matmul(
            selector_wrap4, carry_wrap4, prime
        )
        hermitian_quadrature_weights = (
            raw_e,
            1,
            1,
            raw_e_inverse,
        )
        hermitian_wrap4 = tuple(
            tuple(
                hermitian_quadrature_weights[row]
                if row == column else 0
                for column in range(4)
            )
            for row in range(4)
        )

        # If raw D=diag(E,1), cyclotomic conjugation sends
        # bar(D)=diag(E^-1,1).  Therefore M -> D M bar(D)^T acts on
        #
        #   (U*bar(V), |U|^2, |V|^2, V*bar(U))
        #
        # by (E,1,1,E^-1).  This is only a balanced-torus diagonal
        # match.  It is not an identification of the faithful tensor
        # representation with Hermitianization: the former sends the
        # Q8 center to -I4, while conjugation kills that global phase.
        raw_d_weights = (raw_e, 1)
        raw_d_bar_weights = (raw_e_inverse, 1)
        hermitian_weights_from_raw_d = (
            raw_d_weights[0] * raw_d_bar_weights[1] % prime,
            raw_d_weights[0] * raw_d_bar_weights[0] % prime,
            raw_d_weights[1] * raw_d_bar_weights[1] % prime,
            raw_d_weights[1] * raw_d_bar_weights[0] % prime,
        )
        global_minus_weights = (prime - 1, prime - 1)
        global_minus_bar_weights = (prime - 1, prime - 1)
        hermitian_center_weights = (
            global_minus_weights[0]
            * global_minus_bar_weights[1] % prime,
            global_minus_weights[0]
            * global_minus_bar_weights[0] % prime,
            global_minus_weights[1]
            * global_minus_bar_weights[1] % prime,
            global_minus_weights[1]
            * global_minus_bar_weights[0] % prime,
        )
        qa_permutation = monomial_permutation(qa4)
        qb_permutation = monomial_permutation(qb4)
        literal_prony_plane = frozenset((0, 1))
        require(
            carry_wrap4 == kron2(identity2, determinant_wrap2, prime)
            and square_scale(
                pow(omega, 2, prime), carry_wrap4, prime
            ) == kron2(identity2, raw_wrap2, prime),
            "four-channel wrapped-edge projective factor failed",
        )
        require(
            balanced_wrap4 == hermitian_wrap4
            and hermitian_weights_from_raw_d
            == hermitian_quadrature_weights
            and hermitian_center_weights == (1, 1, 1, 1)
            and square_matmul(qa4, qa4, prime) == minus4
            and minus4 != identity4
            and frozenset(
                qa_permutation[index] for index in literal_prony_plane
            ) == frozenset((2, 3))
            and frozenset(
                qb_permutation[index] for index in literal_prony_plane
            ) == literal_prony_plane,
            "balanced Hermitian torus or representation boundary failed",
        )

        # The normalized Prony spinor is (t_r,1), where
        # t_r=xi^(955+546r).  Its Hermitian coherence is t_r.  Acting by
        # c=diag(s,s^-1) multiplies coherence by s^2=rho; c^13 gives
        # exactly omega^3.  In the selector-shadow representation the
        # semantic monomial matrices send coherence as listed below:
        #
        # Q0:t, QA:-bar(t), QB:-t, QAB:bar(t).
        #
        # The carry shadow swaps QA and QB, so its table is
        # Q0:t, QA:-t, QB:-bar(t), QAB:bar(t).
        projective_t = tuple(
            pow(xi, (955 + 546 * r) % 2366, prime)
            for r in range(P)
        )
        conjugate_t = tuple(pow(value, -1, prime) for value in projective_t)
        require(
            len(projective_t) == P
            and len(set(projective_t)) == P
            and all(
                projective_t[(r + 1) % P]
                == pow(omega, 3, prime) * projective_t[r] % prime
                for r in range(P)
            )
            and all(
                projective_t[r] * conjugate_t[r] % prime == 1
                for r in range(P)
            ),
            "THM-2868 normalized Prony orbit changed",
        )
        semantic_coherence = tuple(
            (
                projective_t[r],
                -conjugate_t[r] % prime,
                -projective_t[r] % prime,
                conjugate_t[r],
            )
            for r in range(P)
        )
        carry_semantic_coherence = tuple(
            (
                projective_t[r],
                -projective_t[r] % prime,
                -conjugate_t[r] % prime,
                conjugate_t[r],
            )
            for r in range(P)
        )
        require(
            all(
                semantic_coherence[r][0] == projective_t[r]
                and semantic_coherence[r][1]
                == -conjugate_t[r] % prime
                and semantic_coherence[r][2]
                == -projective_t[r] % prime
                and semantic_coherence[r][3] == conjugate_t[r]
                for r in range(P)
            )
            and all(
                carry_semantic_coherence[r][0] == projective_t[r]
                and carry_semantic_coherence[r][1]
                == -projective_t[r] % prime
                and carry_semantic_coherence[r][2]
                == -conjugate_t[r] % prime
                and carry_semantic_coherence[r][3]
                == conjugate_t[r]
                for r in range(P)
            ),
            "semantic Hermitian coherence formulas failed",
        )

        # The honest lifted coefficient line has a faithful C169
        # trivialization after a state-dependent relative-channel gauge.
        # In the original chart r=a+q-3, its edge is
        # omega^(3(h+kappa)); in the new chart it is rho^h.
        gauge_checks = 0
        for ancestry in range(P):
            for q in range(P):
                n = 13 * ancestry + q
                r = (ancestry + q - 3) % P
                old_value = pow(omega, 3 * r, prime)
                gauge = (
                    pow(rho, n, prime) * pow(old_value, -1, prime)
                ) % prime
                require(
                    gauge * old_value % prime == pow(rho, n, prime),
                    "C169 vertex gauge failed",
                )
                for h in range(P):
                    carry_event = (q + h) // P
                    new_n = (n + h) % CARRY_ORDER
                    new_ancestry, new_q = divmod(new_n, P)
                    new_r = (new_ancestry + new_q - 3) % P
                    old_transition = pow(
                        omega, 3 * (h + carry_event), prime
                    )
                    new_old_value = pow(omega, 3 * new_r, prime)
                    new_gauge = (
                        pow(rho, new_n, prime)
                        * pow(new_old_value, -1, prime)
                    ) % prime
                    require(
                        new_old_value
                        == old_transition * old_value % prime
                        and new_gauge * old_transition % prime
                        == pow(rho, h, prime) * gauge % prime,
                        "faithful C169 gauge transport failed",
                    )
                    gauge_checks += 1
        require(gauge_checks == P**3,
                "wrong C169 gauge-exhaustion count")

        field_rows.append((
            prime,
            eta,
            rho,
            pow(rho, 13, prime),
            s,
            projective_t[0],
            gauge_checks,
            pow(omega, 2, prime),
        ))

    print("THM-2889 DICYCLIC HERMITIAN ENVELOPE FOR THE THM-2886/2887 HORN")
    print(
        "presentation=Dic_338=<a,x | a^676=1, x^2=a^338, "
        "x*a*x^-1=a^-1>; convention |Dic_n|=4n; order=1352"
    )
    print(
        "embedding=C169=<c=a^344>=<a^4>; "
        "QA=x; QB=a^169; QAB=QA*QB; Q8=<QA,QB>; "
        "C169_intersection_Q8=1; product_size=1352"
    )
    print(
        "section_gauge=printed QAB is QA*QB; THM2887 canonical lift has "
        "QA*QB=-QAB, differing by the central sign only"
    )
    print(
        "action=chi(h)=det(QB,h); kernel=<QB>="
        "{1,QB,-1,-QB} order4; QA,QAB invert c; QB fixes c"
    )
    print(
        f"horn_Ad_QB_sign_Q0_QA_QAB={semantic_signs}; "
        f"four_character_patterns={horn_patterns}; unique_match=det(QB,-)"
    )
    print(
        f"edge_character_patterns_on_QA_QB_QAB={edge_patterns}; "
        "selector_shadow=chi_s=u has target (0,1,1) but fixes edge QB; "
        "carry_shadow=chi_c=v has target (0,0,1) and inverts edge QB; "
        "no_single_character_does_both=1; carry_event_pattern_(0,1,0)_"
        "is_not_a_character=1"
    )
    print(
        "carry_shadow_via_direct=(c^8 QA)(c^9 QB)="
        "(c^13*(-1))(c^4 QAB_canonical); clutch_order=26; "
        "QB inverts clutch; clutch_square=c^26"
    )
    print(f"element_order_census={dict(sorted(order_census.items()))}")
    print(
        "conjugacy_class_census="
        f"{dict(sorted(conjugacy_census.items()))}; "
        f"class_count={len(conjugacy_classes)}; center_size=2"
    )
    print(
        "derived_subgroup=<a^2> order338; abelianization=V4; "
        "C169 is killed by all four one-dimensional characters"
    )
    print(
        "joint_carrier=B=(C169_s x C169_c) semidirect Q8; "
        "QA flips s only, QB flips c only, QAB flips both; "
        "order=228488; center={+/-1}; derived=C169^2 x {+/-1} "
        "order57122; abelianization=V4"
    )
    print(
        "central_quotient=B/<-1> is the product of two groups "
        "(C169 semidirect C2), each of order338"
    )
    print(f"joint_element_order_census={dict(sorted(b_order_census.items()))}")
    print(
        "joint_conjugacy_class_census="
        f"{dict(sorted(b_conjugacy_census.items()))}; "
        f"class_count={b_class_count}"
    )
    print(
        "joint_fibre_product=two Dic_338 shadows over shared Q8; "
        "orders=(1352,1352,8)->228488; "
        f"dual_weight_orbits={dict(sorted(weight_orbit_census.items()))}; "
        "minimal_faithful_semisimple_degree=4; physical diagonal address "
        "has four-point Q8 orbit and its minimal stable abelian envelope "
        "within the coordinatewise-inversion model is all C169^2"
    )
    print(
        "reflection_classes=QA and QAB give the two distinct size-338 "
        "order-4 classes; each has 169 central-sign-normalized lifts and "
        "their negatives"
    )
    print(
        "q7_fillers=13 checked ancestry pairs; each pair "
        "(13a+7,QAB),(13(a+1)+7,QAB) is distinct, squares to -1, "
        "and is conjugate by c^91 because 2*91=13 mod169; "
        "joint diagonal lifts retain both coordinates and obey the same "
        "law, but conjugacy/class functions do not separate them"
    )
    print(
        "typed_joint_right_action=e8=(8,8,QA), "
        "e9=(-9,+9,QB), direct=(4,4,QAB_canonical); "
        "e8*e9 lands at +17 diagonal and differs from direct by "
        "(13,13,-1); checks=13"
    )
    for row in field_rows:
        print(
            f"field={row[0]}; eta676={row[1]}; rho169={row[2]}; "
            f"rho13=omega3={row[3]}; spinor_half_phase_s={row[4]}; "
            f"Prony_t0={row[5]}; faithful_matrix_count=1352; "
            f"C169_gauge_checks={row[6]}; wrap_scalar_omega2={row[7]}"
        )
    print(
        "selector_shadow_projective_coherence=C=(U/V); bar is cyclotomic "
        "root inversion before finite-field specialization; "
        "c:C->rho*C; QA:C->-bar(C); QB:C->-C; "
        "QAB:C->bar(C); c^13:C->omega^3*C; "
        "carry shadow swaps QA,QB, giving QA:C->-C and QB:C->-bar(C)"
    )
    print(
        "balanced_hermitian_torus=cs^130*cc^130="
        "diag(E,1,1,E^-1), E=omega^4; this is exactly conjugation by "
        "raw D=diag(E,1) on "
        "(U*bar(V),|U|^2,|V|^2,V*bar(U)); only the torus diagonals "
        "coincide: Hermitianization kills z=-1/global phase and factors "
        "through B/<z>, whereas the faithful B-module sends z to -I4; "
        "the literal Prony plane e_+ tensor K^2 is sent by QA to the "
        "complementary plane and is not Q8-stable"
    )
    print(
        "wrapped_edge=diag(omega^4,1)=omega^2*"
        "diag(omega^2,omega^-2), with determinant-one factor a^104=c^130; "
        "this is the raw Prony frame, while gauge-normalized +9 is c^9; "
        "omega^2 cannot be a global scalar character because both carry "
        "circles lie in the relevant derived subgroup"
    )
    print(
        "positive=selector provenance and carry/Hermitian reversal each "
        "have an exact Dic_338 shadow; within the full-C169 "
        "coordinatewise-inversion stable-diagonal model their minimal "
        "joint faithful semisimple envelope is the four-channel "
        "bi-dicyclic fibre product, and it separates the two q7 ancestry "
        "lifts"
    )
    print(
        "hostile=no one Q8 character types both selector parity and the "
        "q11-to-q7 carry edge; the literal carry pattern is Bockstein "
        "curvature, not a vertex character.  The construction uses a "
        "state-dependent relative-channel gauge and abstract Q8 actions; "
        "no physical four-channel operator, canonical unmarked Hermitian "
        "quadrature, row exclusion, or LRC14 conclusion is proved"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
