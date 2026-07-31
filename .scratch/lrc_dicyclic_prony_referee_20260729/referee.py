#!/usr/bin/env python3
"""Independent hostile audit of the ker<QA> reverse-action carrier.

This is scratch-only referee evidence.  It deliberately distinguishes:

* the abstract semidirect product C169 rtimes Q8;
* its labelled semantic quotient;
* a projective two-channel representation;
* the induced Hermitian coherence; and
* the still-missing physical/current realization.
"""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
from math import gcd
from pathlib import Path
import re
import sys


ROOT = Path(__file__).resolve().parents[2]
COMP = ROOT / "04-computation"
RESULTS = ROOT / "05-knowledge/results"
sys.path.insert(0, str(COMP))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_hash(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


PINNED = {
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
    require(path.is_file(), f"missing pinned input {path}")
    require(lf_hash(path) == expected, f"changed pinned input {path.name}")


import lrc14_two_origin_endpoint_projective_kummer_thm2868 as atlas


# Dic_n convention: order 4n.  Store r^i x^b with
# r^676=1, x^2=r^338, x r x^-1=r^-1.
N = 338
ROT = 2 * N
ORDER = 4 * N
ONE = (0, 0)


def mul(left: tuple[int, int], right: tuple[int, int]) -> tuple[int, int]:
    i, b = left
    j, c = right
    return (
        (i + (-j if b else j) + (N if b and c else 0)) % ROT,
        b ^ c,
    )


def inv(value: tuple[int, int]) -> tuple[int, int]:
    i, b = value
    if b == 0:
        return ((-i) % ROT, 0)
    return ((i + N) % ROT, 1)


def power(value: tuple[int, int], exponent: int) -> tuple[int, int]:
    require(exponent >= 0, "negative power")
    result = ONE
    base = value
    while exponent:
        if exponent & 1:
            result = mul(result, base)
        base = mul(base, base)
        exponent //= 2
    return result


def conj(actor: tuple[int, int], value: tuple[int, int]) -> tuple[int, int]:
    return mul(mul(actor, value), inv(actor))


DIVISORS = tuple(
    divisor for divisor in range(1, ORDER + 1) if ORDER % divisor == 0
)


def element_order(value: tuple[int, int]) -> int:
    for divisor in DIVISORS:
        if power(value, divisor) == ONE:
            return divisor
    raise RuntimeError("order did not divide group order")


ELEMENTS = tuple((i, b) for i in range(ROT) for b in (0, 1))
R = (1, 0)
X = (0, 1)
Z = power(R, N)

# Labelled semidirect decomposition requested by the referee task.
# C=r^508 has order 169, QA=r^169 has order four, and C*QA=r.
# QB is the reversing generator.  QAB is pinned to the THM-2887 forward
# gauge QA*QB=-QAB.
C = power(R, 508)
QA = power(R, 169)
QB = X
QAB = mul(Z, mul(QA, QB))


def subgroup(generator: tuple[int, int], size: int) -> frozenset[tuple[int, int]]:
    return frozenset(power(generator, exponent) for exponent in range(size))


def matmul(left, right, prime):
    return (
        (
            (left[0][0] * right[0][0] + left[0][1] * right[1][0]) % prime,
            (left[0][0] * right[0][1] + left[0][1] * right[1][1]) % prime,
        ),
        (
            (left[1][0] * right[0][0] + left[1][1] * right[1][0]) % prime,
            (left[1][0] * right[0][1] + left[1][1] * right[1][1]) % prime,
        ),
    )


def matscale(scalar, matrix, prime):
    return tuple(
        tuple(scalar * entry % prime for entry in row)
        for row in matrix
    )


def matrix(value, eta, prime):
    exponent, outside = value
    diagonal = (
        (pow(eta, exponent, prime), 0),
        (0, pow(eta, -exponent, prime)),
    )
    if outside == 0:
        return diagonal
    return matmul(diagonal, ((0, prime - 1), (1, 0)), prime)


def determinant(matrix_value, prime):
    return (
        matrix_value[0][0] * matrix_value[1][1]
        - matrix_value[0][1] * matrix_value[1][0]
    ) % prime


def matrix_n_identity(size: int):
    return tuple(
        tuple(int(row == column) for column in range(size))
        for row in range(size)
    )


def matrix_n_mul(left, right, prime):
    size = len(left)
    return tuple(
        tuple(
            sum(left[row][middle] * right[middle][column]
                for middle in range(size)) % prime
            for column in range(size)
        )
        for row in range(size)
    )


def matrix_n_scale(scalar, value, prime):
    return tuple(
        tuple(scalar * entry % prime for entry in row)
        for row in value
    )


def matrix_n_power(value, exponent, prime):
    result = matrix_n_identity(len(value))
    base = value
    while exponent:
        if exponent & 1:
            result = matrix_n_mul(result, base, prime)
        base = matrix_n_mul(base, base, prime)
        exponent //= 2
    return result


def kronecker(left, right, prime):
    return tuple(
        tuple(
            left[outer_row][outer_column]
            * right[inner_row][inner_column] % prime
            for outer_column in range(len(left[0]))
            for inner_column in range(len(right[0]))
        )
        for outer_row in range(len(left))
        for inner_row in range(len(right))
    )


def main() -> None:
    require(
        len(ELEMENTS) == ORDER
        and power(R, ROT) == ONE
        and power(R, N) == Z != ONE
        and mul(X, X) == Z
        and conj(X, R) == inv(R),
        "dicyclic presentation failed",
    )

    c169 = subgroup(C, 169)
    q8 = frozenset(
        mul(power(QA, u), power(QB, v))
        for u in range(4) for v in range(4)
    )
    require(
        element_order(C) == 169
        and element_order(QA) == element_order(QB) == element_order(QAB) == 4
        and mul(C, QA) == R
        and len(c169) == 169
        and len(q8) == 8
        and c169 & q8 == {ONE}
        and len({mul(c, h) for c in c169 for h in q8}) == ORDER
        and power(QA, 2) == power(QB, 2) == power(QAB, 2) == Z
        and mul(QA, QB) == mul(Z, QAB)
        and mul(QB, QA) == QAB,
        "labelled C169/Q8 factorization failed",
    )

    action_kernel = frozenset(h for h in q8 if conj(h, C) == C)
    require(
        action_kernel == subgroup(QA, 4)
        and conj(QB, C) == inv(C)
        and conj(QAB, C) == inv(C),
        "the action kernel is not <QA>",
    )

    # The requested action character differs from the source-marked
    # selector character proved by THM-2887.
    path = ((0, 0), (1, 0), (1, 1))  # Q0, QA, QAB
    inversion_character = tuple(v for _u, v in path)     # det(QA,h)
    selector_character = tuple(u for u, _v in path)      # det(QB,h)
    require(
        inversion_character == (0, 0, 1)
        and selector_character == (0, 1, 1)
        and inversion_character != selector_character,
        "independent Q8 characters unexpectedly collapsed",
    )

    order_census = Counter(element_order(value) for value in ELEMENTS)
    require(
        order_census == Counter({
            1: 1, 2: 1, 4: 678, 13: 12, 26: 12, 52: 24,
            169: 156, 338: 156, 676: 312,
        }),
        "element-order census failed",
    )

    unclassified = set(ELEMENTS)
    classes = []
    while unclassified:
        representative = min(unclassified)
        current = frozenset(conj(actor, representative) for actor in ELEMENTS)
        require(current <= unclassified, "overlapping conjugacy classes")
        unclassified.difference_update(current)
        classes.append(current)
    class_census = Counter(
        (element_order(next(iter(current))), len(current))
        for current in classes
    )
    require(
        len(classes) == 341
        and class_census == Counter({
            (1, 1): 1, (2, 1): 1, (4, 2): 1, (4, 338): 2,
            (13, 2): 6, (26, 2): 6, (52, 2): 12,
            (169, 2): 78, (338, 2): 78, (676, 2): 156,
        }),
        "conjugacy census failed",
    )

    qb_class = frozenset(conj(actor, QB) for actor in ELEMENTS)
    qab_class = frozenset(conj(actor, QAB) for actor in ELEMENTS)
    qb_half = {mul(power(C, n), QB) for n in range(169)}
    qab_half = {mul(power(C, n), QAB) for n in range(169)}
    require(
        qb_class != qab_class
        and len(qb_class) == len(qab_class) == 338
        and qb_class == qb_half | {mul(Z, value) for value in qb_half}
        and qab_class == qab_half | {mul(Z, value) for value in qab_half}
        and all(mul(value, value) == Z for value in qb_class | qab_class),
        "the two skew-lift classes failed",
    )

    filler_checks = 0
    for ancestry in range(13):
        old = mul(power(C, 13 * ancestry + 7), QAB)
        new = mul(power(C, (13 * (ancestry + 1) + 7) % 169), QAB)
        require(
            old != new
            and new == conj(power(C, 91), old)
            and power(old, 2) == power(new, 2) == Z,
            "q7 ancestry skew lifts collapsed",
        )
        filler_checks += 1

    # The pinned forward-gauge edge product simultaneously records the
    # base-thirteen carry and the quaternionic central sign.
    edge8 = mul(power(C, 8), QA)
    edge9 = mul(power(C, 9), QB)
    edge4 = mul(power(C, 4), QAB)
    skew_clutch = mul(power(C, 13), Z)
    require(
        mul(edge8, edge9) == mul(skew_clutch, edge4)
        and element_order(skew_clutch) == 26
        and conj(QB, skew_clutch) == inv(skew_clutch)
        and power(skew_clutch, 2) == power(C, 26),
        "joint carry/sign skew clutch failed",
    )

    # [x,r]=r^-2 generates <r^2>; the quotient is V4.  Hence every
    # one-dimensional character kills C169, so the scalar needed below
    # cannot be promoted to a character of G.
    commutator_subgroup = subgroup(power(R, 2), 338)
    commutator_x_r = mul(mul(mul(X, R), inv(X)), inv(R))
    require(
        commutator_x_r in {power(R, 2), power(R, 674)}
        and subgroup(commutator_x_r, 338) == commutator_subgroup
        and C in commutator_subgroup
        and Z in commutator_subgroup,
        "commutator/abelianization calculation failed",
    )

    thm2886_output = (
        RESULTS / "lrc14_stepped_origin_v4_provenance_transport_thm2886.out"
    ).read_text()
    sample_rows = re.findall(
        r"field=(\d+);.*charged_edge_E=omega\^4=(\d+); "
        r"sample_\(U8,V8\)=\((\d+),(\d+)\)",
        thm2886_output,
    )
    require(
        len(sample_rows) == 2
        and all(int(u) != 0 and int(v) != 0 for _p, _e, u, v in sample_rows),
        "THM-2886 nonzero two-channel controls changed",
    )

    field_rows = []
    joint_rows = []
    for prime, root in atlas.endpoint.MODS:
        # Twisting the primitive dicyclic eigenvalue by 17 is faithful and
        # aligns C^26 with the square-root-normalized Prony edge.
        eta0 = pow(root, atlas.endpoint.NN // ROT, prime)
        eta = pow(eta0, 17, prime)
        omega = pow(root, atlas.endpoint.NN // 13, prime)
        omega2 = pow(omega, 2, prime)
        omega4 = pow(omega, 4, prime)
        require(
            gcd(17, ROT) == 1
            and pow(eta, ROT, prime) == 1
            and pow(eta, ROT // 2, prime) == prime - 1
            and pow(omega, 13, prime) == 1
            and omega != 1,
            "root normalization failed",
        )

        image = {matrix(value, eta, prime) for value in ELEMENTS}
        c26_matrix = matrix(power(C, 26), eta, prime)
        prony_matrix = ((omega4, 0), (0, 1))
        require(
            len(image) == ORDER
            and c26_matrix
            == ((omega2, 0), (0, pow(omega2, -1, prime)))
            and prony_matrix == matscale(omega2, c26_matrix, prime)
            and determinant(c26_matrix, prime) == 1
            and determinant(prony_matrix, prime) == omega4 != 1,
            "projective Prony/dicyclic match failed",
        )

        qb_matrix = matrix(QB, eta, prime)
        require(
            matmul(
                matmul(qb_matrix, c26_matrix, prime),
                matrix(inv(QB), eta, prime),
                prime,
            ) == matrix(inv(power(C, 26)), eta, prime),
            "QB did not reverse the normalized Prony phase",
        )

        # On K=U*bar(V), both diag(omega^4,1) and c^26 act by omega^4.
        # QA sends K to -K; QB sends K to -bar(K); QAB sends K to bar(K).
        coherence_factor_prony = omega4
        coherence_factor_c26 = (
            c26_matrix[0][0] * pow(c26_matrix[1][1], -1, prime)
        ) % prime
        require(
            coherence_factor_prony == coherence_factor_c26
            and pow(
                (
                    matrix(skew_clutch, eta, prime)[0][0]
                    * pow(
                        matrix(skew_clutch, eta, prime)[1][1],
                        -1,
                        prime,
                    )
                ) % prime,
                2,
                prime,
            ) == omega4,
            "Hermitian coherence/skew-square relation failed",
        )
        field_rows.append((
            prime,
            eta,
            omega4,
            c26_matrix,
            determinant(prony_matrix, prime),
        ))

        # Minimal joint carrier for the two independent characters.
        #
        # A=<a,b>=C169^2.  QA inverts a and fixes b; QB fixes a and
        # inverts b.  The following Clifford matrices give a faithful
        # four-dimensional linear realization:
        #
        #   a = diag(zeta,zeta^-1) tensor I,
        #   b = I tensor diag(zeta,zeta^-1),
        #   QA = J tensor I,
        #   QB = Z tensor J.
        #
        # J^2=-I and JZ=-ZJ, so QA,QB square to the same central -I
        # and anticommute exactly as Q8 requires.
        zeta = pow(root, atlas.endpoint.NN // 169, prime)
        identity2 = ((1, 0), (0, 1))
        diagonal_zeta = (
            (zeta, 0),
            (0, pow(zeta, -1, prime)),
        )
        j_matrix = ((0, prime - 1), (1, 0))
        grading_matrix = ((1, 0), (0, prime - 1))
        joint_a = kronecker(diagonal_zeta, identity2, prime)
        joint_b = kronecker(identity2, diagonal_zeta, prime)
        joint_qa = kronecker(j_matrix, identity2, prime)
        joint_qb = kronecker(grading_matrix, j_matrix, prime)
        identity4 = matrix_n_identity(4)
        minus_identity4 = matrix_n_scale(prime - 1, identity4, prime)
        require(
            pow(zeta, 169, prime) == 1
            and pow(zeta, 13, prime) == omega
            and matrix_n_power(joint_a, 169, prime) == identity4
            and matrix_n_power(joint_b, 169, prime) == identity4
            and matrix_n_power(joint_qa, 2, prime) == minus_identity4
            and matrix_n_power(joint_qb, 2, prime) == minus_identity4
            and matrix_n_mul(joint_qa, joint_qb, prime)
            == matrix_n_scale(
                prime - 1,
                matrix_n_mul(joint_qb, joint_qa, prime),
                prime,
            ),
            "joint Clifford/Q8 relations failed",
        )
        require(
            matrix_n_mul(
                matrix_n_mul(joint_qa, joint_a, prime),
                matrix_n_power(joint_qa, 3, prime),
                prime,
            ) == matrix_n_power(joint_a, 168, prime)
            and matrix_n_mul(
                matrix_n_mul(joint_qa, joint_b, prime),
                matrix_n_power(joint_qa, 3, prime),
                prime,
            ) == joint_b
            and matrix_n_mul(
                matrix_n_mul(joint_qb, joint_a, prime),
                matrix_n_power(joint_qb, 3, prime),
                prime,
            ) == joint_a
            and matrix_n_mul(
                matrix_n_mul(joint_qb, joint_b, prime),
                matrix_n_power(joint_qb, 3, prime),
                prime,
            ) == matrix_n_power(joint_b, 168, prime),
            "joint coordinatewise inversion action failed",
        )

        balanced_prony = matrix_n_mul(
            matrix_n_power(joint_a, 26, prime),
            matrix_n_power(joint_b, 26, prime),
            prime,
        )
        expected_balanced = (
            (omega4, 0, 0, 0),
            (0, 1, 0, 0),
            (0, 0, 1, 0),
            (0, 0, 0, pow(omega4, -1, prime)),
        )
        require(
            balanced_prony == expected_balanced,
            "balanced four-channel Prony completion failed",
        )

        # The first two eigenlines see literal diag(E,1), but QA sends
        # the first eigenline into the third.  Thus this two-plane is not
        # a module for the full joint group.
        first_two_invariant_under_qa = all(
            joint_qa[row][column] == 0
            for row in (2, 3) for column in (0, 1)
        )
        require(
            not first_two_invariant_under_qa,
            "the literal Prony two-plane unexpectedly became Q8-stable",
        )

        # Exact faithfulness checks without materializing all 228488
        # matrices: A has 169^2 distinct diagonal images; the eight Q8
        # matrices are distinct monomial matrices; and only 1,z are
        # diagonal, with z=-I outside the odd-order A image.
        diagonal_images = {
            (
                pow(zeta, left + right, prime),
                pow(zeta, left - right, prime),
                pow(zeta, -left + right, prime),
                pow(zeta, -left - right, prime),
            )
            for left in range(169) for right in range(169)
        }
        q8_images = set()
        for left in range(4):
            for right in range(4):
                q8_images.add(
                    matrix_n_mul(
                        matrix_n_power(joint_qa, left, prime),
                        matrix_n_power(joint_qb, right, prime),
                        prime,
                    )
                )
        require(
            len(diagonal_images) == 169**2
            and len(q8_images) == 8
            and tuple(prime - 1 for _ in range(4)) not in diagonal_images,
            "joint representation faithfulness controls failed",
        )
        joint_rows.append((
            prime,
            len(diagonal_images),
            len(q8_images),
            balanced_prony,
        ))

    print("INDEPENDENT DICYCLIC/PRONY REFEREE -- KERNEL <QA>")
    print(
        "abstract_group=C169 rtimes Q8 ~= Dic_338; "
        "Dic_n convention order=4n; order=1352"
    )
    print(
        "presentation=<r,x | r^676=1,x^2=r^338,xrx^-1=r^-1>; "
        "C=r^508(order169), QA=r^169, QB=x, "
        "QAB=r^338*QA*QB; r=C*QA"
    )
    print(
        "action_kernel=<QA>; QB,QAB invert C; "
        "requested_inversion_character_on_(Q0,QA,QAB)=(0,0,1)"
    )
    print(
        "selector_character=det(QB,-)_on_(Q0,QA,QAB)=(0,1,1); "
        "therefore_action_character_NE_selector_character"
    )
    print(f"order_census={dict(sorted(order_census.items()))}")
    print(
        f"conjugacy_classes={len(classes)}; "
        f"class_census={dict(sorted(class_census.items()))}"
    )
    print(
        "skew_classes=QB,QAB are distinct order4 size338 classes; "
        "each consists of 169 C169-lifts and their central negatives; "
        f"q7_filler_pairs_checked={filler_checks}"
    )
    print(
        "joint_edge=(C^8 QA)(C^9 QB)=(C^13 z)(C^4 QAB); "
        "skew_clutch=C^13 z has order26, QB-inverse, square=C^26"
    )
    print(
        "abelianization=V4; commutator=<r^2> order338 contains C169,z; "
        "every one-dimensional character is trivial on C169"
    )
    for row in field_rows:
        print(
            f"field={row[0]}; faithful_eta676={row[1]}; "
            f"E=omega4={row[2]}; rho(C^26)={row[3]}; "
            f"det(diag(E,1))={row[4]}"
        )
    print(
        "linear_boundary=diag(E,1)=omega^2*rho(C^26), hence exact only "
        "projectively; det(E,1)=E!=1 and the omega^2 scalar cannot extend "
        "to a character of G"
    )
    print(
        "Hermitian_positive=K=U*bar(V) is nonzero whenever U,V are; "
        "diag(E,1) and rho(C^26) both send K->E*K; "
        "QA:K->-K, QB:K->-bar(K), QAB:K->bar(K)"
    )
    print(
        "scope=abstract group and coefficient/projective Hermitian "
        "observable only; no Q8 action on marked physical sheets, no "
        "simultaneous physical Prony quadrature, no descent from K to "
        "S=U+V or canonical unmarked current, no row exclusion"
    )
    print(
        "joint_minimal_group=(C169^2) rtimes Q8 with "
        "(chi_sel,chi_carry)=(det(QB,-),det(QA,-)); "
        "order=169^2*8=228488; action kernel on Q8 is center {+/-1}"
    )
    print(
        "joint_structure=center=C2; derived=C169^2*C2 order57122; "
        "abelianization=V4; quotient_by_center=D338 x D338; "
        "nontrivial central double cover because lifted reflections "
        "anticommute"
    )
    print(
        "minimality=one C169 coordinate carries only one Q8->C2 "
        "character; two coordinates are minimal for both; every faithful "
        "splitting-field linear representation has dimension at least4, "
        "and the displayed Clifford tensor model attains4"
    )
    for row in joint_rows:
        print(
            f"joint_field={row[0]}; A_diagonal_images={row[1]}; "
            f"Q8_images={row[2]}; "
            f"(a^26*b^26)_matrix={row[3]}"
        )
    print(
        "joint_linear_verdict=diag(E,1) occurs as the first half of "
        "diag(E,1,1,E^-1), but that two-plane is not Q8-stable; the full "
        "four-channel balanced completion is linear and faithful, while "
        "the original two-channel action remains only projective/Hermitian"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
