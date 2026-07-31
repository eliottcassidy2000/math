#!/usr/bin/env python3
"""Exact finite companion for THM-2743.

Enumerate the eight affine lifts from THM-2595 and compute the rational
C3 Reynolds projector, its C2 off-diagonal block, the projector commutator,
and the S3 relation defect.  All truth checks use explicit exceptions.
"""

from fractions import Fraction
from itertools import product


V = ((0, 0), (1, 0), (0, 1), (1, 1))
ZERO = V[0]
S2 = ((0, 1), (1, 0))
T2 = ((0, 1), (1, 1))
ID2 = ((1, 0), (0, 1))
ID4_PERM = tuple(range(4))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def vadd(x, y):
    return (x[0] ^ y[0], x[1] ^ y[1])


def mact(a, x):
    return (
        (a[0][0] * x[0]) ^ (a[0][1] * x[1]),
        (a[1][0] * x[0]) ^ (a[1][1] * x[1]),
    )


def mmul2(a, b):
    return tuple(
        tuple(
            (a[i][0] * b[0][j]) ^ (a[i][1] * b[1][j])
            for j in range(2)
        )
        for i in range(2)
    )


def pcompose(p, q):
    """Permutation p after q."""
    return tuple(p[q[i]] for i in range(4))


def affine_perm(a, c):
    return tuple(V.index(vadd(mact(a, x), c)) for x in V)


def generated_group(generators):
    group = {ID4_PERM}
    frontier = list(generators)
    while frontier:
        x = frontier.pop()
        if x in group:
            continue
        group.add(x)
        current = tuple(group)
        frontier.extend(pcompose(x, g) for g in current)
        frontier.extend(pcompose(g, x) for g in current)
    return frozenset(group)


def eye(n):
    return tuple(
        tuple(Fraction(int(i == j)) for j in range(n)) for i in range(n)
    )


def madd(a, b):
    return tuple(
        tuple(a[i][j] + b[i][j] for j in range(len(a[0])))
        for i in range(len(a))
    )


def msub(a, b):
    return tuple(
        tuple(a[i][j] - b[i][j] for j in range(len(a[0])))
        for i in range(len(a))
    )


def mscale(c, a):
    return tuple(tuple(c * entry for entry in row) for row in a)


def mmul(a, b):
    return tuple(
        tuple(
            sum((a[i][k] * b[k][j] for k in range(len(b))), Fraction(0))
            for j in range(len(b[0]))
        )
        for i in range(len(a))
    )


def perm_matrix(p):
    answer = [[Fraction(0) for _ in range(4)] for _ in range(4)]
    for source, target in enumerate(p):
        answer[target][source] = Fraction(1)
    return tuple(tuple(row) for row in answer)


def matrix_rank(a):
    rows = [list(row) for row in a]
    m = len(rows)
    n = len(rows[0])
    rank = 0
    for col in range(n):
        pivot = next((i for i in range(rank, m) if rows[i][col]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        scale = rows[rank][col]
        rows[rank] = [entry / scale for entry in rows[rank]]
        for i in range(m):
            if i == rank or not rows[i][col]:
                continue
            factor = rows[i][col]
            rows[i] = [rows[i][j] - factor * rows[rank][j] for j in range(n)]
        rank += 1
        if rank == m:
            break
    return rank


def frobenius_squared(a):
    return sum((entry * entry for row in a for entry in row), Fraction(0))


def main() -> None:
    t2 = mmul2(T2, T2)
    require(mmul2(S2, S2) == ID2, "linear C2 generator")
    require(mmul2(t2, T2) == ID2, "linear C3 generator")

    cocycles = tuple(
        (a, b)
        for a, b in product(V, repeat=2)
        if vadd(a, mact(S2, a)) == ZERO
        and vadd(vadd(b, mact(T2, b)), mact(t2, b)) == ZERO
    )
    coboundaries = frozenset(
        (vadd(v, mact(S2, v)), vadd(v, mact(T2, v))) for v in V
    )
    require(len(cocycles) == 8, "eight affine cocycles")
    require(len(coboundaries) == 4, "four affine coboundaries")

    zero_off = []
    nonzero_off = []
    zero_comm = []
    nonzero_comm = []
    zero_relation = []
    nonzero_relation = []
    nonzero_frob = []
    zero_images = []
    nonzero_images = []

    for a, b in cocycles:
        sigma_perm = affine_perm(S2, a)
        tau_perm = affine_perm(T2, b)
        require(pcompose(sigma_perm, sigma_perm) == ID4_PERM,
                "affine involution")
        require(pcompose(pcompose(tau_perm, tau_perm), tau_perm) == ID4_PERM,
                "affine order-three generator")
        require(len(generated_group((sigma_perm,))) == 2,
                "affine involution subgroup order")
        require(len(generated_group((tau_perm,))) == 3,
                "affine order-three subgroup order")

        sigma = perm_matrix(sigma_perm)
        tau = perm_matrix(tau_perm)
        identity = eye(4)
        tau_sq = mmul(tau, tau)
        projector = mscale(Fraction(1, 3), madd(madd(identity, tau), tau_sq))
        charged = msub(identity, projector)
        off = mmul(mmul(charged, sigma), projector)
        commutator = msub(mmul(sigma, projector), mmul(projector, sigma))
        relation = msub(mmul(mmul(sigma, tau), sigma), tau_sq)

        require(mmul(projector, projector) == projector, "Reynolds idempotent")
        require(mmul(tau, projector) == projector, "projector image is C3 fixed")
        require(mmul(charged, off) == off, "off-diagonal output is charged")
        require(mmul(off, projector) == off, "off-diagonal input is fixed")

        image_order = len(generated_group((sigma_perm, tau_perm)))
        is_zero = (a, b) in coboundaries
        if is_zero:
            zero_off.append(matrix_rank(off))
            zero_comm.append(matrix_rank(commutator))
            zero_relation.append(matrix_rank(relation))
            zero_images.append(image_order)
            require(frobenius_squared(off) == 0, "zero-class off block")
        else:
            nonzero_off.append(matrix_rank(off))
            nonzero_comm.append(matrix_rank(commutator))
            nonzero_relation.append(matrix_rank(relation))
            nonzero_frob.append(frobenius_squared(off))
            nonzero_images.append(image_order)

    require(zero_off == [0, 0, 0, 0], "zero-class off ranks")
    require(nonzero_off == [1, 1, 1, 1], "nonzero-class off ranks")
    require(zero_comm == [0, 0, 0, 0], "zero-class commutator ranks")
    require(nonzero_comm == [2, 2, 2, 2], "nonzero-class commutator ranks")
    require(zero_relation == [0, 0, 0, 0], "zero-class S3 relation")
    require(nonzero_relation == [2, 2, 2, 2], "nonzero-class relation defect")
    require(nonzero_frob == [Fraction(8, 9)] * 4,
            "nonzero-class Hilbert-Schmidt norm")
    require(zero_images == [6, 6, 6, 6], "zero-class S3 images")
    require(nonzero_images == [24, 24, 24, 24], "nonzero-class S4 images")

    print("THM-2743 C2-C3 OFF-DIAGONAL PROJECTOR AUDIT")
    print("affine_lifts=8 coboundaries=4 nonzero_H1=4")
    print(f"zero_class_offdiag_ranks={zero_off}")
    print(f"nonzero_class_offdiag_ranks={nonzero_off}")
    print(f"zero_class_commutator_ranks={zero_comm}")
    print(f"nonzero_class_commutator_ranks={nonzero_comm}")
    print(f"zero_class_S3_relation_defect_ranks={zero_relation}")
    print(f"nonzero_class_S3_relation_defect_ranks={nonzero_relation}")
    print("nonzero_class_offdiag_frobenius_squared=[8/9, 8/9, 8/9, 8/9]")
    print(f"image_orders=zero:{zero_images} nonzero:{nonzero_images}")
    print("equal_arm_charged_component=(I-Pi3)*S*Pi3")
    print("S3_normalizer_boundary=offdiag_zero")
    print("LRC_scope=no_common_physical_C2_action_supplied")
    print("JC_scope=finite_marked_resolvent_detector_not_Keller_exclusion")
    print("PASS")


if __name__ == "__main__":
    main()
