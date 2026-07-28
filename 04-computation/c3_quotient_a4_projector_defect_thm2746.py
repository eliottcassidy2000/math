#!/usr/bin/env python3
"""Exact finite companion for THM-2746.

Classify affine lifts of C2*C3 -> C3 on V=F2^2 and compute the
equal-arm-to-charged Reynolds block.  All truth checks use explicit
exceptions and remain active under python -O.
"""

from fractions import Fraction
from itertools import product


V = ((0, 0), (1, 0), (0, 1), (1, 1))
ZERO = V[0]
I2 = ((1, 0), (0, 1))
T2 = ((0, 1), (1, 1))
I4P = tuple(range(4))


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
    group = {I4P}
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


def permutation_orbits(group):
    unseen = set(range(4))
    answer = []
    while unseen:
        x = min(unseen)
        orbit = {g[x] for g in group}
        answer.append(tuple(sorted(orbit)))
        unseen -= orbit
    return tuple(sorted((len(o) for o in answer)))


def cycle_type(p):
    unseen = set(range(4))
    sizes = []
    while unseen:
        x = min(unseen)
        y = x
        size = 0
        while y in unseen:
            unseen.remove(y)
            size += 1
            y = p[y]
        sizes.append(size)
    return tuple(sorted(sizes))


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
    rank = 0
    for col in range(len(rows[0])):
        pivot = next(
            (i for i in range(rank, len(rows)) if rows[i][col]), None
        )
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        scale = rows[rank][col]
        rows[rank] = [entry / scale for entry in rows[rank]]
        for i in range(len(rows)):
            if i == rank or not rows[i][col]:
                continue
            factor = rows[i][col]
            rows[i] = [
                rows[i][j] - factor * rows[rank][j]
                for j in range(len(rows[0]))
            ]
        rank += 1
    return rank


def frobenius_squared(a):
    return sum((entry * entry for row in a for entry in row), Fraction(0))


def reduce_mod_two(a):
    answer = []
    for row in a:
        reduced = []
        for entry in row:
            require(entry.denominator % 2 == 1, "odd projector denominator")
            reduced.append(entry.numerator % 2)
        answer.append(tuple(reduced))
    return tuple(answer)


def matrix_rank_mod_two(a):
    rows = [list(row) for row in a]
    rank = 0
    for col in range(len(rows[0])):
        pivot = next(
            (i for i in range(rank, len(rows)) if rows[i][col]), None
        )
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for i in range(len(rows)):
            if i != rank and rows[i][col]:
                rows[i] = [x ^ y for x, y in zip(rows[i], rows[rank])]
        rank += 1
    return rank


def matrix_vector_mod_two(a, v):
    return tuple(
        sum((a[i][j] * v[j] for j in range(len(v))), 0) % 2
        for i in range(len(a))
    )


def xor_vectors(*vectors):
    return tuple(sum(entries) % 2 for entries in zip(*vectors))


def matrix_vector(a, v):
    return tuple(
        sum((a[i][j] * v[j] for j in range(len(v))), Fraction(0))
        for i in range(len(a))
    )


def permute_vector(p, v):
    answer = [0] * len(v)
    for source, target in enumerate(p):
        answer[target] = v[source]
    return tuple(answer)


def determinant_pairing(x, y):
    """The canonical GL(2,F2)-invariant alternating pairing on V."""
    return (x[0] * y[1]) ^ (x[1] * y[0])


def cyclic_orbits(p):
    group = generated_group((p,))
    unseen = set(range(4))
    answer = []
    while unseen:
        x = min(unseen)
        orbit = frozenset(g[x] for g in group)
        answer.append(orbit)
        unseen -= orbit
    return tuple(sorted(answer, key=lambda o: (len(o), min(o))))


def contingency(sigma, tau):
    blocks = cyclic_orbits(tau)
    require(tuple(map(len, blocks)) == (1, 3), "C3 orbit partition 1+3")
    return tuple(
        tuple(len({sigma[x] for x in oi} & oj) for oj in blocks)
        for oi in blocks
    )


def main() -> None:
    t2 = mmul2(T2, T2)
    require(mmul2(t2, T2) == I2, "linear C3 generator")
    require({vadd(v, mact(T2, v)) for v in V} == set(V),
            "I+T is invertible")

    cocycles = tuple(product(V, repeat=2))
    coboundaries = frozenset(
        (ZERO, vadd(v, mact(T2, v))) for v in V
    )
    require(len(cocycles) == 16, "sixteen affine cocycles")
    require(len(coboundaries) == 4, "four affine coboundaries")
    require({b for _, b in coboundaries} == set(V),
            "every C3 translation is gauged away")

    zero_orders = []
    nonzero_orders = []
    zero_orbits = []
    nonzero_orbits = []
    zero_off = []
    nonzero_off = []
    zero_comm = []
    nonzero_comm = []
    nonzero_frob = []
    zero_binary_types = []
    nonzero_binary_types = []
    nonzero_contingencies = []
    nonzero_mod_two_ranks = []
    nonzero_mod_two_orbit_spans = []
    canonical_mod_two_lines = []
    leakage_characters = {b: {} for b in V}
    scaled_leakage_rows = 0

    for a, b in cocycles:
        sigma_perm = affine_perm(I2, a)
        tau_perm = affine_perm(T2, b)
        require(pcompose(sigma_perm, sigma_perm) == I4P,
                "affine C2 relation")
        require(pcompose(pcompose(tau_perm, tau_perm), tau_perm) == I4P,
                "affine C3 relation")
        require(len(generated_group((tau_perm,))) == 3,
                "isolated C3 subgroup order")
        st = pcompose(sigma_perm, tau_perm)
        require(pcompose(pcompose(st, st), st) == I4P,
                "A4 triangle relation (sigma tau)^3")

        sigma = perm_matrix(sigma_perm)
        tau = perm_matrix(tau_perm)
        identity = eye(4)
        projector = mscale(
            Fraction(1, 3),
            madd(madd(identity, tau), mmul(tau, tau)),
        )
        charged = msub(identity, projector)
        off = mmul(mmul(charged, sigma), projector)
        commutator = msub(mmul(sigma, projector), mmul(projector, sigma))
        require(mmul(projector, projector) == projector,
                "Reynolds idempotent")

        fixed_index = next(i for i in range(4) if tau_perm[i] == i)
        fixed_point = V[fixed_index]
        fixed_contrast = tuple(
            Fraction(3 if x == fixed_point else -1) for x in V
        )
        scaled_leakage = tuple(
            Fraction(3, 4) * entry
            for entry in matrix_vector(off, fixed_contrast)
        )
        expected_leakage = [Fraction(0)] * 4
        expected_leakage[V.index(vadd(fixed_point, a))] += 2
        expected_leakage[
            V.index(vadd(fixed_point, mact(T2, a)))
        ] -= 1
        expected_leakage[
            V.index(vadd(fixed_point, mact(t2, a)))
        ] -= 1
        require(scaled_leakage == tuple(expected_leakage),
                "integral scaled leakage formula")
        require(all(entry.denominator == 1 for entry in scaled_leakage),
                "scaled leakage is integral")

        leakage_mod_two = tuple(int(entry) % 2 for entry in scaled_leakage)
        character_table = tuple(
            determinant_pairing(a, vadd(x, fixed_point)) for x in V
        )
        require(leakage_mod_two == character_table,
                "scaled leakage is the determinant character")
        require(leakage_mod_two[fixed_index] == 0,
                "character vanishes at the marked origin")
        require(sum(leakage_mod_two) % 2 == 0,
                "character lies in the three-orbit augmentation plane")
        leakage_characters[b][a] = leakage_mod_two
        scaled_leakage_rows += 1

        group = generated_group((sigma_perm, tau_perm))
        if a == ZERO:
            zero_orders.append(len(group))
            zero_orbits.append(permutation_orbits(group))
            zero_off.append(matrix_rank(off))
            zero_comm.append(matrix_rank(commutator))
            zero_binary_types.append(cycle_type(sigma_perm))
        else:
            nonzero_orders.append(len(group))
            nonzero_orbits.append(permutation_orbits(group))
            nonzero_off.append(matrix_rank(off))
            nonzero_comm.append(matrix_rank(commutator))
            nonzero_frob.append(frobenius_squared(off))
            nonzero_binary_types.append(cycle_type(sigma_perm))
            nonzero_contingencies.append(contingency(sigma_perm, tau_perm))

            off_two = reduce_mod_two(off)
            tau_two = reduce_mod_two(tau)
            rank_two = matrix_rank_mod_two(off_two)
            nonzero_mod_two_ranks.append(rank_two)
            require(rank_two == 1, "mod-two projector image line")
            columns = tuple(
                tuple(off_two[i][j] for i in range(4)) for j in range(4)
            )
            image_vector = next(v for v in columns if any(v))
            image_orbit = (
                image_vector,
                matrix_vector_mod_two(tau_two, image_vector),
                matrix_vector_mod_two(
                    tau_two, matrix_vector_mod_two(tau_two, image_vector)
                ),
            )
            require(len(set(image_orbit)) == 3, "three charged image lines")
            require(xor_vectors(*image_orbit) == (0, 0, 0, 0),
                    "standard-plane orbit relation")
            orbit_span = matrix_rank_mod_two(tuple(zip(*image_orbit)))
            nonzero_mod_two_orbit_spans.append(orbit_span)
            require(orbit_span == 2, "standard-plane orbit span")
            if b == ZERO:
                canonical_mod_two_lines.append(image_vector)

    for b in V:
        tau_perm = affine_perm(T2, b)
        fixed_index = next(i for i in range(4) if tau_perm[i] == i)
        characters = leakage_characters[b]
        for a, c in product(V, repeat=2):
            require(
                characters[vadd(a, c)]
                == tuple(x ^ y for x, y in zip(characters[a], characters[c])),
                "leakage-character map is F2-linear",
            )
        for a in V:
            require(
                characters[mact(T2, a)]
                == permute_vector(tau_perm, characters[a]),
                "leakage-character map is C3-equivariant",
            )
        augmentation_plane = {
            bits
            for bits in product((0, 1), repeat=4)
            if bits[fixed_index] == 0 and sum(bits) % 2 == 0
        }
        require(set(characters.values()) == augmentation_plane,
                "leakage characters fill the augmentation plane")

    require(zero_orders == [3] * 4, "zero-class C3 images")
    require(nonzero_orders == [12] * 12, "nonzero-class A4 images")
    require(zero_orbits == [(1, 3)] * 4, "zero-class orbit type")
    require(nonzero_orbits == [(4,)] * 12, "nonzero-class transitivity")
    require(zero_binary_types == [(1, 1, 1, 1)] * 4,
            "zero-class binary generator is identity")
    require(nonzero_binary_types == [(2, 2)] * 12,
            "A4 binary generator is a double transposition")
    require(zero_off == [0] * 4, "zero-class off-diagonal ranks")
    require(nonzero_off == [1] * 12, "A4 off-diagonal ranks")
    require(zero_comm == [0] * 4, "zero-class commutator ranks")
    require(nonzero_comm == [2] * 12, "A4 commutator ranks")
    require(nonzero_frob == [Fraction(8, 9)] * 12,
            "A4 Hilbert--Schmidt defect")
    require(nonzero_contingencies == [((0, 1), (1, 2))] * 12,
            "A4 C3-orbit contingency")
    require(nonzero_mod_two_ranks == [1] * 12,
            "all A4 mod-two image ranks")
    require(nonzero_mod_two_orbit_spans == [2] * 12,
            "all A4 standard-plane orbit spans")
    require(canonical_mod_two_lines
            == [(0, 0, 1, 1), (0, 1, 0, 1), (0, 1, 1, 0)],
            "canonical three standard-plane lines")
    require(scaled_leakage_rows == 16, "all scaled leakage rows")

    invertible_linear_maps = tuple(
        matrix
        for entries in product((0, 1), repeat=4)
        for matrix in ((entries[:2], entries[2:]),)
        if len({mact(matrix, x) for x in V}) == 4
    )
    linear_reflections = tuple(
        matrix
        for matrix in invertible_linear_maps
        if mmul2(matrix, matrix) == I2
        and mmul2(mmul2(matrix, T2), matrix) == t2
    )
    require(len(linear_reflections) == 3,
            "three C3-inverting linear reflections")

    reflection_rows = 0
    reflection_completion_orders = []
    for b in V:
        tau_perm = affine_perm(T2, b)
        tau_inverse = pcompose(tau_perm, tau_perm)
        a4_group = generated_group(
            tuple(affine_perm(I2, a) for a in V) + (tau_perm,)
        )
        require(len(a4_group) == 12, "marked A4 subgroup order")
        compatible = []
        for reflection in linear_reflections:
            for translation in V:
                candidate = affine_perm(reflection, translation)
                if (
                    pcompose(candidate, candidate) == I4P
                    and pcompose(
                        pcompose(candidate, tau_perm), candidate
                    ) == tau_inverse
                ):
                    compatible.append(candidate)
                    require(cycle_type(candidate) == (1, 1, 2),
                            "compatible reflection is a transposition")
                    require(candidate not in a4_group,
                            "compatible reflection lies outside A4")
                    completion_order = len(
                        generated_group(tuple(a4_group) + (candidate,))
                    )
                    require(completion_order == 24,
                            "compatible reflection completes A4 to S4")
                    reflection_completion_orders.append(completion_order)
        require(len(compatible) == 3,
                "three affine reflections for each marked C3 lift")
        reflection_rows += len(compatible)

    require(reflection_rows == 12, "twelve compatible reflection rows")
    require(reflection_completion_orders == [24] * 12,
            "every compatible reflection gives S4")

    print("THM-2746 C3-QUOTIENT A4 PROJECTOR AUDIT")
    print("affine_lifts=16 coboundaries=4 H1_classes=4")
    print("restriction_H1=C2:V4 C3:0")
    print(f"zero_class_image_orders={zero_orders}")
    print(f"nonzero_class_image_orders={nonzero_orders}")
    print(f"zero_class_orbits={zero_orbits}")
    print(f"nonzero_class_orbits={nonzero_orbits}")
    print(f"zero_class_offdiag_ranks={zero_off}")
    print(f"nonzero_class_offdiag_ranks={nonzero_off}")
    print(f"zero_class_commutator_ranks={zero_comm}")
    print(f"nonzero_class_commutator_ranks={nonzero_comm}")
    print("nonzero_class_offdiag_frobenius_squared=[8/9]*12")
    print("nonzero_class_contingency=[[0,1],[1,2]]")
    print("nonzero_class_mod2_offdiag_ranks=[1]*12")
    print("nonzero_class_mod2_standard_orbit_spans=[2]*12")
    print("canonical_mod2_image_lines=[0011,0101,0110]")
    print("binary_cycle_types=zero:1^4 nonzero:2^2")
    print("triangle_relation=(sigma*tau)^3=1 for all 16 lifts")
    print("scaled_leakage_formula=(3/4)M_a*w_p=2e_(p+a)-e_(p+Ta)-e_(p+T2a)")
    print("scaled_leakage_mod2=det(a,x-p) on all 16 lifts")
    print("leakage_character_plane=F2-linear C3-equivariant isomorphism x4 gauges")
    print("compatible_reflections=12 per_C3_lift=3 cycle_type=1^2+2")
    print("reflection_completion_orders=[24]*12")
    print("boundary_scope=character_plane_to_ker(delta_mod2);unit_vs_saturation_unresolved")
    print("quartic_scope=marked_A4_detector_not_Keller_exclusion")
    print("LRC_scope=no_common_physical_involution_or_endpoint_current")
    print("PASS")


if __name__ == "__main__":
    main()
