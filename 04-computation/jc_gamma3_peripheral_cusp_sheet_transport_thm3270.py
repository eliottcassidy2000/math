#!/usr/bin/env python3
"""Exact finite companion for THM-3270.

The infinite-group statement in THM-3270 is proved by the double-coset
identity

    Gamma(3)\PSL_2(Z)/<T> = A_4/<T mod 3>.

This companion checks every finite consequence used by the theorem: the
four cusp/sheet cosets, their stabilizers, the four oriented opposite-face
cycles and their integral homology relation, equivariant transport, the
completion with the three matching directions to the twelve flag torsor,
and the exact holonomy obstruction inside A_4.
"""

from itertools import combinations, permutations


POINTS = (0, 1, 2, 3)  # 3 denotes infinity in P^1(F_3)
INF = 3
IDENTITY = POINTS
EDGES = tuple(combinations(POINTS, 2))
CHORDS = ((0, 1), (0, 2), (1, 2))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def compose(p, q):
    """Permutation product p*q: apply q, then p."""
    return tuple(p[q[i]] for i in POINTS)


def inverse(p):
    return tuple(p.index(i) for i in POINTS)


def parity(p):
    return sum(
        p[i] > p[j]
        for i in POINTS
        for j in POINTS
        if i < j
    ) % 2


def perm_order(p):
    x = IDENTITY
    for n in range(1, 25):
        x = compose(x, p)
        if x == IDENTITY:
            return n
    raise RuntimeError("permutation order exceeded group bound")


def generated_group(generators):
    group = {IDENTITY}
    changed = True
    while changed:
        old_size = len(group)
        current = tuple(group)
        for a in current:
            for b in tuple(generators) + current:
                group.add(compose(a, b))
                group.add(compose(b, a))
        changed = len(group) != old_size
    return frozenset(group)


def conjugate(a, x):
    return compose(compose(a, x), inverse(a))


def right_cosets(group, subgroup):
    unseen = set(group)
    cosets = []
    while unseen:
        g = min(unseen)
        coset = frozenset(compose(g, h) for h in subgroup)
        require(coset <= set(group), "right coset escaped ambient group")
        cosets.append(coset)
        unseen -= coset
    return tuple(cosets)


def oriented_face_chain(cycle):
    fixed = [x for x in POINTS if cycle[x] == x]
    require(len(fixed) == 1, "face cycle must fix exactly one sheet")
    x = min(y for y in POINTS if y != fixed[0])
    chain = [0] * len(EDGES)
    for _ in range(3):
        y = cycle[x]
        edge = tuple(sorted((x, y)))
        chain[EDGES.index(edge)] += 1 if x < y else -1
        x = y
    return tuple(chain)


def transport_chain(a, chain):
    out = [0] * len(EDGES)
    for coefficient, (x, y) in zip(chain, EDGES):
        if coefficient == 0:
            continue
        ax, ay = a[x], a[y]
        edge = tuple(sorted((ax, ay)))
        out[EDGES.index(edge)] += coefficient if ax < ay else -coefficient
    return tuple(out)


def det3(rows):
    a, b, c = rows
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def matmul(a, b):
    return tuple(
        tuple(sum(a[i][k] * b[k][j] for k in range(2)) for j in range(2))
        for i in range(2)
    )


def matpow(a, n):
    out = ((1, 0), (0, 1))
    for _ in range(n):
        out = matmul(out, a)
    return out


def main():
    # The compatible modular generators from THM-3141 on P^1(F_3).
    s = (3, 2, 1, 0)  # z -> -1/z
    r = (2, 1, 3, 0)  # z -> -1/(z+1)
    t = compose(s, r)  # z -> z+1 = (0 1 2), fixing infinity

    a4 = generated_group((s, r))
    require(len(a4) == 12, "mod-three image is not A4")
    require(all(parity(a) == 0 for a in a4), "image escaped A4")
    require((perm_order(s), perm_order(r), perm_order(t)) == (2, 3, 3),
            "modular generator orders are wrong")

    # T^3 is a nontrivial integral parabolic but vanishes modulo three.
    s_matrix = ((0, -1), (1, 0))
    r_matrix = ((0, -1), (1, 1))
    t_matrix = matmul(s_matrix, r_matrix)
    t3_matrix = matpow(t_matrix, 3)
    require(t_matrix in (((1, 1), (0, 1)), ((-1, -1), (0, -1))),
            "SR is not the standard parabolic up to central sign")
    if t_matrix[0][0] == -1:
        t_matrix = tuple(tuple(-entry for entry in row) for row in t_matrix)
        t3_matrix = matpow(t_matrix, 3)
    require(t3_matrix == ((1, 3), (0, 1)), "T^3 matrix is wrong")
    require(t3_matrix != ((1, 0), (0, 1)), "T^3 collapsed integrally")
    require(tuple(tuple(entry % 3 for entry in row) for row in t3_matrix)
            == ((1, 0), (0, 1)), "T^3 does not vanish modulo three")

    # Cusps of barGamma(3) are the four right cosets A4/<Tbar>.
    h = generated_group((t,))
    cosets = right_cosets(a4, h)
    require(len(h) == 3 and len(cosets) == 4, "cusp double-coset count failed")
    cusp_sheet = {}
    for coset in cosets:
        images = {g[INF] for g in coset}
        require(len(images) == 1, "cusp coset has no well-defined fixed sheet")
        sheet = next(iter(images))
        require(sheet not in cusp_sheet.values(), "two cusp cosets share a sheet")
        cusp_sheet[coset] = sheet
        stabilizer = frozenset(a for a in a4 if frozenset(compose(a, g) for g in coset) == coset)
        require(len(stabilizer) == 3, "cusp stabilizer is not C3")
        require(all(a[sheet] == sheet for a in stabilizer),
                "cusp stabilizer does not fix its sheet")

    # The positive parabolic conjugacy class gives the four coherently
    # oriented opposite faces, one for each fixed sheet.
    positive_cycles = {conjugate(a, t) for a in a4}
    require(len(positive_cycles) == 4, "positive C3 conjugacy class size failed")
    face_by_sheet = {}
    cycle_by_sheet = {}
    for cycle in positive_cycles:
        fixed = [x for x in POINTS if cycle[x] == x]
        require(len(fixed) == 1, "positive cycle lacks unique fixed sheet")
        sheet = fixed[0]
        require(sheet not in face_by_sheet, "two positive faces share a fixed sheet")
        face_by_sheet[sheet] = oriented_face_chain(cycle)
        cycle_by_sheet[sheet] = cycle

    require(set(face_by_sheet) == set(POINTS), "not every sheet has a face")
    relation = tuple(sum(face_by_sheet[v][i] for v in POINTS) for i in range(len(EDGES)))
    require(relation == (0,) * len(EDGES), "four peripheral face cycles do not sum to zero")
    chord_indices = tuple(EDGES.index(edge) for edge in CHORDS)
    determinants = []
    for triple in combinations(POINTS, 3):
        rows = [tuple(face_by_sheet[v][i] for i in chord_indices) for v in triple]
        determinants.append(det3(rows))
    require(all(abs(value) == 1 for value in determinants),
            "three peripheral cycles fail to form an integral H1 basis")
    for v, w in combinations(POINTS, 2):
        require(face_by_sheet[v] != face_by_sheet[w], "peripheral cycles collided")
        require(face_by_sheet[v] != tuple(-x for x in face_by_sheet[w]),
                "unoriented peripheral rays collided")

    transport_checks = 0
    for a in a4:
        for v in POINTS:
            require(transport_chain(a, face_by_sheet[v]) == face_by_sheet[a[v]],
                    "cusp/sheet transport is not A4-equivariant")
            require(conjugate(a, cycle_by_sheet[v]) == cycle_by_sheet[a[v]],
                    "positive peripheral generator transport failed")
            transport_checks += 1

    # The determinant-minus-one modular reflection extends A4 to
    # PGL_2(F_3)=S4.  It preserves the peripheral subgroup/ray but reverses
    # the positive face orientation.  Thus S4 transport is available only
    # after the generator orientation is forgotten.
    reflection = (0, 2, 1, 3)  # z -> -z
    s4 = generated_group(tuple(a4) + (reflection,))
    require(len(s4) == 24 and parity(reflection) == 1,
            "odd modular reflection did not generate S4")
    unoriented_transport_checks = 0
    for a in s4:
        sign = 1 if parity(a) == 0 else -1
        for v in POINTS:
            transported = transport_chain(a, face_by_sheet[v])
            target = tuple(sign * value for value in face_by_sheet[a[v]])
            require(transported == target,
                    "PGL transport failed on an unoriented peripheral ray")
            unoriented_transport_checks += 1

    # Cusp owner plus one of the three matching directions is exactly the
    # twelve-point flag torsor used by THM-3072.
    directions = tuple(a for a in a4 if perm_order(a) == 2)
    require(len(directions) == 3, "nonzero V4 direction count failed")
    flags = {(v, d) for v in POINTS for d in directions}
    base = (INF, min(directions))
    orbit = {
        (a[base[0]], conjugate(a, base[1]))
        for a in a4
    }
    require(len(flags) == 12 and orbit == flags,
            "cusp-owner times direction is not the regular A4 flag torsor")

    # Exhaust the A4 holonomy criterion.  Every subgroup is generated by at
    # most two elements here; there are exactly ten of them.
    subgroups = {
        generated_group(gens)
        for size in range(3)
        for gens in combinations(tuple(a4), size)
    }
    require(len(subgroups) == 10, "A4 subgroup census failed")
    c3_stabilizers = {
        frozenset(a for a in a4 if a[v] == v)
        for v in POINTS
    }
    require(all(len(stab) == 3 for stab in c3_stabilizers),
            "point stabilizer is not C3")
    section_subgroups = 0
    for subgroup in subgroups:
        fixed = {v for v in POINTS if all(a[v] == v for a in subgroup)}
        contained = any(subgroup <= stab for stab in c3_stabilizers)
        require(bool(fixed) == contained,
                "fixed-sheet criterion differs from containment in C3")
        section_subgroups += bool(fixed)

    distinct_c3_pairs = 0
    for first, second in combinations(c3_stabilizers, 2):
        require(generated_group(tuple(first | second)) == a4,
                "two distinct pure-C3 stabilizers did not generate A4")
        distinct_c3_pairs += 1

    require(not any(s[v] == v for v in POINTS),
            "V4 hostile unexpectedly fixes a sheet")
    require(sum(t[v] == v for v in POINTS) == 1,
            "parabolic quotient should fix exactly one cusp")
    require(compose(compose(t, t), t) == IDENTITY,
            "kernel parabolic does not act trivially on cusp classes")

    print("THM3270 gamma3 peripheral cusp / fixed-sheet transport exact companion")
    print(f"A4_order={len(a4)} cusp_cosets={len(cosets)} C3_stabilizers={len(c3_stabilizers)}")
    print(f"T3_matrix={t3_matrix} nontrivial_integrally=True identity_mod3=True")
    print("cusp_sheet_bijection=4 stabilizers_match=True")
    print("opposite_face_cycles=4 H1_rank=3 sum_zero=True any_three_unimodular=True")
    print(f"equivariant_transport_checks={transport_checks}")
    print(f"unoriented_PGL_transport_checks={unoriented_transport_checks} odd_reverses_orientation=True")
    print(f"cusp_owner_x_direction_flags={len(flags)} regular_A4_torsor=True")
    print(f"A4_subgroups={len(subgroups)} section_subgroups={section_subgroups} section_iff_in_C3=True")
    print(f"distinct_C3_pairs={distinct_c3_pairs} all_generate_A4=True")
    print("hostiles: V4_loop_fixed_cusps=0 kernel_T3_fixed_cusps=4")


if __name__ == "__main__":
    main()
