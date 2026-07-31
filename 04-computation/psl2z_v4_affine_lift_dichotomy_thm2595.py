#!/usr/bin/env python3
"""Exact finite audit for THM-2595.

The calculation uses only tuples of bits.  It classifies affine lifts of the
standard quotient C2 * C3 -> GL(2,F2) inside AGL(2,F2), computes cocycles and
coboundaries, checks the local fixed-locus criterion, and verifies that S3
cannot carry a left-translation-invariant tournament.
"""

from itertools import permutations, product


V = ((0, 0), (1, 0), (0, 1), (1, 1))
ZERO = V[0]
ID2 = ((1, 0), (0, 1))
S = ((0, 1), (1, 0))
T = ((0, 1), (1, 1))
ID4 = tuple(range(4))


def vadd(x, y):
    return (x[0] ^ y[0], x[1] ^ y[1])


def mact(a, x):
    return (
        (a[0][0] * x[0]) ^ (a[0][1] * x[1]),
        (a[1][0] * x[0]) ^ (a[1][1] * x[1]),
    )


def mmul(a, b):
    return tuple(
        tuple(
            (a[i][0] * b[0][j]) ^ (a[i][1] * b[1][j])
            for j in range(2)
        )
        for i in range(2)
    )


def pcompose(p, q):
    """Permutation product p after q."""
    return tuple(p[q[i]] for i in range(4))


def affine_perm(a, c):
    return tuple(V.index(vadd(mact(a, x), c)) for x in V)


def generated_group(generators):
    group = {ID4}
    frontier = list(generators)
    while frontier:
        x = frontier.pop()
        if x in group:
            continue
        old = tuple(group)
        group.add(x)
        frontier.extend(pcompose(x, g) for g in old)
        frontier.extend(pcompose(g, x) for g in old)
    return frozenset(group)


def orbit_partition(group):
    unseen = set(range(4))
    answer = []
    while unseen:
        x = min(unseen)
        orbit = frozenset(g[x] for g in group)
        unseen.difference_update(orbit)
        answer.append(tuple(sorted(orbit)))
    return tuple(sorted(answer, key=lambda part: (len(part), part)))


def fixed_points(p):
    return frozenset(i for i in range(4) if p[i] == i)


def determinant(a):
    return (a[0][0] * a[1][1]) ^ (a[0][1] * a[1][0])


def matrix_order(a):
    x = ID2
    for order in range(1, 7):
        x = mmul(x, a)
        if x == ID2:
            return order
    raise AssertionError("matrix order exceeded |GL(2,F2)|")


def main():
    gl2 = tuple(
        a
        for entries in product((0, 1), repeat=4)
        for a in (((entries[0], entries[1]), (entries[2], entries[3])),)
        if determinant(a) == 1
    )
    assert len(gl2) == 6
    assert sorted(matrix_order(a) for a in gl2) == [1, 2, 2, 2, 3, 3]
    assert mmul(S, S) == ID2
    assert mmul(mmul(T, T), T) == ID2
    assert generated_group((affine_perm(S, ZERO), affine_perm(T, ZERO))) == frozenset(
        affine_perm(a, ZERO) for a in gl2
    )

    # A cocycle on C2*C3 is freely determined by a=c(s), b=c(t), subject
    # only to (1+S)a=0 and (1+T+T^2)b=0.
    t2 = mmul(T, T)
    cocycles = tuple(
        (a, b)
        for a, b in product(V, repeat=2)
        if vadd(a, mact(S, a)) == ZERO
        and vadd(vadd(b, mact(T, b)), mact(t2, b)) == ZERO
    )
    coboundaries = frozenset(
        (vadd(v, mact(S, v)), vadd(v, mact(T, v))) for v in V
    )
    assert len(cocycles) == 8
    assert len(coboundaries) == 4

    # The restrictions to both free factors have zero H^1.
    z1_c2 = frozenset(a for a in V if vadd(a, mact(S, a)) == ZERO)
    b1_c2 = frozenset(vadd(v, mact(S, v)) for v in V)
    z1_c3 = frozenset(
        b for b in V if vadd(vadd(b, mact(T, b)), mact(t2, b)) == ZERO
    )
    b1_c3 = frozenset(vadd(v, mact(T, v)) for v in V)
    assert z1_c2 == b1_c2 and len(z1_c2) == 2
    assert z1_c3 == b1_c3 and len(z1_c3) == 4

    rows = []
    for a, b in cocycles:
        ps = affine_perm(S, a)
        pt = affine_perm(T, b)
        group = generated_group((ps, pt))
        fs = fixed_points(ps)
        ft = fixed_points(pt)
        is_boundary = (a, b) in coboundaries
        assert len(fs) == 2
        assert len(ft) == 1
        assert is_boundary == bool(fs & ft)
        if is_boundary:
            assert len(group) == 6
            assert tuple(sorted(map(len, orbit_partition(group)))) == (1, 3)
        else:
            assert len(group) == 24
            assert orbit_partition(group) == ((0, 1, 2, 3),)
        rows.append(
            (
                a,
                b,
                "zero" if is_boundary else "nonzero",
                tuple(V[i] for i in sorted(fs)),
                V[next(iter(ft))],
                len(group),
                tuple(sorted(map(len, orbit_partition(group)))),
            )
        )

    # AGL(2,F2) is the full permutation group on its four affine points.
    agl = frozenset(affine_perm(a, c) for a in gl2 for c in V)
    assert len(agl) == 24
    assert agl == frozenset(permutations(range(4)))

    # A left-invariant tournament on S3 is impossible.  It already fails on
    # each involution u: e -> u would translate by u to u -> e, and conversely.
    linear_s3 = generated_group((affine_perm(S, ZERO), affine_perm(T, ZERO)))
    involutions = tuple(
        p for p in linear_s3 if p != ID4 and pcompose(p, p) == ID4
    )
    assert len(linear_s3) == 6
    assert len(involutions) == 3
    left_invariant_tournament_possible = len(involutions) == 0
    assert not left_invariant_tournament_possible

    print("THM-2595 exact audit")
    print("GL(2,F2): order=6, element_orders=[1,2,2,2,3,3]")
    print("presentation: S^2=T^3=1, <S,T>=GL(2,F2)")
    print(f"Z1(C2*C3,V): {len(cocycles)} cocycles")
    print(f"B1(C2*C3,V): {len(coboundaries)} coboundaries")
    print("H1(C2*C3,V): order=2")
    print(f"restrictions: |Z1(C2)|=|B1(C2)|={len(z1_c2)}")
    print(f"restrictions: |Z1(C3)|=|B1(C3)|={len(z1_c3)}")
    print("rows: c(s), c(t), H1-class, Fix(s-affine), Fix(t-affine), image, orbits")
    for row in rows:
        print(row)
    print("class census: zero -> 4 lifts of image 6 and orbit sizes 1+3")
    print("class census: nonzero -> 4 lifts of image 24 and orbit size 4")
    print("AGL(2,F2)=Sym(V): order=24")
    print("left-invariant tournament on S3: impossible (3 involution reversals)")
    print("FAILED CHECKS: NONE")


if __name__ == "__main__":
    main()
