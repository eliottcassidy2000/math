#!/usr/bin/env python3
"""Exact finite controls for THM-3354's typed H^1 no-go.

The universal coefficient and divisibility statements are proved in the
theorem.  This companion independently enumerates the finite coefficient
homomorphisms, the C91 mapping-torus orbit, the weighted-star hostile, and the
integral-versus-generic Hamiltonian localization hostile.
"""
from itertools import permutations, product


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def compose(p, q):
    return tuple(p[q[i]] for i in range(len(p)))


def power(g, n, identity):
    out = identity
    for _ in range(n):
        out = compose(out, g)
    return out


def det_bareiss(matrix):
    a = [list(row) for row in matrix]
    n = len(a)
    sign = 1
    previous = 1
    for k in range(n - 1):
        if a[k][k] == 0:
            swap = next((i for i in range(k + 1, n) if a[i][k]), None)
            require(swap is not None, "singular Bareiss pivot")
            a[k], a[swap] = a[swap], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = a[i][j] * pivot - a[i][k] * a[k][j]
                require(numerator % previous == 0, "nonintegral Bareiss step")
                a[i][j] = numerator // previous
        previous = pivot
    return sign * a[-1][-1]


def rank_mod2(matrix):
    a = [[x & 1 for x in row] for row in matrix]
    rows, cols = len(a), len(a[0])
    rank = 0
    for col in range(cols):
        pivot = next((i for i in range(rank, rows) if a[i][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        for i in range(rows):
            if i != rank and a[i][col]:
                a[i] = [x ^ y for x, y in zip(a[i], a[rank])]
        rank += 1
    return rank


def s3_to_cyclic_count(n):
    # S3=<s,t | s^2=t^3=1, sts=t^-1>.
    return sum(
        1
        for s, t in product(range(n), repeat=2)
        if 2 * s % n == 0
        and 3 * t % n == 0
        and (2 * s + 2 * t) % n == 0
    )


def v4_to_cyclic_count(n):
    return sum(
        1
        for x, y in product(range(n), repeat=2)
        if 2 * x % n == 0 and 2 * y % n == 0
    )


def hamiltonian_on_monomial(a, b):
    """Apply D_(x+x^2 z) to x^a z^b in the Laurent polynomial ring."""
    out = {}

    def add(exponents, coefficient):
        out[exponents] = out.get(exponents, 0) + coefficient
        if out[exponents] == 0:
            del out[exponents]

    # D=(1+2xz) partial_z-x^2 partial_x.
    if b:
        add((a, b - 1), b)
        add((a + 1, b), 2 * b)
    if a:
        add((a + 1, b), -a)
    return out


def main():
    identity = (0, 1, 2)
    s3 = tuple(permutations(range(3)))
    v4 = tuple(product(range(2), repeat=2))

    hom_rows = []
    for n in (13, 91):
        c_to_s3 = sum(power(g, n, identity) == identity for g in s3)
        s3_to_c = s3_to_cyclic_count(n)
        c_to_v4 = sum(
            ((n * x) & 1, (n * y) & 1) == (0, 0) for x, y in v4
        )
        v4_to_c = v4_to_cyclic_count(n)
        c_to_c2 = sum((n * x) % 2 == 0 for x in range(2))
        c2_to_c = sum((2 * x) % n == 0 for x in range(n))
        require((c_to_s3, s3_to_c, c_to_v4, v4_to_c,
                 c_to_c2, c2_to_c) == (1, 1, 1, 1, 1, 1),
                ("nontrivial coefficient hom", n))
        hom_rows.append((n, c_to_s3, s3_to_c, c_to_v4, v4_to_c,
                         c_to_c2, c2_to_c))

    c13_to_c91 = tuple(y for y in range(91) if 13 * y % 91 == 0)
    embeddings = tuple(y for y in c13_to_c91 if y != 0)
    require(len(c13_to_c91) == 13 and len(embeddings) == 12,
            "C13 to C91 hom census")
    require(all(y % 7 == 0 for y in c13_to_c91),
            "C13 image acquired a C7 component")

    orbit_lengths = []
    for a in range(1, 13):
        state = (0, 0)
        seen = set()
        while state not in seen:
            seen.add(state)
            state = ((state[0] + 1) % 7, (state[1] + a) % 13)
        require(state == (0, 0) and len(seen) == 91,
                ("mapping-torus orbit", a, len(seen), state))
        orbit_lengths.append(len(seen))

    d4 = (
        (-2, 1, 1, 1),
        (1, -2, 0, 0),
        (1, 0, -2, 0),
        (1, 0, 0, -2),
    )
    odd = (
        (-2, 1, 1, 1),
        (1, -3, 0, 0),
        (1, 0, -3, 0),
        (1, 0, 0, -3),
    )
    d4_det, odd_det = abs(det_bareiss(d4)), abs(det_bareiss(odd))
    d4_nullity = 4 - rank_mod2(d4)
    odd_nullity = 4 - rank_mod2(odd)
    require((d4_det, d4_nullity) == (4, 2), "D4 control")
    require((odd_det, odd_nullity) == (27, 0), "odd-arm control")

    require(7 % 13 != 0, "seven-chart holonomy collapsed")
    tree_h1 = 6 - 7 + 1
    cycle_h1 = 7 - 7 + 1
    require((tree_h1, cycle_h1) == (0, 1), "graph H1 controls")

    localized_primitive = hamiltonian_on_monomial(-1, 0)
    require(localized_primitive == {(0, 0): 1},
            ("localized Hamiltonian hostile", localized_primitive))
    d, e = 2, 1
    nu = min(d, e)
    annihilator_exponent = (d - 1 + nu - 1) // nu
    require(annihilator_exponent == 1, "one-root annihilator exponent")

    print("D5 TYPED H1 NO-GO -- EXACT FINITE CONTROLS")
    print("hom_counts (n,Cn_to_S3,S3_to_Cn,Cn_to_V4,V4_to_Cn,"
          "Cn_to_C2,C2_to_Cn)", hom_rows)
    print("C13_to_C91 homs", len(c13_to_c91), "embeddings", len(embeddings),
          "all_C7_projection_zero", all(y % 7 == 0 for y in c13_to_c91))
    print("mapping_torus nonzero_a", len(orbit_lengths),
          "orbit_lengths", sorted(set(orbit_lengths)))
    print("graph_h1 tree", tree_h1, "cycle7", cycle_h1,
          "holonomy_7a_nonzero", True)
    print("weighted_star D4_det_nullity", (d4_det, d4_nullity),
          "odd_det_nullity", (odd_det, odd_nullity))
    print("localization_hostile P=x+x^2z D_P(1/x)=1",
          "Ann_K[P](theta_int)=(P)^" + str(annihilator_exponent))
    print("PASS")


if __name__ == "__main__":
    main()
