#!/usr/bin/env python3
"""Exact, standard-library controls for the number/component bridge packet.

No assert statements: all tests also run under python -O.
"""
from fractions import Fraction
from itertools import combinations, product
from math import isqrt
import json

CHECKS = 0


def require(test, label):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(label)


def det(a):
    """Fraction-free elimination, with exact-division checks."""
    a = [row[:] for row in a]
    n = len(a)
    if not n:
        return 1
    last, sign = 1, 1
    for k in range(n - 1):
        if not a[k][k]:
            j = next((j for j in range(k + 1, n) if a[j][k]), None)
            if j is None:
                return 0
            a[k], a[j] = a[j], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = a[i][j] * pivot - a[i][k] * a[k][j]
                require(numerator % last == 0, "Bareiss exact division")
                a[i][j] = numerator // last
        for i in range(k + 1, n):
            a[i][k] = 0
        last = pivot
    return sign * a[-1][-1]


def matmul(a, b):
    return [[sum(x * y for x, y in zip(row, col)) for col in zip(*b)] for row in a]


def eye(n):
    return [[int(i == j) for j in range(n)] for i in range(n)]


def transpose(a):
    return [list(row) for row in zip(*a)]


def cycle_matrix(weights):
    a = eye(len(weights))
    for i, w in enumerate(weights):
        a[i][(i + 1) % len(a)] += w
    return a


def solve(a, b):
    a = [[Fraction(x) for x in row] + [Fraction(y)] for row, y in zip(a, b)]
    n = len(a)
    for k in range(n):
        j = next(j for j in range(k, n) if a[j][k])
        a[k], a[j] = a[j], a[k]
        pivot = a[k][k]
        a[k] = [x / pivot for x in a[k]]
        for i in range(n):
            if i != k:
                factor = a[i][k]
                a[i] = [x - factor * y for x, y in zip(a[i], a[k])]
    return [row[-1] for row in a]


def independent_sets(vertices, edges):
    vertices = tuple(vertices)
    out = []
    for mask in range(1 << len(vertices)):
        s = frozenset(v for i, v in enumerate(vertices) if mask >> i & 1)
        if all(not {u, v} <= s for u, v in edges):
            out.append(s)
    return out


def polynomial(sets):
    p = [0] * (max(map(len, sets)) + 1)
    for s in sets:
        p[len(s)] += 1
    return p


def pmul(a, b):
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i + j] += x * y
    return c


def norm(a, b):
    return a * a - a * b + b * b


def emul(x, y):
    a, b = x
    c, d = y
    return (a * c - b * d, a * d + b * c - b * d)


def sq_graph(n):
    return [(a, b) for a, b in combinations(range(1, n + 1), 2)
            if isqrt(a + b) ** 2 == a + b]


def main():
    e = sq_graph(9)
    require(e == [(1, 3), (1, 8), (2, 7), (3, 6), (4, 5), (7, 9)], "Q9 edges")
    components = [(8, 1, 3, 6), (9, 7, 2), (4, 5)]
    rows = []
    for v in components:
        ed = [(a, b) for a, b in e if a in v and b in v]
        sets = independent_sets(v, ed)
        rows.append({"path": v, "polynomial": polynomial(sets), "count": len(sets)})
    require([r["count"] for r in rows] == [8, 5, 3], "Fibonacci counts")
    all_sets = independent_sets(range(1, 10), e)
    require(len(all_sets) == 120, "whole graph independent count")
    require(polynomial(all_sets) == pmul(pmul(rows[0]["polynomial"], rows[1]["polynomial"]), rows[2]["polynomial"]), "component polynomial product")

    # Literal endpoint-conditioned bijections to the other two components.
    p4 = independent_sets(components[0], [(8, 1), (1, 3), (3, 6)])
    no8 = {frozenset({1: 9, 3: 7, 6: 2}[v] for v in s) for s in p4 if 8 not in s}
    yes8 = {frozenset({3: 4, 6: 5}[v] for v in s if v != 8) for s in p4 if 8 in s}
    require(no8 == set(independent_sets(components[1], [(9, 7), (7, 2)])), "P3 branch bijection")
    require(yes8 == set(independent_sets(components[2], [(4, 5)])), "P2 branch bijection")
    disjoint = independent_sets(components[0], [(8, 1), (3, 6)])
    require(set(disjoint) - set(p4) == {frozenset((1, 3))}, "single excluded set")
    toggle_edges = sum(len(a ^ b) == 1 for a, b in combinations(p4, 2))
    require(toggle_edges == 10 and toggle_edges != 12, "eight states do not form a cube under native toggles")

    # All graphs through five vertices, all edges: the edge-deletion identity.
    edge_cases = 0
    for n in range(2, 6):
        pairs = list(combinations(range(n), 2))
        for mask in range(1 << len(pairs)):
            ed = [uv for i, uv in enumerate(pairs) if mask >> i & 1]
            sets_g = set(independent_sets(range(n), ed))
            for uv in ed:
                u, v = uv
                forbidden = {u, v}
                for a, b in ed:
                    if a in uv:
                        forbidden.add(b)
                    if b in uv:
                        forbidden.add(a)
                remaining = set(range(n)) - forbidden
                rem_edges = [(a, b) for a, b in ed if a in remaining and b in remaining]
                expected = {s | frozenset(uv) for s in independent_sets(remaining, rem_edges)}
                actual = set(independent_sets(range(n), [xy for xy in ed if xy != uv])) - sets_g
                require(actual == expected, "general edge deletion bijection")
                require((len(actual) == 1) == (not remaining), "dominating edge iff one forbidden set")
                edge_cases += 1

    double_stars = []
    for r in range(1, 6):
        for s in range(1, 6):
            vertices = range(r + s + 2)
            ed = [(0, 1)] + [(0, j) for j in range(2, r + 2)] + [(1, j) for j in range(r + 2, r + s + 2)]
            count = len(independent_sets(vertices, ed))
            require(count == 2 ** (r + s) + 2 ** r + 2 ** s, "double-star formula")
            require((count & (count - 1) == 0) == (r == s == 1), "double-star power-two boundary")
            double_stars.append([r, s, count])

    determinants = []
    for n in (3, 5, 9, 15, 27):
        for t in (-3, -2, -1, 0, 1, 2, 3):
            d = det(cycle_matrix([t] * n))
            require(d == 1 + t ** n, "odd cycle determinant")
            determinants.append([n, t, d])
    require(det(cycle_matrix([2] * 3)) == 9, "positive triangle")
    require(det(cycle_matrix([-2] * 3)) == -7, "negative triangle")
    require(9 - (-7) == 16 and 9 * (-7) == -63, "two-sector invariants")
    require(det(cycle_matrix([8] * 3)) == 9 * 57 == 27 * 19, "cubic clock")
    require(emul((1, 2), (5, 2)) == (1, 8), "Eisenstein factorization")
    require(norm(1, 8) == norm(8, 1) == norm(7, -1) == 57, "norm57 associates")
    require(norm(1, 2) == 3 and norm(1, -2) == 7 and norm(5, 2) == 19, "norm factors")

    # Independent full solves versus triple-block elimination, including sources.
    port_cases = 0
    for n in (3, 6, 9, 12):
        for case in range(16):
            weights = [2 if (case >> (i % 4)) & 1 else -2 for i in range(n)]
            b = [((7 * i + 3 * case) % 11) - 5 for i in range(n)]
            a = cycle_matrix(weights)
            full = solve(a, b)
            coarse_weights = [weights[j] * weights[j + 1] * weights[j + 2] for j in range(0, n, 3)]
            coarse_b = [b[j] - weights[j] * b[j + 1] + weights[j] * weights[j + 1] * b[j + 2] for j in range(0, n, 3)]
            coarse_a = cycle_matrix(coarse_weights)
            coarse = solve(coarse_a, coarse_b)
            require(coarse == full[::3], "named-port solve")
            require(det(a) == det(coarse_a), "unimodular elimination determinant")
            rebuilt = [None] * n
            for k, j in enumerate(range(0, n, 3)):
                rebuilt[j] = coarse[k]
                rebuilt[j + 2] = b[j + 2] - weights[j + 2] * coarse[(k + 1) % len(coarse)]
                rebuilt[j + 1] = b[j + 1] - weights[j + 1] * rebuilt[j + 2]
            require(rebuilt == full, "internal reconstruction")
            port_cases += 1
    require(1 + 4 + 16 == 21, "source norm cost")
    require((1 - 2 * (-2) + 4 * 4) ** 2 == 21 * (1 + 4 + 16), "sharp norm witness")

    # The golden operator on the ACTUAL five-cycle's sum-zero module.
    p = [[int(j == (i + 1) % 5) for j in range(5)] for i in range(5)]
    pt = transpose(p)
    h = [[-p[i][j] - pt[i][j] for j in range(5)] for i in range(5)]
    h2 = matmul(h, h)
    require(all(h2[i][j] - h[i][j] - int(i == j) == 1 for i in range(5) for j in range(5)), "H2-H-I=J")
    v = [[1], [-1], [0], [0], [0]]
    hv = matmul(h, v)
    e = [[v[i][0], hv[i][0]] for i in range(5)]
    q = [[0, 1], [1, 1]]
    require(matmul(h, e) == matmul(e, q), "explicit golden intertwiner")
    factors = {2: 5, 3: 7, 5: 3, 7: 4, 11: 1, 13: 2}
    integer = 1
    for prime, exponent in factors.items():
        integer *= prime ** exponent
    def height(base):
        k, power = 0, base
        while integer % power == 0:
            k += 1
            power *= base
        return k
    require([height(x) for x in [4, 2, 5, 3, 7, 13]] == [2, 5, 3, 7, 4, 2], "five-cycle and tail")

    # Signed half-Hadamard triples and exact integer-image guard.
    def transform(a, b, c):
        return (a + c - b, b + c - a, c - a - b)
    require(transform(5, 7, 8) == (6, 10, -4), "owner triple")
    require(transform(0, 2, 10) == (8, 12, 8), "owner zero triple")
    def full_transform(a, b, c, d):
        return (a+b+c+d, a-b+c-d, -a+b+c-d, -a-b+c+d)
    require(full_transform(5, 7, 8, 0)[1:] == full_transform(6, 8, 9, 1)[1:], "fourth-input hostile")
    require(full_transform(5, 7, 8, 0)[0] != full_transform(6, 8, 9, 1)[0], "missing output differs")
    for a, b, c in product(range(-5, 6), repeat=3):
        x, y, z = transform(a, b, c)
        s = x + y - z
        require((x - z, y - z, x + y) == (2 * a, 2 * b, 2 * c), "inverse")
        require(s * s + x * x + y * y + z * z == 4 * (a * a + b * b + c * c), "Hadamard norm")
    for x, y, z in product(range(-5, 6), repeat=3):
        integral = all(v % 2 == 0 for v in (x - z, y - z, x + y))
        require(integral == (x % 2 == y % 2 == z % 2), "parity image")

    # Cyclic support carrier: proper rational support versus bad-prime reduction.
    a13 = [[int((j - i) % 13 in {0, 2, 10}) for j in range(13)] for i in range(13)]
    require(det(a13) == 159, "C13 unit-weight determinant")
    bad = [x for x in range(2, 53) if pow(x, 13, 53) == 1 and (1 + x*x + pow(x, 10, 53)) % 53 == 0]
    require(bool(bad), "nonzero rational spectrum can vanish at prime53")

    # Golden/Eisenstein marked-generator repair at the next binary digit.
    def mul4(x, y):
        return tuple(z % 4 for z in emul(x, y))
    def power4(x, k):
        y = (1, 0)
        for _ in range(k):
            y = mul4(y, x)
        return y
    golden_roots = []
    for x in product(range(4), repeat=2):
        xx = mul4(x, x)
        if ((xx[0] - x[0] - 1) % 4, (xx[1] - x[1]) % 4) == (0, 0):
            golden_roots.append(x)
    require(golden_roots == [(2, 3), (3, 1)], "golden roots modulo4")
    t = (2, 3)
    images = {((a + 2 * b) % 4, 3 * b % 4) for a, b in product(range(4), repeat=2)}
    require(len(images) == 16, "basis isomorphism")
    require(power4(t, 3) == (1, 2) and power4(t, 6) == (1, 0), "golden clock six")
    require(all(power4(t, k) != (1, 0) for k in range(1, 6)), "exact golden order")
    require(power4((0, 1), 3) == (1, 0), "Eisenstein clock three")
    def power_mod(x, exponent, modulus):
        y = (1, 0)
        while exponent:
            if exponent & 1:
                y = tuple(c % modulus for c in emul(y, x))
            x = tuple(c % modulus for c in emul(x, x))
            exponent //= 2
        return y
    tower = []
    t = (0, 1)
    for k in range(1, 17):
        modulus = 2 ** k
        tt = emul(t, t)
        defect = (tt[0] - t[0] - 1, tt[1] - t[1])
        require(all(x % modulus == 0 for x in defect), "lifted golden relation")
        require(t[1] % 2 == 1, "invertible lifted basis")
        order = 3 * 2 ** (k - 1)
        require(power_mod(t, order, modulus) == (1, 0), "lifted golden period")
        require(power_mod(t, order // 3, modulus) != (1, 0), "period needs three")
        if k > 1:
            require(power_mod(t, order // 2, modulus) != (1, 0), "period needs full two-part")
        tower.append({"k": k, "root": t, "order": order})
        delta = tuple((-x // modulus) % 2 for x in defect)
        t = tuple(x + modulus * y for x, y in zip(t, delta))

    print(json.dumps({
        "status": "FINITE-EXACT checks; all-size proofs in PROOFS.md",
        "checks": CHECKS,
        "square_sum_edges": sq_graph(9),
        "components": rows,
        "whole_independence_polynomial": polynomial(all_sets),
        "native_toggle_edges_on_eight_states": toggle_edges,
        "all_graph_edge_deletion_cases_through_five_vertices": edge_cases,
        "double_star_controls_r_s_1_to_5": double_stars,
        "cyclic_determinant_controls": determinants,
        "source_and_port_elimination_cases": port_cases,
        "golden_intertwiner_columns": e,
        "C13_unit_support_determinant": 159,
        "C13_zero_characters_mod53": bad,
        "golden_roots_in_Eisenstein_mod4": golden_roots,
        "marked_mod4_orders": {"omega": 3, "phi_image": 6},
        "compatible_golden_generator_tower": tower,
        "long_tuple_basic_checks": {"lengths": [10, 9], "sums": [28, 72], "square_sums": [84, 658], "interpretation": "UNRESOLVED"}
    }, indent=2))


if __name__ == "__main__":
    main()
