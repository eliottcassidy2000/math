#!/usr/bin/env python3
"""Exact consequence audit for THM-3333.

Universe
--------
* primitive first-octant Gaussian rays (m,n), 0 <= n <= m <= 48;
* every unordered pair of those rays;
* every Euclid parameter pair 1 <= n < m <= 30;
* explicit graph and clique-pair controls for five small right triangles.

All truth-bearing arithmetic is integral.  The checks use explicit exceptions
so that ``python`` and ``python -O`` execute the same validation path.
"""

from itertools import combinations
from math import gcd


RAY_BOUND = 48
TRIPLE_BOUND = 30
SIGN_REPRESENTATIVES = (
    (1, 1, 1),
    (-1, 1, 1),
    (1, -1, 1),
    (1, 1, -1),
)
OWNER_CONTROLS = ((1, 0), (5, 3), (-4, 7), (12, -13))
GL2_CONTROLS = (
    ((1, 1), (0, 1)),
    ((1, 0), (1, 1)),
    ((0, 1), (1, 0)),
    ((0, -1), (1, 1)),
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def triangular(n):
    return n * (n + 1) // 2


def polygonal(sides, n):
    return ((sides - 2) * n * n - (sides - 4) * n) // 2


def det2(u, v):
    return u[0] * v[1] - u[1] * v[0]


def dot2(u, v):
    return u[0] * v[0] + u[1] * v[1]


def gaussian_mul(z, w):
    return (z[0] * w[0] - z[1] * w[1],
            z[0] * w[1] + z[1] * w[0])


def phi(u):
    m, n = u
    return (m * m - n * n, 2 * m * n, m * m + n * n)


def phi_independent(u):
    real, imag = gaussian_mul(u, u)
    norm = dot2(u, u)
    return (real, imag, norm)


def symmetric_matrix(x):
    return ((x[2] + x[0], x[1]), (x[1], x[2] - x[0]))


def matrix_vector(matrix, u):
    return (
        matrix[0][0] * u[0] + matrix[0][1] * u[1],
        matrix[1][0] * u[0] + matrix[1][1] * u[1],
    )


def congruence2(matrix, symmetric):
    return tuple(
        tuple(
            sum(
                matrix[i][a] * symmetric[a][b] * matrix[j][b]
                for a in range(2)
                for b in range(2)
            )
            for j in range(2)
        )
        for i in range(2)
    )


def rank_one_matrix(u):
    return (
        (2 * u[0] * u[0], 2 * u[0] * u[1]),
        (2 * u[0] * u[1], 2 * u[1] * u[1]),
    )


def lorentz(x, y):
    return x[2] * y[2] - x[0] * y[0] - x[1] * y[1]


def owner_defect(owner, u):
    return dot2(u, u) - 91 * dot2(owner, u)


def tri_polar(x, y):
    return triangular(x + y) - triangular(x) - triangular(y)


def triangular_lorentz(x, y):
    return (tri_polar(x[2], y[2])
            - tri_polar(x[0], y[0])
            - tri_polar(x[1], y[1]))


def det3_columns(x, y, z):
    return (
        x[0] * (y[1] * z[2] - y[2] * z[1])
        - y[0] * (x[1] * z[2] - x[2] * z[1])
        + z[0] * (x[1] * y[2] - x[2] * y[1])
    )


def signed_radius_potential(x, signs=(1, 1, 1)):
    return (
        triangular(signs[0] * x[0])
        + triangular(signs[1] * x[1])
        - triangular(signs[2] * x[2])
    )


def polygonal_radius_potential(x, sides):
    return polygonal(sides, x[0]) + polygonal(sides, x[1]) \
        - polygonal(sides, x[2])


def content3(x):
    return gcd(gcd(abs(x[0]), abs(x[1])), abs(x[2]))


def primitive_normalize3(x):
    common = content3(x)
    return tuple(value // common for value in x)


def add2(u, v):
    return (u[0] + v[0], u[1] + v[1])


def sub2(u, v):
    return (u[0] - v[0], u[1] - v[1])


def add3(x, y):
    return tuple(a + b for a, b in zip(x, y))


def scale3(k, x):
    return tuple(k * a for a in x)


def complete_edges(vertices):
    ordered = sorted(vertices)
    return {pair for pair in combinations(ordered, 2)}


def cycle_rank(vertices, edges):
    parent = {v: v for v in vertices}

    def find(v):
        while parent[v] != v:
            parent[v] = parent[parent[v]]
            v = parent[v]
        return v

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for a, b in edges:
        union(a, b)
    components = len({find(v) for v in vertices})
    return len(edges) - len(vertices) + components


def gf2_rank(columns):
    pivots = {}
    for value in columns:
        x = value
        while x:
            pivot = x.bit_length() - 1
            if pivot in pivots:
                x ^= pivots[pivot]
            else:
                pivots[pivot] = x
                break
    return len(pivots)


def clique_pair_relative_h1_rank(vertex_count, y_edges, missing_edges):
    """Compute H_1(Cl(K_V), Cl(Y); F_2) from the relative boundary."""
    missing_index = {edge: i for i, edge in enumerate(sorted(missing_edges))}
    columns = []
    for triple in combinations(range(vertex_count), 3):
        boundary = list(combinations(triple, 2))
        if all(edge in y_edges for edge in boundary):
            continue
        vector = 0
        for edge in boundary:
            index = missing_index.get(edge)
            if index is not None:
                vector ^= 1 << index
        columns.append(vector)
    # C_0(X,Y)=0, so every relative 1-chain is a cycle.
    return len(missing_edges) - gf2_rank(columns)


def pythagorean_graph_invoice(m, n, check_clique=False):
    a = m * m - n * n
    b = 2 * m * n
    c = m * m + n * n
    r = n * (m - n)
    require(a > 0 and b > 0 and c > 0, (m, n, "nonpositive triple"))
    require(a * a + b * b == c * c, (m, n, "not Pythagorean"))
    require(
        triangular(a) + triangular(b) - triangular(c) == r,
        (m, n, "triangular inradius"),
    )
    require((a + b - c) // 2 == r and a + b - c == 2 * r,
            (m, n, "boundary inradius"))
    require(triangular(m) - triangular(n) - triangular(m - n) == r,
            (m, n, "parameter polarization"))

    r_a = m * (m - n)
    r_b = n * (m + n)
    r_c = m * (m + n)
    signed_values = tuple(
        signed_radius_potential((a, b, c), signs)
        for signs in SIGN_REPRESENTATIVES
    )
    require(signed_values == (r, -r_a, -r_b, r_c),
            (m, n, "signed in/exradius potentials", signed_values))
    for sides in range(3, 13):
        require(
            polygonal_radius_potential((a, b, c), sides)
            == (4 - sides) * r,
            (m, n, sides, "polygonal radius potential"),
        )

    edge_a = triangular(a)
    edge_b = triangular(b)
    edge_c = triangular(c)
    beta_a = triangular(a - 1)
    beta_b = triangular(b - 1)
    beta_c = triangular(c - 1)
    require(edge_a + edge_b - edge_c == r,
            (m, n, "edge inradius"))
    require(beta_c - beta_a - beta_b == r,
            (m, n, "cycle inradius"))

    overlap = a + b - c + 1
    require(overlap == 2 * r + 1, (m, n, "overlap"))
    left_only = c - b
    right_only = c - a
    require(left_only == (m - n) ** 2, (m, n, "left-only square"))
    require(right_only == 2 * n * n, (m, n, "right-only double square"))
    missing_count = left_only * right_only
    require(missing_count == 2 * r * r, (m, n, "relative rank formula"))

    parameter_clique_relative_h1 = None
    clique_relative_h1 = None
    if check_clique:
        parameter_vertices = set(range(m + 1))
        parameter_left = set(range(n + 1))
        parameter_right = set(range(n, m + 1))
        require(parameter_left | parameter_right == parameter_vertices,
                (m, n, "parameter cover union"))
        require(parameter_left & parameter_right == {n},
                (m, n, "parameter cover point overlap"))
        parameter_x_edges = complete_edges(parameter_vertices)
        parameter_y_edges = (
            complete_edges(parameter_left) | complete_edges(parameter_right)
        )
        parameter_missing = parameter_x_edges - parameter_y_edges
        require(len(parameter_missing) == r,
                (m, n, "parameter missing biclique"))
        require(
            cycle_rank(parameter_vertices, parameter_x_edges)
            - cycle_rank(parameter_vertices, parameter_y_edges) == r,
            (m, n, "parameter graph relative H1"),
        )
        parameter_clique_relative_h1 = clique_pair_relative_h1_rank(
            m + 1, parameter_y_edges, parameter_missing
        )
        require(parameter_clique_relative_h1 == 0,
                (m, n, "parameter clique filling"))

        vertices = set(range(c + 1))
        set_a = set(range(a + 1))
        set_b = set(range(c - b, c + 1))
        require(set_a | set_b == vertices, (m, n, "cover union"))
        require(len(set_a & set_b) == overlap, (m, n, "cover overlap"))
        x_edges = complete_edges(vertices)
        y_edges = complete_edges(set_a) | complete_edges(set_b)
        missing_edges = x_edges - y_edges
        require(len(missing_edges) == missing_count,
                (m, n, "explicit missing biclique"))
        require(cycle_rank(vertices, x_edges) == beta_c,
                (m, n, "explicit beta X"))
        expected_y_beta = beta_a + beta_b - triangular(2 * r - 1)
        require(cycle_rank(vertices, y_edges) == expected_y_beta,
                (m, n, "explicit beta Y"))
        require(cycle_rank(vertices, x_edges) - cycle_rank(vertices, y_edges)
                == missing_count, (m, n, "graph relative H1"))
        clique_relative_h1 = clique_pair_relative_h1_rank(
            c + 1, y_edges, missing_edges
        )
        require(clique_relative_h1 == 0,
                (m, n, "clique filling should kill relative H1"))

    return {
        "parameters": (m, n),
        "triple": (a, b, c),
        "r": r,
        "overlap": overlap,
        "parameter_missing": r,
        "parameter_clique_relative_h1": parameter_clique_relative_h1,
        "missing": missing_count,
        "clique_relative_h1": clique_relative_h1,
    }


def main():
    rays = [
        (m, n)
        for m in range(1, RAY_BOUND + 1)
        for n in range(0, m + 1)
        if gcd(m, n) == 1
    ]

    image_to_ray = {}
    odd_odd = 0
    for u in rays:
        x = phi(u)
        require(x == phi_independent(u), (u, "Gaussian square mismatch"))
        require(symmetric_matrix(x) == rank_one_matrix(u),
                (u, "rank-one symmetric-square mismatch"))
        for matrix in GL2_CONTROLS:
            transformed = matrix_vector(matrix, u)
            require(
                symmetric_matrix(phi(transformed))
                == rank_one_matrix(transformed)
                == congruence2(matrix, symmetric_matrix(x)),
                (u, matrix, "symmetric-square covariance"),
            )
        require(x[0] * x[0] + x[1] * x[1] == x[2] * x[2],
                (u, "light-cone mismatch"))
        common = content3(x)
        expected_common = 2 if (u[0] % 2 == 1 and u[1] % 2 == 1) else 1
        require(common == expected_common, (u, common, expected_common))
        odd_odd += expected_common == 2
        require(gcd(u[0], x[2]) == 1 and gcd(u[1], x[2]) == 1,
                (u, "Kelvin denominator not reduced"))
        require(x[1] * u[0] == u[1] * (x[2] + x[0]),
                (u, "inverse slope identity"))
        require(x not in image_to_ray, (u, image_to_ray.get(x), "Phi collision"))
        image_to_ray[x] = u

    pair_count = 0
    farey_edges = 0
    farey_faces = 0
    fixed_owner_updates = 0
    for u, v in combinations(rays, 2):
        pair_count += 1
        x, y = phi(u), phi(v)
        delta = det2(u, v)
        gram_cross = dot2(u, v)
        require(x[2] * y[2] == gram_cross * gram_cross + delta * delta,
                (u, v, "two-square product identity"))
        require(lorentz(x, y) == 2 * delta * delta,
                (u, v, "Lorentz determinant identity"))
        require(triangular_lorentz(x, y) == 2 * delta * delta,
                (u, v, "triangular polarization identity"))
        for signs in SIGN_REPRESENTATIVES:
            require(
                signed_radius_potential(x, signs)
                + signed_radius_potential(y, signs)
                - signed_radius_potential(add3(x, y), signs)
                == 2 * delta * delta,
                (u, v, signs, "signed triangular polarization"),
            )
        for sides in range(3, 9):
            require(
                polygonal_radius_potential(x, sides)
                + polygonal_radius_potential(y, sides)
                - polygonal_radius_potential(add3(x, y), sides)
                == 2 * (sides - 2) * delta * delta,
                (u, v, sides, "polygonal Lorentz polarization"),
            )
        require((abs(delta) == 1) == (lorentz(x, y) == 2),
                (u, v, "Farey edge equivalence"))
        if abs(delta) != 1:
            continue

        farey_edges += 1
        plus, minus = add2(u, v), sub2(u, v)
        xp, xm = phi(plus), phi(minus)
        diamond_left = add3(xp, xm)
        diamond_right = scale3(2, add3(x, y))
        require(diamond_left == diamond_right, (u, v, "Farey diamond"))
        require(lorentz(x, xp) == 2 and lorentz(y, xp) == 2,
                (u, v, "mediant face adjacency"))
        require(det3_columns(x, y, xp) == -4 * delta ** 3,
                (u, v, "oriented face determinant"))
        require(xp[2] - xm[2] == 4 * gram_cross,
                (u, v, "two-face Gram recovery"))
        require(xp[2] + xm[2] == 2 * (x[2] + y[2]),
                (u, v, "Vieta face sum"))
        require(xp[2] * xm[2] == (x[2] - y[2]) ** 2 + 4,
                (u, v, "Vieta face product"))
        for face_norm in (xp[2], xm[2]):
            require(
                x[2] * x[2] + y[2] * y[2] + face_norm * face_norm
                - 2 * x[2] * y[2] - 2 * x[2] * face_norm
                - 2 * y[2] * face_norm == -4,
                (u, v, face_norm, "Vieta face surface"),
            )
        if gram_cross >= 0:
            require(gram_cross * gram_cross == x[2] * y[2] - 1,
                    (u, v, "acute Pell edge"))
            require(xp[2] == x[2] + y[2] + 2 * gram_cross,
                    (u, v, "mediant norm recurrence"))
            for owner in OWNER_CONTROLS:
                f_u = owner_defect(owner, u)
                f_v = owner_defect(owner, v)
                f_plus = owner_defect(owner, plus)
                require(x[2] - f_u == 91 * dot2(owner, u)
                        and y[2] - f_v == 91 * dot2(owner, v),
                        (owner, u, v, "owner-coordinate recovery"))
                require(f_plus == f_u + f_v + 2 * gram_cross,
                        (owner, u, v, "fixed-owner scalar recurrence"))
                fixed_owner_updates += 1
        face_images = (x, y, xp)
        require(sorted(content3(point) for point in face_images) == [1, 1, 2],
                (u, v, "Farey face parity contents"))
        normalized_face = tuple(primitive_normalize3(point)
                                for point in face_images)
        normalized_weights = sorted(
            lorentz(left, right)
            for left, right in combinations(normalized_face, 2)
        )
        require(normalized_weights == [1, 1, 2],
                (u, v, "normalized Farey face weights"))
        farey_faces += 1

    # The THM-2056 gate becomes a constant Lorentz cap after Kelvin-square
    # rescaling.  This is checked without floating point by clearing squares.
    landmarks = ((1, 0), (0, 1), (1, 1))
    gate_probe_rays = ((1, 0), (90, 1), (91, 1), (92, 1), (144, 55))
    kelvin_gate_passes = 0
    for u in gate_probe_rays:
        q_u = phi(u)[2]
        determinant_norm = max(abs(det2(u, column)) for column in landmarks)
        max_lorentz = max(lorentz(phi(u), phi(column))
                          for column in landmarks)
        require(max_lorentz == 2 * determinant_norm * determinant_norm,
                (u, "Kelvin landmark pairing"))
        original_gate = 91 * determinant_norm <= q_u
        lorentz_cap = 91 * 91 * max_lorentz <= 2 * q_u * q_u
        require(original_gate == lorentz_cap,
                (u, "Kelvin Lorentz-cap equivalence"))
        kelvin_gate_passes += original_gate
    require(0 < kelvin_gate_passes < len(gate_probe_rays),
            "Kelvin gate controls need positive and hostile cases")

    triple_rows = 0
    for m in range(2, TRIPLE_BOUND + 1):
        for n in range(1, m):
            pythagorean_graph_invoice(m, n)
            triple_rows += 1

    graph_controls = []
    for parameters in [(2, 1), (3, 2), (4, 1), (4, 3), (5, 2)]:
        graph_controls.append(pythagorean_graph_invoice(*parameters,
                                                        check_clique=True))

    # Hostile 1: one sum-of-two-squares value can have differently adjacent
    # representations.  Hypotenuse values alone are not a Farey graph.
    u, u_alt, v = (8, 1), (7, 4), (1, 0)
    require(phi(u)[2] == phi(u_alt)[2] == 65, "norm-65 hostile setup")
    require(abs(det2(u, v)) == 1 and abs(det2(u_alt, v)) == 4,
            "norm-65 hostile determinants")
    require(lorentz(phi(u), phi(v)) == 2,
            "norm-65 positive control")
    require(lorentz(phi(u_alt), phi(v)) == 32,
            "norm-65 hostile control")

    # Hostile 2: primitive-normalizing an odd/odd spinor divides its light-cone
    # point by two and changes the uniform Farey pairing from 2 to 1.
    odd_spinor, endpoint = (3, 1), (2, 1)
    raw = phi(odd_spinor)
    normalized = tuple(value // 2 for value in raw)
    require(raw == (8, 6, 10) and normalized == (4, 3, 5),
            "parity hostile setup")
    require(lorentz(raw, phi(endpoint)) == 2,
            "raw parity positive control")
    require(lorentz(normalized, phi(endpoint)) == 1,
            "primitive-normalization hostile")
    require(sorted(normalized[:2]) == sorted(phi(endpoint)[:2])
            and normalized[2] == phi(endpoint)[2],
            "unordered primitive-triple edge collapse")

    # Hostile 3: the scalar edge pairing squares determinant orientation.
    orient_u, orient_v = (2, 1), (1, 1)
    require(det2(orient_u, orient_v) == 1,
            "orientation hostile setup")
    require(lorentz(phi(orient_u), phi(orient_v))
            == lorentz(phi(orient_v), phi(orient_u)) == 2,
            "orientation scalar hostile")
    require(
        det3_columns(phi(orient_u), phi(orient_v),
                     phi(add2(orient_u, orient_v))) == -4,
        "orientation face positive control",
    )
    require(
        det3_columns(phi(orient_v), phi(orient_u),
                     phi(add2(orient_u, orient_v))) == 4,
        "orientation face reversal control",
    )

    # Hostile 4: determinant two is the first non-Farey quantized shell.
    require(lorentz(phi((1, 0)), phi((1, 2))) == 8,
            "determinant-two hostile")

    # Hostile 5: coordinate addition in the polarization is not spinor
    # mediant addition, even on a Farey edge.
    ambient_u, ambient_v = (2, 1), (3, 2)
    require(abs(det2(ambient_u, ambient_v)) == 1,
            "ambient-addition hostile setup")
    require(add3(phi(ambient_u), phi(ambient_v)) == (8, 16, 18),
            "ambient-addition positive control")
    require(phi(add2(ambient_u, ambient_v)) == (16, 30, 34),
            "spinor-mediant positive control")
    require(add3(phi(ambient_u), phi(ambient_v))
            != phi(add2(ambient_u, ambient_v)),
            "ambient addition must not equal spinor mediant")

    # Hostile 6: leg order and raw-spinor membership are load-bearing.
    triple_x, triple_y = (3, 4, 5), (5, 12, 13)
    require(lorentz(triple_x, triple_y) == 2,
            "ordered-leg positive control")
    require(lorentz((4, 3, 5), triple_y) == 9,
            "leg-swap hostile")
    require(lorentz((4, 3, 5), (6, 8, 10)) == 2,
            "arbitrary-null pairing-two hostile")
    require(phi((2, 1)) == triple_x and phi((4, 1)) == (15, 8, 17),
            "Berggren/Farey hostile setup")
    require(abs(det2((2, 1), (4, 1))) == 2
            and lorentz(phi((2, 1)), phi((4, 1))) == 8,
            "Berggren child need not be a Farey edge")

    print("THM-3333 exact consequence audit")
    print(f"ray_bound={RAY_BOUND}")
    print(f"primitive_projective_rays={len(rays)}")
    print(f"odd_odd_parity_vertices={odd_odd}")
    print(f"unordered_ray_pairs={pair_count}")
    print(f"farey_edges={farey_edges}")
    print(f"oriented_mediant_faces_checked={farey_faces}")
    print("farey_face_raw_contents=1,1,2; normalized_edge_weights=1,1,2")
    print(f"fixed_owner_scalar_updates={fixed_owner_updates}")
    print(
        "kelvin_gate_controls={} pass={} fail={}".format(
            len(gate_probe_rays), kelvin_gate_passes,
            len(gate_probe_rays) - kelvin_gate_passes,
        )
    )
    print(f"euclid_parameter_rows={triple_rows}")
    print("graph_clique_controls:")
    for row in graph_controls:
        print(
            "  params={} triple={} r={} parameter_graph_h1={} "
            "parameter_clique_h1={} overlap={} graph_relative_h1={} "
            "clique_relative_h1={}".format(
                row["parameters"], row["triple"], row["r"],
                row["parameter_missing"],
                row["parameter_clique_relative_h1"], row["overlap"],
                row["missing"],
                row["clique_relative_h1"]
            )
        )
    print("signed_in_exradius_and_polygonal_rows=ALL")
    print("hostile_norm_collision_65: Lorentz pairings 2 versus 32")
    print("hostile_primitive_normalization: Lorentz pairing 2 -> 1")
    print("hostile_unordered_triple: adjacent spinors -> same 3-4-5")
    print("hostile_edge_orientation: scalar 2 both ways; face det -4/+4")
    print("hostile_det_two_shell: Lorentz pairing 8")
    print("hostile_ambient_addition: (8,16,18) != (16,30,34)")
    print("hostile_leg_swap: Lorentz pairing 2 -> 9")
    print("hostile_arbitrary_nulls: pairing 2 without raw spinor edge")
    print("hostile_berggren_edge: determinant 2; Lorentz pairing 8")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
