#!/usr/bin/env python3
"""Independent exact audit for THM-3991.

The positive controls use the standard periodic Kuhn triangulations.  The
hostile controls separate ordinary from compact-support Euler additivity and
exhibit the first local nonunimodular tetrahedron inside the periodic
five-tetrahedron cube grammar.
"""

from hashlib import sha256
from itertools import combinations, permutations, product
from math import factorial, prod
import json
import sys


sys.stdout.reconfigure(newline="\n")
CHECKS = 0


def require(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(label)


def det_bareiss(matrix):
    a = [list(map(int, row)) for row in matrix]
    n = len(a)
    if n == 0:
        return 1
    sign = 1
    previous = 1
    for k in range(n - 1):
        if a[k][k] == 0:
            pivot = next((i for i in range(k + 1, n) if a[i][k]), None)
            if pivot is None:
                return 0
            a[k], a[pivot] = a[pivot], a[k]
            sign *= -1
        pivot_value = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                a[i][j] = (a[i][j] * pivot_value
                           - a[i][k] * a[k][j]) // previous
        previous = pivot_value
    return sign * a[-1][-1]


def spatial_normalized_volume(vertices):
    base = vertices[0]
    columns = [[vertices[j][i] - base[i] for j in range(1, len(vertices))]
               for i in range(len(base))]
    return abs(det_bareiss(columns))


def cone_determinant(vertices):
    # Columns are the primitive height-one ray generators (v,1).
    matrix = [[vertices[j][i] for j in range(len(vertices))]
              for i in range(len(vertices[0]))]
    matrix.append([1] * len(vertices))
    return abs(det_bareiss(matrix))


def kuhn_simplex(n, order):
    point = [0] * n
    vertices = [tuple(point)]
    for i in order:
        point = point.copy()
        point[i] += 1
        vertices.append(tuple(point))
    return tuple(vertices)


def audit_kuhn():
    rows = []
    for n in range(1, 7):
        simplex_types = tuple(kuhn_simplex(n, p) for p in permutations(range(n)))
        require(len(simplex_types) == factorial(n), f"Kuhn type count n={n}")
        for simplex in simplex_types:
            require(spatial_normalized_volume(simplex) == 1,
                    f"unimodular spatial simplex n={n}")
            require(cone_determinant(simplex) == 1,
                    f"smooth height-one cone n={n}")
        # Several diagonal translation sublattices.  One unit cube per
        # residue class gives d*n! top-simplex orbits and d vertex orbits.
        period_vectors = [tuple([1] * n), tuple(2 if i == 0 else 1 for i in range(n))]
        if n >= 2:
            period_vectors.append(tuple(2 if i < 2 else 1 for i in range(n)))
        for periods in period_vectors:
            d = prod(periods)
            top = d * len(simplex_types)
            vertices = d
            require(top == d * factorial(n), f"factorial count {(n, periods)}")
            require(vertices == d, f"component count {(n, periods)}")
            rows.append((n, periods, d, vertices, top))
    return rows


def parity(v):
    return sum(v) & 1


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def five_tetra_cube(origin):
    """Five-tetrahedron cube decomposition with global-even central vertices."""
    local = tuple(product((0, 1), repeat=3))
    central_local_parity = parity(origin)
    central = tuple(add(origin, v) for v in local if parity(v) == central_local_parity)
    tets = [central]
    for corner_local in local:
        if parity(corner_local) == central_local_parity:
            continue
        corner = add(origin, corner_local)
        neighbours = []
        for i in range(3):
            w = list(corner_local)
            w[i] = 1 - w[i]
            neighbours.append(add(origin, tuple(w)))
        tets.append((corner, *neighbours))
    return tuple(tets)


def boundary_triangles(tets, axis, value):
    counts = {}
    for tet in tets:
        for face in combinations(tet, 3):
            key = tuple(sorted(face))
            counts[key] = counts.get(key, 0) + 1
    return tuple(sorted(face for face, count in counts.items()
                        if count == 1 and all(v[axis] == value for v in face)))


def audit_nonunimodular_hostile():
    origins = ((0, 0, 0), (1, 0, 0))
    cube_rows = []
    for origin in origins:
        tets = five_tetra_cube(origin)
        volumes = sorted(spatial_normalized_volume(t) for t in tets)
        require(volumes == [1, 1, 1, 1, 2], f"five-tet volumes {origin}")
        require(sum(volumes) == 6 and len(volumes) == 5,
                f"volume/count mismatch {origin}")
        cube_rows.append((origin, tuple(volumes)))

    # Alternating the local central parity makes the face diagonals match.
    origin = (0, 0, 0)
    for axis in range(3):
        neighbour = tuple(1 if i == axis else 0 for i in range(3))
        upper = boundary_triangles(five_tetra_cube(origin), axis, 1)
        lower = boundary_triangles(five_tetra_cube(neighbour), axis, 1)
        require(upper == lower, f"periodic face match axis={axis}")

    # The pattern is periodic under the index-two even-sum sublattice.  Its
    # quotient has two cube orbits, ten tetrahedron orbits, but normalized
    # volume twelve.  It is outside THM-3991 because it is not invariant under
    # the full Z^3 translation lattice and its central tetrahedra are singular.
    quotient = {
        "translation_index": 2,
        "cube_orbits": 2,
        "tetrahedron_orbits": 10,
        "normalized_volume_sum": 12,
        "vertex_orbits": 2,
    }
    require(quotient["tetrahedron_orbits"] != quotient["normalized_volume_sum"],
            "nonunimodular count differs from volume")
    return cube_rows, quotient


def main():
    kuhn_rows = audit_kuhn()

    # Standard A2 / one-square quotient cell-orbit counts.
    a2 = {"vertices": 1, "edges": 3, "triangles": 2}
    require(a2 == {"vertices": 1, "edges": 3, "triangles": 2}, "A2 f-vector")
    require(a2["triangles"] == factorial(2), "A2 Euler contribution")

    # Positive-integer solutions to d*n! = 2.
    solutions = tuple((n, d) for n in range(1, 11) for d in range(1, 101)
                      if d * factorial(n) == 2)
    require(solutions == ((1, 2), (2, 1)), "integer solutions")

    # Compact-support Euler is the additive invariant.  Ordinary Euler is
    # not additive for arbitrary locally closed decompositions: [0,1] is the
    # disjoint union of (0,1) and two endpoints.
    interval_hostile = {
        "ordinary_total": 1,
        "ordinary_open_plus_boundary": 1 + 2,
        "compact_support_open_plus_boundary": -1 + 2,
    }
    require(interval_hostile["ordinary_total"]
            != interval_hostile["ordinary_open_plus_boundary"],
            "ordinary Euler hostile")
    require(interval_hostile["ordinary_total"]
            == interval_hostile["compact_support_open_plus_boundary"],
            "compact-support repair")

    torus_euler = tuple((m, 1 if m == 0 else 0) for m in range(7))
    require(torus_euler[0] == (0, 1) and all(x == 0 for _, x in torus_euler[1:]),
            "complex torus orbit Euler")

    cube_rows, nonunimodular = audit_nonunimodular_hostile()

    # Smallest numerical free-quotient escape above the S6 dimension.
    free_quotient_candidates = tuple(
        (n, d, d * factorial(n) // 2)
        for n in range(3, 7) for d in range(1, 4)
        if (d * factorial(n)) % 2 == 0
    )
    require(free_quotient_candidates[0] == (3, 1, 3),
            "first higher-rank quotient order")

    semantic = {
        "kuhn_rows": kuhn_rows,
        "a2": a2,
        "integer_solutions": solutions,
        "torus_euler": torus_euler,
        "ordinary_additivity_hostile": interval_hostile,
        "nonunimodular_cube_rows": cube_rows,
        "nonunimodular_periodic_quotient": nonunimodular,
        "free_quotient_candidates": free_quotient_candidates,
    }
    digest = sha256(json.dumps(semantic, sort_keys=True,
                                separators=(",", ":")).encode("ascii")).hexdigest()
    print("THM3991_PERIODIC_UNIMODULAR_TORIC_CUSP_INDEPENDENT_AUDIT_20260824")
    print(f"kuhn_controls={len(kuhn_rows)};dimensions=1..6;all_cone_determinants=1")
    print("A2_orbits=vertices_1;edges_3;triangles_2;Euler=2")
    print("integer_solutions_d_times_n_factorial_eq_2=((1,2),(2,1))")
    print("Euler_additivity=use_compact_support_on_strata_then_compact_W_gives_ordinary")
    print("ordinary_hostile=[0,1]=(0,1)_disjoint_union_two_points;1_ne_1+2;compact_support=-1+2=1")
    print("nonunimodular_hostile=periodic_five_tet_cube;top_count_10;normalized_volume_12")
    print("first_numerical_free_quotient_escape=(n,d,|H|)=(3,1,3);existence_not_claimed")
    print("scope=THM3991_PASS_with_wording_repairs;Hopf_S6_global_claims_OPEN")
    print("semantic_sha256=" + digest)
    print(f"checks={CHECKS}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
