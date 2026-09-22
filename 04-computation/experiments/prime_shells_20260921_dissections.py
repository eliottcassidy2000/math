"""Exact, optimized-safe controls for the dissection/finite-conjugacy note.

Normal invocation writes the sibling JSON. No external packages are required.
The all-dimension and irrationality statements are proved in the note; the
declared finite universes here audit their coordinates and boundaries.
"""
from fractions import Fraction as F
from itertools import combinations, permutations, product
from math import factorial
from pathlib import Path
import hashlib
import json


CHECKS = 0


def require(test, label):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(label)


def solve(rows, rhs):
    n = len(rows)
    a = [[F(x) for x in row] + [F(y)] for row, y in zip(rows, rhs)]
    for j in range(n):
        pivot = next((i for i in range(j, n) if a[i][j]), None)
        if pivot is None:
            return None
        a[j], a[pivot] = a[pivot], a[j]
        z = a[j][j]
        a[j] = [x / z for x in a[j]]
        for i in range(n):
            if i != j:
                z = a[i][j]
                a[i] = [x - z * y for x, y in zip(a[i], a[j])]
    return tuple(a[i][-1] for i in range(n))


def determinant(rows):
    n = len(rows)
    a = [[F(x) for x in row] for row in rows]
    result = F(1)
    for j in range(n):
        pivot = next((i for i in range(j, n) if a[i][j]), None)
        if pivot is None:
            return F(0)
        if pivot != j:
            a[j], a[pivot] = a[pivot], a[j]
            result = -result
        z = a[j][j]
        result *= z
        for i in range(j + 1, n):
            scale = a[i][j] / z
            for k in range(j + 1, n):
                a[i][k] -= scale * a[j][k]
    return result


def volume(vertices):
    n = len(vertices) - 1
    return abs(determinant([[v[j] - vertices[0][j] for j in range(n)]
                            for v in vertices[1:]])) / factorial(n)


def squared_edges(vertices):
    return sorted(sum((x - y) ** 2 for x, y in zip(a, b))
                  for a, b in combinations(vertices, 2))


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def distance1(a, b):
    return sum(abs(x - y) for x, y in zip(a, b))


def barycentric_functions(vertices):
    # Each returned affine functional is one barycentric coordinate.
    rows = [tuple(v) + (1,) for v in vertices]
    return [solve(rows, [int(i == j) for j in range(len(vertices))])
            for i in range(len(vertices))]


def demicube_constraints(n):
    constraints = []
    for i in range(n):
        constraints.append((tuple(-int(i == j) for j in range(n)), 0))
        constraints.append((tuple(int(i == j) for j in range(n)), 1))
    for v in product((0, 1), repeat=n):
        if sum(v) % 2:
            # distance1(x,v)>=1, expressed as a.x<=b on the cube.
            constraints.append((tuple(2 * x - 1 for x in v), sum(v) - 1))
    return constraints


def check_demicube_vertices():
    counts = {}
    for n in range(2, 5):
        constraints = demicube_constraints(n)
        found = set()
        for active in combinations(constraints, n):
            p = solve([a for a, _ in active], [b for _, b in active])
            if p is not None and all(dot(a, p) <= b for a, b in constraints):
                found.add(p)
        expected = {v for v in product((0, 1), repeat=n) if sum(v) % 2 == 0}
        require(found == expected, f"demicube vertices n={n}")
        counts[str(n)] = len(found)
    return counts


def check_cube():
    vertices = tuple(product((0, 1), repeat=3))
    central = tuple(v for v in vertices if sum(v) % 2 == 0)
    odd = tuple(v for v in vertices if sum(v) % 2)
    corners = [tuple([v] + [u for u in central if distance1(u, v) == 1])
               for v in odd]
    cells = [central] + corners
    require([volume(t) for t in cells] == [F(1, 3)] + [F(1, 6)] * 4,
            "five-cell exact volumes")
    require(sum(volume(t) for t in cells) == 1, "volume sum")
    require(squared_edges(central) == [2] * 6, "regular central tetrahedron")
    require(all(squared_edges(t) == [1, 1, 1, 2, 2, 2] for t in corners),
            "corner metric")
    bary = [barycentric_functions(t) for t in cells]
    grid = 0
    for denominator in range(1, 13):
        for raw in product(range(denominator + 1), repeat=3):
            p = tuple(F(x, denominator) for x in raw)
            values = [[dot(row, p + (1,)) for row in fs] for fs in bary]
            membership = [all(x >= 0 for x in vals) for vals in values]
            inequalities = [all(distance1(p, v) >= 1 for v in odd)]
            inequalities += [distance1(p, v) <= 1 for v in odd]
            require(membership == inequalities, "barycentric/inequality agreement")
            require(any(membership), "cube covered")
            require(sum(all(x > 0 for x in vals) for vals in values) <= 1,
                    "disjoint cell interiors")
            grid += 1
    midpoint = tuple(F(x + y, 2) for x, y in zip(central[2], central[3]))
    true_half = (central[0], central[1], central[2], midpoint)
    require(volume(true_half) == F(1, 6), "actual regular-tetrahedron half volume")
    require(squared_edges(true_half) == [F(1, 2), F(3, 2), F(3, 2), 2, 2, 2],
            "actual half-cut edge spectrum")
    require(squared_edges(true_half) != squared_edges(corners[0]),
            "equal volume does not give congruence")
    center = (F(1, 2),) * 3
    octa_volume = F(0)
    for signs in product((-1, 1), repeat=3):
        axes = [tuple(center[j] + (F(signs[i], 2) if i == j else 0)
                      for j in range(3)) for i in range(3)]
        octa_volume += volume([center] + axes)
    require(octa_volume == F(1, 6), "dual-tetrahedron intersection octahedron")
    require(2 * F(1, 3) - octa_volume == F(1, 2), "dual-tetrahedron union")
    hostile = (F(1, 2), F(0), F(0))
    require(not all(dot(row, hostile + (1,)) >= 0 for row in bary[0]),
            "edge midpoint outside even tetrahedron")
    odd_bary = barycentric_functions(odd)
    require(not all(dot(row, hostile + (1,)) >= 0 for row in odd_bary),
            "edge midpoint outside odd tetrahedron")
    scaled = [tuple((v[0], 2 * v[1], 3 * v[2]) for v in t) for t in cells]
    require([volume(t) for t in scaled] == [2, 1, 1, 1, 1], "cuboid scaling")
    require(squared_edges(scaled[0]) == [5, 5, 10, 10, 13, 13],
            "cuboid destroys regularity")
    sheared = [tuple((v[0] + v[1], v[1], v[2]) for v in t) for t in cells]
    require(len({tuple(squared_edges(t)) for t in sheared[1:]}) > 1,
            "general affine map destroys corner congruence")
    return {"grid_points_with_repetition": grid,
            "normalized_cell_volumes": [2, 1, 1, 1, 1],
            "intersection_volume": str(octa_volume), "union_volume": "1/2"}


def check_rectangles():
    count = 0
    for w, h in product(range(1, 9), repeat=2):
        a, b, m = (F(0), F(0)), (F(w), F(0)), (F(w, 2), F(h))
        left_half = (a, (F(w, 2), 0), m)
        left_ear = (a, (0, F(h)), m)
        right_ear = (b, (F(w), F(h)), m)
        require(volume((a, b, m)) == F(w * h, 2), "central triangle area")
        require(volume(left_half) == volume(left_ear) == volume(right_ear)
                == F(w * h, 4), "ear areas")
        require(squared_edges(left_half) == squared_edges(left_ear)
                == squared_edges(right_ear), "ears are congruent half-copies")
        require((4 * h * h == 3 * w * w) is False,
                "no positive rational rectangle aspect is exactly equilateral")
        count += 1
    # In Q(sqrt(3)), (r+t sqrt(3))^2 has rational part r^2+3t^2.
    # These are finite witnesses; the note proves the all-rational obstruction.
    witness_count = 0
    for r, t in product(range(-12, 13), repeat=2):
        if r or t:
            require(r * r + 3 * t * t > 0, "quadratic-field square obstruction")
            witness_count += 1
    return {"rectangles": count, "quadratic_field_witnesses": witness_count}


def check_general_dimensions():
    rows = []
    for n in range(2, 11):
        central_fraction = 1 - F(2 ** (n - 1), factorial(n))
        corner_fraction = F(1, factorial(n))
        require(central_fraction + 2 ** (n - 1) * corner_fraction == 1,
                "all-dimension volume identity")
        if n >= 3:
            require((corner_fraction == central_fraction / 2) == (n == 3),
                    "half-volume boundary")
            require((2 ** (n - 1) == n + 1) == (n == 3),
                    "simplex vertex-count boundary")
        rows.append({"dimension": n, "parity_vertices": 2 ** (n - 1),
                     "central_volume": str(central_fraction),
                     "one_corner_volume": str(corner_fraction)})
    return rows


def check_finite_conjugacy():
    vertices = tuple(product((0, 1), repeat=3))
    g = lambda v: (1 - v[1], v[0], 1 - v[2])
    c = lambda v: (v[0], v[1], 1 - v[2])
    cycle_a, cycle_b = (1, 15, 13, 9), (3, 11, 5, 7)
    phi = {}
    for roots, start in ((cycle_a, (0, 0, 0)), (cycle_b, (0, 0, 1))):
        p = start
        for root in roots:
            phi[root] = p
            p = g(p)
        require(p == start, "four-cycle closes")
    inv = {v: k for k, v in phi.items()}
    require(len(inv) == 8, "conjugacy is bijective")
    switch = {}
    for root, v in phi.items():
        require(phi[abs(17 - 2 * root)] == g(v), "D17/G conjugacy")
        require(pow(root, 8, 17) == (1 if sum(v) % 2 == 0 else 16),
                "quadratic character matches cube parity")
        require(g(c(v)) == c(g(v)), "commuting sheet swap")
        switch[root] = inv[c(v)]
    require(switch == {1: 3, 15: 11, 13: 5, 9: 7,
                       3: 1, 11: 15, 5: 13, 7: 9}, "transported involution")
    oddrep = lambda x: x % 17 if (x % 17) % 2 else 17 - (x % 17)
    multipliers = tuple(sorted(phi))
    require(not any(all(oddrep(a * v) == switch[v] for v in phi)
                    for a in multipliers), "sheet swap is not a scalar multiplier")
    involutions = [a for a in multipliers if oddrep(a * a) == 1]
    require(involutions == [1, 13], "scalar involutions stay within sectors")
    d = lambda v: abs(17 - 2 * v)
    require(all(oddrep(9 * v) == d(d(d(v))) for v in phi),
            "M3 squared equals D cubed")
    require(sum((x - y) ** 2 for x, y in zip(phi[1], phi[7])) == 1
            and sum((x - y) ** 2 for x, y in zip(phi[3], phi[13])) == 3,
            "M3 transport destroys the cube metric")
    commute, preserve = 0, 0
    for image in permutations(vertices):
        h = dict(zip(vertices, image))
        if all(h[g(v)] == g(h[v]) for v in vertices):
            commute += 1
            preserve += int(all(sum(h[v]) % 2 == sum(v) % 2 for v in vertices))
    require((commute, preserve) == (32, 16), "phase choices / finite centralizer")
    isometry_counts = {"total": 0, "parity_preserving": 0,
                       "parity_and_orientation_preserving": 0}
    isometry_orders = set()
    for perm, flips in product(permutations(range(3)), product((0, 1), repeat=3)):
        inv_count = sum(perm[i] > perm[j] for i in range(3) for j in range(i + 1, 3))
        parity = sum(flips) % 2
        det_sign = (-1) ** (inv_count + sum(flips))
        images = {v: tuple(v[perm[i]] ^ flips[i] for i in range(3)) for v in vertices}
        current = {v: v for v in vertices}
        order = 0
        while True:
            order += 1
            current = {v: images[current[v]] for v in vertices}
            if all(current[v] == v for v in vertices):
                break
        isometry_orders.add(order)
        require(all((sum(images[v]) - sum(v)) % 2 == parity for v in vertices),
                "cube parity character under every isometry")
        isometry_counts["total"] += 1
        isometry_counts["parity_preserving"] += int(parity == 0)
        isometry_counts["parity_and_orientation_preserving"] += int(parity == 0 and det_sign == 1)
    require(list(isometry_counts.values()) == [48, 24, 12], "cube symmetry counts")
    require(isometry_orders == {1, 2, 3, 4, 6}, "no order-eight cube isometry")
    canonical_matching = lambda pairs: tuple(sorted(tuple(sorted(p)) for p in pairs))
    matchings = tuple(canonical_matching(p) for p in
                      (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))))
    quotient = {}
    for p in permutations(range(4)):
        action = tuple(matchings.index(canonical_matching([(p[a], p[b]) for a, b in m]))
                       for m in matchings)
        quotient[action] = quotient.get(action, 0) + 1
    require(len(quotient) == 6 and set(quotient.values()) == {4},
            "tetrahedral S4/V4 quotient S3")
    return {"phi": {str(k): list(v) for k, v in sorted(phi.items())},
            "transported_switch": {str(k): v for k, v in sorted(switch.items())},
            "scalar_involutions": involutions, "centralizer_order": commute,
            "character_preserving_centralizer_order": preserve,
            "cube_isometries": isometry_counts,
            "cube_isometry_orders": sorted(isometry_orders),
            "matching_quotient_fiber_sizes": sorted(quotient.values())}


def main():
    global CHECKS
    CHECKS = 0
    data = {"status": "FINITE-EXACT; all-quantifier proofs are in the note",
            "demicube_vertex_enumeration": check_demicube_vertices(),
            "cube": check_cube(), "rectangles": check_rectangles(),
            "dimension_table": check_general_dimensions(),
            "finite_conjugacy": check_finite_conjugacy()}
    data["checks"] = CHECKS
    source = Path(__file__).read_bytes().replace(b"\r\n", b"\n")
    data["source_sha256_lf"] = hashlib.sha256(source).hexdigest()
    return data


if __name__ == "__main__":
    result = main()
    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"PASS: {result['checks']} exact checks; wrote {output.name}")
