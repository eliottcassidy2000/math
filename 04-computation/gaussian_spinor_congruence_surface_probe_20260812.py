#!/usr/bin/env python3
"""Exact finite-field audit of the Gaussian-spinor congruence surfaces.

For an odd prime p the vertices are nonzero vectors in F_p^2 modulo sign,
and two vertices are adjacent when their determinant has square one.  All
graph triangles are filled.  The probe checks the resulting closed orientable
surface, its projection to P^1(F_p), and the Gaussian-square light-cone map.

The calculation uses only exact integer and finite-field arithmetic.  Every
truth-bearing check uses require(), so python3 and python3 -O execute the same
validation path.
"""

from collections import Counter, defaultdict, deque
from itertools import combinations, combinations_with_replacement
from math import comb, gcd


PRIMES = (3, 5, 7, 11, 13)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def canon_pm(vector, prime):
    vector = tuple(coordinate % prime for coordinate in vector)
    negative = tuple((-coordinate) % prime for coordinate in vector)
    return min(vector, negative)


def vertices_mod_sign(prime):
    return tuple(sorted({
        canon_pm((x, y), prime)
        for x in range(prime)
        for y in range(prime)
        if (x, y) != (0, 0)
    }))


def det2(left, right, prime):
    return (left[0] * right[1] - left[1] * right[0]) % prime


def edge_key(left, right):
    return tuple(sorted((left, right)))


def gaussian_square(vector, prime):
    m, n = vector
    return (
        (m * m - n * n) % prime,
        (2 * m * n) % prime,
        (m * m + n * n) % prime,
    )


def lorentz(left, right, prime):
    return (
        left[2] * right[2]
        - left[0] * right[0]
        - left[1] * right[1]
    ) % prime


def projective_slope(vector, prime):
    x, y = vector
    if x % prime:
        inverse = pow(x, -1, prime)
        return (1, y * inverse % prime)
    return (0, 1)


def quadratic_character(value, prime):
    value %= prime
    require(value != 0, (prime, value, "quadratic character at zero"))
    residue = pow(value, (prime - 1) // 2, prime)
    require(residue in (1, prime - 1), (prime, value, residue))
    return 1 if residue == 1 else -1


def graph_and_faces(prime):
    vertices = vertices_mod_sign(prime)
    edges = frozenset(
        edge_key(left, right)
        for left, right in combinations(vertices, 2)
        if det2(left, right, prime) ** 2 % prime == 1
    )
    adjacency = {vertex: set() for vertex in vertices}
    for left, right in edges:
        adjacency[left].add(right)
        adjacency[right].add(left)
    faces = tuple(
        triple
        for triple in combinations(vertices, 3)
        if all(edge_key(left, right) in edges
               for left, right in combinations(triple, 2))
    )
    return vertices, edges, adjacency, faces


def require_connected(vertices, adjacency, label):
    reached = {vertices[0]}
    queue = deque((vertices[0],))
    while queue:
        vertex = queue.popleft()
        for neighbor in sorted(adjacency[vertex]):
            if neighbor not in reached:
                reached.add(neighbor)
                queue.append(neighbor)
    require(len(reached) == len(vertices), (label, "disconnected", len(reached)))


def oriented_edge_sign(face, edge):
    directed_boundary = (
        (face[0], face[1]),
        (face[1], face[2]),
        (face[2], face[0]),
    )
    for start, end in directed_boundary:
        if edge_key(start, end) == edge:
            return 1 if (start, end) == edge else -1
    raise RuntimeError((face, edge, "edge absent from face"))


def require_orientable(faces, edge_faces, label):
    face_sign = {faces[0]: 1}
    queue = deque((faces[0],))
    while queue:
        face = queue.popleft()
        for edge in combinations(face, 2):
            edge = edge_key(*edge)
            incident = edge_faces[edge]
            require(len(incident) == 2, (label, edge, "not two incident faces"))
            other = incident[0] if incident[1] == face else incident[1]
            forced = (
                -face_sign[face]
                * oriented_edge_sign(face, edge)
                * oriented_edge_sign(other, edge)
            )
            if other in face_sign:
                require(face_sign[other] == forced,
                        (label, edge, "orientation conflict"))
            else:
                face_sign[other] = forced
                queue.append(other)
    require(len(face_sign) == len(faces),
            (label, "face dual disconnected", len(face_sign), len(faces)))


def require_link_cycle(vertex, adjacency, edges, prime):
    neighbors = tuple(sorted(adjacency[vertex]))
    link_edges = frozenset(
        edge_key(left, right)
        for left, right in combinations(neighbors, 2)
        if edge_key(left, right) in edges
    )
    link_adjacency = {neighbor: set() for neighbor in neighbors}
    for left, right in link_edges:
        link_adjacency[left].add(right)
        link_adjacency[right].add(left)
    require(len(neighbors) == prime, (prime, vertex, "link vertex count"))
    require(len(link_edges) == prime, (prime, vertex, "link edge count"))
    require(all(len(link_adjacency[node]) == 2 for node in neighbors),
            (prime, vertex, "link is not 2-regular"))
    require_connected(neighbors, link_adjacency, (prime, vertex, "link"))


def base_triangle_lifts(vertices, faces, prime):
    slopes = tuple(sorted({projective_slope(vertex, prime)
                           for vertex in vertices}))
    require(len(slopes) == prime + 1, (prime, "slope count"))
    lift_count = Counter()
    for face in faces:
        base = tuple(sorted(projective_slope(vertex, prime) for vertex in face))
        require(len(set(base)) == 3, (prime, face, "collapsed face projection"))
        lift_count[base] += 1

    histogram = Counter()
    for base in combinations(slopes, 3):
        actual = lift_count[base]
        if prime % 4 == 3:
            expected = 1
        else:
            determinant_product = 1
            for left, right in ((base[0], base[1]),
                                (base[1], base[2]),
                                (base[2], base[0])):
                determinant_product *= det2(left, right, prime)
            expected = 2 if quadratic_character(
                determinant_product, prime
            ) == 1 else 0
        require(actual == expected,
                (prime, base, actual, expected, "base-triangle lift law"))
        histogram[actual] += 1
    require(sum(histogram.values()) == comb(prime + 1, 3),
            (prime, "base-triangle histogram size"))
    return slopes, histogram


def audit_prime(prime):
    vertices, edges, adjacency, faces = graph_and_faces(prime)
    expected_vertices = (prime * prime - 1) // 2
    expected_edges = prime * (prime * prime - 1) // 4
    expected_faces = prime * (prime * prime - 1) // 6
    require(len(vertices) == expected_vertices, (prime, "vertex formula"))
    require(len(edges) == expected_edges, (prime, "edge formula"))
    require(len(faces) == expected_faces, (prime, "face formula"))
    require(all(len(adjacency[vertex]) == prime for vertex in vertices),
            (prime, "degree formula"))
    require_connected(vertices, adjacency, (prime, "surface graph"))

    edge_faces = defaultdict(list)
    for face in faces:
        for edge in combinations(face, 2):
            edge_faces[edge_key(*edge)].append(face)
    require(set(edge_faces) == set(edges), (prime, "edge-face support"))
    require(all(len(incident) == 2 for incident in edge_faces.values()),
            (prime, "two faces per edge"))
    for vertex in vertices:
        require_link_cycle(vertex, adjacency, edges, prime)
    require_orientable(faces, edge_faces, prime)

    spinor_images = {vertex: gaussian_square(vertex, prime)
                     for vertex in vertices}
    require(len(set(spinor_images.values())) == len(vertices),
            (prime, "Gaussian-square map not injective modulo sign"))
    for vertex, image in spinor_images.items():
        require(lorentz(image, image, prime) == 0,
                (prime, vertex, image, "non-null spinor image"))
    for left, right in combinations_with_replacement(vertices, 2):
        expected_pairing = 2 * det2(left, right, prime) ** 2 % prime
        require(lorentz(spinor_images[left], spinor_images[right], prime)
                == expected_pairing,
                (prime, left, right, "Lorentz determinant identity"))

    slopes, lift_histogram = base_triangle_lifts(vertices, faces, prime)
    fibre_size = (prime - 1) // 2
    fibre_counts = Counter(projective_slope(vertex, prime)
                           for vertex in vertices)
    require(set(fibre_counts) == set(slopes), (prime, "projection support"))
    require(set(fibre_counts.values()) == {fibre_size},
            (prime, "projection fibre size"))
    for vertex in vertices:
        own_slope = projective_slope(vertex, prime)
        neighbor_slopes = Counter(projective_slope(neighbor, prime)
                                  for neighbor in adjacency[vertex])
        require(own_slope not in neighbor_slopes,
                (prime, vertex, "neighbor in own fibre"))
        require(set(neighbor_slopes) == set(slopes) - {own_slope},
                (prime, vertex, "missing neighbor slope"))
        require(set(neighbor_slopes.values()) == {1},
                (prime, vertex, "nonunique neighbor over slope"))

    euler_characteristic = len(vertices) - len(edges) + len(faces)
    genus_numerator = 2 - euler_characteristic
    require(genus_numerator >= 0 and genus_numerator % 2 == 0,
            (prime, euler_characteristic, "orientable genus not integral"))
    genus = genus_numerator // 2
    formula_numerator = 24 + (prime * prime - 1) * (prime - 6)
    require(formula_numerator % 24 == 0,
            (prime, formula_numerator, "genus formula not integral"))
    require(genus == formula_numerator // 24,
            (prime, genus, formula_numerator, "genus formula"))

    return {
        "prime": prime,
        "vertices": vertices,
        "edges": edges,
        "adjacency": adjacency,
        "faces": faces,
        "V": len(vertices),
        "E": len(edges),
        "F": len(faces),
        "chi": euler_characteristic,
        "genus": genus,
        "fibre": fibre_size,
        "lifts": lift_histogram,
    }


def frame_coordinate(vertex, prime):
    x, y = vertex
    if y % prime == 0:
        return ("inf", x * x % prime)
    inverse_y = pow(y, -1, prime)
    return (
        "aff",
        x * inverse_y % prime,
        inverse_y * inverse_y % prime,
    )


def frame_edge_law(left, right, prime):
    if left[0] == "inf" and right[0] == "inf":
        return False
    if left[0] == "inf":
        return left[1] == right[2]
    if right[0] == "inf":
        return right[1] == left[2]
    return (left[1] - right[1]) ** 2 % prime == left[2] * right[2] % prime


def audit_p13_frame(row):
    prime = 13
    vertices = row["vertices"]
    edges = row["edges"]
    adjacency = row["adjacency"]
    faces = set(row["faces"])
    squares = tuple(value for value in range(1, prime)
                    if quadratic_character(value, prime) == 1)
    coordinates = {vertex: frame_coordinate(vertex, prime)
                   for vertex in vertices}
    require(len(set(coordinates.values())) == len(vertices),
            "p=13 frame coordinates not bijective")
    expected_coordinates = {
        *(('inf', eta) for eta in squares),
        *(('aff', x, eta) for x in range(prime) for eta in squares),
    }
    require(set(coordinates.values()) == expected_coordinates,
            "p=13 frame coordinate support")
    inverse_coordinates = {coordinate: vertex
                           for vertex, coordinate in coordinates.items()}
    for left, right in combinations(vertices, 2):
        actual = edge_key(left, right) in edges
        expected = frame_edge_law(coordinates[left], coordinates[right], prime)
        require(actual == expected,
                (coordinates[left], coordinates[right], "p=13 frame edge law"))

    wheel_vertices = []
    wheel_faces = set()
    paley_arcs = set()
    for eta in squares:
        center = inverse_coordinates[("inf", eta)]
        rim = {inverse_coordinates[("aff", x, eta)] for x in range(prime)}
        wheel_vertices.append({center} | rim)
        require(adjacency[center] == rim,
                (eta, "p=13 wheel spoke set"))
        expected_rim_edges = {
            edge_key(inverse_coordinates[("aff", x, eta)],
                     inverse_coordinates[("aff", (x + eta) % prime, eta)])
            for x in range(prime)
        }
        actual_rim_edges = {
            edge_key(left, right)
            for left, right in combinations(rim, 2)
            if edge_key(left, right) in edges
        }
        require(actual_rim_edges == expected_rim_edges,
                (eta, "p=13 wheel rim cycle"))
        for x in range(prime):
            face = tuple(sorted((
                center,
                inverse_coordinates[("aff", x, eta)],
                inverse_coordinates[("aff", (x + eta) % prime, eta)],
            )))
            require(face in faces, (eta, x, "p=13 wheel face"))
            wheel_faces.add(face)
            paley_arcs.add((x, (x + eta) % prime))

    require(sum(len(wheel) for wheel in wheel_vertices) == len(vertices),
            "p=13 wheel vertex total")
    require(set().union(*wheel_vertices) == set(vertices),
            "p=13 wheel vertex support")
    require(sum(len(left & right) for left, right
                in combinations(wheel_vertices, 2)) == 0,
            "p=13 wheels not vertex-disjoint")
    require(len(wheel_faces) == 6 * 13,
            "p=13 wheel face count")
    infinity_vertices = {inverse_coordinates[("inf", eta)] for eta in squares}
    incident_to_infinity = {
        face for face in faces if set(face) & infinity_vertices
    }
    require(incident_to_infinity == wheel_faces,
            "p=13 wheel faces are not exactly infinity-incident faces")
    expected_paley_arcs = {
        (start, end)
        for start in range(prime)
        for end in range(prime)
        if start != end
        and quadratic_character(end - start, prime) == 1
    }
    require(paley_arcs == expected_paley_arcs,
            "p=13 wheel/Paley arc correspondence")
    return len(expected_coordinates), len(wheel_faces), len(paley_arcs)


def audit_p2_hostile():
    prime = 2
    vertices, edges, adjacency, faces = graph_and_faces(prime)
    images = {vertex: gaussian_square(vertex, prime) for vertex in vertices}
    require((len(vertices), len(edges), len(faces)) == (3, 3, 1),
            "p=2 graph hostile counts")
    require(set(len(adjacency[vertex]) for vertex in vertices) == {2},
            "p=2 graph hostile degrees")
    require(images[canon_pm((1, 0), prime)]
            == images[canon_pm((0, 1), prime)] == (1, 0, 1),
            "p=2 Gaussian-square collision")
    require(images[canon_pm((1, 1), prime)] == (0, 0, 0),
            "p=2 zero spinor image")
    return len(vertices), len(edges), len(faces)


def integer_phi(vector):
    m, n = vector
    return (m * m - n * n, 2 * m * n, m * m + n * n)


def lrc_gate_invoice(vector):
    a, b = vector
    q = a * a + b * b
    denominator = max(13 * abs(b), abs(a - 12 * b))
    return q, denominator, 91 * denominator <= q


def audit_lrc_hostiles():
    low = (1, 0)
    high = (1, 91 * 13)
    low_q, low_d, low_pass = lrc_gate_invoice(low)
    high_q, high_d, high_pass = lrc_gate_invoice(high)
    require(tuple(value % 13 for value in low)
            == tuple(value % 13 for value in high),
            "LRC same-residue hostile setup")
    require((low_q, low_d, low_pass) == (1, 1, False),
            "LRC same-residue hostile low invoice")
    require((high_q, high_d, high_pass)
            == (1_399_490, 15_379, True),
            "LRC same-residue hostile high invoice")
    require(high_q - 91 * high_d == 1,
            "LRC same-residue hostile boundary slack")

    kelvin_spinor = (5, 1)
    kelvin_image = integer_phi(kelvin_spinor)
    require(gcd(*kelvin_spinor) == 1,
            "Kelvin singular-residue spinor not primitive")
    require(kelvin_image == (24, 10, 26),
            "Kelvin singular-residue integer image")
    require(tuple(value % 13 for value in kelvin_image) == (11, 10, 0),
            "Kelvin singular-residue reduction")
    require(kelvin_image[2] % 13 == 0,
            "Kelvin singular-residue denominator unexpectedly invertible")

    columns = tuple((index, 0) if index != 12 else (12, 1)
                    for index in range(1, 14))
    residues = tuple((a % 13, b % 13) for a, b in columns)
    zero_count = residues.count((0, 0))
    nonzero_slopes = Counter(projective_slope(residue, 13)
                             for residue in residues if residue != (0, 0))
    require(zero_count == 1, "one-tail residue hostile zero count")
    require(nonzero_slopes == Counter({(1, 0): 11, (1, 12): 1}),
            (nonzero_slopes, "one-tail residue hostile fibre collapse"))
    return {
        "low": (low_q, low_d, low_pass),
        "high": (high_q, high_d, high_pass),
        "kelvin": kelvin_image,
        "zero": zero_count,
    }


def histogram_text(histogram):
    return ",".join(f"{key}:{histogram[key]}" for key in sorted(histogram))


def main():
    rows = [audit_prime(prime) for prime in PRIMES]
    frame_vertices, wheel_faces, paley_arcs = audit_p13_frame(rows[-1])
    p2_v, p2_e, p2_f = audit_p2_hostile()
    hostiles = audit_lrc_hostiles()

    print("GAUSSIAN SPINOR CONGRUENCE SURFACE PROBE 2026-08-12")
    print("universe=Vp=(F_p^2-{0})/{+-1}; det^2=1 edges; all graph triangles")
    print("spinor=Phi(m,n)=(m^2-n^2,2mn,m^2+n^2); exact arithmetic")
    for row in rows:
        print(
            "p={prime} pmod4={pmod4} pmod12={pmod12} V={V} E={E} F={F} degree={prime} "
            "edge_faces=2 links=C{prime} orientable=yes chi={chi} genus={genus} "
            "cover=K{base} fibre={fibre} one_neighbor_each_other_slope=yes "
            "base_triangle_lifts={lifts} spinor_injective=yes".format(
                prime=row["prime"],
                pmod4=row["prime"] % 4,
                pmod12=row["prime"] % 12,
                V=row["V"],
                E=row["E"],
                F=row["F"],
                chi=row["chi"],
                genus=row["genus"],
                base=row["prime"] + 1,
                fibre=row["fibre"],
                lifts=histogram_text(row["lifts"]),
            )
        )
    print(
        "p13_frame: coordinates={} law=[no inf-inf; (inf,e)~(x,t) iff e=t; "
        "(x,e)~(y,t) iff (x-y)^2=e*t]".format(frame_vertices)
    )
    print(
        "p13_wheels: 6 vertex-disjoint W14 disks; each V=14 E=26 F=13; "
        f"infinity_faces={wheel_faces}; directed_Paley_arcs={paley_arcs}"
    )
    print(
        f"p2_hostile: V={p2_v} E={p2_e} F={p2_f}; edge_faces=1; "
        "Phi(1,0)=Phi(0,1)=(1,0,1); Phi(1,1)=(0,0,0)"
    )
    print(
        "lrc_same_residue_mod13: (1,0) q={} D={} gate={}; "
        "(1,1183) q={} D={} gate={} slack=1".format(
            hostiles["low"][0], hostiles["low"][1], hostiles["low"][2],
            hostiles["high"][0], hostiles["high"][1], hostiles["high"][2],
        )
    )
    print(
        "lrc_kelvin_hostile: primitive_spinor=(5,1) Phi={} mod13=(11,10,0); "
        "q_inverse=no".format(hostiles["kelvin"])
    )
    print(
        "lrc_one_tail_residue_hostile: 11x[1:0]+1x[12:1]+1xzero; "
        f"zero_count={hostiles['zero']}"
    )
    print("ALL REQUIRE CHECKS PASSED")


if __name__ == "__main__":
    main()
