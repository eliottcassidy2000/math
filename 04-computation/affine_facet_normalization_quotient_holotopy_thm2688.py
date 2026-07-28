#!/usr/bin/env python3
"""Exact hostile audit for the abstract THM-2672 facet-torsor calculus.

This is the exact THM-2688 companion.  It checks the simplicial image, the
diagonal cyclic quotient, and mapping-cone chain complexes in small prime
controls, while printing the closed n=13 formulas.
"""

from itertools import combinations
from math import comb, gcd


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def rank_mod(matrix, prime):
    """Rank of an integer matrix over F_prime."""
    if not matrix:
        return 0
    a = [[entry % prime for entry in row] for row in matrix]
    rows, cols = len(a), len(a[0]) if a else 0
    rank = 0
    for col in range(cols):
        pivot = next((r for r in range(rank, rows) if a[r][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        inv = pow(a[rank][col], -1, prime)
        a[rank] = [(inv * x) % prime for x in a[rank]]
        for r in range(rows):
            if r != rank and a[r][col]:
                scalar = a[r][col]
                a[r] = [(x - scalar * y) % prime
                        for x, y in zip(a[r], a[rank])]
        rank += 1
        if rank == rows:
            break
    return rank


def boundary_terms(face):
    for i in range(len(face)):
        yield face[:i] + face[i + 1 :], -1 if i % 2 else 1


def boundary_faces(n):
    """Faces of boundary Delta^(n-1), grouped by dimension."""
    return {
        k: tuple(combinations(range(n), k + 1))
        for k in range(n - 1)
    }


def carry_faces(n, missing):
    """Faces of the disjoint facet resolution, delta-oriented."""
    return {
        k: tuple(
            (c, ds)
            for c in range(n)
            for ds in combinations(
                tuple(d for d in range(n) if d != missing(c)), k + 1
            )
        )
        for k in range(n - 1)
    }


def quotient_faces(n):
    """Faces of Delta^(n-2) on nonzero invariant t."""
    return {
        k: tuple(combinations(range(1, n), k + 1))
        for k in range(n - 1)
    }


def mapping_cone_betti(n, target_faces, source_faces, image, prime):
    """Betti numbers of Cone(source -> target) over F_prime.

    The cone differential on C_k(Y) + C_(k-1)(X) is
      (y,x) -> (d_Y y + f x, -d_X x).
    """
    max_target = max(target_faces, default=-1)
    max_source = max(source_faces, default=-1)
    max_cone = max(max_target, max_source + 1)
    cone_bases = {}
    for k in range(max_cone + 1):
        cone_bases[k] = (
            tuple(("Y", face) for face in target_faces.get(k, ()))
            + tuple(("X", face) for face in source_faces.get(k - 1, ()))
        )
    # The topological cone CX has one common apex, including when X is
    # disconnected.  The unaugmented algebraic-cone formula omits this vertex;
    # include it explicitly so the degree-one cone edges have boundary
    # f(v)-apex.
    cone_bases[0] += (("A", ()),)

    ranks = {0: 0, max_cone + 1: 0}
    for k in range(1, max_cone + 1):
        domain = cone_bases[k]
        codomain = cone_bases[k - 1]
        row_index = {basis: i for i, basis in enumerate(codomain)}
        matrix = [[0] * len(domain) for _ in codomain]
        for j, (kind, face) in enumerate(domain):
            if kind == "Y":
                for subface, sign in boundary_terms(face):
                    matrix[row_index[("Y", subface)]][j] += sign
            else:
                # The source face lives one degree lower than its cone cell.
                target_face, orient = image(face)
                matrix[row_index[("Y", target_face)]][j] += orient
                if k == 1:
                    matrix[row_index[("A", ())]][j] -= 1
                if len(face[1] if isinstance(face[0], int) and
                           isinstance(face[1], tuple) else face) > 1:
                    if isinstance(face[0], int) and isinstance(face[1], tuple):
                        c, vertices = face
                        for subface, sign in boundary_terms(vertices):
                            matrix[row_index[("X", (c, subface))]][j] -= sign
                    else:
                        for subface, sign in boundary_terms(face):
                            matrix[row_index[("X", subface)]][j] -= sign
        ranks[k] = rank_mod(matrix, prime)

    return tuple(
        len(cone_bases[k]) - ranks.get(k, 0) - ranks.get(k + 1, 0)
        for k in range(max_cone + 1)
    )


def vertical_image_faces(n, missing):
    image = set()
    for c in range(n):
        active = tuple(d for d in range(n) if d != missing(c))
        for size in range(n):
            image.update(combinations(active, size))
    return image


def audit_prime(n):
    require(n in (3, 5, 7, 13), "unexpected prime control")
    # The n=13 theorem has m(c)=12-2c.  The uniform affine model below is
    # m(c)=-1-2c mod n, which specializes to it.
    missing = lambda c: (-1 - 2 * c) % n
    require(sorted(missing(c) for c in range(n)) == list(range(n)),
            "affine missing map lost bijectivity")

    K = carry_faces(n, missing)
    B = boundary_faces(n)
    Q = quotient_faces(n)
    image = vertical_image_faces(n, missing)
    expected_boundary = {
        face
        for size in range(n)
        for face in combinations(range(n), size)
    }
    require(image == expected_boundary,
            "vertical simplicial image is not the full boundary")

    # Diagonal generator: for n=13, s=7 and r=-1.  For the other prime
    # controls choose s=1 and its compatible r=-2.
    s = 7 if n == 13 else 1
    r = (-2 * s) % n
    vertices = tuple(
        (c, d) for c in range(n) for d in range(n) if d != missing(c)
    )
    g = lambda v: ((v[0] + s) % n, (v[1] + r) % n)
    require(all(g(v)[1] != missing(g(v)[0]) for v in vertices),
            "diagonal generator left the refined complex")
    require(all((g(v)[1] - missing(g(v)[0])) % n ==
                (v[1] - missing(v[0])) % n for v in vertices),
            "diagonal invariant changed")
    for v in vertices:
        orbit = []
        w = v
        for _ in range(n):
            orbit.append(w)
            w = g(w)
        require(w == v and len(set(orbit)) == n,
                "diagonal action is not free of order n")
    vertex_orbits = {
        (d - missing(c)) % n for c, d in vertices
    }
    require(vertex_orbits == set(range(1, n)),
            "diagonal quotient lost a nonzero offset")

    # The compatible label rotation is free on the geometric boundary when n
    # is prime: every nonidentity power is one n-cycle, whose only fixed point
    # in the full simplex is the (interior) barycenter.
    for power in range(1, n):
        step = (power * r) % n
        require(gcd(step, n) == 1,
                "nonidentity label rotation is not transitive")
        orbit = []
        d = 0
        for _ in range(n):
            orbit.append(d)
            d = (d + step) % n
        require(d == 0 and len(set(orbit)) == n,
                "label rotation gained a boundary fixed orbit")

    fK = tuple(len(K[k]) for k in range(n - 1))
    fB = tuple(len(B[k]) for k in range(n - 1))
    fQ = tuple(len(Q[k]) for k in range(n - 1))
    require(fK == tuple(n * comb(n - 1, k + 1) for k in range(n - 1)),
            "refined f-vector changed")
    require(fB == tuple(comb(n, k + 1) for k in range(n - 1)),
            "boundary f-vector changed")
    require(fQ == tuple(comb(n - 1, k + 1) for k in range(n - 1)),
            "diagonal quotient f-vector changed")

    if n <= 7:
        # Chain-level positive controls for both maps.  Several characteristics
        # rule out accidental rational-only rank behavior in these controls.
        image_pi = lambda face: (face[1], 1)

        # Reorient K by invariant t for the orbit quotient chain map.
        Kt = {
            k: tuple(
                (c, ts)
                for c in range(n)
                for ts in combinations(range(1, n), k + 1)
            )
            for k in range(n - 1)
        }
        image_q = lambda face: (face[1], 1)
        for prime in (2, 3, 5, 7, 11):
            beta_pi = mapping_cone_betti(n, B, K, image_pi, prime)
            expected_pi = [0] * n
            expected_pi[0] = 1
            expected_pi[1] += n - 1
            expected_pi[n - 2] += 1
            require(beta_pi == tuple(expected_pi),
                    f"vertical mapping-cone Betti mismatch {(n, prime, beta_pi)}")

            beta_q = mapping_cone_betti(n, Q, Kt, image_q, prime)
            expected_q = [0] * n
            expected_q[0] = 1
            expected_q[1] = n - 1
            require(beta_q == tuple(expected_q),
                    f"orbit mapping-cone Betti mismatch {(n, prime, beta_q)}")

    return fK, fB, fQ


def hostile_controls():
    # Non-surjective missing-label map: the vertical image is a filled
    # simplex on labels 1,...,4, not boundary Delta^4.
    n = 5
    constant = lambda _c: 0
    image = vertical_image_faces(n, constant)
    require((1, 2, 3, 4) in image, "constant-map hostile lost its top face")
    require((0,) not in image and (0, 1) not in image,
            "constant-map hostile unexpectedly gained label zero")
    require(len(image) == 2 ** (n - 1),
            "constant-map hostile is not one filled simplex")

    # A non-generator component step on Z/6 has two component orbits, so its
    # quotient is two disjoint Delta^4 rather than one Delta^4.
    n = 6
    s = 2
    require(gcd(s, n) == 2, "nongenerator hostile changed")
    component_orbits = []
    unseen = set(range(n))
    while unseen:
        c = min(unseen)
        orbit = set()
        x = c
        while x not in orbit:
            orbit.add(x)
            x = (x + s) % n
        unseen -= orbit
        component_orbits.append(orbit)
    require(len(component_orbits) == 2,
            "nongenerator hostile did not split the quotient")


def lens_and_bockstein_controls():
    """Exact p=13 character, cellular-homology, and extension checks."""
    p = 13
    weights = tuple(range(1, (p - 1) // 2 + 1))
    paired_characters = {
        residue
        for weight in weights
        for residue in (weight, (-weight) % p)
    }
    require(weights == (1, 2, 3, 4, 5, 6)
            and paired_characters == set(range(1, p)),
            "reduced regular representation lost a Fourier weight")

    # One cellular generator in each degree 0,...,11, with differential p in
    # positive even degrees and zero in odd degrees, gives the standard lens
    # homology pattern.
    homology = []
    for degree in range(12):
        incoming = p if degree + 1 < 12 and (degree + 1) % 2 == 0 else 0
        outgoing = p if degree > 0 and degree % 2 == 0 else 0
        if outgoing:
            homology.append("0")
        elif incoming:
            homology.append("Z/13")
        else:
            homology.append("Z")
    require(tuple(homology) == (
        "Z", "Z/13", "0", "Z/13", "0", "Z/13",
        "0", "Z/13", "0", "Z/13", "0", "Z",
    ), "lens cellular homology changed")

    # THM-2657 gauges: i(a)=13a, pi(k)=2k mod13, alpha(g)=1, s(1)=7.
    total_modulus = p**6
    kernel_modulus = p**5
    lift = 7
    require((2 * lift) % p == 1, "physical lift no longer covers alpha(g)=1")
    defect = p * lift
    require(defect == 91 and defect % p == 0,
            "thirteen-turn extension defect changed")
    kernel_coordinate = (defect // p) % kernel_modulus
    require(kernel_coordinate == 7 and kernel_coordinate % p != 0,
            "coefficient Bockstein is no longer the unit class 7")
    require(((-kernel_coordinate) % p) == 6,
            "inverse-generator gauge hostile changed")
    require(total_modulus == p * kernel_modulus,
            "cyclic extension sizes changed")
    return weights, tuple(homology), kernel_coordinate


def main():
    weights, lens_homology, bockstein = lens_and_bockstein_controls()
    for n in (3, 5, 7, 13):
        fK, fB, fQ = audit_prime(n)
        if n == 13:
            print("n=13")
            print(f"refined_K_fvector={fK}")
            print(f"vertical_image_boundary_fvector={fB}")
            print(f"diagonal_quotient_fvector={fQ}")
            print("K_homology=reduced_H0:12")
            print("vertical_image_homology=reduced_H11:1")
            print("diagonal_quotient_homology=0")
            print("Cone(vertical_projection)=wedge(S11,12*S1)")
            print("Cone(vertical_projection)_reduced_H1=Z^12 H11=Z")
            print("Cone(diagonal_quotient)=12*S1")
            print("mapping_torus(diagonal_generator)=Delta11_x_S1")
            print("label_action_on_boundary=free_C13")
            print("boundary_label_quotient=L11(13;1,2,3,4,5,6)")
            print("lens_reduced_homology=H1,H3,H5,H7,H9:Z/13 H11:Z")
            print(
                f"lens_exact_controls=weights={weights}:"
                f"homology={lens_homology}:bockstein={bockstein}"
            )
    hostile_controls()
    print("positive_controls=n3,n5,n7 chain complexes over F2,F3,F5,F7,F11")
    print("hostiles=constant_missing_map; non-generator_component_step")
    print("physical_boundary=g^13_is_not_identity_after_lift: 91/13^6=7/13^5")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()

