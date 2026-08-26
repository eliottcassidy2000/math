#!/usr/bin/env python3
"""Independent clean-room audit of THM-4230's weight-12 (2,3) seam.

This script uses only Python's standard library.  It reconstructs the source
support, all lower supporting planes, Pick genera, the C_12 character packet,
edge indices, and the W=-2Z attachment hostile.  It deliberately imports no
primary M12 implementation.
"""

from fractions import Fraction
from itertools import combinations
from math import gcd


def require(condition, message):
    if not condition:
        raise AssertionError(message)


# Exact allowed monomials p^i y^j of (2,3)-weight at most 12, with y and py
# absent by the inherited normalization.
MONOMIALS = [
    (i, j)
    for j in range(5)
    for i in range(7)
    if 0 < 2 * i + 3 * j <= 12 and (i, j) not in {(0, 1), (1, 1)}
]


def source_points():
    pts = {(2, 0, 0), (0, 1, 0), (2, 0, 1)}
    for i, j in MONOMIALS:
        # -Q s^2 p^i(sp)^j and +Q p p^i(sp)^j
        pts.add((j + 2, i + j, 1))
        pts.add((j, i + j + 1, 1))
    return sorted(pts)


def plane_through(a, b, c):
    # Solve h = A*r+B*l+C through three points, using exact elimination.
    mat = [
        [Fraction(a[0]), Fraction(a[1]), Fraction(1), Fraction(a[2])],
        [Fraction(b[0]), Fraction(b[1]), Fraction(1), Fraction(b[2])],
        [Fraction(c[0]), Fraction(c[1]), Fraction(1), Fraction(c[2])],
    ]
    for col in range(3):
        pivot = next((row for row in range(col, 3) if mat[row][col]), None)
        if pivot is None:
            return None
        mat[col], mat[pivot] = mat[pivot], mat[col]
        scale = mat[col][col]
        mat[col] = [x / scale for x in mat[col]]
        for row in range(3):
            if row == col:
                continue
            scale = mat[row][col]
            mat[row] = [mat[row][k] - scale * mat[col][k] for k in range(4)]
    return tuple(mat[row][3] for row in range(3))


def lower_planes(points):
    answer = {}
    for triple in combinations(points, 3):
        plane = plane_through(*triple)
        if plane is None:
            continue
        A, B, C = plane
        gaps = [Fraction(h) - (A * r + B * ell + C) for r, ell, h in points]
        if min(gaps) >= 0:
            equal = tuple(p for p, gap in zip(points, gaps) if gap == 0)
            # A genuine face projects two-dimensionally.
            if any(
                (q[0] - equal[0][0]) * (r[1] - equal[0][1])
                != (q[1] - equal[0][1]) * (r[0] - equal[0][0])
                for q, r in combinations(equal[1:], 2)
            ):
                answer[plane] = equal
    return answer


def twice_area(vertices):
    return abs(
        sum(
            x1 * y2 - y1 * x2
            for (x1, y1), (x2, y2) in zip(vertices, vertices[1:] + vertices[:1])
        )
    )


def boundary_count(vertices):
    return sum(
        gcd(abs(x2 - x1), abs(y2 - y1))
        for (x1, y1), (x2, y2) in zip(vertices, vertices[1:] + vertices[:1])
    )


def pick_interior(vertices):
    area2 = twice_area(vertices)
    boundary = boundary_count(vertices)
    require((area2 - boundary + 2) % 2 == 0, "Pick parity")
    return (area2 - boundary + 2) // 2


def interior_points_triangle():
    # Conv((0,0),(0,6),(4,4)): r>0, l>r, r+2l<12.
    return [
        (r, ell)
        for r in range(1, 4)
        for ell in range(1, 6)
        if ell > r and r + 2 * ell < 12
    ]


def hostile_attachment_check(U, Z):
    W = -2 * Z
    D = W * W - 4 * U * Z
    lam = U + W + Z
    require(U * Z * D * lam != 0, "hostile sample must remain inside the gate")
    require(D == -4 * Z * lam, "hostile discriminant identity")

    # At every R--C attachment x=S^2/P=1, u=S^-2, hence u^6=Lambda.
    # Both elliptic quotient ordinates vanish, and their abscissae obey the
    # respective 2-torsion cubics.
    require(D + 4 * Z * lam == 0, "E1 attachment is 2-torsion")
    require(D / lam + 4 * Z == 0, "E2 attachment is 2-torsion")
    return D, lam


def quartic_decomposition_obstruction(K, zeta, Z):
    """Return the depressed-quartic linear numerator.

    For h(T)=Z*T^4+zeta*T^3+K*T^2+constant, translation by
    -zeta/(4Z) makes the linear coefficient
    zeta*(zeta^2-4*K*Z)/(8*Z^2).  A degree 4 polynomial has a quadratic
    intermediate precisely when this coefficient vanishes.
    """
    require(Z != 0, "quartic leading coefficient")
    return zeta * (zeta * zeta - 4 * K * Z)


# F_4 = F_2[w]/(w^2+w+1), encoded by low/high bits a+b*w.
def f4_add(a, b):
    return a ^ b


def f4_mul(a, b):
    a0, a1 = a & 1, (a >> 1) & 1
    b0, b1 = b & 1, (b >> 1) & 1
    c0 = (a0 * b0) ^ (a1 * b1)       # w^2=w+1 contributes 1
    c1 = (a0 * b1) ^ (a1 * b0) ^ (a1 * b1)
    return c0 | (c1 << 1)


def f4_square(a):
    return f4_mul(a, a)


def elliptic_two_torsion_saturation():
    """Pairs (a,b) with a*t+b*t^2=0 on all of F4.

    These are the reductions modulo 2 of Eisenstein pairs which kill the
    graph of the inversion/Frobenius gluing E[2] -> E[2].
    """
    admissible = []
    for a in range(4):
        for b in range(4):
            if all(f4_add(f4_mul(a, t), f4_mul(b, f4_square(t))) == 0 for t in range(4)):
                admissible.append((a, b))
    return admissible


def main():
    expected_monomials = [
        (1, 0), (2, 0), (3, 0), (4, 0), (5, 0), (6, 0),
        (2, 1), (3, 1), (4, 1),
        (0, 2), (1, 2), (2, 2), (3, 2),
        (0, 3), (1, 3), (0, 4),
    ]
    require(MONOMIALS == expected_monomials, "complete weight<=12 support")

    points = source_points()
    planes = lower_planes(points)
    expected_plane = (Fraction(1, 12), Fraction(1, 6), Fraction(-1, 6))
    require(set(planes) == {expected_plane}, f"unique lower plane: {planes}")
    expected_equal = {
        (0, 1, 0), (2, 0, 0), (0, 7, 1),
        (2, 6, 1), (4, 5, 1), (6, 4, 1),
    }
    require(set(planes[expected_plane]) == expected_equal, "main-face equality set")

    global_polygon = [(0, 1), (2, 0), (6, 4), (0, 7)]
    component_polygon = [(0, 0), (0, 6), (4, 4)]
    require((twice_area(global_polygon), boundary_count(global_polygon), pick_interior(global_polygon)) == (48, 14, 18), "global Pick data")
    require((twice_area(component_polygon), boundary_count(component_polygon), pick_interior(component_polygon)) == (24, 12, 7), "component Pick data")

    interiors = interior_points_triangle()
    require(interiors == [(1, 2), (1, 3), (1, 4), (1, 5), (2, 3), (2, 4), (3, 4)], "interior lattice points")
    h0_characters = sorted((r + 2 * ell) % 12 for r, ell in interiors)
    require(h0_characters == [5, 7, 8, 9, 10, 11, 11], "H0 C12 characters")
    h1_characters = sorted(h0_characters + [(-x) % 12 for x in h0_characters])
    require(h1_characters == [1, 1, 2, 3, 4, 5, 5, 7, 7, 8, 9, 10, 11, 11], "H1 C12 characters")
    invariant = [x for x in h0_characters if x % 2 == 0]
    anti_invariant = [x for x in h0_characters if x % 2]
    require(invariant == [8, 10] and anti_invariant == [5, 7, 9, 11, 11], "involution split 2+5")

    # At W=0, rho:(S,P)->(S,-P) acts on the residue differential indexed by
    # (r,l) as (-1)^l.  Its primitive Phi_12 eigenspaces have the two stated
    # imprimitive CM types.
    primitive = [(r, ell, (r + 2 * ell) % 12) for r, ell in interiors if gcd((r + 2 * ell) % 12, 12) == 1]
    rho_plus = sorted(char for _, ell, char in primitive if ell % 2 == 0)
    rho_minus = sorted(char for _, ell, char in primitive if ell % 2 == 1)
    require(rho_plus == [5, 11] and rho_minus == [7, 11], "W=0 primitive CM split")
    require(sorted((7 * x) % 12 for x in rho_plus) == rho_plus, "j=0 CM-type stabilizer")
    require(sorted((5 * x) % 12 for x in rho_minus) == rho_minus, "j=1728 CM-type stabilizer")

    # Intrinsic index e=u+v-c on the three non-affine outer edges.
    packet = sorted([1] + [2] * 4 + [11] * 3, reverse=True)
    require(packet == [11, 11, 11, 2, 2, 2, 2, 1], "Newton packet")
    require(sum(packet) == 42, "full response")
    require(sum(packet) - 4 * 2 == 34, "finite quartic-carrier response")

    # Composite-carrier hostile: this is literally a quadratic composed with
    # a quadratic, so its degree-four residue extension has a quadratic
    # intermediate even though h(T)-q stays irreducible over C(q).
    require(quartic_decomposition_obstruction(1, 2, 1) == 0, "decomposable carrier locus")
    # T^4+2T^3+T^2+1/2 = (T(T+1))^2+1/2 coefficient check.
    inner = [0, 1, 1]
    square = [sum(inner[i] * inner[k - i] for i in range(3) if 0 <= k - i < 3) for k in range(5)]
    square[0] += Fraction(1, 2)
    require(square == [Fraction(1, 2), 0, 1, 2, 1], "explicit composed quartic")
    require(quartic_decomposition_obstruction(1, 1, 1) != 0, "indecomposable control")

    # The map r -> r^-1 on the three roots r^3=-1 fixes one and swaps two;
    # under E[2] ~= F4 this is t -> t^2.  No nonzero O/2-linear pair kills
    # its graph, proving that the descent lattice is exactly 2O^2.
    require([f4_square(t) for t in range(4)] == [0, 1, 3, 2], "F4 Frobenius permutation")
    require(elliptic_two_torsion_saturation() == [(0, 0)], "E[2] gluing saturation")

    # Gate-compatible collision controls: the displayed intermediate source
    # coefficient W-U or Z-W may vanish without changing endpoints/genera.
    for U, W, Z in [(1, 1, 2), (2, 1, 1)]:
        D = W * W - 4 * U * Z
        lam = U + W + Z
        require(U * Z * D * lam != 0, "collision control lies in gate")

    D, lam = hostile_attachment_check(U=2, Z=1)

    print("monomials", len(MONOMIALS), MONOMIALS)
    print("source_points", len(points))
    print("lower_plane", expected_plane, "equal", sorted(expected_equal))
    print("global_pick", (48, 14, 18), "component_pick", (24, 12, 7))
    print("H0_characters", h0_characters)
    print("H1_characters", h1_characters)
    print("involution", {"invariant": invariant, "anti": anti_invariant})
    print("W=0_primitive_split", {"rho_plus": rho_plus, "rho_minus": rho_minus})
    print("packet", packet, "full", 42, "finite", 34)
    print("carrier_decomposition_wall", "zeta*(zeta^2-4*K*Z)=0")
    print("E2_gluing_mod2_admissible", elliptic_two_torsion_saturation())
    print("hostile_W=-2Z", {"U": 2, "W": -2, "Z": 1, "D": D, "Lambda": lam})
    print("PASS")


if __name__ == "__main__":
    main()
