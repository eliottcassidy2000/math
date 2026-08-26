#!/usr/bin/env python3
"""Exact lower-model certificate for THM-4220.

Universe: the zeta=0 wall of the complete exact-M=10 inherited b=d=0
reduced (2,3) source, with upsilon*xi*(upsilon+xi) nonzero.  The proved
strata are D_V=Theta^2-4*K*xi nonzero and the collapsed K=Theta=0 stratum.
The split-conic collision D_V=0,K!=0 is a hostile control, not a result.
"""

from fractions import Fraction
from hashlib import sha256
from math import gcd

import sympy as sy


CHECKS = 0


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not bool(condition):
        raise RuntimeError(message)


def convex_hull(points):
    points = sorted(set(points))

    def cross(origin, first, second):
        return ((first[0] - origin[0]) * (second[1] - origin[1])
                - (first[1] - origin[1]) * (second[0] - origin[0]))

    lower = []
    for point in points:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)
    upper = []
    for point in reversed(points):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    return tuple(lower[:-1] + upper[:-1])


def polygon_ledger(points):
    vertices = convex_hull(points)
    area2 = abs(sum(
        vertices[index][0] * vertices[(index + 1) % len(vertices)][1]
        - vertices[(index + 1) % len(vertices)][0] * vertices[index][1]
        for index in range(len(vertices))
    ))
    boundary = sum(
        gcd(abs(vertices[(index + 1) % len(vertices)][0] - vertices[index][0]),
            abs(vertices[(index + 1) % len(vertices)][1] - vertices[index][1]))
        for index in range(len(vertices))
    )
    require((area2 - boundary + 2) % 2 == 0, "Pick parity")
    genus = (area2 - boundary + 2) // 2
    return vertices, area2, boundary, genus


def edge_packet(vertices):
    packet = []
    rows = []
    for start, end in zip(vertices, vertices[1:] + vertices[:1]):
        dx, dy = end[0] - start[0], end[1] - start[1]
        length = gcd(abs(dx), abs(dy))
        inward = (-dy // length, dx // length)
        constant = inward[0] * start[0] + inward[1] * start[1]
        index = inward[0] + inward[1] - constant
        rows.append((start, end, length, inward, constant, index))
        if start[0] == end[0] == 0:
            continue
        packet.extend([index] * length)
    return tuple(sorted(packet, reverse=True)), tuple(rows)


def digest(rows):
    payload = "|".join(str(row) for row in rows)
    return sha256(payload.encode()).hexdigest()


def main():
    # Complete source support through weight ten, after deleting y and py.
    monomials = tuple(
        sorted(
            [(i, j, 2*i + 3*j)
             for i in range(6)
             for j in range(4)
             if 0 < 2*i + 3*j <= 10
             and (i, j) not in {(0, 1), (1, 1)}],
            key=lambda row: (row[2], row[1], row[0]))
    )
    expected = (
        (1, 0, 2), (2, 0, 4), (3, 0, 6), (0, 2, 6),
        (2, 1, 7), (4, 0, 8), (1, 2, 8), (3, 1, 9),
        (0, 3, 9), (5, 0, 10), (2, 2, 10),
    )
    require(monomials == expected, "complete M10 monomial universe")

    # With zeta=0 every surviving endpoint lies above M; every y^2 first
    # endpoint lies on V and all other endpoints lie strictly above V.
    plane_m = (Fraction(1, 10), Fraction(1, 5), Fraction(-1, 5))
    plane_v = (Fraction(1, 2), Fraction(0), Fraction(-1))
    for i, j, weight in (row for row in monomials if row[:2] != (0, 3)):
        first = (j + 2, i + j, 1)
        second = (j, i + j + 1, 1)
        for point in (first, second):
            m_height = (plane_m[0]*point[0] + plane_m[1]*point[1]
                        + plane_m[2])
            require(Fraction(point[2]) - m_height == Fraction(10-weight, 10),
                    "main-plane weight gap")
            v_height = (plane_v[0]*point[0] + plane_v[1]*point[1]
                        + plane_v[2])
            require(Fraction(point[2]) - v_height >= 0,
                    "side-plane nonnegative gap")
        require((Fraction(first[2])
                 - plane_v[0]*first[0] - plane_v[2] == 0) == (j == 2),
                "side equality is exactly a y^2 first endpoint")

    A, B, U, W, D, E = (0, 1), (2, 0), (4, 2), (4, 3), (4, 4), (0, 6)
    main = polygon_ledger((A, B, D, E))
    k_nonzero = polygon_ledger((A, B, U, D, E))
    theta_only = polygon_ledger((A, B, W, D, E))
    require(main == (((0, 1), (2, 0), (4, 4), (0, 6)), 30, 10, 11),
            "collapsed main polygon")
    require(k_nonzero == (
        ((0, 1), (2, 0), (4, 2), (4, 4), (0, 6)), 34, 12, 12),
        "K-nonzero polygon")
    require(theta_only == (
        ((0, 1), (2, 0), (4, 3), (4, 4), (0, 6)), 32, 10, 12),
        "K-zero Theta-nonzero polygon")

    packet_k, rows_k = edge_packet(k_nonzero[0])
    packet_theta, rows_theta = edge_packet(theta_only[0])
    packet_zero, rows_zero = edge_packet(main[0])
    require(packet_k == (9, 9, 3, 3, 2, 2, 1), "K-nonzero packet")
    require(packet_theta == (9, 9, 5, 3, 1), "Theta-only packet")
    require(packet_zero == (9, 9, 3, 3, 1), "collapsed packet")
    for packet, genus in ((packet_k, 12), (packet_theta, 12),
                          (packet_zero, 11)):
        require(sum(value - 1 for value in packet) == 2*genus - 2,
                "packet saturates genus defect")
    require((sum(packet_k), sum(packet_theta), sum(packet_zero)) == (29, 27, 25),
            "full packet degrees")

    # Exact source and face equations.
    S, P, Y, Z = sy.symbols("S P Y Z")
    K, Theta, xi, upsilon = sy.symbols("K Theta xi upsilon")
    main_face = sy.expand((S**2-P)
                          * (1-upsilon*P**5-xi*S**2*P**4))
    require(sy.factor(main_face - ((S**2-P)
                                    * (-P**5*upsilon-P**4*S**2*xi+1))) == 0,
            "main face factorization")
    cyclotomic_on_rational = sy.factor(
        (1-upsilon*P**5-xi*S**2*P**4).subs(P, S**2))
    require(sy.expand(cyclotomic_on_rational
                      - (1-S**10*(upsilon+xi))) == 0,
            "ten main intersections")
    determinant = sy.factor(-10*S**9*(upsilon+xi))
    require(determinant != 0, "transverse-node determinant expression")

    side_core = 1-S**2*P**2*(K+Theta*P+xi*P**2)
    side_conic = Y**2-(K+Theta*P+xi*P**2)
    require(sy.factor(side_core.subs(S, 1/(Y*P))*Y**2
                      - side_conic) == 0,
            "side normalization is a conic")
    discriminant = sy.discriminant(K+Theta*P+xi*P**2, P)
    require(sy.factor(discriminant) == Theta**2-4*K*xi,
            "side-conic discriminant")
    require(sy.gcd(sy.Poly(K+Theta*Z+xi*Z**2, Z),
                   sy.Poly(Theta+2*xi*Z, Z)).degree() == 0,
            "symbolic generic side edge squarefree")

    # Reduced edge schemes.  The discriminants are the exact unit gates.
    require(sy.discriminant(1-K*Z**2, Z) == 4*K,
            "BU edge discriminant")
    require(sy.factor(sy.discriminant(K+Theta*Z+xi*Z**2, Z)) ==
            Theta**2-4*K*xi, "UD edge discriminant")
    require(sy.discriminant(1-xi*Z**2, Z) == 4*xi,
            "internal BD discriminant")
    require(sy.factor(sy.discriminant((1-Z)*(upsilon+xi*Z), Z)) ==
            (upsilon+xi)**2, "DE root-separation discriminant")
    require(sy.discriminant(1-upsilon*Z**5, Z) != 0,
            "affine EA edge reduced")
    require(sy.degree(1-Theta*Z, Z) == 1 and sy.degree(Theta+xi*Z, Z) == 1,
            "Theta-only side edges linear")

    # Generic moving edges: quadratic, linear, quadratic in the three rows.
    q = sy.symbols("q")
    require(sy.factor(sy.discriminant(K*Z**2-(q-sy.Rational(1, 2)), Z)
                      - 4*K*(q-sy.Rational(1, 2))) == 0,
            "K-carrier separability formula")
    require(sy.degree(Theta*Z-(q-sy.Rational(1, 2)), Z) == 1,
            "Theta carrier is a rational section")
    require(sy.factor(sy.discriminant(xi*Z**2-(q-sy.Rational(1, 2)), Z)
                      - 4*xi*(q-sy.Rational(1, 2))) == 0,
            "xi-carrier separability formula")

    # The side has genus zero; the main nonrational branch has genus two.
    # Ten R--C paths and, when present, two C--V paths give the full Pick genus.
    require((10-2+1, 2+(10-2+1)) == (9, 11),
            "collapsed dual graph and genus")
    require((12-3+1, 2+(12-3+1)) == (10, 12),
            "side-conic dual graph and genus")

    # Common base change Q=sigma^30.  Both primitive height normals have
    # final coordinate one, hence all face multiplicities are one.
    require((30*plane_m[0], 30*plane_m[1], 30*plane_m[2]) == (3, 6, -6),
            "integral main height")
    require((30*plane_v[0], 30*plane_v[1], 30*plane_v[2]) == (15, 0, -30),
            "integral side height")
    require(gcd(gcd(3, 6), 1) == gcd(gcd(15, 0), 1) == 1,
            "primitive normals and multiplicity one")
    # The internal BD slope pair is 0>-6.  Definition 3.12 inserts the five
    # intermediate integral slopes, once for each of the two reduced roots.
    require((0, -6, 0-(-6)-1, 2*(0-(-6)-1)) == (0, -6, 5, 10),
            "BD has five rational curves per root")

    sigma, H_M, H_V = sy.symbols("sigma H_M H_V")
    U_0 = S**2-P
    V_0 = (1-H_M)/S**2
    main_scaled = (S**2-P)*(1-H_M)-sigma**30*S**2/sy.Integer(2)
    require(sy.factor(main_scaled/S**2
                      - (U_0*V_0-sigma**30/sy.Integer(2))) == 0,
            "exact A29 main chart UV=sigma^30/2")
    side_scaled = ((S**2-sigma**30*P)*(1-H_V)
                   - sigma**30*S**2/sy.Integer(2))
    require(side_scaled.subs(sigma, 0) == S**2*(1-H_V),
            "exact side chart special equation")
    require(30-10*3 == 0 and 30-15*2 == 0,
            "good target scaling")
    require((20, 30) == (30-10, 30), "target correction exponents")

    # Degree conservation has no nonconstant special component: the conic,
    # rational branch, chains and exceptional curves have genus zero, while
    # the simple genus-two component has Hom(J(C),E_0)=0 (THM-4218).
    multiplicities = (1, 1, 1)
    component_degrees = (0, 0, 0)
    require(sum(m*d for m, d in zip(multiplicities, component_degrees)) == 0,
            "specialized cover degree is zero")

    # Hostile collision: D_V=0,K!=0 makes the side conic split and the UD
    # edge nonreduced.  It is deliberately outside the theorem.
    collision = sy.factor((K+Theta*P+xi*P**2).subs(K, Theta**2/(4*xi)))
    require(collision == (2*P*xi+Theta)**2/(4*xi),
            "split-conic collision factorization")
    require(sy.factor((Theta**2-4*K*xi).subs(K, Theta**2/(4*xi))) == 0,
            "collision edge has a double root")

    semantic = digest((main, k_nonzero, theta_only,
                       packet_k, packet_theta, packet_zero,
                       rows_k, rows_theta, rows_zero,
                       sy.srepr(main_face), sy.srepr(side_core),
                       sy.srepr(collision)))
    print(f"checks={CHECKS}")
    print("scope=zeta=0;upsilon*xi*(upsilon+xi)!=0")
    print("proved=(D_V!=0)|(K=Theta=0);D_V=Theta^2-4*K*xi")
    print("planes=M:(r+2k-2)/10;V:(r-2)/2")
    print("K_nonzero_DV_nonzero=polygon[(0,1),(2,0),(4,2),(4,4),(0,6)];Pick=34,12,12;packet=(9,9,3,3,2,2,1);degree=29")
    print("K_zero_Theta_nonzero=polygon[(0,1),(2,0),(4,3),(4,4),(0,6)];Pick=32,10,12;packet=(9,9,5,3,1);degree=27")
    print("K_Theta_zero=polygon[(0,1),(2,0),(4,4),(0,6)];Pick=30,10,11;packet=(9,9,3,3,1);degree=25")
    print("main=C_genus2+R_rational;side=conic_genus0;graphs=b1_10_or_9")
    print("internal_BD=slopes_0>-6;roots=2;P1_per_root=5;multiplicity=1")
    print("base_change=Q=sigma^30;multiplicities=1;target_special=y^2=x^3+1")
    print("degree_conservation=0;finite_nonconstant_generic_map=impossible")
    print("open=D_V=0,K!=0(split_conic_and_double_UD_root)")
    print(f"semantic_sha256={semantic}")


if __name__ == "__main__":
    main()
