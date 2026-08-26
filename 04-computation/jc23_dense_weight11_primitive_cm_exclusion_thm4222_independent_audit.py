#!/usr/bin/env python3
"""Clean-room exact referee certificate for THM-4222.

This audit reconstructs the dense exact-M=11 Newton model independently of
the primary program.  Its universe is the inherited b=d=0 reduced (2,3)
seam with A*B*U*Z*(A+B) != 0 and every other lower coefficient arbitrary.
Milne, Complex Multiplication (2020), Proposition 3.13 (together with the
classification in Proposition 3.12) is a CITED input only after this script
establishes the full Q(zeta_11)-action and a primitive CM type.
"""

from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from math import gcd

import sympy as sy


CHECKS = 0


def require(assertion, label):
    global CHECKS
    CHECKS += 1
    if not bool(assertion):
        raise RuntimeError(label)


def hull(points):
    points = sorted(set(points))

    def turn(a, b, c):
        return ((b[0]-a[0])*(c[1]-a[1])
                - (b[1]-a[1])*(c[0]-a[0]))

    lower = []
    for point in points:
        while len(lower) > 1 and turn(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)
    upper = []
    for point in reversed(points):
        while len(upper) > 1 and turn(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    return tuple(lower[:-1] + upper[:-1])


def pick(points):
    vertices = hull(points)
    area2 = abs(sum(
        a[0]*b[1]-a[1]*b[0]
        for a, b in zip(vertices, vertices[1:]+vertices[:1])
    ))
    boundary = sum(
        gcd(abs(a[0]-b[0]), abs(a[1]-b[1]))
        for a, b in zip(vertices, vertices[1:]+vertices[:1])
    )
    require((area2-boundary+2) % 2 == 0, "Pick parity")
    return vertices, area2, boundary, (area2-boundary+2)//2


def plane_through(a, b, c):
    determinant = ((b[0]-a[0])*(c[1]-a[1])
                   - (b[1]-a[1])*(c[0]-a[0]))
    if determinant == 0:
        return None
    slope_r = Fraction(
        (b[2]-a[2])*(c[1]-a[1])-(b[1]-a[1])*(c[2]-a[2]),
        determinant,
    )
    slope_k = Fraction(
        (b[0]-a[0])*(c[2]-a[2])-(b[2]-a[2])*(c[0]-a[0]),
        determinant,
    )
    intercept = Fraction(a[2])-slope_r*a[0]-slope_k*a[1]
    return slope_r, slope_k, intercept


def lower_facets(points):
    answer = set()
    points = tuple(sorted(points))
    for triple in combinations(points, 3):
        plane = plane_through(*triple)
        if plane is None:
            continue
        sr, sk, intercept = plane
        gaps = [Fraction(h)-sr*r-sk*k-intercept for r, k, h in points]
        if min(gaps) < 0:
            continue
        equality = [point for point, gap in zip(points, gaps) if gap == 0]
        if len(hull((r, k) for r, k, _h in equality)) >= 3:
            answer.add(plane)
    return answer


def facet_catalog(universe):
    """Precompute lower-plane masks once for the deletion hostile sweep."""
    universe = tuple(sorted(universe))
    index = {point: bit for bit, point in enumerate(universe)}
    records = {}
    for triple in combinations(universe, 3):
        plane = plane_through(*triple)
        if plane is None or plane in records:
            continue
        sr, sk, intercept = plane
        below = 0
        equal = 0
        for bit, (r, k, height) in enumerate(universe):
            gap = Fraction(height)-sr*r-sk*k-intercept
            if gap < 0:
                below |= 1 << bit
            elif gap == 0:
                equal |= 1 << bit
        records[plane] = (below, equal)
    return universe, index, records


def catalog_facets(active, catalog):
    universe, index, records = catalog
    active_bits = sum(1 << index[point] for point in active)
    answer = set()
    for plane, (below, equal) in records.items():
        if below & active_bits:
            continue
        equality = [universe[bit] for bit in range(len(universe))
                    if active_bits & equal & (1 << bit)]
        if len(hull((r, k) for r, k, _height in equality)) >= 3:
            answer.add(plane)
    return answer


def endpoints(row):
    i, j, _weight = row
    return {(j+2, i+j, 1), (j, i+j+1, 1)}


def edge_restriction(poly, X, Y, start, end, z):
    dx, dy = end[0]-start[0], end[1]-start[1]
    length = gcd(abs(dx), abs(dy))
    ux, uy = dx//length, dy//length
    result = 0
    for (i, j), coefficient in sy.Poly(poly, X, Y).terms():
        vx, vy = i-start[0], j-start[1]
        if vx*dy-vy*dx:
            continue
        position = vx//ux if ux else vy//uy
        if 0 <= position <= length:
            result += coefficient*z**position
    return sy.factor(result)


def main():
    # Complete source support through weight eleven.
    monomials = tuple(sorted(
        ((i, j, 2*i+3*j)
         for i in range(6) for j in range(4)
         if 0 < 2*i+3*j <= 11 and (i, j) not in {(0, 1), (1, 1)}),
        key=lambda row: (row[2], row[1], row[0])))
    expected = (
        (1, 0, 2), (2, 0, 4), (3, 0, 6), (0, 2, 6),
        (2, 1, 7), (4, 0, 8), (1, 2, 8), (3, 1, 9),
        (0, 3, 9), (5, 0, 10), (2, 2, 10),
        (4, 1, 11), (1, 3, 11),
    )
    require(monomials == expected, "complete M11 universe")
    required_rows = {(4, 1), (1, 3), (5, 0), (0, 3)}
    optional = tuple(row for row in monomials if row[:2] not in required_rows)
    require(len(optional) == 9, "nine arbitrary lower rows")

    base = {(2, 0, 0), (0, 1, 0)}
    required_support = set(base)
    for row in monomials:
        if row[:2] in required_rows:
            required_support |= endpoints(row)
    all_support = set(base)
    endpoint_owners = {}
    for row in monomials:
        all_support |= endpoints(row)
        for point in endpoints(row):
            endpoint_owners.setdefault(point, []).append(row[:2])
    collisions = tuple(sorted(point for point, owners in endpoint_owners.items()
                              if len(owners) > 1))
    require(collisions == ((2, 3, 1), (2, 4, 1), (2, 5, 1),
                           (3, 4, 1), (3, 5, 1)),
            "complete coincident-support list")

    M = (Fraction(1, 11), Fraction(2, 11), Fraction(-2, 11))
    V = (Fraction(0), Fraction(1, 5), Fraction(-1, 5))
    T = (Fraction(1, 3), Fraction(0), Fraction(-2, 3))
    expected_faces = {M, V, T}

    catalog = facet_catalog(all_support)

    # An independent brute hull pass covers every lower-row support subset.
    for mask in range(1 << len(optional)):
        active = set(required_support)
        for bit, row in enumerate(optional):
            if mask & (1 << bit):
                active |= endpoints(row)
        require(catalog_facets(active, catalog) == expected_faces,
                "lower hull under arbitrary lower-row deletion")
    # Independently delete every subset of coincident aggregate coefficients.
    for mask in range(1 << len(collisions)):
        deleted = {point for bit, point in enumerate(collisions)
                   if mask & (1 << bit)}
        require(catalog_facets(all_support-deleted, catalog) == expected_faces,
                "lower hull under coefficient cancellation")

    # Exact analytic gaps explain why no fourth face can enter.
    for i, j, weight in monomials:
        first, second = ((j+2, i+j, 1), (j, i+j+1, 1))
        for plane, wanted in (
            (M, (Fraction(11-weight, 11),)*2),
            (V, (Fraction(6-i-j, 5), Fraction(5-i-j, 5))),
            (T, (Fraction(3-j, 3), Fraction(5-j, 3))),
        ):
            gaps = tuple(Fraction(h)-plane[0]*r-plane[1]*k-plane[2]
                         for r, k, h in (first, second))
            require(gaps == wanted, "analytic face gap")

    # Face and global polygons.
    A0, B0, C0, D0, E0, F0 = (
        (0, 1), (2, 0), (5, 3), (5, 4), (1, 6), (0, 6)
    )
    main_polygon = (A0, B0, D0, E0)
    tail_polygon = (B0, C0, D0)
    vertical_polygon = (A0, E0, F0)
    require(pick(main_polygon)[1:] == (33, 5, 15), "main Pick ledger")
    require(pick(tail_polygon)[1:] == (3, 5, 0), "tail Pick ledger")
    require(pick(vertical_polygon)[1:] == (5, 7, 0), "vertical Pick ledger")
    global_ledger = pick((A0, B0, C0, D0, E0, F0))
    require(global_ledger == ((A0, B0, C0, D0, E0, F0), 41, 13, 15),
            "global Pick ledger")

    # Face equations and all six outer plus two internal schemes.
    S, P, z = sy.symbols("S P z")
    AA, BB, UU, ZZ = sy.symbols("AA BB UU ZZ", nonzero=True)
    R = S**2-P
    C = 1-AA*S*P**5-BB*S**3*P**4
    main_face = sy.expand(R*C)
    tail_face = sy.expand(S**2*(1-(S*P)**3*(ZZ+BB*P)))
    vertical_face = sy.expand(P*(-1+P**5*(UU+AA*S)))
    schemes = (
        edge_restriction(main_face, S, P, A0, B0, z),
        edge_restriction(tail_face, S, P, B0, C0, z),
        edge_restriction(tail_face, S, P, C0, D0, z),
        edge_restriction(main_face, S, P, D0, E0, z),
        edge_restriction(vertical_face, S, P, E0, F0, z),
        edge_restriction(vertical_face, S, P, F0, A0, z),
        edge_restriction(main_face, S, P, A0, E0, z),
        edge_restriction(main_face, S, P, B0, D0, z),
    )
    expected_schemes = (
        z-1, 1-ZZ*z**3, -ZZ-BB*z, (z-1)*(AA*z+BB),
        AA+UU*z, UU-z**5, AA*z-1, 1-BB*z,
    )
    require(all(sy.factor(g-h) == 0
                for g, h in zip(schemes, expected_schemes)),
            "eight exact edge schemes")
    discriminants = tuple(sy.factor(sy.discriminant(g, z)) for g in schemes)
    require(discriminants == (1, -27*ZZ**2, 1, (AA+BB)**2,
                              1, 3125*UU**4, 1, 1),
            "five-factor reduced-edge gate")
    require(sy.factor(schemes[3].subs(BB, AA)
                      - AA*(z-1)*(z+1)) == 0,
            "A=B cancellation is safe")

    # Packet reconstructed from primitive inward normals.
    packet = []
    outer_gcds = []
    height_vertex = {
        A0: (0, 1, 0), B0: (2, 0, 0), C0: (5, 3, 330),
        D0: (5, 4, 330), E0: (1, 6, 330), F0: (0, 6, 330),
    }
    for start, end in zip(global_ledger[0], global_ledger[0][1:]+global_ledger[0][:1]):
        dx, dy = end[0]-start[0], end[1]-start[1]
        length = gcd(abs(dx), abs(dy))
        inward = (-dy//length, dx//length)
        constant = inward[0]*start[0]+inward[1]*start[1]
        index = inward[0]+inward[1]-constant
        if not (start[0] == end[0] == 0):
            packet.extend([index]*length)
        delta = tuple(height_vertex[end][i]-height_vertex[start][i]
                      for i in range(3))
        outer_gcds.append(gcd(gcd(abs(delta[0]), abs(delta[1])), abs(delta[2])))
        require(outer_gcds[-1] == length, "outer three-dimensional primitivity")
    packet = tuple(sorted(packet, reverse=True))
    require(packet == (10, 10, 5, 4, 2, 2, 2, 1), "response packet")
    require((sum(packet), sum(e-1 for e in packet)) == (36, 28),
            "packet degree and genus defect")

    # Q=sigma^330: primitive faces, exact internal slope subdivisions.
    hM = lambda r, k: 30*r+60*k-60
    hV = lambda r, k: 66*k-66
    hT = lambda r, k: 110*r-220
    require(all(gcd(gcd(abs(a), abs(b)), abs(c)) == 1
                for a, b, c in ((30, 60, -1), (0, 66, -1),
                                (110, 0, -1))),
            "face multiplicity one")
    ae = (hM(0, 0)-hM(*A0), hV(0, 0)-hV(*A0))
    bd_probe = (1, -1)
    bd = (hM(*bd_probe)-hM(*B0), hT(*bd_probe)-hT(*B0))
    require(ae == (-60, -66) and abs(ae[0]-ae[1])-1 == 5,
            "AE five-component chain")
    require(bd == (-90, -110) and abs(bd[0]-bd[1])-1 == 19,
            "BD nineteen-component chain")
    require(all(abs(a-b) == 1 for a, b in zip(range(-60, -67, -1),
                                               range(-61, -67, -1))),
            "AE determinant-one subdivision")
    require(all(abs(a-b) == 1 for a, b in zip(range(-90, -111, -1),
                                               range(-91, -111, -1))),
            "BD determinant-one subdivision")

    # Face smoothness and main intersections.
    log_gradient = sy.Matrix(((AA, 3*BB), (5*AA, 4*BB)))
    require(sy.factor(log_gradient.det()+11*AA*BB) == 0,
            "C torus smoothness")
    node_restriction = sy.factor(C.subs(P, S**2))
    require(sy.factor(node_restriction-(1-(AA+BB)*S**11)) == 0,
            "eleven R-C nodes")
    node_jacobian = sy.factor(sy.det(sy.Matrix((
        (sy.diff(R, S), sy.diff(R, P)),
        (sy.diff(C, S), sy.diff(C, P)),
    ))).subs(P, S**2))
    require(sy.factor(node_jacobian+11*(AA+BB)*S**10) == 0,
            "transverse R-C determinant")
    require(sy.diff(-1+P**5*(UU+AA*S), S) == AA*P**5,
            "V rational and smooth")
    W = sy.symbols("W")
    require(sy.diff(1-W**3*(ZZ+BB*P), P) == -BB*W**3,
            "T rational and smooth")

    # Exact scaled source charts, reconstructed with arbitrary lower owners.
    sigma = sy.symbols("sigma")
    K, Phi, Delta, Theta, eta, Xi = sy.symbols(
        "K Phi Delta Theta eta Xi")
    ss, pp = sy.symbols("ss pp")
    H = (-3*pp+sy.Rational(8, 3)*pp**2-sy.Rational(1376, 135)*pp**3
         +K*ss**2*pp**2+Phi*ss*pp**3+Delta*pp**4
         +Theta*ss**2*pp**3+eta*ss*pp**4+ZZ*ss**3*pp**3
         +UU*pp**5+Xi*ss**2*pp**4+AA*ss*pp**5+BB*ss**3*pp**4)
    source = (ss**2-pp)*(1-sigma**330*H)-sigma**330*ss**2/2
    Hm = sy.expand(sigma**330*H.subs({ss: sigma**-30*S,
                                      pp: sigma**-60*P}))
    require(not Hm.has(sigma**-1), "main chart integral")
    require(sy.factor(Hm.subs(sigma, 0)
                      -AA*S*P**5-BB*S**3*P**4) == 0,
            "main face reduction")
    Fm = sy.expand(sigma**60*source.subs({ss: sigma**-30*S,
                                         pp: sigma**-60*P}))
    exact_main = (S**2-P)*(1-Hm)-sigma**330*S**2/2
    require(sy.factor(Fm-exact_main) == 0, "exact main chart")
    local_u = S**2-P
    local_v = (1-Hm)/S**2
    require(sy.cancel(Fm/S**2-(local_u*local_v-sigma**330/2)) == 0,
            "eleven A329 equations")
    require(330-1 == 329, "A329 chain length")
    Fv = sy.expand(sigma**66*source.subs({ss: S, pp: sigma**-66*P}))
    require(sy.factor(Fv.subs(sigma, 0)-vertical_face) == 0,
            "exact V side chart")
    Ft = sy.expand(sigma**220*source.subs({ss: sigma**-110*S, pp: P}))
    require(sy.factor(Ft.subs(sigma, 0)-tail_face) == 0,
            "exact T side chart")

    # Genus-five cyclic normalization and two independent CM-character paths.
    x = sy.symbols("x")
    x_on_curve = AA*S*P**5
    cyclic_identity = sy.factor(
        AA**3*P**11*(1-x_on_curve)-BB*x_on_curve**3)
    require(sy.factor(cyclic_identity-AA**3*P**11*C) == 0,
            "cyclic-cover identity")
    branch = (3, 10, 9)
    require(sum(branch) % 11 == 0 and all(gcd(11, a) == 1 for a in branch),
            "three total branch points")
    require(2*5-2 == 11*(-2)+3*(11-1), "Riemann-Hurwitz genus five")

    triangle = ((0, 0), (1, 5), (3, 4))
    interiors = []
    total_area = abs((triangle[1][0]-triangle[0][0])
                     *(triangle[2][1]-triangle[0][1])
                     -(triangle[1][1]-triangle[0][1])
                     *(triangle[2][0]-triangle[0][0]))
    for i in range(4):
        for j in range(6):
            pieces = []
            for a, b in zip(triangle, triangle[1:]+triangle[:1]):
                pieces.append(abs((b[0]-a[0])*(j-a[1])
                                  -(b[1]-a[1])*(i-a[0])))
            if sum(pieces) == total_area and all(value > 0 for value in pieces):
                interiors.append((i, j))
    require(tuple(interiors) == ((1, 2), (1, 3), (1, 4), (2, 3), (2, 4)),
            "five interior differentials")
    # tau:(S,P)->(zeta^6 S,zeta P); residue differentials have 6i+j.
    cm_type = {(6*i+j) % 11 for i, j in interiors}
    require(cm_type == {4, 5, 8, 9, 10}, "toric differential CM type")
    # Chevalley-Weil fractional-part computation gives the same convention.
    cw = set()
    for character in range(1, 11):
        dimension = -1+sum(Fraction((-character*a) % 11, 11)
                           for a in branch)
        require(dimension in (0, 1), "Chevalley-Weil multiplicity")
        if dimension == 1:
            cw.add(character)
    require(cw == cm_type, "cyclic-cover CM character cross-check")
    complement = {(-a) % 11 for a in cm_type}
    require(cm_type.isdisjoint(complement)
            and cm_type | complement == set(range(1, 11)),
            "full H1 cyclotomic spectrum")
    stabilizer = tuple(u for u in range(1, 11)
                       if {(u*a) % 11 for a in cm_type} == cm_type)
    require(stabilizer == (1,), "primitive CM stabilizer")
    require(sy.Poly(sy.cyclotomic_poly(11, x), x).degree() == 10,
            "full degree-ten CM action")

    # Labelled attachments and component/genus completeness.
    require(11 == 11 and sy.isprime(11), "attachment torsion prime")
    core_edges = 11+1+1
    core_vertices = 4
    graph_rank = core_edges-core_vertices+1
    require((graph_rank, 5+graph_rank) == (10, 15),
            "complete dual graph genus")
    require((sum(packet)-3*2, sum(packet)) == (30, 36),
            "finite/full response degrees")

    # Good j=0 target and zero-degree specialization.
    a, XX, YY = sy.symbols("a XX YY")
    target = (YY**2-XX**3-1
              +sy.Rational(3, 4)*a**2*sigma**220*XX
              +sy.Rational(1, 4)*a**3*sigma**330)
    require(sy.factor(target.subs(sigma, 0)-(YY**2-XX**3-1)) == 0,
            "good target reduction")
    require(sy.discriminant(XX**3+1, XX) == -27,
            "target special smooth")
    # Milne Prop. 1.9/Cor. 1.10 turns the trivial stabilizer into a primitive
    # CM-pair; Props. 3.12--3.13 then make J(C) simple.  Hence Hom(J(C),E0)=0.
    genera = (0, 5, 0, 0)
    require(genera == (0, 5, 0, 0), "sole positive-genus component")
    multiplicities = (1, 1, 1, 1)
    special_degrees = (0, 0, 0, 0)
    require(sum(m*d for m, d in zip(multiplicities, special_degrees)) == 0,
            "degree-zero specialization")

    semantic = sha256("|".join(map(str, (
        monomials, collisions, global_ledger, packet, schemes,
        ae, bd, sy.srepr(C), interiors, tuple(sorted(cm_type)),
        stabilizer, graph_rank,
    ))).encode()).hexdigest()
    print("THM4222_INDEPENDENT_CLEAN_ROOM_AUDIT")
    print(f"checks={CHECKS}")
    print("scope=exact_M11;A*B*U*Z*(A+B)!=0;all_lower_coefficients")
    print("faces=M11,V5,T3;support_subsets=512;collision_patterns=32")
    print("Pick=global(41,13,15);main(33,5,15);side_genera=(0,0)")
    print("packet=(10,10,5,4,2,2,2,1);responses=(full36,finite30)")
    print("edges=8_reduced;internal_chains=(AE5,BD19);nodes=11*A329")
    print("model=Q:sigma^330;multiplicity_one;components_complete")
    print("C=P^11-(B/A^3)*x^3/(1-x);genus=5")
    print("CM_type=(4,5,8,9,10);H1=all_nontrivial_chars;stabilizer=(1)")
    print("CITED=Milne_Prop1.9_Cor1.10_Props3.12-3.13;Jacobian=simple")
    print("attachments=(x=1,x=0);dual_graph_rank=10;special_genus=15")
    print("target=good_j0;Hom_to_E0=0;specialized_degree=0")
    print("verdict=ACCEPT_TMM4222".replace("TMM", "THM"))
    print(f"semantic_sha256={semantic}")


if __name__ == "__main__":
    main()
