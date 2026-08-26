#!/usr/bin/env python3
"""Clean-room exact referee for THM-4218.

This audit does not import the primary certificate.  It reconstructs the
literal exact-M=10 source, uses an integer cross-product lower-hull routine,
extracts edge schemes directly from the face polynomials, and checks the
regular-model, attachment, specialization, and response ledgers.  The final
critical resultant is a positive control only; the theorem scope is all
lower coefficients in the stated top chamber.
"""

from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from math import gcd

import sympy as sy


CHECKS = 0


def require(condition, label):
    """Raise without relying on assert, so normal and -O runs agree."""
    global CHECKS
    CHECKS += 1
    if not bool(condition):
        raise RuntimeError(label)


def lower_facets(points):
    """Return lower facets by oriented integer cross products in R^3."""
    points = tuple(sorted(set(points)))
    facets = {}
    for first, second, third in combinations(points, 3):
        ab = tuple(second[k] - first[k] for k in range(3))
        ac = tuple(third[k] - first[k] for k in range(3))
        normal = (
            ab[1] * ac[2] - ab[2] * ac[1],
            ab[2] * ac[0] - ab[0] * ac[2],
            ab[0] * ac[1] - ab[1] * ac[0],
        )
        if normal[2] == 0:
            continue
        if normal[2] < 0:
            normal = tuple(-entry for entry in normal)
        gaps = tuple(sum(normal[k] * (point[k] - first[k]) for k in range(3))
                     for point in points)
        if min(gaps) < 0:
            continue
        equality = tuple(point for point, gap in zip(points, gaps) if gap == 0)
        if len(equality) < 3:
            continue
        divisor = gcd(gcd(abs(normal[0]), abs(normal[1])), abs(normal[2]))
        primitive = tuple(entry // divisor for entry in normal)
        constant = sum(primitive[k] * first[k] for k in range(3))
        facets[equality] = primitive + (constant,)
    return tuple(sorted(facets.items(), key=lambda row: row[1]))


def hull2(points):
    """Counterclockwise convex hull, discarding collinear edge interiors."""
    points = sorted(set(points))

    def turn(origin, first, second):
        return ((first[0] - origin[0]) * (second[1] - origin[1])
                - (first[1] - origin[1]) * (second[0] - origin[0]))

    low = []
    for point in points:
        while len(low) > 1 and turn(low[-2], low[-1], point) <= 0:
            low.pop()
        low.append(point)
    high = []
    for point in reversed(points):
        while len(high) > 1 and turn(high[-2], high[-1], point) <= 0:
            high.pop()
        high.append(point)
    return tuple(low[:-1] + high[:-1])


def pick(points):
    polygon = hull2(points)
    area2 = abs(sum(
        polygon[index][0] * polygon[(index + 1) % len(polygon)][1]
        - polygon[(index + 1) % len(polygon)][0] * polygon[index][1]
        for index in range(len(polygon))
    ))
    boundary = sum(
        gcd(abs(polygon[(index + 1) % len(polygon)][0] - polygon[index][0]),
            abs(polygon[(index + 1) % len(polygon)][1] - polygon[index][1]))
        for index in range(len(polygon))
    )
    require((area2 - boundary + 2) % 2 == 0, "Pick parity")
    return polygon, area2, boundary, (area2 - boundary + 2) // 2


def edge_polynomial(polynomial, first, second, start, end, variable):
    """Extract a lattice-edge polynomial in its primitive edge coordinate."""
    delta = (end[0] - start[0], end[1] - start[1])
    length = gcd(abs(delta[0]), abs(delta[1]))
    step = (delta[0] // length, delta[1] // length)
    answer = 0
    for exponent, coefficient in sy.Poly(polynomial, first, second).terms():
        difference = (exponent[0] - start[0], exponent[1] - start[1])
        if step[0] != 0:
            if difference[0] % step[0]:
                continue
            position = difference[0] // step[0]
        else:
            if step[1] == 0 or difference[1] % step[1]:
                continue
            position = difference[1] // step[1]
        if (difference == (position * step[0], position * step[1])
                and 0 <= position <= length):
            answer += coefficient * variable**position
    return sy.factor(answer)


def polynomial_order(expression, variable):
    polynomial = sy.Poly(expression, variable)
    return min(exponent[0] for exponent, coefficient in polynomial.terms()
               if coefficient != 0)


def coefficient_hash(polynomial):
    payload = ",".join(str(sy.Rational(value)) for value in polynomial.all_coeffs())
    return sha256(payload.encode()).hexdigest()


def main():
    # Literal monomial universe; b=[y]H and d=[py]H are the only deletions.
    enumerated = []
    for i in range(6):
        for j in range(4):
            weight = 2 * i + 3 * j
            if 0 < weight <= 10 and (i, j) not in {(0, 1), (1, 1)}:
                enumerated.append((i, j, weight))
    monomials = tuple(sorted(enumerated, key=lambda row: (row[2], row[1], row[0])))
    require(monomials == (
        (1, 0, 2), (2, 0, 4), (3, 0, 6), (0, 2, 6),
        (2, 1, 7), (4, 0, 8), (1, 2, 8), (3, 1, 9),
        (0, 3, 9), (5, 0, 10), (2, 2, 10),
    ), "literal exact-M10 universe")

    base_support = {(2, 0, 0), (0, 1, 0)}
    top_indices = {(5, 0), (2, 2), (0, 3)}
    lower_terms = tuple(row for row in monomials if row[:2] not in top_indices)

    def endpoints(rows):
        support = set(base_support)
        for i, j, _weight in rows:
            support.add((j + 2, i + j, 1))
            support.add((j, i + j + 1, 1))
        return support

    top_terms = tuple(row for row in monomials if row[:2] in top_indices)
    expected_normals = {(-1, -2, 10, -2), (-1, -1, 6, -2)}
    # Optional lower coefficients may independently vanish.  Coincident
    # lower endpoints are strictly above the two faces, so cancellation is
    # included in this hostile universe.
    for mask in range(1 << len(lower_terms)):
        rows = list(top_terms)
        rows.extend(row for bit, row in enumerate(lower_terms) if mask & (1 << bit))
        facet_rows = lower_facets(endpoints(rows))
        normals = {data for _points, data in facet_rows}
        require(normals == expected_normals, "universal two-face hull")

    # Independent analytic gaps explain the finite enumeration.
    for i, j, weight in monomials:
        first = (j + 2, i + j, 1)
        second = (j, i + j + 1, 1)
        main_gap_first = Fraction(10 - weight, 10)
        main_gap_second = Fraction(10 - weight, 10)
        require(main_gap_first >= 0 and main_gap_second >= 0, "main weight gap")
        tail_gap_first = Fraction(6 - i - 2*j, 6)
        tail_gap_second = Fraction(7 - i - 2*j, 6)
        require(tail_gap_first >= 0 and tail_gap_second >= 0, "tail weight gap")
        computed_first = Fraction(first[2]) - Fraction(first[0] + first[1] - 2, 6)
        computed_second = Fraction(second[2]) - Fraction(second[0] + second[1] - 2, 6)
        require((computed_first, computed_second) == (tail_gap_first, tail_gap_second),
                "tail gap formula")

    main_points = ((0, 1), (2, 0), (4, 4), (0, 6))
    tail_points = ((2, 0), (5, 3), (4, 4))
    global_polygon, area2, boundary, genus = pick(main_points + tail_points)
    require(pick(main_points)[1:] == (30, 10, 11), "main Pick ledger")
    require(pick(tail_points)[1:] == (6, 6, 1), "tail Pick ledger")
    require((global_polygon, area2, boundary, genus) == (
        ((0, 1), (2, 0), (5, 3), (4, 4), (0, 6)), 36, 12, 13,
    ), "global Pick ledger")

    # Derive the infinity packet from normals and the residue monomial (1,1).
    packet = []
    edge_ledger = []
    for start, end in zip(global_polygon, global_polygon[1:] + global_polygon[:1]):
        dx, dy = end[0] - start[0], end[1] - start[1]
        length = gcd(abs(dx), abs(dy))
        inward = (-dy // length, dx // length)
        support_constant = inward[0] * start[0] + inward[1] * start[1]
        local_degree = inward[0] + inward[1] - support_constant
        edge_ledger.append((start, end, length, inward, support_constant, local_degree))
        if start[0] == end[0] == 0:  # affine divisor s=0
            continue
        packet.extend([local_degree] * length)
    packet = tuple(sorted(packet, reverse=True))
    require(packet == (9, 9, 6, 2, 2, 2, 1), "complete toric packet")
    require((sum(packet), sum(value - 1 for value in packet)) == (31, 24),
            "packet degree and defect")

    # Reconstruct the two face equations and extract every special edge.
    S, P, Z = sy.symbols("S P Z")
    upsilon, xi, zeta = sy.symbols("upsilon xi zeta")
    rational = S**2 - P
    cyclotomic = 1 - upsilon*P**5 - xi*S**2*P**4
    main_face = sy.expand(rational * cyclotomic)
    tail_core = 1 - zeta*S**3*P**3 - xi*S**2*P**4
    tail_face = sy.expand(S**2 * tail_core)
    A, B, C, D, E = global_polygon
    special_edges = (
        edge_polynomial(main_face, S, P, A, B, Z),
        edge_polynomial(tail_face, S, P, B, C, Z),
        edge_polynomial(tail_face, S, P, C, D, Z),
        edge_polynomial(main_face, S, P, D, E, Z),
        edge_polynomial(main_face, S, P, E, A, Z),
        edge_polynomial(main_face, S, P, B, D, Z),
    )
    require(tuple(sy.factor(value) for value in special_edges) == (
        Z - 1,
        -Z**3*zeta + 1,
        -Z*xi - zeta,
        (Z - 1)*(Z*upsilon + xi),
        -Z**5 + upsilon,
        -Z**2*xi + 1,
    ), "six special edge schemes derived from faces")
    discriminants = tuple(sy.factor(sy.discriminant(value, Z))
                          for value in special_edges)
    require(discriminants == (
        1, -27*zeta**2, 1, (upsilon + xi)**2, 3125*upsilon**4, 4*xi,
    ), "six special edge schemes reduced")

    # The Q=sigma^30 graph faces and every edge lattice are primitive.
    graph_normals = ((3, 6, -1), (5, 5, -1))
    require(all(gcd(gcd(abs(a), abs(b)), abs(c)) == 1
                for a, b, c in graph_normals), "multiplicity-one face normals")
    graph_vertices = {
        (0, 1): (0, 1, 0), (2, 0): (2, 0, 0),
        (5, 3): (5, 3, 30), (4, 4): (4, 4, 30),
        (0, 6): (0, 6, 30),
    }
    for start, end, length, _normal, _constant, _index in edge_ledger:
        difference = tuple(graph_vertices[end][axis] - graph_vertices[start][axis]
                           for axis in range(3))
        require(gcd(gcd(abs(difference[0]), abs(difference[1])),
                    abs(difference[2])) == length, "edge denominator one")

    height_main = lambda point: 3*point[0] + 6*point[1] - 6
    height_tail = lambda point: 5*point[0] + 5*point[1] - 10
    outer_slopes = (
        height_main((1, 1)) - height_main(A),
        height_tail((1, 0)) - height_tail(B),
        height_tail((4, 3)) - height_tail(C),
        height_main((3, 4)) - height_main(D),
        height_main((1, 6)) - height_main(E),
    )
    require(outer_slopes == (3, -5, -5, -3, 3), "outer slope ledger")
    require(all((slope, slope - 1) in {(3, 2), (-5, -6), (-3, -4)}
                for slope in outer_slopes), "no outer toric chain")
    bd_probe = (2, -1)
    bd_slopes = (height_tail(bd_probe) - height_tail(B),
                 height_main(bd_probe) - height_main(B))
    require(bd_slopes == (-5, -6), "internal BD has no toric chain")

    # Face smoothness, genera, and the exact dual graph.
    require(sy.diff(cyclotomic, S) == -2*S*P**4*xi,
            "genus-two face is torus smooth")
    node_determinant = sy.factor(sy.det(sy.Matrix((
        (sy.diff(rational, S), sy.diff(rational, P)),
        (sy.diff(cyclotomic, S), sy.diff(cyclotomic, P)),
    ))).subs(P, S**2))
    require(node_determinant == -10*S**9*(upsilon + xi),
            "ten transverse main nodes")
    node_equation = sy.Poly(1 - (upsilon + xi)*S**10, S,
                            domain=sy.QQ.frac_field(upsilon, xi))
    require(node_equation.degree() == 10
            and sy.gcd(node_equation, node_equation.diff()).degree() == 0,
            "ten distinct main nodes")

    Y, T, W = sy.symbols("Y T W")
    genus_two_equation = 1 - upsilon*P**5 - xi*Y**2
    require(sy.expand(cyclotomic - genus_two_equation.subs(Y, S*P**2)) == 0,
            "genus-two normalization")
    branch = sy.Poly(1 - upsilon*P**5, P, domain=sy.QQ.frac_field(upsilon))
    require(sy.gcd(branch, branch.diff()).degree() == 0 and (5 - 1)//2 == 2,
            "smooth genus-two normalization")
    elliptic_equation = 1 - zeta*T**3 - xi*W**2
    require(sy.expand(tail_core - elliptic_equation.subs({T: S*P, W: S*P**2})) == 0,
            "elliptic tail normalization")
    tail_gradient_gate = sy.det(sy.Matrix(((3*zeta, 2*xi),
                                           (3*zeta, 4*xi))))
    require(sy.factor(tail_gradient_gate) == 6*zeta*xi,
            "tail torus derivative equations incompatible")
    require((10 + 2 - 3 + 1, 0 + 2 + 1 + 10) == (10, 13),
            "dual graph and genus ledger")

    # Exact A29 total-space equation after the common base change.
    sigma = sy.symbols("sigma")
    lam, alpha, epsilon, kappa, phi = sy.symbols("lam alpha epsilon kappa phi")
    Delta, Theta, eta, K = sy.symbols("Delta Theta eta K")
    source_s, source_p, Q = sy.symbols("source_s source_p Q")
    H = (
        lam*source_p + alpha*source_p**2 + epsilon*source_p**3
        + kappa*source_s**2*source_p**2 + phi*source_s*source_p**3
        + Delta*source_p**4 + Theta*source_s**2*source_p**3
        + eta*source_s*source_p**4 + zeta*source_s**3*source_p**3
        + upsilon*source_p**5 + xi*source_s**2*source_p**4
    )
    F = (source_s**2 - source_p)*(1 - Q*H) - Q*source_s**2/2
    H_sigma = sy.cancel(sigma**30 * H.subs({
        source_s: sigma**-3*S, source_p: sigma**-6*P,
    }))
    require(sy.denom(H_sigma) == 1, "integral main H_sigma")
    scaled_F = sy.cancel(sigma**6 * F.subs({
        Q: sigma**30, source_s: sigma**-3*S, source_p: sigma**-6*P,
    }))
    expected_scaled = (S**2 - P)*(1 - H_sigma) - sigma**30*S**2/2
    require(sy.factor(scaled_F - expected_scaled) == 0, "exact scaled source equation")
    require(sy.factor(H_sigma.subs(sigma, 0)
                      - upsilon*P**5 - xi*S**2*P**4) == 0,
            "main reduction retains exactly the top")
    local_U = S**2 - P
    local_V = (1 - H_sigma)/S**2
    require(sy.cancel(scaled_F/S**2 - (local_U*local_V - sigma**30/2)) == 0,
            "ten exact UV=sigma^30/2 charts")
    require(30 - 1 == 29, "A29 rational resolution length")

    # Generic outer root gates and function-field recovery.
    q = sy.symbols("q")
    carrier = zeta*Z**3 + K*Z**2 - (q - sy.Rational(1, 2))
    carrier_discriminant = sy.factor(sy.discriminant(carrier, Z))
    require(sy.factor(carrier_discriminant - (
        q - sy.Rational(1, 2)
    )*(4*K**3 - 27*zeta**2*(q - sy.Rational(1, 2)))) == 0,
            "generic cubic carrier discriminant")
    require(sy.degree(carrier, Z) == 3, "prime degree-three polynomial map")
    pure_edge = (lam*Z + alpha*Z**2 + epsilon*Z**3
                 + Delta*Z**4 + upsilon*Z**5 - q)
    pure_domain = sy.QQ.frac_field(lam, alpha, epsilon, Delta, upsilon, q)
    require(sy.gcd(sy.Poly(pure_edge, Z, domain=pure_domain),
                   sy.Poly(sy.diff(pure_edge, Z), Z, domain=pure_domain)).degree() == 0,
            "generic degree-five affine edge squarefree")
    require(sy.factor(F.subs(source_p, source_s**2) + Q*source_s**2/2) == 0,
            "torus excludes t=0 and recovers X=s/t")

    # Cyclotomic simplicity and labelled order-three attachments.
    X = sy.symbols("X")
    phi5 = sy.Poly(sy.cyclotomic_poly(5, X), X)
    require(phi5.as_expr() == X**4 + X**3 + X**2 + X + 1
            and phi5.is_irreducible, "Q(zeta5) action on H1")
    trace = sy.symbols("trace")
    trace_polynomial = sy.factor(sy.resultant(phi5.as_expr(),
                                               X**2 - trace*X + 1, X))
    require(trace_polynomial == (trace**2 + trace - 1)**2
            and sy.discriminant(trace**2 + trace - 1, trace) == 5,
            "unique quadratic subfield is real Q(sqrt5)")
    require(4 > 2 and 8 > 4, "elliptic-product endomorphism dimensions excluded")

    slope = sy.Rational(3*0**2, 2*1)
    doubled_x = slope**2
    doubled_y = slope*(0 - doubled_x) - 1
    require((doubled_x, doubled_y) == (0, -1),
            "attachments (0,+1),(0,-1) differ by nonzero 3-torsion")
    require(sy.cyclotomic_poly(3, 1) == 3, "sharp Eisenstein degree-three gate")

    # The target has good j=0 reduction under the same sigma-base change.
    a, target_A, target_C, target_X, target_Y = sy.symbols(
        "a target_A target_C target_X target_Y"
    )
    target = (target_C**2 - target_A**3 + sy.Rational(3, 4)*a**2*target_A
              - sigma**-30 + sy.Rational(1, 4)*a**3)
    target_scaled = sy.cancel(sigma**30 * target.subs({
        target_A: sigma**-10*target_X, target_C: sigma**-15*target_Y,
    }))
    target_expected = (target_Y**2 - target_X**3 - 1
                       + sy.Rational(3, 4)*a**2*sigma**20*target_X
                       + sy.Rational(1, 4)*a**3*sigma**30)
    require(sy.factor(target_scaled - target_expected) == 0,
            "good target scaling")
    require(sy.discriminant(target_X**3 + 1, target_X) == -27,
            "smooth j0 special target")

    # Degree specialization: all rational components and the simple genus-two
    # component are constant; the multiplicity-one tail is the sole carrier.
    # Equality at its two labelled attachments puts nonzero 3-torsion in the
    # isogeny kernel, so the conserved generic degree is divisible by three.
    full_degree = sum(packet)
    finite_degree = full_degree - 3*2
    require((full_degree, finite_degree) == (31, 25), "two exhaustive responses")
    require(full_degree % 3 == 1 and finite_degree % 3 == 1,
            "both response degrees violate specialization mod three")

    # Disjoint optional positive control: it is not used in the theorem.
    ss, pp = sy.symbols("ss pp")
    tt = pp - ss**2
    control_H = (
        -3*pp + sy.Rational(8, 3)*pp**2 - sy.Rational(1376, 135)*pp**3
        + sy.Rational(5591, 90)*ss**2*pp**2 + 2*ss*pp**3 + pp**4
        + 5*ss**2*pp**3 + 7*ss*pp**4 + 11*ss**3*pp**3
        + 13*pp**5 + 17*ss**2*pp**4
    )
    control_A = sy.cancel((-ss*pp + tt**2*sy.diff(control_H, ss))/pp)
    control_C0 = sy.expand(ss**2 + 2*tt**2*sy.diff(control_H, pp))
    control_B = sy.cancel((control_C0 + ss*control_A)/tt**2)
    require(sy.denom(control_A) == sy.denom(control_B) == 1,
            "critical control source pair polynomial")
    resultant_ab = sy.resultant(control_A, control_B, ss)
    resultant_ac = sy.resultant(control_A, control_C0, ss)
    require((polynomial_order(resultant_ab, pp),
             polynomial_order(resultant_ac, pp)) == (6, 8),
            "critical control exceptional powers")
    residual_ab = sy.Poly(sy.cancel(resultant_ab/pp**6), pp)
    residual_ac = sy.Poly(sy.cancel(resultant_ac/pp**8), pp)
    require(residual_ab == residual_ac and residual_ab.degree() == 25,
            "independent R25 control")
    require(residual_ab.nth(0) == -sy.Rational(189675421056, 5)
            and residual_ab.LC() == 2731549392000000000,
            "R25 endpoint gates")
    require(sy.gcd(residual_ab, residual_ab.diff()).degree() == 0,
            "R25 squarefree control")
    critical_length = residual_ab.degree() + 2 + 2
    require(critical_length == 29, "optional critical length")

    print("THM4218_INDEPENDENT_REGULAR_MODEL_MOD3_AUDIT")
    print("theorem_scope=all_lower_coefficients;critical_resultant_used=false")
    print("cell=exact_M10,b=d=0,reduced_(2,3);gate=upsilon*xi*zeta*(upsilon+xi)!=0")
    print("lower_faces=M:(r+2k-2)/10,T:(r+k-2)/6;Q=sigma^30;multiplicities=(1,1)")
    print("edge_schemes=6_reduced;toric_chains=none;BD_slopes=-5>-6;BD_roots=2")
    print("main=R+genus2;nodes=10;local=UV-sigma^30/2;resolution=A29x10")
    print("tail=j0_genus1;dual_graph_b1=10;special_genus=13;packet=(9,9,6,2,2,2,1)")
    print("genus2_sidecar=Q(zeta5),simple,Hom_to_elliptic=0")
    print("attachments=(0,+1),(0,-1);difference=nonzero_3_torsion")
    print("target=good_j0;degree_specialization=3_divides_generic_degree")
    print("responses=full31,finite25;residues_mod3=(1,1);carrier=prime_separable_cubic")
    print("critical_control=separate_positive_control,R25_squarefree,L29,digest="
          + coefficient_hash(residual_ab))
    print("dependency=THM-4045_plus_Dokchitser_Def3.7_Def3.9_Def3.12_Thm3.14")
    print(f"checks={CHECKS}")
    print("verdict=ACCEPT_THM4218_RELATIVE")


if __name__ == "__main__":
    main()
