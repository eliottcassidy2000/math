#!/usr/bin/env python3
"""Independent exact audit for the THM-4147 packet.

This program was written without importing the primary implementation.  It
reconstructs the valued support from the displayed source equation itself,
computes the projected Newton polygon and all toric residue indices, uses a
new rational coefficient control, recomputes the affine critical eliminants,
and attacks the commutator/overlap inequality through S_6.

It is deliberately an audit of the algebraic ledgers, not a substitute for
the geometric carrier and fixed-sheet arguments discussed in the report.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations
from math import gcd

import sympy as sy


def check(test: bool, label: str) -> None:
    if not test:
        raise ArithmeticError(label)


def q_valued_support(poly, s, p, qinv):
    """Return (s exponent,p exponent,q valuation,leading coefficient).

    Unlike the primary certificate's endpoint-construction routine, this
    starts with the fully expanded source polynomial and combines every
    coincident (s,p) coefficient before taking its q-adic valuation.
    """
    projected = sy.Poly(sy.expand(poly), s, p)
    answer = []
    for (es, ep), coefficient in projected.terms():
        coefficient = sy.Poly(coefficient, qinv)
        nonzero = [(power[0], value) for power, value in coefficient.terms()
                   if value != 0]
        check(nonzero, "zero projected coefficient survived Poly")
        valuation = min(power for power, _ in nonzero)
        leading = coefficient.nth(valuation)
        answer.append((es, ep, valuation, sy.factor(leading)))
    return tuple(sorted(answer))


def lower_planes(points):
    """Enumerate rational lower supporting planes z=ax+by+c."""
    xyz = [(row[0], row[1], row[2]) for row in points]
    planes = set()
    for a0, a1, a2 in combinations(xyz, 3):
        matrix = sy.Matrix([
            [a0[0], a0[1], 1],
            [a1[0], a1[1], 1],
            [a2[0], a2[1], 1],
        ])
        if matrix.det() == 0:
            continue
        aa, bb, cc = matrix.inv() * sy.Matrix([a0[2], a1[2], a2[2]])
        gaps = [sy.factor(z - aa*x - bb*y - cc) for x, y, z in xyz]
        if all(gap >= 0 for gap in gaps):
            incident = tuple(sorted(point for point, gap in zip(xyz, gaps)
                                    if gap == 0))
            planes.add((sy.Rational(aa), sy.Rational(bb), sy.Rational(cc),
                        incident))
    # discard planes supported on only a collinear segment after duplicates
    return tuple(sorted(planes, key=str))


def cross(o, a, b):
    return (a[0]-o[0])*(b[1]-o[1]) - (a[1]-o[1])*(b[0]-o[0])


def polygon(points):
    pts = sorted(set((row[0], row[1]) for row in points))
    low = []
    for pt in pts:
        while len(low) >= 2 and cross(low[-2], low[-1], pt) <= 0:
            low.pop()
        low.append(pt)
    high = []
    for pt in reversed(pts):
        while len(high) >= 2 and cross(high[-2], high[-1], pt) <= 0:
            high.pop()
        high.append(pt)
    return tuple(low[:-1] + high[:-1])


def polygon_and_edges(vertices):
    area2_signed = sum(
        vertices[i][0]*vertices[(i+1) % len(vertices)][1]
        - vertices[(i+1) % len(vertices)][0]*vertices[i][1]
        for i in range(len(vertices))
    )
    check(area2_signed > 0, "polygon orientation is not counterclockwise")
    boundary = 0
    rows = []
    for start, stop in zip(vertices, vertices[1:] + vertices[:1]):
        dx, dy = stop[0]-start[0], stop[1]-start[1]
        length = gcd(abs(dx), abs(dy))
        boundary += length
        inward = (-dy//length, dx//length)
        constant = inward[0]*start[0] + inward[1]*start[1]
        residue_index = inward[0] + inward[1] - constant
        rows.append((start, stop, length, inward, constant, residue_index))
    check((area2_signed-boundary+2) % 2 == 0, "Pick parity")
    genus = (area2_signed-boundary+2)//2
    return area2_signed, boundary, genus, tuple(rows)


def edge_polynomial(poly, vertices, s, p, w):
    """Coefficient polynomial along each oriented projected boundary edge."""
    pp = sy.Poly(sy.expand(poly), s, p)
    answer = []
    for start, stop in zip(vertices, vertices[1:] + vertices[:1]):
        dx, dy = stop[0]-start[0], stop[1]-start[1]
        length = gcd(abs(dx), abs(dy))
        step = (dx//length, dy//length)
        coeffs = []
        for k in range(length+1):
            exp = (start[0]+k*step[0], start[1]+k*step[1])
            coeffs.append(pp.coeff_monomial(s**exp[0]*p**exp[1]))
        initial = sy.factor(sum(coeffs[k]*w**k for k in range(length+1)))
        answer.append((start, stop, sy.factor(initial)))
    return tuple(answer)


def t_valuation(poly, t):
    terms = sy.Poly(poly, t).terms()
    check(terms, "valuation of zero")
    return min(exp[0] for exp, coefficient in terms if coefficient != 0)


def monic_digest(poly):
    p = sy.Poly(poly).monic()
    payload = ";".join(str(sy.factor(c)) for c in p.all_coeffs())
    return sha256(payload.encode("utf-8")).hexdigest()


def compose(a, b):
    return tuple(a[b[i]] for i in range(len(a)))


def invert(a):
    result = [0]*len(a)
    for i, value in enumerate(a):
        result[value] = i
    return tuple(result)


def perm_index(a):
    seen = set()
    cycles = 0
    for i in range(len(a)):
        if i in seen:
            continue
        cycles += 1
        j = i
        while j not in seen:
            seen.add(j)
            j = a[j]
    return len(a)-cycles


def commutator_hostile(limit=6):
    tested = 0
    equality = 0
    for n in range(1, limit+1):
        group = tuple(permutations(range(n)))
        inverses = {a: invert(a) for a in group}
        supports = {a: frozenset(i for i, value in enumerate(a) if i != value)
                    for a in group}
        for x in group:
            xi = inverses[x]
            for y in group:
                yi = inverses[y]
                comm = compose(compose(compose(x, y), xi), yi)
                overlap = len(supports[x] & supports[y])
                lhs = perm_index(comm)
                check(lhs <= 2*overlap,
                      f"commutator overlap counterexample in S_{n}: {x},{y}")
                equality += int(lhs == 2*overlap and overlap > 0)
                tested += 1
    return tested, equality


def source_polynomial(s, p, qinv, K, phi, delta, theta, eta, zeta):
    y = s*p
    H = (-3*p + sy.Rational(8, 3)*p**2 - sy.Rational(1376, 135)*p**3
         + K*y**2 + phi*p**2*y + delta*p**4 + theta*p*y**2
         + eta*p**3*y + zeta*y**3)
    return sy.expand((s**2-p)*(1-qinv*H) - qinv*s**2/2)


def critical_control(name, delta, phi, theta, eta, zeta, expected_degree):
    X, T = sy.symbols("X T")
    K = sy.Rational(2848, 45)-sy.Rational(7, 6)*delta
    P = T+X**2*T**2
    Y = X*T*P
    G = (-X**2*T/2 - 3*P + sy.Rational(8, 3)*P**2
         - sy.Rational(1376, 135)*P**3 + K*Y**2 + phi*P**2*Y
         + delta*P**4 + theta*P*Y**2 + eta*P**3*Y + zeta*Y**3)
    gx_over_t = sy.cancel(sy.diff(G, X)/T)
    gt = sy.diff(G, T)
    check(sy.denom(gx_over_t) == 1, f"{name}: G_X/T is not polynomial")

    # Degree-at-infinity audit: for T != 0, the two X leading rows cannot
    # simultaneously drop inside any admitted chamber.
    gxpoly, gtpoly = sy.Poly(gx_over_t, X), sy.Poly(gt, X)
    top = sy.factor(eta+zeta)
    check(gxpoly.degree() == 8 and gtpoly.degree() == 9,
          f"{name}: unexpected generic X degrees")
    check(sy.factor(gxpoly.LC() - 9*top*T**8) == 0,
          f"{name}: G_X/T leading row")
    check(sy.factor(gtpoly.LC() - 9*top*T**8) == 0,
          f"{name}: G_T leading row")
    check(top != 0, f"{name}: planted top collision")

    resultant = sy.factor(sy.resultant(gx_over_t, gt, X))
    valuation = t_valuation(resultant, T)
    check(valuation == 56, f"{name}: resultant T valuation {valuation}")
    quotient = sy.cancel(resultant/(T**56*(6*T+1)**2))
    check(sy.denom(quotient) == 1, f"{name}: nonpolynomial strict quotient")
    residual = sy.Poly(quotient, T)
    check(residual.degree() == expected_degree,
          f"{name}: residual degree {residual.degree()}")
    check(residual.nth(0) != 0 and residual.LC() != 0,
          f"{name}: endpoint failure")
    check(sy.gcd(residual, residual.diff()).degree() == 0,
          f"{name}: residual is not squarefree")

    # Restore the two coordinate strata directly, and audit that the
    # universal -1/6 factor has precisely the advertised common roots.
    check(sy.factor(gx_over_t.subs(T, 0)+X) == 0,
          f"{name}: T=0 first specialization")
    check(sy.factor(gt.subs(T, 0)+(X**2+6)/2) == 0,
          f"{name}: T=0 second specialization")
    special_f = sy.Poly(gx_over_t.subs(T, -sy.Rational(1, 6)), X)
    special_h = sy.Poly(gt.subs(T, -sy.Rational(1, 6)), X)
    common = sy.gcd(special_f, special_h).monic()
    check(sy.expand(common.as_expr()-(X**2-6)) == 0,
          f"{name}: universal -1/6 gcd is {common.as_expr()}")
    remainder = sy.rem(sy.Poly(G.subs(T, -sy.Rational(1, 6))-sy.Rational(1,2), X),
                       common)
    check(remainder.is_zero, f"{name}: universal pair has wrong value")

    # A second chart is obtained directly from the rational function in
    # (s,p), not by substituting into the first resultant.  Saturating the
    # identically vanishing factor p in G_s gives an unrelated eliminant.
    ss, pp = sy.symbols("ss pp")
    tt = pp-ss**2
    yy = ss*pp
    HH = (-3*pp + sy.Rational(8,3)*pp**2
          - sy.Rational(1376,135)*pp**3 + K*yy**2
          + phi*pp**2*yy + delta*pp**4 + theta*pp*yy**2
          + eta*pp**3*yy + zeta*yy**3)
    Gsp = -ss**2/(2*tt)+HH
    raw_s = sy.cancel(tt**2*sy.diff(Gsp,ss))
    raw_p = sy.cancel(2*tt**2*sy.diff(Gsp,pp))
    check(sy.denom(raw_s) == sy.denom(raw_p) == 1,
          f"{name}: second-chart clearing")
    sat_s = sy.cancel(raw_s/pp)
    check(sy.denom(sat_s) == 1,
          f"{name}: G_s numerator lacks the p saturation factor")
    second_resultant = sy.factor(sy.resultant(sat_s,raw_p,ss))
    pval = t_valuation(second_resultant,pp)
    expected_pval = 4 if name == "P" else 8
    check(pval == expected_pval,
          f"{name}: second-chart p valuation {pval}")
    second_residual = sy.Poly(sy.cancel(second_resultant/pp**pval),pp)
    check(second_residual.degree() == expected_degree,
          f"{name}: second-chart residual degree")
    check(second_residual.nth(0) != 0 and second_residual.LC() != 0
          and sy.gcd(second_residual,second_residual.diff()).degree() == 0,
          f"{name}: second-chart residual endpoint/squarefree gate")

    # The p=0 line was removed by saturation.  Directly on the rational
    # chart G_s vanishes there and G_p restores s^2=1/6, value 1/2.
    check(sy.factor(sy.diff(Gsp,ss).subs(pp,0)) == 0,
          f"{name}: collapsed p=0 G_s")
    gp_at_zero = sy.factor(sy.diff(Gsp,pp).subs(pp,0))
    check(sy.factor(gp_at_zero - (1-6*ss**2)/(2*ss**2)) == 0,
          f"{name}: collapsed p=0 G_p is {gp_at_zero}")
    check(sy.factor(Gsp.subs({pp:0,ss**2:sy.Rational(1,6)})
                    - sy.Rational(1,2)) == 0,
          f"{name}: collapsed pair value")

    return {
        "name": name,
        "K": K,
        "valuation": valuation,
        "degree": residual.degree(),
        "constant": residual.nth(0),
        "leading": residual.LC(),
        "digest": monic_digest(residual),
        "second_pval": pval,
        "second_digest": monic_digest(second_residual),
        "critical_length": residual.degree()+4,
    }


def main():
    s, p, Q, w = sy.symbols("s p Q w")

    carrier_K, carrier_zeta, carrier_q = sy.symbols(
        "carrier_K carrier_zeta carrier_q"
    )
    carrier_polynomial = (carrier_zeta*w**3 + carrier_K*w**2
                          - (carrier_q-sy.Rational(1, 2)))
    carrier_expected = ((carrier_q-sy.Rational(1, 2))
                        * (4*carrier_K**3
                           - 27*carrier_zeta**2
                           * (carrier_q-sy.Rational(1, 2))))
    check(sy.factor(sy.discriminant(carrier_polynomial,w)
                    - carrier_expected) == 0,
          "cubic carrier Hurwitz discriminant")

    # This control is disjoint from the primary script's rational point.
    # Delta=-64/105 makes the forced K exactly 64.
    delta = -sy.Rational(64, 105)
    K = sy.Rational(64)
    phi = sy.Rational(1)
    theta = sy.Rational(2)
    check(K == sy.Rational(2848,45)-sy.Rational(7,6)*delta,
          "forced K/Delta row")

    cases = [
        ("P", sy.Rational(1), sy.Rational(0),
         ((0,1),(2,0),(4,2),(4,3),(3,4),(1,5),(0,5)),
         (29,11,10), (8,5,4,3,2,2,1), 20, 21, 2, 25),
        ("Y", sy.Rational(0), sy.Rational(2),
         ((0,1),(2,0),(5,3),(3,4),(0,5)),
         (30,10,11), (11,8,2,2,2,1), 20, 20, 3, 26),
        ("B", sy.Rational(1), sy.Rational(2),
         ((0,1),(2,0),(5,3),(1,5),(0,5)),
         (31,11,11), (8,8,4,2,2,2,1), 21, 21, 3, 27),
    ]

    geometry_rows = []
    critical_rows = []
    for (name, eta, zeta, expected_vertices, expected_pick, expected_packet,
         expected_residual, finite_n, beta, full_n) in cases:
        F = source_polynomial(s,p,Q,K,phi,delta,theta,eta,zeta)
        support = q_valued_support(F,s,p,Q)
        vertices = polygon(support)
        check(vertices == expected_vertices, f"{name}: polygon {vertices}")
        area2, boundary, genus, edges = polygon_and_edges(vertices)
        check((area2,boundary,genus) == expected_pick,
              f"{name}: Pick ledger {(area2,boundary,genus)}")
        packet = []
        for start, stop, length, inward, constant, index in edges:
            if start[0] == stop[0] == 0:
                check(length == 4 and index == 1,
                      f"{name}: vertical affine edge")
                continue
            packet.extend([index]*length)
        packet = tuple(sorted(packet,reverse=True))
        check(packet == expected_packet, f"{name}: packet {packet}")
        defect = sum(index-1 for index in packet)
        check(defect == 2*genus-2, f"{name}: RH saturation")
        planes = lower_planes(support)
        initials = edge_polynomial(F,vertices,s,p,w)
        geometry_rows.append((name, support, vertices, planes, edges, initials,
                              (area2,boundary,genus),packet,defect))

        control = critical_control(name,delta,phi,theta,eta,zeta,
                                   expected_residual)
        L = control["critical_length"]
        check(L > finite_n+beta, f"{name}: finite index gate")
        check(defect > 2*(full_n-L), f"{name}: full overlap gate")
        control.update(finite_n=finite_n,beta=beta,full_n=full_n,
                       defect=defect,finite_margin=L-finite_n-beta,
                       full_margin=defect-2*(full_n-L))
        critical_rows.append(control)

    # Hostile/cancellation walls.  These use exact coefficients different
    # from the main B control and check which cancellations are geometric.
    epsilon = -sy.Rational(1376,135)
    delta_keqeps = sy.factor(sy.Rational(6,7)
                             *(sy.Rational(2848,45)-epsilon))
    support_keqeps = q_valued_support(
        source_polynomial(s,p,Q,epsilon,phi,delta_keqeps,theta,1,0),s,p,Q)
    check(polygon(support_keqeps) == cases[0][3],
          "K=epsilon changed the P polygon")
    support_thetaeqdelta = q_valued_support(
        source_polynomial(s,p,Q,K,phi,delta,delta,1,0),s,p,Q)
    check(polygon(support_thetaeqdelta) == cases[0][3],
          "Theta=Delta changed the P polygon")
    support_etaeqzeta = q_valued_support(
        source_polynomial(s,p,Q,K,phi,delta,theta,1,1),s,p,Q)
    check(polygon(support_etaeqzeta) == cases[2][3],
          "eta=zeta changed the B polygon")
    midpoint_control = critical_control("B_midpoint",delta,phi,theta,1,1,21)

    Xc,Tc = sy.symbols("Xc Tc")
    Pc = Tc+Xc**2*Tc**2
    Yc = Xc*Tc*Pc
    Gc = (-Xc**2*Tc/2 - 3*Pc + sy.Rational(8,3)*Pc**2
          - sy.Rational(1376,135)*Pc**3 + K*Yc**2 + phi*Pc**2*Yc
          + delta*Pc**4 + theta*Pc*Yc**2 + Pc**3*Yc-Yc**3)
    collision = sy.factor(sy.resultant(sy.cancel(sy.diff(Gc,Xc)/Tc),
                                       sy.diff(Gc,Tc),Xc))
    collision_val = t_valuation(collision,Tc)
    collision_strict = sy.cancel(collision/(Tc**collision_val*(6*Tc+1)**2))
    check(collision_val == 42 and sy.denom(collision_strict) == 1
          and sy.degree(collision_strict,Tc) == 19,
          "eta+zeta=0 collision ledger")

    hostile_pairs, equality_pairs = commutator_hostile(6)

    print("THM-4147 INDEPENDENT EXACT AUDIT")
    print("control: Delta=-64/105, K=64, Phi=1, Theta=2, eta/zeta by chamber")
    for row in geometry_rows:
        name,support,vertices,planes,edges,initials,pick,packet,defect=row
        print(f"case={name}")
        print(f"  valued_support={support}")
        print(f"  polygon={vertices}; Pick={pick}")
        print(f"  lower_planes={tuple((a,b,c) for a,b,c,_ in planes)}")
        print(f"  edges={edges}")
        print(f"  edge_initials={initials}")
        print(f"  packet={packet}; defect={defect}")
    for row in critical_rows:
        print(f"critical={row['name']}: T^{row['valuation']}*(6T+1)^2*Q{row['degree']}; "
              f"second=p^{row['second_pval']}*R{row['degree']}; "
              f"L={row['critical_length']}; residual_sha256={row['digest']}; "
              f"second_sha256={row['second_digest']}")
        print(f"  constant={row['constant']}; leading={row['leading']}")
        print(f"  finite=(n={row['finite_n']},beta={row['beta']},"
              f"margin={row['finite_margin']}); full=(n={row['full_n']},"
              f"defect={row['defect']},margin={row['full_margin']})")
    print(f"commutator_overlap: pairs={hostile_pairs}; sharp_nonzero_pairs={equality_pairs}; PASS")
    print("cubic_carrier_discriminant=(q-1/2)*(4K^3-27zeta^2*(q-1/2)): PASS")
    print("support_hostiles: K=epsilon,Theta=Delta,eta=zeta preserve polygons: PASS")
    print(f"midpoint_eta_equals_zeta: T^{midpoint_control['valuation']}*(6T+1)^2*Q{midpoint_control['degree']}: PASS")
    print(f"top_collision_eta_plus_zeta_zero: T^{collision_val}*(6T+1)^2*Q{sy.degree(collision_strict,Tc)}: CHANGES_STRICT_TRANSFORM")
    print("verdict=ALGEBRAIC_LEDGERS_PASS")


if __name__ == "__main__":
    main()
