#!/usr/bin/env python3
"""Exact certificates for AP anchors, common neighbors, and square scaling.

Standard library only. All checks survive python -O. The finite parameter
universe below checks identities and examples, not the infinite theorem.
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd, isqrt, lcm
import json

CHECKS = 0


def check(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise ArithmeticError(message)


def square_root(q):
    q = F(q)
    if q < 0:
        return None
    n, d = isqrt(q.numerator), isqrt(q.denominator)
    return F(n, d) if n*n == q.numerator and d*d == q.denominator else None


def add(P, Q, d):
    if P is None:
        return Q
    if Q is None:
        return P
    x, y = P
    X, Y = Q
    if x == X and y == -Y:
        return None
    s = (3*x*x-d*d)/(2*y) if P == Q else (Y-y)/(X-x)
    z = s*s-x-X
    R = z, s*(x-z)-y
    check(R[1]**2 == R[0]**3-d*d*R[0], "group-law curve membership")
    return R


def area_point(A, B, C, d):
    A, B, C, d = map(F, (A, B, C, d))
    check(A*A+B*B == C*C and A*B/2 == d, "input right triangle")
    X = d*(A+C)/B
    Y = 2*d*X/B
    check(Y*Y == X**3-d*d*X and X > d, "input curve point")
    return X, Y


def neighbor(P, a, d):
    X, Y = P
    check(Y != 0 and X > d, "positive component")
    A = abs((X*X-d*d)/Y)
    B = abs(2*d*X/Y)
    C = abs((X*X+d*d)/Y)
    check(A*A+B*B == C*C and A*B/2 == d, "area identity")
    roots = abs(A-B)/2, C/2, (A+B)/2
    x = roots[0]**2-a
    check(roots[1]**2-roots[0]**2 == d, "first AP difference")
    check(roots[2]**2-roots[1]**2 == d, "second AP difference")
    check(x == add(P, P, d)[0]-(a+d), "duplication coordinate")
    for j, r in enumerate(roots):
        check(x+a+j*d == r*r, "common-neighbor equation")
    return x, roots


def integer_packet(a, d, chosen):
    L = lcm(*(r.denominator for _, roots in chosen for r in roots))
    core = tuple(L*L*(a+j*d) for j in range(3))
    leaves = tuple(L*L*x for x, _ in chosen)
    labels = core+leaves
    check(all(F(x).denominator == 1 and x > 0 for x in labels), "positive integers")
    labels = tuple(int(x) for x in labels)
    check(len(set(labels)) == len(labels), "all vertices distinct")
    check(labels[1]-labels[0] == labels[2]-labels[1], "AP core survives")
    roots_out = []
    edges = list(combinations(range(3), 2))
    edges += [(i, j) for i in range(3) for j in range(3, len(labels))]
    for i, j in edges:
        root = isqrt(labels[i]+labels[j])
        check(root*root == labels[i]+labels[j], "integer edge certificate")
        roots_out.append({"edge": [i, j], "root": str(root)})
    return L, labels, roots_out


def fixed_integer_neighbors(a, d):
    solutions = []
    pairs = 0
    for u in range(1, isqrt(d)+1):
        if d % u:
            continue
        v = d//u
        if (u+v) % 2:
            continue
        pairs += 1
        p, q = (v-u)//2, (v+u)//2
        check(q*q-p*p == d, "factor-pair difference")
        r = isqrt(q*q+d)
        if r*r == q*q+d:
            solutions.append((p*p-a, p, q, r))
    return solutions, pairs


def main():
    # Complete primitive Euclidean-parameter universe m <= 30, both leg
    # orientations collapsed; exactly the strict positive AP cores survive.
    admissible = 0
    root_samples = 0
    selected_neighbor_counts = {}
    for m in range(2, 31):
        for n in range(1, m):
            if gcd(m, n) != 1 or (m-n) % 2 == 0:
                continue
            A, B, C = m*m-n*n, 2*m*n, m*m+n*n
            a, d = 2*(C*C-4*A*B), 8*A*B
            if a <= 0:
                continue
            admissible += 1
            core = (a, a+d, a+2*d)
            for i, j in combinations(range(3), 2):
                check(square_root(core[i]+core[j]) is not None, "PPT core edge")
            P = area_point(4*A, 4*B, 4*C, d)
            Q = None
            chosen = []
            for multiple in range(1, 9):
                Q = add(Q, P, d)
                x, roots = neighbor(Q, a, d)
                root_samples += 1
                if multiple == 1:
                    check(x == a+d, "initial equation-solution is middle core vertex")
                if len(chosen) < 3 and x > 0 and x not in core and x not in [z for z, _ in chosen]:
                    chosen.append((x, roots))
            selected_neighbor_counts[str(len(chosen))] = selected_neighbor_counts.get(str(len(chosen)), 0)+1
            integer_packet(a, d, chosen)

    # Displayed certificate: E5 point, scaled by (24^2,24^3) to E2880.
    P = (F(25, 4), F(75, 8))
    Q = None
    records = []
    chosen = []
    for multiple in range(1, 6):
        Q = add(Q, P, 5)
        R = 576*Q[0], 13824*Q[1]
        x, roots = neighbor(R, 482, 2880)
        if multiple == 1:
            check(x == 3362, "excluded core collision")
        if multiple == 2:
            check(x == F(-1856688659521, 3878673841), "negative double control")
        if multiple >= 3:
            chosen.append((x, roots))
            records.append({
                "multiple": multiple,
                "point_on_E5": list(map(str, Q)),
                "neighbor": str(x),
                "roots": list(map(str, roots)),
                "root_denominator": str(lcm(*(r.denominator for r in roots))),
            })
    L, labels, edge_roots = integer_packet(482, 2880, chosen)
    check(L == 13851986720466207228947843753473356093372173219780057654419750125010709636370,
          "displayed common scale")
    for x, y in combinations(labels[3:], 2):
        check(isqrt(x+y)**2 != x+y, "certificate leaf pairs are nonedges")

    # Complete divisor enumeration for the original fixed core. Independent
    # exhaustive p-scan uses (q-p)(q+p)=2880 and hence p <= 719.
    fixed, pairs = fixed_integer_neighbors(482, 2880)
    brute = []
    for p in range(720):
        q = isqrt(p*p+2880)
        r = isqrt(p*p+5760)
        if q*q == p*p+2880 and r*r == p*p+5760:
            brute.append((p*p-482, p, q, r))
    check(fixed == brute == [(3362, 62, 82, 98)], "complete unscaled hostile")
    check(not [s for s in fixed if s[0] > 0 and s[0] not in (482, 3362, 6242)],
          "unscaled core has zero graph common neighbors")

    # Explicit three-anchor success without a clique: gap 6 from (3,4,5).
    # (1,7,13) is not a triangle; x=1/4-1 is negative initially, but later
    # multiples supply a positive common neighbor. Gap controls existence.
    P6 = area_point(3, 4, 5, 6)
    Q6 = add(P6, P6, 6)
    x6, roots6 = neighbor(Q6, 1, 6)
    check(x6 > 0 and x6 not in (1, 7, 13), "nonclique triple has common neighbor")
    check(square_root(1+7) is None, "nonclique positive control")

    # Four-anchor hostile: these three consecutive AP square equations work,
    # but adding the fourth anchor fails; the universal impossibility is the
    # cited Fermat theorem, not this bounded example.
    check(square_root(F(1, 4)) == F(1, 2), "three-square start")
    check(square_root(F(25, 4)) == F(5, 2), "three-square middle")
    check(square_root(F(49, 4)) == F(7, 2), "three-square end")
    check(square_root(F(73, 4)) is None, "fourth square fails")

    print(json.dumps({
        "status": "FINITE-EXACT; universal claims have separate proofs",
        "checks": CHECKS,
        "universe": {"euclidean_m_max": 30, "positive_AP_cores": admissible,
                     "elliptic_neighbor_samples": root_samples,
                     "multipliers_checked_per_core": 8,
                     "fresh_positive_neighbors_selected_up_to_three": selected_neighbor_counts},
        "certificate": {"core_before_scaling": [482, 3362, 6242],
                        "common_scale_L": str(L), "label_multiplier_L_squared": str(L*L),
                        "rational_neighbors": records, "integer_labels": list(map(str, labels)),
                        "required_edges": edge_roots, "induced_graph_on_these_six_vertices": "K3 join I3"},
        "fixed_core_hostile": {"factor_pairs_tested": pairs,
                               "all_equation_solutions": fixed, "external_integer_neighbors": []},
        "nonclique_positive_control": {"anchors": [1, 7, 13], "neighbor": str(x6),
                                       "roots": list(map(str, roots6))},
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
