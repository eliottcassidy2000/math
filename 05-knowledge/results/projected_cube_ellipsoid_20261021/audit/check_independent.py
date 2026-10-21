#!/usr/bin/env python3
"""Independent exact audit; no producer imports and no third-party packages.

Uses ambient projection matrices, exact linear-system ranks, and a direct
facet-basis covariance construction rather than the producer implementations.
Universal conclusions require the separately retained written proof audit.
"""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement, product
from math import prod
import json

C = Counter()


def check(value, label):
    C[label] += 1
    if not value:
        raise RuntimeError(label)


def identity(n):
    return [[F(i == j) for j in range(n)] for i in range(n)]


def transpose(a):
    return [list(x) for x in zip(*a)]


def mul(a, b):
    return [[sum(x*y for x, y in zip(row, col)) for col in zip(*b)] for row in a]


def plus(a, b):
    return [[x+y for x, y in zip(row, other)] for row, other in zip(a, b)]


def scale(a, s):
    return [[s*x for x in row] for row in a]


def outer(a, b):
    return [[x*y for y in b] for x in a]


def diagonal(q):
    return [[x if i == j else F(0) for j in range(len(q))] for i, x in enumerate(q)]


def determinant(a):
    a = [list(map(F, row)) for row in a]
    result = F(1)
    for i in range(len(a)):
        j = next((j for j in range(i, len(a)) if a[j][i]), None)
        if j is None:
            return F(0)
        if j != i:
            a[j], a[i] = a[i], a[j]
            result = -result
        pivot = a[i][i]
        result *= pivot
        a[i] = [x/pivot for x in a[i]]
        for j in range(i+1, len(a)):
            t = a[j][i]
            a[j] = [x-t*y for x, y in zip(a[j], a[i])]
    return result


def inverse(a):
    n = len(a)
    rows = [list(map(F, row))+unit for row, unit in zip(a, identity(n))]
    for i in range(n):
        j = next(j for j in range(i, n) if rows[j][i])
        rows[j], rows[i] = rows[i], rows[j]
        p = rows[i][i]
        rows[i] = [x/p for x in rows[i]]
        for j in range(n):
            if j != i:
                t = rows[j][i]
                rows[j] = [x-t*y for x, y in zip(rows[j], rows[i])]
    return [row[n:] for row in rows]


def rank(a):
    a = [list(map(F, row)) for row in a]
    r = 0
    for c in range(len(a[0])):
        pivot = next((i for i in range(r, len(a)) if a[i][c]), None)
        if pivot is None:
            continue
        a[pivot], a[r] = a[r], a[pivot]
        p = a[r][c]
        a[r] = [x/p for x in a[r]]
        for i in range(r+1, len(a)):
            t = a[i][c]
            a[i] = [x-t*y for x, y in zip(a[i], a[r])]
        r += 1
        if r == len(a):
            break
    return r


def quad(a, w):
    return sum(w[i]*a[i][j]*w[j] for i in range(len(w)) for j in range(len(w)))


def trace(a):
    return sum(a[i][i] for i in range(len(a)))


def geometry(v):
    n = len(v)
    N = scale(outer(v, v), 1/sum(x*x for x in v))
    P = plus(identity(n), scale(N, -1))
    normals = {}
    for i, j in combinations(range(n), 2):
        w = [F(0)]*n
        w[i], w[j] = v[j]/(v[i]+v[j]), -v[i]/(v[i]+v[j])
        normals[i, j] = w
    L = [normals[i, n-1] for i in range(n-1)]
    J = mul(transpose(L), inverse(mul(L, transpose(L))))
    return N, P, normals, L, J


def covariance(P, q):
    return mul(mul(P, diagonal(q)), P)


def feasible(v, q):
    return min(q) >= 0 and all(q[i]*v[j]**2+q[j]*v[i]**2 <= (v[i]+v[j])**2
                              for i, j in combinations(range(len(v)), 2))


def run_shape(v):
    n, m, M = len(v), len(v)-1, v[-1]
    N, P, normals, L, J = geometry(v)
    I = identity(m)
    check(mul(L, J) == I and mul(J, L) == P, "exact_chart_inverse_and_projection")
    q = [(1+a/M)**2 for a in v[:-1]]+[F(0)]
    star = covariance(P, q)
    check(star == mul(J, transpose(J)), "singular_star_ambient_covariance")
    delta = determinant(plus(star, N))
    expected = M*M*prod(q[:-1])/sum(x*x for x in v)
    check(delta == expected and delta > 0, "ambient_rank_one_determinant_normalization")
    bad = [ij for ij, w in normals.items() if quad(star, w) > 1]
    condition = M*M >= M*(v[-2]+v[-3])+v[-2]*v[-3]
    check((not bad) == condition, "exact_threshold_from_ambient_supports")
    src = [1+2*a/M-(a/M)**2 for a in v]
    source_cov = covariance(P, src)
    check(all(quad(source_cov, w) <= 1 for w in normals.values()), "source_containment_ambient")
    if not bad:
        check(delta > determinant(plus(source_cov, N)), "dominant_strict_source_volume")
    if len(bad) == 1:
        i, j = bad[0]
        w = normals[i, j]
        # Derive chart normal by pullback rather than use the stated A,B formulas.
        z = mul([w], J)[0]
        c = (1-sum(t*t for t in z))/(2*z[i]*z[j])
        S = identity(m)
        S[i][j] = S[j][i] = c
        repaired = mul(mul(J, S), transpose(J))
        check(0 < c < 1, "repair_positive_definite")
        check(all(quad(repaired, x) <= 1 for x in normals.values()), "repair_all_ambient_facets")
        Sinv = inverse(S)
        lam = Sinv[i][j]/(z[i]*z[j])
        axis = [Sinv[h][h]-lam*z[h]**2 for h in range(m)]
        check(min(axis+[lam]) > 0, "derived_positive_contact_weights")
        check(plus(diagonal(axis), scale(outer(z, z), lam)) == Sinv,
              "derived_inverse_contact_decomposition")
        check(sum(axis)+lam == m and quad(S, z) == 1, "contact_trace_and_pair_tangency")
        measured = determinant(plus(repaired, N))
        simple = 4*(M*(v[i]+v[j])+v[i]*v[j])/sum(x*x for x in v)*prod(
            (1+v[h]/M)**2 for h in range(m) if h not in (i, j))
        check(measured == simple and measured == delta*determinant(S),
              "repaired_euclidean_determinant")
        equations = [[P[r][h]*P[s][h] for h in range(n)]
                     for r in range(n) for s in range(r, n)]
        rhs = [repaired[r][s] for r in range(n) for s in range(r, n)]
        augmented = [row+[b] for row, b in zip(equations, rhs)]
        check(rank(equations) == n, "diagonal_restriction_injectivity")
        if n >= 4:
            check(rank(augmented) == n+1, "no_real_diagonal_lift_linear_system")
        else:
            a, b = v[:2]
            t = (-M*M+M*(a+b)+a*b)/(a*b)
            lift = [(1+a/M)**2-a*a*t/M**2, (1+b/M)**2-b*b*t/M**2, t]
            check(min(lift) > 0 and covariance(P, lift) == repaired,
                  "triple_positive_diagonal_representation")
            check(feasible(v, lift) and rank(augmented) == n, "triple_feasibility_and_consistency")
    return "star" if not bad else "one_pair" if len(bad) == 1 else "multiple_pairs"


def boundary_checks():
    for n in (3, 4):
        for raw in combinations_with_replacement((F(1, 2), F(1), F(2)), n):
            v = list(raw)
            N, P, _, _, _ = geometry(v)
            for p in range(n):
                for vals in product((F(1, 2), F(1)), repeat=n-1):
                    q, it = [], iter(vals)
                    for j in range(n):
                        q.append(F(0) if j == p else next(it))
                    if not feasible(v, q):
                        continue
                    A = plus(covariance(P, q), N)
                    for i in range(n):
                        if i == p:
                            continue
                        slack = (v[p]+v[i])**2-q[i]*v[p]**2
                        if slack <= 0:
                            continue
                        d = [v[j]**2/v[p]**2*(1 if j == i else -1) for j in range(n)]
                        d[p] = F(1)
                        epsilon = min([slack/(4*v[i]**2)]+
                                      [q[j]/(-2*d[j]) for j in range(n) if d[j] < 0])
                        check(feasible(v, [x+epsilon*y for x, y in zip(q, d)]),
                              "boundary_direction_exact_feasibility")
                        derivative = determinant(A)*trace(mul(inverse(A), covariance(P, d)))
                        expected = 2*v[i]**2*prod(q[j] for j in range(n) if j not in (p, i))/sum(x*x for x in v)
                        check(derivative == expected and derivative > 0,
                              "boundary_direction_matrix_derivative")


def tree_checks():
    edges = list(combinations(range(4), 2))
    for raw in product((F(1, 2), F(1), F(3)), repeat=4):
        v = list(raw)
        _, _, normals, _, _ = geometry(v)
        bases = 0
        for choice in combinations(edges, 3):
            R = [normals[e] for e in choice]
            if rank(R) != 3:
                continue
            bases += 1
            J = mul(transpose(R), inverse(mul(R, transpose(R))))
            cov = mul(J, transpose(J))
            valid = all(quad(cov, w) <= 1 for w in normals.values())
            deg = Counter(i for e in choice for i in e)
            centers = [i for i in range(4) if deg[i] == 3]
            if not centers:
                check(not valid, "nonstar_facet_basis_direct_matrix_failure")
            else:
                p = centers[0]
                a, b = sorted(v[i] for i in range(4) if i != p)[-2:]
                check(valid == (v[p]**2 >= v[p]*(a+b)+a*b), "star_facet_basis_exact_criterion")
        check(bases == 16, "labelled_four_vertex_basis_count")


def controls():
    for raw, desired in [((1, 2, 4), F(75, 28)), ((1, 2, 3), F(22, 7)),
                         ((1, 2, 3, 5), F(1488, 325)), ((2, 3, 6), F(144, 49))]:
        v = list(map(F, raw))
        N, P, normals, L, J = geometry(v)
        S = identity(len(v)-1)
        for e, w in normals.items():
            z = mul([w], J)[0]
            if sum(x*x for x in z) > 1:
                i, j = e
                S[i][j] = S[j][i] = (1-sum(x*x for x in z))/(2*z[i]*z[j])
        cov = mul(mul(J, S), transpose(J))
        check(determinant(plus(cov, N)) == desired, "explicit_optimum_determinants")
    v = list(map(F, (1, 2, 4)))
    N, P, _, L, _ = geometry(v)
    qsrc = [F(23, 16), F(7, 4), F(2)]
    scov = covariance(P, qsrc)
    check(determinant(plus(scov, N)) == F(221, 84), "source_example_determinant")
    chart = mul(mul(L, scov), transpose(L))
    check(quad(chart, [1, 1]) > 2 and quad(chart, [1, -1]) < 2, "volume_not_containment_hostile")
    for raw in combinations_with_replacement((F(1, 2), F(1), F(3)), 2):
        v = list(raw)
        N, P, _, _, _ = geometry(v)
        a, M = v
        for u in (F(0), F(1, 5), F(1, 2), F(1)):
            q = [(1-u)*(a+M)**2/M**2, u*(a+M)**2/a**2]
            check(determinant(plus(covariance(P, q), N)) == (a+M)**2/(a*a+M*M),
                  "k2_continuum_same_projected_interval")
    # Matrix implementation of the inverse/compression hostile, without cancellation.
    v = list(map(F, (1, 2)))
    N, P, _, _, _ = geometry(v)
    forward = covariance(P, list(map(F, (2, 1))))
    inv_restricted = plus(inverse(plus(forward, N)), scale(N, -1))
    wrong = covariance(P, [F(1, 2), F(1)])
    check(trace(inv_restricted) == F(5, 9) and trace(wrong) == F(3, 5),
          "compressed_inverse_hostile_matrix")
    for k in range(4, 11):
        check(run_shape(list(map(F, [1]*(k-3)+[2, 2, 4]))) == "one_pair",
              "dimension_independent_nonempty_regime")


def translated_checks():
    # Covariances, not only linear image matrices: the exact triangle optimum
    # itself is present even when no rational square-root matrix exists.
    for c, normal in [(F(0), None), (F(1, 10), [F(8, 9), -F(5, 9)])]:
        S = [[F(1), c], [c, F(1)]]
        normals = [[F(1), F(0)], [F(0), F(1)]] + ([] if normal is None else [normal])
        for a, b, d in product((F(0), F(1, 4), F(1, 2), F(3, 4), F(1)),
                               (F(0), F(1, 4), F(1, 2), F(3, 4), F(1)),
                               (F(-1), F(-1, 2), F(0), F(1, 10), F(1, 2), F(1))):
            T = [[a, d], [d, b]]
            if determinant(T) < 0:
                continue
            for h in product((F(-1, 2), F(0), F(1, 2)), repeat=2):
                if not all(abs(sum(x*y for x, y in zip(w, h))) <= 1 and
                           quad(T, w) <= (1-abs(sum(x*y for x, y in zip(w, h))))**2
                           for w in normals):
                    continue
                check(trace(mul(inverse(S), T)) <= 2 and determinant(T) <= determinant(S),
                      "translated_covariance_trace_and_volume_bound")
                if determinant(T) == determinant(S):
                    check(T == S and h == (0, 0), "translated_covariance_unique_equality")


def main():
    regimes = {}
    pool = (F(1, 2), F(2, 3), F(1), F(3, 2), F(2), F(3))
    for k in range(3, 6):
        count = Counter(run_shape(list(v)) for v in combinations_with_replacement(pool, k))
        regimes[str(k)] = dict(sorted(count.items()))
    boundary_checks()
    tree_checks()
    controls()
    translated_checks()
    print(json.dumps({"status": "PASS", "arithmetic": "Python standard-library Fraction only",
                      "independence": "No producer code imported; ambient projection and exact linear-system checks",
                      "universe": "nondecreasing tuples from {1/2,2/3,1,3/2,2,3}, k=3..5",
                      "regime_counts": regimes,
                      "additional_universes": ["all 16 labelled four-vertex facet bases; speeds in {1/2,1,3}^4",
                                               "boundary controls k=3,4; speeds in {1/2,1,2}; positive q in {1/2,1}",
                                               "nonempty one-pair family k=4..10 and exact named examples"],
                      "checks_by_label": dict(sorted(C.items())),
                      "total_exact_gates": sum(C.values()),
                      "scope": "Finite companion to the retained analytic audit, not a substitute for real-parameter proofs"},
                     indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
