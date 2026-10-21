#!/usr/bin/env python3
"""Exact standalone single-violating-pair certificate replay."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement
from math import prod
import json

COUNTS = Counter()


def need(test, label):
    COUNTS[label] += 1
    if not test:
        raise RuntimeError("FAILED: "+label)


def eye(n):
    return [[F(i == j) for j in range(n)] for i in range(n)]


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def quadratic(S, n):
    return sum(n[i]*S[i][j]*n[j] for i in range(len(n)) for j in range(len(n)))


def delta(v, q):
    return sum(v[i]**2*prod(q[j] for j in range(len(v)) if j != i)
               for i in range(len(v)))/sum(a*a for a in v)


def chart_covariance(v, q):
    M = v[-1]
    return [[(M*M*q[i]+v[i]**2*q[-1])/(M+v[i])**2 if i == j else
             v[i]*v[j]*q[-1]/((M+v[i])*(M+v[j]))
             for j in range(len(v)-1)] for i in range(len(v)-1)]


def inspect(v):
    m, M = len(v)-1, v[-1]
    pairs = {}
    for i, j in combinations(range(m), 2):
        a, b = v[i], v[j]
        A, B = b*(M+a)/(M*(a+b)), a*(M+b)/(M*(a+b))
        need(0 < A <= 1 and 0 < B <= 1, "positive_coefficients_at_most_one")
        pairs[i, j] = A, B
    bad = [ij for ij, (A, B) in pairs.items() if A*A+B*B > 1]
    if m >= 3:
        a, b, d = v[-2], v[-3], v[-4]
        stated = M*M < M*(a+b)+a*b and M*M >= M*(a+d)+a*d
        need((len(bad) == 1) == stated, "ordered_single_pair_criterion")
    star_delta = M*M*prod((1+a/M)**2 for a in v[:-1])/sum(a*a for a in v)
    if not bad:
        q = [(1+a/M)**2 for a in v[:-1]]+[F(0)]
        need(chart_covariance(v, q) == eye(m), "star_lift_identity")
        if m == 2:
            a, b = v[:2]
            need(star_delta == (M+a)**2*(M+b)**2/(M*M*(M*M+a*a+b*b)),
                 "triple_star_piece")
            D = -M*M+M*(a+b)+a*b
            if D == 0:
                q_cont = [(a+b)/b+a*(b-a)/(b*M), (a+b)/a+b*(a-b)/(a*M), F(0)]
                need(q_cont == q, "triple_boundary_continuity")
        return 0
    if len(bad) != 1:
        return len(bad)
    i, j = bad[0]
    a, b = v[i], v[j]
    A, B = pairs[i, j]
    c = (A*A+B*B-1)/(2*A*B)
    D = -M*M+M*(a+b)+a*b
    need(c == D/((M+a)*(M+b)) and 0 < c < 1, "exact_correlation_positive_definite")
    S = eye(m)
    S[i][j] = S[j][i] = c
    lp = c/(A*B*(1-c*c))
    weights = [F(1)]*m
    weights[i] = (1-c*A/B)/(1-c*c)
    weights[j] = (1-c*B/A)/(1-c*c)
    need(min(weights+[lp]) > 0, "strictly_positive_contact_weights")
    n = [F(0)]*m
    n[i], n[j] = A, -B
    need(quadratic(S, n) == 1, "repaired_pair_contact")
    cert = [[(weights[u] if u == w else F(0))+lp*n[u]*n[w]
             for w in range(m)] for u in range(m)]
    for u in range(m):
        for w in range(m):
            need(sum(cert[u][h]*S[h][w] for h in range(m)) == (u == w),
                 "exact_inverse_certificate")
    need(sum(weights)+lp == m, "trace_weight_normalization")
    for (u, w), (C, E) in pairs.items():
        nn = [F(0)]*m
        nn[u], nn[w] = C, -E
        need(quadratic(S, nn) <= 1, "all_leaf_facets_contain_repair")
    repaired = star_delta*(1-c*c)
    simplified = 4*(M*(a+b)+a*b)/sum(x*x for x in v)*prod(
        (1+v[h]/M)**2 for h in range(m) if h not in (i, j))
    need(repaired == simplified, "simplified_repaired_determinant")
    need(((M+a)*(M+b))**2-D*D == 4*M*M*(M*(a+b)+a*b),
         "difference_of_squares_identity")
    if m == 2:
        q = [(a+b)/b+a*(b-a)/(b*M), (a+b)/a+b*(a-b)/(a*M), D/(a*b)]
        need(min(q) > 0, "triple_lift_strict_positivity")
        need(chart_covariance(v, q) == S, "triple_lift_represents_repair")
        need(delta(v, q) == repaired, "triple_piece_determinant")
        for u, w in combinations(range(3), 2):
            need(q[u]*v[w]**2+q[w]*v[u]**2 <= (v[u]+v[w])**2,
                 "triple_lift_feasibility")
    else:
        other = next(h for h in range(m) if h not in (i, j))
        need(S[i][j] > 0 and S[i][other] == 0 and
             v[i]*v[other]/((M+v[i])*(M+v[other])) > 0,
             "higher_dimension_diagonal_obstruction")
    return 1


def main():
    regime_counts = {}
    for k in range(3, 7):
        regimes = Counter()
        for raw in combinations_with_replacement(range(1, 9), k):
            nbad = inspect(list(map(F, raw)))
            regimes["star" if nbad == 0 else "one_pair" if nbad == 1 else "multiple_pairs"] += 1
        regime_counts[str(k)] = dict(sorted(regimes.items()))
    v = list(map(F, (1, 2, 3)))
    need(delta(v, [F(5, 3), F(7, 3), F(1)]) == F(22, 7), "triple_example_optimum")
    src = [1+2*a/v[-1]-(a/v[-1])**2 for a in v]
    need(delta(v, src) == F(64, 21), "triple_example_source")
    need(F(22, 7)/F(64, 21) == F(33, 32), "triple_example_squared_volume_ratio")
    v = list(map(F, (1, 2, 3, 5)))
    need(4*(5*(2+3)+2*3)/sum(a*a for a in v)*(1+v[0]/v[-1])**2 == F(1488, 325),
         "quadruple_example_optimum")
    need(F(-25+5*(2+3)+2*3, (5+2)*(5+3)) == F(3, 28),
         "quadruple_example_correlation")
    for k in range(4, 11):
        family = list(map(F, [1]*(k-3)+[2, 2, 4]))
        need(inspect(family) == 1, "arbitrary_dimension_family_single_pair")
    need(1-F(1, 9)**2 == F(80, 81), "arbitrary_dimension_family_chart_determinant")
    print(json.dumps({"status": "PASS", "arithmetic": "fractions.Fraction; no floating point",
                      "universe": "nondecreasing positive integer tuples in 1..8, k=3..6",
                      "additional_controls": "family (1,...,1,2,2,4) in k=4..10",
                      "regime_counts": regime_counts,
                      "scope": "complete k=3 classification; only one-violating-pair repair for k>=4",
                      "checks_by_label": dict(sorted(COUNTS.items())),
                      "total_exact_gates": sum(COUNTS.values())}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
