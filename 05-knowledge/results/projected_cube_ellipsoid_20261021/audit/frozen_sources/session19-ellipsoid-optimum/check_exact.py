#!/usr/bin/env python3
"""Portable exact companion. Standard library only; no floating arithmetic.

Finite universes are declared in main(). Universal statements are proved in
README.md, not inferred from these tests. Every gate survives python -O.
"""
from collections import Counter, deque
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement, product
from math import prod
import json

COUNTS = Counter()


def require(test, label):
    COUNTS[label] += 1
    if not test:
        raise RuntimeError("FAILED: " + label)


def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def padd(a, b):
    return trim([(a[i] if i < len(a) else 0) +
                 (b[i] if i < len(b) else 0)
                 for i in range(max(len(a), len(b)))])


def pscale(a, s):
    return trim([s*x for x in a])


def pmul(a, b):
    c = [F(0)] * (len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i+j] += x*y
    return trim(c)


def pprod(ps):
    out = [F(1)]
    for p in ps:
        out = pmul(out, p)
    return out


def deriv(p):
    return trim([i*p[i] for i in range(1, len(p))] or [F(0)])


def peval(p, t):
    out = F(0)
    for x in reversed(p):
        out = out*t+x
    return out


def det(a):
    a = [list(row) for row in a]
    out = F(1)
    for i in range(len(a)):
        pivot = next((j for j in range(i, len(a)) if a[j][i]), None)
        if pivot is None:
            return F(0)
        if pivot != i:
            a[pivot], a[i] = a[i], a[pivot]
            out = -out
        x = a[i][i]
        out *= x
        for j in range(i+1, len(a)):
            ratio = a[j][i]/x
            for h in range(i+1, len(a)):
                a[j][h] -= ratio*a[i][h]
    return out


def num(v, q):
    return sum(v[i]**2*prod(q[j] for j in range(len(v)) if j != i)
               for i in range(len(v)))


def delta(v, q):
    return num(v, q)/sum(x*x for x in v)


def delta_basis(v, q):
    # Columns e_i - (v_i/v_last)e_last give a rational basis of H.
    m = len(v)-1
    z = [a/v[-1] for a in v[:-1]]
    gram = [[F(i == j)+z[i]*z[j] for j in range(m)] for i in range(m)]
    form = [[(q[i] if i == j else 0)+q[-1]*z[i]*z[j]
             for j in range(m)] for i in range(m)]
    return det(form)/det(gram)


def slacks(v, q):
    return {(i, j): (v[i]+v[j])**2-q[i]*v[j]**2-q[j]*v[i]**2
            for i, j in combinations(range(len(v)), 2)}


def feasible(v, q):
    return min(q) >= 0 and min(slacks(v, q).values()) >= 0


def star(v):
    return [(1+a/v[-1])**2 for a in v[:-1]]+[F(0)]


def dominant(v):
    a, b = sorted(v[:-1])[-2:]
    return v[-1]**2 >= v[-1]*(a+b)+a*b


def source(v):
    return [1+2*a/v[-1]-(a/v[-1])**2 for a in v]


def inspect_speed(v):
    k, M = len(v), v[-1]
    q = star(v)
    require(feasible(v, q) == dominant(v), "star_iff_threshold")
    require(delta(v, q) == delta_basis(v, q), "star_det_independent_basis")
    expected = M*M*prod(q[:-1])/sum(a*a for a in v)
    require(delta(v, q) == expected, "star_volume_formula")
    require(delta(v, q) > 0, "singular_lift_positive_projection")
    R = [[(M/(M+v[i]) if j == i else
           -v[i]/(M+v[i]) if j == k-1 else F(0))
          for j in range(k)] for i in range(k-1)]
    for i in range(k-1):
        require(sum(R[i][h]*v[h] for h in range(k)) == 0,
                "star_rows_annihilate_speed")
        for j in range(k-1):
            require(sum(R[i][h]*q[h]*R[j][h] for h in range(k)) == (i == j),
                    "star_covariance_is_identity")
    gram_R = [[sum(R[i][h]*R[j][h] for h in range(k))
               for j in range(k-1)] for i in range(k-1)]
    require(det(gram_R)*expected == 1, "inverse_star_jacobian")
    for i, j in combinations(range(k-1), 2):
        rhs = 2*v[i]*v[j]/M**2*(M**2-M*(v[i]+v[j])-v[i]*v[j])
        require(slacks(v, q)[i, j] == rhs, "exact_small_pair_slack")
    src = source(v)
    require(feasible(v, src), "source_inclusion")
    require(delta(v, src) == delta_basis(v, src), "source_det_independent_basis")
    if dominant(v):
        require(delta(v, q) > delta(v, src), "strict_source_volume_gain")
    A = q[:-1]
    b = [(a/M)**2 for a in v[:-1]]
    x = [bi/Ai for bi, Ai in zip(b, A)]
    caps = [[Ai, -bi] for Ai, bi in zip(A, b)] + [[F(0), F(1)]]
    poly_num = [F(0)]
    for i in range(k):
        poly_num = padd(poly_num, pscale(pprod(caps[j] for j in range(k) if i != j), v[i]**2))
    f = pprod([F(1), -xi] for xi in x)
    expected_poly = pscale(padd(f, pscale([F(0)]+deriv(f), -1)), M*M*prod(A))
    require(poly_num == expected_poly, "cap_polynomial_identity")
    rhs_derivative = pscale([F(0)]+deriv(deriv(f)), -M*M*prod(A))
    require(deriv(poly_num) == trim(rhs_derivative), "cap_derivative_identity")
    T = min(1/xi for xi in x)
    for t in (F(0), T/3, T/2, T):
        f2sum = 2*sum(x[i]*x[j]*prod(1-x[h]*t for h in range(k-1) if h not in (i, j))
                     for i, j in combinations(range(k-1), 2))
        require(peval(deriv(deriv(f)), t) == f2sum, "factored_second_derivative")
        require(f2sum >= 0 and (t == T or f2sum > 0), "second_derivative_sign")
        capped = [Ai-bi*t for Ai, bi in zip(A, b)]+[t]
        require(num(v, capped) == peval(poly_num, t), "cap_value_endpoint_safe")
        if t > 0:
            require(num(v, capped) < num(v, q), "strict_loss_including_endpoint")
        if dominant(v):
            require(feasible(v, capped), "dominant_caps_feasible")


def tree_from_prufer(seq):
    n = len(seq)+2
    degree = [1+seq.count(i) for i in range(n)]
    edges = []
    for x in seq:
        leaf = next(i for i in range(n) if degree[i] == 1)
        edges.append((leaf, x))
        degree[leaf] -= 1
        degree[x] -= 1
    leaves = [i for i in range(n) if degree[i] == 1]
    edges.append(tuple(leaves))
    return edges


def path(edges, a, b):
    adj = [[] for _ in range(len(edges)+1)]
    for i, j in edges:
        adj[i].append(j)
        adj[j].append(i)
    queue = deque([(a, [])])
    seen = {a}
    while queue:
        u, history = queue.popleft()
        if u == b:
            return history
        for w in adj[u]:
            if w not in seen:
                seen.add(w)
                queue.append((w, history+[(u, w)]))
    raise RuntimeError("Disconnected tree")


def tree_tests():
    trees = [tree_from_prufer(list(seq)) for seq in product(range(4), repeat=2)]
    require(len({tuple(sorted(tuple(sorted(e)) for e in T)) for T in trees}) == 16,
            "all_labelled_four_vertex_trees")
    successful, successful_nonstar = 0, 0
    for raw in product(range(1, 6), repeat=4):
        v = list(map(F, raw))
        r = [1/x for x in v]
        for T in trees:
            degrees = Counter(i for e in T for i in e)
            center = next((i for i, degree in degrees.items() if degree == 3), None)
            ok = all(sum((r[u]+r[w])**2 for u, w in path(T, i, j)) <= (r[i]+r[j])**2
                     for i, j in combinations(range(4), 2))
            if center is None:
                require(not ok, "nonstar_tree_impossible")
                successful_nonstar += int(ok)
            else:
                remaining = sorted(v[i] for i in range(4) if i != center)
                a, b = remaining[-2:]
                require(ok == (v[center]**2 >= v[center]*(a+b)+a*b),
                        "tree_star_iff_dominance")
            successful += int(ok)
    # Triangle plus isolated vertex: its weighted rows have rank two.
    R = [[F(1), -F(1), F(0), F(0)],
         [F(0), F(1), -F(1), F(0)],
         [F(1), F(0), -F(1), F(0)]]
    require(all(R[2][j] == R[0][j]+R[1][j] for j in range(4)), "triangle_rank_defect")
    return {"tree_speed_assignments": 625, "trees_per_assignment": 16,
            "successful_tree_certificates": successful,
            "successful_nonstar_certificates": successful_nonstar}


def boundary_direction_tests():
    # Every q_p=0 face control on sorted triples in 1..5 with other q in
    # {1/2,1,3/2}; check each slack-pivot direction with an explicit step.
    tested = 0
    for raw in combinations_with_replacement(range(1, 6), 3):
        v = list(map(F, raw))
        for p in range(3):
            other = [j for j in range(3) if j != p]
            for vals in product((F(1, 2), F(1), F(3, 2)), repeat=2):
                q = [F(0)]*3
                for j, val in zip(other, vals):
                    q[j] = val
                if not feasible(v, q):
                    continue
                for i in other:
                    slack = slacks(v, q)[tuple(sorted((p, i)))]
                    if slack <= 0:
                        continue
                    d = [v[j]**2/v[p]**2*(1 if j == i else -1) for j in range(3)]
                    d[p] = F(1)
                    steps = [slack/(4*v[i]**2)]
                    steps += [q[j]/(-2*d[j]) for j in range(3) if d[j] < 0]
                    epsilon = min(steps)
                    require(feasible(v, [qj+epsilon*dj for qj, dj in zip(q, d)]),
                            "boundary_direction_feasible")
                    poly = [F(0)]
                    for h in range(3):
                        poly = padd(poly, pscale(pprod([q[j], d[j]] for j in range(3) if j != h), v[h]**2))
                    derivative = poly[1] if len(poly) > 1 else F(0)
                    expected = 2*v[i]**2*prod(q[j] for j in range(3) if j not in (p, i))
                    require(derivative == expected and derivative > 0,
                            "boundary_direction_positive_derivative")
                    tested += 1
    return tested


def translated_cube_tests():
    grid = (F(-1), F(-1, 2), F(0), F(1, 2), F(1))
    feasible_count, equality_count = 0, 0
    for c in product(grid, repeat=2):
        for flat in product(grid, repeat=4):
            T = [list(flat[:2]), list(flat[2:])]
            rownorms = [sum(x*x for x in row) for row in T]
            if any(rownorms[i] > (1-abs(c[i]))**2 for i in range(2)):
                continue
            feasible_count += 1
            d2 = det(T)**2
            require(d2 <= 1, "translated_cube_volume_bound")
            if d2 == 1:
                require(c == (0, 0) and rownorms == [1, 1] and
                        sum(T[0][j]*T[1][j] for j in range(2)) == 0,
                        "translated_cube_unique_ellipse")
                equality_count += 1
    return {"translated_cube_feasible_controls": feasible_count,
            "translated_cube_equality_representatives": equality_count}


def explicit_controls():
    v = list(map(F, (1, 2, 4)))
    require(delta(v, star(v)) == F(75, 28), "example_optimal_determinant")
    require(delta(v, source(v)) == F(221, 84), "example_source_determinant")
    require(delta(v, star(v))-delta(v, source(v)) == F(1, 21), "example_difference")
    require(delta(v, star(v))/delta(v, source(v)) == F(225, 221), "example_squared_ratio")
    R = [[F(4, 5), F(0), -F(1, 5)], [F(0), F(2, 3), -F(1, 3)]]
    S = [[sum(R[i][h]*source(v)[h]*R[j][h] for h in range(3))
          for j in range(2)] for i in range(2)]
    require(S == [[F(1), F(2, 15)], [F(2, 15), F(1)]], "source_star_chart_covariance")
    require(S[0][0]+2*S[0][1]+S[1][1] > 2 and
            S[0][0]-2*S[0][1]+S[1][1] < 2, "volume_gain_not_containment")
    v = list(map(F, (1, 2, 3)))
    q = star(v)
    require((q[0]*v[1]**2+q[1]*v[0]**2)/(v[0]+v[1])**2 == F(89, 81),
            "outside_cone_hostile")
    v = list(map(F, (2, 3, 6)))
    require(dominant(v) and slacks(v, star(v))[0, 1] == 0, "equality_threshold_contact")
    # Two zero lift coordinates have zero projected determinant at k>=3.
    require(delta(list(map(F, (1, 2, 4))), list(map(F, (0, 0, 1)))) == 0,
            "two_zero_coordinates_degenerate")
    for a, M in combinations_with_replacement(range(1, 8), 2):
        v = list(map(F, (a, M)))
        tmax = F((a+M)**2, a*a)
        for fraction in (F(0), F(1, 4), F(1, 2), F(1)):
            t = fraction*tmax
            q = [F((a+M)**2, M*M)-F(a*a, M*M)*t, t]
            require(feasible(v, q) and delta(v, q) == F((a+M)**2, a*a+M*M),
                    "k2_nonunique_lift_constant")
    require(F(1, 1)/F(9, 5) == F(5, 9) and F(5, 9) != F(3, 5),
            "compressed_inverse_not_ambient_inverse")
    for M in range(2, 31):
        v = list(map(F, (1, 2, M)))
        require(dominant(v) == (M >= 4), "integer_family_threshold")
        require(delta(v, star(v)) == F((M+1)**2*(M+2)**2, M*M*(M*M+5)),
                "integer_family_determinant")


def main():
    speed_counts = {}
    for k in range(3, 7):
        count = 0
        for raw in combinations_with_replacement(range(1, 7), k):
            inspect_speed(list(map(F, raw)))
            count += 1
        speed_counts[str(k)] = count
    feasible_grid = 0
    for raw in combinations_with_replacement(range(1, 7), 3):
        v = list(map(F, raw))
        for q in product((F(0), F(1, 2), F(1), F(3, 2), F(2)), repeat=3):
            require(delta(v, q) == delta_basis(v, q), "grid_det_independent_basis")
            if feasible(v, q):
                feasible_grid += 1
                if dominant(v):
                    require(delta(v, q) <= delta(v, star(v)), "grid_global_bound")
                    if delta(v, q) == delta(v, star(v)):
                        require(list(q) == star(v), "grid_unique_diagonal_lift")
    explicit_controls()
    tree_report = tree_tests()
    cube_report = translated_cube_tests()
    boundary_count = boundary_direction_tests()
    report = {"status": "PASS", "arithmetic": "fractions.Fraction; no floating point",
              "speed_universe": "nondecreasing positive integer tuples in 1..6, k=3..6",
              "speed_tuple_counts": speed_counts,
              "grid_universe": "all sorted triples in 1..6; q in {0,1/2,1,3/2,2}^3",
              "feasible_grid_controls": feasible_grid,
              "boundary_improvement_controls": boundary_count,
              "tree_universe": "all 16 labelled trees on 4 vertices; v in {1,..,5}^4",
              "source_example": {"v": [1, 2, 4], "optimal_delta": "75/28",
                                 "source_delta": "221/84", "squared_volume_ratio": "225/221"},
              **tree_report, **cube_report,
              "checks_by_label": dict(sorted(COUNTS.items())),
              "total_exact_gates": sum(COUNTS.values())}
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
