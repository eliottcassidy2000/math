#!/usr/bin/env python3
"""Exact finite controls for the marked infinity-plumbing obstruction.

Universe: every assignment of seven vertex fibres in A_6 to the declared
tree presentation, with arrow a single 3-cycle. No transitivity filter.
The complete analytic argument, topology and scope live in the matching note.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
from itertools import permutations
import json

checks = 0


def check(condition, message):
    global checks
    checks += 1
    if not condition:
        raise RuntimeError(message)


ONE = tuple(range(6))


def mul(p, q):
    return tuple(p[q[i]] for i in range(6))


def inv(p):
    return tuple(p.index(i) for i in range(6))


def power(p, n):
    out = ONE
    for _ in range(n):
        out = mul(out, p)
    return out


def even(p):
    return sum(p[i] > p[j] for i in range(6) for j in range(i + 1, 6)) % 2 == 0


def cycle_type(p):
    seen, sizes = set(), []
    for i in range(6):
        if i in seen:
            continue
        j, size = i, 0
        while j not in seen:
            seen.add(j)
            size += 1
            j = p[j]
        sizes.append(size)
    return tuple(sorted(sizes, reverse=True))


def generated(gens):
    seen, todo = {ONE}, [ONE]
    while todo:
        p = todo.pop()
        for q in gens:
            z = mul(p, q)
            if z not in seen:
                seen.add(z)
                todo.append(z)
    return seen


def determinant(matrix):
    a = [[Fraction(x) for x in row] for row in matrix]
    ans = Fraction(1)
    for k in range(len(a)):
        p = next((i for i in range(k, len(a)) if a[i][k]), None)
        if p is None:
            return Fraction(0)
        if p != k:
            a[k], a[p] = a[p], a[k]
            ans = -ans
        pivot = a[k][k]
        ans *= pivot
        for i in range(k + 1, len(a)):
            r = a[i][k] / pivot
            for j in range(k + 1, len(a)):
                a[i][j] -= r * a[k][j]
    return ans


# Independently build the compact boundary from the six recorded centres.
weights, edges = {"L": 1}, set()


def blow(name, parents):
    for p in parents:
        weights[p] -= 1
    if len(parents) == 2:
        edges.remove(tuple(sorted(parents)))
    weights[name] = -1
    for p in parents:
        edges.add(tuple(sorted((name, p))))


for name, parents in [
    ("E1", ("L",)), ("E2", ("L", "E1")),
    ("E3", ("L", "E2")), ("E4", ("E3",)),
    ("E5", ("E4",)), ("E6", ("E4", "E5")),
]:
    blow(name, parents)
names = ["L", "E1", "E2", "E3", "E4", "E5", "E6"]
check([weights[n] for n in names] == [-2, -2, -2, -2, -3, -2, -1], "weights")
expected_edges = {tuple(sorted(e)) for e in [
    ("L", "E3"), ("E1", "E2"), ("E2", "E3"),
    ("E3", "E4"), ("E4", "E6"), ("E5", "E6"),
]}
check(edges == expected_edges, "tree edges")
matrix = [[-weights[a] if a == b else -int(tuple(sorted((a, b))) in edges)
           for b in names] for a in names]
check(determinant(matrix) == -1, "total boundary determinant")
meridian_degrees = [-6, -4, -8, -12, -10, -9, -18]
check([sum(matrix[i][j] * meridian_degrees[j] for j in range(7))
       for i in range(7)] == [0, 0, 0, 0, 0, 0, 1], "marked H1 sign control")
arm_data = []
for arm, expected in [(["L"], 2), (["E1", "E2"], 3),
                      (["E4", "E5", "E6"], 1),
                      (["L", "E1", "E2", "E3", "E4"], 9),
                      (["E5"], 2)]:
    ii = [names.index(n) for n in arm]
    value = determinant([[matrix[i][j] for j in ii] for i in ii])
    check(value == expected, "splice arm")
    arm_data.append([arm, int(value)])

A6 = [p for p in permutations(range(6)) if even(p)]
check(len(A6) == 360, "A6 universe")
squares, cubes = defaultdict(list), defaultdict(list)
for p in A6:
    squares[power(p, 2)].append(p)
    cubes[power(p, 3)].append(p)
common_types = {cycle_type(p) for p in set(squares) & set(cubes)}
check(common_types == {(1, 1, 1, 1, 1, 1), (2, 2, 1, 1), (5, 1)}, "square/cube types")


def relations(row):
    l, a, b, c, d, e, f, mu = row
    equations = [
        power(l, 2) == c, power(a, 2) == b, power(b, 2) == mul(a, c),
        power(c, 2) == mul(mul(l, b), d), power(d, 3) == mul(c, f),
        power(e, 2) == f, f == mul(mul(d, e), mu),
    ]
    equations += [mul(p, q) == mul(q, p) for p, q in
                  [(l, c), (a, b), (b, c), (c, d), (d, f), (e, f), (mu, f)]]
    return all(equations)


# Elimination is complete: c=l²=a³; b=a²; d=(lb)^(-1)c²;
# f=c^(-1)d³; e ranges over ALL square roots of f; mu=(de)^(-1)f.
rows, by_fixed, by_arrow = [], Counter(), Counter()
trivial_arrow = 0
for c in sorted(set(squares) & set(cubes)):
    for l in squares[c]:
        for a in cubes[c]:
            b = power(a, 2)
            d = mul(inv(mul(l, b)), power(c, 2))
            if mul(c, d) != mul(d, c):
                continue
            f = mul(inv(c), power(d, 3))
            for e in squares[f]:
                mu = mul(inv(mul(d, e)), f)
                row = (l, a, b, c, d, e, f, mu)
                if not relations(row):
                    continue
                if mu == ONE:
                    trivial_arrow += 1
                    check(all(p == ONE for p in row), "filled-arrow S3 control")
                if cycle_type(mu) != (3, 1, 1, 1):
                    continue
                fixed = [i for i in range(6) if all(p[i] == i for p in row)]
                check(c == f == ONE, "both central fibres vanish")
                check(len(fixed) >= 1, "no transitive marked image")
                check(cycle_type(d) == (3, 1, 1, 1), "shared three-cycle")
                support_left = {i for i in range(6) if l[i] != i or a[i] != i}
                support_right = {i for i in range(6) if d[i] != i or e[i] != i}
                check(len(support_left) <= 4 and len(support_right) <= 4, "triangle supports")
                check(len(support_left | support_right) <= 5, "joint support")
                by_fixed[len(fixed)] += 1
                by_arrow[mu] += 1
                rows.append(row)
check(trivial_arrow == 1, "one identity-meridian assignment")
check(len(rows) == 4000, "all labelled assignment count")
check(dict(by_fixed) == {3: 40, 2: 1800, 1: 2160}, "fixed-label distribution")
check(len(by_arrow) == 40 and set(by_arrow.values()) == {100}, "conjugacy completeness")

# A positive control prevents the false strengthening to an abelian image.
control = (
    (0, 1, 3, 2, 5, 4), (0, 1, 2, 4, 5, 3),
    (0, 1, 2, 5, 3, 4), ONE, (0, 1, 4, 2, 3, 5),
    (0, 2, 1, 4, 3, 5), ONE, (0, 2, 4, 3, 1, 5),
)
check(control in rows and relations(control), "named marked plumbing control")
control_group = generated((control[0], control[1], control[5]))
check(len(control_group) == 60, "A5 positive control")
check([i for i in range(6) if all(p[i] == i for p in control_group)] == [0], "exactly one fixed label")

# Independent finite corroboration of the triangle-action support lemma.
triangle_rows = 0
for x in squares[ONE]:
    for y in cubes[ONE]:
        xy = mul(x, y)
        if power(xy, 3) != ONE or cycle_type(xy) != (3, 1, 1, 1):
            continue
        triangle_rows += 1
        check(len({i for i in range(6) if x[i] != i or y[i] != i}) <= 4, "triangle action lemma")

trace = json.dumps(rows, separators=(",", ":"))
report = {
    "status": "FINITE-EXACT; analytic and actual-complement scope in matching note",
    "universe": "all A6 assignments to seven fibres, arrow a single 3-cycle",
    "graph_weights": weights,
    "graph_edges": sorted(edges),
    "vertex_multiples_of_arrow_in_H1": meridian_degrees,
    "arm_determinants": arm_data,
    "labelled_assignments": len(rows),
    "fixed_label_distribution": dict(sorted(by_fixed.items())),
    "arrows": len(by_arrow),
    "assignments_per_arrow": sorted(set(by_arrow.values())),
    "trivial_arrow_assignments": trivial_arrow,
    "triangle_controls": triangle_rows,
    "positive_control": control,
    "positive_control_group_order": len(control_group),
    "semantic_sha256": sha256(trace.encode()).hexdigest(),
    "checks": checks,
}
print(json.dumps(report, indent=2, sort_keys=True))
