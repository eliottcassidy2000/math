#!/usr/bin/env python3
"""Realizable unit-ear boundary tomography: exact arithmetic, no packages.

All squared-distance exclusions use rational enclosures of explicit real
radicals. Prescribed unit edges have independent algebraic identities. stdout
is a deterministic receipt; --write additionally saves coordinates and tables.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations, product
from math import factorial, gcd, isqrt
from pathlib import Path
import hashlib
import json
import sys


def need(ok, detail):
    if not ok:
        raise RuntimeError(detail)


@dataclass(frozen=True)
class Q3:
    a: F = F(0)
    b: F = F(0)

    def __add__(self, y):
        if not isinstance(y, Q3):
            y = Q3(F(y))
        return Q3(self.a + y.a, self.b + y.b)

    __radd__ = __add__

    def __neg__(self):
        return Q3(-self.a, -self.b)

    def __sub__(self, y):
        return self + (-y if isinstance(y, Q3) else -F(y))

    def __mul__(self, y):
        if not isinstance(y, Q3):
            y = Q3(F(y))
        return Q3(self.a * y.a + 3 * self.b * y.b,
                  self.a * y.b + self.b * y.a)

    __rmul__ = __mul__

    def inv(self):
        d = self.a * self.a - 3 * self.b * self.b
        need(d != 0, ("zero denominator", self))
        return Q3(self.a / d, -self.b / d)

    def __truediv__(self, y):
        if not isinstance(y, Q3):
            y = Q3(F(y))
        return self * y.inv()

    def serial(self):
        return [str(self.a), str(self.b)]


BITS = 160
SCALE = 1 << BITS


def plus(A, B):
    return A[0] + B[0], A[1] + B[1]


def times(A, B):
    values = [x * y for x in A for y in B]
    return min(values), max(values)


def minus(A, B):
    return A[0] - B[1], A[1] - B[0]


def square(A):
    lo = 0 if A[0] <= 0 <= A[1] else min(x * x for x in A)
    return F(lo), max(x * x for x in A)


def root(A):
    need(A[0] >= 0, ("negative root interval", A))
    lo = isqrt(A[0].numerator * SCALE * SCALE // A[0].denominator)
    hi = isqrt(A[1].numerator * SCALE * SCALE // A[1].denominator) + 1
    return F(lo, SCALE), F(hi, SCALE)


ROOT3 = root((F(3), F(3)))


def interval(z):
    return plus((z.a, z.a), times((z.b, z.b), ROOT3))


def qnorm(v):
    return v[0] * v[0] + v[1] * v[1]


@dataclass(frozen=True)
class Point:
    # (base_x,base_y)+(coeff_x,coeff_y)*sqrt(rad), all entries in Q(sqrt3).
    base: tuple
    coeff: tuple = (Q3(), Q3())
    rad: Q3 = Q3()

    def boxes(self):
        R = root(interval(self.rad)) if self.rad != Q3() else (F(0), F(0))
        return tuple(plus(interval(a), times(interval(b), R))
                     for a, b in zip(self.base, self.coeff))

    def serial(self):
        return dict(base=[a.serial() for a in self.base],
                    coefficient=[a.serial() for a in self.coeff],
                    radicand=self.rad.serial())


def distance_box(A, B):
    return plus(*(square(minus(a, b)) for a, b in zip(A.boxes(), B.boxes())))


def rational_parameters():
    yield F(0)
    for h in range(1, 65):
        for denominator in range(1, h + 1):
            for numerator in range(-h, h + 1):
                if max(abs(numerator), denominator) == h and gcd(numerator, denominator) == 1:
                    yield F(numerator, denominator)


def ear(points, u, v):
    """Find a certified induced 3-edge path between two Q(sqrt3) old points."""
    U, V = points[u], points[v]
    need(U.coeff == V.coeff == (Q3(), Q3()), "old terminals must be Q(sqrt3)")
    attempts = 0
    for t in rational_parameters():
        denominator = 1 + t * t
        shift = (Q3((1 - t * t) / denominator), Q3(2 * t / denominator))
        need(qnorm(shift) == Q3(F(1)), ("unit parameter", t))
        X = tuple(a + b for a, b in zip(U.base, shift))
        delta = tuple(a - b for a, b in zip(V.base, X))
        d2 = qnorm(delta)
        dbox = interval(d2)
        if not (dbox[0] > 0 and dbox[1] < 4):
            continue
        D = d2.inv() - F(1, 4)
        need(interval(D)[0] > 0, "intersection radicand")
        midpoint = tuple((a + b) / 2 for a, b in zip(X, V.base))
        perp = (-delta[1], delta[0])
        # (delta/2) is perpendicular to perp and their squared norms are d2/4,d2.
        need(delta[0] * perp[0] + delta[1] * perp[1] == Q3(), "perpendicular")
        need(d2 * (D + F(1, 4)) == Q3(F(1)), "two unit intersection edges")
        for sign in (1, -1):
            attempts += 1
            x = Point(X)
            y = Point(midpoint, tuple(sign * z for z in perp), D)
            candidates = (x, y)
            good = True
            for j, z in enumerate(candidates):
                for i, old in enumerate(points):
                    db = distance_box(z, old)
                    if db[0] <= 0:
                        good = False
                        break
                    if i != (u if j == 0 else v) and db[0] <= 1 <= db[1]:
                        good = False
                        break
                if not good:
                    break
            if good:
                return candidates, dict(terminals=[u, v], t=str(t), sign=sign,
                                        attempts=attempts)
    raise RuntimeError(("bounded rational search found no ear", u, v, len(points)))


def unit_graph(points, required):
    edges = set()
    separated = 0
    for i, j in combinations(range(len(points)), 2):
        db = distance_box(points[i], points[j])
        need(db[0] > 0, ("collision not excluded", i, j))
        if (i, j) in required:
            need(db[0] <= 1 <= db[1], ("unit edge interval", i, j, db))
            edges.add((i, j))
        else:
            need(not db[0] <= 1 <= db[1], ("unexpected or unresolved unit contact", i, j))
            separated += 1
    return edges, separated


def partition(words):
    names = {}
    return tuple(names.setdefault(x, len(names)) for x in words)


def partitions(n):
    if n == 0:
        return [()]
    result = [(0,)]
    for _ in range(1, n):
        result = [a + (j,) for a in result for j in range(max(a) + 2)]
    return result


def forest(pi):
    roots = {}
    edges = []
    for i, b in enumerate(pi):
        if b in roots:
            edges.append((roots[b], i))
        else:
            roots[b] = i
    return tuple(edges)


def component_partition(n, edges):
    parent = list(range(n))
    def find(i):
        while parent[i] != i:
            i = parent[i]
        return i
    for i, j in edges:
        parent[find(j)] = find(i)
    return partition([find(i) for i in range(n)])


def coarser(beta, pi):
    return all(pi[i] != pi[j] or beta[i] == beta[j]
               for i, j in combinations(range(len(pi)), 2))


def falling(q, k):
    return factorial(q) // factorial(q - k) if k <= q else 0


def recover(P, queries, q):
    K = q * q - 3 * q + 3
    moments = {}
    for pi in P:
        E = forest(pi)
        M = 0
        for mask in range(1 << len(E)):
            selected = tuple(e for j, e in enumerate(E) if mask >> j & 1)
            sub = component_partition(len(pi), selected)
            need(forest(sub) == selected, ("canonical subforest", pi, selected, sub))
            M += (-1) ** len(selected) * K ** (len(E) - len(selected)) * queries[sub]
        moments[pi] = M
    masses = {}
    for pi in sorted(P, key=lambda a: max(a)):
        masses[pi] = moments[pi] - sum(mass for beta, mass in masses.items() if coarser(beta, pi))
    return moments, masses


def proper_count(n, edges, q):
    # Independent labelled backtracking for small hostile/control graphs.
    neighbors = [set() for _ in range(n)]
    for u, v in edges:
        neighbors[u].add(v)
        neighbors[v].add(u)
    order = sorted(range(n), key=lambda v: -len(neighbors[v]))
    colours = [-1] * n
    def visit(i):
        if i == n:
            return 1
        u = order[i]
        excluded = {colours[v] for v in neighbors[u] if colours[v] >= 0}
        result = 0
        for colour in range(q):
            if colour not in excluded:
                colours[u] = colour
                result += visit(i + 1)
        colours[u] = -1
        return result
    return visit(0)


def main():
    source = Path(__file__).read_bytes()
    root_path = Path(__file__).resolve().parents[4]
    source_pins = json.loads(Path("SOURCE_PINS.json").read_text(encoding="utf-8"))
    inherited_bytes = {}
    for pin in source_pins:
        body = (root_path / pin["path"]).read_bytes()
        need(hashlib.sha256(body).hexdigest() == pin["sha256"], ("inherited source pin", pin["path"]))
        inherited_bytes[pin["path"]] = body
    half = F(1, 2)
    points = [Point((Q3(F(1)), Q3())),
              Point((Q3(half), Q3(F(0), half))),
              Point((Q3(-half), Q3(F(0), half))),
              Point((Q3(F(-1)), Q3())),
              Point((Q3(-half), Q3(F(0), -half))),
              Point((Q3(half), Q3(F(0), -half))),
              Point((Q3(), Q3()))]
    cycle = {tuple(sorted((i, (i + 1) % 6))) for i in range(6)}
    required = cycle | {(i, 6) for i in range(6)}
    for i, j in required:
        need(qnorm(tuple(a - b for a, b in zip(points[i].base, points[j].base))) == Q3(F(1)),
             ("base edge identity", i, j))
    records = []
    ear_vertices = {}
    for u, v in combinations(range(6), 2):
        new, record = ear(points, u, v)
        x, y = len(points), len(points) + 1
        points.extend(new)
        required.update({(u, x), (x, y), (v, y)})
        ear_vertices[u, v] = (x, y)
        record["vertices"] = [x, y]
        records.append(record)
    edges, nonedges = unit_graph(points, required)
    need(len(points) == 37 and len(edges) == 57 and nonedges == 609, "wheel-ear universe")

    # A literal 3-colouring, hence also 5- and 6-colourings, of the full37-point graph.
    colouring = [i % 2 for i in range(6)] + [2] + [-1] * 30
    for (u, v), (x, y) in ear_vertices.items():
        colouring[x], colouring[y] = next((a, b) for a, b in product(range(3), repeat=2)
                                          if a != colouring[u] and b != colouring[v] and a != b)
    need(all(colouring[u] != colouring[v] for u, v in edges), "literal 3-colouring")

    # Full six-colour wheel boundary, no discarded equality states.
    q, n = 6, 6
    P = partitions(n)
    need(len(P) == 203, "Bell6")
    weights = {pi: q - len(set(pi)) if all(pi[u] != pi[v] for u, v in cycle) else 0 for pi in P}
    masses = {pi: falling(q, len(set(pi))) * w for pi, w in weights.items()}
    K = q * q - 3 * q + 3
    queries = {}
    for pi in P:
        queries[pi] = sum(masses[beta] * prod(K - int(beta[u] == beta[v]) for u, v in forest(pi))
                          for beta in P)
    moments, recovered = recover(P, queries, q)
    need(masses == recovered, "full partition reconstruction")
    need(sum(masses.values()) == 24600, "wheel total")

    # An independent full46656-assignment route with literal colour-pair loops.
    direct_weights = {pi: 0 for pi in P}
    direct_queries = {pi: 0 for pi in P}
    ear_kernel = {(a, b): sum(x != a and y != b and x != y for x, y in product(range(q), repeat=2))
                  for a, b in product(range(q), repeat=2)}
    for beta in product(range(q), repeat=n):
        if any(beta[u] == beta[v] for u, v in cycle):
            continue
        completions = sum(all(c != b for b in beta) for c in range(q))
        canonical = partition(beta)
        direct_weights[canonical] += completions
        for pi in P:
            direct_queries[pi] += completions * prod(ear_kernel[beta[u], beta[v]] for u, v in forest(pi))
    need(direct_weights == masses and direct_queries == queries, "literal colour-table route")
    support_assignments = sum(falling(q, len(set(pi))) for pi in P if weights[pi])
    need(support_assignments == 14910, "wheel boundary support")

    # Sharp palette failure at q2: induced P3 plus a realizable 3-edge ear is C5.
    hostile = [Point((Q3(F(-1)), Q3())), Point((Q3(F(1)), Q3())), Point((Q3(), Q3()))]
    new, hostile_record = ear(hostile, 0, 1)
    hostile.extend(new)
    hostile_edges = {(0, 2), (1, 2), (0, 3), (3, 4), (1, 4)}
    unit_graph(hostile, hostile_edges)
    need(proper_count(3, {(0, 2), (1, 2)}, 2) == 2 and proper_count(5, hostile_edges, 2) == 0,
         "q2 forcing hostile")
    need(proper_count(5, hostile_edges, 3) == 30, "q3 positive hostile")

    # Literal independent total counts on one wheel ear, three metric pair types.
    small = []
    for uv in ((0, 1), (0, 2), (0, 3)):
        es = cycle | {(i, 6) for i in range(6)} | {(uv[0], 7), (7, 8), (uv[1], 8)}
        actual = proper_count(9, es, 6)
        expected = sum(masses[b] * (K - int(b[uv[0]] == b[uv[1]])) for b in P)
        need(actual == expected, ("one-ear backtracking", uv, actual, expected))
        small.append(dict(terminals=uv, count=actual))

    tables = [dict(partition=list(pi), weight=weights[pi], mass=masses[pi],
                   query=queries[pi], equality_moment=moments[pi], recovered_mass=recovered[pi]) for pi in P]
    coordinates = dict(field="Q(sqrt(3),sqrt(D_e) for each ear); D_e in Q(sqrt(3))",
                       point_format="base + coefficient*sqrt(radicand); [a,b] means a+b*sqrt(3)",
                       interval_precision_bits=BITS, points=[p.serial() for p in points],
                       unit_edges=sorted(map(list, edges)), ears=records, colouring3=colouring,
                       q2_hostile=dict(points=[p.serial() for p in hostile],
                                       edges=sorted(map(list, hostile_edges)), construction=hostile_record))
    receipt = dict(points=37, unit_edges=57, certified_nonunit_pairs=609,
                   independent_literal_assignments=6 ** 6, partitions=203,
                   positive_boundary_partitions=sum(w > 0 for w in weights.values()),
                   full_boundary_support=support_assignments, base_colourings6=sum(masses.values()),
                   all203_queries_recover_exact_table=True, literal_colour_loop_agreement=True,
                   query_min=min(queries.values()), query_max=max(queries.values()),
                   witness_colours=3, one_ear_counts=small,
                   q2_hostile_base=2, q2_hostile_ear=0, q3_hostile_ear=30,
                   max_construction_attempts=max(r["attempts"] for r in records),
                   inherited_source_pins=len(source_pins),
                   source_sha256=hashlib.sha256(source).hexdigest())
    need(Path(__file__).read_bytes() == source, "source changed during verification")
    need(all((root_path / name).read_bytes() == body for name, body in inherited_bytes.items()),
         "inherited source changed during verification")
    if "--write" in sys.argv:
        Path("coordinates.json").write_text(json.dumps(coordinates, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
        Path("boundary_tables.json").write_text(json.dumps(tables, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    else:
        need(json.loads(Path("coordinates.json").read_text(encoding="utf-8")) == coordinates,
             "saved exact coordinates disagree")
        need(json.loads(Path("boundary_tables.json").read_text(encoding="utf-8")) == tables,
             "saved complete boundary table disagrees")
    print(json.dumps(receipt, indent=2, sort_keys=True))


def prod(items):
    result = 1
    for item in items:
        result *= item
    return result


if __name__ == "__main__":
    main()
