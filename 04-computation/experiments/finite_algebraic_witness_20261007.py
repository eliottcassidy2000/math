"""Exact algebraic geometry, colouring refutations and triangle-packing certificates.

This verifies small controls, not a six-chromatic plane graph.  No accepted
external theorem's proof is tested here; no floating tolerances are used.
"""
from fractions import Fraction as F
from itertools import combinations, permutations
from collections import deque
import json


def need(ok, why):
    if not ok:
        raise ValueError(why)


def poly(a):
    need(type(a) is tuple and a and all(type(x) in (int, F) for x in a),
         "nonempty tuple of exact rational coefficients required")
    a = tuple(F(x) for x in a)
    while len(a) > 1 and a[-1] == 0:
        a = a[:-1]
    return a


ZERO, ONE = poly((0,)), poly((1,))


def add(a, b):
    n = max(len(a), len(b))
    return poly(tuple((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
                      for i in range(n)))


def scale(a, c):
    need(type(c) in (int, F), "exact rational scale required")
    return poly(tuple(c * x for x in a))


def sub(a, b):
    return add(a, scale(b, -1))


def mul(a, b):
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return poly(tuple(out))


def divide(a, b):
    need(b != ZERO, "zero polynomial divisor")
    r = list(a)
    q = [F(0)] * max(1, len(a) - len(b) + 1)
    while len(r) >= len(b) and poly(tuple(r)) != ZERO:
        k = len(r) - len(b)
        c = r[-1] / b[-1]
        q[k] += c
        for j, x in enumerate(b):
            r[k + j] -= c * x
        r = list(poly(tuple(r)))
    return poly(tuple(q)), poly(tuple(r))


def derivative(a):
    return poly(tuple(i * a[i] for i in range(1, len(a))) or (0,))


def gcd(a, b):
    while b != ZERO:
        a, b = b, divide(a, b)[1]
    return scale(a, 1 / a[-1]) if a != ZERO else ZERO


def evaluate(a, x):
    out = F(0)
    for c in reversed(a):
        out = out * x + c
    return out


def sturm(a):
    seq = [a, derivative(a)]
    if seq[-1] == ZERO:
        return (a,)
    while True:
        r = scale(divide(seq[-2], seq[-1])[1], -1)
        if r == ZERO:
            return tuple(seq)
        seq.append(r)


def changes(values):
    signs = [1 if x > 0 else -1 for x in values if x]
    return sum(x != y for x, y in zip(signs, signs[1:]))


def root_count(a, lo, hi):
    need(lo < hi and evaluate(a, lo) != 0 and evaluate(a, hi) != 0,
         "nonroot rational endpoints required")
    seq = sturm(a)
    return changes([evaluate(p, lo) for p in seq]) - changes([evaluate(p, hi) for p in seq])


class AlgebraicReal:
    """A selected real root, with exact zero tests even for reducible f."""
    def __init__(self, coefficients, lo, hi):
        self.f = poly(coefficients)
        need(len(self.f) >= 2 and self.f[-1] == 1, "monic nonconstant polynomial required")
        need(type(lo) in (int, F) and type(hi) in (int, F), "rational interval required")
        self.lo, self.hi = F(lo), F(hi)
        need(len(gcd(self.f, derivative(self.f))) == 1, "squarefree polynomial required")
        need(root_count(self.f, self.lo, self.hi) == 1, "interval must select exactly one real root")

    def reduce(self, p):
        return divide(poly(p), self.f)[1]

    def zero(self, p):
        p = self.reduce(p)
        if p == ZERO:
            return True
        g = gcd(self.f, p)
        return len(g) > 1 and root_count(g, self.lo, self.hi) == 1

    def equal(self, p, q):
        return self.zero(sub(p, q))

    def times(self, p, q):
        return self.reduce(mul(p, q))


def cmul(a, b, alg):
    return (sub(alg.times(a[0], b[0]), alg.times(a[1], b[1])),
            add(alg.times(a[0], b[1]), alg.times(a[1], b[0])))


def same_point(a, b, alg):
    return alg.equal(a[0], b[0]) and alg.equal(a[1], b[1])


def distance_squared(a, b, alg):
    x, y = sub(a[0], b[0]), sub(a[1], b[1])
    return add(alg.times(x, x), alg.times(y, y))


def unit_graph(points, alg):
    need(type(points) is tuple, "tuple of points required")
    need(all(type(p) is tuple and len(p) == 2 for p in points), "each point must have exactly two coordinates")
    pts = tuple((poly(p[0]), poly(p[1])) for p in points)
    edges = set()
    for i, j in combinations(range(len(pts)), 2):
        need(not same_point(pts[i], pts[j], alg), "coincident vertices rejected")
        if alg.equal(distance_squared(pts[i], pts[j], alg), ONE):
            edges.add((i, j))
    return frozenset(edges)


def graph_checked(n, edges):
    need(type(n) is int and n >= 0, "natural vertex count required")
    need(type(edges) in (set, frozenset), "edge set required")
    need(all(type(e) is tuple and len(e) == 2 and all(type(x) is int for x in e)
             and 0 <= e[0] < e[1] < n for e in edges), "ordered distinct vertex indices required")
    return tuple(frozenset(v if u == i else u for u, v in edges if i in (u, v)) for i in range(n))


def colour_or_refute(n, edges, colours):
    """Finite search; returned refutation is independently checkable below."""
    adj = graph_checked(n, edges)
    need(type(colours) is int and colours >= 1, "positive exact colour count required")

    def visit(a):
        for u, v in sorted(edges):
            if u in a and v in a and a[u] == a[v]:
                return None, ("edge", u, v)
        if len(a) == n:
            return tuple(a[i] for i in range(n)), None
        v = max((i for i in range(n) if i not in a),
                key=lambda i: (len({a[j] for j in adj[i] if j in a}), len(adj[i]), -i))
        children = []
        for c in range(colours):
            witness, proof = visit(a | {v: c})
            if witness is not None:
                return witness, None
            children.append(proof)
        return None, ("split", v, tuple(children))

    return visit({})


def verify_refutation(n, edges, colours, proof):
    graph_checked(n, edges)
    need(type(colours) is int and colours >= 1, "positive colour count required")

    def visit(node, assigned):
        need(type(node) is tuple and len(node) == 3, "proof node required")
        tag, a, b = node
        if tag == "edge":
            need(type(a) is int and type(b) is int and (a, b) in edges,
                 "leaf must name an actual edge")
            need(a in assigned and b in assigned and assigned[a] == assigned[b],
                 "leaf must be an assigned monochromatic edge")
            return 1
        need(tag == "split" and type(a) is int and 0 <= a < n and a not in assigned,
             "split must select a fresh vertex")
        need(type(b) is tuple and len(b) == colours, "all colours must be covered")
        return 1 + sum(visit(child, assigned | {a: c}) for c, child in enumerate(b))

    return visit(proof, {})


def verify_colouring(n, edges, colours, assignment):
    graph_checked(n, edges)
    need(type(colours) is int and colours >= 1, "positive colour count required")
    need(type(assignment) is tuple and len(assignment) == n and
         all(type(c) is int and 0 <= c < colours for c in assignment), "colouring shape invalid")
    need(all(assignment[u] != assignment[v] for u, v in edges), "monochromatic edge")
    return True


def triangles_of(n, edges):
    graph_checked(n, edges)
    return tuple(t for t in combinations(range(n), 3)
                 if all(e in edges for e in combinations(t, 2)))


def conflict_graph(triangles):
    need(type(triangles) is tuple, "tuple of triangles required")
    need(all(type(t) is tuple and len(t) == 3 and all(type(v) is int and v >= 0 for v in t)
             and t[0] < t[1] < t[2] for t in triangles), "ordered triangle indices required")
    need(len(set(triangles)) == len(triangles), "distinct triangles required")
    return frozenset((i, j) for i, j in combinations(range(len(triangles)), 2)
                     if len(set(triangles[i]) & set(triangles[j])) == 2)


def centred_cubic(triangle, points, alg):
    center = tuple(scale(add(add(points[triangle[0]][j], points[triangle[1]][j]),
                             points[triangle[2]][j]), F(1, 3)) for j in (0, 1))
    result = (ONE, ZERO)
    for v in triangle:
        result = cmul(result, (sub(points[v][0], center[0]), sub(points[v][1], center[1])), alg)
    return result


def packing_certificate(n, edges):
    triangles = triangles_of(n, edges)
    conflicts = conflict_graph(triangles)
    adj = graph_checked(len(triangles), conflicts)
    side = {}
    for start in range(len(triangles)):
        if start in side:
            continue
        side[start] = 0
        queue = deque([start])
        while queue:
            v = queue.popleft()
            for u in sorted(adj[v]):
                if u not in side:
                    side[u] = 1 - side[v]
                    queue.append(u)
                need(side[u] != side[v], "triangle conflicts are not bipartite")
    left = {v for v in side if side[v] == 0}
    right_match = {}

    def augment(v, seen):
        for u in sorted(adj[v]):
            if u in seen:
                continue
            seen.add(u)
            if u not in right_match or augment(right_match[u], seen):
                right_match[u] = v
                return True
        return False

    for v in sorted(left):
        augment(v, set())
    matching = frozenset(tuple(sorted((v, u))) for u, v in right_match.items())
    matched_vertices = {x for e in matching for x in e}
    reached = set(left - matched_vertices)
    queue = deque(sorted(reached))
    while queue:
        v = queue.popleft()
        for u in sorted(adj[v]):
            edge = tuple(sorted((v, u)))
            allowed = (side[v] == 0 and edge not in matching) or (side[v] == 1 and edge in matching)
            if allowed and u not in reached:
                reached.add(u)
                queue.append(u)
    cover = (left - reached) | (reached - left)
    packing = set(range(len(triangles))) - cover
    return triangles, conflicts, frozenset(packing), matching, frozenset(cover)


def verify_packing_certificate(vertices, edges, packet):
    graph_checked(vertices, edges)
    need(type(packet) is tuple and len(packet) == 5, "packing packet required")
    triangles, conflicts, packing, matching, cover = packet
    computed_conflicts = conflict_graph(triangles)
    need(all(v < vertices for t in triangles for v in t), "triangle vertices out of range")
    need(triangles == triangles_of(vertices, edges), "triangle list must be complete for the source graph")
    need(all(type(s) is frozenset for s in (conflicts, packing, matching, cover)), "frozen incidence sets required")
    n = len(triangles)
    graph_checked(n, conflicts)
    graph_checked(n, matching)
    need(conflicts == computed_conflicts, "wrong conflict incidence")
    need(all(type(i) is int and 0 <= i < n for i in packing | cover), "bad triangle index")
    need(packing == frozenset(range(n)) - cover, "packing and cover must be complementary")
    need(all(u in cover or v in cover for u, v in conflicts), "cover misses a conflict")
    need(matching <= conflicts and len({v for e in matching for v in e}) == 2 * len(matching),
         "not a matching of actual conflicts")
    need(len(cover) == len(matching), "primal/dual sizes differ")
    used = []
    for i in packing:
        used.extend(combinations(triangles[i], 2))
    need(len(used) == len(set(used)), "packing repeats a geometric edge")
    return True


def main():
    checks = 0
    def check(ok, label):
        nonlocal checks
        need(ok, label)
        checks += 1

    alg = AlgebraicReal((64, 0, -28, 0, 1), 5, 6)
    root = poly((0, 1))
    s3 = scale(poly((0, -20, 0, 1)), F(1, 16))
    s11 = scale(poly((0, 36, 0, -1)), F(1, 16))
    s33 = alg.times(s3, s11)
    check(alg.equal(add(s3, s11), root), "selected primitive element")
    check(alg.equal(alg.times(s3, s3), poly((3,))), "square root three")
    check(alg.equal(alg.times(s11, s11), poly((11,))), "square root eleven")
    points = ((ZERO, ZERO), (s3, ZERO), (scale(s3, F(1, 2)), poly((F(1, 2),))),
              (scale(s3, F(1, 2)), poly((F(-1, 2),))),
              (scale(s3, F(5, 6)), scale(s33, F(1, 6))),
              (scale(sub(scale(s3, 5), s11), F(1, 12)), scale(add(s33, poly((5,))), F(1, 12))),
              (scale(add(scale(s3, 5), s11), F(1, 12)), scale(sub(s33, poly((5,))), F(1, 12))))
    edges = unit_graph(points, alg)
    check(len(edges) == 11, "Moser spindle exact edge count")
    colouring, proof = colour_or_refute(7, edges, 3)
    check(colouring is None, "Moser spindle rejects three colours")
    proof_nodes = verify_refutation(7, edges, 3, proof)
    check(proof_nodes > 1, "independent full branching refutation")
    four, no_proof = colour_or_refute(7, edges, 4)
    check(no_proof is None and verify_colouring(7, edges, 4, four), "four-colour witness")
    for deleted in range(7):
        keep = [i for i in range(7) if i != deleted]
        smaller = frozenset((keep.index(u), keep.index(v)) for u, v in edges if deleted not in (u, v))
        witness, _ = colour_or_refute(6, smaller, 3)
        check(witness is not None and verify_colouring(6, smaller, 3, witness), "vertex-critical control")

    # All subsets of a 3-by-3 triangular-lattice patch: exact geometric and
    # independent combinatorial packing controls, not a six-chromatic search.
    patch = tuple((poly((F(i) + F(j, 2),)), scale(s3, F(j, 2)))
                  for i in range(3) for j in range(3))
    full_edges = unit_graph(patch, alg)
    patch_controls, total_triangles = 0, 0
    for mask in range(1 << len(patch)):
        keep = [i for i in range(len(patch)) if mask & (1 << i)]
        renaming = {v: i for i, v in enumerate(keep)}
        e = frozenset((renaming[u], renaming[v]) for u, v in full_edges if u in renaming and v in renaming)
        packet = packing_certificate(len(keep), e)
        t, c, pack, matching, cover = packet
        check(verify_packing_certificate(len(keep), e, packet), "matching/cover exact certificate")
        # Independent exhaustive packing sizes (the full patch has 8 triangles).
        best = max((m.bit_count() for m in range(1 << len(t))
                    if all(not (m & (1 << u) and m & (1 << v)) for u, v in c)), default=0)
        check(len(pack) == best, "maximum versus independently enumerated packing")
        check(max((sum(v in edge for edge in c) for v in range(len(t))), default=0) <= 3,
              "geometric conflict degree bound")
        patch_controls += 1
        total_triangles += len(t)

    for pts in (points, patch):
        e = unit_graph(pts, alg)
        packet = packing_certificate(len(pts), e)
        check(verify_packing_certificate(len(pts), e, packet), "full geometry packing")
        t, conflicts, *_ = packet
        etas = [centred_cubic(tri, pts, alg) for tri in t]
        for eta in etas:
            check(not (alg.zero(eta[0]) and alg.zero(eta[1])), "nonzero centered cubic")
        for u, v in conflicts:
            check(alg.zero(add(etas[u][0], etas[v][0])) and alg.zero(add(etas[u][1], etas[v][1])),
                  "shared edge negates unordered complex cubic")
        for tri in t:
            for perm in permutations(tri):
                check(same_point(centred_cubic(perm, pts, alg), centred_cubic(tri, pts, alg), alg),
                      "cubic independent of vertex ordering")

    # Three triangles in a chain: greedy maximal need not be maximum.
    omega = (poly((F(1, 2),)), scale(s3, F(1, 2)))
    chain = ((ZERO, ZERO), (ONE, ZERO), omega,
             (add(ONE, omega[0]), omega[1]), (poly((2,)), ZERO))
    chain_edges = unit_graph(chain, alg)
    packet = packing_certificate(5, chain_edges)
    triangles, conflicts, pack, matching, cover = packet
    check(len(chain_edges) == 7 and len(triangles) == 3 and len(conflicts) == 2, "three-triangle chain")
    middle = next(v for v in range(3) if sum(v in edge for edge in conflicts) == 2)
    check(all(middle in edge for edge in conflicts) and len(pack) == 2, "greedy one versus optimal two")
    check(len(chain_edges) - 3 == 4 and len(chain_edges) - 3 * len(pack) == 1, "terminal leaves four versus one")

    # One convex hull and volume product cannot determine chromatic number.
    square = tuple((poly((x,)), poly((y,))) for x, y in ((-10,-10),(-10,10),(10,-10),(10,10)))
    square_edges = unit_graph(square, alg)
    with_spindle = unit_graph(square + points, alg)
    check(not square_edges and len(with_spindle) == 11, "same square hull, chromatic one versus four")
    check(F(400) * F(1, 50) == 8 > F(27, 4), "square/polar product and planar Mahler bound")

    # A zero test must use the selected root, not a nonzero quotient remainder.
    reducible = AlgebraicReal((2, 0, -3, 0, 1), F(5, 4), F(3, 2))
    check(reducible.zero(poly((-2, 0, 1))), "selected-root zero with nonzero polynomial remainder")
    two_edges = frozenset(((0, 1), (0, 2), (1, 2), (0, 3), (1, 3)))
    tt, tc, tp, tm, tv = packing_certificate(4, two_edges)
    hostiles = [lambda: AlgebraicReal((64, 0, -28, 0, 1), -6, 6),
                lambda: AlgebraicReal((1, -2, 1), 0, 2),
                lambda: AlgebraicReal((True, 1), -2, 0),
                lambda: unit_graph((points[0], points[0]), alg),
                lambda: verify_refutation(7, edges, 3, ("split", 0, (proof,))),
                lambda: verify_refutation(7, edges, 3, ("edge", 0, 1)),
                lambda: verify_colouring(7, edges, 4, (0,) * 7),
                lambda: packing_certificate(4, frozenset(combinations(range(4), 2))),
                lambda: verify_packing_certificate(5, chain_edges,
                                                  (triangles, conflicts, frozenset(range(3)),
                                                   frozenset(), frozenset())),
                lambda: verify_packing_certificate(5, chain_edges,
                                                  (triangles[:-1], frozenset(), frozenset((0, 1)),
                                                   frozenset(), frozenset())),
                lambda: unit_graph((((0,), (0,), (7,)), ((1,), (0,), (8,))), alg),
                lambda: verify_packing_certificate(4, two_edges, (tt, tc, tp, frozenset(((False, True),)), tv)),
                lambda: verify_packing_certificate(4, two_edges, (tt, tc, tp, frozenset(((0.0, 1.0),)), tv)),
                lambda: verify_packing_certificate(4, two_edges, (((False, 1, 2), (0, 1, 3)), tc, tp, tm, tv)),
                lambda: verify_packing_certificate(4, two_edges, (tt, frozenset(((False, True),)), tp, tm, tv))]
    for action in hostiles:
        try:
            action()
        except ValueError:
            check(True, "malformed geometry/proof/packing rejected")
        else:
            raise ValueError("hostile unexpectedly accepted")

    print(json.dumps({"status": "Exact certificate controls; no six-chromatic witness claimed",
                      "checks": checks, "primitive_polynomial": [64, 0, -28, 0, 1],
                      "isolating_interval": [5, 6], "moser_vertices": 7, "moser_edges": 11,
                      "moser_chromatic_number": 4, "refutation_nodes": proof_nodes,
                      "moser_four_colouring": four, "lattice_subsets": patch_controls,
                      "triangle_occurrences_in_subsets": total_triangles,
                      "greedy_chain_leave": 4, "optimal_chain_leave": 1,
                      "typed_geometry_hostiles": len(hostiles)}, sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
