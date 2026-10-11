"""Checked vertex-critical extraction. No six-chromatic example is supplied."""
from dataclasses import dataclass, replace
from fractions import Fraction as F
from itertools import combinations, product
import json

import finite_algebraic_witness_20261007 as exact


def need(ok, message):
    if not ok:
        raise ValueError(message)


def induced(n, edges, keep):
    exact.graph_checked(n, edges)
    need(type(keep) is tuple and all(type(v) is int and 0 <= v < n for v in keep),
         "exact vertex labels required")
    need(tuple(sorted(set(keep))) == keep, "increasing distinct labels required")
    index = {v: i for i, v in enumerate(keep)}
    return frozenset((index[u], index[v]) for u, v in edges if u in index and v in index)


@dataclass(frozen=True)
class Critical:
    n: int
    edges: frozenset
    colours: int
    refutation: tuple
    colouring: tuple
    deletions: tuple


def verify_critical(packet):
    need(type(packet) is Critical, "critical packet required")
    n, edges, k = packet.n, packet.edges, packet.colours
    exact.graph_checked(n, edges)
    need(type(k) is int and k >= 2 and n >= 2, "nontrivial exact palette required")
    need(type(edges) is frozenset, "frozen edges required")
    nodes = exact.verify_refutation(n, edges, k - 1, packet.refutation)
    exact.verify_colouring(n, edges, k, packet.colouring)
    need(type(packet.deletions) is tuple and len(packet.deletions) == n,
         "one marked colouring for every vertex deletion required")
    for v, colouring in enumerate(packet.deletions):
        keep = tuple(u for u in range(n) if u != v)
        exact.verify_colouring(n - 1, induced(n, edges, keep), k - 1, colouring)
    return nodes


def extract_critical(n, edges, k):
    """Return original vertex labels and an induced vertex-critical packet."""
    exact.graph_checked(n, edges)
    need(type(k) is int and k >= 2, "exact palette at least two required")
    colour, _ = exact.colour_or_refute(n, edges, k - 1)
    need(colour is None, "input already admits the smaller palette")
    keep = tuple(range(n))
    # A single pass suffices: colourability persists under further deletion.
    for v in range(n):
        trial = tuple(u for u in keep if u != v)
        colour, _ = exact.colour_or_refute(len(trial), induced(n, edges, trial), k - 1)
        if colour is None:
            keep = trial
    kept_edges = induced(n, edges, keep)
    colour, proof = exact.colour_or_refute(len(keep), kept_edges, k - 1)
    need(colour is None, "retained obstruction")
    upper, _ = exact.colour_or_refute(len(keep), kept_edges, k)
    need(upper is not None, "vertex minimality gives one extra colour")
    deletions = []
    for v in range(len(keep)):
        rest = tuple(u for u in range(len(keep)) if u != v)
        c, _ = exact.colour_or_refute(len(rest), induced(len(keep), kept_edges, rest), k - 1)
        need(c is not None, "every deletion must be colourable")
        deletions.append(c)
    result = Critical(len(keep), kept_edges, k, proof, upper, tuple(deletions))
    verify_critical(result)
    return keep, result


@dataclass(frozen=True)
class AlgebraicWitness:
    polynomial: tuple
    lower: F
    upper: F
    points: tuple
    critical: Critical


def verify_algebraic(packet):
    need(type(packet) is AlgebraicWitness, "algebraic witness required")
    alg = exact.AlgebraicReal(packet.polynomial, packet.lower, packet.upper)
    edges = exact.unit_graph(packet.points, alg)
    nodes = verify_critical(packet.critical)
    need(len(packet.points) == packet.critical.n and edges == packet.critical.edges,
         "complete induced unit graph must agree with the colouring certificate")
    return nodes


def extract_algebraic(polynomial, lower, upper, points, k):
    alg = exact.AlgebraicReal(polynomial, lower, upper)
    edges = exact.unit_graph(points, alg)
    keep, critical = extract_critical(len(points), edges, k)
    packet = AlgebraicWitness(alg.f, alg.lo, alg.hi,
                              tuple(points[v] for v in keep), critical)
    verify_algebraic(packet)
    return keep, packet


def first_witness(candidates, k):
    """Consume a caller-supplied stream; completeness/fairness is its premise."""
    need(type(k) is int and k >= 2, "exact palette at least two required")
    for polynomial, lo, hi, points in candidates:
        alg = exact.AlgebraicReal(polynomial, lo, hi)
        edges = exact.unit_graph(points, alg)
        c, _ = exact.colour_or_refute(len(points), edges, k - 1)
        if c is None:
            return extract_algebraic(polynomial, lo, hi, points, k)
    return None


def colouring_cnf(n, edges, k):
    exact.graph_checked(n, edges)
    need(type(k) is int and k >= 1, "positive exact palette required")
    var = lambda v, c: k * v + c + 1
    clauses = []
    for v in range(n):
        clauses.append(tuple(var(v, c) for c in range(k)))
        clauses.extend((-var(v, c), -var(v, d)) for c, d in combinations(range(k), 2))
    for u, v in sorted(edges):
        clauses.extend((-var(u, c), -var(v, c)) for c in range(k))
    return tuple(clauses)


def moser_points():
    z = exact.ZERO
    s = exact.scale(exact.poly((0, -20, 0, 1)), F(1, 16))
    t = exact.scale(exact.poly((0, 36, 0, -1)), F(1, 16))
    alg = exact.AlgebraicReal((64, 0, -28, 0, 1), 5, 6)
    st = alg.times(s, t)
    return (
        (z, z), (s, z),
        (exact.scale(s, F(1, 2)), exact.poly((F(1, 2),))),
        (exact.scale(s, F(1, 2)), exact.poly((F(-1, 2),))),
        (exact.scale(s, F(5, 6)), exact.scale(st, F(1, 6))),
        (exact.scale(exact.sub(exact.scale(s, 5), t), F(1, 12)),
         exact.scale(exact.add(st, exact.poly((5,))), F(1, 12))),
        (exact.scale(exact.add(exact.scale(s, 5), t), F(1, 12)),
         exact.scale(exact.sub(st, exact.poly((5,))), F(1, 12))),
    )


def main():
    checks = 0
    def check(ok, label):
        nonlocal checks
        need(ok, label)
        checks += 1
    def rejects(fn):
        nonlocal checks
        try:
            fn()
        except (ValueError, TypeError):
            checks += 1
            return
        raise ValueError("malformed packet accepted")

    points = moser_points()
    keep, spindle = extract_algebraic((64, 0, -28, 0, 1), 5, 6, points, 4)
    check(keep == tuple(range(7)), "all spindle vertices essential")
    nodes = verify_algebraic(spindle)
    check(nodes == 139 and len(spindle.critical.edges) == 11, "inherited exact spindle")
    extra = points + ((exact.poly((10,)), exact.poly((10,))),)
    keep2, spindle2 = extract_algebraic((64, 0, -28, 0, 1), 5, 6, extra, 4)
    check(keep2 == keep and spindle2 == spindle, "irrelevant isolated point removed")
    triangle = points[:1] + (points[2], points[3])
    stream = (((64, 0, -28, 0, 1), 5, 6, triangle),
              ((64, 0, -28, 0, 1), 5, 6, extra))
    check(first_witness(stream, 4) == (keep, spindle), "finite stream search and source indices")
    check(first_witness(stream[:1], 4) is None, "finite miss is not global impossibility")

    # Complete graph universe on one through five labelled vertices, at k=3.
    graph_count = obstruction_count = 0
    critical_sizes = set()
    for n in range(1, 6):
        possible = tuple(combinations(range(n), 2))
        for mask in range(1 << len(possible)):
            edges = frozenset(e for i, e in enumerate(possible) if mask >> i & 1)
            graph_count += 1
            c, _ = exact.colour_or_refute(n, edges, 2)
            if c is not None:
                continue
            obstruction_count += 1
            retained, packet = extract_critical(n, edges, 3)
            check(packet.edges == induced(n, edges, retained), "induced source sidecar")
            check(verify_critical(packet) > 0, "independent critical proof check")
            critical_sizes.add(packet.n)
            # Independent full truth table on every proper induced subgraph.
            for subset in range((1 << packet.n) - 1):
                rest = tuple(v for v in range(packet.n) if subset >> v & 1)
                re = induced(packet.n, packet.edges, rest)
                found = any(all(a[u] != a[v] for u, v in re)
                            for a in product(range(2), repeat=len(rest)))
                check(found, "every proper induced subgraph is two-colourable")
    check(graph_count == 1099 and critical_sizes == {3, 5}, "complete small critical universe")

    cnf = colouring_cnf(7, spindle.critical.edges, 3)
    check(len(cnf) == 61, "explicit one-hot CNF size")
    cnf_models = 0
    for a in product(range(3), repeat=7):
        true_vars = {3 * v + a[v] + 1 for v in range(7)}
        sat = all(any((lit in true_vars) if lit > 0 else (-lit not in true_vars)
                      for lit in clause) for clause in cnf)
        proper = all(a[u] != a[v] for u, v in spindle.critical.edges)
        check(sat == proper, "CNF versus direct colouring truth")
        cnf_models += sat
    check(cnf_models == 0, "independent 2187-assignment spindle refutation")

    p = spindle.critical
    rejects(lambda: verify_critical(replace(p, colours=True)))
    rejects(lambda: verify_critical(replace(p, n=7.0)))
    rejects(lambda: verify_critical(replace(p, colouring=(False,) + p.colouring[1:])))
    rejects(lambda: verify_critical(replace(p, deletions=p.deletions[:-1])))
    rejects(lambda: verify_critical(replace(p, deletions=((0,) * 6,) + p.deletions[1:])))
    rejects(lambda: verify_critical(replace(p, refutation=("split", 0, ()))))
    rejects(lambda: verify_critical(replace(p, edges=frozenset({(False, 1)}) | p.edges)))
    rejects(lambda: verify_algebraic(replace(spindle, points=points[:-1])))
    rejects(lambda: verify_algebraic(replace(spindle, points=(points[0],) + points[:-1])))
    rejects(lambda: verify_algebraic(replace(spindle, points=((exact.ZERO, exact.ZERO, exact.ONE),) + points[1:])))
    rejects(lambda: verify_algebraic(replace(spindle, lower=6, upper=7)))
    scaled = tuple((exact.scale(x, 2), exact.scale(y, 2)) for x, y in points)
    rejects(lambda: verify_algebraic(replace(spindle, points=scaled)))
    rejects(lambda: induced(7, p.edges, (False, 1)))
    rejects(lambda: induced(7, p.edges, (1, 0)))
    rejects(lambda: extract_critical(3, frozenset(), 3))
    rejects(lambda: colouring_cnf(7, p.edges, 3.0))
    rejects(lambda: first_witness((), True))

    print(json.dumps({
        "scope": "finite critical extraction and exact algebraic controls; no six-chromatic graph supplied",
        "checks": checks, "small_graphs": graph_count,
        "nonbipartite_small_graphs": obstruction_count,
        "critical_sizes": sorted(critical_sizes),
        "moser": {"vertices": 7, "edges": 11, "chromatic": 4,
                  "refutation_nodes": nodes, "deletion_colourings": 7,
                  "cnf_variables": 21, "cnf_clauses": len(cnf),
                  "one_hot_assignments": 3 ** 7},
        "malformed_controls": 17,
    }, indent=2))


if __name__ == "__main__":
    main()
