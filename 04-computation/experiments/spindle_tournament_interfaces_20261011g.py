"""Marked orientations, geometric equalities, and complete multiport interfaces.
The cited 509-vertex result is not loaded or reverified by this small packet.
"""
from dataclasses import replace
from fractions import Fraction as F
from itertools import combinations, product
import json

import boundary_palette_compression_20261011 as boundary
import partial_four_tournament_20261011c as tournament
import unit_equality_candidates_20261011e as geometry
import finite_algebraic_witness_20261007 as proof

CHECKS = 0


def need(ok, why):
    if not ok:
        raise ValueError(why)


def check(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(why)


def lift_witness(table, assignment):
    """Relabel a retained witness to a specified boundary colouring."""
    key = boundary.canonical(tuple(assignment[v] for v in table.boundary))
    witness = dict(table.entries)[key]
    old = dict(zip(table.graph.vertices, witness))
    transport = {old[v]: assignment[v] for v in table.boundary}
    free = iter(c for c in range(table.palette) if c not in transport.values())
    for c in range(table.palette):
        if c not in transport:
            transport[c] = next(free)
    return {v: transport[c] for v, c in old.items()}


def compose_interfaces(left, right, cross):
    """Complete relation for two marked pieces and declared boundary edges.

    This is an abstract graph operation. Geometric edges are checked separately.
    The inherited audit recomputes completeness, rather than trusting witnesses.
    """
    boundary.audit_table(left)
    boundary.audit_table(right)
    need(left.palette == right.palette, "same palette")
    shared = set(left.graph.vertices) & set(right.graph.vertices)
    need(shared <= set(left.boundary) & set(right.boundary),
         "every shared vertex must be retained by both interfaces")
    vertices = tuple(sorted(set(left.graph.vertices) | set(right.graph.vertices)))
    ports = tuple(sorted(set(left.boundary) | set(right.boundary)))
    boundary.audit_graph(boundary.Graph(ports, cross))
    graph = boundary.Graph(vertices, left.graph.edges | right.graph.edges | cross)
    lk, rk = boundary.keys(left), boundary.keys(right)
    entries = []
    for key in boundary.patterns(len(ports), left.palette):
        assignment = dict(zip(ports, key))
        if boundary.canonical(tuple(assignment[v] for v in left.boundary)) not in lk:
            continue
        if boundary.canonical(tuple(assignment[v] for v in right.boundary)) not in rk:
            continue
        if any(assignment[u] == assignment[v] for u, v in cross):
            continue
        a, b = lift_witness(left, assignment), lift_witness(right, assignment)
        need(all(a[v] == b[v] for v in shared), "matched shared colours")
        both = a | b
        witness = tuple(both[v] for v in vertices)
        need(boundary.valid_colouring(graph, left.palette, witness),
             "source graph witness survives all retained edges")
        entries.append((key, witness))
    return boundary.Table(graph, ports, left.palette, tuple(entries))


def pair_projection(table):
    """Only pairwise allowed equality bits; deliberately loses joint support."""
    boundary.audit_table(table)
    return tuple((table.boundary[i], table.boundary[j],
                  frozenset(int(key[i] == key[j]) for key, _ in table.entries))
                 for i, j in combinations(range(len(table.boundary)), 2))


def geometric_graph(points, edges=None):
    """Same ordered coordinate source; omitted unit edges may be explicit."""
    full = geometry.unit_edges(points)
    n = len(points)
    if edges is None:
        edges = full
    graph = boundary.Graph(tuple(range(n)), edges)
    boundary.audit_graph(graph)
    need(edges <= full, "every retained edge is an actual unit segment")
    return graph


def induced_piece(graph, vertices):
    boundary.audit_graph(graph)
    need(type(vertices) is tuple and all(type(v) is int for v in vertices)
         and tuple(sorted(set(vertices))) == vertices
         and set(vertices) <= set(graph.vertices), "actual ordered vertex subset")
    return boundary.Graph(vertices, frozenset(e for e in graph.edges if set(e) <= set(vertices)))


def nae_triangle():
    """Induced six-point unit graph, with ports 3,4,5."""
    z, one, a = geometry.ZERO, geometry.ONE, geometry.A
    points = (z, one, a, geometry.scale(-1, a),
              geometry.scale(2, one), geometry.scale(2, a))
    return points, geometric_graph(points)


def count_fixed(graph, k, ports, values):
    answer = 0
    for word in product(range(k), repeat=len(graph.vertices)):
        at = dict(zip(graph.vertices, word))
        if tuple(at[v] for v in ports) == values and boundary.valid_colouring(graph, k, word):
            answer += 1
    return answer


def rejects(call):
    try:
        call()
    except (ValueError, TypeError):
        return True
    return False


def main():
    z, one, a = geometry.ZERO, geometry.ONE, geometry.A
    rhombus = (z, one, geometry.add(one, a), a)
    diamond = geometric_graph(rhombus)
    square = geometric_graph(rhombus, tournament.edges(tournament.BASE))
    check(diamond.edges == square.edges | {(1, 3)}, "one actual unit diagonal")
    check(geometry.norm(geometry.sub(rhombus[0], rhombus[2])) == geometry.scale(3, one),
          "other diagonal has length sqrt(3)")
    orientation_rows = []
    for bits in product((0, 1), repeat=2):
        full = tournament.complete(tournament.BASE, bits)
        check(tournament.read(tournament.BASE, full) == bits, "same marked orientation")
        check(tournament.classify(full) == tournament.TYPES[bits], "four types")
        check(tournament.edges(full) == frozenset(combinations(range(4), 2)),
              "every tournament completion forgets to K4")
        check(rejects(lambda: geometric_graph(rhombus, tournament.edges(full))),
              "orientation completion is not a geometric unit completion")
        orientation_rows.append((bits, tournament.classify(full)))
    counts = []
    for k in range(3, 7):
        eq = count_fixed(diamond, k, (0, 2), (0, 0))
        diff = count_fixed(diamond, k, (0, 2), (0, 1))
        check(eq == (k-1)*(k-2), "diamond fixed-equal count")
        check(diff == (k-2)*(k-3), "diamond fixed-different count")
        counts.append((k, eq, diff))
    check(boundary.keys(boundary.table(square, (0, 2), 3)) == {(0, 0), (0, 1)},
          "square alone does not force its tips equal")
    check(boundary.keys(boundary.table(diamond, (0, 2), 3)) == {(0, 0)},
          "unit diagonal creates the three-colour equality")

    points = geometry.spindle()
    spindle = geometric_graph(points)
    left = induced_piece(spindle, (0, 1, 2, 3))
    right = induced_piece(spindle, (0, 4, 5, 6))
    check(len(spindle.vertices) == 7 and len(spindle.edges) == 11, "exact Moser graph")
    check(spindle.edges == left.edges | right.edges | {(3, 6)}, "complete cross-edge list")
    joined = {}
    for k in (3, 4, 6):
        lt, rt = boundary.table(left, (0, 3), k), boundary.table(right, (0, 6), k)
        result = compose_interfaces(lt, rt, frozenset({(3, 6)}))
        direct = boundary.table(spindle, (0, 3, 6), k)
        check(boundary.keys(result) == boundary.keys(direct), "spindle interface is exact")
        boundary.audit_table(result)
        joined[k] = sorted(boundary.keys(result))
        for _, witness in result.entries:
            check(boundary.valid_colouring(spindle, k, witness), "lift to exact spindle")
    check(not joined[3] and joined[4] and joined[6], "palette-specific obstruction")
    check(F(7, 8)**2 + F(15, 64) == 1, "cited Type M rotation is a unit rotation")
    check(4*(2-2*F(7, 8)) == 1, "Type M rotation closes radius-two tips")
    colouring, refutation = proof.colour_or_refute(7, spindle.edges, 3)
    check(colouring is None and proof.verify_refutation(7, spindle.edges, 3, refutation) > 0,
          "source-bound full non-three-colour proof")
    for edge in sorted(spindle.edges):
        reduced = spindle.edges - {edge}
        colouring, _ = proof.colour_or_refute(7, reduced, 3)
        check(proof.verify_colouring(7, reduced, 3, colouring), "all eleven edge deletions colour")
        check(colouring[edge[0]] == colouring[edge[1]], "deleted-edge equality witness")

    _, nae = nae_triangle()
    expected = frozenset({(0, 1), (0, 2), (1, 2), (0, 3), (1, 4), (2, 5)})
    check(nae.edges == expected, "NAE toy is fully induced unit graph")
    nt = boundary.table(nae, (3, 4, 5), 3)
    check(boundary.keys(nt) == set(boundary.patterns(3, 3)) - {(0, 0, 0)},
          "complete non-monochromatic triple support")
    for port_values in product(range(3), repeat=3):
        check(count_fixed(nae, 3, (3, 4, 5), port_values)
              == (0 if len(set(port_values)) == 1 else 2), "every labelled NAE extension count")
    empty_ports = boundary.Graph((3, 4, 5), frozenset())
    unconstrained = boundary.table(empty_ports, (3, 4, 5), 3)
    check(pair_projection(nt) == pair_projection(unconstrained), "all pair projections collide")
    check(boundary.keys(nt) != boundary.keys(unconstrained), "joint support does not collide")
    complete_witnesses = replace(nt, entries=nt.entries[:-1])
    check(rejects(lambda: compose_interfaces(complete_witnesses, unconstrained, frozenset())),
          "a list of positive witnesses is not a complete relation")

    # Independent finite composition universe, including overlapping ports.
    contexts = 0
    lp, rp = tuple(combinations((0, 1, 2), 2)), tuple(combinations((2, 3, 4), 2))
    for lm, rm, crossing, k in product(range(8), range(8), (False, True), range(2, 5)):
        l = boundary.Graph((0, 1, 2), frozenset(e for i, e in enumerate(lp) if lm >> i & 1))
        r = boundary.Graph((2, 3, 4), frozenset(e for i, e in enumerate(rp) if rm >> i & 1))
        result = compose_interfaces(boundary.table(l, (1, 2), k),
                                    boundary.table(r, (2, 3), k),
                                    frozenset({(1, 3)}) if crossing else frozenset())
        direct = boundary.table(result.graph, result.boundary, k)
        check(boundary.keys(result) == boundary.keys(direct), "complete five-vertex context")
        for _, witness in result.entries:
            check(boundary.valid_colouring(result.graph, k, witness), "context witness")
        contexts += 1
    lt, rt = boundary.table(left, (0, 3), 3), boundary.table(right, (0, 6), 3)
    hostiles = (
        lambda: compose_interfaces(lt, rt, frozenset({(False, 3)})),
        lambda: compose_interfaces(lt, rt, frozenset({(0.0, 3)})),
        lambda: compose_interfaces(lt, rt, frozenset({(1, 6)})),
        lambda: compose_interfaces(lt, boundary.table(right, (6,), 3), frozenset()),
        lambda: compose_interfaces(lt, boundary.table(right, (0, 6), 4), frozenset()),
        lambda: geometric_graph(rhombus, frozenset({(0, 2)})),
    )
    for hostile in hostiles:
        check(rejects(hostile), "malformed or unmarked interface rejected")
    print(json.dumps({
        "checks": CHECKS,
        "orientation_completions": orientation_rows,
        "diamond_fixed_colour_counts": counts,
        "spindle_vertices_edges": [7, 11],
        "spindle_boundary_support": joined,
        "nae_labelled_port_assignments": 24,
        "nae_extensions_each": 2,
        "nae_pairwise_summary_equals_unconstrained": True,
        "complete_join_contexts": contexts,
        "malformed_controls": len(hostiles)+1,
        "published_509_result_replayed": False,
        "new_six_colour_equality": False,
    }, indent=2))


if __name__ == "__main__":
    main()
