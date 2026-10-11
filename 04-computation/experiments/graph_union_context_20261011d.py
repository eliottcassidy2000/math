"""Exact graph-union receipts; no geometric realization or plane bound is inferred.

Run with python -B [-O] -X utf8.  All verification gates survive optimization.
"""
from collections import deque
from dataclasses import dataclass, replace
from itertools import combinations, product
from math import prod


CHECKS = 0


def need(ok, message):
    if not ok:
        raise ValueError(message)


def check(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)


def natural(n):
    need(type(n) is int and n >= 0, "exact natural number required")


@dataclass(frozen=True)
class Graph:
    n: int
    edges: frozenset


def validate_graph(g):
    need(type(g) is Graph, "Graph required")
    natural(g.n)
    need(type(g.edges) is frozenset, "canonical edge set required")
    for e in g.edges:
        need(type(e) is tuple and len(e) == 2, "edge must be a pair")
        u, v = e
        need(type(u) is int and type(v) is int and 0 <= u < v < g.n,
             "exact canonical endpoints required")
    return g


def graph(n, edges):
    natural(n)
    # Validate before set conversion: aliases such as False/0 must not disappear.
    records = tuple(edges)
    for e in records:
        need(type(e) is tuple and len(e) == 2, "edge must be a pair")
        u, v = e
        need(type(u) is int and type(v) is int and 0 <= u < v < n,
             "exact canonical endpoints required")
    return validate_graph(Graph(n, frozenset(records)))


def union_graph(graphs):
    need(type(graphs) is tuple and graphs, "nonempty tuple of graphs required")
    for g in graphs:
        validate_graph(g)
    n = graphs[0].n
    need(all(g.n == n for g in graphs), "shared labelled vertex set required")
    return Graph(n, frozenset().union(*(g.edges for g in graphs)))


def validate_colouring(g, colours, k):
    validate_graph(g)
    need(type(k) is int and k >= 1, "positive exact palette required")
    need(type(colours) is tuple and len(colours) == g.n, "one colour per vertex")
    need(all(type(c) is int and 0 <= c < k for c in colours), "colour type/range")
    need(all(colours[u] != colours[v] for u, v in g.edges), "improper colouring")


def colouring(g, k):
    """Complete finite search, returning one labelled colouring or None."""
    validate_graph(g)
    need(type(k) is int and k >= 1, "positive exact palette required")
    adj = [set() for _ in range(g.n)]
    for u, v in g.edges:
        adj[u].add(v)
        adj[v].add(u)
    order = sorted(range(g.n), key=lambda v: (-len(adj[v]), v))
    colours = [-1] * g.n

    def visit(i, largest):
        if i == g.n:
            return tuple(colours)
        v = order[i]
        forbidden = {colours[u] for u in adj[v] if colours[u] >= 0}
        for c in range(min(k - 1, largest + 1) + 1):
            if c not in forbidden:
                colours[v] = c
                result = visit(i + 1, max(largest, c))
                if result is not None:
                    return result
        colours[v] = -1
        return None

    return visit(0, -1)


def chromatic(g):
    if g.n == 0:
        return 0
    for k in range(1, g.n + 1):
        if colouring(g, k) is not None:
            return k
    raise RuntimeError("unreachable")


def clique_number(g):
    validate_graph(g)
    for size in range(g.n, 0, -1):
        if any(all(e in g.edges for e in combinations(s, 2))
               for s in combinations(range(g.n), size)):
            return size
    return 0


def girth(g):
    """None denotes a forest (infinite girth)."""
    validate_graph(g)
    adj = [set() for _ in range(g.n)]
    for u, v in g.edges:
        adj[u].add(v)
        adj[v].add(u)
    best = None
    for start in range(g.n):
        dist, parent = {start: 0}, {start: -1}
        queue = deque([start])
        while queue:
            u = queue.popleft()
            for v in adj[u]:
                if v not in dist:
                    dist[v], parent[v] = dist[u] + 1, u
                    queue.append(v)
                elif parent[u] != v:
                    length = dist[u] + dist[v] + 1
                    if best is None or length < best:
                        best = length
    return best


def owner(a, b):
    for r, (x, y) in enumerate(zip(a, b)):
        if x != y:
            return r
    return None


def validate_frame(g, palettes, labels):
    validate_graph(g)
    need(type(palettes) is tuple and palettes, "nonempty palette tuple required")
    need(all(type(k) is int and k >= 2 for k in palettes), "palette sizes >=2")
    need(type(labels) is tuple and len(labels) == g.n, "one product label per vertex")
    for label in labels:
        need(type(label) is tuple and len(label) == len(palettes), "product label arity")
        need(all(type(c) is int and 0 <= c < k for c, k in zip(label, palettes)),
             "product label type/range")
    need(all(labels[u] != labels[v] for u, v in g.edges), "source product colouring")


def required_girth(palettes):
    need(type(palettes) is tuple and palettes, "nonempty palette tuple required")
    need(all(type(k) is int and k >= 2 for k in palettes), "palette sizes >=2")
    if len(palettes) == 1:
        return None  # no mixed triangle exists
    return 1 + max(2 * (palettes[r] - 1) + palettes[s] - 1
                   for r in range(len(palettes)) for s in range(r + 1, len(palettes)))


@dataclass(frozen=True)
class UnionPacket:
    source: Graph
    palettes: tuple
    labels: tuple
    layers: tuple
    paths: tuple


def compile_union(g, palettes, labels):
    """First-difference orders with an original-edge path for every closure edge."""
    validate_frame(g, palettes, labels)
    layers, witnesses = [], []
    for r in range(len(palettes)):
        reach = {}
        for u, v in sorted(g.edges):
            if owner(labels[u], labels[v]) == r:
                if labels[u][r] > labels[v][r]:
                    u, v = v, u
                reach[u, v] = (u, v)
        # Floyd's recurrence also retains an actual path, not just reachability.
        for mid in range(g.n):
            for u in range(g.n):
                for v in range(g.n):
                    if (u, v) not in reach and (u, mid) in reach and (mid, v) in reach:
                        reach[u, v] = reach[u, mid] + reach[mid, v][1:]
        records = tuple((u, v, path) for (u, v), path in sorted(reach.items()))
        layers.append(Graph(g.n, frozenset(tuple(sorted((u, v))) for u, v in reach)))
        witnesses.append(records)
    result = UnionPacket(g, palettes, labels, tuple(layers), tuple(witnesses))
    verify_packet(result)
    return result


def verify_packet(packet):
    need(type(packet) is UnionPacket, "UnionPacket required")
    g, sizes, labels = packet.source, packet.palettes, packet.labels
    validate_frame(g, sizes, labels)
    need(type(packet.layers) is tuple and len(packet.layers) == len(sizes), "layer arity")
    need(type(packet.paths) is tuple and len(packet.paths) == len(sizes), "path arity")
    for r, (layer, records) in enumerate(zip(packet.layers, packet.paths)):
        validate_graph(layer)
        need(layer.n == g.n, "layer source labels")
        need(type(records) is tuple, "canonical path records required")
        seen = set()
        for row in records:
            need(type(row) is tuple and len(row) == 3, "path record arity")
            u, v, path = row
            need(type(u) is int and type(v) is int and 0 <= u < g.n and 0 <= v < g.n,
                 "path endpoint type/range")
            need(type(path) is tuple and len(path) >= 2, "nonempty edge path required")
            need(all(type(x) is int and 0 <= x < g.n for x in path), "path vertex type")
            need(path[0] == u and path[-1] == v and len(set(path)) == len(path), "path endpoints")
            need((u, v) not in seen, "duplicate closure record")
            seen.add((u, v))
            for x, y in zip(path, path[1:]):
                need(tuple(sorted((x, y))) in g.edges, "path must use original edges")
                need(owner(labels[x], labels[y]) == r and labels[x][r] < labels[y][r],
                     "path must keep owner and orientation")
            need(len(path) <= sizes[r], "height budget")
        expected = frozenset(tuple(sorted(e)) for e in seen)
        need(expected == layer.edges, "closure/path incidence mismatch")
        for u, v in g.edges:
            if owner(labels[u], labels[v]) == r:
                if labels[u][r] > labels[v][r]:
                    u, v = v, u
                need((u, v) in seen, "missing original edge")
        for u, v in seen:
            for x, y in seen:
                if v == x:
                    need((u, y) in seen, "closure is not transitive")
        validate_colouring(layer, tuple(label[r] for label in labels), sizes[r])
    return union_graph(packet.layers)


@dataclass(frozen=True)
class Quotient:
    source: Graph
    joint: tuple
    states: tuple
    projection: tuple
    conflict: Graph


def joint_quotient(graphs, colourings, palettes):
    original = union_graph(graphs)
    need(type(colourings) is tuple and type(palettes) is tuple,
         "colouring and palette tuples required")
    need(len(graphs) == len(colourings) == len(palettes), "colouring arity")
    for g, c, k in zip(graphs, colourings, palettes):
        validate_colouring(g, c, k)
    joint = tuple(tuple(c[v] for c in colourings) for v in range(original.n))
    states = tuple(sorted(set(joint)))
    index = {s: i for i, s in enumerate(states)}
    projection = tuple(index[s] for s in joint)
    edges = frozenset(tuple(sorted((projection[u], projection[v]))) for u, v in original.edges)
    conflict = graph(len(states), edges)
    return Quotient(original, joint, states, projection, conflict)


def lift(quotient, state_colours, k):
    need(type(quotient) is Quotient, "Quotient required")
    validate_graph(quotient.source)
    need(type(quotient.joint) is tuple and len(quotient.joint) == quotient.source.n,
         "joint-source arity")
    need(all(type(s) is tuple and all(type(x) is int and x >= 0 for x in s)
             for s in quotient.joint), "exact product states")
    states = tuple(sorted(set(quotient.joint)))
    need(type(quotient.states) is tuple and quotient.states == states, "state labels")
    need(all(type(s) is tuple and all(type(x) is int for x in s) for s in quotient.states),
         "state types")
    need(type(quotient.projection) is tuple and
         all(type(i) is int for i in quotient.projection), "projection types")
    index = {s: i for i, s in enumerate(states)}
    projection = tuple(index[s] for s in quotient.joint)
    need(quotient.projection == projection, "projection mismatch")
    edges = frozenset(tuple(sorted((projection[u], projection[v])))
                      for u, v in quotient.source.edges)
    expected = graph(len(states), edges)
    validate_graph(quotient.conflict)
    need(quotient.conflict == expected, "missing or forged context edge")
    validate_colouring(expected, state_colours, k)
    result = tuple(state_colours[i] for i in projection)
    validate_colouring(quotient.source, result, k)
    return result


def complete_graph(n):
    return graph(n, combinations(range(n), 2))


def cycle(n):
    need(type(n) is int and n >= 3, "cycle size >=3")
    return graph(n, (tuple(sorted((i, (i + 1) % n))) for i in range(n)))


def mycielski(g):
    validate_graph(g)
    n, edges = g.n, set(g.edges)
    for u, v in g.edges:
        edges.add((u, v + n))
        edges.add((v, u + n))
    edges.update((u + n, 2 * n) for u in range(n))
    return graph(2 * n + 1, (tuple(sorted(e)) for e in edges))


def rejects(f):
    try:
        f()
    except (ValueError, TypeError):
        return True
    return False


def main():
    # Complete small source/label universes, including all non-girth-safe graphs.
    compiled = safe = 0
    for n in range(5):
        potential = tuple(combinations(range(n), 2))
        labels_universe = tuple(product(range(2), range(3)))
        # Canonical equality patterns use at most four of six labels; all actual
        # labels are still tested, so ordered-coordinate comparisons are retained.
        for labels in product(labels_universe, repeat=n):
            allowed = tuple(e for e in potential if labels[e[0]] != labels[e[1]])
            for mask in range(1 << len(allowed)):
                g = graph(n, (e for j, e in enumerate(allowed) if mask >> j & 1))
                packet = compile_union(g, (2, 3), labels)
                out = verify_packet(packet)
                check(g.edges <= out.edges, "source edges retained")
                check(all(owner(labels[u], labels[v]) == r
                          for r, h in enumerate(packet.layers) for u, v in h.edges),
                      "unique owner survives closure")
                gg = girth(g)
                if gg is None or gg >= required_girth((2, 3)):
                    check(clique_number(out) <= 3, "mixed palette clique bound")
                    safe += 1
                compiled += 1
    print("universe: palettes=(2,3), all labelled graphs/compatible labels n<=4",
          "packets=", compiled, "girth-safe=", safe)

    # Cyclic, nonforest controls for the sharpened 2-by-3 threshold five.
    c5 = cycle(5)
    cycles = 0
    for labels in product(tuple(product(range(2), range(3))), repeat=5):
        if all(labels[u] != labels[v] for u, v in c5.edges):
            p = compile_union(c5, (2, 3), labels)
            out = verify_packet(p)
            check(clique_number(out) <= 3, "girth-five closure")
            check(all(sum(owner(labels[u], labels[v]) == r for r in range(2)) == 1
                      for u, v in out.edges), "edge owners distinct")
            cycles += 1
    print("C5 complete product-colouring controls:", cycles,
          "required girth (2,3)/(3,2)/(3,3)=",
          tuple(required_girth(k) for k in ((2,3),(3,2),(3,3))))

    # A fully checked tight instance of the paper's k=2,d=2 construction.
    grotzsch = mycielski(c5)
    c4 = colouring(grotzsch, 4)
    check(c4 is not None and colouring(grotzsch, 3) is None, "Grotzsch chromatic four")
    lab = tuple((c // 2, c % 2) for c in c4)
    packet = compile_union(grotzsch, (2,2), lab)
    out = verify_packet(packet)
    quotient = joint_quotient(packet.layers,
                              tuple(tuple(s[r] for s in lab) for r in range(2)), (2,2))
    check(out == grotzsch and girth(grotzsch) == 4 and clique_number(out) == 2,
          "height-two closure and exact source")
    check(quotient.conflict == complete_graph(4), "tight joint-state saturation")
    check(lift(quotient, tuple(range(4)), 4) == c4, "source-owned lift")
    print("tight concrete control: Mycielski(C5), vertices/edges/omega/chi=",
          (grotzsch.n, len(grotzsch.edges), clique_number(grotzsch), chromatic(grotzsch)),
          "joint conflict=K4")

    # Colour-fibre restrictions can cost colours even on an even cycle.
    sequence = (0,1,0,2,1,3,2,3)
    c8 = cycle(8)
    first = graph(8, (e for e in c8.edges if sequence[e[0]]//2 != sequence[e[1]]//2))
    second = graph(8, c8.edges-first.edges)
    q = joint_quotient((first, second),
                       (tuple(x//2 for x in sequence), tuple(x%2 for x in sequence)), (2,2))
    check(chromatic(c8) == 2 and q.conflict == complete_graph(4), "fixed-fibre hostile")
    check(rejects(lambda: lift(q, (0,1,0,1), 2)), "cannot identify conflicting states")
    print("fixed-context hostile: C8 chi=2, chosen two-bit quotient chi=4")

    # Complete pair-state iff, independently tested for every candidate recolouring.
    maps = 0
    for palettes in ((2,2),(2,3),(3,3)):
        states = tuple(product(*(range(k) for k in palettes)))
        n = len(states)
        layers = tuple(graph(n, (e for e in combinations(range(n),2)
                                 if states[e[0]][r] != states[e[1]][r]))
                       for r in range(len(palettes)))
        qq = joint_quotient(layers, tuple(tuple(s[r] for s in states)
                                        for r in range(len(palettes))), palettes)
        check(qq.source == complete_graph(n) and qq.conflict == complete_graph(n),
              "product bound attained on complete multipartite layers")
        check(all(chromatic(g) == k for g, k in zip(layers, palettes)), "component palette exact")
    for k in range(1,5):
        for state_c in product(range(k), repeat=4):
            original = tuple(state_c[i] for i in q.projection)
            expected = all(original[u] != original[v] for u,v in c8.edges)
            check(rejects(lambda: lift(q,state_c,k)) == (not expected), "context-colouring iff")
            maps += 1
    print("edge-union sharp product controls: 2x2=4, 2x3=6, 3x3=9>3+3; maps=", maps)

    # Two bipartite fragments sharing an induced edgeless boundary need not glue.
    even_path = graph(3, ((0,2),(1,2)))
    odd_path = graph(4, ((0,2),(2,3),(1,3)))
    check(colouring(even_path,2) is not None and colouring(odd_path,2) is not None,
          "local colourability")
    same_even = {c[0] == c[1] for c in product(range(2),repeat=3)
                 if all(c[u]!=c[v] for u,v in even_path.edges)}
    same_odd = {c[0] == c[1] for c in product(range(2),repeat=4)
                if all(c[u]!=c[v] for u,v in odd_path.edges)}
    check(same_even == {True} and same_odd == {False}, "overlap parity sidecar")
    check(chromatic(c5) == 3, "glued fragments need three")
    print("overlap hostile: length-2 and length-3 paths share only their ends; union C5")

    # A unit path closed transitively acquires a nonunit edge: no metric transfer.
    line = graph(3, ((0,1),(1,2)))
    lp = compile_union(line, (3,), ((0,),(1,),(2,)))
    check(verify_packet(lp) == complete_graph(3), "transitive completion adds edge02")
    check((2-0)**2 == 4, "added edge has length two in retained realization")
    # Without girth, even the component heights do not bound the union clique.
    bad = compile_union(complete_graph(4),(2,3),((0,0),(0,1),(1,0),(1,1)))
    check(clique_number(verify_packet(bad)) == 4 and girth(bad.source) == 3,
          "girth hypothesis necessary to this proof")

    # Typed packet controls; forged records must not exploit Python aliases.
    invalid = [
        lambda: graph(True, ()),
        lambda: graph(2, ((False,1),)),
        lambda: graph(2, ((0,1),(False,1))),
        lambda: compile_union(line,(True,),((0,),(1,),(0,))),
        lambda: compile_union(line,(3,),((False,),(1,),(2,))),
        lambda: compile_union(line,(3,),((0,),(0,),(2,))),
        lambda: verify_packet(replace(lp,layers=(line,))),
        lambda: verify_packet(replace(lp,paths=(((0,1,(False,1)),),))),
        lambda: verify_packet(replace(lp,paths=(((0,1,(0,1)),(1,2,(1,2))),))),
        lambda: lift(replace(q,projection=(False,)+q.projection[1:]),tuple(range(4)),4),
        lambda: lift(replace(q,conflict=graph(4,())),tuple(range(4)),4),
        lambda: lift(q,(False,1,2,3),4),
        lambda: joint_quotient((first,second),((0,)*8,(0,)*8),(2,2)),
    ]
    for f in invalid:
        check(rejects(f), "malformed packet rejected")
    print("typed malformed packets rejected:",len(invalid))
    print("exact checks:",CHECKS)
    print("scope: mixed-palette closure and fixed-state compression; no plane six-colouring")


if __name__ == '__main__':
    main()
