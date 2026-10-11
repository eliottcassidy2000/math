"""Exact fixed-ternary boundary parity, gluing, and lifted odd-cycle receipts.

The full six-colour question still quantifies over every ternary assignment.
This is a finite graph instrument, not a new bound on the plane chromatic number.
Run: python -B -X utf8 04-computation/experiments/plane_six_reframes_parity_20261011b.py
"""

from collections import deque
from dataclasses import dataclass, replace
from fractions import Fraction
from itertools import combinations, product


@dataclass(frozen=True)
class Graph:
    vertices: tuple
    edges: tuple


@dataclass(frozen=True)
class Packet:
    graph: Graph
    labels: tuple
    boundary: tuple
    roots: tuple
    parity: tuple
    parents: tuple
    signature: tuple
    odd_cycle: tuple


@dataclass(frozen=True)
class Constraint:
    start: int
    end: int
    parity: int
    path: tuple


@dataclass(frozen=True)
class Gluing:
    graph: Graph
    labels: tuple
    boundary: tuple
    seam_bits: tuple
    colours: tuple
    odd_cycle: tuple


def exact_ints(values, label):
    if type(values) is not tuple or any(type(x) is not int for x in values):
        raise ValueError(label + ': exact integer tuple required')


def audit_graph(graph):
    if type(graph) is not Graph:
        raise ValueError('exact Graph required')
    exact_ints(graph.vertices, 'vertices')
    if graph.vertices != tuple(sorted(set(graph.vertices))) or any(v < 0 for v in graph.vertices):
        raise ValueError('vertices must be distinct sorted nonnegative integers')
    if type(graph.edges) is not tuple:
        raise ValueError('edge tuple required')
    vertices = set(graph.vertices)
    for edge in graph.edges:
        exact_ints(edge, 'edge')
        if len(edge) != 2 or not edge[0] < edge[1] or not set(edge) <= vertices:
            raise ValueError('simple edge with declared endpoints required')
    if graph.edges != tuple(sorted(set(graph.edges))):
        raise ValueError('edges must be distinct and sorted')


def audit_input(graph, labels, boundary):
    audit_graph(graph)
    exact_ints(labels, 'ternary labels')
    if len(labels) != len(graph.vertices) or any(q not in (0, 1, 2) for q in labels):
        raise ValueError('one ternary label per vertex required')
    exact_ints(boundary, 'boundary')
    if len(set(boundary)) != len(boundary) or not set(boundary) <= set(graph.vertices):
        raise ValueError('distinct declared boundary vertices required')


def tree_path(u, v, parents):
    """Unique path in one rooted forest component, retaining actual vertices."""
    left = [u]
    while parents[left[-1]] != -1:
        left.append(parents[left[-1]])
    right = [v]
    while parents[right[-1]] != -1:
        right.append(parents[right[-1]])
    if left[-1] != right[-1]:
        raise ValueError('different forest components')
    i, j = len(left) - 1, len(right) - 1
    while i >= 0 and j >= 0 and left[i] == right[j]:
        i, j = i - 1, j - 1
    return tuple(left[:i + 2] + list(reversed(right[:j + 1])))


def audit_odd_cycle(graph, labels, cycle):
    audit_input(graph, labels, ())
    exact_ints(cycle, 'cycle')
    if len(cycle) < 4 or cycle[0] != cycle[-1] or (len(cycle) - 1) % 2 != 1:
        raise ValueError('closed odd cycle required')
    if len(set(cycle[:-1])) != len(cycle) - 1:
        raise ValueError('simple cycle required')
    q, edges = dict(zip(graph.vertices, labels)), set(graph.edges)
    if any(v not in q for v in cycle) or len({q[v] for v in cycle}) != 1:
        raise ValueError('cycle must be declared and monochromatic')
    if any(tuple(sorted((u, v))) not in edges for u, v in zip(cycle, cycle[1:])):
        raise ValueError('cycle edge absent')
    return True


def compile_boundary(graph, labels, boundary):
    audit_input(graph, labels, boundary)
    q = dict(zip(graph.vertices, labels))
    adj = {v: [] for v in graph.vertices}
    for u, v in graph.edges:
        if q[u] == q[v]:
            adj[u].append(v)
            adj[v].append(u)
    roots, parity, parents, obstruction = {}, {}, {}, ()
    for root in graph.vertices:
        if root in parity:
            continue
        roots[root], parity[root], parents[root] = root, 0, -1
        queue = deque((root,))
        while queue:
            u = queue.popleft()
            for v in adj[u]:
                if v not in parity:
                    roots[v], parity[v], parents[v] = root, parity[u] ^ 1, u
                    queue.append(v)
                elif parity[v] == parity[u] and not obstruction:
                    obstruction = tree_path(u, v, parents) + (u,)
    signature, boundary_roots = [], {}
    if not obstruction:
        for v in boundary:
            representative = boundary_roots.setdefault(roots[v], v)
            signature.append((representative, parity[v] ^ parity[representative]))
    return Packet(graph, labels, boundary,
                  tuple(roots[v] for v in graph.vertices),
                  tuple(parity[v] for v in graph.vertices),
                  tuple(parents[v] for v in graph.vertices),
                  tuple(signature), obstruction)


def audit_packet(packet):
    """Reauthenticate an exported packet; typed guards precede numeric equality."""
    if type(packet) is not Packet:
        raise ValueError('exact Packet required')
    audit_input(packet.graph, packet.labels, packet.boundary)
    for name in ('roots', 'parity', 'parents', 'odd_cycle'):
        exact_ints(getattr(packet, name), name)
    if type(packet.signature) is not tuple:
        raise ValueError('signature tuple required')
    for pair in packet.signature:
        exact_ints(pair, 'signature pair')
        if len(pair) != 2:
            raise ValueError('signature pair arity')
    if packet != compile_boundary(packet.graph, packet.labels, packet.boundary):
        raise ValueError('packet differs from canonical checked compilation')
    return True


def _accepts(packet, bits):
    if packet.odd_cycle:
        return False
    assignment = dict(zip(packet.boundary, bits))
    return all(assignment[v] ^ assignment[root] == p
               for v, (root, p) in zip(packet.boundary, packet.signature))


def accepts(packet, bits):
    audit_packet(packet)
    exact_ints(bits, 'boundary bits')
    if len(bits) != len(packet.boundary) or any(b not in (0, 1) for b in bits):
        raise ValueError('one binary value per port required')
    return _accepts(packet, bits)


def audit_colouring(graph, labels, colours):
    audit_input(graph, labels, ())
    exact_ints(colours, 'colours')
    if len(colours) != len(graph.vertices) or any(c < 0 or c >= 6 for c in colours):
        raise ValueError('six-colour tuple required')
    if tuple(c // 2 for c in colours) != labels:
        raise ValueError('colour does not retain supplied ternary label')
    colour = dict(zip(graph.vertices, colours))
    if any(colour[u] == colour[v] for u, v in graph.edges):
        raise ValueError('improper colouring')
    return True


def lift(packet, bits):
    if not accepts(packet, bits):
        raise ValueError('boundary bits do not extend')
    roots = dict(zip(packet.graph.vertices, packet.roots))
    parity = dict(zip(packet.graph.vertices, packet.parity))
    flips = {}
    for v, bit in zip(packet.boundary, bits):
        flips[roots[v]] = bit ^ parity[v]
    result = tuple(2 * q + (parity[v] ^ flips.get(roots[v], 0))
                   for v, q in zip(packet.graph.vertices, packet.labels))
    audit_colouring(packet.graph, packet.labels, result)
    return result


def constraints(packet):
    audit_packet(packet)
    if packet.odd_cycle:
        raise ValueError('inconsistent fragment has no satisfiable interface')
    parents = dict(zip(packet.graph.vertices, packet.parents))
    return tuple(Constraint(root, v, p, tree_path(root, v, parents))
                 for v, (root, p) in zip(packet.boundary, packet.signature) if root != v)


def simplify_odd_walk(walk):
    """Delete even subloops, or retain an odd subloop, until a simple odd cycle."""
    exact_ints(walk, 'walk')
    if len(walk) < 4 or walk[0] != walk[-1] or (len(walk) - 1) % 2 != 1:
        raise ValueError('closed odd walk required')
    while True:
        first, repeated = {}, None
        for j, v in enumerate(walk[:-1]):
            if v in first:
                repeated = (first[v], j)
                break
            first[v] = j
        if repeated is None:
            return walk
        i, j = repeated
        if (j - i) % 2:
            walk = walk[i:j + 1]
        else:
            walk = walk[:i] + walk[j:]


def glue(left, right):
    """Solve only signed boundary equations, then lift to a colour or real cycle."""
    audit_packet(left)
    audit_packet(right)
    if left.boundary != right.boundary:
        raise ValueError('same ordered boundary required')
    boundary = left.boundary
    if set(left.graph.vertices) & set(right.graph.vertices) != set(boundary):
        raise ValueError('interiors must be disjoint; intersection exactly the boundary')
    ql = dict(zip(left.graph.vertices, left.labels))
    qr = dict(zip(right.graph.vertices, right.labels))
    if any(ql[v] != qr[v] for v in boundary):
        raise ValueError('ternary labels disagree on boundary')
    graph = Graph(tuple(sorted(set(left.graph.vertices) | set(right.graph.vertices))),
                  tuple(sorted(set(left.graph.edges) | set(right.graph.edges))))
    q = ql | qr
    labels = tuple(q[v] for v in graph.vertices)
    for packet in (left, right):
        if packet.odd_cycle:
            audit_odd_cycle(graph, labels, packet.odd_cycle)
            return Gluing(graph, labels, boundary, (), (), packet.odd_cycle)
    equations = constraints(left) + constraints(right)
    adjacency = {v: [] for v in boundary}
    for i, equation in enumerate(equations):
        adjacency[equation.start].append((equation.end, i))
        adjacency[equation.end].append((equation.start, i))
    bits, parents, edge_to_parent = {}, {}, {}
    for root in boundary:
        if root in bits:
            continue
        bits[root], parents[root] = 0, -1
        queue = deque((root,))
        while queue:
            u = queue.popleft()
            for v, edge_id in adjacency[u]:
                equation = equations[edge_id]
                candidate = bits[u] ^ equation.parity
                if v not in bits:
                    bits[v], parents[v], edge_to_parent[v] = candidate, u, edge_id
                    queue.append(v)
                elif bits[v] != candidate:
                    ports = tree_path(u, v, parents)
                    port_edges = []
                    for a, b in zip(ports, ports[1:]):
                        node = a if parents[a] == b else b
                        port_edges.append((a, b, edge_to_parent[node]))
                    port_edges.append((v, u, edge_id))
                    walk = (u,)
                    for a, b, i in port_edges:
                        witness = equations[i].path
                        if witness[0] != a:
                            witness = tuple(reversed(witness))
                        if witness[0] != a or witness[-1] != b:
                            raise ValueError('internal path orientation defect')
                        walk += witness[1:]
                    cycle = simplify_odd_walk(walk)
                    audit_odd_cycle(graph, labels, cycle)
                    return Gluing(graph, labels, boundary, (), (), cycle)
    seam = tuple(bits[v] for v in boundary)
    lc = dict(zip(left.graph.vertices, lift(left, seam)))
    rc = dict(zip(right.graph.vertices, lift(right, seam)))
    if any(lc[v] != rc[v] for v in boundary):
        raise ValueError('internal seam mismatch')
    colour_map = lc | rc
    colours = tuple(colour_map[v] for v in graph.vertices)
    audit_colouring(graph, labels, colours)
    return Gluing(graph, labels, boundary, seam, colours, ())


def ternary_partitions(n):
    """Complete ternary assignments modulo global colour-name permutations."""
    if type(n) is not int or n < 0:
        raise ValueError('nonnegative exact order required')
    if n == 0:
        yield ()
        return
    def visit(prefix):
        if len(prefix) == n:
            yield prefix
        else:
            for q in range(min(2, max(prefix) + 1) + 1):
                yield from visit(prefix + (q,))
    yield from visit((0,))


def all_graphs(vertices):
    pairs = tuple(combinations(vertices, 2))
    for mask in range(1 << len(pairs)):
        yield Graph(vertices, tuple(edge for i, edge in enumerate(pairs) if mask >> i & 1))


def binary_solutions(graph, labels):
    """Independent exhaustive reference, with no compiled forests/signatures."""
    index = {v: i for i, v in enumerate(graph.vertices)}
    return tuple(bits for bits in product((0, 1), repeat=len(graph.vertices))
                 if all(labels[index[u]] != labels[index[v]] or bits[index[u]] != bits[index[v]]
                        for u, v in graph.edges))


def main():
    checks = 0
    def check(test, message):
        nonlocal checks
        checks += 1
        if not test:
            raise RuntimeError(message)
    def reject(call, message):
        nonlocal checks
        checks += 1
        try:
            call()
        except ValueError:
            return
        raise RuntimeError('hostile accepted: ' + message)

    graph_count = assignment_count = positive = negative = signature_count = lifted = 0
    for n in range(1, 6):
        vertices = tuple(range(n))
        boundaries = tuple(dict.fromkeys(((), vertices[:min(3, n)], vertices)))
        for graph in all_graphs(vertices):
            graph_count += 1
            for q in ternary_partitions(n):
                assignment_count += 1
                solutions = binary_solutions(graph, q)
                positive += bool(solutions)
                negative += not solutions
                for boundary in boundaries:
                    signature_count += 1
                    packet = compile_boundary(graph, q, boundary)
                    check(audit_packet(packet), 'canonical compilation')
                    check(bool(solutions) == (not packet.odd_cycle), 'bipartite iff binary lift')
                    expected = {tuple(bits[v] for v in boundary) for bits in solutions}
                    actual = {bits for bits in product((0, 1), repeat=len(boundary)) if accepts(packet, bits)}
                    check(expected == actual, 'complete boundary relation')
                    if packet.odd_cycle:
                        check(audit_odd_cycle(graph, q, packet.odd_cycle), 'lifted odd cycle')
                    else:
                        for bits in actual:
                            colour = lift(packet, bits)
                            check(audit_colouring(graph, q, colour), 'complete actual lift')
                            check(tuple(colour[v] % 2 for v in boundary) == bits, 'boundary preserved')
                            lifted += 1

    glue_count = glue_yes = glue_no = 0
    for left_graph in all_graphs((0, 1, 2)):
        for right_graph in all_graphs((0, 1, 3)):
            for a, b, c, d in product(range(3), repeat=4):
                left = compile_boundary(left_graph, (a, b, c), (0, 1))
                right = compile_boundary(right_graph, (a, b, d), (0, 1))
                result = glue(left, right)
                direct = binary_solutions(result.graph, result.labels)
                check(bool(direct) == (not result.odd_cycle), 'glue iff direct solution')
                if result.odd_cycle:
                    check(audit_odd_cycle(result.graph, result.labels, result.odd_cycle), 'seam path odd-cycle lift')
                    glue_no += 1
                else:
                    check(audit_colouring(result.graph, result.labels, result.colours), 'glue full colour lift')
                    glue_yes += 1
                glue_count += 1

    even = compile_boundary(Graph((0, 1, 2), ((0, 2), (1, 2))), (0, 0, 0), (0, 1))
    odd = compile_boundary(Graph((0, 1), ((0, 1),)), (0, 0), (0, 1))
    triangle = glue(even, odd)
    check(even.signature == ((0, 0), (0, 0)), 'even path equals')
    check(odd.signature == ((0, 0), (0, 1)), 'odd path differs')
    check(bool(triangle.odd_cycle), 'connectivity alone loses odd seam')
    check(not compile_boundary(triangle.graph, (0, 0, 1), (0, 1)).odd_cycle,
          'fixed-ternary obstruction is not a six-colour obstruction')
    # Coordinates x+y*sqrt(3)*i make all three squared distances exact rational.
    points = ((Fraction(0), Fraction(0)), (Fraction(1), Fraction(0)),
              (Fraction(1, 2), Fraction(1, 2)))
    for u, v in combinations(range(3), 2):
        check((points[u][0] - points[v][0]) ** 2 + 3 * (points[u][1] - points[v][1]) ** 2 == 1,
              'unit triangle edge')

    # Signed seam cycles can use several constraints and repeated real vertices.
    l4 = compile_boundary(Graph((0, 1, 2, 3, 4), ((0, 4), (1, 4), (2, 4), (3, 4))), (0,) * 5, (0, 1, 2, 3))
    r4 = compile_boundary(Graph((0, 1, 2, 3), ((0, 1), (1, 2), (2, 3))), (0,) * 4, (0, 1, 2, 3))
    multi = glue(l4, r4)
    check(audit_odd_cycle(multi.graph, multi.labels, multi.odd_cycle), 'path-provenance collision')
    walk = (0, 1, 2, 1, 3, 0)
    check(simplify_odd_walk(walk) == (0, 1, 3, 0), 'delete even detour')
    check(simplify_odd_walk((0, 1, 2, 3, 1, 0)) == (1, 2, 3, 1), 'retain odd subloop')

    for length in (3, 5, 7, 9, 11):
        vertices = tuple(range(length))
        edges = tuple(sorted({tuple(sorted((i, (i + 1) % length))) for i in vertices}))
        cycle_graph = Graph(vertices, edges)
        packet = compile_boundary(cycle_graph, (0,) * length, ())
        check(len(packet.odd_cycle) - 1 == length, 'no bounded odd-cycle-length completeness')
        for omitted in range(length):
            path_graph = Graph(vertices, tuple(e for i, e in enumerate(edges) if i != omitted))
            check(not compile_boundary(path_graph, (0,) * length, ()).odd_cycle, 'every proper cycle-edge subgraph bipartite')

    pairs = tuple(product(range(3), range(2)))
    categorical = tuple((i, j) for i, j in combinations(range(6), 2)
                        if pairs[i][0] != pairs[j][0] and pairs[i][1] != pairs[j][1])
    lexicographic = tuple((i, j) for i, j in combinations(range(6), 2)
                          if pairs[i][0] != pairs[j][0] or pairs[i][1] != pairs[j][1])
    check(len(categorical) == 6 and all(sum(v in edge for edge in categorical) == 2 for v in range(6)),
          'categorical product is six-cycle')
    check(len(lexicographic) == 15, 'lexicographic product is K6')
    check(bool(binary_solutions(Graph(tuple(range(6)), categorical), (0,) * 6)), 'categorical product bipartite')
    check(not any(all(tuple(sorted((image[u], image[v]))) in categorical for u, v in ((0, 1), (1, 2), (0, 2)))
                  for image in product(range(6), repeat=3)), 'triangle cannot map to categorical product')
    check(graph_count == sum(2 ** (n * (n - 1) // 2) for n in range(1, 6)) == 1099, 'labelled graph count')
    check(10 ** 3 + 10 ** 2 - 10 == 1090 != 1099, 'literal decimal hostile')
    check(10 ** 3 + 10 ** 2 - 1 == 1099, 'nearby true decimal identity')

    hostile_count = 0
    hostiles = (
        lambda: compile_boundary(Graph((False, 1), ((0, 1),)), (0, 0), ()),
        lambda: compile_boundary(Graph((0, 1), ((0, 1.0),)), (0, 0), ()),
        lambda: compile_boundary(Graph((0, 1), ((0, 0),)), (0, 0), ()),
        lambda: compile_boundary(Graph((0, 1), ((0, 1), (0, 1))), (0, 0), ()),
        lambda: compile_boundary(Graph((0, 1), ((0, 1),)), (False, 0), ()),
        lambda: compile_boundary(Graph((0, 1), ((0, 1),)), (0, 3), ()),
        lambda: compile_boundary(Graph((0, 1), ((0, 1),)), (0, 0), (0, 0)),
        lambda: accepts(odd, (False, 1)),
        lambda: lift(odd, (0, 0)),
        lambda: audit_packet(replace(odd, parity=(False, 1))),
        lambda: audit_packet(replace(odd, signature=((0, 0), (0, 0)))),
        lambda: audit_packet(replace(odd, parents=(-1, 1))),
        lambda: glue(even, compile_boundary(odd.graph, (0, 0), (1, 0))),
        lambda: glue(even, compile_boundary(odd.graph, (0, 1), (0, 1))),
        lambda: glue(even, even),
        lambda: audit_odd_cycle(triangle.graph, (0, 0, 0), (0, 2, 0)),
        lambda: simplify_odd_walk((0, 1, 0)),
        lambda: tuple(ternary_partitions(True)),
    )
    for call in hostiles:
        reject(call, 'typed or forged boundary')
        hostile_count += 1
    empty = compile_boundary(Graph((), ()), (), ())
    check(lift(empty, ()) == (), 'empty graph and boundary')
    check(glue(empty, empty).colours == (), 'empty gluing')
    check(not accepts(compile_boundary(triangle.graph, (0, 0, 0), ()), ()),
          'empty boundary does not erase internal obstruction')

    print('Plane six reframes: signed boundary parity and exact odd-cycle lifts')
    print(f'Labelled simple graphs, orders 1..5: {graph_count}')
    print(f'Ternary assignments modulo global S3: {assignment_count}; extendible={positive}, obstructed={negative}')
    print(f'Complete boundary relations: {signature_count}; authenticated colouring lifts: {lifted}')
    print(f'Two-port gluings: {glue_count}; compatible={glue_yes}, obstructed={glue_no}')
    print('Unit-triangle hostile: same boundary connectivity, opposite parity; lifted odd cycle', triangle.odd_cycle)
    print('Graph products: K3[K2]=K6; categorical K3 x K2=C6; a triangle separates them')
    print('Arithmetic: 10^3+10^2-10=1090; sum(n=1..5) 2^binom(n,2)=1099')
    print(f'Typed/provenance hostiles rejected: {hostile_count}')
    print(f'Checks: {checks}')
    print('Scope: complete fixed-label gluing; no universal six-colouring of unit-distance graphs proved')


if __name__ == '__main__':
    main()
