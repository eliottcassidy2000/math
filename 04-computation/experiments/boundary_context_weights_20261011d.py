"""Exact context-relative colour weights; no plane chromatic decision.

Universe and controls are printed by main. Public inputs are checked before
enumeration; gates stay active under python -O. Uses the inherited graph and
exact Moser spindle constructor, not a numerical geometric approximation.
"""
from collections import Counter, deque
from itertools import combinations, product
from math import factorial
import boundary_palette_compression_20261011 as p

CHECKS = 0


def check(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)


def falling(k, j):
    return factorial(k) // factorial(k-j)


def conditional_count(graph, boundary, k, pattern):
    p.check_domain(graph, boundary, k)
    p.need(type(pattern) is tuple and len(pattern) == len(boundary)
           and p.canonical(pattern) == pattern and len(set(pattern)) <= k,
           'canonical boundary partition in this palette')
    assignment = dict(zip(boundary, pattern))
    neighbours = {v: set() for v in graph.vertices}
    for u, v in graph.edges:
        neighbours[u].add(v)
        neighbours[v].add(u)
        if u in assignment and v in assignment and assignment[u] == assignment[v]:
            return 0
    interior = tuple(v for v in graph.vertices if v not in assignment)

    def visit(i):
        if i == len(interior):
            return 1
        v = interior[i]
        used = {assignment[u] for u in neighbours[v] if u in assignment}
        result = 0
        for colour in range(k):
            if colour not in used:
                assignment[v] = colour
                result += visit(i+1)
        assignment.pop(v, None)
        return result
    return visit(0)


def weights(graph, boundary, k):
    p.check_domain(graph, boundary, k)
    answer = {}
    for pattern in p.patterns(len(boundary), k):
        count = conditional_count(graph, boundary, k, pattern)
        if count:
            answer[pattern] = count
    return answer


def literal_counts(graph, boundary, k):
    """Independent enumeration: count all colour names, then quotient."""
    p.check_domain(graph, boundary, k)
    result = Counter()
    for cs in product(range(k), repeat=len(graph.vertices)):
        assignment = dict(zip(graph.vertices, cs))
        if all(assignment[u] != assignment[v] for u, v in graph.edges):
            result[p.canonical(tuple(assignment[v] for v in boundary))] += 1
    return dict(result)


def count_from_weights(w, k):
    return sum(falling(k, len(set(pi))) * value for pi, value in w.items())


def projected_weights(graph, boundary, k, vertex):
    p.check_domain(graph, boundary, k)
    p.need(type(vertex) is int and vertex in boundary, 'actual retained port')
    at = boundary.index(vertex)
    reduced = boundary[:at] + boundary[at+1:]
    old = weights(graph, boundary, k)
    result = {}
    for pi in p.patterns(len(reduced), k):
        j = len(set(pi))
        total = 0
        for colour in range(j+1):
            multiplicity = 1 if colour < j else k-j
            full = p.canonical(pi[:at] + (colour,) + pi[at:])
            total += multiplicity * old.get(full, 0)
        if total:
            result[pi] = total
    return reduced, result


def singleton_context(pattern, k):
    p.need(type(k) is int and k >= 1 and type(pattern) is tuple
           and p.canonical(pattern) == pattern and len(set(pattern)) <= k,
           'canonical pattern and exact positive palette')
    b = len(pattern)
    anchors = tuple(range(b, b+k))
    edges = set(combinations(anchors, 2))
    for port, colour in enumerate(pattern):
        edges.update((port, b+c) for c in range(k) if c != colour)
    return p.Graph(tuple(range(b+k)), frozenset(edges)), tuple(range(b))


def graph_from_mask(n, mask):
    pairs = tuple(combinations(range(n), 2))
    return p.Graph(tuple(range(n)), frozenset(e for i,e in enumerate(pairs) if mask >> i & 1))


def add_hidden(graph, r, leaves):
    start = max(graph.vertices, default=-1)+1
    hidden = tuple(range(start, start+r))
    edges = graph.edges | (frozenset((0,v) for v in hidden) if leaves else frozenset())
    return p.Graph(graph.vertices+hidden, edges)


def kempe_separate(graph, ports, k, colours):
    """Return a separated colouring, or paths blocking only single swaps.

    A path result is NOT a forced-equality or uncolourability certificate.
    This combinatorial routine does not authenticate unit coordinates.
    """
    p.check_domain(graph, ports, k)
    p.need(len(ports) == 2 and tuple(sorted(ports)) not in graph.edges,
           'two nonadjacent retained ports')
    p.need(p.valid_colouring(graph,k,colours), 'proper source colouring')
    assignment = dict(zip(graph.vertices,colours))
    u,v = ports
    if assignment[u] != assignment[v]:
        return colours, ()
    a = assignment[u]
    adjacency = {w:set() for w in graph.vertices}
    for x,y in graph.edges:
        adjacency[x].add(y)
        adjacency[y].add(x)
    paths = []
    for b in range(k):
        if b == a:
            continue
        parent, queue = {u:None}, deque([u])
        while queue:
            x = queue.popleft()
            for y in sorted(adjacency[x]):
                if assignment[y] in (a,b) and y not in parent:
                    parent[y] = x
                    queue.append(y)
        if v not in parent:
            result = tuple(a+b-assignment[w] if w in parent else assignment[w]
                           for w in graph.vertices)
            return result, ()
        path, at = [], v
        while at is not None:
            path.append(at)
            at = parent[at]
        paths.append(tuple(reversed(path)))
    return None, tuple(paths)


def rejects(fn):
    try:
        fn()
    except (ValueError, TypeError):
        check(True, 'hostile rejected')
    else:
        check(False, 'hostile accepted')


def main():
    tables = 0
    for n in range(5):
        for mask in range(1 << (n*(n-1)//2)):
            g = graph_from_mask(n, mask)
            for k in range(1,5):
                for bmask in range(1 << n):
                    boundary = tuple(v for v in range(n) if bmask >> v & 1)
                    w = weights(g, boundary, k)
                    literal = literal_counts(g, boundary, k)
                    check(literal == {pi: falling(k,len(set(pi)))*c for pi,c in w.items()},
                          'conditional weights versus independently named colours')
                    check(set(w) == set(p.keys(p.table(g,boundary,k))), 'inherited Boolean support')
                    for v in boundary:
                        reduced, projected = projected_weights(g,boundary,k,v)
                        check(projected == weights(g,reduced,k), 'weighted port projection')
                    tables += 1

    selectors = 0
    for k in range(1,7):
        for b in range(5):
            for pi in p.patterns(b,k):
                g,boundary = singleton_context(pi,k)
                check(weights(g,boundary,k) == {pi:factorial(k-len(set(pi)))},
                      'abstract singleton selector and exact stabilizer weight')
                selectors += 1

    # Every pair of labelled two-port, one-interior fragments; interiors disjoint.
    glued_pairs = 0
    for lm, rm, k in product(range(8), range(8), range(1,5)):
        left = graph_from_mask(3,lm)
        raw = graph_from_mask(3,rm)
        rename = {0:0,1:1,2:3}
        right = p.Graph((0,1,3),frozenset(tuple(sorted((rename[u],rename[v]))) for u,v in raw.edges))
        joined = p.Graph((0,1,2,3),left.edges | right.edges)
        lw,rw = weights(left,(0,1),k),weights(right,(0,1),k)
        expected = {pi:c*rw[pi] for pi,c in lw.items() if pi in rw}
        check(weights(joined,(0,1),k) == expected, 'conditional fibre product')
        glued_pairs += 1

    # All 16 contexts on four vertices containing the marked P3 core.
    # The fourth context vertex may touch any of the three core ports.
    core = frozenset({(0,1),(0,2)})
    optional = ((1,2),(0,3),(1,3),(2,3))
    hidden_tests = 0
    for mask,k,r in product(range(16),range(2,7),range(4)):
        context = p.Graph((0,1,2,3),core | frozenset(e for i,e in enumerate(optional) if mask >> i & 1))
        star,isolates = add_hidden(context,r,True),add_hidden(context,r,False)
        sw,iw = weights(star,(0,1,2),k),weights(isolates,(0,1,2),k)
        check(set(sw) == set(iw), 'same future-observable boundary support')
        check({pi:k**r*c for pi,c in sw.items()} == {pi:(k-1)**r*c for pi,c in iw.items()},
              'weighted relation remains valid for empty supports')
        hidden_tests += 1

    # A forgotten leaf cannot subsequently be used as an attachment port.
    path = p.Graph((0,1,2,3),core)
    star = p.Graph((0,1,2,3),core | {(0,3)})
    check(weights(path,(0,1,2),3).keys() == weights(star,(0,1,2),3).keys(), 'hidden leaf support')
    check(weights(path,(0,1,2,3),3).keys() != weights(star,(0,1,2,3),3).keys(), 'exposed leaf distinguishes')

    kempe_tests = 0
    for n in range(2,5):
        for mask in range(1 << (n*(n-1)//2)):
            g = graph_from_mask(n,mask)
            nonedges = tuple(e for e in combinations(range(n),2) if e not in g.edges)
            for k in range(2,5):
                for cs in product(range(k),repeat=n):
                    if not p.valid_colouring(g,k,cs):
                        continue
                    for ports in nonedges:
                        repaired, paths = kempe_separate(g,ports,k,cs)
                        if repaired is not None:
                            check(p.valid_colouring(g,k,repaired) and repaired[ports[0]] != repaired[ports[1]],
                                  'returned separated colouring is literal')
                        else:
                            a = cs[ports[0]]
                            check(cs[ports[1]] == a and len(paths) == k-1, 'all alternative channels recorded')
                            channel_colours = set()
                            for route in paths:
                                check(route[0] == ports[0] and route[-1] == ports[1]
                                      and len(route) % 2 == 1 and len(route) == len(set(route))
                                      and all(tuple(sorted(e)) in g.edges for e in zip(route,route[1:])),
                                      'literal simple even Kempe path')
                                used = {cs[w] for w in route}
                                check(len(used) == 2 and a in used, 'one alternating channel')
                                channel_colours.update(used-{a})
                            check(channel_colours == set(range(k))-{a}, 'all alternative colours')
                        kempe_tests += 1
    # Unit rhombus with unit-separated omitted diagonal; a blocked single
    # swap does not imply forced equality in all colourings.
    rhombus = p.Graph((0,1,2,3),frozenset({(0,2),(0,3),(1,2),(1,3)}))
    blocked, routes = kempe_separate(rhombus,(0,1),3,(0,0,1,2))
    check(blocked is None and len(routes) == 2, 'both single swaps blocked')
    check(p.valid_colouring(rhombus,3,(0,1,2,2)), 'a different colouring still separates')
    for k in range(4,8):
        pairs = frozenset(e for e in combinations(range(k+1),2) if e != (0,1))
        abstract = p.Graph(tuple(range(k+1)),pairs)
        repaired,routes = kempe_separate(abstract,(0,1),k,(0,0)+tuple(range(1,k)))
        check(repaired is None and len(routes) == k-1 and all(len(route) == 3 for route in routes),
              'abstract complete-minus-edge violates geometric common-neighbour cap')
    # Unit triangle: deleting one edge forces equal endpoint colours at k=2.
    triangle_path = p.Graph((0,1,2),frozenset({(0,2),(1,2)}))
    check(set(weights(triangle_path,(0,1),2)) == {(0,0)}, 'two-colour equality tester')
    # Inherited exact coordinates prove that every spindle edge is unit.
    spindle = p.spindle()
    check(count_from_weights(weights(spindle,(),3),3) == 0, 'spindle rejects three colours')
    for edge in sorted(spindle.edges):
        deleted = p.Graph(spindle.vertices,spindle.edges-{edge})
        check(set(weights(deleted,edge,3)) == {(0,0)}, 'unit-separation equality after spindle edge deletion')
    # Six-colour abstract selector includes K6: not a planar unit gadget.
    selector,boundary = singleton_context((0,0),6)
    check(len(selector.vertices) == 8 and set(weights(selector,boundary,6)) == {(0,0)},
          'abstract equality requires its geometric sidecar')
    one = p.Graph((0,),frozenset())
    check(weights(one,(),6) == {():6} and weights(one,(0,),6) == {(0,):1}, 'six names versus one orbit')
    rejects(lambda: weights(one,(False,),6))
    rejects(lambda: weights(one,(0,),True))
    rejects(lambda: singleton_context((0,2),3))
    rejects(lambda: conditional_count(one,(0,),2,(False,)))
    rejects(lambda: projected_weights(one,(0,),2,False))
    rejects(lambda: kempe_separate(path,(0,False),3,(0,1,1,2)))
    rejects(lambda: kempe_separate(path,(0,1),3,(0,1,1,2)))
    print('Boundary contexts, weighted lifts, and unit-separation testers')
    print('Complete weighted tables:',tables,'(all labelled graphs n=0..4, all boundaries, k=1..4)')
    print('Abstract singleton selectors:',selectors,'(boundary size0..4, k=1..6)')
    print('Two-port fibre products:',glued_pairs)
    print('Hidden-leaf versus isolate contexts:',hidden_tests,'(16 contexts, k=2..6, r=0..3)')
    print('Single-Kempe-swap separator inputs:',kempe_tests,'(all proper k=2..4 colourings, n=2..4, every nonedge)')
    print('Geometric positive controls: unit triangle and every one of 11 exact spindle edge deletions.')
    print('PROVED scope: Boolean context equivalence, weighted gluing/projection, abstract minimality.')
    print('OPEN: a six-colour equality tester whose distinct ports have unit separation.')
    print('Always-active exact checks:',CHECKS)


if __name__ == '__main__':
    main()
