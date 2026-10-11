"""Exact boundary-colouring relations with witness-preserving gluing.

This is a finite prototype and independent completeness checker, not a
six-chromatic plane construction or a general fast colouring algorithm.
"""
from dataclasses import dataclass, replace
from functools import lru_cache
from fractions import Fraction as F
from itertools import combinations, product
import finite_algebraic_witness_20261007 as algebraic

CHECKS = 0


def need(ok, why):
    if not ok:
        raise ValueError(why)


def check(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(why)


@dataclass(frozen=True)
class Graph:
    vertices: tuple
    edges: frozenset


def audit_graph(g):
    need(type(g) is Graph and type(g.vertices) is tuple, 'typed graph and vertex tuple')
    need(all(type(v) is int and v >= 0 for v in g.vertices), 'exact natural vertex labels')
    need(tuple(sorted(set(g.vertices))) == g.vertices, 'sorted distinct vertices')
    need(type(g.edges) is frozenset, 'frozen edge set')
    need(all(type(e) is tuple and len(e) == 2
             and all(type(v) is int and v in g.vertices for v in e)
             and e[0] < e[1] for e in g.edges), 'actual ordered edge pairs')


def canonical(values):
    need(type(values) is tuple and all(type(c) is int and c >= 0 for c in values),
         'exact nonnegative colour labels')
    labels = {}
    return tuple(labels.setdefault(c, len(labels)) for c in values)


def patterns(size, k):
    need(type(size) is int and size >= 0 and type(k) is int and k >= 1,
         'exact size and positive palette')
    def rec(prefix, used):
        if len(prefix) == size:
            yield prefix
        else:
            for c in range(min(k, used+1)):
                yield from rec(prefix+(c,), max(used, c+1))
    return tuple(rec((), 0))


def valid_colouring(g, k, witness):
    audit_graph(g)
    need(type(k) is int and k >= 1, 'positive exact palette')
    need(type(witness) is tuple and len(witness) == len(g.vertices)
         and all(type(c) is int and 0 <= c < k for c in witness), 'full literal colouring')
    a = dict(zip(g.vertices, witness))
    return all(a[u] != a[v] for u, v in g.edges)


@dataclass(frozen=True)
class Table:
    graph: Graph
    boundary: tuple
    palette: int
    entries: tuple


def check_domain(g, boundary, k):
    audit_graph(g)
    need(type(boundary) is tuple and all(type(v) is int and v in g.vertices for v in boundary)
         and len(set(boundary)) == len(boundary), 'distinct named boundary vertices')
    need(type(k) is int and k >= 1, 'positive exact palette')


def table(g, boundary, k):
    check_domain(g, boundary, k)
    return _table(g, boundary, k)


@lru_cache(maxsize=None)
def _table(g, boundary, k):
    # Public validation must precede lookup: bool/int hash and equality alias.
    entries = {}
    for colouring in patterns(len(g.vertices), k):
        if valid_colouring(g, k, colouring):
            a = dict(zip(g.vertices, colouring))
            key = canonical(tuple(a[v] for v in boundary))
            entries.setdefault(key, colouring)
    return Table(g, boundary, k, tuple(sorted(entries.items())))


def keys(t):
    return frozenset(key for key, _ in t.entries)


def audit_table(t):
    """Finite independent-of-certificate completeness check by recomputation."""
    need(type(t) is Table, 'typed boundary table')
    check_domain(t.graph, t.boundary, t.palette)
    need(type(t.entries) is tuple, 'ordered entry tuple')
    found = set()
    for entry in t.entries:
        need(type(entry) is tuple and len(entry) == 2, 'partition/witness pair')
        key, witness = entry
        need(type(key) is tuple and len(key) == len(t.boundary)
             and canonical(key) == key and len(set(key)) <= t.palette, 'canonical partition')
        need(key not in found and valid_colouring(t.graph, t.palette, witness),
             'unique key and proper witness')
        a = dict(zip(t.graph.vertices, witness))
        need(canonical(tuple(a[v] for v in t.boundary)) == key, 'witness has indicated boundary')
        found.add(key)
    need(found == keys(table(t.graph, t.boundary, t.palette)),
         'complete table required; positive witnesses cannot justify omitted states')
    return t


def forget(t, vertex):
    audit_table(t)
    need(type(vertex) is int and vertex in t.boundary, 'retained vertex to forget')
    at = t.boundary.index(vertex)
    boundary = t.boundary[:at]+t.boundary[at+1:]
    entries = {}
    for key, witness in t.entries:
        reduced = canonical(key[:at]+key[at+1:])
        entries.setdefault(reduced, witness)
    return Table(t.graph, boundary, t.palette, tuple(sorted(entries.items())))


def add_boundary_edge(t, u, v):
    audit_table(t)
    need(type(u) is int and type(v) is int and u != v
         and u in t.boundary and v in t.boundary, 'both edge endpoints must be retained')
    i, j = t.boundary.index(u), t.boundary.index(v)
    graph = Graph(t.graph.vertices, t.graph.edges | {tuple(sorted((u,v)))})
    return Table(graph, t.boundary, t.palette,
                 tuple((key, w) for key, w in t.entries if key[i] != key[j]))


def glue(left, right):
    """Shared named boundary, disjoint interiors, and no omitted cross-edge."""
    audit_table(left); audit_table(right)
    need(left.palette == right.palette and left.boundary == right.boundary,
         'same palette and ordered geometric/combinatorial ports')
    need(set(left.graph.vertices) & set(right.graph.vertices) == set(left.boundary),
         'only declared boundary may be shared')
    k, boundary = left.palette, left.boundary
    graph = Graph(tuple(sorted(set(left.graph.vertices) | set(right.graph.vertices))),
                  left.graph.edges | right.graph.edges)
    entries = []
    ldict, rdict = dict(left.entries), dict(right.entries)
    for key in sorted(ldict.keys() & rdict.keys()):
        a = dict(zip(left.graph.vertices, ldict[key]))
        b = dict(zip(right.graph.vertices, rdict[key]))
        permutation = {b[v]: a[v] for v in boundary}
        unused = iter(c for c in range(k) if c not in permutation.values())
        for c in range(k):
            if c not in permutation:
                permutation[c] = next(unused)
        combined = a | {v: permutation[c] for v, c in b.items()}
        witness = tuple(combined[v] for v in graph.vertices)
        need(valid_colouring(graph, k, witness), 'glued witness must be proper')
        entries.append((key, witness))
    return Table(graph, boundary, k, tuple(entries))


def literal_keys(g, boundary, k):
    """Independent full labelled-palette enumeration, without RGS pruning."""
    result = set()
    for word in product(range(k), repeat=len(g.vertices)):
        a = dict(zip(g.vertices, word))
        if all(a[u] != a[v] for u, v in g.edges):
            result.add(canonical(tuple(a[v] for v in boundary)))
    return result


def spindle():
    a = algebraic
    field = a.AlgebraicReal((64,0,-28,0,1), 5, 6)
    s = a.scale(a.poly((0,-20,0,1)), F(1,16))
    t = a.scale(a.poly((0,36,0,-1)), F(1,16))
    st = field.times(s,t)
    points = ((a.ZERO,a.ZERO),(s,a.ZERO),(a.scale(s,F(1,2)),a.poly((F(1,2),))),
              (a.scale(s,F(1,2)),a.poly((F(-1,2),))),
              (a.scale(s,F(5,6)),a.scale(st,F(1,6))),
              (a.scale(a.sub(a.scale(s,5),t),F(1,12)),a.scale(a.add(st,a.poly((5,))),F(1,12))),
              (a.scale(a.add(a.scale(s,5),t),F(1,12)),a.scale(a.sub(st,a.poly((5,))),F(1,12))))
    return Graph(tuple(range(7)),a.unit_graph(points,field))


def rejects(call):
    try:
        call()
    except ValueError:
        check(True, 'hostile correctly rejected')
    else:
        raise RuntimeError('hostile accepted')


def main():
    graph_count = table_count = 0
    for size in range(5):
        vertices = tuple(range(size))
        possible = tuple(combinations(vertices,2))
        for mask in range(1 << len(possible)):
            g = Graph(vertices,frozenset(e for i,e in enumerate(possible) if mask >> i & 1))
            graph_count += 1
            for bmask in range(1 << size):
                boundary = tuple(v for v in vertices if bmask >> v & 1)
                for k in (1,2,3):
                    current = table(g,boundary,k)
                    check(keys(current) == literal_keys(g,boundary,k), 'complete literal-palette agreement')
                    table_count += 1
                    for v in boundary:
                        out = forget(current,v)
                        check(keys(out) == keys(table(g,out.boundary,k)), 'existential boundary forgetting')
                    for u,v in combinations(boundary,2):
                        out = add_boundary_edge(current,u,v)
                        check(keys(out) == keys(table(out.graph,boundary,k)), 'boundary edge filter')

    # 64 pairs with shared ports and private interiors, in each of 3 palettes.
    left_edges = ((0,1),(0,2),(1,2))
    right_edges = ((0,1),(0,3),(1,3))
    for lm,rm,k in product(range(8),range(8),(1,2,3)):
        left = Graph((0,1,2),frozenset(e for i,e in enumerate(left_edges) if lm >> i & 1))
        right = Graph((0,1,3),frozenset(e for i,e in enumerate(right_edges) if rm >> i & 1))
        joined = glue(table(left,(0,1),k),table(right,(0,1),k))
        check(keys(joined) == literal_keys(joined.graph,(0,1),k), 'complete gluing, not merely witness overlap')
        for _,w in joined.entries:
            check(valid_colouring(joined.graph,k,w), 'lifted concrete colouring')

    # Both pieces have exactly two boundary colours; their partitions disagree.
    left = Graph((0,1,2),frozenset(((0,2),(1,2))))
    right = Graph((0,1,2),frozenset(((0,1),(1,2))))
    lt,rt = table(left,(0,1,2),2),table(right,(0,1,2),2)
    check(keys(lt) == {(0,0,1)} and keys(rt) == {(0,1,0)}, 'named partition hostile')
    check({len(set(p)) for p in keys(lt)} == {len(set(p)) for p in keys(rt)} == {2},
          'block counts lose the necessary coordinate')
    check(not glue(lt,rt).entries, 'two paths glue into non-bipartite triangle')

    full = spindle()
    check(len(full.edges) == 11, 'exact algebraic spindle incidence')
    profile = []
    for v in full.vertices:
        boundary = tuple(u for u in full.vertices if tuple(sorted((u,v))) in full.edges)
        rest = Graph(tuple(u for u in full.vertices if u != v),
                     frozenset(e for e in full.edges if v not in e))
        star = Graph(tuple(sorted(boundary+(v,))),
                     frozenset(tuple(sorted((u,v))) for u in boundary))
        outside = table(rest,boundary,3)
        check(bool(outside.entries) and all(len(set(p)) == 3 for p in keys(outside)),
              'outside forces palette saturation on the circle')
        check(not glue(outside,table(star,boundary,3)).entries, 'saturated boundary blocks centre extension')
        local = Graph(boundary,frozenset(e for e in rest.edges if set(e) <= set(boundary)))
        check(bool(table(local,boundary,2).entries), 'local neighbour graph is bipartite')
        four = glue(table(rest,boundary,4),table(star,boundary,4))
        check(bool(four.entries), 'four-colour reconstruction remains possible')
        profile.append((v,len(boundary),len(outside.entries)))
    triangle = Graph((0,1,2),frozenset(combinations(range(3),2)))
    local_edge = Graph((1,2),frozenset(((1,2),)))
    check(keys(table(local_edge,(1,2),2)) == {(0,1)},
          'critical order three: the local edge alone forces its full two-colour palette')
    check(not table(triangle,(),2).entries and not table(full,(),2).entries, 'same empty two-colour interface')
    check(bool(table(triangle,(),3).entries) and not table(full,(),3).entries,
          'changing the palette requires more information')

    counts = {(b,k):sum(len(set(p)) == k for p in patterns(b,k))
              for b,k in ((5,5),(6,5),(7,5),(6,6),(7,6))}
    check(tuple(counts.values()) == (1,15,140,1,21), 'exact palette-saturation partition counts')
    good = table(Graph((0,1),frozenset()),(0,1),2)
    rejects(lambda:audit_table(replace(good,entries=good.entries[:-1])))
    rejects(lambda:add_boundary_edge(forget(good,1),0,1))
    rejects(lambda:table(Graph((False,1),frozenset()),(0,),2))
    rejects(lambda:table(good.graph,(False,1),2))
    rejects(lambda:table(good.graph,(0,0),2))
    rejects(lambda:table(good.graph,(0,1),True))
    rejects(lambda:glue(table(left,(0,),2),table(right,(0,),2)))
    rejects(lambda:glue(good,table(good.graph,(1,0),2)))
    check(audit_table(good) == good, 'positive type/control packet')
    print('PROVED/FINITE-EXACT: exact boundary palette relations and liftable gluing')
    print('Graphs on0..4 vertices:',graph_count,'; all boundary subsets, palettes1..3:',table_count)
    print('Shared two-port gluing controls:192; count-only, local-two-colour, missing-state and forgotten-port hostiles retained.')
    print('Exact Moser circle profile (centre,degree,three-colour boundary partitions):',profile)
    print('Saturated partitions (boundary,palette):',counts)
    print('No six-chromatic graph or universal Collatz coverage is claimed.')
    print('PASS:',CHECKS,'exact checks')


if __name__ == '__main__':
    main()
