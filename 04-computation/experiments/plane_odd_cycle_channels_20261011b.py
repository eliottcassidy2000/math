"""Exact ternary-label/binary-lift and odd-cycle obstruction certificates.

Complete controls on graphs of orders 1..5, abstract K6/K7 controls, and
algebraic unit-distance C5/Moser controls. No plane six/seven decision.
"""
from dataclasses import dataclass, replace
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product

import boundary_palette_compression_20261011 as palette
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


def graph(n, edges):
    need(type(n) is int and n >= 0 and type(edges) is frozenset, 'exact graph domain')
    need(all(type(e) is tuple and len(e) == 2 and
             all(type(v) is int for v in e) and 0 <= e[0] < e[1] < n
             for e in edges), 'ordered exact edges')


def edge(u, v):
    return tuple(sorted((u, v)))


def cycle_checked(n, edges, cycle):
    graph(n, edges)
    need(type(cycle) is tuple and len(cycle) >= 3 and len(cycle) % 2 == 1,
         'odd simple cycle required')
    need(all(type(v) is int and 0 <= v < n for v in cycle)
         and len(set(cycle)) == len(cycle), 'distinct actual cycle vertices')
    need(all(edge(cycle[i], cycle[(i+1) % len(cycle)]) in edges
             for i in range(len(cycle))), 'cycle edges are present')


def chordless(edges, cycle):
    cycle = list(cycle)
    while True:
        for i in range(len(cycle)):
            for j in range(i+2, len(cycle)):
                if i == 0 and j == len(cycle)-1:
                    continue
                if edge(cycle[i], cycle[j]) in edges:
                    one = cycle[i:j+1]
                    two = cycle[j:]+cycle[:i+1]
                    cycle = one if len(one) % 2 else two
                    break
            else:
                continue
            break
        else:
            return tuple(cycle)


def binary(n, edges):
    graph(n, edges)  # Type guards precede cache lookup.
    return _binary(n, edges)


@lru_cache(maxsize=None)
def _binary(n, edges):
    adjacent = [[] for _ in range(n)]
    for u, v in sorted(edges):
        adjacent[u].append(v); adjacent[v].append(u)
    bits, parent, components = [-1]*n, [-1]*n, 0
    for root in range(n):
        if bits[root] != -1:
            continue
        components += 1
        bits[root] = 0
        queue = [root]
        for u in queue:
            for v in adjacent[u]:
                if bits[v] == -1:
                    bits[v], parent[v] = 1-bits[u], u
                    queue.append(v)
                elif bits[v] == bits[u]:
                    left, at = [], u
                    while at != -1:
                        left.append(at); at = parent[at]
                    right, at = [], v
                    while at not in left:
                        right.append(at); at = parent[at]
                    cyc = tuple(left[:left.index(at)+1]+list(reversed(right)))
                    cyc = chordless(edges, cyc)
                    cycle_checked(n, edges, cyc)
                    return None, cyc, None
    return tuple(bits), None, components


def mono_edges(n, edges, labels, r):
    graph(n, edges)
    need(type(r) is int and r >= 1 and type(labels) is tuple and len(labels) == n
         and all(type(t) is int and 0 <= t < r for t in labels), 'exact retained labels')
    return frozenset((u, v) for u, v in edges if labels[u] == labels[v])


def lift(n, edges, labels, r):
    bits, cycle, _ = binary(n, mono_edges(n, edges, labels, r))
    if bits is None:
        cycle_checked(n, edges, cycle)
        need(len({labels[v] for v in cycle}) == 1, 'monochromatic obstruction')
        return None, cycle
    colours = tuple(2*t+b for t, b in zip(labels, bits))
    need(all(colours[u] != colours[v] for u, v in edges), 'full proper colour lift')
    return colours, None


def even_character(n, edges):
    """Independent full edge-subset sum; small declared universe only."""
    graph(n, edges)
    return _even_character(n, edges)


@lru_cache(maxsize=None)
def _even_character(n, edges):
    es = tuple(sorted(edges))
    total = 0
    for mask in range(1 << len(es)):
        boundary = 0
        for i, (u, v) in enumerate(es):
            if mask >> i & 1:
                boundary ^= (1 << u) | (1 << v)
        if boundary == 0:
            total += -1 if mask.bit_count() % 2 else 1
    return total


def count_grouped(n, edges, r):
    graph(n, edges)
    need(type(r) is int and r >= 1, 'positive exact channel count')
    total = 0
    for labels in product(range(r), repeat=n):
        mono = mono_edges(n, edges, labels, r)
        bits, _, components = binary(n, mono)
        total += 0 if bits is None else 2**components
    return total


def count_unlabelled_partitions(n, edges, k):
    """Independent positive count via full-colour equality partitions."""
    total = 0
    for assignment in palette.patterns(n, k):
        if all(assignment[u] != assignment[v] for u, v in edges):
            used = len(set(assignment))
            ways = 1
            for j in range(used):
                ways *= k-j
            total += ways
    return total


@dataclass(frozen=True)
class SearchCertificate:
    channels: int
    cycles: tuple
    labels: tuple | None
    colours: tuple | None


def verify(n, edges, cert):
    graph(n, edges)
    need(type(cert) is SearchCertificate and type(cert.channels) is int and cert.channels >= 1,
         'typed certificate and palette')
    need(type(cert.cycles) is tuple, 'retained cycle witnesses')
    for cycle in cert.cycles:
        cycle_checked(n, edges, cycle)
    if cert.labels is not None:
        mono_edges(n, edges, cert.labels, cert.channels)
        need(type(cert.colours) is tuple and len(cert.colours) == n
             and all(type(c) is int and 0 <= c < 2*cert.channels for c in cert.colours),
             'typed full colour witness')
        need(all(c//2 == t for c, t in zip(cert.colours, cert.labels))
             and all(cert.colours[u] != cert.colours[v] for u, v in edges), 'literal colour lift')
    else:
        need(cert.colours is None, 'negative certificate has no invented lift')
        # Exhaustive certificate checking in the declared finite test sizes.
        for assignment in palette.patterns(n, cert.channels):
            need(any(len({assignment[v] for v in c}) == 1 for c in cert.cycles),
                 'cycle constraints must cover every assignment, not just tested ones')
    return cert


def search(n, edges, r):
    graph(n, edges)
    need(type(r) is int and r >= 1, 'exact channel count')
    cycles = []
    while True:
        candidate = next((a for a in palette.patterns(n, r)
                          if all(len({a[v] for v in c}) > 1 for c in cycles)), None)
        if candidate is None:
            return verify(n, edges, SearchCertificate(r, tuple(cycles), None, None))
        colours, bad = lift(n, edges, candidate, r)
        if colours is not None:
            return verify(n, edges, SearchCertificate(r, tuple(cycles), candidate, colours))
        need(frozenset(bad) not in {frozenset(c) for c in cycles}, 'new violated constraint')
        cycles.append(bad)


def unit_pentagon():
    # theta=2sin(pi/5), the unique root in(1,3/2) of t^4-5t^2+5.
    a = algebraic
    field = a.AlgebraicReal((5,0,-5,0,1), F(1), F(3,2))
    height = a.poly((0,F(3,2),0,F(-1,2)))
    points = ((a.ZERO,a.ZERO),(a.ONE,a.ZERO),
              (a.poly((2,0,F(-1,2))),height),
              (a.poly((F(1,2),)),a.poly((0,2,0,F(-1,2)))),
              (a.poly((-1,0,F(1,2))),height))
    return a.unit_graph(points,field)


def rejects(fn):
    try:
        fn()
    except ValueError:
        check(True, 'hostile rejected')
    else:
        raise RuntimeError('hostile accepted')


def main():
    counts, tables = [], 0
    for n in range(1,6):
        es = tuple(combinations(range(n),2)); counts.append(1 << len(es))
        for mask in range(1 << len(es)):
            edges = frozenset(e for i,e in enumerate(es) if mask >> i & 1)
            bits, cycle, components = binary(n, edges)
            value = F(2**n, 2**len(edges))*even_character(n, edges)
            literal = sum(all(a[u] != a[v] for u,v in edges)
                          for a in product(range(2), repeat=n))
            check(value == literal == (0 if bits is None else 2**components),
                  'even-character, literal bit assignments, and BFS agree')
            if cycle is not None:
                cycle_checked(n, edges, cycle)
                check(sum(edge(u,v) in edges for u,v in combinations(cycle,2)) == len(cycle),
                      'learned cycle is chordless')
            for r in (1,2,3):
                exact_count = count_grouped(n,edges,r)
                check(exact_count == count_unlabelled_partitions(n,edges,2*r),
                      'grouped binary channels recover complete labelled-colour count')
                cert = search(n,edges,r)
                check((cert.labels is not None) == (exact_count > 0), 'lazy cycle certificate is complete')
                tables += 1
    check(counts == [1,2,8,64,1024] and sum(counts) == 1099, 'triangular pair-bit census')
    check(10**3+10**2-10 == 1090 and 10**3+10**2-1 == 1099,
          'decimal identity corrected rather than silently affirmed')
    k4 = frozenset(combinations(range(4),2))
    labels = (1,0,0,0)
    star_basis = ((0,1,2),(0,1,3),(0,2,3))
    check(all(len({labels[v] for v in cycle}) > 1 for cycle in star_basis),
          'original-graph fundamental cycles all pass the label constraints')
    colours, obstruction = lift(4,k4,labels,2)
    check(colours is None and frozenset(obstruction) == frozenset((1,2,3)),
          'residual odd cycle survives the original fixed cycle-basis tests')
    pentagon = unit_pentagon()
    expected = frozenset(edge(i,(i+1)%5) for i in range(5))
    check(pentagon == expected, 'algebraic induced regular unit pentagon')
    check(not any(all(edge(u,v) in pentagon for u,v in combinations(triple,2))
                  for triple in combinations(range(5),3)), 'triangle tests see no defect')
    bad = search(5,pentagon,1)
    check(bad.labels is None and len(bad.cycles) == 1 and len(bad.cycles[0]) == 5,
          'unit odd hole defeats triangle-only palette test')
    spindle = palette.spindle()
    check(search(7,spindle.edges,1).labels is None and search(7,spindle.edges,2).labels is not None,
          'exact algebraic Moser two/four-colour controls')
    k6, k7 = frozenset(combinations(range(6),2)),frozenset(combinations(range(7),2))
    check(search(6,k6,3).labels is not None, 'abstract six-colour positive control')
    negative = search(7,k7,3)
    check(negative.labels is None and len(negative.cycles) == 35,
          'abstract K7 obstruction needs every triple; no unit realization asserted')
    rejects(lambda:verify(7,k7,replace(negative,cycles=negative.cycles[:-1])))
    rejects(lambda:search(True,frozenset(),3))
    rejects(lambda:search(3,frozenset(),True))
    rejects(lambda:binary(2,frozenset(((False,1),))))
    rejects(lambda:lift(2,frozenset(((0,1),)),(False,0),1))
    rejects(lambda:cycle_checked(5,pentagon,(0,1,2)))
    check(verify(7,k7,negative) == negative, 'retained complete obstruction positive check')
    print('PROVED elementary: 2r-colours iff r labels with bipartite monochromatic graph.')
    print('Graph census orders1..5:',counts,'sum',sum(counts),'; decimal1000+100-10=1090.')
    print('Complete palette/channel tables:',tables,'; even-subgraph parity checked independently.')
    print('Basis hostile: original K4 star-cycle NAE clauses miss its monochromatic residual triangle.')
    print('Geometry controls: exact induced unit C5 and Moser spindle; triangle-only test fails on C5.')
    print('Abstract K7 negative control:35 learned odd-cycle constraints; no unit-distance claim.')
    print('No plane-six/seven decision, no universal Collatz result, no speed benchmark.')
    print('PASS:',CHECKS,'exact checks')


if __name__ == '__main__':
    main()
