"""Exact colour-orbit quanta and necessary guards for unit equality testers.

Small finite checks authenticate the count identities, not a plane colouring
decision. The 27-vertex guard separately uses two explicitly cited theorems.
"""
from collections import Counter
from fractions import Fraction
from itertools import combinations
from math import prod
import boundary_palette_compression_20261011 as p
import boundary_context_weights_20261011d as w

CHECKS = 0


def check(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(why)


def falling(a, b):
    p.need(type(a) is int and type(b) is int and 0 <= b <= a,
           'falling factorial exact nonnegative domain')
    return prod(range(a-b+1, a+1))


def profile(g, ports):
    """Complete partitions, independent of a chosen labelled palette."""
    p.check_domain(g, ports, 1)
    p.need(len(ports) == 2 and tuple(sorted(ports)) not in g.edges,
           'two nonadjacent ports')
    pos = {v: i for i, v in enumerate(g.vertices)}
    eq, different = Counter(), Counter()
    for cs in p.patterns(len(g.vertices), len(g.vertices)):
        if all(cs[pos[a]] != cs[pos[b]] for a, b in g.edges):
            target = eq if cs[pos[ports[0]]] == cs[pos[ports[1]]] else different
            target[len(set(cs))] += 1
    return dict(eq), dict(different)


def labelled_counts(eq, different, k):
    p.need(type(k) is int and k >= 2, 'at least two exact colours')
    e = sum(falling(k-1, j-1)*a for j, a in eq.items() if j <= k)
    d = sum(falling(k-2, j-2)*a for j, a in different.items() if j <= k)
    return e, d


def quantum(k, s):
    p.need(type(k) is int and type(s) is int and 2 <= s <= k,
           'independent chromatic lower premise in palette')
    return falling(k-1, s-1), falling(k-2, s-2)


def authenticate_small_receipt(g, ports, k, s, upper, witness):
    """Finite prototype: recompute the premise AND the claimed upper bound.

    For large candidates these two enumerations need real proof certificates;
    a bare integer claimed as an upper bound is not trusted evidence.
    This routine authenticates no Euclidean embedding.
    """
    ge, gd = quantum(k, s)
    p.need(type(upper) is int and 0 <= upper < gd, 'upper below orbit quantum')
    p.need(p.valid_colouring(g, k, witness), 'actual proper witness')
    eq, different = profile(g, ports)
    chi = min(set(eq) | set(different))
    p.need(chi >= s, 'authenticated chromatic lower premise')
    e, d = labelled_counts(eq, different, k)
    a = dict(zip(g.vertices, witness))
    p.need(a[ports[0]] == a[ports[1]], 'positive equal-port witness')
    p.need(d <= upper, 'authenticated absolute separating count bound')
    p.need(e % ge == 0 and d % gd == 0, 'orbit divisibility')
    return d == 0


def rejects(fn):
    try:
        fn()
    except ValueError:
        return True
    return False


def neighbours(g):
    out = {v: set() for v in g.vertices}
    for a, b in g.edges:
        out[a].add(b)
        out[b].add(a)
    return out


def shrikhande():
    pts = tuple((a, b) for a in range(4) for b in range(4))
    steps = {(1, 0), (3, 0), (0, 1), (0, 3), (1, 1), (3, 3)}
    edges = frozenset((i, j) for i, j in combinations(range(16), 2)
                      if ((pts[j][0]-pts[i][0]) % 4,
                          (pts[j][1]-pts[i][1]) % 4) in steps)
    return p.Graph(tuple(range(16)), edges)


def fixed_universe_trace(target, ports, k, witness):
    """Greedily delete bad assignments by adding edges on fixed vertices.

    The actual target edges and a surviving equal-port colouring are retained.
    Every factor is measured on the preceding survivor, never assumed1-1/k.
    """
    p.need(p.valid_colouring(target, k, witness), 'retained target witness')
    a = dict(zip(target.vertices, witness))
    p.need(a[ports[0]] == a[ports[1]], 'equal ports in retained witness')
    current = p.Graph(target.vertices, frozenset())
    old = w.conditional_count(current, ports, k, (0, 1))
    initial, factor, trace = old, Fraction(1), []
    remaining = set(target.edges)
    while remaining and old:
        choices = []
        for edge in remaining:
            trial = p.Graph(current.vertices, current.edges | {edge})
            choices.append((w.conditional_count(trial, ports, k, (0, 1)), edge))
        new, edge = min(choices)
        current = p.Graph(current.vertices, current.edges | {edge})
        remaining.remove(edge)
        rho = Fraction(new, old)
        factor *= rho
        check(initial*factor == new, 'conditional counting chain on one universe')
        check(p.valid_colouring(current, k, witness), 'specific equal witness survives')
        trace.append((edge, new, str(rho)))
        old = new
    return initial, old, trace


def main():
    print('Universe: every labelled simple graph on 2..5 vertices; every nonedge port pair.')
    print('Partition counts in palettes2..7; independent literal fixed-port enumeration on2..4.')
    graph_cases = port_cases = 0
    for n in range(2, 6):
        for mask in range(1 << (n*(n-1)//2)):
            g = w.graph_from_mask(n, mask)
            graph_cases += 1
            for ports in combinations(range(n), 2):
                if ports in g.edges:
                    continue
                port_cases += 1
                eq, different = profile(g, ports)
                chi = min(set(eq) | set(different))
                for k in range(2, 8):
                    e, d = labelled_counts(eq, different, k)
                    for s in range(2, min(chi, k)+1):
                        ge, gd = quantum(k, s)
                        check(e % ge == 0 and d % gd == 0, 'all independent chromatic premises')
                    if n <= 4:
                        check(e == w.conditional_count(g, ports, k, (0, 0)), 'equal literal count')
                        check(d == w.conditional_count(g, ports, k, (0, 1)), 'different literal count')
    print(f'Enumerated {graph_cases} graphs and {port_cases} marked nonedge cases.')
    print('k6 equal coefficients:', [falling(5, j-1) for j in range(1, 7)])
    print('k6 different coefficients:', [falling(4, j-2) for j in range(2, 7)])
    check(quantum(6, 5) == (120, 24), 'six-colour five-chromatic orbit floor')

    # Abstract positive controls. The diamond apex distance is sqrt(3), not1;
    # K7 minus an edge has no planar unit-distance embedding.
    diamond = p.Graph(tuple(range(4)), frozenset(combinations(range(4), 2))-{(0, 1)})
    complete_minus = p.Graph(tuple(range(7)), frozenset(combinations(range(7), 2))-{(0, 1)})
    check(authenticate_small_receipt(diamond, (0, 1), 3, 3, 0, (0, 0, 1, 2)), 'diamond receipt')
    check(authenticate_small_receipt(complete_minus, (0, 1), 6, 5, 23,
                                      (0, 0, 1, 2, 3, 4, 5)), 'abstract six-palette receipt')
    print('Abstract K7-minus-edge fixed-port counts:', labelled_counts(*profile(complete_minus, (0, 1)), 6))
    rhombus = p.Graph(tuple(range(4)), diamond.edges-{(2, 3)})
    check(rejects(lambda: authenticate_small_receipt(rhombus, (0, 1), 3, 3, 0,
                                                     (0, 0, 1, 1))), 'forged chromatic premise')
    check(rejects(lambda: authenticate_small_receipt(diamond, (0, 1), 6, 3, 3,
                                                     (0, 0, 1, 2))), 'forged numerical upper bound')
    check(rejects(lambda: quantum(True, 2)), 'bool palette alias')
    check(rejects(lambda: profile(diamond, (False, 1))), 'bool vertex alias')
    check(rejects(lambda: authenticate_small_receipt(diamond, (0, 1), 3, 3, 0,
                                                     (0, 0, 0, 2))), 'improper claimed witness')

    spindle = p.spindle()
    ports = min(spindle.edges)
    fragment = p.Graph(spindle.vertices, spindle.edges-{ports})
    witness = dict(p.table(fragment, ports, 3).entries)[(0, 0)]
    for k in (3, 6):
        initial, last, trace = fixed_universe_trace(fragment, ports, k, witness)
        check(initial == k**5, 'fixed seven-vertex ambient universe')
        check(last == w.conditional_count(fragment, ports, k, (0, 1)), 'complete residual count')
        check((last == 0) == (k == 3), 'positive3-palette and hostile6-palette control')
        print(f'Fixed-universe Moser-minus-{ports}, palette{k}: {initial}->{last}; (edge,residual,factor)={trace}')

    # Pure counting is saturated by an abstract object without unit geometry.
    g = shrikhande()
    adj = neighbours(g)
    check(len(g.edges) == 48 and all(len(adj[v]) == 6 for v in g.vertices), 'Shrikhande degrees')
    check(all(len(adj[a] & adj[b]) == 2 for a, b in combinations(g.vertices, 2)), 'all codegrees2')
    for v in g.vertices:
        check(all(len(adj[x] & adj[v]) == 2 for x in adj[v]), 'cyclic neighbourhood boundary obstruction')
    check(sum(len(adj[v])*(len(adj[v])-1)//2 for v in g.vertices) == 16*15,
          'cherry inequality exact equality')
    print('Shrikhande16 saturates cherry count; extreme-point geometry excludes its unit embedding.')

    # Cited inputs, not rediscovered by this finite program. The elementary
    # extreme-point proof already excludes critical cores on fewer than17.
    unit_upper = dict(zip(range(16, 31), (41, 43, 46, 50, 54, 57, 61, 66, 72, 78, 84, 90, 96, 103, 110)))
    required = lambda n: (10*n-7+2)//3
    for n in range(17, 27):
        check(required(n) > unit_upper[n], 'cited edge-bound incompatibility')
    check(required(27) <= unit_upper[27], 'bound stops at27; no existence inference')
    print('CITED: Kostochka-Yancey7-critical edge lower bound ceil((10n-7)/3).')
    print('CITED: Alexeev-Mixon-Parshall2412.11914v2 Theorem1 unit-edge upper bounds.')
    print('n25/26/27 (critical requirement, unit upper):', [(n, required(n), unit_upper[n]) for n in (25, 26, 27)])
    print('Necessary six-colour unit-separated equality tester order >=27; no tester constructed.')
    print(f'PASS: {CHECKS} exact checks. No unproved geometry is certified by count receipts.')


if __name__ == '__main__':
    main()
