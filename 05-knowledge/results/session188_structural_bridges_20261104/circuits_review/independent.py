#!/usr/bin/env python3
"""Independent incidence check including zero-, one-, and two-pin gates.
Does not import the producer checker.
"""
from itertools import combinations_with_replacement, product
import json


def check(ok, why):
    if not ok:
        raise RuntimeError(why)


def examine(parents, C):
    g = set(range(2, 2 + len(parents)))
    out = max(g)
    ps = {v: set(parents[v-2]) for v in g}
    arcs = {(u, v) for v in g for u in ps[v] if u in g}
    border = {v for u, v in arcs if u in C and v not in C}
    ind = {v: sum(u in C for u in ps[v]) for v in C}
    deg = {u: sum(v == u for v, _ in arcs) for u in C}
    r = sum(d == 0 for d in ind.values())
    f = sum(d == 1 for d in ind.values())
    surplus = sum(d-1 for d in deg.values())
    check(len(C) + len(border) == 2*r + f + surplus, 'identity')
    check(surplus >= 0, 'negative surplus')
    eq = len(C) + len(border) == 2*r
    topology_eq = (all(d == 1 for d in deg.values())
                   and all(d in (0, 2) for d in ind.values()))
    check(eq == topology_eq, 'equality degree criterion')
    if eq:
        unseen = set(C)
        exits_seen = set()
        while unseen:
            seed = min(unseen)
            component, frontier = {seed}, [seed]
            while frontier:
                x = frontier.pop()
                adj = {v for u, v in arcs if u == x and v in C}
                adj |= {u for u, v in arcs if v == x and u in C}
                new = adj - component
                component |= new
                frontier.extend(new)
            unseen -= component
            edges = {(u, v) for u, v in arcs if u in component and v in component}
            exit_edges = {(u, v) for u, v in arcs if u in component and v not in C}
            leaves = sum(ind[v] == 0 for v in component)
            check(len(edges) == len(component)-1, 'underlying component not a tree')
            check(len(component) == 2*leaves-1, 'not full binary')
            check(len(exit_edges) == 1, 'exit nonunique')
            child = next(iter(exit_edges))[1]
            check(child not in exits_seen, 'components share a boundary child')
            exits_seen.add(child)
        check(exits_seen == border, 'missing boundary')
    return eq


def main():
    rows = []
    for length in range(1, 5):
        slots = []
        for k in range(length):
            vertices = range(k + 2)
            slots.append([()] + [(v,) for v in vertices]
                         + list(combinations_with_replacement(vertices, 2)))
        total = live = regions = equality = 0
        for parents in product(*slots):
            total += 1
            g = set(range(2, 2 + length))
            seen, frontier = {max(g)}, [max(g)]
            while frontier:
                v = frontier.pop()
                for u in parents[v-2]:
                    if u in g and u not in seen:
                        seen.add(u)
                        frontier.append(u)
            if seen != g:
                continue
            live += 1
            for mask in range(1, 1 << (length-1)):
                C = {2+i for i in range(length-1) if mask & (1 << i)}
                if any(set(parents[v-2]) <= C and v not in C for v in g):
                    continue
                equality += examine(parents, C)
                regions += 1
        rows.append(dict(gates=length, presentations=total, live=live,
                         closed_regions=regions, equality_cases=equality))
    # C need not be all semantic constants. Boundary output may be constant.
    check(examine(((0, 1), (2, 0)), {2}), 'constant-output structural control')
    # A=a AND b; O=A AND x; restrict a=0. Every surviving assignment outputs zero.
    check(all(((0 & b) & x) == 0 for b, x in product((0, 1), repeat=2)),
          'constant output semantic control')
    print(json.dumps(dict(status='PASS', independent_universe='all live at-most-two-pin DAGs, 2 external vertices, through 4 gates',
                          rows=rows, constant_boundary_output=True,
                          endpoint_note='If terminal constants are charged, at most one replacement constant-or-NOT gate is needed.'),
                     indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
