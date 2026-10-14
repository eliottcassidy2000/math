#!/usr/bin/env python3
"""Vertex-critical 5-chromatic subgraphs of type-M graphs by randomized core-guided deletion.

Graph spec: G = A u w B, where A, B are Minkowski sums of Parts' pieces (V25, V31, V37, V49, H1, H2) and w is one of
the rotations in rotation_probe.w_of (t = squared reference radius).  For each seed: shrink to an UNSAT core over the
vertex selectors, then repeat deletion passes in random order (each failed deletion attempt keeps the vertex; each
successful one replaces the active set by the returned core) until a full pass removes nothing.  The result is
vertex-critical: non-4-colourable, and every single-vertex deletion is 4-colourable (checked at the end).

Each best set is written to best_<tag>.json: exact coordinates (field coefficients as strings in the basis of
hn_udg.PRIMES), the induced edge list, and a 4-colouring certificate of every vertex-deleted subgraph is NOT stored
(re-run verify_critical.py to check).

Usage: python3 -u minimize.py A B t seeds [solver]
  e.g. python3 -u minimize.py V49+V37 V31+V25 4 8
"""
import sys
import json
import time
import random
from fractions import Fraction as Fr

from hn_udg import (PRIMES, H, H_rot, union, rotate, minkowski, unit_graph, triangles_through, Colourer,
                    is_vertex_critical, I_OVER_SQ3, I_SQ3)
from rotation_probe import w_of


def pieces():
    H1 = H_rot([-1, 0, 1])
    H2 = H_rot([-2, -1, 0, 1, 2])
    return {
        'H': H, 'H1': H1, 'H2': H2,
        'V25': union(H1, rotate(H, I_OVER_SQ3)),
        'V31': union(H1, rotate(H_rot([-1, 1]), I_SQ3)),
        'V37': union(H1, rotate(H1, I_OVER_SQ3)),
        'V49': union(H2, rotate(H1, I_OVER_SQ3)),
    }


def build(spec, P):
    names = spec.split('+')
    B = P[names[0]]
    for nm in names[1:]:
        B = minkowski(B, P[nm])
    return B


def log(*a):
    print(*a, flush=True)


def minimize_once(col, n, start, rng, solver='cadical195'):
    """Core-guided deletion passes with the colourer's symmetry-breaking triangle fixed (valid: a colouring exists iff
    one exists with that triangle coloured 0, 1, 2, as long as the triangle is present).  The triangle's own vertices
    are then retested with another fixed triangle inside the remaining set; if one goes, the passes restart."""
    ok, act = col.colourable(start)
    assert not ok
    active = set(act)
    passes = 0
    cur = col
    while True:
        while True:
            passes += 1
            removed = 0
            order = sorted(active)
            rng.shuffle(order)
            for v in order:
                if v not in active or v in cur.fixed:
                    continue
                ok, core = cur.colourable(sorted(active - {v}))
                if not ok:
                    removed += len(active) - len(core)
                    active = set(core)
            if removed == 0:
                break
        dropped = False
        for v in list(cur.fixed):
            rest = active - {v}
            tris = [t for t in triangles_through(n, [e for e in cur.edges if e[0] in rest and e[1] in rest])]
            if not tris:
                continue
            alt = Colourer(n, cur.edges, k=cur.k, solver=solver, fixed_triangle=tris[rng.randrange(len(tris))])
            ok, core = alt.colourable(sorted(rest))
            if not ok:
                active = set(core)
                if cur is not col:
                    cur.close()
                cur = alt
                dropped = True
                break
            alt.close()
        if not dropped:
            if cur is not col:
                cur.close()
            return sorted(active), passes


def main():
    A, B, t = sys.argv[1], sys.argv[2], Fr(sys.argv[3])
    seeds = int(sys.argv[4]) if len(sys.argv) > 4 else 4
    solver = sys.argv[5] if len(sys.argv) > 5 else 'cadical195'
    P = pieces()
    t0 = time.time()
    GA, GB = build(A, P), build(B, P)
    w = w_of(t)
    G = union(GA, rotate(GB, w))
    pts, xy, E = unit_graph(G)
    n = len(pts)
    tris = triangles_through(n, E)
    tag = f"{A}_{B}_t{str(t).replace('/', 'o')}"
    log(f"G = {A} u w_{t} {B}: |V|={n} |E|={len(E)}  [{time.time() - t0:.1f}s]")
    best = None
    for seed in range(seeds):
        rng = random.Random(seed)
        tri = tris[rng.randrange(len(tris))]
        col = Colourer(n, E, k=4, solver=solver, fixed_triangle=tri)
        ok, core = col.colourable(range(n))
        if ok:
            log("   4-colourable: nothing to minimize")
            return
        t1 = time.time()
        crit, passes = minimize_once(col, n, core, rng)
        free = Colourer(n, E, k=4, solver=solver)
        crit_ok = is_vertex_critical(free, crit)
        free.close()
        nA = sum(1 for v in crit if pts[v] in GA)
        log(f"   seed {seed}: vertex-critical {crit_ok}, {len(crit)} vertices ({nA} from A, {len(crit) - nA} from wB), "
            f"{passes} passes  [{time.time() - t1:.1f}s]")
        col.close()
        if crit_ok:
            keep = set(crit)
            sub_edges = [(i, j) for (i, j) in E if i in keep and j in keep]
            out = {
                'spec': {'A': A, 'B': B, 't': str(t), 'primes': list(PRIMES), 'seed': seed},
                'n': len(crit),
                'edges': len(sub_edges),
                'vertices': [[[str(c) for c in pts[v][0]], [str(c) for c in pts[v][1]]] for v in crit],
                'edge_list': [[crit.index(i), crit.index(j)] for (i, j) in sub_edges],
            }
            with open(f"crit_{tag}_seed{seed}.json", 'w') as f:
                json.dump(out, f)
            if best is None or len(crit) < len(best):
                best = crit
                with open(f"best_{tag}.json", 'w') as f:
                    json.dump(out, f)
    log(f"   best: {len(best) if best else None}  [total {time.time() - t0:.1f}s]")


if __name__ == '__main__':
    main()
