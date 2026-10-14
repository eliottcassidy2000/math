#!/usr/bin/env python3
"""C6-symmetric (orbit-level) minimisation of Parts' type-M graph G = L727 u rho S361, then vertex-level refinement.

Rotation by zeta6 = e^(i pi/3) preserves L727 and S361 (both built from 60-degree-symmetric pieces) and commutes with
rho, so G is C6-invariant.  Phase 1: failed-assumption core over orbit selectors.  Phase 2: random single-orbit deletion
until orbit-critical.  Phase 3: single-vertex deletion (triangle fixed, its vertices retested under another triangle)
until vertex-critical.  Several seeds.

Usage: python3 -u orbit_core_typeM.py [seeds=6]
"""
import sys
import json
import time
import random

from pysat.solvers import Solver

from hn_udg import (H, H_rot, union, rotate, minkowski, unit_graph, triangles_through, Colourer, is_vertex_critical,
                    I_OVER_SQ3, I_SQ3, RHO, el, cmul)
from fractions import Fraction as Fr


def log(*a):
    print(*a, flush=True)


def main():
    seeds = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    t0 = time.time()
    H1 = H_rot([-1, 0, 1]); H2 = H_rot([-2, -1, 0, 1, 2])
    V25 = union(H1, rotate(H, I_OVER_SQ3)); V37 = union(H1, rotate(H1, I_OVER_SQ3))
    V49 = union(H2, rotate(H1, I_OVER_SQ3)); V31 = union(H1, rotate(H_rot([-1, 1]), I_SQ3))
    L = minkowski(V49, V37); S = minkowski(V31, V25)
    G = union(L, rotate(S, RHO))
    pts, xy, E = unit_graph(G)
    n = len(pts)
    idx = {p: i for i, p in enumerate(pts)}
    z6 = (el(r1=Fr(1, 2)), el(r3=Fr(1, 2)))
    rot = [idx[cmul(p, z6)] for p in pts]
    orbit_of = [-1] * n
    orbits = []
    for v in range(n):
        if orbit_of[v] < 0:
            o, w = [], v
            while orbit_of[w] < 0:
                orbit_of[w] = len(orbits); o.append(w); w = rot[w]
            orbits.append(o)
    sz = {}
    for o in orbits:
        sz[len(o)] = sz.get(len(o), 0) + 1
    log(f"G = L727 u rho S361: {n} vertices, {len(E)} edges, C6-invariant, {len(orbits)} orbits {sz}  [{time.time() - t0:.1f}s]")
    best = None
    for seed in range(seeds):
        rng = random.Random(seed)
        tris = triangles_through(n, E)
        tri = tris[rng.randrange(len(tris))]
        fixed = set(orbit_of[v] for v in tri)
        s = Solver(name='cadical195')
        x = lambda v, c: v * 4 + c + 1
        sel = lambda o: n * 4 + o + 1
        for v in range(n):
            s.add_clause([-sel(orbit_of[v])] + [x(v, c) for c in range(4)])
        for i, j in E:
            for c in range(4):
                s.add_clause([-x(i, c), -x(j, c)])
        for c, v in enumerate(tri):
            s.add_clause([x(v, c)])
        for o in fixed:
            s.add_clause([sel(o)])

        def core(active):
            if s.solve(assumptions=[sel(o) for o in sorted(active)]):
                return None
            return set(l - n * 4 - 1 for l in (s.get_core() or []) if l > n * 4) | fixed

        active = core(set(range(len(orbits))))
        while True:
            removed = 0
            order = sorted(active - fixed); rng.shuffle(order)
            for o in order:
                if o in active:
                    c = core(active - {o})
                    if c is not None:
                        removed += len(active) - len(c); active = c
            if removed == 0:
                break
        orb_vertices = sorted(v for o in active for v in orbits[o])
        s.delete()
        # phase 3: vertex-level
        col = Colourer(n, E, k=4, solver='cadical195', fixed_triangle=tri)
        act = set(orb_vertices)
        ok, c0 = col.colourable(sorted(act))
        assert not ok
        act = set(c0)
        cur = col
        while True:
            while True:
                rem = 0
                order = sorted(act - set(cur.fixed)); rng.shuffle(order)
                for v in order:
                    if v in act:
                        ok, cc = cur.colourable(sorted(act - {v}))
                        if not ok:
                            rem += len(act) - len(cc); act = set(cc)
                if rem == 0:
                    break
            dropped = False
            for v in list(cur.fixed):
                rest = act - {v}
                st = triangles_through(n, [e for e in E if e[0] in rest and e[1] in rest])
                if st:
                    alt = Colourer(n, E, k=4, solver='cadical195', fixed_triangle=st[rng.randrange(len(st))])
                    ok, cc = alt.colourable(sorted(rest))
                    if not ok:
                        act, cur, dropped = set(cc), alt, True
                        break
            if not dropped:
                break
        crit = sorted(act)
        free = Colourer(n, E, k=4, solver='cadical195')
        crit_ok = is_vertex_critical(free, crit)
        nL = sum(1 for v in crit if pts[v] in L)
        log(f"seed {seed}: orbit-critical {len(active)} orbits = {len(orb_vertices)} vertices; vertex-critical {crit_ok}: "
            f"{len(crit)} vertices ({nL} from L, {len(crit) - nL} from rho S)  [{time.time() - t0:.1f}s]")
        if crit_ok and (best is None or len(crit) < len(best)):
            best = crit
            keep = set(best)
            json.dump({'n': len(best), 'vertices': [[[str(c) for c in pts[v][0]], [str(c) for c in pts[v][1]]] for v in best],
                       'edge_list': [[best.index(i), best.index(j)] for i, j in E if i in keep and j in keep]},
                      open('best_orbit_typeM.json', 'w'))
    log(f"best: {len(best) if best else None}")


if __name__ == '__main__':
    main()
