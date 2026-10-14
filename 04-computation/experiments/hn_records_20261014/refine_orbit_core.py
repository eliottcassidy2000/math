#!/usr/bin/env python3
"""Vertex-level refinement of an orbit-critical spindle-free subgraph (output of orbit_core.py).

Rebuilds G = M3 u psi P1 exactly, induces the subgraph on the saved vertex_ids, then runs core-guided single-vertex
deletion (symmetry-breaking triangle fixed; the triangle's own vertices retested under another triangle) until
vertex-critical.  Writes the result in verify_graph.py's format (integer numerators X, Y over den in the basis of
Q(sqrt2, sqrt3, sqrt5); induced edges), which verify_graph.py re-checks from scratch.

Usage: HN_PRIMES=2,3,5 python3 -u refine_orbit_core.py orbitcore.json out.json [seed=0]
"""
import os
import sys
import json
import time
import random

from pysat.solvers import Solver

assert os.environ.get('HN_PRIMES') == '2,3,5'
from hn_udg import cmul, cfloat, union, minkowski, triangles_through
import spindle_free_build as sfb
from spindle_free_pieces import gens, clip


def log(*a):
    print(*a, flush=True)


class Oracle:
    def __init__(self, n, E, fixed=()):
        self.k, self.fixed = 4, tuple(fixed)
        self.n = n
        self.s = Solver(name='cadical195')
        x = lambda v, c: v * 4 + c + 1
        for v in range(n):
            self.s.add_clause([-(n * 4 + v + 1)] + [x(v, c) for c in range(4)])
        for i, j in E:
            for c in range(4):
                self.s.add_clause([-x(i, c), -x(j, c)])
        for c, v in enumerate(self.fixed):
            self.s.add_clause([x(v, c)])

    def core(self, active):
        act = set(active) | set(self.fixed)
        if self.s.solve(assumptions=[self.n * 4 + v + 1 for v in sorted(act)]):
            return None
        return set(l - self.n * 4 - 1 for l in (self.s.get_core() or []) if l > self.n * 4) | set(self.fixed)


def main():
    src, out = sys.argv[1], sys.argv[2]
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    t0 = time.time()
    phi0, phi1, psi, M1, M1c = gens()
    M2 = clip(minkowski(M1, M1), 1.0)
    G = union(minkowski(M2, M1), frozenset(cmul(p, psi) for p in minkowski(M2, M1c)))
    pts = sorted(G, key=cfloat)
    keep = sorted(json.load(open(src))['vertex_ids'])
    sub = [pts[v] for v in keep]
    X, Y, den = sfb.to_int(sub)
    I, J = sfb.edges_of(sub, X, Y, den)
    n = len(sub)
    E = list(zip(I.tolist(), J.tolist()))
    log(f"orbit-critical input: {n} vertices, {len(E)} edges  [{time.time() - t0:.1f}s]")
    rng = random.Random(seed)
    tris = triangles_through(n, E)
    orc = Oracle(n, E, fixed=tris[rng.randrange(len(tris))])
    active = orc.core(range(n))
    assert active is not None, "input is 4-colourable"
    while True:
        while True:
            removed = 0
            order = sorted(active - set(orc.fixed))
            rng.shuffle(order)
            t1 = time.time()
            for v in order:
                if v not in active:
                    continue
                c = orc.core(active - {v})
                if c is not None:
                    removed += len(active) - len(c)
                    active = c
            log(f"pass: active {len(active)}, removed {removed}  [{time.time() - t1:.1f}s]")
            if removed == 0:
                break
        dropped = False
        for v in list(orc.fixed):
            rest = active - {v}
            sub_tris = triangles_through(n, [(i, j) for i, j in E if i in rest and j in rest])
            if not sub_tris:
                continue
            alt = Oracle(n, E, fixed=sub_tris[rng.randrange(len(sub_tris))])
            c = alt.core(rest)
            if c is not None:
                active, orc, dropped = c, alt, True
                break
        if not dropped:
            break
    crit = sorted(active)
    idx = {v: i for i, v in enumerate(crit)}
    cs = set(crit)
    sub_E = [[idx[i], idx[j]] for i, j in E if i in cs and j in cs]
    json.dump({'den': den, 'X': X[crit].tolist(), 'Y': Y[crit].tolist(), 'edges': sub_E,
               'source': src, 'spindle_free_field': 'Q(i, sqrt2, sqrt3, sqrt5)'}, open(out, 'w'))
    log(f"RESULT: {len(crit)} vertices, {len(sub_E)} edges (vertex-critical by construction; re-check with verify_graph.py)  "
        f"[total {time.time() - t0:.1f}s]")


if __name__ == '__main__':
    main()
