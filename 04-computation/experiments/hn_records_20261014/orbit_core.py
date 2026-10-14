#!/usr/bin/env python3
"""Orbit-level (C24-symmetric) minimisation of the spindle-free 5-chromatic graph G = M3 u psi P1 (48,121 vertices).

G is invariant under rotation by phi0 = e^(i pi/12) (both pieces are 24-fold symmetric and rotations commute with psi), so
its vertex set splits into the centre and orbits of size 24 (or divisors).  One SAT selector per orbit: an orbit is
present iff its selector is true.  Phase 1 shrinks the active orbit set to the solver's failed-assumption core; phase 2
deletes orbits one at a time (random order, core replacement) until every single-orbit deletion is 4-colourable.
The result is a C24-symmetric non-4-colourable graph that is orbit-critical (not necessarily vertex-critical).
Note: Heule's spindle-free record 1,441 = 1 + 24*60.

Usage: HN_PRIMES=2,3,5 python3 -u orbit_core.py [seed=0] [out_prefix=orbitcore]
"""
import os
import sys
import time
import json
import random

import numpy as np
from pysat.solvers import Solver

assert os.environ.get('HN_PRIMES') == '2,3,5'
from hn_udg import cmul, cfloat, union, minkowski, triangles_through
import spindle_free_build as sfb
from spindle_free_pieces import gens, clip


def log(*a):
    print(*a, flush=True)


def main():
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    prefix = sys.argv[2] if len(sys.argv) > 2 else 'orbitcore'
    t0 = time.time()
    phi0, phi1, psi, M1, M1c = gens()
    M2 = clip(minkowski(M1, M1), 1.0)
    M3 = minkowski(M2, M1)
    P1 = minkowski(M2, M1c)
    G = union(M3, frozenset(cmul(p, psi) for p in P1))
    pts = sorted(G, key=cfloat)
    idx = {p: i for i, p in enumerate(pts)}
    n = len(pts)
    rot = [idx[cmul(p, phi0)] for p in pts]          # exact: KeyError if not invariant
    orbit_of = [-1] * n
    orbits = []
    for v in range(n):
        if orbit_of[v] >= 0:
            continue
        o, w = [], v
        while orbit_of[w] < 0:
            orbit_of[w] = len(orbits)
            o.append(w)
            w = rot[w]
        orbits.append(o)
    sizes = {}
    for o in orbits:
        sizes[len(o)] = sizes.get(len(o), 0) + 1
    log(f"|G| = {n}; C24-invariant; {len(orbits)} orbits, sizes {sizes}  [{time.time() - t0:.1f}s]")
    X, Y, den = sfb.to_int(pts)
    I, J = sfb.edges_of(pts, X, Y, den)
    E = list(zip(I.tolist(), J.tolist()))
    log(f"|E| = {len(E)}  [{time.time() - t0:.1f}s]")
    tri = triangles_through(n, E, limit=1)[0]
    fixed_orbits = sorted(set(orbit_of[v] for v in tri))
    k = 4
    s = Solver(name='cadical195')
    x = lambda v, c: v * k + c + 1
    sel = lambda o: n * k + o + 1
    for v in range(n):
        s.add_clause([-sel(orbit_of[v])] + [x(v, c) for c in range(k)])
    for i, j in E:
        for c in range(k):
            s.add_clause([-x(i, c), -x(j, c)])
    for c, v in enumerate(tri):
        s.add_clause([x(v, c)])
    for o in fixed_orbits:
        s.add_clause([sel(o)])

    def unsat_core(active):
        ok = s.solve(assumptions=[sel(o) for o in sorted(active)])
        if ok:
            return None
        core = s.get_core() or []
        return set(l - n * k - 1 for l in core if l > n * k) | set(fixed_orbits)

    def size(active):
        return sum(len(orbits[o]) for o in active)

    active = set(range(len(orbits)))
    for rnd in range(30):
        t1 = time.time()
        core = unsat_core(active)
        assert core is not None, "colourable?"
        log(f"phase 1 round {rnd}: core {len(core)} orbits, {size(core)} vertices  [{time.time() - t1:.1f}s]")
        if len(core) >= len(active):
            break
        active = core
    rng = random.Random(seed)
    while True:
        removed = 0
        order = sorted(active - set(fixed_orbits))
        rng.shuffle(order)
        t1 = time.time()
        for o in order:
            if o not in active:
                continue
            core = unsat_core(active - {o})
            if core is not None:
                removed += len(active) - len(core)
                active = core
        log(f"phase 2 pass: {len(active)} orbits, {size(active)} vertices, removed {removed}  [{time.time() - t1:.1f}s]")
        with open(f"{prefix}_seed{seed}.json", 'w') as f:
            json.dump({'den': den, 'orbits': [orbits[o] for o in sorted(active)],
                       'X': X.tolist(), 'Y': Y.tolist()} if False else
                      {'den': den, 'n_vertices': size(active), 'n_orbits': len(active),
                       'vertex_ids': sorted(v for o in active for v in orbits[o])}, f)
        if removed == 0:
            break
    log(f"RESULT seed {seed}: {len(active)} orbits, {size(active)} vertices (orbit-critical under the fixed triangle's orbits)  "
        f"[total {time.time() - t0:.1f}s]")


if __name__ == '__main__':
    main()
