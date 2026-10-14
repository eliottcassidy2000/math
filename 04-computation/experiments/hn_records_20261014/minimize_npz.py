#!/usr/bin/env python3
"""Core-guided minimisation of a large exact unit-distance graph saved by spindle_free_build.py (npz: X, Y integer
numerators over den in the basis of Q(sqrt2, sqrt3, sqrt5); I, J exact unit edges).

Phase 1 (core shrinking): solve 4-colourability under assumptions 'every active vertex is present'; replace the
active set by the solver's failed-assumption core; repeat while it shrinks.
Phase 2 (deletion): randomized single-vertex deletion passes with core replacement, then a final pass without the
symmetry-breaking triangle, until vertex-critical.  The best set is saved as JSON (integer numerators, den, edges).

Usage: python3 -u minimize_npz.py graph.npz out_prefix [seed=0] [phase2_threshold=6000]
"""
import sys
import json
import time
import random

import numpy as np
from pysat.solvers import Solver


def log(*a):
    print(*a, flush=True)


class Oracle:
    def __init__(self, n, E, k=4, fixed=(), solver='cadical195'):
        self.k, self.fixed = k, tuple(fixed)
        self.s = Solver(name=solver)
        self.sel = np.arange(n * k + 1, n * k + n + 1)
        x = lambda v, c: v * k + c + 1
        for v in range(n):
            self.s.add_clause([-int(self.sel[v])] + [x(v, c) for c in range(k)])
        for i, j in E:
            for c in range(k):
                self.s.add_clause([-x(i, c), -x(j, c)])
        for c, v in enumerate(self.fixed):
            self.s.add_clause([x(v, c)])
        self.base = n * k

    def colourable(self, active):
        act = set(active) | set(self.fixed)
        ok = self.s.solve(assumptions=[int(self.sel[v]) for v in sorted(act)])
        if ok:
            return True, None
        core = self.s.get_core() or []
        return False, sorted(set(int(l) - self.base - 1 for l in core if l > self.base) | set(self.fixed))


def triangle(n, E):
    adj = [set() for _ in range(n)]
    for i, j in E:
        adj[i].add(j)
        adj[j].add(i)
    for i, j in E:
        common = adj[i] & adj[j]
        if common:
            return (i, j, min(common))
    return None


def main():
    path, prefix = sys.argv[1], sys.argv[2]
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    thr = int(sys.argv[4]) if len(sys.argv) > 4 else 6000
    d = np.load(path)
    X, Y, den, I, J = d['X'], d['Y'], int(d['den'][0]), d['I'], d['J']
    n = len(X)
    E = list(zip(I.tolist(), J.tolist()))
    tri = triangle(n, E)
    log(f"graph: n={n} |E|={len(E)} triangle={tri}")
    t0 = time.time()
    orc = Oracle(n, E, fixed=tri)
    act = list(range(n))
    rnd = 0
    while True:
        t1 = time.time()
        ok, core = orc.colourable(act)
        if ok:
            log("colourable: stop")
            return
        log(f"phase 1 round {rnd}: core {len(core)}  [{time.time() - t1:.1f}s]")
        rnd += 1
        if len(core) >= len(act) or rnd > 40:
            break
        act = core
    if len(act) > thr:
        log(f"core {len(act)} above threshold {thr}: stop after phase 1")
        np.save(f"{prefix}_phase1_core.npy", np.array(sorted(act), dtype=np.int64))
        return
    rng = random.Random(seed)
    # phase 2 on the induced subgraph, renumbered
    keep = sorted(act)
    pos = {v: i for i, v in enumerate(keep)}
    E2 = [(pos[i], pos[j]) for i, j in E if i in pos and j in pos]
    tri2 = tuple(pos[v] for v in tri)
    orc2 = Oracle(len(keep), E2, fixed=tri2)
    active = set(range(len(keep)))
    while True:
        removed = 0
        order = sorted(active)
        rng.shuffle(order)
        t1 = time.time()
        for v in order:
            if v not in active or v in tri2:
                continue
            ok, core = orc2.colourable(sorted(active - {v}))
            if not ok:
                removed += len(active) - len(core)
                active = set(core)
        log(f"phase 2 pass: active {len(active)}, removed {removed}  [{time.time() - t1:.1f}s]")
        if removed == 0:
            break
    free = Oracle(len(keep), E2)
    changed = True
    while changed:
        changed = False
        for v in sorted(active):
            if v not in active:
                continue
            ok, core = free.colourable(sorted(active - {v}))
            if not ok:
                active = set(core)
                changed = True
    crit = sorted(active)
    ok, _ = free.colourable(crit)
    assert not ok
    critical = all(free.colourable([u for u in crit if u != v])[0] for v in crit)
    verts = [keep[v] for v in crit]
    vset = set(verts)
    sub = [(i, j) for i, j in E if i in vset and j in vset]
    idx = {v: k for k, v in enumerate(verts)}
    with open(f"{prefix}_seed{seed}.json", 'w') as f:
        json.dump({'den': den, 'basis': ['1', 'sqrt2', 'sqrt3', 'sqrt6', 'sqrt5', 'sqrt10', 'sqrt15', 'sqrt30'],
                   'X': X[verts].tolist(), 'Y': Y[verts].tolist(),
                   'edges': [[idx[i], idx[j]] for i, j in sub]}, f)
    log(f"RESULT seed {seed}: {len(crit)} vertices, {len(sub)} edges, vertex-critical {critical}  [total {time.time() - t0:.1f}s]")


if __name__ == '__main__':
    main()
