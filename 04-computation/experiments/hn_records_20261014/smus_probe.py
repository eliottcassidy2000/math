#!/usr/bin/env python3
"""Smallest non-k-colourable induced subgraph (smallest MUS over vertex clauses) via pysat's OptUx (implicit hitting
sets with RC2).  Soft clauses: one 'vertex v gets a colour' clause per candidate vertex; hard clauses: the edge
constraints.  A smallest MUS is then a smallest vertex set whose induced subgraph is not k-colourable.

Usage:
  python3 -u smus_probe.py SPEC k                      # whole Minkowski-sum graph SPEC (e.g. V25, H1+H1)
  python3 -u smus_probe.py typeM A B t k best.json X   # type-M graph; candidates = best.json's vertices + X random others
"""
import sys
import json
import time
import random
from fractions import Fraction as Fr

from pysat.formula import WCNF
from pysat.examples.optux import OptUx

from hn_udg import union, rotate, unit_graph
from minimize import pieces, build


def smus(n, E, k, cand, verbose=0):
    cand = sorted(cand)
    pos = {v: i for i, v in enumerate(cand)}
    x = lambda v, c: pos[v] * k + c + 1
    w = WCNF()
    cs = set(cand)
    for i, j in E:
        if i in cs and j in cs:
            for c in range(k):
                w.append([-x(i, c), -x(j, c)])
    for v in cand:
        w.append([x(v, c) for c in range(k)], weight=1)
    t0 = time.time()
    with OptUx(w, solver='g4', verbose=verbose) as ox:
        mus = ox.compute()
    # soft clauses are numbered 1..len(cand) in insertion order
    return sorted(cand[i - 1] for i in mus), time.time() - t0


def main():
    if sys.argv[1] != 'typeM':
        spec, k = sys.argv[1], int(sys.argv[2])
        P = pieces()
        pts, xy, E = unit_graph(build(spec, P))
        m, dt = smus(len(pts), E, k, range(len(pts)))
        print(f"{spec}: |V|={len(pts)} |E|={len(E)}; smallest non-{k}-colourable induced subgraph: {len(m)} vertices "
              f"[{dt:.1f}s]", flush=True)
        return
    A, B, t, k = sys.argv[2], sys.argv[3], Fr(sys.argv[4]), int(sys.argv[5])
    best = json.load(open(sys.argv[6]))
    extra = int(sys.argv[7]) if len(sys.argv) > 7 else 0
    from rotation_probe import w_of
    P = pieces()
    G = union(build(A, P), rotate(build(B, P), w_of(t)))
    pts, xy, E = unit_graph(G)
    key = {(tuple(str(c) for c in p[0]), tuple(str(c) for c in p[1])): i for i, p in enumerate(pts)}
    base = [key[(tuple(v[0]), tuple(v[1]))] for v in best['vertices']]
    rng = random.Random(0)
    others = [v for v in range(len(pts)) if v not in set(base)]
    cand = base + rng.sample(others, min(extra, len(others)))
    print(f"type M {A} u w_{t} {B}: |V|={len(pts)}; candidates {len(cand)} ({len(base)} from {sys.argv[6]})", flush=True)
    m, dt = smus(len(pts), E, k, cand, verbose=1)
    print(f"   smallest non-{k}-colourable induced subgraph among candidates: {len(m)} vertices [{dt:.1f}s]", flush=True)


if __name__ == '__main__':
    main()
