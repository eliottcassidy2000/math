#!/usr/bin/env python3
"""Type-M probe: G = L727 u w S361 for a chosen rotation w about the shared centre 0 (Parts, arXiv:2010.12665, section 3).

Base graphs (exact, Q(sqrt3, sqrt5, sqrt11)), with eta = exp((i/2) arccos(5/6)) and H the hexagonal wheel:
  H1 = U_{a=-1..1} eta^a H,  H2 = U_{a=-2..2} eta^a H,
  V25 = H1 u (i/sqrt3) H,  V37 = H1 u (i/sqrt3) H1,  V49 = H2 u (i/sqrt3) H1,  V31 = H1 u (i sqrt3) (eta^-1 H u eta H),
  L727 = V49 + V37 (Minkowski),  S361 = V31 + V25.
Rotations w_t = exp(i theta_t), cos theta_t = 1 - 1/(2t): a point at squared radius t and its image are at unit distance.
  rho  = w_4   = (7 + i sqrt15)/8      (Parts / Heule; field adds sqrt15)
  w73  = w_7/3 = (11 + 5 i sqrt3)/14   (radius^2 7/3 orbit present in L727 and S361; stays in Q(sqrt-3))
  eta2 = w_3   = (5 + i sqrt11)/6      (internal to the Moser ring)

Stages: A build + exact edges; B non-4-colourability (plain solve, triangle fixed); C core shrinking (selectors);
D vertex-critical deletion with progress.  Output is unbuffered.

Usage: python3 -u type_m_probe.py ROT [solver=cadical195] [stage=D] [seed=1]
"""
import sys
import time
import random
from fractions import Fraction as Fr

from hn_udg import (H, H_rot, union, rotate, minkowski, unit_graph, check_float_gap, triangles_through,
                    Colourer, shrink_to_core, is_vertex_critical, el, ZERO, is_unit, cnorm,
                    I_OVER_SQ3, I_SQ3, RHO, ETA, cmul)

ROTS = {
    'rho': RHO,
    'w73': (el(r1=Fr(11, 14)), el(r3=Fr(5, 14))),
    'eta2': cmul(ETA, ETA),
}


def log(*a):
    print(*a, flush=True)


def build(rot):
    H1 = H_rot([-1, 0, 1])
    H2 = H_rot([-2, -1, 0, 1, 2])
    V25 = union(H1, rotate(H, I_OVER_SQ3))
    V37 = union(H1, rotate(H1, I_OVER_SQ3))
    V49 = union(H2, rotate(H1, I_OVER_SQ3))
    V31 = union(H1, rotate(H_rot([-1, 1]), I_SQ3))
    L = minkowski(V49, V37)
    S = minkowski(V31, V25)
    assert len(L) == 727 and len(S) == 361
    w = ROTS[rot]
    assert is_unit(w)
    G = union(L, rotate(S, w))
    return L, S, G


def main():
    rot = sys.argv[1] if len(sys.argv) > 1 else 'rho'
    solver = sys.argv[2] if len(sys.argv) > 2 else 'cadical195'
    stage = sys.argv[3] if len(sys.argv) > 3 else 'D'
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    t0 = time.time()
    L, S, G = build(rot)
    pts, xy, E = unit_graph(G)
    n = len(pts)
    deg = [0] * n
    for i, j in E:
        deg[i] += 1
        deg[j] += 1
    tris = triangles_through(n, E)
    log(f"A  rot={rot}: |L|=727 |S|=361 |G|={n} |E|={len(E)} maxdeg={max(deg)} mindeg={min(deg)} "
        f"triangles={len(tris)} float-gap violations={check_float_gap(xy, E)}  [{time.time() - t0:.1f}s]")
    if stage == 'A':
        return
    tri = tris[0]
    col = Colourer(n, E, k=4, solver=solver, fixed_triangle=tri)
    t = time.time()
    ok, core = col.colourable(range(n))
    log(f"B  4-colourable? {ok}" + ("" if ok else f"; first core {len(core)}") + f"  [{time.time() - t:.1f}s]")
    if ok or stage == 'B':
        return
    t = time.time()
    act = core
    for r in range(30):
        ok, c2 = col.colourable(act)
        assert not ok
        log(f"C  round {r}: core {len(c2)}  [{time.time() - t:.1f}s]")
        if len(c2) >= len(act):
            break
        act = c2
    if stage == 'C':
        return
    rng = random.Random(seed)
    active = set(act)
    order = sorted(active, key=lambda v: (deg[v], rng.random()))
    t = time.time()
    tried = 0
    for v in order:
        if v not in active or v in col.fixed:
            continue
        ok, c2 = col.colourable(sorted(active - {v}))
        tried += 1
        if not ok:
            active = set(c2)
        if tried % 25 == 0:
            log(f"D  tried {tried}, active {len(active)}  [{time.time() - t:.1f}s]")
    crit = sorted(active)
    log(f"D  vertex-critical candidate: {len(crit)} vertices, {sum(1 for i, j in E if i in active and j in active)} edges "
        f"[{time.time() - t:.1f}s]; verifying criticality...")
    log(f"D  vertex-critical: {is_vertex_critical(col, crit)}; from L: {sum(1 for v in crit if pts[v] in L)}, "
        f"from wS: {sum(1 for v in crit if pts[v] not in L)}  [total {time.time() - t0:.1f}s]")


if __name__ == '__main__':
    main()
