#!/usr/bin/env python3
"""Knight-torus tour-blocking numbers on Z_n x Z_n, n = 5..10 (six-seven session, 2026-10-06).

G_n = knight graph on the n x n torus (8-regular for n >= 5, edge-transitive under
translations x| D4).  beta(n) = fewest edges whose deletion leaves no Hamiltonian cycle.
The 6x6 case was settled exhaustively in chessboard_weave_20261006_knight_torus.py:
beta(6) = 7 and the only minimum blocking sets strip a square to degree 1.

Method (lazy cuts; the pool of cycles is the certificate):
  master  : min sum x_e  s.t.  every pool cycle C has sum_{e in C} x_e >= 1,  x_{e0} = 1
            (edge-transitivity: any blocking set can be moved to contain e0).
  oracle  : CP-SAT AddCircuit finds a Hamiltonian cycle of G - B, or proves none exists.
  Every cycle found is added together with its images under the point group D4 (all
  translates are implied by re-centring; we add the D4 images and the translates of the
  cycle that keep e0 relevant, i.e. the full orbit restricted to cycles through e0's
  neighbourhood is not needed -- we simply add the whole Z_n^2 x| D4 orbit, deduplicated).
  Phase A ends when the master optimum is >= 7: then every 6-set misses a pool cycle, so
  beta(n) >= 7 (each pool cycle is re-verified as a Hamiltonian cycle of G_n).
  Phase B adds "no star" rows sum_{e at v} x_e <= 6 (a 7-set with 7 edges at one vertex
  IS a star) and asks for a 7-set: either the master becomes infeasible at size 7 (no
  non-star 7-set hits every pool cycle => the stars are the only minimum blocking sets),
  or the oracle certifies a non-star blocking 7-set (an exotic minimum).
Run: python3 sixseven_20261006_knight_torus_block.py [n ...]
"""
import sys, time, random
from ortools.sat.python import cp_model

MOVES = [(1, 2), (2, 1), (-1, 2), (-2, 1), (1, -2), (2, -1), (-1, -2), (-2, -1)]


def knight_torus(n):
    eid, E = {}, []
    for i in range(n):
        for j in range(n):
            u = n * i + j
            for di, dj in MOVES:
                v = n * ((i + di) % n) + (j + dj) % n
                e = (min(u, v), max(u, v))
                if e not in eid:
                    eid[e] = len(E); E.append(e)
    return eid, E


def point_group(n):
    """the 8 elements of D4 acting on Z_n^2 coordinates, plus translations generated separately"""
    return [lambda p: (p[0], p[1]), lambda p: (p[1], -p[0]), lambda p: (-p[0], -p[1]), lambda p: (-p[1], p[0]),
            lambda p: (p[1], p[0]), lambda p: (-p[0], p[1]), lambda p: (-p[1], -p[0]), lambda p: (p[0], -p[1])]


def cycle_orbit(n, eid, cyc_edges):
    """all images of an edge set under Z_n^2 x| D4, as frozensets of edge ids"""
    out = set()
    pts = [(divmod(a, n), divmod(b, n)) for a, b in cyc_edges]
    for g in point_group(n):
        for ti in range(n):
            for tj in range(n):
                img = []
                for (p, q) in pts:
                    gp, gq = g(p), g(q)
                    a = n * ((gp[0] + ti) % n) + (gp[1] + tj) % n
                    b = n * ((gq[0] + ti) % n) + (gq[1] + tj) % n
                    img.append(eid[(min(a, b), max(a, b))])
                out.add(frozenset(img))
    return out


def find_hc(n, E, removed, seed=0, timeout=60.0):
    """Hamiltonian cycle of G - removed via AddCircuit; returns list of edge ids, None if none, 'UNKNOWN' on timeout"""
    m = cp_model.CpModel()
    arcs, lit = [], {}
    for k, (u, v) in enumerate(E):
        if k in removed:
            continue
        a = m.NewBoolVar(""); b = m.NewBoolVar("")
        arcs.append((u, v, a)); arcs.append((v, u, b)); lit[k] = (a, b)
    m.AddCircuit(arcs)
    rng = random.Random(seed)
    if seed:
        m.Maximize(sum(rng.randint(0, 3) * (a + b) for a, b in lit.values()))
    s = cp_model.CpSolver()
    s.parameters.num_workers = 2
    s.parameters.max_time_in_seconds = timeout if not seed else 2.0
    s.parameters.random_seed = seed
    st = s.Solve(m)
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return [k for k, (a, b) in lit.items() if s.Value(a) or s.Value(b)]
    if st == cp_model.INFEASIBLE:
        return None
    return "UNKNOWN"


def is_hc(n, E, cyc):
    V = n * n
    if len(cyc) != V:
        return False
    adj = {v: [] for v in range(V)}
    for k in cyc:
        u, v = E[k]; adj[u].append(v); adj[v].append(u)
    if any(len(adj[v]) != 2 for v in adj):
        return False
    seen, prev, cur = {0}, None, 0
    for _ in range(V - 1):
        nxt = adj[cur][0] if adj[cur][0] != prev else adj[cur][1]
        if nxt in seen:
            return False
        seen.add(nxt); prev, cur = cur, nxt
    return len(seen) == V


def master(nE, pool, e0, stars_excluded=None, size=None, timeout=600.0):
    m = cp_model.CpModel()
    x = [m.NewBoolVar("") for _ in range(nE)]
    m.Add(x[e0] == 1)
    for C in pool:
        m.AddBoolOr([x[k] for k in C])
    if stars_excluded is not None:
        for inc in stars_excluded:
            m.Add(sum(x[k] for k in inc) <= 6)
    if size is not None:
        m.Add(sum(x) <= size)
    else:
        m.Minimize(sum(x))
    s = cp_model.CpSolver()
    s.parameters.num_workers = 4
    s.parameters.max_time_in_seconds = timeout
    st = s.Solve(m)
    if st == cp_model.INFEASIBLE:
        return "INFEASIBLE", None
    if st == cp_model.OPTIMAL or (size is not None and st == cp_model.FEASIBLE):
        return "OK", sorted(k for k in range(nE) if s.Value(x[k]))
    return "UNKNOWN", None


def run(n):
    t0 = time.time()
    eid, E = knight_torus(n)
    V, nE = n * n, len(E)
    assert nE == 4 * V, (n, nE)
    inc = [[] for _ in range(V)]
    for k, (u, v) in enumerate(E):
        inc[u].append(k); inc[v].append(k)
    assert all(len(r) == 8 for r in inc)
    e0 = 0
    pool, added = [], set()

    def add_cycle(c):
        assert is_hc(n, E, c)
        new = 0
        for img in cycle_orbit(n, eid, [E[k] for k in c]):
            if img not in added:
                added.add(img); pool.append(sorted(img)); new += 1
        return new

    for sd in range(1, 6):
        c = find_hc(n, E, set(), seed=sd)
        if isinstance(c, list):
            add_cycle(c)
    print(f"n={n}: V={V} E={nE}; initial pool {len(pool)} cycles", flush=True)
    # Phase A: lower bound 7
    it = 0
    while True:
        it += 1
        st, B = master(nE, pool, e0)
        assert st == "OK", st
        k = len(B)
        if k >= 7:
            print(f"  phase A done after {it} master solves: every 6-set containing e0 misses one of "
                  f"{len(pool)} pool cycles => beta >= 7  ({time.time()-t0:.0f}s)", flush=True)
            break
        c = find_hc(n, E, set(B))
        if c is None:
            print(f"  BLOCKING SET OF SIZE {k}: {[E[j] for j in B]}  => beta <= {k}", flush=True)
            return
        assert c != "UNKNOWN"
        add_cycle(c)
    # Phase B: are stars the only minimum blocking 7-sets?
    star = [i for i in inc[E[e0][0]] if i != inc[E[e0][0]][1]][:7]
    assert e0 in star and len(star) == 7
    assert find_hc(n, E, set(star)) is None, "star should block"
    it, exotic = 0, []
    excl = [inc[v] for v in range(V)]
    while True:
        it += 1
        st, B = master(nE, pool, e0, stars_excluded=excl, size=7)
        if st == "INFEASIBLE":
            print(f"  phase B done after {it} solves: no non-star 7-set containing e0 hits all "
                  f"{len(pool)} pool cycles => the minimum blocking sets are exactly the stars "
                  f"({V}*8 = {8*V} of them)  ({time.time()-t0:.0f}s)", flush=True)
            break
        assert st == "OK", st
        c = find_hc(n, E, set(B))
        if c is None:
            print(f"  EXOTIC minimum blocking 7-set: {[E[j] for j in B]}", flush=True)
            exotic.append(B)
            # exclude this exact set and continue
            pool.append([k for k in range(nE) if k not in B])  # forces a different choice (no-good)
            continue
        assert c != "UNKNOWN"
        add_cycle(c)
    # verify every pool cycle once more (non-no-good entries are Hamiltonian cycles)
    nh = sum(1 for C in pool if len(C) == V and is_hc(n, E, C))
    print(f"  pool re-verified: {nh} Hamiltonian cycles; exotic minima: {len(exotic)}", flush=True)


if __name__ == "__main__":
    ns = [int(a) for a in sys.argv[1:]] or [5, 7, 8]
    for n in ns:
        run(n)
