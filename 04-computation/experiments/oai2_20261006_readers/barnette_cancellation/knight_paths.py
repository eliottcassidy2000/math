#!/usr/bin/env python3
"""Prescribed-path tour blocking on the knight torus G_n = Cay(Z_n^2, {(+-1,+-2),(+-2,+-1)}).

Motivation (Barnette preprint, openai/math #180): for cubic (k=3) Barnette graphs the preprint's
Corollary (edge avoidance) says beta_0 >= 2 = k-1, and its Pfaffian corollary (P4-Hamiltonicity)
says beta_{P4} >= 1 = k-2, where beta_P = fewest edges (disjoint from P) whose deletion leaves no
Hamiltonian cycle containing the path P.  For the knight torus (k=8), HYP-9211 says beta_0 = 7.
Here we compute beta_P for prescribed paths P with 1, 2, 3 (and 4) edges and compare with the
LOCAL bound lambda(P) (starve one vertex, or force a premature closing at one vertex).

Method: lazy cuts.  master = min hitting set of the pool of Hamiltonian cycles containing P
(CP-SAT); oracle = CP-SAT AddCircuit with P forced and the hitting set deleted.  At termination
the optimal hitting set blocks, so beta_P = its size (the pool is the certificate of >=).
Phase B: forbid 'local' sets (all deleted edges at one vertex): if infeasible at size beta_P,
every minimum blocking set is local.
"""
import sys, time, random, itertools
from ortools.sat.python import cp_model

MOVES = [(1, 2), (2, 1), (-1, 2), (-2, 1), (1, -2), (2, -1), (-1, -2), (-2, -1)]
D4 = [lambda p: (p[0], p[1]), lambda p: (p[1], -p[0]), lambda p: (-p[0], -p[1]), lambda p: (-p[1], p[0]),
      lambda p: (p[1], p[0]), lambda p: (-p[0], p[1]), lambda p: (-p[1], -p[0]), lambda p: (p[0], -p[1])]


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
    V = n * n
    adj = [set() for _ in range(V)]
    inc = [[] for _ in range(V)]
    for k, (u, v) in enumerate(E):
        adj[u].add(v); adj[v].add(u); inc[u].append(k); inc[v].append(k)
    return eid, E, adj, inc


def find_hc(E, nV, removed, forced, seed=0, timeout=120.0, workers=2):
    m = cp_model.CpModel()
    arcs, lit = [], {}
    for k, (u, v) in enumerate(E):
        if k in removed:
            continue
        a = m.NewBoolVar(""); b = m.NewBoolVar("")
        arcs.append((u, v, a)); arcs.append((v, u, b)); lit[k] = (a, b)
    for k in forced:
        if k not in lit:
            return None
        a, b = lit[k]; m.AddBoolOr([a, b])
    m.AddCircuit(arcs)
    if seed:
        rng = random.Random(seed)
        m.Maximize(sum(rng.randint(0, 5) * (a + b) for a, b in lit.values()))
    s = cp_model.CpSolver()
    s.parameters.num_workers = workers
    s.parameters.max_time_in_seconds = timeout if not seed else 3.0
    s.parameters.random_seed = seed
    st = s.Solve(m)
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return [k for k, (a, b) in lit.items() if s.Value(a) or s.Value(b)]
    if st == cp_model.INFEASIBLE:
        return None
    return "UNKNOWN"


def is_hc(E, nV, cyc):
    if len(cyc) != nV:
        return False
    adj = {v: [] for v in range(nV)}
    for k in cyc:
        u, v = E[k]; adj[u].append(v); adj[v].append(u)
    if any(len(adj[v]) != 2 for v in adj):
        return False
    seen, prev, cur = {0}, None, 0
    for _ in range(nV - 1):
        nxt = adj[cur][0] if adj[cur][0] != prev else adj[cur][1]
        if nxt in seen:
            return False
        seen.add(nxt); prev, cur = cur, nxt
    return len(seen) == nV


def master(cand, pool, size=None, nonlocal_rows=None, timeout=600.0):
    m = cp_model.CpModel()
    x = {k: m.NewBoolVar("") for k in cand}
    for C in pool:
        m.AddBoolOr([x[k] for k in C])
    if nonlocal_rows is not None:
        for row in nonlocal_rows:
            r = [x[k] for k in row if k in x]
            if len(r) > size - 1:
                m.Add(sum(r) <= size - 1)
    if size is not None:
        m.Add(sum(x.values()) <= size)
    else:
        m.Minimize(sum(x.values()))
    s = cp_model.CpSolver()
    s.parameters.num_workers = 2
    s.parameters.max_time_in_seconds = timeout
    st = s.Solve(m)
    if st == cp_model.INFEASIBLE:
        return "INFEASIBLE", None
    if st == cp_model.OPTIMAL or (size is not None and st == cp_model.FEASIBLE):
        return "OK", sorted(k for k in cand if s.Value(x[k]))
    return "UNKNOWN", None


def local_bound(n, adj, path):
    """lambda(P): cheapest single-vertex obstruction (starvation or premature closing)."""
    V = n * n
    I = set(path[1:-1]); a, b = path[0], path[-1]; P = set(path)
    best, why = 10 ** 9, None
    m = len(path) - 1
    for v in range(V):
        if v in P:
            continue
        usable = len(adj[v] - I)
        if usable - 1 < best:
            best, why = usable - 1, ("starve", v)
        if a in adj[v] and b in adj[v] and m >= 1:
            c = usable - 2
            if c < best:
                best, why = c, ("close", v)
    for e in (a, b):
        other = b if e == a else a
        usable = len(adj[e] - I - {other}) if m >= 2 else len(adj[e] - {other})
        if usable < best:
            best, why = usable, ("end", e)
    return best, why


def path_edges(eid, path):
    return [eid[(min(u, v), max(u, v))] for u, v in zip(path, path[1:])]


def beta_path(n, eid, E, adj, inc, path, verbose=False, phaseB=True):
    nV = n * n
    EP = path_edges(eid, path)
    EPs = set(EP)
    cand = [k for k in range(len(E)) if k not in EPs]
    pool, seen = [], set()

    def add(c):
        assert is_hc(E, nV, c) and EPs <= set(c)
        key = frozenset(c)
        if key in seen:
            return
        seen.add(key)
        pool.append(sorted(set(c) - EPs))

    for sd in range(1, 9):
        c = find_hc(E, nV, set(), EP, seed=sd)
        if isinstance(c, list):
            add(c)
    if not pool:
        c = find_hc(E, nV, set(), EP)
        if c is None:
            return 0, None, 0, None
        add(c)
    it = 0
    t0 = time.time()
    while True:
        it += 1
        st, B = master(cand, pool)
        assert st == "OK", st
        c = find_hc(E, nV, set(B), EP)
        if c is None:
            beta = len(B); Bblock = B
            break
        assert c != "UNKNOWN"
        add(c)
        # diversify: a second random cycle avoiding B
        c2 = find_hc(E, nV, set(B), EP, seed=it + 100)
        if isinstance(c2, list):
            add(c2)
    res_local = None
    if phaseB:
        rows = [inc[v] for v in range(nV)]
        itb = 0
        while True:
            itb += 1
            st, B = master(cand, pool, size=beta, nonlocal_rows=rows)
            if st == "INFEASIBLE":
                res_local = True
                break
            assert st == "OK", st
            c = find_hc(E, nV, set(B), EP)
            if c is None:
                res_local = ("EXOTIC", [E[k] for k in B])
                break
            assert c != "UNKNOWN"
            add(c)
    return beta, [E[k] for k in Bblock], len(pool), res_local


def path_reps(n, m):
    """representatives of paths with m edges starting at 0 with first step (1,2), up to reversal+D4."""
    reps, seenkeys = [], set()
    for steps in itertools.product(range(8), repeat=m - 1):
        st = [(1, 2)] + [MOVES[s] for s in steps]
        pts = [(0, 0)]
        for d in st:
            pts.append(((pts[-1][0] + d[0]) % n, (pts[-1][1] + d[1]) % n))
        if len(set(pts)) != len(pts):
            continue
        # canonical key under reversal + D4 + translation: normalise both directions
        keys = []
        for seq in (pts, pts[::-1]):
            dsteps = [((q[0] - p[0]) % n, (q[1] - p[1]) % n) for p, q in zip(seq, seq[1:])]
            for g in D4:
                gs = [tuple(c % n for c in g(d)) for d in dsteps]
                keys.append(tuple(gs))
        key = min(keys)
        if key in seenkeys:
            continue
        seenkeys.add(key)
        reps.append(pts)
    return reps


def main():
    n = int(sys.argv[1]); mmax = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    eid, E, adj, inc = knight_torus(n)
    nV = n * n
    tri = sum(1 for (u, v) in E for w in adj[u] & adj[v]) // 3
    print(f"n={n}: V={nV} E={len(E)} triangles={tri}", flush=True)
    for m in range(1, mmax + 1):
        reps = path_reps(n, m)
        print(f"-- paths with {m} edges: {len(reps)} orbit representatives (translations x| D4, reversal)", flush=True)
        summary = {}
        for pts in reps:
            path = [n * p[0] + p[1] for p in pts]
            lam, why = local_bound(n, adj, path)
            t0 = time.time()
            beta, B, psize, loc = beta_path(n, eid, E, adj, inc, path)
            steps = [((q[0] - p[0]) % n, (q[1] - p[1]) % n) for p, q in zip(pts, pts[1:])]
            flag = "OK" if beta == lam else "DIFF"
            summary.setdefault((beta, lam, str(loc)), 0)
            summary[(beta, lam, str(loc))] += 1
            print(f"   P={pts} beta={beta} lambda={lam} ({why[0]}) {flag}; minimum sets all local: {loc}; "
                  f"pool={psize} ({time.time()-t0:.1f}s)", flush=True)
        print(f"   summary (beta, lambda, all-local) -> count: {summary}", flush=True)


if __name__ == "__main__":
    main()
