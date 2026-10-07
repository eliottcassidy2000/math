#!/usr/bin/env python3
"""8x8 board: which knight paths with m <= 3 moves lie in a closed knight tour?  Compare with
LOCAL FORCING CLOSURE (degree-2 forcing, saturation, premature-closing exclusion), iterated.
The preprint's P4-Hamiltonicity is for cubic graphs; on the board the degree-2 corners make it
fail, and the question is whether every failure is local.
Usage: python3 board_paths.py N mmax"""
import sys, time, itertools
from ortools.sat.python import cp_model

MOVES = [(1, 2), (2, 1), (-1, 2), (-2, 1), (1, -2), (2, -1), (-1, -2), (-2, -1)]


def board(N):
    eid, E = {}, []
    adj = [set() for _ in range(N * N)]
    for i in range(N):
        for j in range(N):
            for di, dj in MOVES:
                a, b = i + di, j + dj
                if 0 <= a < N and 0 <= b < N:
                    u, v = N * i + j, N * a + b
                    e = (min(u, v), max(u, v))
                    if e not in eid:
                        eid[e] = len(E); E.append(e); adj[u].add(v); adj[v].add(u)
    return eid, E, adj


def find_hc(E, nV, forced):
    m = cp_model.CpModel()
    arcs, lit = [], {}
    for k, (u, v) in enumerate(E):
        a = m.NewBoolVar(""); b = m.NewBoolVar("")
        arcs.append((u, v, a)); arcs.append((v, u, b)); lit[k] = (a, b)
    for k in forced:
        m.AddBoolOr(list(lit[k]))
    m.AddCircuit(arcs)
    s = cp_model.CpSolver(); s.parameters.num_workers = 2; s.parameters.max_time_in_seconds = 120
    st = s.Solve(m)
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return True
    if st == cp_model.INFEASIBLE:
        return False
    return None


def closure_ok(nV, E, adj, eid, forced_edges):
    """iterated local forcing; False = local contradiction"""
    A = set(range(len(E)))
    Fo = set(forced_edges)
    inc = [[] for _ in range(nV)]
    for k, (u, v) in enumerate(E):
        inc[u].append(k); inc[v].append(k)
    changed = True
    while changed:
        changed = False
        # premature closing: forced edges form paths; forbid edge between endpoints of a forced path
        par = list(range(nV))
        def find(x):
            while par[x] != x:
                par[x] = par[par[x]]; x = par[x]
            return x
        for k in Fo:
            u, v = E[k]
            ru, rv = find(u), find(v)
            if ru == rv:
                # cycle among forced edges: fine only if Hamiltonian
                if len(Fo) != nV:
                    return False
            par[ru] = rv
        degF = [0] * nV
        for k in Fo:
            u, v = E[k]; degF[u] += 1; degF[v] += 1
        if max(degF) > 2:
            return False
        for k in list(A - Fo):
            u, v = E[k]
            if degF[u] == 2 or degF[v] == 2:
                A.discard(k); changed = True
            elif find(u) == find(v) and len(Fo) + 1 < nV:
                A.discard(k); changed = True
        for v in range(nV):
            av = [k for k in inc[v] if k in A]
            if len(av) < 2:
                return False
            if len(av) == 2:
                for k in av:
                    if k not in Fo:
                        Fo.add(k); changed = True
    return True


def main():
    N, mmax = int(sys.argv[1]), int(sys.argv[2])
    eid, E, adj = board(N)
    nV = N * N
    syms = [lambda i, j: (i, j), lambda i, j: (j, N - 1 - i), lambda i, j: (N - 1 - i, N - 1 - j),
            lambda i, j: (N - 1 - j, i), lambda i, j: (j, i), lambda i, j: (N - 1 - i, j),
            lambda i, j: (N - 1 - j, N - 1 - i), lambda i, j: (i, N - 1 - j)]
    T0 = time.time()
    for m in range(1, mmax + 1):
        seen, stats = set(), {"ext": 0, "nonext_local": 0, "nonext_NONLOCAL": 0, "local_says_bad_but_ext": 0}
        examples = []
        for v0 in range(nV):
            def extend(path):
                if len(path) == m + 1:
                    yield list(path); return
                for w in adj[path[-1]]:
                    if w not in path:
                        path.append(w); yield from extend(path); path.pop()
            for path in extend([v0]):
                keys = []
                for seq in (path, path[::-1]):
                    for g in syms:
                        keys.append(tuple(g(*divmod(x, N)) for x in seq))
                key = min(keys)
                if key in seen:
                    continue
                seen.add(key)
                EP = [eid[(min(a, b), max(a, b))] for a, b in zip(path, path[1:])]
                loc = closure_ok(nV, E, adj, eid, EP)
                if not loc:
                    stats["nonext_local"] += 1
                    continue
                r = find_hc(E, nV, EP)
                assert r is not None
                if r:
                    stats["ext"] += 1
                else:
                    stats["nonext_NONLOCAL"] += 1
                    examples.append([divmod(x, N) for x in path])
        print(f"N={N} m={m}: {len(seen)} path classes (D4 x reversal): {stats}; nonlocal examples {examples[:5]} "
              f"({time.time()-T0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
