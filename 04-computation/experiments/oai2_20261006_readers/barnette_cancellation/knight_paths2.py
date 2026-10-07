#!/usr/bin/env python3
"""Prescribed-path tour blocking on the knight torus G_n (k = 8), with the symmetry trick.

beta_P = fewest edges (disjoint from the path P) whose deletion leaves no Hamiltonian cycle
containing P.  lambda_P = cheapest single-vertex ('local') obstruction: starve a vertex (leave
<= 1 usable edge; edges into the interior of P are unusable), force a premature closing at a
common neighbour of the two ends, or starve an end.

Barnette dictionary (k = 3): the preprint's edge-avoidance corollary is beta_empty >= 2 = k-1;
its Pfaffian corollary (P4-Hamiltonicity) is beta_{P4} >= 1 = k-2.  For k = 8, HYP-9211 is
beta_empty = 7 = k-1 (stars only).  We test beta_P = lambda_P for paths with <= 3 or 4 edges.

Lazy cuts with orbit expansion: every Hamiltonian cycle C found by the oracle is stored once;
for a target path P, the pool consists of all g(C) with g in Z_n^2 x| D4 and g^{-1}(P) a
subpath of C.  master = min hitting set of the pool (CP-SAT); oracle = CP-SAT AddCircuit with
P forced and the master's set deleted.  Termination certificate: pool optimum = size of a set
the oracle proves blocking.  Phase B forbids sets with all edges at one vertex; INFEASIBLE at
size beta means every minimum blocking set is local.
Usage: python3 knight_paths2.py n mmax [mmin]
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


def group_perms(n, eid, E):
    perms = []
    for g in D4:
        for ti in range(n):
            for tj in range(n):
                def mp(a):
                    i, j = divmod(a, n); gi, gj = g((i, j))
                    return n * ((gi + ti) % n) + (gj + tj) % n
                pe = []
                for (a, b) in E:
                    x, y = mp(a), mp(b)
                    pe.append(eid[(min(x, y), max(x, y))])
                inv = [0] * len(E)
                for k, kk in enumerate(pe):
                    inv[kk] = k
                perms.append((pe, inv))
    return perms


def find_hc(E, nV, removed, forced, seed=0, timeout=300.0, workers=2):
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


def master(cand, pool, size=None, rows=None, timeout=900.0):
    m = cp_model.CpModel()
    x = {k: m.NewBoolVar("") for k in cand}
    for C in pool:
        m.AddBoolOr([x[k] for k in C])
    if rows is not None:
        for row in rows:
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
    V = n * n
    I = set(path[1:-1]); a, b = path[0], path[-1]; P = set(path)
    m = len(path) - 1
    best, why = 10 ** 9, None
    for v in range(V):
        if v in P:
            continue
        usable = len(adj[v] - I)
        if usable - 1 < best:
            best, why = usable - 1, ("starve", v)
        if a in adj[v] and b in adj[v] and usable - 2 < best:
            best, why = usable - 2, ("close", v)
    for e in (a, b):
        other = b if e == a else a
        usable = len(adj[e] - I - {other})
        if usable < best:
            best, why = usable, ("end", e)
    return best, why


class Engine:
    def __init__(self, n):
        self.n = n
        self.eid, self.E, self.adj, self.inc = knight_torus(n)
        self.nV = n * n
        self.perms = group_perms(n, self.eid, self.E)
        self.cycles = []        # list of (frozenset of edge ids)
        self.cycset = set()
        self.ncalls = 0

    def add_cycle(self, c):
        assert is_hc(self.E, self.nV, c)
        key = frozenset(c)
        if key not in self.cycset:
            self.cycset.add(key); self.cycles.append(key)
            return True
        return False

    def seed(self, k=20):
        for sd in range(1, k + 1):
            c = find_hc(self.E, self.nV, set(), [], seed=sd)
            if isinstance(c, list):
                self.add_cycle(c)

    def pool_for(self, EP, start=0, pool=None, seen=None):
        """images g(C) containing EP, for stored cycles C[start:]"""
        if pool is None:
            pool, seen = [], set()
        EPs = set(EP)
        pre = [[inv[e] for e in EP] for (pe, inv) in self.perms]
        for C in self.cycles[start:]:
            for gi, (pe, inv) in enumerate(self.perms):
                if all(k in C for k in pre[gi]):
                    img = frozenset(pe[k] for k in C)
                    if img not in seen:
                        seen.add(img)
                        pool.append(sorted(img - EPs))
        return pool, seen

    def beta(self, path, phaseB=True):
        E, nV, eid = self.E, self.nV, self.eid
        EP = [eid[(min(u, v), max(u, v))] for u, v in zip(path, path[1:])]
        EPs = set(EP)
        cand = [k for k in range(len(E)) if k not in EPs]
        pool, seen = self.pool_for(EP)
        upto = len(self.cycles)
        if not pool:
            c = find_hc(E, nV, set(), EP)
            if c is None:
                return 0, None, None, 0
            self.add_cycle(c)
            pool, seen = self.pool_for(EP, upto, pool, seen); upto = len(self.cycles)
        it, lb, tm, to = 0, 1, 0.0, 0.0
        while True:
            it += 1
            t1 = time.time()
            st, B = master(cand, pool, size=lb)
            tm += time.time() - t1
            if st == "INFEASIBLE":
                lb += 1
                continue
            assert st == "OK", st
            t1 = time.time()
            c = find_hc(E, nV, set(B), EP); self.ncalls += 1
            to += time.time() - t1
            if c is None:
                beta, Bblock = len(B), B
                assert beta == lb, (beta, lb)
                break
            assert c != "UNKNOWN"
            self.add_cycle(c)
            pool, seen = self.pool_for(EP, upto, pool, seen); upto = len(self.cycles)
            if it % 25 == 0:
                print(f"      it={it} lb={lb} pool={len(pool)} master {tm:.0f}s oracle {to:.0f}s", flush=True)
        loc = None
        if phaseB:
            rows = [self.inc[v] for v in range(nV)]
            exotic = []
            while True:
                st, B = master(cand, pool, size=beta, rows=rows)
                if st == "INFEASIBLE":
                    loc = "all-local" if not exotic else ("exotic", exotic)
                    break
                assert st == "OK", st
                c = find_hc(E, nV, set(B), EP); self.ncalls += 1
                if c is None:
                    exotic.append([E[k] for k in B])
                    pool.append([k for k in cand if k not in B])   # no-good: exclude this exact set
                    if len(exotic) > 20:
                        loc = ("exotic(>20)", exotic[:3]); break
                    continue
                assert c != "UNKNOWN"
                self.add_cycle(c)
                pool, seen = self.pool_for(EP, upto, pool, seen); upto = len(self.cycles)
        return beta, [E[k] for k in Bblock], loc, len(pool)


def path_reps(n, m):
    reps, seenkeys = [], set()
    for steps in itertools.product(range(8), repeat=m - 1):
        st = [(1, 2)] + [MOVES[s] for s in steps]
        pts = [(0, 0)]
        for d in st:
            pts.append(((pts[-1][0] + d[0]) % n, (pts[-1][1] + d[1]) % n))
        if len(set(pts)) != len(pts):
            continue
        keys = []
        for seq in (pts, pts[::-1]):
            dsteps = [((q[0] - p[0]) % n, (q[1] - p[1]) % n) for p, q in zip(seq, seq[1:])]
            for g in D4:
                keys.append(tuple(tuple(c % n for c in g(d)) for d in dsteps))
        key = min(keys)
        if key in seenkeys:
            continue
        seenkeys.add(key)
        reps.append(pts)
    return reps


def main():
    n = int(sys.argv[1]); mmax = int(sys.argv[2]); mmin = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    T0 = time.time()
    eng = Engine(n)
    tri = sum(1 for (u, v) in eng.E for w in eng.adj[u] & eng.adj[v]) // 3
    print(f"n={n}: V={eng.nV} E={len(eng.E)} triangles={tri} |group used|={len(eng.perms)}", flush=True)
    eng.seed(20)
    print(f"  seeded {len(eng.cycles)} Hamiltonian cycles ({time.time()-T0:.0f}s)", flush=True)
    for m in range(mmin, mmax + 1):
        reps = path_reps(n, m)
        print(f"-- prescribed paths with {m} edges: {len(reps)} representatives", flush=True)
        summary = {}
        for pts in reps:
            path = [n * p[0] + p[1] for p in pts]
            lam, why = local_bound(n, eng.adj, path)
            t0 = time.time()
            beta, B, loc, psize = eng.beta(path)
            flag = "beta=lambda" if beta == lam else "DIFF"
            key = (beta, lam, loc if isinstance(loc, str) else "exotic")
            summary[key] = summary.get(key, 0) + 1
            print(f"   P={pts}: beta={beta} lambda={lam} [{why[0]}] {flag}; minimum sets: {loc}; "
                  f"pool={psize}, stored HCs={len(eng.cycles)} ({time.time()-t0:.0f}s)", flush=True)
        print(f"   SUMMARY m={m}: (beta, lambda, minimum-set type) -> #reps: {summary}", flush=True)
    print(f"done ({time.time()-T0:.0f}s); every stored cycle re-verified as an HC on insertion", flush=True)


if __name__ == "__main__":
    main()
