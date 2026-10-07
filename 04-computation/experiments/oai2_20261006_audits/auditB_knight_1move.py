#!/usr/bin/env python3
"""Audit B, item 6 (HYP-9215 spot check): n x n knight torus G_n, P = one move e0 = (0,0)-(1,2).
 (i)   e0 lies in a Hamiltonian cycle (so beta_P >= 1 is meaningful);
 (ii)  lambda_P <= 7: deleting the 7 other moves at (0,0) leaves (0,0) with degree 1 (explicit local set);
       no 'early closing' for a 1-move path when G_n is triangle-free (n = 6: bipartite);
 (iii) CEGAR lower bound: no set D of <= 6 moves (e0 not in D) meets every Hamiltonian cycle through e0.
       Hitting-set side: CaDiCaL with a totalizer cardinality constraint; cycle side: randomized
       Warnsdorff DFS with degree pruning, falling back to an exact SAT Hamiltonicity check.
Independent code (no use of the reader's kp_* scripts)."""
import itertools, random, sys, time
from pysat.solvers import Solver
from pysat.card import CardEnc, EncType

n = int(sys.argv[1]) if len(sys.argv) > 1 else 6
K = int(sys.argv[2]) if len(sys.argv) > 2 else 6        # test whether a K-set can block
TLIM = float(sys.argv[3]) if len(sys.argv) > 3 else 540
rng = random.Random(20261006)
V = [(x, y) for x in range(n) for y in range(n)]
vid = {v: i for i, v in enumerate(V)}
moves = [(1, 2), (2, 1), (-1, 2), (-2, 1), (1, -2), (2, -1), (-1, -2), (-2, -1)]
E = {}
nbr = [set() for _ in V]
for v in V:
    for dx, dy in moves:
        w = ((v[0] + dx) % n, (v[1] + dy) % n)
        a, b = sorted((vid[v], vid[w]))
        if a != b and (a, b) not in E:
            E[(a, b)] = len(E)
        nbr[vid[v]].add(vid[w])
NE = len(E)
assert all(len(s) == 8 for s in nbr), "not 8-regular"
u0, v0 = vid[(0, 0)], vid[(1, 2)]
e0 = E[tuple(sorted((u0, v0)))]
# triangle check
tri = sum(1 for a in range(len(V)) for b in nbr[a] for c in nbr[b] if a in nbr[c]) // 6
print(f"G_{n}: {len(V)} vertices, {NE} edges, 8-regular, triangles = {tri}")


def eid(a, b):
    return E[(a, b) if a < b else (b, a)]


def ham_dfs(deleted, budget=200000):
    """Hamiltonian cycle containing e0 avoiding 'deleted' (set of edge ids): a Ham path v0 -> ... -> u0."""
    N = len(V)
    adj = [[w for w in nbr[a] if eid(a, w) not in deleted] for a in range(N)]
    if any(len(adj[a]) < 2 for a in range(N)):
        return None, True       # some vertex has degree < 2: certainly no Ham cycle
    visited = [False] * N
    path = [v0]
    visited[v0] = True
    avail = [len(adj[a]) for a in range(N)]   # number of unvisited neighbours (+ current end) approx
    nodes = [0]

    def free_deg(a):
        return sum(1 for w in adj[a] if not visited[w])

    def rec(cur):
        nodes[0] += 1
        if nodes[0] > budget:
            raise TimeoutError
        if len(path) == N:
            return cur == u0
        cand = [w for w in adj[cur] if not visited[w] and (w != u0 or len(path) == N - 1)]
        rng.shuffle(cand)
        cand.sort(key=free_deg)
        for w in cand:
            visited[w] = True
            path.append(w)
            # pruning: every unvisited vertex other than u0 needs >= 2 options (unvisited nbrs or w);
            # u0 needs >= 1 unvisited neighbour or adjacency to w
            ok = True
            for z in adj[w]:
                if not visited[z]:
                    if z == u0:
                        if len(path) < N - 1 and free_deg(z) < 1:
                            ok = False
                            break
                    elif free_deg(z) + 1 < 2:
                        ok = False
                        break
            if ok and rec(w):
                return True
            visited[w] = False
            path.pop()
        return False

    try:
        if rec(v0):
            return list(path), True
        return None, True          # exhaustive search finished: no cycle
    except TimeoutError:
        return None, False


def ham_sat(deleted):
    """exact: Hamiltonian cycle through e0 avoiding deleted, via position SAT encoding."""
    N = len(V)
    var = lambda a, p: a * N + p + 1
    s = Solver(name='cadical153')
    for a in range(N):
        s.append_formula(CardEnc.equals([var(a, p) for p in range(N)], 1, top_id=10**6 + a * 10**4, encoding=EncType.pairwise).clauses)
    for p in range(N):
        s.append_formula(CardEnc.equals([var(a, p) for a in range(N)], 1, top_id=2 * 10**6 + p * 10**4, encoding=EncType.pairwise).clauses)
    s.add_clause([var(v0, 0)])
    s.add_clause([var(u0, N - 1)])
    for p in range(N - 1):
        for a in range(N):
            ok = [var(b, p + 1) for b in nbr[a] if eid(a, b) not in deleted]
            s.add_clause([-var(a, p)] + ok)
    r = s.solve()
    path = None
    if r:
        m = set(l for l in s.get_model() if l > 0)
        path = [None] * N
        for a in range(N):
            for p in range(N):
                if var(a, p) in m:
                    path[p] = a
    s.delete()
    return path


def cyc_edges(path):
    es = {eid(path[i], path[i + 1]) for i in range(len(path) - 1)}
    es.add(eid(path[-1], path[0]))
    return es


# (i)
p, done = ham_dfs(set())
assert p is not None and e0 in cyc_edges(p)
print("(i) e0 lies in a Hamiltonian cycle: yes")
# (ii)
D7 = {eid(u0, w) for w in nbr[u0] if w != v0}
p, done = ham_dfs(D7)
print(f"(ii) local set: the 7 other moves at (0,0) -> Hamiltonian cycle through e0 exists: {p is not None} "
      f"(degree of (0,0) after deletion = {sum(1 for w in nbr[u0] if eid(u0, w) not in D7)})")
# (iii) CEGAR
t0 = time.time()
others = [e for e in range(NE) if e != e0]
xv = {e: i + 1 for i, e in enumerate(others)}
top = len(others)
card = CardEnc.atmost([xv[e] for e in others], bound=K, top_id=top, encoding=EncType.totalizer)
hs = Solver(name='cadical153', bootstrap_with=card.clauses)
pool = 0
# seed with random cycles
for _ in range(400):
    p, done = ham_dfs(set(), budget=50000)
    if p:
        hs.add_clause([xv[e] for e in cyc_edges(p) if e != e0]); pool += 1
it = 0
status = 'TIMEOUT'
fallback = 0
while time.time() - t0 < TLIM:
    it += 1
    if not hs.solve():
        status = 'UNSAT'
        break
    m = set(l for l in hs.get_model() if l > 0)
    D = {e for e in others if xv[e] in m}
    p, done = ham_dfs(D)
    if p is None and not done:
        fallback += 1
        p = ham_sat(D)
    if p is None:
        status = f'BLOCKER FOUND {sorted(D)}'
        break
    hs.add_clause([xv[e] for e in cyc_edges(p) if e != e0]); pool += 1
    # a few more diverse cycles avoiding D
    for _ in range(3):
        q, d2 = ham_dfs(D, budget=20000)
        if q:
            hs.add_clause([xv[e] for e in cyc_edges(q) if e != e0]); pool += 1
    if it % 500 == 0:
        print(f"   ... iteration {it}, pool {pool}, {time.time()-t0:.0f}s", flush=True)
print(f"(iii) n={n}: does a set of <= {K} moves (avoiding e0) block every tour through e0? -> {status}; "
      f"iterations {it}, cycle pool {pool}, SAT fallbacks {fallback}, {time.time()-t0:.0f}s")
if status == 'UNSAT':
    print(f"     => beta_P >= {K+1} for the 1-move path at n = {n}; with (ii), beta_P = lambda_P = {K+1}" if K == 6 else "")
