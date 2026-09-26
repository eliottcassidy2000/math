"""procgen_smallgraph_20260926_lib.py -- shared helpers for the smallgraph lane.

Session collatz-procgen-20260922, lane "smallgraph" (2026-09-26).

Target-sum graphs: for a target set S of positive integers and n >= 1,
G_S(n) has vertex set {1..n} and an edge {x, y} iff x != y and x + y in S.
Q_n = G_S(n) for S = squares.

Everything here is deterministic.  Hamiltonian searches are delegated to the C
program procgen_smallgraph_20260926_ham.c (compiled into scratch/ on first use);
every witness it returns is re-verified here in pure Python before use.
"""
import hashlib
import os
import subprocess
from math import isqrt

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
SCRATCH = os.path.join(REPO, 'scratch', 'procgen_smallgraph')
HAM_SRC = os.path.join(HERE, 'procgen_smallgraph_20260926_ham.c')
HAM_BIN = os.path.join(SCRATCH, 'bin', 'ham')


class CheckFailed(Exception):
    pass


N_CHECKS = [0]


def check(name, cond, detail=''):
    """Every printed claim goes through here; a failed claim raises."""
    if not cond:
        raise CheckFailed(f'{name}: {detail}')
    N_CHECKS[0] += 1
    print(f'[OK] {name}: {detail}', flush=True)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        h.update(f.read())
    return h.hexdigest()


# ----------------------------------------------------------------------------------
# target sets (all return a sorted list of targets <= limit)
# ----------------------------------------------------------------------------------
def is_square(v):
    return v >= 0 and isqrt(v) ** 2 == v


def squares(limit):
    return [k * k for k in range(1, isqrt(limit) + 1)]


def icbrt(v):
    r = int(round(v ** (1.0 / 3))) if v > 0 else 0
    while r ** 3 > v:
        r -= 1
    while (r + 1) ** 3 <= v:
        r += 1
    return r


def cubes(limit):
    return [k ** 3 for k in range(1, icbrt(limit) + 1)]


def triangular(limit):
    out, k = [], 1
    while k * (k + 1) // 2 <= limit:
        out.append(k * (k + 1) // 2)
        k += 1
    return out


def powers(base, limit):
    out, v = [], 1
    while v <= limit:
        out.append(v)
        v *= base
    return out


def pow23(limit):
    return sorted(set(powers(2, limit)) | set(powers(3, limit)))


def perfect_powers(limit):
    out = {1}
    a = 2
    while a * a <= limit:
        v = a * a
        while v <= limit:
            out.add(v)
            v *= a
        a += 1
    return sorted(out)


def shifted_squares(c, limit):
    """{k^2 + c : k >= 1} intersected with [1, limit]."""
    out = []
    k = 1
    while k * k + c <= limit:
        if k * k + c >= 1:
            out.append(k * k + c)
        k += 1
    return sorted(set(out))


def pentagonal(limit):
    out, k = [], 1
    while k * (3 * k - 1) // 2 <= limit:
        out.append(k * (3 * k - 1) // 2)
        k += 1
    return out


FAMILIES = {
    'squares': squares,
    'cubes': cubes,
    'triangular': triangular,
    'pow2': lambda L: powers(2, L),
    'pow3': lambda L: powers(3, L),
    'pow2or3': pow23,
    'perfect_powers': perfect_powers,
    'sq_plus_1': lambda L: shifted_squares(1, L),
    'sq_minus_1': lambda L: shifted_squares(-1, L),
    'pentagonal': pentagonal,
}


def family(name):
    if name in FAMILIES:
        return FAMILIES[name]
    if name.startswith('sq_c'):  # sq_c<int>, e.g. sq_c-5
        c = int(name[4:])
        return lambda L, c=c: shifted_squares(c, L)
    raise KeyError(name)


# ----------------------------------------------------------------------------------
# graphs
# ----------------------------------------------------------------------------------
def sum_graph(n, targets):
    """Edge list of G_S(n) as sorted pairs (x, y), x < y."""
    tset = sorted(t for t in set(targets) if 3 <= t <= 2 * n - 1)
    edges = []
    for s in tset:
        lo = max(1, s - n)
        hi = (s - 1) // 2
        for x in range(lo, hi + 1):
            y = s - x
            if x < y <= n:
                edges.append((x, y))
    edges.sort()
    return edges


def adjacency(n, edges):
    adj = {v: [] for v in range(1, n + 1)}
    for x, y in edges:
        adj[x].append(y)
        adj[y].append(x)
    for v in adj:
        adj[v].sort()
    return adj


def degree_formula(n, x, tset):
    """|S cap [x+1, x+n]| - [2x in S]  (Lemma 1)."""
    c = sum(1 for s in tset if x + 1 <= s <= x + n)
    return c - (1 if (2 * x) in tset else 0)


def components(n, edges):
    parent = list(range(n + 1))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    for x, y in edges:
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry
    comps = {}
    for v in range(1, n + 1):
        comps.setdefault(find(v), []).append(v)
    return sorted(comps.values(), key=lambda c: c[0])


def leaves(n, edges):
    """Vertices of degree <= 1."""
    deg = [0] * (n + 1)
    for x, y in edges:
        deg[x] += 1
        deg[y] += 1
    return [v for v in range(1, n + 1) if deg[v] <= 1], deg


def verify_path(n, seq, tset, cycle=False):
    if sorted(seq) != list(range(1, n + 1)):
        return False
    for i in range(n - 1):
        if seq[i] + seq[i + 1] not in tset:
            return False
    if cycle and n >= 3 and seq[0] + seq[-1] not in tset:
        return False
    return True


# ----------------------------------------------------------------------------------
# C solver interface
# ----------------------------------------------------------------------------------
def ensure_ham():
    os.makedirs(os.path.dirname(HAM_BIN), exist_ok=True)
    need = (not os.path.exists(HAM_BIN)) or os.path.getmtime(HAM_BIN) < os.path.getmtime(HAM_SRC)
    if need:
        subprocess.run(['gcc', '-O2', '-o', HAM_BIN, HAM_SRC], check=True)
    return HAM_BIN


def run_ham(mode, n, edges, args=()):
    ensure_ham()
    inp = f'{n} {len(edges)}\n' + ''.join(f'{x} {y}\n' for x, y in edges)
    r = subprocess.run([HAM_BIN, mode] + [str(a) for a in args], input=inp,
                       capture_output=True, text=True, check=True)
    return r.stdout.strip().splitlines()


def parse_result(lines):
    """Returns (status, seq or None, info)."""
    first = lines[0].split()
    if first[0] != 'RESULT':
        raise ValueError(lines)
    st = first[1]
    if st in ('PATH', 'CYCLE'):
        nn = int(first[2])
        seq = [int(t) for t in first[3:3 + nn]]
        return st, seq, ' '.join(first[3 + nn:])
    return st, None, ' '.join(first[2:])


def ham_exist(n, edges, cycle=False, limit=0):
    lines = run_ham('exist_cycle' if cycle else 'exist_path', n, edges, (['-L', limit] if limit else []))
    return parse_result(lines)


def ham_count(n, edges, cycle=False, limit=0):
    lines = run_ham('count_cycles' if cycle else 'count_paths', n, edges, (['-L', limit] if limit else []))
    head = lines[0].split()
    assert head[0] == 'COUNT'
    cnt = int(head[1])
    complete = head[2] == 'COMPLETE'
    eps = {}
    for ln in lines[1:]:
        p = ln.split()
        if p[0] == 'EP':
            eps[(int(p[1]), int(p[2]))] = int(p[3])
    return cnt, complete, eps


def ham_count_capped(n, edges, cap, limit=0):
    """Path counting that stops once `cap` raw paths are seen.  Returns (count, status, eps) with
    status COMPLETE (exact count), CAPPED (at least 2 distinct paths when cap >= 3) or ABORTED."""
    args = ['-K', cap] + (['-L', limit] if limit else [])
    lines = run_ham('count_paths', n, edges, args)
    head = lines[0].split()
    eps = {}
    for ln in lines[1:]:
        p = ln.split()
        if p[0] == 'EP':
            eps[(int(p[1]), int(p[2]))] = int(p[3])
    return int(head[1]), head[2], eps


def ham_heur(n, edges, cycle=False, seed=1, restarts=200, budget=200000, start=None, interior=()):
    args = ['-S', seed, '-R', restarts, '-B', budget]
    if start is not None:
        args += ['-s', start]
    if interior:
        args += ['-f', ','.join(str(v) for v in interior)]
    lines = run_ham('heur_cycle' if cycle else 'heur_path', n, edges, args)
    return parse_result(lines)


# ----------------------------------------------------------------------------------
# Independent exact method 2: Hamiltonian-cycle search by edge forcing + branching
# (Python; trail-based undo).  Sound rules (every vertex of a Hamiltonian cycle has
# exactly two cycle edges):
#   F1  a vertex with exactly two available edges forces both;
#   F2  a vertex with two forced edges loses all its other edges;
#   F3  the edge joining the two ends of a forced segment that does not contain
#       every vertex is deleted (it would close a short cycle);
#   F4  fewer than two available edges, or three forced edges, is a contradiction.
# Branching forces or deletes one edge at a vertex of least available degree, so
# the search is exhaustive.  A Hamiltonian PATH of G is a Hamiltonian cycle of G + z
# (a new vertex z = 0 joined to every vertex).
# ----------------------------------------------------------------------------------
class _Contra(Exception):
    pass


class HCForcing:
    def __init__(self, nv, edges):
        self.nv = nv
        self.avail = [set() for _ in range(nv)]
        for u, v in edges:
            self.avail[u].add(v)
            self.avail[v].add(u)
        self.fdeg = [0] * nv
        self.fadj = [[] for _ in range(nv)]
        self.other = list(range(nv))
        self.size = [1] * nv
        self.trail = []
        self.nodes = 0
        self.closed = False

    # -- primitive modifications (all trailed) --
    def _delete(self, u, v, queue):
        self.avail[u].discard(v)
        self.avail[v].discard(u)
        self.trail.append(('d', u, v))
        queue.append(u)
        queue.append(v)

    def _force(self, u, v, queue):
        if self.fdeg[u] >= 2 or self.fdeg[v] >= 2:
            raise _Contra()
        a, b = self.other[u], self.other[v]
        if a == v:  # u, v are the two ends of one segment
            if self.size[u] != self.nv:
                raise _Contra()
            self.trail.append(('c',))
            self.closed = True
        else:
            ns = self.size[u] + self.size[v]
            self.trail.append(('m', a, b, self.other[a], self.other[b], self.size[a], self.size[b]))
            self.other[a], self.other[b] = b, a
            self.size[a] = self.size[b] = ns
        self.fdeg[u] += 1
        self.fdeg[v] += 1
        self.fadj[u].append(v)
        self.fadj[v].append(u)
        self.trail.append(('f', u, v))
        queue.extend((u, v, a, b))
        if not self.closed and a != v and self.size[a] < self.nv and b in self.avail[a] and b not in self.fadj[a]:
            self._delete(a, b, queue)

    def _undo(self, mark):
        while len(self.trail) > mark:
            op = self.trail.pop()
            if op[0] == 'd':
                _, u, v = op
                self.avail[u].add(v)
                self.avail[v].add(u)
            elif op[0] == 'f':
                _, u, v = op
                self.fdeg[u] -= 1
                self.fdeg[v] -= 1
                self.fadj[u].pop()
                self.fadj[v].pop()
            elif op[0] == 'm':
                _, a, b, oa, ob, sa, sb = op
                self.other[a], self.other[b] = oa, ob
                self.size[a], self.size[b] = sa, sb
            elif op[0] == 'c':
                self.closed = False

    def _propagate(self, queue):
        while queue:
            w = queue.pop()
            av = self.avail[w]
            if len(av) < 2 or self.fdeg[w] > 2:
                raise _Contra()
            if self.fdeg[w] == 2 and len(av) > 2:
                keep = set(self.fadj[w])
                for x in list(av):
                    if x not in keep:
                        self._delete(w, x, queue)
            elif len(av) == 2 and self.fdeg[w] < 2:
                for x in list(av):
                    if x not in self.fadj[w]:
                        self._force(w, x, queue)

    def solve(self, node_limit=0):
        """Returns ('FOUND', successor-structure) / ('NONE', None) / ('UNKNOWN', None)."""
        self.aborted = False
        try:
            self._propagate(list(range(self.nv)))
        except _Contra:
            return 'NONE', None
        r = self._rec(node_limit)
        if r:
            return 'FOUND', [list(x) for x in self.fadj]
        return ('UNKNOWN' if self.aborted else 'NONE'), None

    def _rec(self, node_limit):
        self.nodes += 1
        if node_limit and self.nodes > node_limit:
            self.aborted = True
            return False
        if self.closed:
            return all(d == 2 for d in self.fdeg)
        best, bv = None, -1
        for v in range(self.nv):
            if self.fdeg[v] < 2:
                k = len(self.avail[v])
                if best is None or k < best:
                    best, bv = k, v
        if bv < 0:
            return False
        x = min(y for y in self.avail[bv] if y not in self.fadj[bv])
        for choice in ('force', 'delete'):
            mark = len(self.trail)
            try:
                q = []
                if choice == 'force':
                    self._force(bv, x, q)
                else:
                    self._delete(bv, x, q)
                self._propagate(q)
                if self._rec(node_limit):
                    return True
            except _Contra:
                pass
            self._undo(mark)
            if self.aborted:
                return False
        return False


def forcing_ham(n, edges, cycle=False, node_limit=0):
    """Exact Hamiltonian path (cycle=False) or cycle search by forcing.
    Returns (status, seq, nodes) with status in PATH/CYCLE/NONE/UNKNOWN."""
    if n == 1:
        return ('PATH', [1], 0) if not cycle else ('NONE', None, 0)
    if cycle:
        if n < 3:
            return 'NONE', None, 0
        solver = HCForcing(n, [(u - 1, v - 1) for u, v in edges])   # relabel v -> v-1
        st, fadj = solver.solve(node_limit)
        if st != 'FOUND':
            return st, None, solver.nodes
        seq = [0]
        prev = fadj[0][0]
        cur = 0
        while len(seq) < n:
            a, b = fadj[cur]
            nxt = b if a == prev else a
            prev, cur = cur, nxt
            seq.append(cur)
        return 'CYCLE', [v + 1 for v in seq], solver.nodes
    E2 = list(edges) + [(0, v) for v in range(1, n + 1)]
    solver = HCForcing(n + 1, E2)
    st, fadj = solver.solve(node_limit)
    if st != 'FOUND':
        return st, None, solver.nodes
    a, b = fadj[0]
    seq = [a]
    prev, cur = 0, a
    while True:
        x, y = fadj[cur]
        nxt = y if x == prev else x
        if nxt == 0:
            break
        prev, cur = cur, nxt
        seq.append(cur)
    return 'PATH', seq, solver.nodes


# ----------------------------------------------------------------------------------
# Independent exact method 3: OR-tools CP-SAT circuit model (cross-check only).
# ----------------------------------------------------------------------------------
def cpsat_ham(n, edges, cycle=False, timeout=120.0, workers=2, mem_mb=400):
    from ortools.sat.python import cp_model
    m = cp_model.CpModel()
    arcs = []
    for (u, v) in edges:
        arcs.append((u, v, m.NewBoolVar('')))
        arcs.append((v, u, m.NewBoolVar('')))
    if not cycle:
        for v in range(1, n + 1):
            arcs.append((0, v, m.NewBoolVar('')))
            arcs.append((v, 0, m.NewBoolVar('')))
    if not arcs:
        return 'NONE', None
    m.AddCircuit(arcs)
    s = cp_model.CpSolver()
    s.parameters.num_workers = workers
    s.parameters.max_time_in_seconds = timeout
    s.parameters.max_memory_in_mb = mem_mb
    s.parameters.random_seed = 1
    st = s.Solve(m)
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        succ = {}
        for (u, v, l) in arcs:
            if s.Value(l):
                succ[u] = v
        if cycle:
            seq = [1]
            while len(seq) < n:
                seq.append(succ[seq[-1]])
            return 'CYCLE', seq
        seq = [succ[0]]
        while succ[seq[-1]] != 0:
            seq.append(succ[seq[-1]])
        return 'PATH', seq
    if st == cp_model.INFEASIBLE:
        return 'NONE', None
    return 'UNKNOWN', None


def path_uniqueness(n, edges, node_limit=0):
    """Exact: returns ('NONE', None, None), ('UNIQUE', P, None) or ('MULTIPLE', P, Q) where P, Q are two
    distinct Hamiltonian paths.  Uses: #paths >= 2  iff  G - e has a Hamiltonian path for some edge e of P."""
    st, P, _ = forcing_ham(n, edges, node_limit=node_limit)
    if st == 'UNKNOWN':
        raise RuntimeError('undecided')
    if st != 'PATH':
        return 'NONE', None, None
    Es = sorted(edges)
    for a, b in zip(P, P[1:]):
        e = (min(a, b), max(a, b))
        E2 = [f for f in Es if f != e]
        st2, Q, _ = forcing_ham(n, E2, node_limit=node_limit)
        if st2 == 'UNKNOWN':
            raise RuntimeError('undecided')
        if st2 == 'PATH':
            return 'MULTIPLE', P, Q
    return 'UNIQUE', P, None
