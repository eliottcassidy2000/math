"""procgen_sumgraph_20260926_solver.py -- an independent exact Hamiltonian-path solver for sum graphs,
with conflict cores.

Session collatz-procgen-20260922, lane "sumgraph" (2026-09-26).  Written from scratch; it does NOT
use the smallgraph library's solvers (HCForcing / ham.c / CP-SAT), so it is a second code path.

Sum graph G_S(n): vertices 1..n, edge {x,y} (x != y) iff x+y in S.  Every edge {x, t-x} lies in the
reflection matching M_t of r_t(x) = t - x.

Method (for Hamiltonian PATHS with prescribed end set):
  need[v] = 1 for a prescribed end, 2 otherwise.  Unit rules
    R1  a vertex whose non-deleted edges number need[v]  -> all of them are path edges;
    R2  a vertex with need[v] path edges                 -> its other edges are deleted;
    R4  the edge joining the two ends of a forced segment is deleted (no cycles);
  contradictions: fewer than need[v] non-deleted edges (LOW), more than need[v] forced (OVER),
  a forced cycle (CYC), a forced segment joining the two prescribed ends that misses a vertex (SHORT).
  Search: DPLL branching (force / delete one edge at a most-constrained vertex), exhaustive, with a
  connectivity cut every 16th node.

Unknown second end (only one vertex of degree <= 1, e.g. the largest power of two for the
powers-of-2-or-3 graph C_n): a BASE run with need = 2 at every other vertex must fail (degree-sum
parity).  Every contradiction it reaches has a CORE: the set of vertices whose need value was used in
its derivation.  If the second end e is outside the core, the same derivation is valid with
need[e] = 1, so e is impossible.  Hence only core vertices are candidates (Lemma K of the note).
When the base run reaches no contradiction the solver falls back to all candidates.

Everything is deterministic.  Every path returned is re-verified by `verify_path`.
"""
import sys

sys.setrecursionlimit(1000000)


def pow23(lim):
    s = set()
    v = 1
    while v <= lim:
        s.add(v)
        v *= 2
    v = 1
    while v <= lim:
        s.add(v)
        v *= 3
    return sorted(s)


def squares(lim):
    out = []
    k = 1
    while k * k <= lim:
        out.append(k * k)
        k += 1
    return out


def targets_for(n, family):
    lim = 2 * n - 1
    T = pow23(lim) if family == 'pow23' else squares(lim)
    return [t for t in T if 3 <= t <= lim]


def edges_of(n, T):
    E = []
    for t in T:
        for x in range(max(1, t - n), (t - 1) // 2 + 1):
            y = t - x
            if x < y <= n:
                E.append((x, y, t))
    return E


def verify_path(n, seq, T):
    Ts = set(T)
    if sorted(seq) != list(range(1, n + 1)):
        return False
    return all(seq[i] + seq[i + 1] in Ts for i in range(n - 1))


class Contra(Exception):
    pass


AV, FO, DE = 0, 1, 2


class State:
    def __init__(self, n, E, ends, track=False):
        self.n = n
        m = len(E)
        self.eu = [e[0] for e in E]
        self.ev = [e[1] for e in E]
        self.st = [AV] * m
        self.tm = [0] * m
        self.rule = [None] * m
        self.cause = [0] * m
        self.adj = [[] for _ in range(n + 1)]
        for i, (x, y, t) in enumerate(E):
            self.adj[x].append(i)
            self.adj[y].append(i)
        self.navail = [len(a) for a in self.adj]
        self.nforced = [0] * (n + 1)
        self.need = [2] * (n + 1)
        self.need[0] = 0
        self.isend = [False] * (n + 1)
        for e in ends:
            self.need[e] = 1
            self.isend[e] = True
        self.nends = len(ends)
        self.other = list(range(n + 1))
        self.size = [1] * (n + 1)
        self.trail = []
        self.clock = 0
        self.track = track
        self.contra = None

    def opp(self, i, v):
        return self.ev[i] if self.eu[i] == v else self.eu[i]

    # ---------------- primitive changes (trailed) ----------------
    def _delete(self, i, rule, cause, q):
        if self.st[i] != AV:
            if self.st[i] == FO:
                self._fail(('DELFORCED', i, rule, cause))
            return
        self.st[i] = DE
        self.clock += 1
        self.tm[i] = self.clock
        self.rule[i] = rule
        self.cause[i] = cause
        u, v = self.eu[i], self.ev[i]
        self.navail[u] -= 1
        self.navail[v] -= 1
        self.trail.append(('s', i))
        q.append(u)
        q.append(v)

    def _force(self, i, rule, cause, q):
        if self.st[i] == FO:
            return
        if self.st[i] == DE:
            self._fail(('FORCEDEL', i, rule, cause))
        u, v = self.eu[i], self.ev[i]
        # record the attempted decision first so that cores can see it
        self.clock += 1
        t_now = self.clock
        if self.nforced[u] >= self.need[u]:
            self._fail(('OVER', u, i, rule, cause, t_now))
        if self.nforced[v] >= self.need[v]:
            self._fail(('OVER', v, i, rule, cause, t_now))
        a, b = self.other[u], self.other[v]
        if a == v:
            self._fail(('CYC', u, v, i, rule, cause, t_now))
        self.st[i] = FO
        self.tm[i] = t_now
        self.rule[i] = rule
        self.cause[i] = cause
        self.nforced[u] += 1
        self.nforced[v] += 1
        ns = self.size[u] + self.size[v]
        self.trail.append(('f', i, a, b, self.other[a], self.other[b], self.size[a], self.size[b]))
        self.other[a], self.other[b] = b, a
        self.size[a] = self.size[b] = ns
        q.append(u)
        q.append(v)
        q.append(a)
        q.append(b)
        # SHORT: a segment from a prescribed end to a prescribed end, both saturated, not spanning
        if (ns < self.n and self.isend[a] and self.isend[b] and self.nforced[a] == self.need[a]
                and self.nforced[b] == self.need[b]):
            self._fail(('SHORT', a, b, t_now))
        # R4: never close a cycle
        if ns < self.n or True:
            for j in self.adj[a]:
                if self.st[j] == AV and self.opp(j, a) == b:
                    self._delete(j, 'R4', (a, b, t_now), q)
                    break

    def _fail(self, info):
        self.contra = info
        raise Contra(info)

    def undo(self, mark):
        while len(self.trail) > mark:
            op = self.trail.pop()
            if op[0] == 's':
                i = op[1]
                self.st[i] = AV
                self.navail[self.eu[i]] += 1
                self.navail[self.ev[i]] += 1
            else:
                _, i, a, b, oa, ob, sa, sb = op
                self.st[i] = AV
                self.nforced[self.eu[i]] -= 1
                self.nforced[self.ev[i]] -= 1
                self.other[a], self.other[b] = oa, ob
                self.size[a], self.size[b] = sa, sb

    def propagate(self, q):
        adj, st = self.adj, self.st
        while q:
            w = q.pop()
            k = self.navail[w]
            f = self.nforced[w]
            nd = self.need[w]
            if k < nd:
                self._fail(('LOW', w, self.clock))
            if f > nd:
                self._fail(('OVERSTATE', w, self.clock))
            if f == nd and k > f:
                for j in adj[w]:
                    if st[j] == AV:
                        self._delete(j, 'R2', w, q)
            elif k == nd and f < nd:
                for j in adj[w]:
                    if st[j] == AV:
                        self._force(j, 'R1', w, q)

    def complete(self):
        return all(self.nforced[v] == self.need[v] for v in range(1, self.n + 1))

    def extract_path(self):
        n = self.n
        fadj = [[] for _ in range(n + 1)]
        for i in range(len(self.st)):
            if self.st[i] == FO:
                fadj[self.eu[i]].append(self.ev[i])
                fadj[self.ev[i]].append(self.eu[i])
        start = next(v for v in range(1, n + 1) if len(fadj[v]) == 1) if n > 1 else 1
        seq = [start]
        prev = 0
        cur = start
        while True:
            nxt = [w for w in fadj[cur] if w != prev]
            if not nxt:
                break
            prev, cur = cur, nxt[0]
            seq.append(cur)
            if len(seq) > n:
                break
        return seq

    # ---------------- cores ----------------
    def core(self):
        """Vertices whose need value is used by the derivation of self.contra (one run, no undo)."""
        info = self.contra
        core = set()
        todo = []   # edge ids whose decisions must be explained
        seen = set()
        kind = info[0]

        def edges_at(w, state, before):
            return [j for j in self.adj[w] if self.st[j] == state and self.tm[j] < before]

        def seg_edges(a, b, before):
            # walk forced edges with time < before from a to b
            out = []
            prev, cur = 0, a
            steps = 0
            while cur != b and steps <= self.n:
                nxt = None
                for j in self.adj[cur]:
                    if self.st[j] == FO and self.tm[j] < before:
                        w = self.opp(j, cur)
                        if w != prev:
                            nxt = (j, w)
                            break
                if nxt is None:
                    break
                out.append(nxt[0])
                prev, cur = cur, nxt[1]
                steps += 1
            return out

        def attempted(rule, cause, t):
            # premises of an attempted (not recorded) R1 force at vertex `cause`
            if rule == 'R1':
                core.add(cause)
                todo.extend(edges_at(cause, DE, t))
            elif rule == 'BR':
                pass

        if kind == 'LOW':
            w, t = info[1], info[2] + 1
            core.add(w)
            todo.extend(edges_at(w, DE, t))
        elif kind == 'OVERSTATE':
            w, t = info[1], info[2] + 1
            core.add(w)
            todo.extend(edges_at(w, FO, t))
        elif kind == 'OVER':
            w, i, rule, cause, t = info[1:]
            core.add(w)
            todo.extend(edges_at(w, FO, t))
            attempted(rule, cause, t)
        elif kind == 'CYC':
            u, v, i, rule, cause, t = info[1:]
            todo.extend(seg_edges(u, v, t))
            attempted(rule, cause, t)
        elif kind == 'SHORT':
            a, b, t = info[1:]
            core.add(a)
            core.add(b)
            todo.extend(seg_edges(a, b, t + 1))
        else:
            return None
        while todo:
            j = todo.pop()
            if j in seen:
                continue
            seen.add(j)
            r, c, t = self.rule[j], self.cause[j], self.tm[j]
            if r == 'R1':
                core.add(c)
                todo.extend(edges_at(c, DE, t))
            elif r == 'R2':
                core.add(c)
                todo.extend(edges_at(c, FO, t))
            elif r == 'R4':
                a, b, t0 = c
                todo.extend(seg_edges(a, b, t0 + 1))
            else:
                return None   # branching decision: no pure-propagation core
        return core


def connected(S):
    """All vertices joined by non-deleted edges (necessary for a Hamiltonian path)."""
    n = S.n
    seen = bytearray(n + 1)
    seen[1] = 1
    stack = [1]
    cnt = 1
    st, eu, ev, adj = S.st, S.eu, S.ev, S.adj
    while stack:
        v = stack.pop()
        for j in adj[v]:
            if st[j] != DE:
                w = ev[j] if eu[j] == v else eu[j]
                if not seen[w]:
                    seen[w] = 1
                    cnt += 1
                    stack.append(w)
    return cnt == n


PICKS = {
    # which available edge to branch on at the chosen vertex bv
    'mindeg': lambda S, bv, c: min(c, key=lambda j: (S.navail[S.opp(j, bv)], S.opp(j, bv))),
    'first': lambda S, bv, c: c[0],
    'large': lambda S, bv, c: max(c, key=lambda j: S.opp(j, bv)),
}


def _dpll(S, node_limit, counter, pick=PICKS['mindeg']):
    counter[0] += 1
    if node_limit and counter[0] > node_limit:
        return 'ABORT'
    if S.complete():
        return 'FOUND'
    if counter[0] % 16 == 1 and not connected(S):   # connectivity cut (sound; checked every 16th node for speed)
        return 'NONE'
    best, bv = None, 0
    for v in range(1, S.n + 1):
        if S.nforced[v] < S.need[v]:
            key = (S.navail[v] - S.nforced[v], S.navail[v])
            if best is None or key < best:
                best, bv = key, v
    cands = [j for j in S.adj[bv] if S.st[j] == AV]
    j = pick(S, bv, cands)
    for choice in (0, 1):
        mark = len(S.trail)
        try:
            q = []
            if choice == 0:
                S._force(j, 'BR', bv, q)
            else:
                S._delete(j, 'BR', bv, q)
            S.propagate(q)
            r = _dpll(S, node_limit, counter, pick)
            if r in ('FOUND', 'ABORT'):
                return r
        except Contra:
            pass
        S.undo(mark)
    return 'NONE'


_ORDER = ['mindeg', 'first', 'large']   # adaptive: the rule that succeeded last is tried first


def portfolio(S, node_limit, counter, quick_only=False):
    """Try the deterministic branching rules with small budgets (the most recently successful rule first), then
    randomized restarts (random edge at the most constrained vertex; fixed seeds 0..7, so runs are reproducible);
    unless quick_only, the last run is exhaustive up to node_limit.  NONE from any complete run is final (each run
    searches the same state exhaustively), so the verdict does not depend on the order; only the path found and the
    node count may."""
    import random
    budgets = [(name, 1500, None) for name in _ORDER] + [('rand', 1500, seed) for seed in range(8)]
    if not quick_only:
        budgets.append((_ORDER[0], node_limit, None))
    for name, lim, seed in budgets:
        if name == 'rand':
            rng = random.Random(seed)
            pick = (lambda r: (lambda S_, bv, c: r.choice(c)))(rng)
        else:
            pick = PICKS[name]
        c = [0]
        mark = len(S.trail)
        r = _dpll(S, lim, c, pick)
        counter[0] += c[0]
        if r == 'FOUND':
            if name != 'rand':
                _ORDER.remove(name)
                _ORDER.insert(0, name)
            return r
        if r == 'NONE':
            return r
        S.undo(mark)
    return 'ABORT'


def solve_with_ends(n, E, ends, node_limit=0):
    """Exact: ('PATH', seq, nodes) / ('NONE', None, nodes) / ('ABORT', None, nodes)."""
    S = State(n, E, ends)
    try:
        S.propagate(list(range(1, n + 1)))
    except Contra:
        return 'NONE', None, 0
    counter = [0]
    r = portfolio(S, node_limit, counter)
    if r == 'FOUND':
        return 'PATH', S.extract_path(), counter[0]
    return r, None, counter[0]


def degrees(n, E):
    d = [0] * (n + 1)
    for x, y, t in E:
        d[x] += 1
        d[y] += 1
    return d


_LAST_FREE_END = []


def _remember(ends, leaves):
    free = [e for e in ends if e not in leaves]
    if free:
        e = free[0]
        _LAST_FREE_END[:] = [e, e + 1, e - 1]


def decide(n, family='pow23', node_limit=200000, want_detail=False):
    """Hamiltonian path of G_S(n)?  Returns dict with keys status in {PATH, NONE, LOCAL, ABORT},
    and details (ends, leaves, candidates, core, survivors, nodes)."""
    T = targets_for(n, family)
    E = edges_of(n, T)
    out = {'n': n}
    if n == 1:
        out.update(status='PATH', seq=[1])
        return out
    d = degrees(n, E)
    leaves = [v for v in range(1, n + 1) if d[v] <= 1]
    out['leaves'] = leaves
    if any(d[v] == 0 for v in range(1, n + 1)) or len(leaves) > 2:
        out['status'] = 'LOCAL'
        return out
    if len(leaves) == 2:
        cand_sets = [tuple(leaves)]
        out['mode'] = 'two-leaf'
    else:
        # base run: prescribed ends = the leaves only (0 or 1 of them), everything else need 2
        B = State(n, E, leaves)
        base_core = None
        try:
            B.propagate(list(range(1, n + 1)))
        except Contra:
            base_core = B.core()
        out['base_core'] = sorted(base_core) if base_core is not None else None
        if len(leaves) == 1:
            p = leaves[0]
            pool = sorted(base_core - {p}) if base_core is not None else [v for v in range(1, n + 1) if v != p]
            # heuristic order only (the candidate set is unchanged): free ends that worked at the previous call first
            hint = [e for e in _LAST_FREE_END if e in set(pool)]
            pool = hint + [e for e in pool if e not in hint]
            cand_sets = [(p, e) for e in pool]
        else:
            # no leaf: both ends free (never happens for C_n, n >= 2: Theorem C1)
            pool = sorted(base_core) if base_core is not None else list(range(1, n + 1))
            cand_sets = [(a, b) for a in pool for b in range(1, n + 1) if b != a]
        out['mode'] = 'one-leaf'
    survivors = []
    total_nodes = 0
    pending = []
    # pass 1: propagation, then a quick search on every surviving candidate (first success wins)
    for ends in cand_sets:
        S = State(n, E, ends)
        try:
            S.propagate(list(range(1, n + 1)))
        except Contra:
            continue
        survivors.append(ends)
        counter = [0]
        r = portfolio(S, node_limit, counter, quick_only=True)
        total_nodes += counter[0]
        if r == 'FOUND':
            seq = S.extract_path()
            assert verify_path(n, seq, T), 'bad path'
            out.update(status='PATH', seq=seq, nodes=total_nodes, ends=ends, survivors=survivors)
            _remember(ends, leaves)
            return out
        if r == 'ABORT':
            pending.append(ends)
    # pass 2: exhaustive search on the candidates not yet decided
    for ends in pending:
        S = State(n, E, ends)
        S.propagate(list(range(1, n + 1)))
        counter = [0]
        r = portfolio(S, node_limit, counter)
        total_nodes += counter[0]
        if r == 'FOUND':
            seq = S.extract_path()
            assert verify_path(n, seq, T), 'bad path'
            out.update(status='PATH', seq=seq, nodes=total_nodes, ends=ends, survivors=survivors)
            _remember(ends, leaves)
            return out
        if r == 'ABORT':
            out.update(status='ABORT', nodes=total_nodes, ends=ends, survivors=survivors)
            return out
    out['survivors'] = survivors
    out.update(status='NONE', nodes=total_nodes)
    return out


if __name__ == '__main__':
    import time
    fam = 'pow23'
    for a in sys.argv[1:]:
        if a in ('pow23', 'squares'):
            fam = a
            continue
        n = int(a)
        t0 = time.time()
        r = decide(n, fam)
        print(n, r['status'], 'leaves', r.get('leaves'), 'core', r.get('base_core'),
              'surv', len(r.get('survivors', [])), 'nodes', r.get('nodes'), round(time.time() - t0, 3))
