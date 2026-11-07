"""S84: SAT checks of sink spectra under partial coherence, any n, with optional acyclicity (CEGAR on directed cycles)."""
import itertools, sys, time
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from perm_relax import edges_and_classes


def build(n, level):
    perms, cls, lits = edges_and_classes(n, level)
    index = {p: k for k, p in enumerate(perms)}
    E = []
    for u, p in enumerate(perms):
        for i in range(n - 1):
            q = list(p); q[i], q[i + 1] = q[i + 1], q[i]
            v = index[tuple(q)]
            if u < v:
                c, b = lits[u][i]
                E.append((u, v, c, b))   # if bit(c) == b the edge points into u, else into v
    return perms, cls, lits, E


def find_cycle(nv, E, model_bits):
    adj = [[] for _ in range(nv)]
    for (u, v, c, b) in E:
        if model_bits[c] == b:
            adj[v].append((u, (u, v, c, b)))   # v -> u
        else:
            adj[u].append((v, (u, v, c, b)))   # u -> v
    color = [0] * nv
    parent = [None] * nv
    for s in range(nv):
        if color[s]:
            continue
        stack = [(s, iter(adj[s]))]
        color[s] = 1
        while stack:
            x, it = stack[-1]
            nxt = next(it, None)
            if nxt is None:
                color[x] = 2
                stack.pop()
                continue
            y, e = nxt
            if color[y] == 0:
                color[y] = 1
                parent[y] = (x, e)
                stack.append((y, iter(adj[y])))
            elif color[y] == 1:
                # cycle y -> ... -> x -> y
                cyc = [e]
                z = x
                while z != y:
                    pz, pe = parent[z]
                    cyc.append(pe)
                    z = pz
                return cyc
    return None


def solve(n, level, k, acyclic=False, mode='eq', max_iter=20000):
    from pysat.formula import CNF, IDPool
    from pysat.card import CardEnc, EncType
    from pysat.solvers import Cadical153
    perms, cls, lits, E = build(n, level)
    nc = len(cls)
    pool = IDPool(start_from=nc + 1)
    cnf = CNF()
    svars = []
    for u, L in enumerate(lits):
        sv = pool.id(('s', u))
        svars.append(sv)
        ll = [(c + 1) if b == 1 else -(c + 1) for c, b in L]
        for l in ll:
            cnf.append([-sv, l])
        cnf.append([sv] + [-l for l in ll])
    if mode == 'eq':
        card = CardEnc.equals(lits=svars, bound=k, vpool=pool, encoding=EncType.totalizer)
    else:
        card = CardEnc.atleast(lits=svars, bound=k, vpool=pool, encoding=EncType.totalizer)
    cnf.extend(card.clauses)
    S = Cadical153(bootstrap_with=cnf.clauses)
    it = 0
    while True:
        it += 1
        if not S.solve():
            S.delete()
            return False, it, None
        model = S.get_model()
        bits = [1 if model[c] > 0 else 0 for c in range(nc)]
        if not acyclic:
            S.delete()
            return True, it, bits
        cyc = find_cycle(len(perms), E, bits)
        if cyc is None:
            S.delete()
            return True, it, bits
        # forbid this cycle: at least one of its edge-classes must flip
        clause = []
        for (u, v, c, b) in cyc:
            clause.append(-(c + 1) if bits[c] == 1 else (c + 1))
        S.add_clause(list(set(clause)))
        if it > max_iter:
            S.delete()
            return None, it, None


if __name__ == '__main__':
    n = int(sys.argv[1]); level = sys.argv[2]; acyc = sys.argv[3] == 'acyclic'
    ks = [int(x) for x in sys.argv[4].split(',')]
    mode = sys.argv[5] if len(sys.argv) > 5 else 'eq'
    for k in ks:
        t0 = time.time()
        ok, it, bits = solve(n, level, k, acyclic=acyc, mode=mode)
        print('n=%d %s %s %s k=%d : %s (iterations %d, %.1fs)' % (n, level, 'acyclic' if acyc else 'any', mode, k, ok, it, time.time() - t0))
        sys.stdout.flush()
