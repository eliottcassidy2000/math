"""S84: acyclicity within the coherence levels at n = 4 (exhaustive), and SAT spectra at n = 5."""
import itertools, sys, time
import numpy as np
from numba import njit
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from perm_relax import edges_and_classes


def edge_list(n, level):
    perms, cls, lits = edges_and_classes(n, level)
    index = {p: k for k, p in enumerate(perms)}
    E = []  # (u, v, class, bit): edge u--v; if class bit == bit then edge points INTO u
    for u, p in enumerate(perms):
        for i in range(n - 1):
            q = list(p); q[i], q[i + 1] = q[i + 1], q[i]
            v = index[tuple(q)]
            if u < v:
                c, bit = lits[u][i]
                E.append((u, v, c, bit))
    return perms, cls, lits, np.array(E, dtype=np.int64)


@njit(cache=True)
def acyclic_and_sinks(masks, E, nv, lit_c, lit_b):
    out_ac = np.zeros(len(masks), dtype=np.bool_)
    out_s = np.zeros(len(masks), dtype=np.int16)
    m = E.shape[0]
    indeg = np.zeros(nv, dtype=np.int64)
    head = np.zeros(m, dtype=np.int64)
    tail = np.zeros(m, dtype=np.int64)
    stack = np.zeros(nv, dtype=np.int64)
    for t in range(len(masks)):
        M = masks[t]
        for k in range(nv):
            indeg[k] = 0
        for e in range(m):
            u = E[e, 0]; v = E[e, 1]; c = E[e, 2]; b = E[e, 3]
            if ((M >> c) & 1) == b:   # points into u
                tail[e] = v; head[e] = u
            else:
                tail[e] = u; head[e] = v
            indeg[head[e]] += 1
        # sinks: vertices with out-degree 0 <=> every incident edge points in; count via literals
        s = 0
        for u in range(nv):
            ok = True
            for j in range(lit_c.shape[1]):
                if ((M >> lit_c[u, j]) & 1) != lit_b[u, j]:
                    ok = False
                    break
            if ok:
                s += 1
        out_s[t] = s
        # Kahn
        top = 0
        for k in range(nv):
            if indeg[k] == 0:
                stack[top] = k; top += 1
        seen = 0
        while top > 0:
            top -= 1
            x = stack[top]
            seen += 1
            for e in range(m):
                if tail[e] == x:
                    h = head[e]
                    indeg[h] -= 1
                    if indeg[h] == 0:
                        stack[top] = h; top += 1
        out_ac[t] = (seen == nv)
    return out_ac, out_s


def n4_acyclic():
    for level in ('L0', 'L2', 'L1'):
        perms, cls, lits, E = edge_list(4, level)
        nc = len(cls)
        lit_c = np.array([[c for c, b in L] for L in lits], dtype=np.int64)
        lit_b = np.array([[b for c, b in L] for L in lits], dtype=np.int64)
        masks = np.arange(1 << nc, dtype=np.int64)
        t0 = time.time()
        ac, s = acyclic_and_sinks(masks, E, len(perms), lit_c, lit_b)
        allspec = dict(zip(*[x.tolist() for x in np.unique(s, return_counts=True)]))
        acspec = dict(zip(*[x.tolist() for x in np.unique(s[ac], return_counts=True)]))
        print('n=4 %s classes=%d  all: %s' % (level, nc, allspec))
        print('            acyclic (%d of %d): %s   [%.1fs]' % (int(ac.sum()), len(masks), acspec, time.time() - t0))


def sat_spectrum(n, level, ks):
    from pysat.formula import CNF, IDPool
    from pysat.card import CardEnc, EncType
    from pysat.solvers import Cadical153
    perms, cls, lits = edges_and_classes(n, level)
    nc = len(cls)
    res = {}
    for k in ks:
        pool = IDPool(start_from=nc + 1)
        cnf = CNF()
        svars = []
        for u, L in enumerate(lits):
            sv = pool.id(('s', u))
            svars.append(sv)
            litlist = [(c + 1) if b == 1 else -(c + 1) for c, b in L]
            for l in litlist:
                cnf.append([-sv, l])
            cnf.append([sv] + [-l for l in litlist])
        card = CardEnc.equals(lits=svars, bound=k, vpool=pool, encoding=EncType.totalizer)
        cnf.extend(card.clauses)
        with Cadical153(bootstrap_with=cnf.clauses) as S:
            ok = S.solve()
            res[k] = ok
    return nc, res


if __name__ == '__main__':
    if sys.argv[1] == 'n4':
        n4_acyclic()
    else:
        level = sys.argv[2]
        ks = list(range(int(sys.argv[3]), int(sys.argv[4]) + 1))
        t0 = time.time()
        nc, res = sat_spectrum(5, level, ks)
        print('n=5 %s classes=%d  achievable: %s' % (level, nc, [k for k in ks if res[k]]))
        print('            not achievable: %s   [%.1fs]' % ([k for k in ks if not res[k]], time.time() - t0))
