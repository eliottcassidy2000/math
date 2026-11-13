"""Constant terms of the tail-law offset c0 and response M over digit-1 clocks, for every out-set O (|O|<=3):
kappa(f) = C when f(j) = a 3^j + C (j = 1,2,3; only ternary-shaped f are recorded)."""
import sys, os, itertools, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tern_check import blowup
from tail_c2_search import clocks
from tail_c2_multi import responses
def kappa(f):
    if f[2] - 3*f[1] != f[1] - 3*f[0]: return None
    return -(f[1] - 3*f[0]) // 2
if __name__ == '__main__':
    C = clocks(); rel = collections.Counter(); single = collections.Counter(); pairsW = collections.Counter()
    for (n, adj, u) in C:
        data = []
        for j in (1, 2, 3):
            A2, n2, block, newid = blowup(adj, n, u, j)
            data.append(responses(A2, n2) + (newid,))
        verts = [y for y in range(n) if y != u]
        for y in verts:
            single[kappa([d[1][d[3][y]] for d in data])] += 1
        for e, f in itertools.combinations(verts, 2):
            pairsW[kappa([d[2][(min(d[3][e], d[3][f]), max(d[3][e], d[3][f]))] for d in data])] += 1
        for size in (1, 2, 3):
            for O in itertools.combinations(verts, size):
                c0 = []; M = []
                for (H, U, W, nid) in data:
                    o = [nid[y] for y in O]
                    c0.append(H - sum(W[(min(a, b), max(a, b))] for a, b in itertools.combinations(o, 2)))
                    M.append(sum(U[y] for y in o))
                kc, km = kappa(c0), kappa(M)
                if kc is not None and km is not None:
                    rel[(size, kc, km)] += 1
    print('kappa(U_y) distribution:', sorted(single.items(), key=lambda t: (t[0] is None, t[0]))[:15])
    print('kappa(W_ef) distribution:', sorted(pairsW.items(), key=lambda t: (t[0] is None, t[0]))[:15])
    for size in (1, 2, 3):
        rows = sorted((kc, km, c) for (s, kc, km), c in rel.items() if s == size)
        print('|O|=%d (kappa(c0), kappa(M), count):' % size, rows[:25])
