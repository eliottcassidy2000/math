# kappa(q) = chi(Cay(F_{q^2}, mu_{q+1})) for prime q.  F_{q^2} = F_q[t]/(t^2 - r), r a non-residue.
import sys, time
from pysat.solvers import Cadical153
from pysat.card import CardEnc
def field(q):
    r = next(x for x in range(2, q) if pow(x, (q-1)//2, q) == q-1)
    els = [(a, b) for a in range(q) for b in range(q)]
    mul = lambda x, y: ((x[0]*y[0] + r*x[1]*y[1]) % q, (x[0]*y[1] + x[1]*y[0]) % q)
    return els, mul, r
def kappa_graph(q):
    els, mul, r = field(q)
    one = (1, 0)
    # mu_{q+1}: z with z^(q+1) = 1
    def pw(x, e):
        res = one
        while e:
            if e & 1: res = mul(res, x)
            x = mul(x, x); e >>= 1
        return res
    mu = [z for z in els if z != (0, 0) and pw(z, q+1) == one]
    assert len(mu) == q + 1
    idx = {e: i for i, e in enumerate(els)}
    edges = set()
    for x in els:
        for z in mu:
            y = ((x[0] + z[0]) % q, (x[1] + z[1]) % q)
            i, j = idx[x], idx[y]
            if i < j: edges.add((i, j))
    return len(els), sorted(edges)
def colourable(n, edges, k, tlimit=None):
    s = Cadical153()
    v = lambda i, c: k*i + c + 1
    for i in range(n): s.add_clause([v(i, c) for c in range(k)])
    for (i, j) in edges:
        for c in range(k): s.add_clause([-v(i, c), -v(j, c)])
    # symmetry breaking: vertex 0 colour 0, and its first neighbour colour 1
    s.add_clause([v(0, 0)])
    nb = [j for (i, j) in edges if i == 0]
    if nb and k > 1: s.add_clause([v(nb[0], 1)])
    t = time.time(); r = s.solve(); return r, time.time() - t
for q, ks in ((3, (2, 3)), (7, (3, 4)), (11, (4, 5)), (19, (5,))):
    n, E = kappa_graph(q)
    for k in ks:
        r, dt = colourable(n, E, k)
        print(f"q={q:2d}: |V|={n}, deg={2*len(E)//n}, {k}-colourable: {r}  ({dt:.1f}s)"); sys.stdout.flush()
