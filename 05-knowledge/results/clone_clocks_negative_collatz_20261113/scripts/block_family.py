"""Block family a, v, c, u, D(p) -> a, X(x), F(q) -> c, Y(y) (transitive order otherwise): two-module tables at
(u, v) and the digits (sigma, sigma_2).  Reproduces the table quoted in CLOCKS.md section 3."""
import sys; sys.path.insert(0, 'scripts')
from twomod import visit_table
def build(order, back):
    n = len(order); idx = {nm: i for i, nm in enumerate(order)}; adj = [0]*n
    for i in range(n):
        for j in range(i+1, n):
            if (order[j], order[i]) in back: adj[j] |= 1 << i
            else: adj[i] |= 1 << j
    return adj, idx
if __name__ == '__main__':
    for p in (1, 2, 3):
        for q in (1, 2, 3):
            for x in (0, 1, 2):
                order = ['a', 'v', 'c', 'u'] + ['d%d' % i for i in range(p)] + ['x%d' % i for i in range(x)] + ['f%d' % i for i in range(q)]
                back = set([('d%d' % i, 'a') for i in range(p)] + [('f%d' % i, 'c') for i in range(q)])
                adj, idx = build(order, back)
                tab = visit_table(adj, len(order), idx['u'], idx['v'], CU=9, CV=9)
                g = lambda i, ip: tab.get((i, ip), 0)
                high = any(c for (i, ip), c in tab.items() if i >= 1 and ip >= 1 and (i >= 4 or ip >= 3))
                A, al = g(3, 1), g(3, 2)
                tern = (not high) and A > 0 and al > 0 and g(2, 1) == 3*A and g(2, 2) == 3*al
                print((p, q, x), 'two-clock' if tern else '-', 'X1=', (g(1, 1), g(2, 1), g(3, 1)), 'X2=', (g(1, 2), g(2, 2), g(3, 2)),
                      ('sigma=%d sig2=%d' % (3*A - g(1, 1), 3*al - g(1, 2))) if tern else '')
