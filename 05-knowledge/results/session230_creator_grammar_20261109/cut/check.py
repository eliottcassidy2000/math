def require(condition, message):
    if not condition:
        raise RuntimeError(message)
from itertools import combinations, permutations
import json

def adj(n, bits):
    a = [[0] * n for _ in range(n)]
    for k, (i, j) in enumerate(combinations(range(n), 2)):
        a[i][j] = bits >> k & 1
        a[j][i] = 1 - a[i][j]
    return a

def counts(a):
    n = len(a)
    N = 1 << n
    f = [[0] * n for _ in range(N)]
    g = [[0] * n for _ in range(N)]
    h = [0] * N
    h[0] = 1
    for S in range(1, N):
        for v in range(n):
            if S >> v & 1:
                R = S ^ 1 << v
                f[S][v] = sum((f[R][u] * a[u][v] for u in range(n))) if R else 1
                g[S][v] = sum((g[R][u] * a[v][u] for u in range(n))) if R else 1
        h[S] = sum(f[S])
    return (h, f, g)

def kernel(a, h):
    n = len(a)
    V = (1 << n) - 1
    K = [[0] * n for _ in range(n)]
    for r in range(2, n + 1, 2):
        for P in permutations(range(n), r):
            if all((a[P[i]][P[i + 1]] for i in range(r - 1))):
                S = sum((1 << i for i in P))
                K[P[0]][P[-1]] += h[V ^ S]
    return K
out = []
witnesses = {}
same_h = {}
for n in range(1, 6):
    stats = {'n': n, 'bases': 0, 'attachments': 0, 'mixed_support_derivatives': 0}
    for bits in range(1 << n * (n - 1) // 2):
        a = adj(n, bits)
        h, f, g = counts(a)
        V = (1 << n) - 1
        K = kernel(a, h)
        for u, v in combinations(range(n), 2):
            Wuv = sum((f[S][u] * g[V ^ S][v] for S in range(1, V) if S >> u & 1 and (not S >> v & 1)))
            Wvu = sum((f[S][v] * g[V ^ S][u] for S in range(1, V) if S >> v & 1 and (not S >> u & 1)))
            require(Wuv + Wvu == 2 * (K[u][v] + K[v][u]), 'check.py:37: validation failed')
            require(K[u][v] + K[v][u] >= h[V ^ 1 << u ^ 1 << v], 'check.py:38: validation failed')
        for u in range(n):
            W = sum((f[S][u] * sum(g[V ^ S]) for S in range(1, V) if S >> u & 1))
            require(f[V][u] - g[V][u] + W == 2 * sum((K[v][u] for v in range(n))), 'check.py:41: validation failed')
        for S in range(1 << n):
            cut = sum((K[u][v] for u in range(n) for v in range(n) if S >> u & 1 and (not S >> v & 1)))
            b = [row + [0 if S >> i & 1 else 1] for i, row in enumerate(a)] + [[S >> i & 1 for i in range(n)] + [0]]
            direct = counts(b)[0][-1]
            require(direct == h[V] + 2 * cut, 'check.py:46: validation failed')
            require(direct not in (7, 21), 'check.py:47: validation failed')
            require((cut == 0) == (not any((a[u][v] for u in range(n) for v in range(n) if S >> u & 1 and (not S >> v & 1)))), 'check.py:49: validation failed')
            stats['attachments'] += 1
        for B in range(1 << n):
            rest = V ^ B
            R = rest
            while R:
                z = 0
                Q = R
                while True:
                    z += (-1) ** (R.bit_count() - Q.bit_count()) * h[B | Q]
                    if not Q:
                        break
                    Q = Q - 1 & R
                require(z >= 0, 'check.py:60: validation failed')
                stats['mixed_support_derivatives'] += 1
                R = R - 1 & rest
        key = (n, h[V])
        sig = tuple(sorted((h[V] + 2 * sum((K[u][v] for u in range(n) for v in range(n) if S >> u & 1 and (not S >> v & 1))) for S in range(1 << n))))
        if key in same_h and same_h[key][1] != sig and ('scalar_loss' not in witnesses):
            witnesses['scalar_loss'] = {'n': n, 'h': h[V], 'bits_A': same_h[key][0], 'bits_B': bits, 'extension_values_A': sorted(set(same_h[key][1])), 'extension_values_B': sorted(set(sig))}
        same_h.setdefault(key, (bits, sig))
        stats['bases'] += 1
    out.append(stats)
print(json.dumps({'universe': 'every labelled tournament of orders 1 through 5, every attachment; all disjoint support/base marked sets', 'checks': out, 'hostile_witnesses': witnesses}, indent=2))
