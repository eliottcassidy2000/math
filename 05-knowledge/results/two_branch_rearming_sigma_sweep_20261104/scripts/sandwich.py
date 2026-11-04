"""Sandwich construction S(X): vertices Y (= X's vertices) u {o, m, r, a}:  o -> Y, Y -> m, m -> o, o -> r, r -> Y, r -> m,
a -> r, every other vertex -> a.  Then a is a sole-exit auxiliary with root r and sigma = e_r(S - a) = h(X).
Compute H(S(X)) for all X of orders 3..6 and fit H against responses of X: h, K = sum_S h(X[S]) h(X[S^c]) (ordered,
empty allowed), and K3 (three ordered slots)."""
import subprocess, itertools, numpy as np
def decode(code, n):
    adj = [[0]*n for _ in range(n)]; idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if code[idx] == '1': adj[i][j] = 1
            else: adj[j][i] = 1
            idx += 1
    return adj
def hp(adj, verts, start=None, end=None):
    n = len(verts)
    if n == 0: return 1 if (start is None and end is None) else 0
    f = [[0]*n for _ in range(1 << n)]
    for i, v in enumerate(verts):
        if start is None or v == start: f[1 << i][i] = 1
    for S in range(1, 1 << n):
        row = f[S]
        for i in range(n):
            c = row[i]
            if not c: continue
            for j in range(n):
                if adj[verts[i]][verts[j]] and not (S >> j & 1): f[S | 1 << j][j] += c
    full = (1 << n) - 1
    return sum(f[full][i] for i, v in enumerate(verts) if end is None or v == end)
def responses(adj, n):
    V = list(range(n)); h = hp(adj, V)
    K = 0
    for r in range(n + 1):
        for S in itertools.combinations(V, r):
            K += hp(adj, list(S)) * hp(adj, [v for v in V if v not in S])
    # K3: ordered triples of slots
    K3 = 0
    subs = {}
    for r in range(n + 1):
        for S in itertools.combinations(V, r): subs[S] = hp(adj, list(S))
    for S1 in subs:
        rest = [v for v in V if v not in S1]
        for r in range(len(rest) + 1):
            for S2 in itertools.combinations(rest, r):
                S3 = tuple(v for v in rest if v not in S2)
                K3 += subs[S1] * subs[S2] * subs[S3]
    return h, K, K3
def sandwich(adj, n):
    # labels: 0..n-1 = Y, n = o, n+1 = m, n+2 = r, n+3 = a
    N = n + 4; o, m, r, a = n, n + 1, n + 2, n + 3
    A = [[0]*N for _ in range(N)]
    for i in range(n):
        for j in range(n): A[i][j] = adj[i][j]
    for y in range(n): A[o][y] = 1; A[y][m] = 1; A[r][y] = 1
    A[m][o] = 1; A[o][r] = 1; A[r][m] = 1
    A[a][r] = 1
    for v in range(N):
        if v not in (a, r): A[v][a] = 1
    return A, N, r, a
rows = []
for n in range(3, 7):
    for code in subprocess.run(['/opt/homebrew/bin/gentourng', '-q', str(n)], capture_output=True, text=True).stdout.split():
        adj = decode(code, n); h, K, K3 = responses(adj, n)
        A, N, r, a = sandwich(adj, n); V = list(range(N))
        H = hp(A, V); sigma = hp(A, [v for v in V if v != a], end=r)
        rows.append((n, h, K, K3, H, sigma))
import collections
print('rows', len(rows), '; sigma == h(X) always:', all(rw[5] == rw[1] for rw in rows))
# fit H = alpha*h + beta*K + gamma*K3 + delta (least squares over exact data), check exactness
M = np.array([[rw[1], rw[2], rw[3], 1.0] for rw in rows]); y = np.array([rw[4] for rw in rows], dtype=float)
coef, *_ = np.linalg.lstsq(M, y, rcond=None); print('lstsq coefficients (h, K, K3, 1):', np.round(coef, 6))
pred = M @ coef; print('max abs residual:', np.max(np.abs(pred - y)))
# try H = alpha*h + beta*K + delta only
M2 = M[:, [0, 1, 3]]; coef2, *_ = np.linalg.lstsq(M2, y, rcond=None); print('(h, K, 1) fit:', np.round(coef2, 6), 'max residual', np.max(np.abs(M2 @ coef2 - y)))
for rw in rows[:8]: print('  n=%d h=%d K=%d K3=%d -> H=%d sigma=%d' % rw)
