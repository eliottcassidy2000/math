"""Word-singletons M_K(N) = 1 for K = 2..8 and N <= X (auditor's check of HYP-9494 (C)); standard library only.
M_K(N) = [N=1] + #{wJ : #J(w) <= K-1, F(t.w) = N} + #{wZ : #J(w) <= K, F(t.w) = N}  (THM-4904 section 0).
The DFS runs over all prefixes with at most 8 letters J and F(t.w) <= X (F increases strictly, so the pruning is exact).
usage: python3 -I k_singletons.py [X]   (default X = 2048)
"""
import sys
from collections import defaultdict

sys.setrecursionlimit(100000)
X = int(sys.argv[1]) if len(sys.argv) > 1 else 2048
KMAX = 8
M = {K: defaultdict(int) for K in range(KMAX + 1)}
pref = defaultdict(list)          # prefixes with >= 6 letters J, by (F, #J)


def dfs(D, E, U, d, j, w):
    F = D - E + U
    if F > X:
        return
    for K in range(KMAX + 1):
        if j <= K:
            M[K][F] += 1          # the word wZ
        if j + 1 <= K:
            M[K][F] += 1          # the word wJ
    if j >= 6:
        pref[(F, j)].append(w)
    dfs(D, 0, 2 * U, d, j, w + 'O')
    dfs(F, d, 2 * U - d, d, j, w + 'Z')
    if j < KMAX:
        dfs(F, D - E, D + U, D, j + 1, w + 'J')


def dagger(w):
    return ''.join({'O': 'Z', 'Z': 'O', 'J': 'J'}[c] for c in reversed(w))


dfs(1, 0, 2, 1, 0, '')
for K in range(2, KMAX + 1):
    s = [N for N in range(3, X + 1, 2) if M[K][N] == 1]
    print('K=%d: odd 1 < N <= %d with M_K(N) = 1: %s' % (K, X, s))
for N in (1877, 1897):
    if N <= X:
        w7 = pref[(N, 7)]
        print('N=%d: M_6 = %d, M_7 = %d, M_8 = %d; prefixes with exactly 7 J: %s; dagger-fixed among them: %s'
              % (N, M[6][N], M[7][N], M[8][N], w7, [w for w in w7 if w == dagger(w)]))
