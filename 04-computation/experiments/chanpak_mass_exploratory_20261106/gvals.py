import sys
K = int(sys.argv[1]) if len(sys.argv) > 1 else 20
X = 1 << K
M = 2 * X + 2
g = [0] * (M + 1)
g[0] = 1
for v in range(1, M + 1):
    if v & 1:
        g[v] = g[v >> 1]
    else:
        u = v >> 1
        g[v] = g[u] + u + (u & 1)
vals = set(x for x in g if x <= X)
for k in range(2, K + 1):
    c = sum(1 for x in vals if x < (1 << k))
    print(k, c, (1 << k) // 2, c / ((1 << k) // 2))
miss = [n for n in range(1, 200, 2) if n not in vals]
print('missing odd < 200:', miss)
# multiplicity: number of even v (or v=0) with g(v)=n, n small
from collections import Counter
cnt = Counter(g[v] for v in range(0, M + 1, 2) if g[v] <= 200)
print('reps (even v) for n<=61:', [(n, cnt[n]) for n in range(1, 62, 2)])
