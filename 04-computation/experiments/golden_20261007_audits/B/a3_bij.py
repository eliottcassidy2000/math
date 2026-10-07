# kappa_n injectivity for 23 <= n <= 29 (words only; points are in bijection with words).
import sys
sys.setrecursionlimit(10000)
from a1_golden import Quot, fibpair, lucas
def run(n):
    g = fibpair(n); Q = Quot((g[0] - 1, g[1]))
    pw = [fibpair(n - 1 - j) for j in range(n)]
    seen = {}
    coll = []
    def dfs(j, first, prev, m, k, word):
        if j == n:
            if prev and first:  # cyclic 11
                return
            c = Q.red((m, k))
            if c in seen: coll.append((seen[c], word))
            else: seen[c] = word
            return
        dfs(j + 1, first if j else 0, 0, m, k, word * 2)
        if not prev:
            p = pw[j]
            dfs(j + 1, first if j else 1, 1, m + p[0], k + p[1], word * 2 + 1)
    dfs(0, 0, 0, 0, 0, 0)
    return len(seen), Q.N, coll
for n in range(23, 30):
    cls, N, coll = run(n)
    desc = [(format(a, f'0{n}b'), format(b, f'0{n}b')) for a, b in coll]
    print(f"n={n} L_n={lucas(n)} |R_n|={N} classes={cls} collisions={desc}")
