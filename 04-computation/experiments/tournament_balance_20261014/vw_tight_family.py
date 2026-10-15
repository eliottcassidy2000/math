# Sharpness of HYP-9279(B): the family W_n = TT_(n-3) => C3 (a -> b -> c -> a) with one back arc b -> 0.
# For each n, compute the vertex-wise constant min_v max_w min(p, 1-p) exactly and the tight vertex.
from fractions import Fraction as Fr
import sys
sys.path.insert(0, '.')
from cutvertex_chain import brute
def W(n):
    A = [[False] * n for _ in range(n)]
    t = n - 3; a, b, c = t, t + 1, t + 2
    for i in range(n):
        for j in range(i + 1, n): A[i][j] = True            # transitive base, everything forward
    A[a][c] = False; A[c][a] = True                          # cyclic triangle a -> b -> c -> a
    A[0][b] = False; A[b][0] = True                          # back arc b -> 0
    return A
for n in range(5, 15):
    A = W(n); H, cnt = brute(A)
    best = {v: max(min(Fr(cnt[v][w], H), 1 - Fr(cnt[v][w], H)) for w in range(n) if w != v) for v in range(n)}
    vmin = min(best.values()); tight = [v for v in best if best[v] == vmin]
    c = n - 1
    print(f"n={n}: H={H}, vertex-wise constant {vmin}, tight vertices {tight} (c = {c}); p(c before a) = {Fr(cnt[c][n-3], H)}, p(c before b) = {Fr(cnt[c][n-2], H)}")
