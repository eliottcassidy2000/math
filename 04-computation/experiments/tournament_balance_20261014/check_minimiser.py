# Independent exact check of a candidate minimiser (subset DP in Python, Fractions) + structure: delta, most balanced pairs,
# cut vertices, minimum feedback arc set (exact DP over orderings), vertex-wise constant.
import sys
from fractions import Fraction as Fr
from cutvertex_chain import brute, components
def dec(s):
    n = int((1 + (1 + 8 * len(s)) ** 0.5) / 2); A = [[False] * n for _ in range(n)]; k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if s[k] == '1': A[i][j] = True
            else: A[j][i] = True
            k += 1
    return A
def fas(A):
    n = len(A); N = 1 << n; INF = 10 ** 9; dp = [INF] * N; dp[0] = 0
    for S in range(N):
        if dp[S] == INF: continue
        for v in range(n):
            if not S >> v & 1:
                c = sum(1 for u in range(n) if S >> u & 1 and A[v][u])
                if dp[S] + c < dp[S | 1 << v]: dp[S | 1 << v] = dp[S] + c
    return dp[N - 1]
if __name__ == '__main__':
  for s in sys.argv[1:]:
      A = dec(s); n = len(A); H, cnt = brute(A)
      d = max(min(Fr(cnt[x][y], H), 1 - Fr(cnt[x][y], H)) for x in range(n) for y in range(x + 1, n))
      vw = min(max(min(Fr(cnt[v][w], H), 1 - Fr(cnt[v][w], H)) for w in range(n) if w != v) for v in range(n))
      cuts = [z for z in range(n) if len(components(A, [x for x in range(n) if x != z])) > 1]
      print(f"n={n} H={H} delta={d} = {float(d):.6f}  vertex-wise={vw}  scores={sorted(sum(r) for r in A)}  cut vertices={cuts}  fas={fas(A) if n <= 14 else '-'}")
