# Second order on random regular tournaments (same 3-cycle-reversal chain and seeds as two_loop_random_regular.py).
# Usage: second_order_random_regular.py n count seed HK.  mac-mini-2026-10-22-sierp.
import sys, random, math, subprocess
import numpy as np
sys.path.insert(0, '.')
from second_order import predict2
n, cnt, seed, HK = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
rng = random.Random(seed); d = (n - 1) // 2
A = [[1 if 0 < (j - i) % n <= d else 0 for j in range(n)] for i in range(n)]
def step():
    while True:
        u = rng.randrange(n); v = rng.randrange(n); w = rng.randrange(n)
        if len({u, v, w}) == 3 and A[u][v] and A[v][w] and A[w][u]:
            A[u][v] = A[v][w] = A[w][u] = 0; A[v][u] = A[w][v] = A[u][w] = 1; return
codes, mats = [], []
for s in range(cnt):
    for _ in range(30 * n * n): step()
    codes.append(''.join('1' if A[i][j] else '0' for i in range(n) for j in range(i + 1, n)))
    mats.append(np.array(A, dtype=float))
Hs = [int(x) for x in subprocess.run([HK], input='\n'.join(codes) + '\n', capture_output=True, text=True).stdout.split()]
r1, r2 = [], []
for M, H in zip(mats, Hs):
    rho, first, second, X2 = predict2(M - M.T)
    ex = H * 2 ** (n - 1) / math.factorial(n) / rho - 1
    r1.append(ex - first); r2.append(ex - first - second)
r1, r2 = np.array(r1), np.array(r2)
print(f"n={n}, {cnt} random regular (seed {seed}): residual 1st order mean {r1.mean():+.3e} sd {r1.std():.3e};  2nd order mean {r2.mean():+.3e} sd {r2.std():.3e}")
