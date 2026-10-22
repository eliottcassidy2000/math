# Two-loop law on random regular tournaments (non-vertex-transitive, larger n).  Sampler: Markov chain of directed
# 3-cycle reversals from R_n (reversals preserve scores and connect all tournaments with one score sequence), 30 n^2
# steps between samples.  Exact H by Held-Karp (hk64.c, stdin codes).  Compares y = n(1 - kappa/rho) with the
# parameter-free y_pred of HYP-9365 (two_loop_regular.two_loop, dumbbells included).
# Usage: two_loop_random_regular.py n count seed HK_BINARY.  mac-mini-2026-10-22-sierp.
import sys, random, math, subprocess
import numpy as np
sys.path.insert(0, '.')
from two_loop_regular import two_loop
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
out = subprocess.run([HK], input='\n'.join(codes) + '\n', capture_output=True, text=True).stdout.split()
Hs = [int(x) for x in out]
Y, P, M4 = [], [], []
for M, H in zip(mats, Hs):
    S = M - M.T
    rho, yp, f8, th, db, m4 = two_loop(S)
    kap = H * 2 ** (n - 1) / math.factorial(n)
    Y.append(n * (1 - kap / rho)); P.append(yp); M4.append(m4)
Y, P, M4 = map(np.array, (Y, P, M4))
rk = lambda v: np.argsort(np.argsort(v))
b, a = np.polyfit(M4, Y, 1)
print(f"n={n}: {cnt} random regular tournaments (3-cycle-reversal chain, seed {seed})")
print(f"  y in [{Y.min():.4f}, {Y.max():.4f}], sd {Y.std():.4f};  y_pred in [{P.min():.4f}, {P.max():.4f}], sd {P.std():.4f}")
print(f"  corr(y, y_pred) = {np.corrcoef(Y, P)[0, 1]:+.4f}, spearman = {np.corrcoef(rk(Y), rk(P))[0, 1]:+.4f}, mean(y - y_pred) = {np.mean(Y - P):+.4f}, rms scatter = {np.std(Y - P):.4f}")
print(f"  empirical a + b M4: a = {a:.3f}, b = {b:.3f}, corr = {np.corrcoef(Y, a + b * M4)[0, 1]:+.4f}")
np.savetxt(f'two_loop_random_regular_{n}.tsv', np.column_stack([Y, P, M4]), header='y y_pred M4', fmt='%.8f')
