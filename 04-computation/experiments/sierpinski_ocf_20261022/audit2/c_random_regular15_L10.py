# Audit2 (C): the first random regular 15-tournament of c_random_regular15.py at z^10 (set partitions, exact).
import sys, random, time
sys.path.insert(0, '.')
from b_partitions import cluster_sums
from common import rot_A, skew, X0_series, X1_series
random.seed(15)
n = 15
A = rot_A(n).copy()
for step in range(3000):
    a, b, c = random.sample(range(n), 3)
    if A[a, b] and A[b, c] and A[c, a]:
        A[a, b] = A[b, c] = A[c, a] = 0; A[b, a] = A[c, b] = A[a, c] = 1
S = skew(A)
t = time.time()
cs = cluster_sums(S, 10, regular=True, emax=1, odd=False, Lmin=10, verbose=False)
x0 = X0_series(S, 10); x1 = X1_series(S, 10, parts=True)
print(f"random regular 15-tournament #0, z^10: X0 exact {cs[10].get(0, 0) == x0[10]}, X1 exact {cs[10].get(1, 0) == x1[10][3]} "
      f"(X1 = {float(x1[10][3]):+.10e}; F8 {float(x1[10][0]):+.3e}, theta {float(x1[10][1]):+.3e}, dumbbell {float(x1[10][2]):+.1e}); {time.time() - t:.0f}s")
