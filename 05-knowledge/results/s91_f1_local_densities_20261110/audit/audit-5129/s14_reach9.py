# audit-5129: first odd N >= 69 that is not D + 2^t U (t >= 1) for a rooted tournament of order <= 9
# (THM-4865 (e): every odd number in [69, 382687] is attained, near-join values with k <= 8)
import numpy as np
p9 = np.load("pairs9.npy")
D = (p9 >> 32).astype(np.int64); U = (p9 & 0xFFFFFFFF).astype(np.int64)
LIM = 1 << 20
for kmax in (1, 8, 40):
    hit = np.zeros(LIM // 2, dtype=bool)
    for t in range(1, kmax + 1):
        N = D + (U << t)
        ok = N < LIM
        if not ok.any():
            break
        hit[(N[ok] - 1) // 2] = True
    Ns = np.arange(69, LIM, 2)
    miss = Ns[~hit[(Ns - 1) // 2]]
    print("t <= %d: least odd N >= 69 not attained: %d ; odd N < 69 not attained: %s" %
          (kmax, int(miss[0]), [n for n in range(1, 69, 2) if not hit[(n - 1) // 2]]))
