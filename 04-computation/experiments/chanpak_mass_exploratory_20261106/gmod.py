import numpy as np
def gtable(n):
    g = np.zeros(n, dtype=np.int64); g[0] = 1
    for v in range(1, n):
        g[v] = g[v >> 1] if v & 1 else g[v >> 1] + (v >> 1) + ((v >> 1) & 1)
    return g
K = 22
G = gtable(1 << (K + 1))
for m in (3, 7, 9, 49, 5):
    print(f"g(v) mod {m}, v in [2^k, 2^(k+1)) (fraction x m; uniform = 1):")
    for k in (10, 14, 18, 22):
        V = 1 << k
        c = np.bincount(G[V:2*V] % m, minlength=m) / V * m
        print(f"  k={k}: " + " ".join(f"{x:.3f}" for x in c[:min(m, 14)]) + (" ..." if m > 14 else ""))
