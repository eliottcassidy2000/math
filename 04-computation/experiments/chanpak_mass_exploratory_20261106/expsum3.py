import numpy as np, math
from fractions import Fraction
def gtable(n):
    g = np.zeros(n, dtype=np.int64); g[0] = 1
    for v in range(1, n):
        g[v] = g[v >> 1] if v & 1 else g[v >> 1] + (v >> 1) + ((v >> 1) & 1)
    return g
K = 21
G = gtable(1 << (K + 1))
print("V | sup over theta = j/L, theta not in {0, 1/2}: value, location, exponent")
for k in range(12, K + 1):
    V = 1 << k; L = 1 << (k + 3)
    hist = np.bincount(G[V:2*V] % L, minlength=L).astype(np.float64)
    S = np.abs(np.fft.fft(hist)) / V
    S[0] = 0; S[L // 2] = 0
    j = int(np.argmax(S)); fr = Fraction(j, L).limit_denominator(500)
    top = np.argsort(S)[-6:][::-1]
    locs = ", ".join(f"{Fraction(int(t), L).limit_denominator(500)}:{S[t]:.3f}" for t in top)
    print(f"2^{k} | {S[j]:.4f} (exponent {math.log(S[j]*V)/math.log(V):.3f}) | top: {locs}")
