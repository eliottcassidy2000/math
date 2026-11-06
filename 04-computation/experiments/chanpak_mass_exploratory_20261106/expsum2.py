import numpy as np, math
from fractions import Fraction
def gtable(n):
    g = np.zeros(n, dtype=np.int64); g[0] = 1
    for v in range(1, n):
        g[v] = g[v >> 1] if v & 1 else g[v >> 1] + (v >> 1) + ((v >> 1) & 1)
    return g
K = 21
G = gtable(1 << (K + 1))
print("V | sup_{j != 0} |S(j/L)|/V and its location | |S(1/7)| | |S(1/3)| | |S(1/21)| | |S(1/49)|")
for k in range(10, K + 1):
    V = 1 << k; L = 1 << (k + 3)
    hist = np.bincount(G[V:2*V] % L, minlength=L).astype(np.float64)
    S = np.abs(np.fft.fft(hist)) / V
    S[0] = 0
    j = int(np.argmax(S)); fr = Fraction(j, L).limit_denominator(200)
    def at(a, q): return S[int(round(a * L / q)) % L]
    print(f"2^{k} | {S[j]:.4f} at ~{fr} (exponent {math.log(S[j]*V)/math.log(V):.3f}) | {at(1,7):.4f} | {at(1,3):.4f} | {at(1,21):.4f} | {at(1,49):.4f}")
