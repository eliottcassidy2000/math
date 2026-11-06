import numpy as np, math, sys
# S(theta) = (1/V) | sum_{v in [V,2V)} e(theta * g(v)) |, evaluated at theta = j / L on a grid via FFT of the histogram of g.
def gtable(n):
    g = np.zeros(n, dtype=np.int64); g[0] = 1
    for v in range(1, n):
        g[v] = g[v >> 1] if v & 1 else g[v >> 1] + (v >> 1) + ((v >> 1) & 1)
    return g
K = int(sys.argv[1]) if len(sys.argv) > 1 else 18
G = gtable(1 << (K + 1))
from fractions import Fraction
for k in range(12, K + 1):
    V = 1 << k
    vals = G[V:2 * V]
    L = 1 << (k + 3)                       # grid of theta = j / L
    hist = np.bincount(vals % L, minlength=L).astype(np.float64)
    S = np.abs(np.fft.fft(hist)) / V       # S[j] = |sum e(-2 pi i j g / L)| / V
    # major arcs: theta within 1/(q*V^{0.8})... use simple rule: exclude j/L within 4/L of a/q with q <= Q
    Q = 12
    mask = np.ones(L, dtype=bool)
    for q in range(1, Q + 1):
        for a in range(q):
            if math.gcd(a, q) != 1: continue
            c = a * L / q
            lo, hi = int(math.floor(c - 8)), int(math.ceil(c + 8))
            for j in range(lo, hi + 1): mask[j % L] = False
    # also exclude dyadic rationals with denominators up to 2^8 (2-adic major arcs)
    for d in range(1, 9):
        q = 1 << d
        for a in range(1, q, 2):
            c = a * L / q
            for j in range(int(c) - 8, int(c) + 9): mask[j % L] = False
    sup_minor = S[mask].max()
    j_arg = np.argmax(np.where(mask, S, -1)); fr = Fraction(int(j_arg), L).limit_denominator(64)
    print(f"V=2^{k}: max over minor grid |S|/V = {sup_minor:.4f} at theta~{j_arg}/{L} (~{fr}); log_V(sup*V) = {math.log(sup_minor*V)/math.log(V):.3f}; "
          f"S(1/3)={S[round(L/3)]:.4f} S(1/5)={S[round(L/5)]:.4f} S(1/7)={S[round(L/7)]:.4f} S(1/2^5)={S[L>>5]:.4f}")
