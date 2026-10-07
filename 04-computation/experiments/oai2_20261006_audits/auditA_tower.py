"""Audit A: independent checks of THM-4561 (tower H_2=[[1,1],[-1,1]], H_2N=[[H,H],[-H^T,H^T]])."""
import numpy as np, itertools, math, random

def tower(k):
    H = np.array([[1, 1], [-1, 1]], dtype=np.int64)
    for _ in range(k - 1):
        H = np.block([[H, H], [-H.T, H.T]])
    return H

def acf(x):  # exact aperiodic autocorrelations C_0..C_{n-1}
    x = np.asarray(x, dtype=np.int64); n = len(x)
    return np.array([int(np.dot(x[:n - u], x[u:])) for u in range(n)], dtype=np.int64)

def merit(x):
    C = acf(x); s = int((C[1:] ** 2).sum())
    return math.inf if s == 0 else len(x) ** 2 / (2 * s)

# (a) closed form, my own derivation: r_i(y) = w_i(y) - 2 sum_{l: i_l=1} [y = i mod 2^l] (-1)^{sum_{m>=l} i_m y_m}
ok = True
for k in range(1, 10):
    H = tower(k); N = 1 << k
    for i in range(N):
        for y in range(N):
            v = (-1) ** bin(i & y).count('1')
            for l in range(k):
                if (i >> l) & 1 and (y - i) % (1 << l) == 0:
                    v -= 2 * (-1) ** bin((i >> l) & (y >> l)).count('1')
            ok &= (v == H[i, y])
print("closed form k<=9:", ok)

# (b) transfer identity by hand-derived formulas, random sequences incl. odd lengths
def n4_S(x):
    C = acf(x); m = len(x)
    n4 = int(C[0] ** 2 + 2 * (C[1:] ** 2).sum())
    S = int(sum(int(C[u]) * int(C[m - u]) for u in range(1, m)))
    return n4, S
rng = random.Random(1)
ok = True
for _ in range(300):
    m = rng.randint(1, 60); x = [rng.choice((1, -1)) for _ in range(m)]
    n4, S = n4_S(x)
    for eps in (1, -1):
        y = x + [eps * v for v in x]
        n4y, Sy = n4_S(y)
        ok &= (n4y == 6 * n4 + 8 * eps * S) and (Sy == n4 + 4 * eps * S)
print("transfer identity (300 random x, both eps):", ok)

# (c) Cauchy-Schwarz consequence F' <= 2F/(F+1)
ok = True
for _ in range(300):
    m = rng.randint(2, 60); x = [rng.choice((1, -1)) for _ in range(m)]
    F = merit(x)
    for eps in (1, -1):
        Fy = merit(x + [eps * v for v in x]); ok &= (Fy <= 2 * F / (F + 1) + 1e-12) and Fy < 2
print("F(x(1+-z^m)) <= 2F/(F+1) < 2:", ok)

# (d) rows F<2, columns max F at k=2,3
for k in range(2, 9):
    H = tower(k)
    rmax = max(merit(H[i]) for i in range(H.shape[0])); cmax = max(merit(H[:, j]) for j in range(H.shape[0]))
    print(f"k={k}: max row F = {rmax:.4f}, max column F = {cmax:.4f}")

# (e) row sums r_i(1): 0 for i != 0
ok = all((tower(k).sum(axis=1)[1:] == 0).all() and tower(k).sum(axis=1)[0] == (1 << k) for k in range(1, 11))
print("row sums: r_0(1)=N, r_i(1)=0 (i>0), k<=10:", ok)

# (f) Rudin-Shapiro switching: exact-ish sup via dense FFT; compare with (3+2sqrt2)sqrt(2N)
for k in range(2, 12):
    N = 1 << k; H = tower(k)
    d = np.array([(-1) ** sum(((y >> m) & 1) * ((y >> (m + 1)) & 1) for m in range(k)) for y in range(N)])
    G = d[:, None] * H * d[None, :]
    M = 16 * N
    sup = np.abs(np.fft.fft(G.astype(float), n=M, axis=1)).max() / math.cos(math.pi * N / (2 * M))
    # also check the decomposition: each spike piece restricted to y = i mod 2^l is Golay (|.|^2 + |partner|^2 = 2 len)
    print(f"k={k}: max_rows sup/sqrtN <= {sup / math.sqrt(N):.3f}  (bound {(3 + 2 * math.sqrt(2)) * math.sqrt(2):.3f})")

# (g) best Walsh row merit factor vs 0.905*(4/(1+sqrt17))^k, and tower best row ratio
def sylvester(k):
    H = np.array([[1]], dtype=np.int64)
    for _ in range(k):
        H = np.block([[H, H], [H, -H]])
    return H
lam = 4 / (1 + math.sqrt(17))
for k in range(6, 14):
    W = sylvester(k); N = 1 << k
    L = 1 << (k + 1)
    Fw = np.fft.rfft(W.astype(float), n=L, axis=1); Cw = np.rint(np.fft.irfft(np.abs(Fw) ** 2, n=L, axis=1)[:, :N])
    Fwalsh = N * N / (2 * (Cw[:, 1:] ** 2).sum(axis=1))
    T = tower(k).astype(float)
    Ft = np.fft.rfft(T, n=L, axis=1); Ct = np.rint(np.fft.irfft(np.abs(Ft) ** 2, n=L, axis=1)[:, :N])
    Ftow = N * N / (2 * (Ct[:, 1:] ** 2).sum(axis=1))
    print(f"k={k}: best Walsh row F = {Fwalsh.max():.5f}, /lam^k = {Fwalsh.max() / lam ** k:.4f}; best tower row F = {Ftow.max():.5f}, ratio tower/Walsh = {Ftow.max() / Fwalsh.max():.3f}")
