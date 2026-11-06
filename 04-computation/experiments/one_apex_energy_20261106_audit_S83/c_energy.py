# Independent audit of THM-4887, part 3: exact energies and the (D4) window integral.  From scratch.
#  E_n(0), E_n(Delta) at Delta = 2^floor(n/3), against (D2) (1.39 n + 19) 2^(4n/3), the exact (D1) value at the
#  (D2) split, and the sharpness lower bound 2^(2n-2)(Delta+1)/(2^n+Delta+1); distinct values against
#  2^(2n/3)/(5.56 n + 76);  (D4): I_n = integral_{|theta| <= 2^(-floor(n/3)-2)} |S_n|^2 computed exactly through the
#  autocorrelation of the value histogram, against 3.92 (1.39 n + 19) 2^n, and the central-peak lower bound
#  integral_{|theta| <= 2^(-n-3)} |S_n|^2 >= 2^(n-5).
import math, sys
import numpy as np

def levels_upto(nmax):
    cur = np.array([1], dtype=np.int64)   # X_1 = {1}, g(1) = 1
    yield 1, cur
    for j in range(1, nmax):
        u = np.arange(1 << (j - 1), 1 << j, dtype=np.int64)
        nxt = np.empty(2 * len(cur), dtype=np.int64)
        nxt[0::2] = cur + u + (u & 1)
        nxt[1::2] = cur
        cur = nxt
        yield j + 1, cur

def D1(n, Delta0):
    t, r = divmod(n, 3); m = t; k = 2 * t + r
    C = 2.0 ** m; P = 2.0 ** (k - 1); lB = math.log(P / C)
    if Delta0:
        return P * C + 8 * P * P + P * C * C * (33.74 + 12 * lB) + 3 * P * C * (1 + lB)
    return P * C + 10 * P * P + P * C * C * (35.74 + 12 * lB) + 3 * P * C * (1 + lB)

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 26
NFFT = int(sys.argv[2]) if len(sys.argv) > 2 else 22
print(" n   E_n/|X|   E_n(0)/2^(4n/3) [D1(0)]   E_n(D)/2^(4n/3)  [D1(D)]  (1.39n+19)  LB/2^(4n/3)  distinct  2^(2n/3)/(5.56n+76)  max g <= 2^n")
for n, G in levels_upto(NMAX):
    if n < 12: continue
    t = n // 3; Delta = 1 << t
    X = 1 << (n - 1); s = 2.0 ** (4 * n / 3)
    x = np.sort(G)
    _, cnt = np.unique(x, return_counts=True)
    E0 = int((cnt.astype(np.int64) ** 2).sum())
    ED = int((np.searchsorted(x, x + Delta, side="right") - np.searchsorted(x, x - Delta, side="left")).sum())
    LB = 2.0 ** (2 * n - 2) * (Delta + 1) / (2 ** n + Delta + 1)
    assert int(x.max()) <= (1 << n) and int(x.min()) >= 1
    assert E0 <= (1.39 * n + 19) * s and ED <= (1.39 * n + 19) * s and ED >= LB
    print(f"{n:2d}  {E0/X:7.4f}   {E0/s:9.4f} [{D1(n,True)/s:6.2f}]       {ED/s:8.4f}  [{D1(n,False)/s:6.2f}]  {1.39*n+19:8.2f}"
          f"   {LB/s:8.4f}   {len(cnt):9d}   {2**(2*n/3)/(5.56*n+76):10.2f}   {int(x.max()) <= (1 << n)}")
    sys.stdout.flush()

print()
print("(D4) I_n = integral over |theta| <= 2^(-floor(n/3)-2) of |S_n|^2 (exact via autocorrelation), vs 3.92(1.39n+19)2^n;")
print("     central peak J_n = integral over |theta| <= 2^(-n-3) (>= 2^(n-5) by Re e(.) >= cos(pi/4))")
print(" n   I_n/2^n   bound/2^n   J_n/2^n   I_n/J_n   Fejer: I'_n <= (pi^2/4)E_n(D)/(D+1) [ratio]")
for n, G in levels_upto(NFFT):
    if n < 12: continue
    t = n // 3; Delta = 1 << t
    M = 1 << (n + 1)
    hist = np.bincount(G, minlength=(1 << n) + 1).astype(np.float64)
    L = 1 << (n + 2)
    F = np.fft.rfft(hist, L)
    R = np.fft.irfft(F * np.conj(F), L)          # R[d] = #{(v,v'): g(v) - g(v') = d mod L}
    R = np.rint(R)
    d = np.arange(L); d = np.where(d >= L // 2, d - L, d).astype(np.float64)
    def win(w):
        K = np.where(d == 0, 2 * w, np.sin(2 * np.pi * w * d) / (np.pi * np.where(d == 0, 1, d)))
        return float((R * K).sum())
    w = 2.0 ** (-t - 2)
    I = win(w); J = win(2.0 ** (-n - 3))
    w2 = 1.0 / (2 * (Delta + 1)); I2 = win(w2)
    x = np.sort(G)
    ED = int((np.searchsorted(x, x + Delta, side="right") - np.searchsorted(x, x - Delta, side="left")).sum())
    fej = (math.pi ** 2 / 4) * ED / (Delta + 1)
    bnd = 3.92 * (1.39 * n + 19) * 2 ** n
    assert I <= bnd and J >= 2 ** (n - 5) and I2 <= fej * (1 + 1e-9)
    print(f"{n:2d}  {I/2**n:8.4f}  {bnd/2**n:9.2f}   {J/2**n:7.4f}   {I/J:7.3f}   {I2:.1f} <= {fej:.1f} [{I2/fej:.3f}]")
    sys.stdout.flush()
