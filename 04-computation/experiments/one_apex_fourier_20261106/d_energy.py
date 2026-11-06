# S83 (D): exact energy E_n and smoothed energy E_n^(Delta) (Delta = 2^floor(n/3)) of the one-apex counts
# on X_n = [2^(n-1), 2^n), against the proved bound (1.39 n + 19) 2^(4n/3), and the counting corollary.
import numpy as np, math
def g_vec(v):
    v = v.astype(np.int64); acc = np.zeros_like(v)
    cur = v.copy()
    while True:
        alive = cur > 0
        if not alive.any(): break
        bit = cur & 1
        w = cur >> 1
        acc += np.where(alive & (bit == 0), w + (w & 1), 0)
        cur = w
    return 1 + acc
def g(v):
    acc = 0
    while v > 0:
        if v & 1: v >>= 1
        else:
            w = v >> 1; acc += w + (w & 1); v = w
    return 1 + acc
t = np.arange(1, 5000); assert (g_vec(t) == np.array([g(int(x)) for x in t])).all()
print(" n  m   E_n/|X|  E^(D)/|X|  E^(D)/2^(4n/3)  bound/2^(4n/3)  distinct  |X|^2/bound  sqrt(2^n)/4")
for n in range(12, 27):
    m = n // 3; C = 1 << m; Delta = C
    x = np.sort(g_vec(np.arange(1 << (n - 1), 1 << n, dtype=np.int64)))
    _, cnt = np.unique(x, return_counts=True)
    E = int((cnt.astype(np.int64) ** 2).sum())
    hi = np.searchsorted(x, x + Delta, side='right'); lo = np.searchsorted(x, x - Delta, side='left')
    Es = int((hi - lo).sum())
    X = 1 << (n - 1); s = 2 ** (4 * n / 3)
    bound = (1.39 * n + 19) * s
    assert x.max() <= (1 << n) + 1
    print(f"{n:2d} {m:2d}  {E/X:8.4f}  {Es/X:9.3f}  {Es/s:13.4f}  {bound/s:13.2f}  {len(cnt):9d}  {X*X/bound:10.1f}  {2**(n/2)/4:10.1f}")
