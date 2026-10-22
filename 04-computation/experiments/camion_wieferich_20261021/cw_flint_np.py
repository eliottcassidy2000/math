# Camion-Wieferich search, O(p log p) per prime: numpy-vectorised coefficients t_k = 2^C(k,2)/k! mod p (log-depth scans and
# vectorised square-and-multiply), FLINT power-series inversion (python-flint nmod_poly.inverse_series_trunc), and
#   q(p) = u_p - e'(p) mod p,  u_p = [x^p] 1/T_<p(x),  e' = (2^((p-1)/2) - eps)/p;   p Camion-Wieferich iff q(p) = 0.
# Validated against cw_fast.c (3, 13 below 10^4; 22409).  mac-mini-2026-10-21-fanopsp.  Usage: python3 cw_flint_np.py P0 P1
import sys, time, numpy as np, flint
P0, P1 = int(sys.argv[1]), int(sys.argv[2])
def primes(a, b):
    s = np.ones(b + 1, dtype=bool); s[:2] = False
    for i in range(2, int(b ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return [int(p) for p in np.nonzero(s)[0] if p >= max(a, 3)]
def prefix_prod(a, p):                     # inclusive prefix products mod p (Hillis-Steele scan)
    a = a.copy(); s = 1
    while s < len(a):
        a[s:] = a[s:] * a[:-s] % p; s <<= 1
    return a
def vpow(base, e, p):                      # base^e mod p, elementwise (base scalar or array), e int64 array >= 0
    r = np.ones_like(e); b = np.full_like(e, base) if np.isscalar(base) else base.copy(); e = e.copy()
    while e.any():
        m = (e & 1).astype(bool); r[m] = r[m] * b[m] % p; b = b * b % p; e >>= 1
    return r
t0 = time.time(); hits = []; last = time.time()
for p in primes(P0, P1):
    k = np.arange(p, dtype=np.int64)
    e = (k * (k - 1) // 2) % (p - 1)
    pw2 = vpow(2, e, p)
    ar = k.copy(); ar[0] = 1
    fact = prefix_prod(ar, p)               # k! mod p
    invfact = vpow(fact, np.full(p, p - 2, dtype=np.int64), p)
    coeffs = (pw2 * invfact % p).tolist()
    U = flint.nmod_poly(coeffs, p).inverse_series_trunc(p + 1)
    up = int(U[p]); eps = 1 if p % 8 in (1, 7) else -1
    h = pow(2, (p - 1) // 2, p * p); e1 = ((h - eps) % (p * p)) // p % p
    if (up - e1) % p == 0: hits.append(p); print("CAMION-WIEFERICH p=%d" % p, flush=True)
    if time.time() - last > 600: print(f"  progress p={p} ({time.time() - t0:.0f}s)", flush=True); last = time.time()
print(f"searched [{P0},{P1}]: Camion-Wieferich primes {hits}  ({time.time() - t0:.1f}s)")
