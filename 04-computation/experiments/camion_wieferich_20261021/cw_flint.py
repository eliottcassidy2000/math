# Camion-Wieferich test in O(p log p) per prime with FLINT power-series inversion (python-flint nmod_poly):
#   u_p = [x^p] 1/T_<p(x) (mod p),  T_<p = sum_{k<p} 2^C(k,2) x^k/k!,   q(p) = u_p - e'(p),  e' = (2^((p-1)/2) - eps)/p mod p;
# p is Camion-Wieferich (p^2 | D(p)) iff q(p) = 0  (same identity as cw_fast.c: s(p) = eps + p u_p mod p^2).
# mac-mini-2026-10-21-fanopsp.  Usage: python3 cw_flint.py P0 P1
import sys, time, flint
P0, P1 = int(sys.argv[1]), int(sys.argv[2])
def primes(a, b):
    sieve = bytearray([1]) * (b + 1); sieve[0:2] = b'\x00\x00'
    for i in range(2, int(b ** 0.5) + 1):
        if sieve[i]: sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
    return [p for p in range(max(a, 3), b + 1) if sieve[p]]
t0 = time.time(); hits = []
for p in primes(P0, P1):
    # coefficients t_k = 2^C(k,2)/k! mod p: 2^C(k+1,2) = 2^C(k,2) * 2^k, 1/(k+1)! = 1/k! * inv(k+1)
    inv = [0, 1] + [0] * (p - 2)
    for i in range(2, p): inv[i] = (p - (p // i) * inv[p % i] % p) % p
    coeffs = [0] * p; c = 1; twok = 1
    for k in range(p):
        coeffs[k] = c
        if k + 1 < p: c = c * twok % p * inv[k + 1] % p; twok = twok * 2 % p
    T = flint.nmod_poly(coeffs, p)
    U = T.inverse_series_trunc(p + 1)
    up = int(U[p])
    eps = 1 if p % 8 in (1, 7) else -1
    h = pow(2, (p - 1) // 2, p * p); e1 = ((h - eps) % (p * p)) // p % p
    if (up - e1) % p == 0: hits.append(p); print("CAMION-WIEFERICH p=%d" % p, flush=True)
print(f"searched [{P0},{P1}]: Camion-Wieferich primes {hits}  ({time.time() - t0:.1f}s)")
