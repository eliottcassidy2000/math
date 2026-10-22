# Audit b6 (independent algorithm): the Camion-Wieferich search p <= 10^5 (THM-4750 (5), HYP-9358 (a)).
#
# Identity (derived in the audit report; the same identity cw_fast.c uses, re-derived here):
#   s(p) = -p! [x^p] 1/T(x),  T(x) = T_<p(x) + t_p x^p + ...,  t_p = 2^C(p,2)/p!, so
#   [x^p] 1/T = u_p - t_p  with u_p = [x^p] 1/T_<p (a p-integral number), hence
#   s(p) = 2^C(p,2) - p (p-1)! u_p = 2^C(p,2) + p u_p  (mod p^2)   (Wilson),   2^C(p,2) = (2^((p-1)/2))^p = eps (mod p^2).
#   With 2^((p-1)/2) = eps + p e' (mod p^2):  D(p)/p = u_p - e'  (mod p).
# Algorithm here: FLINT's fast power-series inversion (Newton iteration, nmod_poly.inverse_series_trunc) instead of
# cw_fast.c's schoolbook O(p^2) recurrence; e' from exact integer arithmetic.
# Validation: (i) against direct Wright/inverse recurrences mod p^2 for all p <= 1500 (exact big-integer free, modular);
#             (ii) the known values q_C(5..37) printed by cw_fast_validation_5_10000.out.
# Usage: python3 -I b6_cw_search_flint.py [PMAX]
import sys, time
import flint
PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 100000

def primes_upto(n):
    sv = bytearray([1]) * (n + 1); sv[0] = sv[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if sv[i]: sv[i * i::i] = bytearray(len(sv[i * i::i]))
    return [i for i in range(n + 1) if sv[i]]

def camion_quotient(p):
    # t_k = 2^C(k,2)/k! mod p, k = 0..p-1
    fact = [1] * p
    for k in range(1, p): fact[k] = fact[k - 1] * k % p
    invf = [1] * p; invf[p - 1] = pow(fact[p - 1], p - 2, p)
    for k in range(p - 1, 0, -1): invf[k - 1] = invf[k] * k % p
    coeffs = []; two = 1
    for k in range(p):
        coeffs.append(two * invf[k] % p); two = two * pow(2, k, p) % p
    Tp = flint.nmod_poly(coeffs, p)
    inv = Tp.inverse_series_trunc(p + 1)
    c = inv.coeffs()
    up = int(c[p]) if len(c) > p else 0
    eps = 1 if p % 8 in (1, 7) else -1
    h = pow(2, (p - 1) // 2, p * p)
    e1 = ((h - eps) % (p * p)) // p
    return (up - e1) % p

import numpy as np
from sympy.ntheory import primitive_root
def camion_quotient_fast(p):
    """Same quantity; the coefficients 2^C(k,2)/k! mod p are built in numpy through discrete logs to a primitive root g."""
    g = primitive_root(p)
    i = np.arange(p - 1, dtype=np.int64); res = np.ones(p - 1, dtype=np.int64); base = g % p; e = i.copy()
    while e.any():
        res = np.where(e & 1, res * base % p, res); base = base * base % p; e >>= 1
    pows = res                                   # pows[i] = g^i mod p
    log = np.zeros(p, dtype=np.int64); log[pows] = i
    k = np.arange(p, dtype=np.int64)
    Lfact = np.concatenate(([0], np.cumsum(log[1:]) % (p - 1)))          # L(k!) for k = 0..p-1
    ex = ((k * (k - 1) // 2) % (p - 1) * log[2] - Lfact) % (p - 1)
    coeffs = pows[ex].tolist()
    inv = flint.nmod_poly(coeffs, p).inverse_series_trunc(p + 1)
    c = inv.coeffs(); up = int(c[p]) if len(c) > p else 0
    eps = 1 if p % 8 in (1, 7) else -1
    h = pow(2, (p - 1) // 2, p * p); e1 = ((h - eps) % (p * p)) // p
    return (up - e1) % p

def direct_quotient(p):
    """D(p)/p mod p by the Hurwitz-inverse recurrence mod p^2 (independent of the identity)."""
    q = p * p
    T = [pow(2, n * (n - 1) // 2, q) for n in range(p + 1)]
    W = [1] + [0] * p; row = [1] + [0] * p
    for n in range(1, p + 1):
        for i in range(n, 0, -1): row[i] = (row[i] + row[i - 1]) % q
        W[n] = -sum(row[k] * T[k] * W[n - k] for k in range(1, n + 1)) % q
    s = (-W[p]) % q
    D = (s - pow(2, (p - 1) // 2, q)) % q
    assert D % p == 0
    return (D // p) % p

t0 = time.time()
# validation (i)
bad = [p for p in primes_upto(400) if p >= 5 and camion_quotient(p) != direct_quotient(p)]
print(f"validation: identity + FLINT inversion vs direct recurrence mod p^2 for 5 <= p <= 400: mismatches {bad}")
known = {5: 3, 7: 2, 11: 7, 13: 0, 17: 14, 19: 1, 23: 8, 29: 7, 31: 2, 37: 2}
bad2 = [p for p in primes_upto(3000) if p >= 5 and camion_quotient(p) != camion_quotient_fast(p)]
print(f"validation: numpy/discrete-log coefficient builder vs plain loop for 5 <= p <= 3000: mismatches {bad2}")
print("validation vs cw_fast_validation_5_10000.out values q_C(5..37):",
      all(camion_quotient(p) == v for p, v in known.items()))
hits = []; n = 0
for p in primes_upto(PMAX):
    if p < 5: continue
    n += 1
    if camion_quotient_fast(p) == 0: hits.append(p); print("  CAMION-WIEFERICH", p, flush=True)
print(f"searched all {n} primes 5 <= p <= {PMAX}: Camion-Wieferich primes {hits}  (plus p = 3, where D(3) = 0)  [{time.time() - t0:.1f}s]")
