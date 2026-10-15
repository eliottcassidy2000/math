#!/usr/bin/env python3
"""Golden 'Fermat' primes: pi(p) = 2^k (k >= 3)  <=>  p | L_(2^(k-2))  (Lucas numbers; L_(2n) = L_n^2 - 2).
Checked directly for all primes p < 10^6; and the Mersenne primes M_q with q = 3 mod 4 are inert with p + 1 = 2^q, so
pi(M_q) divides 2^(q+1) (checked for q = 3, 7, 19, 31, 107, 127 by computing pi via the Lucas-Lehmer-type recursion)."""
from sympy import primerange
import sys
sys.path.insert(0, '.')
from wss_orbits import pi_p

def lucas_pow2_mod(j, p):          # L_(2^j) mod p, j >= 1: L_2 = 3, then L_(2n) = L_n^2 - 2 for even n
    L = 3 % p
    for _ in range(j - 1):
        L = (L * L - 2) % p
    return L

bad = 0; gf = []
for p in primerange(3, 10 ** 6):
    if p == 5: continue
    pi, _ = pi_p(1, p)
    is_pow2 = pi & (pi - 1) == 0 and pi >= 8
    lucas = any(lucas_pow2_mod(j, p) == 0 for j in range(1, 40))
    if is_pow2: gf.append((p, pi))
    if is_pow2 != lucas: bad += 1
    if is_pow2 and lucas_pow2_mod(pi.bit_length() - 3, p) != 0: bad += 1
print(f"pi(p) = 2^k  <=>  p | some L_(2^j), for all primes p < 10^6: mismatches {bad}; golden-Fermat primes: {gf}")
for q in (3, 7, 19, 31, 107, 127):
    p = 2 ** q - 1
    # pi(M_q) divides 2^(q+1): test Q^(2^(q+1)) = I and Q^(2^q) != I via the L_(2^j) criterion: pi = 2^(j+2) iff p | L_(2^j)
    j = next((j for j in range(1, q + 2) if lucas_pow2_mod(j, p) == 0), None)
    print(f"M_{q} = 2^{q}-1 (q mod 4 = {q % 4}): first j with M | L_(2^j): {j} -> pi = 2^{j + 2 if j else '?'}")
