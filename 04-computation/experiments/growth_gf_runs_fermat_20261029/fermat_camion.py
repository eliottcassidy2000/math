"""Fermat numbers against the rotation (Camion-type) tests of THM-4746.

(a) tournament-count test: 2^C(N,2) = c(N) = 2^((N-1)/2) (mod N) iff ord_N(2) | (N-1)^2/2 -- every F_n, n >= 2, passes.
(b) strong tournament test: s(N) = c(N) (mod N); Korselt form: p | D(r_p), r_p = N mod (p-1) in 1..p-1.
(c) strong semicomplete test (3 states per pair, total 3^C(m,2)): at a prime p >= 5, s3(p) = (3/p) (mod p)
    (THM-4746 (8) with m = 1); for a Fermat prime p >= 5 this is -1. For tournaments s(p) = (2/p), which is -1 at
    p = 5 and +1 at 17, 257 (audit correction).
For a composite N the residue modulo each prime factor p is eps^q * s(r) with N = q(p-1) + r (Theorem 1 and its
semicomplete analogue); it is evaluated here at the small known factors of F_5 and F_6.
"""
from math import comb

def strong_counts_mod(rmax, p, base):
    """s(0..rmax) mod p for strong labelled 'base-state' families with total base^C(m,2) (Wright's recurrence)."""
    T = [pow(base, comb(m, 2), p) for m in range(rmax + 1)]
    # binomials mod p via factorials (rmax < p)
    fact = [1] * (rmax + 1)
    for i in range(1, rmax + 1):
        fact[i] = fact[i - 1] * i % p
    inv = [1] * (rmax + 1)
    inv[rmax] = pow(fact[rmax], p - 2, p)
    for i in range(rmax, 0, -1):
        inv[i - 1] = inv[i] * i % p
    s = [0] * (rmax + 1)
    for m in range(1, rmax + 1):
        acc = T[m]
        for k in range(1, m):
            acc -= fact[m] * inv[k] % p * inv[m - k] % p * s[k] % p * T[m - k]
        s[m] = acc % p
    return s

def legendre(a, p):
    v = pow(a, (p - 1) // 2, p)
    return -1 if v == p - 1 else v

F = lambda n: 2 ** (2 ** n) + 1
out = []
for n in range(2, 12):
    N = F(n)
    ok = ((N - 1) ** 2 // 2) % (2 ** (n + 1)) == 0 and pow(2, 2 ** (n + 1), N) == 1
    out.append('F_%d count test: ord(2)=2^%d divides (N-1)^2/2: %s' % (n, n + 1, ok))
for n, p in ((5, 641), (5, 6700417), (6, 274177)):
    N = F(n)
    assert N % p == 0
    r = N % (p - 1)
    if r == 0:
        r = p - 1
    q = (N - r) // (p - 1)
    line = 'F_%d at p=%d: r=%d, q parity=%d' % (n, p, r, q % 2)
    if r <= 20000:
        s2 = strong_counts_mod(r, p, 2)
        c = pow(2, (r - 1) // 2, p) if r % 2 else 0
        D = (s2[r] - c) % p
        eps2 = legendre(2, p)
        s3 = strong_counts_mod(r, p, 3)
        eps3 = legendre(3, p)
        sN3 = pow(eps3, q, p) * s3[r] % p
        line += ('; strong tournament: D(r) mod p = %d (pass iff 0); (2/p)=%d'
                 '; strong semicomplete: s3(N) mod p = %d (Pepin-type pass iff p-1), (3/p)=%d'
                 % (D, eps2, sN3, eps3))
    else:
        line += '; r too large for the O(r^2) recurrence here'
    out.append(line)
# primes: strong counts s_b(p) mod p computed DIRECTLY by Wright's recurrence modulo p^2 with exact binomials
# (no binomial-vanishing assumption; the earlier draft only evaluated 3^C(p,2), which was circular -- audit).
def strong_direct(p, b):
    mod = p * p
    s = [0] * (p + 1)
    T = [pow(b, comb(m, 2), mod) for m in range(p + 1)]
    for m in range(1, p + 1):
        acc = T[m]
        for k in range(1, m):
            acc -= comb(m, k) % mod * s[k] % mod * T[m - k]
        s[m] = acc % mod
    v = s[p] % p
    return -1 if v == p - 1 else v
for p in (3, 5, 17, 257):
    out.append('p=%d: s_2(p) mod p = %d, (2/p) = %d; s_3(p) mod p = %d, (3/p) = %d (direct Wright mod p^2)'
               % (p, strong_direct(p, 2), legendre(2, p), strong_direct(p, 3), legendre(3, p)))
print('\n'.join(out))
