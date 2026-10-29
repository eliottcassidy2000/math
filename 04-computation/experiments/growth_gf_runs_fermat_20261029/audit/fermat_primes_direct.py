"""Direct (non-circular) check of THM-4847 (C3) at the Fermat primes 5, 17, 257: s_b(p) mod p computed by Wright's
recurrence with exact binomials, working modulo p^2 (so the vanishing of C(p,k) mod p is not used)."""
from math import comb
for p in (3, 5, 17, 257):
    mod = p * p
    for b in (2, 3):
        T = [pow(b, comb(m, 2), mod) for m in range(p + 1)]
        s = [0] * (p + 1)
        for m in range(1, p + 1):
            s[m] = (T[m] - sum(comb(m, k) % mod * s[k] * T[m - k] for k in range(1, m))) % mod
        v = s[p] % p
        print('p=%d b=%d: s_b(p) mod p = %d (%s)' % (p, b, v, '-1' if v == p - 1 else ('+1' if v == 1 else str(v))))
