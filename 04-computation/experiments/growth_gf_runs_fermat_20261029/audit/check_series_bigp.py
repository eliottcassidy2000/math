"""Cross-check the FLINT series-inversion strong counts against Wright's recurrence mod a 47-bit prime (r = 1500)."""
from math import comb
from audit_fermat import strong_mod_p
p = 67280421310721; r = 1500
for b in (2, 3):
    T = [pow(b, comb(m, 2), p) for m in range(r + 1)]
    fact = [1] * (r + 1)
    for i in range(1, r + 1): fact[i] = fact[i-1] * i % p
    inv = [pow(f, p - 2, p) for f in fact]
    s = [0] * (r + 1)
    for m in range(1, r + 1):
        acc = T[m]
        for k in range(1, m):
            acc -= fact[m] * inv[k] % p * inv[m-k] % p * s[k] % p * T[m-k]
        s[m] = acc % p
    assert s == strong_mod_p(r, p, b), b
print('series inversion == Wright recurrence mod 67280421310721 for r <= 1500, b = 2, 3')
