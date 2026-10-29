"""Extension (audit upgrade): strong-tournament and strong-semicomplete Camion tests for F_9 and F_11 at small known
prime factors (one failing factor suffices for N to fail).  Same series-inversion method as audit_fermat.py."""
import sympy
from audit_fermat import strong_mod_p
F = lambda n: 2 ** (2 ** n) + 1
for n, p in ((9, 2424833), (11, 319489), (11, 974849)):
    N = F(n)
    assert N % p == 0 and sympy.isprime(p)
    r = N % (p - 1) or (p - 1)
    q = (N - r) // (p - 1)
    s2 = strong_mod_p(r, p, 2)[r]
    s3 = strong_mod_p(r, p, 3)[r]
    D = (s2 - pow(2, (r - 1) // 2, p)) % p
    e3 = sympy.jacobi_symbol(3, p)
    s3N = pow(e3 % p, q, p) * s3 % p
    print('F_%d at p=%d: r=%d, q parity %d; D(r) mod p = %d (strong test passes iff 0); (3/p)=%d; '
          's_3(N) mod p = %d (Pepin-type pass iff %d)' % (n, p, r, q % 2, D, e3, s3N, p - 1), flush=True)
