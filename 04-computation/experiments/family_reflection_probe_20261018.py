#!/usr/bin/env python3
"""S43 probe: which of the four families has a mod-p reflection F(m) == eps (m/c)^d F(c^2/m)? Only A1 (c = 64, d = (p-1)/2,
eps = (-1/p)); none for A2 (c = 108), A3 (256), A4 (1728), p <= 37. So THM-4686 is specific to C(2k,k)^3.
Usage: python family_reflection_probe_20261018.py. Session opus-2026-10-17-S43."""
from math import comb
# Does F^A(m) = sum_{k<p} A(k) m^-k (mod p) satisfy F(m) == eps * (m/c)^d * F(c^2/m)  for some c, d, eps?
# A1 control: c = 64, A = C(2k,k)^3, expect d = (p-1)/2 (palindrome).  A2: c = 108?  A3: c = 256?  A4: c = 1728?
FAM = {
    'A1': (lambda k: comb(2*k, k)**3, 64),
    'A2': (lambda k: comb(2*k, k)**2*comb(3*k, k), 108),
    'A3': (lambda k: comb(2*k, k)**2*comb(4*k, 2*k), 256),
    'A4': (lambda k: comb(2*k, k)*comb(3*k, k)*comb(6*k, 3*k), 1728),
}
for name, (Af, c) in FAM.items():
    print(name, 'c =', c)
    for p in [13, 17, 19, 23, 29, 31, 37]:
        A = [Af(k) % p for k in range(p)]
        deg = max(k for k in range(p) if A[k])
        found = None
        for d in range(0, p - 1):
            for eps in (1, p - 1):
                ok = True
                for m in range(1, p):
                    if m % p == 0 or (c * c) % p == 0:
                        continue
                    mi = pow(m, p - 2, p)
                    m2 = c * c * mi % p
                    m2i = pow(m2, p - 2, p)
                    F1 = sum(A[k] * pow(mi, k, p) for k in range(p)) % p
                    F2 = sum(A[k] * pow(m2i, k, p) for k in range(p)) % p
                    rhs = eps * pow(m * pow(c, p - 2, p) % p, d, p) * F2 % p
                    if (F1 - rhs) % p:
                        ok = False
                        break
                if ok:
                    found = (d, 'eps=+1' if eps == 1 else 'eps=-1')
                    break
            if found:
                break
        print('  p=%d: degree of A(k) mod p = %d; reflection m -> c^2/m: %s' % (p, deg, found))
