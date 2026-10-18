#!/usr/bin/env python3
"""Is s7 == (-1)^n A at n = (p-1)/2 (mod p^2) an instance of a general swap C(n+k,k) <-> C(2k,n)? (S43, HYP-9291)

For weights F(n,k) in a list, compare U(n) = sum_k F(n,k) C(n+k,k) with V(n) = sum_k F(n,k) C(2k,n) at n = (p-1)/2:
test U == eps V (mod p^2) and (mod p) for eps = (-1)^n, +1, and record which F admit the swap.
Session opus-2026-10-17-S43.
"""
from math import comb
F = {
    'C(n,k)^2 C(n+k,k)  [A vs s7]': lambda n, k: comb(n, k) ** 2 * comb(n + k, k),
    'C(n,k)^2           [Apery b vs sum C(n,k)^2 C(2k,n)]': lambda n, k: comb(n, k) ** 2,
    'C(n,k) C(n+k,k)': lambda n, k: comb(n, k) * comb(n + k, k),
    'C(n,k)': lambda n, k: comb(n, k),
    'C(n,k)^3': lambda n, k: comb(n, k) ** 3,
    'C(n,k)^2 C(2k,k)': lambda n, k: comb(n, k) ** 2 * comb(2 * k, k),
    'C(n,k)^2 C(n+k,k)^2': lambda n, k: comb(n, k) ** 2 * comb(n + k, k) ** 2,
    'C(n,k)^3 C(n+k,k)': lambda n, k: comb(n, k) ** 3 * comb(n + k, k),
    'C(n,k)^4': lambda n, k: comb(n, k) ** 4,
    '1': lambda n, k: 1,
}
primes = [p for p in range(5, 400) if all(p % d for d in range(2, int(p ** 0.5) + 1))]
for name, f in F.items():
    res = []
    for mod in (2, 1):
        for epsname in ('(-1)^n', '+1'):
            ok = True
            for p in primes:
                n = (p - 1) // 2
                U = sum(f(n, k) * comb(n + k, k) for k in range(n + 1))
                V = sum(f(n, k) * comb(2 * k, n) for k in range(n + 1))
                e = (-1) ** n if epsname == '(-1)^n' else 1
                if (U - e * V) % p ** mod:
                    ok = False
                    break
            if ok:
                res.append('mod p^%d with eps %s' % (mod, epsname))
    print('%-58s %s' % (name, '; '.join(res) if res else 'no swap'))
