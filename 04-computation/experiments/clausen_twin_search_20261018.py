#!/usr/bin/env python3
"""Search for new 'Clausen twins' of the Apery numbers at the half-period (S43, around THM-4669).

THM-4669: s_7((p-1)/2) == (-1)^((p-1)/2) A((p-1)/2) mod p^2. We search sums u(n) = sum_k prod of binomial factors
(products of 3 or 4 factors chosen from a library of binomials in (n, k)) for congruences at n = (p-1)/2:
    u(n) == eps * A(n)  or  eps * b(n)  (mod p^2),    eps in {+1, -1, (-1)^n},
where A = Apery zeta(3) numbers (M24 class 4A at the half-period) and b = Apery zeta(2) numbers (class 4C).
Screening on primes 5..97, confirmation on all primes 5..300. Trivial hits (u identically equal to A or b, or to
s_7) are reported separately.
Usage: python clausen_twin_search_20261018.py
Session opus-2026-10-17-S43.
"""
import os, sys, json, itertools
from math import comb


def C(a, b):
    if b < 0 or a < 0 or b > a:
        return 0
    return comb(a, b)


LIB = {
    'C(n,k)': lambda n, k: C(n, k),
    'C(n+k,k)': lambda n, k: C(n + k, k),
    'C(2k,k)': lambda n, k: C(2 * k, k),
    'C(2k,n)': lambda n, k: C(2 * k, n),
    'C(2n-2k,n-k)': lambda n, k: C(2 * n - 2 * k, n - k),
    'C(k,n-k)': lambda n, k: C(k, n - k),
    'C(n+k,2k)': lambda n, k: C(n + k, 2 * k),
    'C(n,2k)': lambda n, k: C(n, 2 * k),
    'C(2k,n-k)': lambda n, k: C(2 * k, n - k),
    'C(n+k,n-k)': lambda n, k: C(n + k, n - k),
}


def A(n): return sum(C(n, k) ** 2 * C(n + k, k) ** 2 for k in range(n + 1))
def b(n): return sum(C(n, k) ** 2 * C(n + k, k) for k in range(n + 1))


def main():
    screen = [p for p in range(5, 98) if all(p % d for d in range(2, int(p ** 0.5) + 1))]
    confirm = [p for p in range(5, 301) if all(p % d for d in range(2, int(p ** 0.5) + 1))]
    targets = {'A': A, 'b': b}
    tv = {t: {p: f((p - 1) // 2) for p in confirm} for t, f in targets.items()}
    names = list(LIB)
    hits = []
    for r in (3, 4):
        for combo in itertools.combinations_with_replacement(names, r):
            fs = [LIB[c] for c in combo]
            def u(n):
                tot = 0
                for k in range(0, n + 1):
                    v = 1
                    for f in fs:
                        v *= f(n, k)
                        if v == 0:
                            break
                    tot += v
                return tot
            uv = {}
            for p in screen:
                uv[p] = u((p - 1) // 2)
            for t in targets:
                for eps in ('+1', '-1', '(-1)^n'):
                    def e(p):
                        n = (p - 1) // 2
                        return 1 if eps == '+1' else (-1 if eps == '-1' else (-1) ** n)
                    if all((uv[p] - e(p) * tv[t][p]) % (p * p) == 0 for p in screen):
                        ok = all((u((p - 1) // 2) - e(p) * tv[t][p]) % (p * p) == 0 for p in confirm)
                        if ok:
                            ident = all(u(n) == targets[t](n) for n in range(12))
                            hits.append(dict(summand=' '.join(combo), target=t, eps=eps, identical=ident))
                            print('HIT  u = sum %-48s == %s * %s (mod p^2), p <= 300%s'
                                  % (' '.join(combo), eps, t, '   [identical sequence]' if ident else ''), flush=True)
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(hits, open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'clausen_twin_search_20261018.json'), 'w'), indent=1)
    print('total hits:', len(hits))


if __name__ == '__main__':
    main()
