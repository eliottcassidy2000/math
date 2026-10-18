#!/usr/bin/env python3
"""S43 referee: rerun the half-period battery (u((p-1)/2) == (d/p) a_g(p) mod p, mod p^2) for the correctly defined
sequences that the session misdefined or omitted (Cooper s18 per Cooper 2012 / OSS 2016, Almkvist-Zudilin (delta),
(epsilon) = sum C(n,k)^2 C(2k,n)^2, Zagier E), against all 21 distinct M24 eta-products and twists d.
Usage: python audit_battery_corrected.py [PMAX]
"""
import sys
from math import comb, factorial, lcm

SHAPES = {'1A': {1: 24}, '2A': {1: 8, 2: 8}, '2B': {2: 12}, '3A': {1: 6, 3: 6}, '3B': {3: 8}, '4A': {2: 4, 4: 4},
          '4B': {1: 4, 2: 2, 4: 4}, '4C': {4: 6}, '5A': {1: 4, 5: 4}, '6A': {1: 2, 2: 2, 3: 2, 6: 2}, '6B': {6: 4},
          '7AB': {1: 3, 7: 3}, '8A': {1: 2, 2: 1, 4: 1, 8: 2}, '10A': {2: 2, 10: 2}, '11A': {1: 2, 11: 2},
          '12A': {2: 1, 4: 1, 6: 1, 12: 1}, '12B': {12: 2}, '14AB': {1: 1, 2: 1, 7: 1, 14: 1},
          '15AB': {1: 1, 3: 1, 5: 1, 15: 1}, '21AB': {3: 1, 21: 1}, '23AB': {1: 1, 23: 1}}

def eta_product(shape, N):
    c = [0] * (N + 1)
    c[0] = 1
    for a, m in shape.items():
        for _ in range(m):
            for n in range(1, N // a + 1):
                st = a * n
                for k in range(N, st - 1, -1):
                    c[k] -= c[k - st]
    return c

def kron(d, p):
    d %= p
    if d == 0:
        return 0
    return 1 if pow(d, (p - 1) // 2, p) == 1 else -1

SEQ = {
    'Cooper s18 (correct)': lambda n: 1 if n == 0 else sum((-1) ** k * comb(n, k) * comb(2 * k, k) * comb(2 * n - 2 * k, n - k) * (comb(2 * n - 3 * k - 1, n) + comb(2 * n - 3 * k, n)) for k in range(n // 3 + 1)),
    'Almkvist-Zudilin (delta)': lambda n: sum((-1) ** k * 3 ** (n - 3 * k) * comb(n, 3 * k) * comb(n + k, n) * factorial(3 * k) // factorial(k) ** 3 for k in range(n // 3 + 1)),
    '(epsilon) sum C(n,k)^2 C(2k,n)^2': lambda n: sum(comb(n, k) ** 2 * comb(2 * k, n) ** 2 for k in range(n + 1)),
    'Zagier E (= session "Zagier D")': lambda n: sum(comb(n, k) * comb(2 * k, k) * comb(2 * n - 2 * k, n - k) for k in range(n + 1)),
    'Zagier C': lambda n: sum(comb(n, k) ** 2 * comb(2 * k, k) for k in range(n + 1)),
    'Zagier F': lambda n: sum((-1) ** k * 8 ** (n - k) * comb(n, k) * sum(comb(k, j) ** 3 for j in range(k + 1)) for k in range(n + 1)),
}
D = [1, -1, 2, -2, 3, -3, 5, -5, 6, -6, 7, -7]

def main():
    PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 400
    N = PMAX + 2
    ps = [p for p in range(5, PMAX + 1) if all(p % q for q in range(2, int(p ** 0.5) + 1))]
    forms = {}
    for g, sh in SHAPES.items():
        c = eta_product(sh, N)
        off = sum(a * m for a, m in sh.items()) // 24
        forms[g] = ([0] * off + c)[:N + 1], lcm(*sh) * min(sh)
    for name, f in SEQ.items():
        vals = {p: f((p - 1) // 2) for p in ps}
        hits = []
        for g, (a, level) in forms.items():
            for d in D:
                tested = [p for p in ps if level % p and 6 % p and d % p]
                if len(tested) < 10:
                    continue
                m1 = all((vals[p] - kron(d, p) * a[p]) % p == 0 for p in tested)
                if m1:
                    m2 = all((vals[p] - kron(d, p) * a[p]) % (p * p) == 0 for p in tested)
                    hits.append('%s*(%d/p) mod %s' % (g, d, 'p^2' if m2 else 'p'))
        print('%-36s half-period matches (primes 5..%d): %s' % (name, PMAX, ', '.join(hits) or 'none'), flush=True)

if __name__ == '__main__':
    main()
