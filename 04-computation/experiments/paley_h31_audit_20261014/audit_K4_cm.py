#!/usr/bin/env python3
"""AUDIT: K_4 sum S(q) = sum_{x,y} chi(x(x-1)y(y-1)(y-x)) = Sigma*(K_4)/(q(q-1)):
   q = 3 mod 4: 0 (THM-4652 part 4);  q = 1 mod 4: claimed 4a^2 - 2q, q = a^2 + b^2, a odd (= 2(a^2 - b^2)).
Own O(q^2) numba loop, all primes < QMAX."""
import sys
import numpy as np
import numba as nb

QMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 1000


@nb.njit(cache=True)
def S(q, chi):
    tot = 0
    for x in range(q):
        fx = chi[(x * (x - 1 + q)) % q]
        if fx == 0:
            continue
        for y in range(q):
            tot += fx * chi[(y * (y - 1 + q)) % q] * chi[(y - x + q) % q]
    return tot


def chi_table(q):
    chi = np.full(q, -1, dtype=np.int64)
    chi[0] = 0
    for x in range(1, q):
        chi[(x * x) % q] = 1
    return chi


def main():
    primes = [p for p in range(3, QMAX) if all(p % d for d in range(2, int(p ** 0.5) + 1))]
    bad3, bad1 = [], []
    n3 = n1 = 0
    for q in primes:
        s = int(S(q, chi_table(q)))
        if q % 4 == 3:
            n3 += 1
            if s != 0:
                bad3.append(q)
        else:
            n1 += 1
            a = next(a for a in range(1, q, 2) if int(round((q - a * a) ** 0.5)) ** 2 == q - a * a)
            if s != 4 * a * a - 2 * q:
                bad1.append((q, s, 4 * a * a - 2 * q))
    print(f'q = 3 mod 4: {n3} primes < {QMAX}, S != 0 at {bad3 if bad3 else "none"}')
    print(f'q = 1 mod 4: {n1} primes < {QMAX}, S != 4a^2 - 2q at {bad1 if bad1 else "none"}')


if __name__ == '__main__':
    main()
