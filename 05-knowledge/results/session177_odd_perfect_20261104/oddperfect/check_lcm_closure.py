#!/usr/bin/env python3
from itertools import combinations
from math import prod
from check_odd_order_support import isprime, order, v2, lcm, sigma_mod

def closure(T, K, b):
    L = b
    for r in T:
        L = lcm(L, K[r])
    return (frozenset((r for r, k in K.items() if L % k == 0)), L)

def main():
    P = [p for p in range(3, 44, 2) if isprime(p)]
    cases = 0
    for p in P:
        for role in (0, 1):
            R = [r for r in P if r != p and v2(order(p, r)) <= role]
            K = {r: r if order(p, r) == 1 else order(p, r) for r in R}
            for mask in range(1 << len(R)):
                T = frozenset((r for j, r in enumerate(R) if mask >> j & 1))
                C, L = closure(T, K, 1 + role)
                if not (T <= C and closure(C, K, 1 + role) == (C, L)):
                    raise AssertionError('Verification failed')
                if C == T:
                    if L == 1:
                        L = max([3] + list(K.values())) + 1
                        while not isprime(L):
                            L += 1
                    A = frozenset((r for r in P if r != p and sigma_mod(p, L - 1, r) == 0))
                    if not A == T:
                        raise AssertionError((p, role, T, L, A))
                cases += 1
    if not isprime(4561):
        raise AssertionError('Verification failed')
    if not [order(3, r) for r in (11, 13, 4561)] == [5, 3, 15]:
        raise AssertionError('Verification failed')
    if not all((sigma_mod(3, 14, r) == 0 for r in (11, 13, 4561))):
        raise AssertionError('Verification failed')
    if not (sigma_mod(3, 4, 11) == 0 and sigma_mod(3, 4, 4561) != 0):
        raise AssertionError('Verification failed')
    if not (sigma_mod(3, 2, 13) == 0 and sigma_mod(3, 2, 4561) != 0):
        raise AssertionError('Verification failed')
    print(f'PASS: {cases} closure subsets; exact fixed-point realization; 11+13 -> 4561 implication')
if __name__ == '__main__':
    main()
