#!/usr/bin/env python3
"""Dual Burnside for the LRC: orbit side (round=A000016, open times) vs fix side
(self-converse regular config, boundary times). opus-2026-06-02-S565."""
from sympy import totient, divisors
from sympy.utilities.iterables import partitions
from math import gcd, factorial
import collections

def A000016(m):  # round / open = ORBIT side
    return sum(int(totient(d)) * 2**(m // d) for d in divisors(m) if d % 2 == 1) // (2 * m)

def A000568(n):  # all tournaments (odd-cycle Burnside)
    tot = 0
    for p in partitions(n):
        cyc = []
        for L, mu in p.items():
            cyc += [L] * mu
        if any(l % 2 == 0 for l in cyc):
            continue
        cc = collections.Counter(cyc); autom = 1
        for L, mu in cc.items():
            autom *= (L**mu) * factorial(mu)
        e = sum((l - 1) // 2 for l in cyc) + sum(gcd(cyc[i], cyc[j])
                for i in range(len(cyc)) for j in range(i + 1, len(cyc)))
        tot += factorial(n) // autom * 2**e
    return tot // factorial(n)

def regular_round(m):
    h = (m - 1) // 2
    return {i: set((i + k) % m for k in range(1, h + 1)) for i in range(m)}

def is_self_converse(T, m):          # T = T^op via i -> -i
    return all(((-i) % m) in T[((-j) % m)] for i in T for j in T[i])

def rot_aut(m):
    T = regular_round(m)
    return sum(1 for r in range(m)
               if all(set((j + r) % m for j in T[i]) == T[(i + r) % m] for i in T))

if __name__ == "__main__":
    print("m   A000016(round)   A000568(all)        round/total   self-conv  rot-aut")
    for m in [3, 5, 7, 9, 11, 13]:
        a, b = A000016(m), A000568(m)
        print(f"{m:2d}  {a:>13d}   {b:>16d}   {a/b:.2e}   "
              f"{str(is_self_converse(regular_round(m), m)):>5s}     {rot_aut(m)}")
    print("\nORBIT side = round/A000016 (open, IGNORE); FIX side = self-converse "
          "regular config (boundary, WORRY; S547 (q,q), Fix=2^(n-1)).")
