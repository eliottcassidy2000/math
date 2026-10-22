#!/usr/bin/env python3
"""AUDIT of THM-4751 (5): covering radius of O_K = Z[(1+sqrt(-q))/2] in C equals (q+1)^2/(16q) (squared), attained at
z* = 1/2 + (q-1)/(4q) sqrt(-q) in K; < 1 iff q < 7 + sqrt(48); among primes q = 3 mod 4 this is q in {3, 7, 11}.
Numerical check: brute-force max-min distance on a fine grid of the fundamental parallelogram, and exact value at z*.
Output: a11_covering.out"""
import os, math
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
out = open(os.path.join(HERE, 'a11_covering.out'), 'w')
def say(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n')
def isprime(n): return n > 1 and all(n % d for d in range(2, int(n ** 0.5) + 1))
for q in (3, 7, 11, 19, 23, 43, 67, 163):
    w = complex(0.5, math.sqrt(q) / 2)
    pts = [a + b * w for a in range(-3, 4) for b in range(-3, 4)]
    best = 0.0; N = 300
    for i in range(N + 1):
        for j in range(N + 1):
            z = (i / N) + (j / N) * w
            d = min(abs(z - p) ** 2 for p in pts)
            best = max(best, d)
    zs = complex(0.5, (q - 1) / (4 * q) * math.sqrt(q))
    dz = min(abs(zs - p) ** 2 for p in pts)
    say(f"q={q}: grid max-min |z-a|^2 = {best:.6f}; at z* = {dz:.6f}; (q+1)^2/(16q) = {(q+1)**2/(16*q):.6f}; norm-Euclidean (<1): {(q+1)**2 < 16*q}")
say("primes q = 3 mod 4 below 10^4 with (q+1)^2 < 16q:", [q for q in range(3, 10000) if isprime(q) and q % 4 == 3 and (q + 1) ** 2 < 16 * q],
    "; 7 + sqrt(48) =", 7 + math.sqrt(48))
out.close()
