#!/usr/bin/env python3
"""primes_near.py: the least prime and the least twin-prime pair above 10^k (k = 1..15), their offsets, the Hardy-Littlewood
expectation for a twin pair within 100 of 10^k; the owner's pairs 10007/10009, 30011/30013, 30021/30027 = 3*(10007/10009),
1e9+7/1e9+9; and the holes 7, 21 in the Mandelbrot counts 2^(p-1) - 1."""
from sympy import nextprime, isprime, primerange
import math
C2 = 0.6601618158
print("k  least prime > 10^k (offset)   least twin pair > 10^k (offsets)   HL-expected pairs in (10^k, 10^k+100]")
near = 0; at79 = []
for k in range(1, 16):
    N = 10**k; p = nextprime(N); q = p
    while not isprime(q + 2): q = nextprime(q)
    exp = 2 * C2 * 100 / (math.log(N) ** 2)
    if q - N <= 100: near += 1
    if (q - N, q + 2 - N) == (7, 9): at79.append(k)
    print(f"{k:2d}  {p} (+{p-N})   ({q}, {q+2}) (+{q-N}, +{q+2-N})   {exp:.3f}")
print(f"twin pair within 100 of 10^k for {near} of 15 values of k; pair at offsets (7, 9) for k in {at79}")
print("30021 = 3*10007:", 30021 == 3*10007, "; 30027 = 3*10009:", 30027 == 3*10009, "; 30011, 30013 prime:", isprime(30011), isprime(30013), "; 1e9+7, 1e9+9 prime:", isprime(10**9+7), isprime(10**9+9))
print("primes p <= 100 with 21 | 2^(p-1)-1:", [p for p in primerange(2, 101) if (2**(p-1) - 1) % 21 == 0], "; p = 1 mod 6:", [p for p in primerange(2, 101) if p % 6 == 1])
print("Mandelbrot hyperbolic components of prime period p = 2^(p-1)-1:", {p: 2**(p-1) - 1 for p in (2, 3, 5, 7, 11, 13)}, "; 63 = 3*21 = 9*7; Zsigmondy exception n = 6")
print("Fermat quotients (2^(p-1)-1)/p:", {p: (2**(p-1) - 1)//p for p in (3, 5, 7, 11, 13, 17, 19)})
print("terms k^(1/k) - 1 of the MRB series, k = 1..8:", [round(k**(1/k) - 1, 4) for k in range(1, 9)], "; maximal at k = 3 (x^(1/x) peaks at e)")
