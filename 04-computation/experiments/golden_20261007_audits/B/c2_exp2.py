#!/usr/bin/env python3
"""THM-4566 addendum: (a) negative discriminants |D| <= 30000 whose form class group has exponent <= 2
(all primitive reduced forms ambiguous); (b) split-prime lemma p^2 >= |D|/4 on all of them, for all
split p (p not dividing D, (D/p) = 1; for p = 2: D = 1 mod 8); (c) idoneal n > 9 never = 2 mod 3."""
from math import gcd, isqrt
X = 30000
bad = set()
for a in range(1, isqrt(X // 3) + 2):
    for B in range(1, a):
        c = a + 1
        while 4 * a * c - B * B <= X:
            if gcd(gcd(a, B), c) == 1:
                bad.add(B * B - 4 * a * c)
            c += 1
discs = [-d for d in range(3, X + 1) if d % 4 in (0, 3)]
good = [D for D in discs if D not in bad]
print("exponent<=2 discriminants |D|<=30000:", len(good), " even:", sum(1 for D in good if D % 4 == 0),
      " odd:", sum(1 for D in good if D % 4 != 0), " largest |D|:", max(-D for D in good))
idoneal = sorted(-D // 4 for D in good if D % 4 == 0)
print("idoneal (from D=-4n):", idoneal)
print("idoneal = 2 mod 4:", [n for n in idoneal if n % 4 == 2], len([n for n in idoneal if n % 4 == 2]))
print("squarefree idoneal = 1 mod 4:", [n for n in idoneal if n % 4 == 1 and all(n % (q * q) for q in range(2, isqrt(n) + 1))])
def primes(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, isqrt(n) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]
P = primes(X)
def kron(D, p):
    if p == 2:
        if D % 2 == 0: return 0
        return 1 if D % 8 in (1, 7) else -1
    r = pow(D % p, (p - 1) // 2, p)
    return 0 if r == 0 else (1 if r == 1 else -1)
viol = []
for D in good:
    for p in P:
        if p * p * 4 >= -D:
            break
        if kron(D, p) == 1:
            viol.append((D, p))
print("split-prime lemma violations (split p with 4p^2 < |D|):", viol)
print("idoneal n > 9 with n = 2 mod 3:", [n for n in idoneal if n > 9 and n % 3 == 2])
print("D (exp<=2) where 2 splits (D = 1 mod 8):", [D for D in good if D % 8 == 1])
