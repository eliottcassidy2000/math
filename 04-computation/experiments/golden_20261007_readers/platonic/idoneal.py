"""Idoneal slices.
(1) All negative discriminants D (D = 0,1 mod 4, |D| <= DMAX) whose primitive form class group has
    exponent <= 2 (every reduced primitive form is ambiguous: b = 0, |b| = a, or a = c).  Expect 101 (A003171).
(2) Idoneal n = |D|/4 for D = 0 mod 4 (65 numbers); slices n = 2 mod 4 and squarefree n = 1 mod 4.
(3) Borwein-Choi: brute-force the n <= NMAX not of the form xy+yz+zx (x,y,z >= 1).
(4) Selling/Conway check: n is NOT represented  <=>  every form (a,2b,c), ac-b^2 = n (primitive or not)
    is equivalent to a diagonal one  <=>  every reduced form of det n has b = 0 or a = c ... (computed directly).
"""
from math import gcd, isqrt
from collections import Counter
from sympy import factorint, n_order, primerange, isprime
DMAX = 30000
def reduced_forms(D):
    """all reduced primitive forms (a,b,c), b^2-4ac = D < 0, |b| <= a <= c, b >= 0 if |b| = a or a = c"""
    out = []
    a = 1
    while 3 * a * a <= -D:
        for b in range(-a + 1, a + 1):
            if (b * b - D) % (4 * a): continue
            c = (b * b - D) // (4 * a)
            if c < a: continue
            if b < 0 and (a == c): continue
            if gcd(gcd(a, abs(b)), c) != 1: continue
            out.append((a, b, c))
        a += 1
    return out
exp2 = []
for N in range(3, DMAX + 1):
    D = -N
    if D % 4 not in (0, 1): continue
    F = reduced_forms(D)
    if all(b == 0 or abs(b) == a or a == c for (a, b, c) in F):
        exp2.append(N)
print("discriminants |D| <= %d with class-group exponent <= 2: count = %d, max = %d" % (DMAX, len(exp2), max(exp2)))
odd = [N for N in exp2 if N % 4 == 3]; even = [N for N in exp2 if N % 4 == 0]
idon = sorted(N // 4 for N in even)
print("  odd |D|: %d ; |D| = 0 mod 4: %d" % (len(odd), len(even)))
print("idoneal numbers (n = |D|/4):", len(idon), idon)
print("idoneal counts mod 8:", dict(sorted(Counter(n % 8 for n in idon).items())))
s2 = [n for n in idon if n % 4 == 2]
s1sf = [n for n in idon if n % 4 == 1 and all(e == 1 for e in factorint(n).values())]
s1nsf = [n for n in idon if n % 4 == 1 and n not in s1sf]
print("idoneal n = 2 mod 4 (%d):" % len(s2), s2)
print("squarefree idoneal n = 1 mod 4 (%d):" % len(s1sf), s1sf, " non-squarefree removed:", s1nsf)
# (3) Borwein-Choi brute force
NMAX = 200000
rep = bytearray(NMAX + 1)
for z in range(1, NMAX):
    if 3 * z * z > NMAX: break
    for y in range(z, NMAX):
        base = y * z
        if base + z * y + z * z > NMAX and z * (2 * y + z) > NMAX: pass
        # n = xy + yz + zx with x >= y >= z: n = x(y+z) + yz
        if y * (y + z) + y * z > NMAX: break
        x = y
        n = x * (y + z) + y * z
        while n <= NMAX:
            rep[n] = 1; x += 1; n += (y + z)
nonrep = [n for n in range(1, NMAX + 1) if not rep[n]]
print("Borwein-Choi non-represented n <= %d (%d):" % (NMAX, len(nonrep)), nonrep)
print("  equals {1,4} U idoneal(2 mod 4)?", nonrep == sorted([1, 4] + s2))
# (5) primes in the 19 planes; 2,3 primitive roots?
P = sorted({p for m in s1sf for p in factorint(m)})
print("primes dividing the 19 planes:", P)
print("   p: ord_p(2), ord_p(3), p-1 :", [(p, n_order(2, p), n_order(3, p), p - 1) for p in P if p > 3])
both = [p for p in primerange(5, 200) if n_order(2, p) == p - 1 and n_order(3, p) == p - 1]
print("primes < 200 with both 2 and 3 primitive roots:", both)
print("planes divisible by 19:", [m for m in s1sf if m % 19 == 0], " idoneal divisible by 19:", [n for n in idon if n % 19 == 0])
print("101 discriminants mod 19 counts:", dict(sorted(Counter(N % 19 for N in exp2).items())))
print("101 discriminants mod 18 counts:", dict(sorted(Counter(N % 18 for N in exp2).items())))
print("65 idoneal mod 18:", dict(sorted(Counter(n % 18 for n in idon).items())))
print("65 idoneal mod 19:", dict(sorted(Counter(n % 19 for n in idon).items())))
