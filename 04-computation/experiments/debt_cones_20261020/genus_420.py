#!/usr/bin/env python3
"""Genus dictionary of the idoneal discriminant -420 = -4*105: reduced primitive forms, genus characters (Legendre symbols
at 3, 5, 7 and the character mod 4), and the residues mod 420 (and the pairs (mod 15, mod 35)) represented by each form.
Idoneal <=> every genus has one class <=> representation of a number coprime to 420 is decided by its residue mod 420;
the deciding characters at 105 are exactly the registers 3, 5, 7 of the mod-15/mod-35 guards."""
from math import gcd
D = -420
forms = []
for a in range(1, 30):
    for b in range(-a + 1, a + 1):
        if (b * b - D) % (4 * a): continue
        c = (b * b - D) // (4 * a)
        if c < a or (c == a and b < 0): continue
        if gcd(gcd(a, abs(b)), c) != 1: continue
        forms.append((a, b, c))
print('reduced primitive forms of discriminant -420:', forms, 'h =', len(forms))
def legendre(n, p):
    r = pow(n % p, (p - 1) // 2, p); return -1 if r == p - 1 else r
reps = {}
for f in forms:
    a, b, c = f; vals = set()
    for x in range(-40, 41):
        for y in range(-40, 41):
            v = a * x * x + b * x * y + c * y * y
            if v > 0 and gcd(v, 420) == 1: vals.add(v % 420)
    chars = {(legendre(v, 3), legendre(v, 5), legendre(v, 7), (1 if v % 4 == 1 else -1)) for v in vals}
    reps[f] = (sorted(vals), chars)
    print(f"form {f}: {len(vals)} residue classes mod 420; genus characters (3,5,7,mod4) = {chars}; classes mod 105 = {len({v % 105 for v in vals})}, (mod 15, mod 35) pairs = {len({(v % 15, v % 35) for v in vals})}")
allres = [set(r[0]) for r in reps.values()]
print('residue classes disjoint across forms:', all(not (allres[i] & allres[j]) for i in range(len(allres)) for j in range(i + 1, len(allres))))
print('union covers all units mod 420:', len(set().union(*allres)), 'of', sum(1 for v in range(420) if gcd(v, 420) == 1))
idoneal = [1,2,3,4,5,6,7,8,9,10,12,13,15,16,18,21,22,24,25,28,30,33,37,40,42,45,48,57,58,60,70,72,78,85,88,93,102,105,112,120,130,133,165,168,177,190,210,232,240,253,273,280,312,330,345,357,385,408,462,520,760,840,1320,1365,1848]
print('65 idoneal numbers:', len(idoneal), '; 15 idoneal', 15 in idoneal, '; 35 idoneal', 35 in idoneal, '; 105 idoneal', 105 in idoneal, '; 21', 21 in idoneal, '; 63', 63 in idoneal, '; 153', 153 in idoneal, '; 420', 420 in idoneal)
print('idoneal multiples of 105:', [m for m in idoneal if m % 105 == 0], '; idoneal products of exactly three distinct odd primes:', [m for m in idoneal if m % 2 and len([p for p in range(3, m + 1, 2) if m % p == 0 and all(p % q for q in range(3, p, 2))]) == 3 and all((m // p) % p for p in range(3, m + 1, 2) if m % p == 0)])
