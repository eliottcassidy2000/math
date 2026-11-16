#!/usr/bin/env python3
"""residues.py: replay of the residue theorem.  (A) products of the explicit counts 3, 5, 11, 13, 17, 257 hit every odd
residue modulo 2^a 3^b 5^c (closure of the multiplicative monoid mod m).  (B) primes <= 63: Frattini check — the vectors
(x mod 4 == 3, Legendre symbols (x/p)) of the realised numbers x generate (Z/2)^(1 + #primes); and realised primitive
roots mod p^2.  (C) with 7 in the modulus: closure including 49, 343 and the L_3 points 1 + 3*2^s (s = 1 mod 3, s >= 4)
hits every odd residue modulo 7^e * m' for small e and m'."""
from itertools import product
def closure(gens, m, odd_only):
    seen = {1 % m}; frontier = [1 % m]
    while frontier:
        nxt = []
        for x in frontier:
            for g in gens:
                y = (x * g) % m
                if y not in seen: seen.add(y); nxt.append(y)
        frontier = nxt
    target = [r for r in range(m) if (r % 2 == 1 or not odd_only)]
    missing = [r for r in target if r not in seen]
    return missing
GA = [3, 5, 11, 13, 17, 257]
print("(A) odd residues modulo 2^a 3^b 5^c hit by products of", GA)
bad = []
for a in range(0, 9):
    for b in range(0, 5):
        for c in range(0, 4):
            m = 2**a * 3**b * 5**c
            if m > 3_000_000 or m == 1: continue
            miss = closure(GA, m, a >= 1)
            if miss: bad.append((a, b, c, miss[:5]))
print("   (a <= 8, b <= 4, c <= 3, m <= 3e6): moduli with a missed residue:", bad if bad else "none")
for m in [2**16, 3**9, 5**7, 2**10*3**3, 2**6*5**4, 3**4*5**4]:
    miss = closure(GA, m, m % 2 == 0); print(f"   m={m}: missed {len(miss)}")
# (B) Frattini check for primes <= 63 except 7, and with 7
def legendre(x, p):
    return 0 if pow(x, (p-1)//2, p) == 1 else 1
primes = [p for p in range(3, 64) if all(p % q for q in range(2, p))]
def frattini_rank(gens, plist):
    vecs = []
    for x in gens:
        v = [1 if x % 4 == 3 else 0] + [legendre(x, p) for p in plist if x % p]
        if any(x % p == 0 for p in plist): continue
        vecs.append(v)
    # rank over GF(2)
    rows = [int(''.join(map(str, v)), 2) for v in vecs]; rank = 0; basis = []
    for r in rows:
        for b in basis: r = min(r, r ^ b)
        if r: basis.append(r); rank += 1
    return rank, 1 + len(plist)
def primitive_root_mod_p2(x, p):
    m = p*p; phi = p*(p-1)
    if x % p == 0: return False
    for q in set(f for f in range(2, phi+1) if phi % f == 0 and all(f % r for r in range(2, f))):
        if pow(x, phi//q, m) == 1: return False
    return True
realised = [x for x in range(3, 200, 2) if x not in (7, 21)]
for plist, name in [([p for p in primes if p != 7], "primes <= 63 except 7"), (primes, "primes <= 63 including 7")]:
    rank, need = frattini_rank(realised, plist)
    roots = {p: [x for x in realised if primitive_root_mod_p2(x, p)][:3] for p in plist}
    print(f"(B) {name}: Frattini rank of the realised numbers 3..199 (not 7, 21) = {rank} of {need} needed ->", "generates Q'" if rank == need else "FAILS")
    print("    realised primitive roots mod p^2 (first three):", {p: roots[p] for p in plist})
# (C) with 7: multiples of 7 from L_3 and 49, 343
L3 = [1 + 3*2**s for s in range(4, 61) if s % 3 == 1]
assert all(x % 7 == 0 for x in L3)
GC = GA + [49, 343] + L3 + [x for x in realised if x % 7 and x <= 61]
bad = []
for e in (1, 2, 3):
    for mp in (1, 2, 4, 8, 16, 3, 9, 5, 25, 15, 45, 11, 13, 33, 65, 120, 360):
        m = 7**e * mp
        miss = closure(GC, m, m % 2 == 0)
        if miss: bad.append((e, mp, miss[:4]))
print("(C) moduli 7^e * m' (e <= 3) with a missed odd residue using", "GA + {49, 343} + L_3 multiples of 7 + realised units <= 61:", bad if bad else "none")
# the s-avoidance: for each target u' mod 7 and each prime p <= 63 (p != 7), the classes s mod ord_p(2) to avoid
print("(C') L_3 multiples of 7: u' = (1+3*2^s)/7 mod 7 by s mod 21 (s = 1 mod 3):", {s: ((1+3*2**s)//7) % 7 for s in range(4, 25) if s % 3 == 1})
