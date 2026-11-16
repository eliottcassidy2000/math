#!/usr/bin/env python3
"""clocks_check.py: exact arithmetic behind the covering-clock and decimal threads (klein-2026-11-15-age).
(1) factorisations of 2^24 - 1 and 2^36 - 1, orders of 2, and which primes are primitive divisors of 2^d - 1 for d | 24
    (Zsigmondy: d = 6 contributes none, 2^6 - 1 = 63 = 3^2 * 7); (2) the Riesel number 509203 and the Sierpinski number
    78557: covering verified over a full period; (3) the ray A_e = 13 * 2^e + 1 (27, 53, 105, 209, ...): divisibility by
    3, 7, 9, 63 and the wheel phases; (4) 7 * 3^j and 13 * 2^e + 1 are disjoint; (5) 2^u = 10^v + 3 * 2^w small search
    and C(n, 3) + 1 = square for n <= 10^6; (6) the five pieces of the (24, 45, 51) sum graph on 1..2024."""
from sympy import factorint, n_order, isprime
from math import gcd, lcm, comb, isqrt
def fac(n): return dict(sorted(factorint(n).items()))
M24 = 2**24 - 1; M36 = 2**36 - 1
print("2^24-1 =", M24, fac(M24)); print("2^36-1 =", M36, fac(M36))
print("orders of 2 mod primes of 2^24-1:", {p: n_order(2, p) for p in fac(M24)})
for d in sorted(k for k in range(1, 25) if 24 % k == 0):
    new = [p for p in fac(2**d - 1) if n_order(2, p) == d]
    print(f"  d = {d:2d}: 2^d-1 = {2**d-1:9d} primitive prime divisors {new}")
# covering checks
def covers(k, sign, primes, period):
    bad = [n for n in range(period) if all((k * 2**n + sign) % p for p in primes)]
    return bad
R = 509203; S = 78557
print("509203*2^n-1 covered by {3,5,7,13,17,241} over n mod 24: uncovered", covers(R, -1, [3,5,7,13,17,241], 24))
print("78557*2^n+1 covered by {3,5,7,13,19,37,73} over n mod 36: uncovered", covers(S, 1, [3,5,7,13,19,37,73], 36))
print("lcm of orders, Riesel set:", lcm(*[n_order(2, p) for p in [3,5,7,13,17,241]]), "; Sierpinski set:", lcm(*[n_order(2, p) for p in [3,5,7,13,19,37,73]]))
print("ord_7(2) =", n_order(2, 7), "(odd: -1 not a power of 2 mod 7, the Brier obstruction of CLOCKS.md)")
# the ray
A = [13 * 2**e + 1 for e in range(1, 13)]
print("13*2^e+1, e=1..12:", A)
print("  divisible by 3 at e:", [e for e in range(1, 61) if (13*2**e+1) % 3 == 0][:8], "(odd e)")
print("  divisible by 7 at e:", [e for e in range(1, 61) if (13*2**e+1) % 7 == 0][:8], "(3 | e)")
print("  divisible by 9 at e:", [e for e in range(1, 61) if (13*2**e+1) % 9 == 0][:8], "(e = 1 mod 6)")
print("  divisible by 63 at e <= 1000:", [e for e in range(1, 1001) if (13*2**e+1) % 63 == 0])
print("  uncovered phases mod 24 by the Riesel set:", covers(13, 1, [3,5,7,13,17,241], 24))
print("  primality of A_e, e=1..12:", [isprime(a) for a in A])
print("27 = (4*21-3)/3 ->", (4*21-3)//3, "; 63 = 3*21; D-images: (3*27+3)/4 =", (3*27+3)//4, "; (13+1)/2 =", (13+1)//2, "; (41+1)/2 =", (41+1)//2)
# Pillai-type search and C(n,3)+1 squares
sols = [(u, v, w) for u in range(1, 200) for v in range(1, 60) for w in range(0, 200) if 2**u == 10**v + 3 * 2**w]
print("2^u = 10^v + 3*2^w, u<200, v<60:", sols)
sq = [n for n in range(3, 10**6 + 1) if isqrt(comb(n, 3) + 1)**2 == comb(n, 3) + 1]
print("n <= 10^6 with C(n,3)+1 a perfect square:", sq, "(C(24,3)+1 = 2025 = 45^2; 2024 = 2^10 + 10^3; 1024 = 10^3 + 24)")
# the five pieces: sums 576, 2025, 2601 on 1..2024
s, m = 24**2, 45**2; V = range(1, m); adj = {v: set() for v in V}
for a in V:
    for t in (s, m, m + s):
        b = t - a
        if 1 <= b < m and b != a: adj[a].add(b)
seen = set(); comps = []
for v in V:
    if v in seen: continue
    stack = [v]; comp = []; seen.add(v)
    while stack:
        x = stack.pop(); comp.append(x)
        for y in adj[x]:
            if y not in seen: seen.add(y); stack.append(y)
    comps.append(comp)
sizes = sorted(len(c) for c in comps)
kinds = []
for c in comps:
    deg = [len(adj[x]) for x in c]; kinds.append(('path' if min(deg) == 1 else 'cycle', len(c), min(c) % 9, max(c)))
print("components of the (576, 2025, 2601)-sum graph on 1..2024:", sizes, sorted(kinds))
print("largest power of two below 2024:", 2**((2024).bit_length() - 1), "; 44*46 =", 44*46, "= C(24,3) =", comb(24, 3))
