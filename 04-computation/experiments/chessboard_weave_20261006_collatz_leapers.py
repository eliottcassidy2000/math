#!/usr/bin/env python3
"""Collatz on the 3-smooth 'exponent board': cycle shapes as chess leapers,
Syracuse steps as (1,v)-leapers, Gersonides' equation as the Fibonacci leapers.
(chessboard-weave session 2026-10-06)"""
from fractions import Fraction
from collections import Counter
import math

def T(n):  # shortcut map on Z
    return n // 2 if n % 2 == 0 else (3 * n + 1) // 2

# 1. all cycles of the shortcut map T with a point of |x| <= 10^5
seen_cycles = {}
for x0 in range(-10**5, 10**5 + 1):
    x, orbit = x0, {}
    for step in range(2000):
        if x in orbit:
            start = x
            cyc = [start]
            y = T(start)
            while y != start:
                cyc.append(y); y = T(y)
            key = min(cyc)
            if key not in seen_cycles:
                seen_cycles[key] = cyc
            break
        orbit[x] = step
        x = T(x)
print("1. cycles of T(n) = n/2 or (3n+1)/2 met from |x0| <= 1e5:")
for key in sorted(seen_cycles):
    cyc = seen_cycles[key]
    p = sum(1 for y in cyc if y % 2)        # odd steps
    A = len(cyc)                            # every T-step halves once: A = T-period
    # in the standard 3n+1 map: halvings = A, odd steps = p
    print(f"   min {key:4d}: T-period {len(cyc):2d}, odd steps p={p}, halvings A={A}, leaper (p,A)=({p},{A}),"
          f" 2^A-3^p={2**A - 3**p}")

# 2. Gersonides: |2^A - 3^p| = 1
sols = [(A, p) for A in range(0, 400) for p in range(0, 260) if abs(2**A - 3**p) == 1]
print("\n2. |2^A - 3^p| = 1 with A<400, p<260:", sols)
F = [0, 1]
while len(F) < 12:
    F.append(F[-1] + F[-2])
fib_leapers = [(F[n], F[n + 1]) for n in range(0, 6)]
print("   Fibonacci leapers Q^n(0,1) = (F_n, F_{n+1}), n=0..5:", fib_leapers)
print("   Gersonides solutions with A>=1 as (p,A):", [(p, A) for (A, p) in sols if A >= 1],
      "== first four Fibonacci leapers:", sorted((p, A) for (A, p) in sols if A >= 1) == fib_leapers[:4])
L = [2, 1]
while len(L) < 12:
    L.append(L[-1] + L[-2])
print("   the -17 cycle (p,A) = (7,11) = (L_4, L_5):", (L[4], L[5]),
      " mediant of (2,3) and (5,8):", (2 + 5, 3 + 8))

# 3. convergents of log2(3) and of phi
def cf(x, k):
    out = []
    for _ in range(k):
        a = math.floor(x); out.append(a); x = 1 / (x - a)
    return out
def convergents(c):
    h0, h1, k0, k1 = 1, c[0], 0, 1
    out = [(h1, k1)]
    for a in c[1:]:
        h0, h1 = h1, a * h1 + h0
        k0, k1 = k1, a * k1 + k0
        out.append((h1, k1))
    return out
l23 = math.log2(3)
c = cf(l23, 10)
print("\n3. log2 3 =", c, " convergents A/p:", convergents(c)[:8])
print("   phi =", cf((1 + 5 ** .5) / 2, 8), " convergents:", convergents(cf((1 + 5 ** .5) / 2, 8))[:7])

# 4. Syracuse step as a (1,v)-leaper on the exponent board (log2, log3): exact Terras densities
print("\n4. Syracuse step S(n) = (3n+1)/2^v for odd n: exponent move (-v,+1)")
names = {1: "bishop (1,1)", 2: "knight (1,2)", 3: "camel (1,3)", 4: "giraffe (1,4)"}
K = 20
cnt = Counter()
for n in range(1, 2**K, 2):
    m = 3 * n + 1
    v = (m & -m).bit_length() - 1
    cnt[v] += 1
tot = 2**(K - 1)
for v in range(1, 7):
    print(f"   v={v}: {names.get(v, '(1,%d)-leaper' % v):14s} density {Fraction(cnt[v], tot)} (exact 2^-v = {Fraction(1, 2**v)})")
print("   mean leap = (-E v, +1) = (-2, +1) = the knight; log-drift per odd step = log(3/4) =", math.log(3 / 4))
print("   critical slope log2 3 =", l23, "lies strictly between bishop (1) and knight (2)")

# 5. 'slide until obstruction': the halving slide n -> n/2^v2(n) stops at the first odd number;
#    3n+1 always lands on the other parity class (the knight's colour switch)
assert all((3 * n + 1) % 2 == 0 for n in range(1, 10**5, 2))
print("\n5. 3n+1 maps every odd n to an even number (colour switch): checked n < 1e5")

# 6. bishop/knight-only words: odd n whose first k Syracuse leaps are all bishops or knights
for k in (1, 2, 4, 8, 12):
    M = 2**(2 * k + 2)
    good = 0
    for n in range(1, M, 2):
        x, ok = n, True
        for _ in range(k):
            m = 3 * x + 1
            v = (m & -m).bit_length() - 1
            if v > 2:
                ok = False; break
            x = m >> v
        good += ok
    print(f"   k={k:2d}: density of odd n (mod 2^{2*k+2}) with k bishop/knight leaps = {Fraction(good, M//2)}"
          f" (3/4)^k = {Fraction(3,4)**k}")
