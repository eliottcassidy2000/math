#!/usr/bin/env python3
"""Wilson / {2,4,8} / Fermat-multiplier dictionary: exact checks, each row typed STRUCTURAL (a proved identity) or
COINCIDENCE (a finite fact with its frequency baseline).  Usage: python3 wilson_248_fermat_dictionary.py"""
import sys, math, itertools, time
from fractions import Fraction as Fr
from sympy import isprime, primerange, n_order, legendre_symbol, factorint
CHECKS = 0
def check(c, msg):
    global CHECKS
    CHECKS += 1
    if not c: print("FAIL:", msg); sys.exit(1)

def wilson_quotient_mod_p(p):
    m = p * p; f = 1
    for i in range(2, p): f = f * i % m
    check((f + 1) % p == 0, f"Wilson's theorem at {p}")
    return ((f + 1) // p) % p

def fermat_q(a, p): return ((pow(a, p - 1, p * p) - 1) // p) % p

t0 = time.time()
print("R1 (STRUCTURAL, Gersonides 1343) powers of 2 adjacent to a power of 3:")
adj = sorted({(a, b) for a in range(0, 200) for b in range(0, 130) if abs(2 ** a - 3 ** b) == 1})
check(adj == [(1, 0), (1, 1), (2, 1), (3, 2)], f"Gersonides pairs {adj}")
print(f"   |2^a - 3^b| = 1 for a < 200, b < 130 exactly at {adj}: the powers of 2 next to a power of 3 are 2, 4, 8")

print("R2 (STRUCTURAL, THM-4640) the three adjacencies are the three generators of the repair algebra at -1:")
f = lambda a, x: (3 * x + 1) / Fr(2 ** a)
check(f(1, Fr(-1)) == -1, "anchor -1 fixed by f_1 (2 = 3 - 1)")
for b in range(0, 30):
    check(3 * 2 ** (b + 1) + 2 ** b == 3 * 2 ** b + 2 ** (b + 2), "unit carry (4 = 3 + 1)")
    check(9 * 2 ** (b + 3) + 2 ** b == 9 * 2 ** b + 2 ** (b + 6), "Catalan carry (8 = 9 - 1)")
for a in range(1, 40):
    c = Fr(1, 2 ** a - 3)
    check(f(a, c) == c, f"anchor c_{a} fixed by f_{a}")
    check(-3 + 2 ** a == 1 / c, "first-generation numerator K((a)) = 2^a - 3 = 1/c_a")
print("   f_1(-1) = -1; 3*2^(b+1) + 2^b = 3*2^b + 2^(b+2); 9*2^(b+3) + 2^b = 9*2^b + 2^(b+6); "
      "K((a)) = 2^a - 3 = 1/c_a for a < 40; anchors for 2^a in {2,4,8}: -1, 1, 1/5")

print("R2b (STRUCTURAL) the Gersonides solutions index exactly the Collatz cycles with unit denominator |2^A - 3^k| = 1:")
def cycle_points(word):
    # fixed point of f_word: x = (3^k x + b)/2^A  ->  x = b / (2^A - 3^k)
    k = len(word); A = sum(word); b = 0; S = 0
    for j, a in enumerate(word):
        b += 3 ** (k - 1 - j) * 2 ** S; S += a
    return Fr(b, 2 ** A - 3 ** k)
unit = {}
for (A, k) in [(1, 1), (2, 1), (3, 2)]:
    pts = set()
    for word in itertools.product(range(1, A + 1), repeat=k):
        if sum(word) == A:
            x = cycle_points(word); check(x.denominator == 1, "unit denominator gives an integer cycle point"); pts.add(int(x))
    unit[(A, k)] = sorted(pts)
check(unit == {(1, 1): [-1], (2, 1): [1], (3, 2): [-7, -5]}, f"unit-denominator cycles {unit}")
def U(n):
    n = 3 * n + 1
    while n % 2 == 0: n //= 2
    return n
cycles = {}
for n0 in range(-100001, 100002, 2):
    seen = []; n = n0
    for _ in range(400):
        if n in seen:
            cyc = seen[seen.index(n):]; key = min(cyc, key=abs)
            if key not in cycles:
                A = 0; m = key
                for _ in range(len(cyc)):
                    t = 3 * m + 1; a = 0
                    while t % 2 == 0: t //= 2; a += 1
                    A += a; m = t
                cycles[key] = (len(cyc), A, abs(2 ** A - 3 ** len(cyc)))
            break
        seen.append(n); n = U(n)
check(sorted(cycles) == [-17, -5, -1, 1], f"integer cycles found {sorted(cycles)}")
check([cycles[c][2] for c in (-17, -5, -1, 1)] == [139, 1, 1, 1], "cycle denominators")
print(f"   (A,k) = (1,1), (2,1), (3,2) give the integer cycles -1, +1, {{-5,-7}} for every word (and (1,0) gives 0); iterating U on all odd "
      f"|n| <= 100001 finds exactly the cycles {{1}}, {{-1}}, {{-5,-7}}, {{-17,...}} with |2^A - 3^k| = 1, 1, 1, {cycles[-17][2]} "
      f"(the -17 cycle: {cycles[-17][0]} odd steps, {cycles[-17][1]} halvings)")

print("R3 (STRUCTURAL) Wilson's -1 is the product of a finite abelian group = product of its involutions:")
for p in primerange(3, 2000):
    prod = 1
    for a in range(1, p): prod = prod * a % p
    check(prod == p - 1, "Wilson")
    d = n_order(2, p); prod2 = pow(2, d * (d - 1) // 2, p)
    check(prod2 == (p - 1 if d % 2 == 0 else 1), f"product of <2> mod {p}")
print("   prod F_p^* = -1 and prod <2> = 2^(d(d-1)/2) = -1 iff d = ord_p(2) is even, all odd primes < 2000")

print("R4 (STRUCTURAL, Lerch 1905) the Wilson quotient is the sum of the Fermat quotients; q is a logarithm:")
for p in primerange(3, 2000):
    W = wilson_quotient_mod_p(p)
    check(W == sum(fermat_q(a, p) for a in range(1, p)) % p, f"Lerch at {p}")
    q2 = fermat_q(2, p)
    check(fermat_q(4, p) == 2 * q2 % p and fermat_q(8, p) == 3 * q2 % p, "q_p(2^j) = j q_p(2) for j = 1, 2, 3")
print("   W_p = sum_a q_p(a) (mod p) and (q_p(2), q_p(4), q_p(8)) = (1, 2, 3) q_p(2) (mod p), all odd primes < 2000")

print("R5 (STRUCTURAL, THM-307) tournaments: T(p) = f_p - (2/p) w_p - A_p (mod p), with the Fermat multiplier (2/p):")
def tournaments_burnside(n):
    # T(n) = sum over partitions of n into odd parts of 2^e(lam) / z_lam   (exact)
    total = Fr(0)
    def parts(n, mx):
        if n == 0: yield []; return
        for k in range(min(n, mx), 0, -1):
            if k % 2:
                for rest in parts(n - k, k): yield [k] + rest
    for lam in parts(n, n):
        mult = {}
        for k in lam: mult[k] = mult.get(k, 0) + 1
        z = 1
        for k, m in mult.items(): z *= k ** m * math.factorial(m)
        e = sum((k - 1) // 2 for k in lam) + sum(math.gcd(lam[i], lam[j]) for i in range(len(lam)) for j in range(i + 1, len(lam)))
        total += Fr(2 ** e, z)
    check(total.denominator == 1, "T(n) integral")
    return total.numerator
known = [1, 1, 1, 2, 4, 12, 56, 456, 6880, 191536, 9733056, 903753248, 154108311168]   # A000568, n = 0..12
for n in range(1, 13): check(tournaments_burnside(n) == known[n], f"T({n}) vs A000568")
for p in primerange(3, 48):
    T = tournaments_burnside(p)
    leg = legendre_symbol(2, p) if p > 2 else 1
    fp = ((pow(2, (p - 1) // 2) - leg) // p) % p
    wp = wilson_quotient_mod_p(p)
    check(fp == leg * fermat_q(2, p) * pow(2, -1, p) % p, f"f_p = (2/p) q_p(2)/2 at {p}")
    # A_p from the Burnside sum: p! T(p) = 2^C(p,2) + (p-1)! 2^((p-1)/2) + p A_p  ->  A_p mod p
    rest = math.factorial(p) * T - 2 ** (p * (p - 1) // 2) - math.factorial(p - 1) * 2 ** ((p - 1) // 2)
    check(rest % p == 0, "remaining Burnside sum divisible by p")
    Ap = (rest // p) % p
    check(T % p == (fp - leg * wp - Ap) % p, f"THM-307 decomposition at p = {p}")
print("   T(n) reproduced from Burnside for n <= 12 (A000568); THM-307's decomposition and f_p = (2/p) q_p(2)/2 hold for p < 48")

print("R6 (COINCIDENCE) the anchor row 2^a - 3 and the Wilson primes:")
row = []
for a in range(3, 25):
    p = 2 ** a - 3
    if isprime(p):
        row.append((a, p, wilson_quotient_mod_p(p) == 0))
check([x for x in row if x[2]] == [(3, 5, True), (4, 13, True)], f"Wilson primes in the anchor row {row}")
print(f"   primes 2^a - 3 (a <= 24): {[(a, p) for a, p, _ in row]}; Wilson among them: 5 = 2^3 - 3, 13 = 2^4 - 3 only. "
      f"CITED (Costa-Gerbicz-Harvey 2014): no Wilson prime other than 5, 13, 563 below 2*10^13, so none of the form "
      f"2^a - 3 with 25 <= a <= 44 either; 563 is not of this form")

print("R7 (COINCIDENCE) 2 is a primitive root modulo 5, 13 and 563:")
for p in (5, 13, 563): check(n_order(2, p) == p - 1, f"2 primitive mod {p}")
for a in range(3, 60):
    p = 2 ** a - 3
    check(p % 8 == 5, "2^a - 3 = 5 mod 8");
ps = list(primerange(5, 10 ** 6)); frac = sum(1 for p in ps if n_order(2, p) == p - 1) / len(ps)
print(f"   ord_p(2) = p - 1 at 5, 13, 563; every 2^a - 3 (a >= 3) is 5 mod 8, so 2 is a non-residue there (4 | ord_p(2)); "
      f"empirical frequency of '2 primitive' among primes < 10^6 = {frac:.4f} (Artin's constant 0.3740); P(all three) ~ {frac ** 3:.3f}")

print("R8 (COINCIDENCE) p - 2 is absolutely Fermat (a^n = a mod n for all a) at the Wilson primes:")
def korselt(n):
    if isprime(n): return True
    fa = factorint(n)
    return n > 2 and all(e == 1 for e in fa.values()) and len(fa) >= 2 and all((n - 1) % (q - 1) == 0 for q in fa)
for p in (5, 13, 563):
    n = p - 2
    check(korselt(n) and all(pow(a, n, n) == a % n for a in range(n)), f"{n} absolutely Fermat")
cnt = sum(1 for p in primerange(5, 10 ** 6) if korselt(p - 2)); tot = len(list(primerange(5, 10 ** 6)))
print(f"   3, 11 prime and 561 = 3*11*17 Carmichael (Korselt); among primes 5 <= p < 10^6, p - 2 is prime or Carmichael "
      f"for {cnt}/{tot} = {cnt / tot:.4f} (twin primes dominate); 563 - 2 = 561 is the least Carmichael number")

print("R9 (STRUCTURAL negative) increasing {2,3}-representations (numerators at the -1 anchor):")
def reps(K, emax):
    out = []
    if K == -1: out.append(())
    for e in range(1, emax):
        if (K - (1 << e)) % 3 == 0:
            P = (K - (1 << e)) // 3
            if P >= -(3 ** (e - 1)):
                for r in reps(P, e): out.append(r + (e,))
    return out
for K, expect in ((5, [(3,), (1, 3)]), (13, [(4,), (1, 4)]), (563, []), (561, []), (1093, []), (3511, [])):
    rs = sorted(reps(K, 64), key=lambda r: (len(r), r))
    check(rs == expect, f"representations of {K}: {rs}")
print("   5 and 13 are numerators only as anchors ((a) and its root twin (1, a)); 563, 561, 1093, 3511 have no "
      "increasing representation with exponents < 64 (pruned search, P >= -3^(e-1) is necessary)")

print("R10 (COINCIDENCE) the Carmichael successor of 561 (results note carmichael_successor_blocks_20261014b):")
N = 684106401
fa = factorint(N); lam = 1
for q in fa: lam = lam * (q - 1) // math.gcd(lam, q - 1)
check(fa == {3: 1, 11: 1, 17: 1, 401: 1, 3041: 1} and korselt(N), "684106401 Carmichael")
check(lam == 15200 and lam + 1 == 27 * 563, "lambda + 1 = 27 * 563")
print(f"   684106401 = 3*11*17*401*3041 is Carmichael, lambda = {lam}, lambda + 1 = 27*563")

print(f"TOTAL exact checks: {CHECKS}  [{time.time() - t0:.0f}s]")
