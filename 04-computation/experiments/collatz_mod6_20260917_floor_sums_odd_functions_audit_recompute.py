#!/usr/bin/env python3
"""
Adversarial recompute audit (lane floor_sums_odd_functions, lens recompute).
INDEPENDENT of the explorer's script: nothing imported from it; cube roots via an
integer Newton root; class numbers via a separately coded reduced-form count AND a
hard-coded literature table (h(-p) for p = 3 mod 4, p < 200); Pillai via direct gcd
sums; cycle counts of 3n+k via a finite signed search.  Every check raises.
"""
import sys, time
from fractions import Fraction as Fr
from math import gcd

T0 = time.time()
def chk(c, m):
    if not c:
        raise RuntimeError("AUDIT CHECK FAILED: " + m)

def factor(n):
    f = {}; d = 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1; n //= d
        d += 1
    if n > 1: f[n] = f.get(n, 0) + 1
    return f

def isprime(n):
    return n >= 2 and all(n % d for d in range(2, int(n ** 0.5) + 1))

def squarefree(n):
    return all(a == 1 for a in factor(n).values())

def icbrt(m):
    """exact floor cube root via Newton + fix-up (independent of the monotone walk)"""
    if m <= 0: return 0
    x = int(round(m ** (1 / 3)))
    while x ** 3 > m: x -= 1
    while (x + 1) ** 3 <= m: x += 1
    return x

# -------------------------------------------------- S1: odd law, direct residue count
def law(f, n):
    Z = sum(1 for k in range(1, n) if f(k) % n == 0)
    return Fr(sum(f(k) for k in range(1, n)), n) - Fr(n - 1, 2) + Fr(Z, 2)

def fsum(f, n):
    return sum(f(k) // n for k in range(1, n))

odd = {"x": lambda k: k, "x^3": lambda k: k ** 3, "x^5": lambda k: k ** 5,
       "x^3-100x (negative values)": lambda k: k ** 3 - 100 * k,
       "-7x^9+x": lambda k: -7 * k ** 9 + k, "3x": lambda k: 3 * k}
for name, f in odd.items():
    bad = [n for n in range(2, 401) if Fr(fsum(f, n)) != law(f, n)]
    chk(not bad, "S1 %s %s" % (name, bad[:3]))
print("S1 odd law: holds for", list(odd), "all 2<=n<=400 (incl. negative-valued odd f)")
# n=1 degenerate: both sides empty -> 0
chk(fsum(lambda k: k ** 3, 1) == 0 and law(lambda k: k ** 3, 1) == 0, "n=1")
# even f minimal failure, hostile n=2 holds for every f
for name, f in {"x^2": lambda k: k * k, "x^4": lambda k: k ** 4, "const 2": lambda k: 2,
                "x^2+x": lambda k: k * k + k}.items():
    chk(Fr(fsum(f, 2)) == law(f, 2), "n=2 %s" % name)
    chk(Fr(fsum(f, 3)) != law(f, 3), "n=3 should fail %s" % name)
print("S1 even f: n=2 holds for all f tested, n=3 minimal failure for x^2,x^4,2,x^2+x; x^2@3: %d vs %s"
      % (fsum(lambda k: k * k, 3), law(lambda k: k * k, 3)))
# a non-polynomial odd function mod n (values in [0,n) lifted) also obeys the law
import random
random.seed(1)
for n in range(2, 60):
    vals = {}
    for k in range(1, n):
        if k in vals: continue
        v = random.randrange(0, 10 * n)
        vals[k] = v
        if n - k != k: vals[n - k] = -v
    if n % 2 == 0:
        vals[n // 2] = random.choice([0, n // 2]) + n * random.randrange(-3, 3)
    g = lambda k, vals=vals: vals[k]
    chk(Fr(fsum(g, n)) == law(g, n), "generic odd function n=%d" % n)
print("S1 generic odd functions (non-polynomial, random lifts, fixed point in {0,n/2} mod n): law holds n<60")

# -------------------------------------------------- S2: Z_e, truth set of (1), squarefree
def n_e(n, e):
    return eval("*".join(["1"] + ["%d**%d" % (p, -(-a // e)) for p, a in factor(n).items()]))
for e in range(1, 10):
    for n in range(2, 1501):
        Z = sum(1 for k in range(1, n) if pow(k, e, n) == 0)
        chk(Z == n // n_e(n, e) - 1, "Z e=%d n=%d" % (e, n))
print("S2 Z_e(n)=n/n_e-1 verified e=1..9, n<=1500")
rhs1 = lambda n: Fr((n - 2) * (n - 1) * (n + 1), 4)
t1 = [n for n in range(2, 2001) if Fr(fsum(lambda k: k ** 3, n)) == rhs1(n)]
sq = [n for n in range(2, 2001) if squarefree(n)]
chk(t1 == sq and len(sq) == 1214, "truth(1) vs squarefree; count %d" % len(sq))
print("S2 truth set of (1) on [2,2000] == squarefree (1214 values)")
for n in (4, 9, 12, 8, 36):
    d = Fr(fsum(lambda k: k ** 3, n)) - rhs1(n)
    Z = sum(1 for k in range(1, n) if k ** 3 % n == 0)
    print("   n=%2d floorsum=%d RHS=%s Z3=%d dev=%s" % (n, fsum(lambda k: k ** 3, n), rhs1(n), Z, d))
    chk(d == Fr(Z, 2), "dev")
# (1) at prime p incl. p=2,3: user's wording
for p in (2, 3, 5, 7):
    chk(Fr(fsum(lambda k: k ** 3, p)) == rhs1(p), "prime %d" % p)
print("S2 user's identity (1) at p=2,3,5,7 holds (p=2: 0=0)")
for e in (3, 5, 7, 9):
    ts = [n for n in range(2, 801) if Fr(fsum(lambda k: k ** e, n)) == Fr(sum(k ** e for k in range(1, n)), n) - Fr(n - 1, 2)]
    chk(ts == [n for n in sq if n <= 800], "odd e=%d" % e)
print("S2 odd e=3,5,7,9 truth sets on [2,800] == squarefree")
# Z_7 = n/rad -1 first failure
def rad(n): return eval("*".join(["1"] + [str(p) for p in factor(n)]))
bad7 = [n for n in range(2, 300) if sum(1 for k in range(1, n) if pow(k, 7, n) == 0) != n // rad(n) - 1]
chk(bad7 == [256], "Z7 %s" % bad7)
print("S2 Z_7(n)=n/rad(n)-1 first failure n=256; Z_7(256)=%d" % sum(1 for k in range(1, 256) if pow(k, 7, 256) == 0))
# squarefree count to 10^6 by independent Moebius sum Q(x)=sum mu(d) floor(x/d^2)
N = 10 ** 6
mu = [1] * 1001
isp = [True] * 1001
for i in range(2, 1001):
    if isp[i]:
        for j in range(i, 1001, i):
            if j > i: isp[j] = False
            mu[j] *= -1
        for j in range(i * i, 1001, i * i): mu[j] = 0
Q = sum(mu[d] * (N // (d * d)) for d in range(1, 1001))
chk(Q == 607926, "Q(10^6)=%d" % Q)
print("S2 #squarefree<=10^6 via Moebius = 607926 (matches)")

# -------------------------------------------------- S3: lattice reciprocity, N2 via Newton cube root
rhs2 = lambda n: Fr((3 * n - 5) * (n - 2) * (n - 1), 4)
t2 = []
uni = list(range(2, 261)) + [289, 300, 343, 361, 400]
for n in uni:
    Y = (n - 1) ** 3 // n
    chk(Y == (n - 1) * (n - 2), "Y n=%d" % n)
    N1 = fsum(lambda k: k ** 3, n)
    N2 = sum(icbrt(n * y) for y in range(1, Y + 1))
    Z3 = sum(1 for k in range(1, n) if k ** 3 % n == 0)
    chk(N1 + N2 == (n - 1) ** 2 * (n - 2) + Z3, "recip n=%d" % n)
    chk(Fr(N2) == rhs2(n) + Fr(Z3, 2), "N2 closed n=%d" % n)
    if Fr(N2) == rhs2(n): t2.append(n)
chk(t2 == [n for n in uni if squarefree(n)], "truth(2)")
print("S3 reciprocity + closed form + truth(2)==squarefree on [2,260]+{289,300,343,361,400}")
n = 8; N1 = fsum(lambda k: k ** 3, 8); N2 = sum(icbrt(8 * y) for y in range(1, 43))
chk((N1, N2) == (96, 201), "n=8 table")
print("S3 n=8: N1=96 N2=201; user's (2) at p=5: %d == 30" % sum(icbrt(5 * y) for y in range(1, 13)))
chk(sum(icbrt(5 * y) for y in range(1, 13)) == 30, "p=5")
# user's wording upper limit (p-1)(p-2) at p=2 gives empty sum: RHS (3*2-5)(0)(1)/4 = 0
chk(rhs2(2) == 0 and rhs2(3) == Fr(4, 4) * 2 / 2 * 1, "p=2,3 (2)")
print("S3 user's (2) at p=2: empty sum = RHS 0; p=3: sum_{y<=2} icbrt(3y) = %d, RHS=%s" % (sum(icbrt(3 * y) for y in (1, 2)), rhs2(3)))
chk(Fr(sum(icbrt(3 * y) for y in (1, 2))) == rhs2(3), "p=3 (2)")

# -------------------------------------------------- S4: bilinear sum, Pillai
P = lambda n: sum(gcd(k, n) for k in range(1, n + 1))
t3 = []
for n in range(2, 301):
    S = sum((i * j) % n for i in range(1, n) for j in range(1, n))
    N3 = sum((i * j) // n for i in range(1, n) for j in range(1, n))
    chk(2 * S == n * (n * n - P(n)), "S n=%d" % n)
    chk(Fr(N3) == Fr((n - 1) ** 2 * (n - 2), 4) + Fr(P(n) - 2 * n + 1, 2), "N3 n=%d" % n)
    if Fr(N3) == Fr((n - 1) ** 2 * (n - 2), 4): t3.append(n)
chk(t3 == [n for n in range(2, 301) if isprime(n)], "truth(3)")
print("S4 S(n), N3(n) closed forms n<=300; truth(3)==primes")
for n in range(2, 5001):
    chk(P(n) >= 2 * n - 1 and ((P(n) == 2 * n - 1) == isprime(n)), "P ineq n=%d" % n)
print("S4 P(n)>=2n-1 with equality iff prime, n<=5000")
tab = {4: (8, 5), 6: (15, 27), 8: (20, 76), 12: (40, 311), 5: (9, 12)}
for n, (Pn, N3n) in tab.items():
    chk(P(n) == Pn and sum((i * j) // n for i in range(1, n) for j in range(1, n)) == N3n, "table n=%d" % n)
print("S4 table values P,N3 at n=4,5,6,8,12 confirmed")
for n in range(2, 120):
    N1 = fsum(lambda k: k ** 3, n); Y = (n - 1) * (n - 2)
    N2 = sum(icbrt(n * y) for y in range(1, Y + 1))
    N3 = sum((i * j) // n for i in range(1, n) for j in range(1, n))
    Z3 = sum(1 for k in range(1, n) if k ** 3 % n == 0)
    chk(4 * N3 - (N1 + N2) == 2 * (P(n) - 2 * n + 1) - Z3, "S4b n=%d" % n)
print("S4b 4N3-(N1+N2)=2(P-2n+1)-Z3 for n<120")

# -------------------------------------------------- S5: class numbers, independently
def h_forms(D):
    """count reduced primitive forms (a,b,c), D=b^2-4ac<0: |b|<=a<=c, b>=0 if |b|==a or a==c"""
    cnt = 0
    for a in range(1, int(((-D) / 3) ** 0.5) + 2):
        for b in range(-a, a + 1):
            if (b * b - D) % (4 * a): continue
            c = (b * b - D) // (4 * a)
            if c < a: continue
            if (abs(b) == a or a == c) and b < 0: continue
            if gcd(gcd(a, abs(b)), c) != 1: continue
            cnt += 1
    return cnt
# literature table (e.g. Cohen, A Course in Comp. Alg. NT, table of h(-p)): p = 3 mod 4, p<200
LIT = {3: 1, 7: 1, 11: 1, 19: 1, 23: 3, 31: 3, 43: 1, 47: 5, 59: 3, 67: 1, 71: 7, 79: 5, 83: 3,
       103: 5, 107: 3, 127: 5, 131: 5, 139: 3, 151: 7, 163: 1, 167: 11, 179: 5, 191: 13, 199: 9}
for p, h in LIT.items():
    chk(h_forms(-p) == h, "lit h(-%d)" % p)
print("S5 reduced-form class numbers match literature table for all p=3 mod 4, p<200")
def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)
primes = [p for p in range(3, 500) if isprime(p)]
for p in primes:
    S2 = fsum(lambda k: k * k, p)
    dev = Fr(S2) - (Fr((p - 1) * (2 * p - 1), 6) - Fr(p - 1, 2))
    if p % 4 == 3:
        h = h_forms(-p)
        w = 6 if p == 3 else 2
        chk(Fr(-sum(leg(a, p) * a for a in range(1, p)) * w, 2 * p) == h, "Dirichlet p=%d" % p)
        chk(dev == Fr(2 * h, w), "dev p=%d %s h=%d" % (p, dev, h))
    else:
        chk(dev == 0, "p=1 mod 4 p=%d" % p)
print("S5 sum floor(k^2/p) = law + 2h(-p)/w for all p=3 mod 4 <500 (p=3: 1/3); 0 for p=1 mod 4")
# hostile p=2: sum_{k=1}^{1} floor(1/2)=0, law = 1/2-1/2 = 0
chk(fsum(lambda k: k * k, 2) == 0, "p=2")
# general even e: the explorer's two verified cases, plus the UNVERIFIED 'only if' direction
counter = []
for p in primes:
    for e in (2, 4, 6, 8, 10, 12):
        g = gcd(e, p - 1)
        dev = Fr(fsum(lambda k: k ** e, p)) - (Fr(sum(k ** e for k in range(1, p)), p) - Fr(p - 1, 2))
        if ((p - 1) // g) % 2 == 0:
            chk(dev == 0, "even case p=%d e=%d" % (p, e))
        elif g == 2 and p > 3:
            chk(dev == h_forms(-p), "g=2 p=%d e=%d" % (p, e))
        elif p % 4 == 3 and p > 3:
            # g in {6,10,...}, (p-1)/g odd: is the deviation ever h(-p)?  (tests the 'iff')
            if dev == h_forms(-p):
                counter.append((p, e, g, dev))
print("S5 even e in {2..12}: both stated cases hold p<500; 'deviation = h(-p) ONLY IF g=2' counterexamples with g>2, p=3 mod 4:", counter)
ex = [(p, e, gcd(e, p - 1), Fr(fsum(lambda k: k ** e, p)) - (Fr(sum(k ** e for k in range(1, p)), p) - Fr(p - 1, 2)), h_forms(-p))
      for p, e in ((7, 6), (19, 6), (31, 6), (11, 10), (43, 6))]
print("   samples (p,e,g,D_e,h(-p)) for g>2, (p-1)/g odd:", ex)

# -------------------------------------------------- quartic deviation: exact character-sum form (extra finding)
import cmath
quart = []
for p in primes:
    if p % 8 != 5: continue
    dev = Fr(fsum(lambda k: k ** 4, p)) - (Fr(sum(k ** 4 for k in range(1, p)), p) - Fr(p - 1, 2))
    # primitive root, index map
    gen = next(g for g in range(2, p) if all(pow(g, (p - 1) // q, p) != 1 for q in factor(p - 1)))
    ind = {}; x = 1
    for i in range(p - 1):
        ind[x] = i; x = x * gen % p
    # Re sum a chi4(a): chi4(a)=i^{ind a}; Re = +a (ind=0 mod 4), -a (ind=2 mod 4), 0 else
    re = sum(a for a in range(1, p) if ind[a] % 4 == 0) - sum(a for a in range(1, p) if ind[a] % 4 == 2)
    im = sum(a for a in range(1, p) if ind[a] % 4 == 1) - sum(a for a in range(1, p) if ind[a] % 4 == 3)
    chk(dev == Fr(-2 * re, p), "D4 = -(2/p) Re sum a chi4(a) p=%d" % p)
    B2 = Fr(re * re + im * im, p * p)   # |B_{1,chi4}|^2
    quart.append((p, dev, Fr(-2 * re, p), B2))
print("S5-extra PROVED-in-range: D_4(p) = -(2/p) Re sum_a a chi_4(a) = -2 Re B_{1,chi_4} for all p=5 mod 8 <500")
print("   (p, D_4, |B_{1,chi4}|^2) first rows:", [(p, str(d), str(b)) for p, d, _, b in quart[:8]])
print("   |B|^2/2 integer? (would be relative class number h^- of the imaginary cyclic quartic field of conductor p, Q=1,w=2):",
      [(p, str(b / 2)) for p, _, _, b in quart[:10]])

# -------------------------------------------------- S6a: negation conjugacy + degenerate 3|k
def T(k, n):
    m = 3 * n + k
    if m == 0: return None
    while m % 2 == 0: m //= 2
    return m
for k in (1, -1, 5, -5, 7, -7, 11, -11, 13, -13, 3, 9):
    for n in range(-4001, 4002, 2):
        a, b = T(-k, -n), T(k, n)
        if a is None or b is None:
            chk(a is None and b is None and 3 * n + k == 0, "degenerate")
            print("   S6a degenerate: k=%d n=%d has 3n+k=0, T undefined (v_2(0) undefined); statement needs 3 !| k or n != -k/3" % (k, n))
            continue
        chk(a == -b, "conj k=%d n=%d" % (k, n))
print("S6a T_{-k}(-n) = -T_k(n) for k in {+-1,+-5,+-7,+-11,+-13,3,9}, odd |n|<=4001 (except the 3n+k=0 point)")
chk(T(-1, 5) == 7 and T(-1, 7) == 5 and T(1, -5) == -7 and T(1, -7) == -5, "(5,7)")
print("S6a (5,7) cycle of 3n-1 <-> (-5,-7) of 3n+1 confirmed")

# -------------------------------------------------- F3 cycle counts: finite signed search (bounds only)
def cycles(k, B, steps=4000):
    found = set()
    for n0 in range(-B, B + 1, 2):
        n = n0; seen = {}
        for s in range(steps):
            if n in seen: break
            seen[n] = s
            n = T(k, n)
            if n is None or abs(n) > 10 ** 12: break
        else:
            continue
        if n is None or n not in seen: continue
        cyc = [n]; m = T(k, n)
        while m != n: cyc.append(m); m = T(k, m)
        found.add(min(cyc))
    return sorted(found)
for k, claimed in ((1, 4), (5, 9), (7, 5), (11, 7), (13, 13)):
    c = cycles(k, 20001)
    print("   F3 finite search 3n+%d, odd |n0|<=20001: %d cycles (claimed %d): mins %s" % (k, len(c), claimed, c))

print("elapsed %.1fs" % (time.time() - T0))
print("AUDIT ALL CHECKS PASSED")
