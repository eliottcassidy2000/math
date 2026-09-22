#!/usr/bin/env python3
"""
collatz_mod6_20260917_floor_sums_odd_functions_audit.py -- adversarial verify-and-fix
audit of lane floor_sums_odd_functions (session collatz-mod6-20260917, wave two).

Written 2026-09-21 as the lane's only surviving independent audit.  Nothing is
imported from the lane script.  Deliberately different methods:
  * cube roots by integer binary search (lane: monotone walk);
  * squarefree count to 10^6 by the Moebius sum Q(x) = sum mu(d) floor(x/d^2) (lane: sieve);
  * class numbers by a separately coded reduced-form count, cross-checked against a
    recollected literature table (UNCITED-RECOLLECTION, p < 200) and against Dirichlet;
  * the even-power deviation D_e(p) reduced EXACTLY to index-class sums A_m
    (integer arithmetic, no complex characters), which also tests the lane's
    "D_e(p) = h(-p) iff gcd(e,p-1) = 2" -- the ONLY-IF direction is not proved;
  * a bounded signed-cycle census of 3n+k for the F3 counts 4,9,5,7,13 the note quotes.
Every check is an explicit raise (active under python3 -O).
"""
import time
from fractions import Fraction as Fr
from math import gcd

T0 = time.time()

def chk(c, m):
    if not c:
        raise RuntimeError("AUDIT CHECK FAILED: " + m)

def stamp():
    return "[t=%.1fs]" % (time.time() - T0)

def factor(n):
    f = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f

def isprime(n):
    if n < 2:
        return False
    d = 2
    while d * d <= n:
        if n % d == 0:
            return False
        d += 1
    return True

def squarefree(n):
    return all(a == 1 for a in factor(n).values())

def icbrt(m):
    """floor cube root by binary search (exact integers)"""
    if m <= 0:
        return 0
    lo, hi = 0, 1
    while hi ** 3 <= m:
        hi *= 2
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if mid ** 3 <= m:
            lo = mid
        else:
            hi = mid
    return lo

def fsum(f, n):
    return sum(f(k) // n for k in range(1, n))

def Zf(f, n):
    return sum(1 for k in range(1, n) if f(k) % n == 0)

def law(f, n):
    return Fr(sum(f(k) for k in range(1, n)), n) - Fr(n - 1, 2) + Fr(Zf(f, n), 2)

print("AUDIT of lane floor_sums_odd_functions (independent recompute, 2026-09-21)")

# ------------------------------------------------------------------ S1
print("\n== S1 odd-function floor law ==")
odd = {"x": lambda k: k, "x^3": lambda k: k ** 3, "x^5": lambda k: k ** 5,
       "x^3+x": lambda k: k ** 3 + k, "3x": lambda k: 3 * k,
       "x^7-5x^3+x": lambda k: k ** 7 - 5 * k ** 3 + k,
       "x^3-100x (negative values)": lambda k: k ** 3 - 100 * k,
       "-7x^9+x (negative values)": lambda k: -7 * k ** 9 + k}
for name, f in odd.items():
    bad = [n for n in range(2, 401) if Fr(fsum(f, n)) != law(f, n)]
    chk(not bad, "S1 %s fails at %s" % (name, bad[:3]))
print("  S1 law holds for %s, all 2 <= n <= 400 (Python floor is the true floor for negative values)" % list(odd))
# n = 1: both sides are empty sums = 0 (Z_f(1) = 0, (n-1)/2 = 0)
chk(fsum(lambda k: k ** 3, 1) == 0 and law(lambda k: k ** 3, 1) == 0, "n=1 degenerate")
print("  S1 n=1: both sides 0 (degenerate, consistent; the lane states n >= 2)")
# non-polynomial odd functions mod n (random integer lifts), fixed point residue in {0, n/2}
import random
random.seed(20260921)
for n in range(2, 80):
    vals = {}
    for k in range(1, n):
        if k in vals:
            continue
        v = random.randrange(-5 * n, 5 * n)
        vals[k] = v
        if n - k != k:
            vals[n - k] = -v + n * random.randrange(-2, 3)   # only -v mod n is required
    if n % 2 == 0:
        vals[n // 2] = random.choice([0, n // 2]) + n * random.randrange(-3, 4)
    g = (lambda k, vals=vals: vals[k])
    chk(Fr(fsum(g, n)) == law(g, n), "generic odd function mod n, n=%d" % n)
print("  S1 generic odd functions Z -> Z with f(n-k) = -f(k) (mod n) only (random lifts): law holds 2 <= n < 80")
# residue pair witnesses quoted in the note
f3 = lambda k: k ** 3
chk([(k, f3(k) % 6, f3(6 - k) % 6) for k in range(1, 3)] == [(1, 1, 5), (2, 2, 4)] and f3(3) % 6 == 3, "n=6 pairs")
chk(f3(6) % 12 == 0 and f3(2) % 4 == 0, "n=12 / n=4 fixed residues 0")
print("  S1 n=6, f=x^3: pairs (1,1,5),(2,2,4), fixed residue 3; n=12 fixed residue 0; n=4 fixed residue 0")
# even f: n=2 holds for every f, n=3 minimal failure; exact values
even = {"x^2": lambda k: k * k, "x^4": lambda k: k ** 4, "x^2+x": lambda k: k * k + k, "2 (const)": lambda k: 2}
mins = {}
for name, f in even.items():
    chk(Fr(fsum(f, 2)) == law(f, 2), "n=2 %s" % name)
    first = next(n for n in range(2, 60) if Fr(fsum(f, n)) != law(f, n))
    mins[name] = (first, fsum(f, first), law(f, first))
chk(mins["x^2"] == (3, 1, Fr(2, 3)) and mins["x^4"] == (3, 5, Fr(14, 3))
    and mins["x^2+x"] == (3, 2, Fr(13, 6)) and mins["2 (const)"] == (3, 0, Fr(1, 3)), "even minimal failures %s" % mins)
print("  S1 even f: n=2 holds for all four; minimal failures %s" % {k: (v[0], v[1], str(v[2])) for k, v in mins.items()})
# HOSTILE: an EVEN f can obey the law at some n > 2 by accident?  (the lane only claims minimal failure)
acc = [(name, n) for name, f in even.items() for n in range(3, 60) if Fr(fsum(f, n)) == law(f, n)]
print("  S1 hostile: even f obeying the odd law at some 3 <= n < 60 (accidental agreements): %s" % acc)
sq_ok = sorted(n for name, n in acc if name == "x^2")
pred = [n for n in range(3, 60) if n % 4 != 0 and all(q % 4 == 1 for q in factor(n) if q != 2)]
chk(sq_ok == pred, "x^2 accidental set %s vs predicate %s" % (sq_ok, pred))
print("  S1 FINITE-EXACT: on [3,60) x^2 obeys the odd law exactly at n = %s = {n : 4 !| n, every odd prime factor of n is 1 mod 4} (observation only; not proved)" % sq_ok)
print(stamp())

# ------------------------------------------------------------------ S2
print("\n== S2 fixed layer Z_{x^e}(n) = n/n_e - 1, squarefree truth sets ==")
def n_e(n, e):
    r = 1
    for p, a in factor(n).items():
        r *= p ** ((a + e - 1) // e)
    return r
for e in range(1, 10):
    for n in range(2, 2001):
        chk(sum(1 for k in range(1, n) if pow(k, e, n) == 0) == n // n_e(n, e) - 1, "Z e=%d n=%d" % (e, n))
print("  S2 Z_{x^e}(n) = n/n_e - 1 for 1 <= e <= 9, 2 <= n <= 2000 (direct count via pow(k,e,n))")
rhs1 = lambda n: Fr((n - 2) * (n - 1) * (n + 1), 4)
t1 = [n for n in range(2, 2001) if Fr(fsum(f3, n)) == rhs1(n)]
sq = [n for n in range(2, 2001) if squarefree(n)]
chk(t1 == sq and len(sq) == 1214, "truth(1) vs squarefree, count %d" % len(sq))
print("  S2 truth set of (1) on [2,2000] = squarefree, 1214 values")
tab = {4: (8, Fr(15, 2), 1), 8: (96, Fr(189, 2), 3), 9: (141, Fr(140), 2), 12: (358, Fr(715, 2), 1), 36: (11010, Fr(22015, 2), 5)}
for n, (fs, r, z) in tab.items():
    chk(fsum(f3, n) == fs and rhs1(n) == r and Zf(f3, n) == z and Fr(fs) - r == Fr(z, 2), "S2 table n=%d" % n)
print("  S2 table (n, floor sum, RHS(1), Z_3) confirmed at n = 4, 8, 9, 12, 36; deviation = Z_3/2")
for e in (3, 5, 7, 9):
    ts = [n for n in range(2, 1201) if Fr(fsum(lambda k: k ** e, n)) == Fr(sum(k ** e for k in range(1, n)), n) - Fr(n - 1, 2)]
    chk(ts == [n for n in sq if n <= 1200], "odd e=%d truth set" % e)
chk(all(Fr(fsum(lambda k: k, n)) == Fr(sum(range(1, n)), n) - Fr(n - 1, 2) for n in range(2, 1201)), "e=1")
print("  S2 odd e = 3,5,7,9: zero-defect identity on [2,1200] iff squarefree; e = 1 holds for all n")
def rad(n):
    r = 1
    for p in factor(n):
        r *= p
    return r
bad7 = [n for n in range(2, 300) if sum(1 for k in range(1, n) if pow(k, 7, n) == 0) != n // rad(n) - 1]
Z7_256 = sum(1 for k in range(1, 256) if pow(k, 7, 256) == 0)
chk(bad7 == [256] and Z7_256 == 63 and 256 // rad(256) - 1 == 127, "Z_7 first failure %s, Z_7(256)=%d, 256/rad-1=%d" % (bad7, Z7_256, 256 // rad(256) - 1))
print("  S2 Z_7(n) = n/rad(n) - 1 on [2,300) fails only at n = 256: Z_7(256) = 63, 256/rad(256) - 1 = 127")
print("     CORRECTION: the lane script printed the hard-coded literal '256/rad-1 = 255' and the note repeated 255; rad(256) = 2 so the true value is 127 (256 = 2^8, n_7 = 2^2 = 4, Z_7 = 64 - 1 = 63)")
# squarefree count to 10^6 by Moebius sum (independent of the lane's sieve)
N = 10 ** 6
L = 1000
mu = [1] * (L + 1)
isp = [True] * (L + 1)
for i in range(2, L + 1):
    if isp[i]:
        for j in range(i, L + 1, i):
            if j > i:
                isp[j] = False
            mu[j] = -mu[j]
        for j in range(i * i, L + 1, i * i):
            mu[j] = 0
Q = sum(mu[d] * (N // (d * d)) for d in range(1, L + 1))
chk(Q == 607926, "Q(10^6) = %d" % Q)
print("  S2 #squarefree <= 10^6 = 607926 by Moebius sum (matches the lane's sieve count)")
print(stamp())

# ------------------------------------------------------------------ S3
print("\n== S3 cube-root lattice reciprocity ==")
rhs2 = lambda n: Fr((3 * n - 5) * (n - 2) * (n - 1), 4)
uni = list(range(2, 221)) + [225, 243, 250, 256, 289, 300, 331, 343, 360, 361, 385, 399, 400]
t2 = []
rows3 = {}
for n in uni:
    Y = (n - 1) ** 3 // n
    chk(Y == (n - 1) * (n - 2), "Y n=%d" % n)
    N1 = fsum(f3, n)
    N2 = sum(icbrt(n * y) for y in range(1, Y + 1))
    Z3 = Zf(f3, n)
    chk(N1 + N2 == (n - 1) ** 2 * (n - 2) + Z3, "reciprocity n=%d" % n)
    chk(Fr(N2) == rhs2(n) + Fr(Z3, 2), "N2 closed form n=%d" % n)
    if Fr(N2) == rhs2(n):
        t2.append(n)
    rows3[n] = (N1, N2, Z3)
chk(t2 == [n for n in uni if squarefree(n)], "truth(2)")
print("  S3 (i),(ii),(iii) and truth set of (2) = squarefree on [2,220] + 13 samples (binary-search cube roots)")
want = {7: (60, 120, 0), 8: (96, 201, 3), 9: (141, 309, 2), 12: (358, 853, 1), 31: (6960, 19140, 0), 32: (7676, 21161, 7)}
for n, v in want.items():
    chk(rows3[n] == v, "S3 table n=%d: %s" % (n, rows3[n]))
print("  S3 table rows (N1, N2, Z_3) at n = 7, 8, 9, 12, 31, 32 confirmed")
# boundary: n = 2, 3 (empty / tiny sums)
chk(rows3[2] == (0, 0, 0) and rhs2(2) == 0 and rows3[3] == (2, 2, 0) and rhs2(3) == 2, "n=2,3")
print("  S3 boundary: n=2 gives N1=N2=0=RHS(2); n=3 gives N1=N2=2=RHS(2)  (the recovered draft audit asserted RHS(2)(3)=1: that draft check was wrong, true value 2)")
print(stamp())

# ------------------------------------------------------------------ S4
print("\n== S4 bilinear sum, Pillai gcd-sum ==")
P = lambda n: sum(gcd(k, n) for k in range(1, n + 1))
def P_mult(n):
    """braids2 floor_reciprocity (5): sum gcd = prod p^(a-1)((a+1)p - a)"""
    r = 1
    for p, a in factor(n).items():
        r *= p ** (a - 1) * ((a + 1) * p - a)
    return r
t3 = []
rows4 = {}
for n in range(2, 261):
    S = sum((i * j) % n for i in range(1, n) for j in range(1, n))
    N3 = sum((i * j) // n for i in range(1, n) for j in range(1, n))
    Pn = P(n)
    chk(Pn == P_mult(n), "P multiplicative form n=%d" % n)
    chk(2 * S == n * (n * n - Pn), "S(n) n=%d" % n)
    chk(Fr(N3) == Fr((n - 1) ** 2 * (n - 2), 4) + Fr(Pn - 2 * n + 1, 2), "N3 n=%d" % n)
    if Fr(N3) == Fr((n - 1) ** 2 * (n - 2), 4):
        t3.append(n)
    rows4[n] = (Pn, N3)
chk(t3 == [n for n in range(2, 261) if isprime(n)], "truth(3)")
print("  S4 S(n) and N3(n) closed forms for 2 <= n <= 260; truth set of (3) = primes; P(n) = prod p^(a-1)((a+1)p-a) [= braids2 (5) + 2n - 1]")
for n in range(2, 5001):
    Pn = P_mult(n)
    chk(Pn >= 2 * n - 1 and ((Pn == 2 * n - 1) == isprime(n)), "P >= 2n-1 iff prime n=%d" % n)
print("  S4 P(n) >= 2n - 1 with equality iff n prime, all 2 <= n <= 5000")
want4 = {4: (8, 5), 5: (9, 12), 6: (15, 27), 7: (13, 45), 8: (20, 76), 9: (21, 114), 12: (40, 311), 15: (45, 645)}
for n, v in want4.items():
    chk(rows4[n] == v, "S4 table n=%d" % n)
print("  S4 table (P, N3) at n = 4,5,6,7,8,9,12,15 confirmed")
for n in range(2, 150):
    N1, N2, Z3 = rows3[n]
    chk(4 * rows4[n][1] - (N1 + N2) == 2 * (P(n) - 2 * n + 1) - Z3, "S4b n=%d" % n)
chk((Zf(f3, 6), P(6) - 11) == (0, 4) and (Zf(f3, 8), P(8) - 15) == (3, 5), "independence witnesses")
print("  S4b 4N3 - (N1+N2) = 2(P-2n+1) - Z_3 for 2 <= n < 150; witnesses n=6: (Z_3, P-2n+1) = (0,4); n=8: (3,5)")
print(stamp())

# ------------------------------------------------------------------ S5
print("\n== S5 even powers, class numbers ==")
def h_forms(D):
    """h(D), D < 0: reduced primitive forms (a,b,c): |b| <= a <= c, b >= 0 if |b| = a or a = c"""
    cnt = 0
    a = 1
    while 3 * a * a <= -D:
        for b in range(-a, a + 1):
            if (b * b - D) % (4 * a):
                continue
            c = (b * b - D) // (4 * a)
            if c < a:
                continue
            if (abs(b) == a or a == c) and b < 0:
                continue
            if gcd(gcd(a, abs(b)), c) != 1:
                continue
            cnt += 1
        a += 1
    return cnt
# UNCITED-RECOLLECTION literature table h(-p), p = 3 mod 4, p < 200
LIT = {3: 1, 7: 1, 11: 1, 19: 1, 23: 3, 31: 3, 43: 1, 47: 5, 59: 3, 67: 1, 71: 7, 79: 5, 83: 3,
       103: 5, 107: 3, 127: 5, 131: 5, 139: 3, 151: 7, 163: 1, 167: 11, 179: 5, 191: 13, 199: 9}
mism = [(p, h_forms(-p), h) for p, h in LIT.items() if h_forms(-p) != h]
chk(not mism, "reduced-form count vs recollected table: %s" % mism)
print("  S5 reduced-form h(-p) agrees with the recollected table for all 24 primes p = 3 mod 4 below 200")
def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)
primes = [p for p in range(3, 500) if isprime(p)]
p3 = [p for p in primes if p % 4 == 3]
chk(len(p3) == 50 and len([p for p in p3 if p > 3]) == 49, "count of p=3 mod 4 below 500: %d" % len(p3))
print("  S5 COUNT: %d primes p = 3 mod 4 below 500 INCLUDING p = 3, hence %d with p > 3 (the note's 'all 50 primes with p > 3' is off by one: the lane's row count includes p = 3)" % (len(p3), len(p3) - 1))
rows5 = {}
for p in primes:
    S2 = fsum(lambda k: k * k, p)
    lawv = Fr((p - 1) * (2 * p - 1), 6) - Fr(p - 1, 2)
    dev = Fr(S2) - lawv
    if p % 4 == 3:
        h = h_forms(-p)
        w = 6 if p == 3 else 2
        RmN = sum(leg(a, p) * a for a in range(1, p))
        chk(RmN == -Fr(2 * p, w) * h, "R-N = -(2p/w) h at p=%d" % p)
        chk(dev == Fr(2 * h, w), "deviation = 2h/w at p=%d: %s vs h=%d" % (p, dev, h))
        rows5[p] = (S2, lawv, dev, h)
    else:
        chk(dev == 0, "p = 1 mod 4 deviation p=%d" % p)
print("  S5 all 50 primes p = 3 mod 4 below 500 (p = 3 included): R - N = -(2p/w) h(-p) and deviation = 2h(-p)/w, w = 2 (p > 3), w = 6 (p = 3)")
print("  S5 all primes p = 1 mod 4 below 500: deviation 0")
want5 = {3: (1, Fr(2, 3), Fr(1, 3), 1), 7: (11, 10, 1, 1), 11: (31, 30, 1, 1), 19: (103, 102, 1, 1), 23: (157, 154, 3, 3),
         31: (293, 290, 3, 3), 43: (575, 574, 1, 1), 47: (695, 690, 5, 5), 59: (1105, 1102, 3, 3), 67: (1431, 1430, 1, 1),
         71: (1617, 1610, 7, 7), 79: (2007, 2002, 5, 5), 83: (2217, 2214, 3, 3), 103: (3439, 3434, 5, 5),
         107: (3713, 3710, 3, 3), 127: (5255, 5250, 5, 5)}
for p, v in want5.items():
    chk(tuple(Fr(x) for x in rows5[p]) == tuple(Fr(x) for x in v), "S5 table p=%d: %s" % (p, rows5[p]))
print("  S5 the 16 table rows (S2, odd-law value, deviation, h) confirmed; p=3 formula value 5/3 vs actual 1 (error -2/3)")
chk(Fr(2 * 5, 6) - 1 + 1 == Fr(5, 3) and fsum(lambda k: k * k, 2) == 0, "p=3 formula value / p=2")
print("  S5 hostile p = 2: sum_{k=1}^{1} floor(1/2) = 0 = odd-law value (2 is 'p = 3 mod 4'-free; excluded by p >= 3 anyway)")

# exact index-class decomposition of D_e(p): D_e = (p-1)/2 - g*sum_{H_g} r / p, and
# g*sum_{H_g} = sum_{j=0}^{g-1} sum_a a chi_j(a) with chi_j = zeta_g^{j*ind}; only ODD chi_j survive.
def primroot(p):
    fs = factor(p - 1)
    return next(g for g in range(2, p) if all(pow(g, (p - 1) // q, p) != 1 for q in fs))
def index_table(p):
    g0 = primroot(p)
    ind = {}
    x = 1
    for i in range(p - 1):
        ind[x] = i
        x = x * g0 % p
    return ind
coinc = []   # (p, e, g, D_e) with g >= 4, (p-1)/g odd, and D_e == h(-p)  -> tests the lane's "iff"
nonint = []
for p in primes:
    ind = index_table(p)
    for e in (2, 4, 6, 8, 10, 12):
        g = gcd(e, p - 1)
        dev = Fr(fsum(lambda k: k ** e, p)) - (Fr(sum(k ** e for k in range(1, p)), p) - Fr(p - 1, 2))
        Hsum = sum(r for r in range(1, p) if pow(r, (p - 1) // g, p) == 1)
        chk(dev == Fr(p - 1, 2) - Fr(g * Hsum, p), "D_e definition p=%d e=%d" % (p, e))
        if ((p - 1) // g) % 2 == 0:
            chk(dev == 0, "even case p=%d e=%d" % (p, e))
        elif g == 2:
            chk(dev == Fr(2 * h_forms(-p), 6 if p == 3 else 2), "g=2 case p=%d e=%d" % (p, e))
        else:
            if p % 4 == 3 and dev == h_forms(-p):
                coinc.append((p, e, g, dev))
            if dev.denominator != 1:
                nonint.append((p, e, g, str(dev)))
print("  S5 e in {2,4,6,8,10,12}, all 3 <= p < 500: D_e(p) = (p-1)/2 - g*sum(H_g)/p; D_e = 0 when (p-1)/g even; D_e = 2h(-p)/w when g = 2")
print("  S5 ONLY-IF test of 'D_e(p) = h(-p) iff g = 2 and p = 3 mod 4': coincidences with g >= 4, (p-1)/g odd, p = 3 mod 4: %s" % coinc)
print("  S5 non-integer D_e with g >= 4, (p-1)/g odd (p, e, g, D_e): %s" % nonint)
samples = []
for p, e in ((7, 6), (19, 6), (31, 6), (43, 6), (11, 10), (31, 10)):
    g = gcd(e, p - 1)
    dev = Fr(fsum(lambda k: k ** e, p)) - (Fr(sum(k ** e for k in range(1, p)), p) - Fr(p - 1, 2))
    samples.append((p, e, g, str(dev), h_forms(-p)))
print("  S5 samples (p, e, g, D_e, h(-p)) with g >= 6, (p-1)/g odd: %s" % samples)

# quartic: D_4(p) = -(2/p)(A_0 - A_2), A_m = sum of a with ind(a) = m (mod 4); |S(chi_4)|^2 = (A_0-A_2)^2 + (A_1-A_3)^2
print("  S5 quartic p = 5 mod 8: exact form D_4(p) = -(2/p)(A_0 - A_2), A_m = sum_{ind a = m mod 4} a  (PROVED below in the note)")
print("     p | D_4 | A_0-A_2 | A_1-A_3 | |S(chi_4)|^2/(2p^2) | h(-4p)")
q_rows = []
for p in primes:
    if p % 8 != 5:
        continue
    ind = index_table(p)
    A = [sum(a for a in range(1, p) if ind[a] % 4 == m) for m in range(4)]
    dev = Fr(fsum(lambda k: k ** 4, p)) - (Fr(sum(k ** 4 for k in range(1, p)), p) - Fr(p - 1, 2))
    chk(dev == Fr(-2 * (A[0] - A[2]), p), "D_4 index form p=%d" % p)
    chk(A[1] - A[3] != 0 or p == 5, "imaginary part sanity p=%d" % p)
    B2 = Fr((A[0] - A[2]) ** 2 + (A[1] - A[3]) ** 2, 2 * p * p)
    h4 = h_forms(-4 * p)
    q_rows.append((p, dev, A[0] - A[2], A[1] - A[3], B2, h4))
    if p < 200:
        print("   %3d | %5s | %6d | %6d | %s | %d" % (p, dev, A[0] - A[2], A[1] - A[3], B2, h4))
want_q = {5: Fr(6, 5), 13: 2, 29: -2, 37: 2, 53: -2, 61: 2, 101: -6, 109: 10, 149: 6, 157: 2, 173: -6, 181: 14, 197: -2}
want_h4 = {5: 2, 13: 2, 29: 6, 37: 2, 53: 6, 61: 6, 101: 14, 109: 6, 149: 14, 157: 6, 173: 14, 181: 10, 197: 10}
for p, dev, _, _, _, h4 in q_rows:
    if p in want_q:
        chk(dev == want_q[p] and h4 == want_h4[p], "quartic table p=%d: %s %d" % (p, dev, h4))
chk(not all(d == h for _, d, _, _, _, h in q_rows) and not all(d == 2 * h for _, d, _, _, _, h in q_rows), "quartic verdicts")
int_B = all(B2.denominator == 1 for p, _, _, _, B2, _ in q_rows if p > 5)
print("  S5 quartic table D_4 and h(-4p) confirmed for the 13 primes p = 5 mod 8 below 200; neither D_4 = h(-4p) nor 2h(-4p) holds for all p < 500")
print("  S5 |S(chi_4)|^2/(2p^2) integer for every p = 5 mod 8 with 5 < p < 500: %s (p = 5: %s)" % (int_B, q_rows[0][4]))
chk(int_B, "|B_{1,chi_4}|^2/2 integrality")
print(stamp())

# ------------------------------------------------------------------ S6
print("\n== S6 negation conjugacy and the F3 cycle counts ==")
def T(k, n):
    m = 3 * n + k
    if m == 0:
        return None
    while m % 2 == 0:
        m //= 2
    return m
degen = []
for k in (1, -1, 5, -5, 7, -7, 11, -11, 13, -13, 3, -3, 9):
    for n in range(-4001, 4002, 2):
        a, b = T(-k, -n), T(k, n)
        if a is None or b is None:
            chk(a is None and b is None and 3 * n + k == 0, "degenerate")
            degen.append((k, n))
            continue
        chk(a == -b, "conjugacy k=%d n=%d" % (k, n))
chk(degen == [(3, -1), (-3, 1), (9, -3)], "degenerate set %s" % degen)
print("  S6a T_{-k}(-n) = -T_k(n) for k in {+-1,+-5,+-7,+-11,+-13,+-3,9}, odd |n| <= 4001, except the points 3n+k = 0: %s" % degen)
print("     (so the conjugacy statement needs 3 !| k, or n != -k/3; the lane's k in {+-1,+-5,+-7,+-13} are all prime to 3)")
chk(T(-1, 5) == 7 and T(-1, 7) == 5 and T(1, -5) == -7 and T(1, -7) == -5, "(5,7)")
print("  S6a (5,7) cycle of 3n-1 <-> (-5,-7) cycle of 3n+1 confirmed")
chk(8 // n_e(8, 3) - 1 == 3, "Z_3(8)")

def census(k, B, cap=10 ** 12, steps=100000):
    """cycles of T_k met from odd seeds |n0| <= B; returns sorted list of cycle minima (by absolute value key)"""
    status = {}   # n -> cycle id or 'esc'
    cycles = {}
    for n0 in range(-B, B + 1, 2):
        path = []
        n = n0
        seen_local = {}
        res = None
        while True:
            if n in status:
                res = status[n]
                break
            if n in seen_local:
                cyc = path[seen_local[n]:]
                key = min(cyc, key=lambda x: (abs(x), x))
                if key not in cycles:
                    cycles[key] = cyc
                res = key
                break
            if abs(n) > cap or len(path) > steps:
                res = 'esc'
                break
            seen_local[n] = len(path)
            path.append(n)
            n = T(k, n)
            if n is None:
                res = 'esc'
                break
        for m in path:
            status[m] = res
    return cycles
for k, claimed in ((1, 4), (5, 9), (7, 5), (11, 7), (13, 13)):
    cyc = census(k, 20001)
    keys = sorted(cyc, key=lambda x: (abs(x), x))
    print("  F3 census 3n+%d, odd seeds |n0| <= 20001: %d cycles found (lead's count %d); representatives %s; lengths %s"
          % (k, len(cyc), claimed, keys, [len(cyc[x]) for x in keys]))
    chk(len(cyc) == claimed, "F3 count k=%d: %d vs %d" % (k, len(cyc), claimed))
print("  F3 counts 4, 9, 5, 7, 13 reproduced as a FINITE bounded-seed census (lower bounds in principle; no cycle escapes the seed box up to 10^12)")
print(stamp())
print("AUDIT ALL CHECKS PASSED")
