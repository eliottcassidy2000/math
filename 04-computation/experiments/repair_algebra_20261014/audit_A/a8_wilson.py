#!/usr/bin/env python3
"""Audit A (independent), THM-4641 rows 1-11 and 2b.  Own implementations; sympy only for isprime.
Usage: python3 a8_wilson.py"""
import math, time, itertools
from fractions import Fraction as Fr
from sympy import isprime

t0 = time.time()
LIM = 10 ** 6
# smallest-prime-factor sieve
spf = list(range(LIM + 1))
for i in range(2, int(LIM ** 0.5) + 1):
    if spf[i] == i:
        for j in range(i * i, LIM + 1, i):
            if spf[j] == j:
                spf[j] = i
primes = [p for p in range(2, LIM) if spf[p] == p]


def factor(n):
    f = {}
    while n > 1:
        p = spf[n]
        while n % p == 0:
            n //= p; f[p] = f.get(p, 0) + 1
    return f


def legendre2(p):
    return 1 if p % 8 in (1, 7) else -1


def wq(p):           # Wilson quotient mod p
    m = p * p; f = 1
    for i in range(2, p):
        f = f * i % m
    assert (f + 1) % p == 0
    return (f + 1) // p % p


def fq(a, p):        # Fermat quotient mod p
    return (pow(a, p - 1, p * p) - 1) // p % p


# ---------------- row 1: Gersonides
G = sorted((a, b) for a in range(1, 300) for b in range(0, 200) if abs(2 ** a - 3 ** b) == 1)
print("row 1: |2^a - 3^b| = 1, 1 <= a < 300, b < 200:", G)

# ---------------- row 3: product of <2> mod p
for p in primes[1:400]:
    d = 1; x = 2 % p
    while x != 1:
        x = x * 2 % p; d += 1
    prod = 1; x = 1
    for _ in range(d):
        prod = prod * x % p; x = x * 2 % p
    assert prod == (p - 1 if d % 2 == 0 else 1)
print("row 3: prod_{j<d} 2^j = -1 iff ord_p(2) even, first 399 odd primes: OK")

# ---------------- rows 4, 5 (part 1): Lerch, log property, Euler quotient
for p in primes[1:400]:
    w = wq(p)
    assert w == sum(fq(a, p) for a in range(1, p)) % p
    q2 = fq(2, p)
    for j in range(1, 8):
        assert fq(2 ** j, p) == j * q2 % p
    fp = ((pow(2, (p - 1) // 2, p * p) - legendre2(p)) % (p * p)) // p % p
    assert fp == legendre2(p) * q2 * pow(2, -1, p) % p
print("rows 4-5: Lerch, q_p(2^j) = j q_p(2) (j <= 7), f_p = (2/p) q_p(2)/2: first 399 odd primes OK", flush=True)


# ---------------- row 5: tournaments, A_p computed directly from the non-trivial partitions
def odd_partitions(n, mx=None):
    if mx is None:
        mx = n
    if n == 0:
        yield []
        return
    for k in range(min(n, mx), 0, -1):
        if k % 2:
            for rest in odd_partitions(n - k, k):
                yield [k] + rest


def z_and_e(lam):
    z = 1
    for k in set(lam):
        c = lam.count(k)
        z *= k ** c * math.factorial(c)
    e = sum((k - 1) // 2 for k in lam) + sum(math.gcd(lam[i], lam[j]) for i in range(len(lam)) for j in range(i + 1, len(lam)))
    return z, e


def T_count(n):
    tot = Fr(0)
    for lam in odd_partitions(n):
        z, e = z_and_e(lam)
        tot += Fr(2 ** e, z)
    assert tot.denominator == 1
    return tot.numerator


A000568 = [1, 1, 1, 2, 4, 12, 56, 456, 6880, 191536, 9733056, 903753248, 154108311168, 48542114686912,
           28401423719122304]
for n in range(len(A000568)):
    assert T_count(n) == A000568[n], n
table307 = {3: (1, 1, 0, 2), 5: (1, 0, 4, 2), 7: (1, 5, 2, 1), 11: (3, 1, 2, 2), 13: (5, 0, 3, 2)}
for p in [p for p in primes if 3 <= p <= 61]:
    fact = math.factorial(p)
    Asum = 0
    for lam in odd_partitions(p):
        if lam == [p] or lam == [1] * p:
            continue
        z, e = z_and_e(lam)
        assert fact % z == 0 and (fact // z) % p == 0
        Asum += (fact // z) * 2 ** e
    assert Asum % p == 0
    Ap = Asum // p % p
    T = T_count(p)
    chi = legendre2(p); q2 = fq(2, p); w = wq(p)
    fp = ((pow(2, (p - 1) // 2, p * p) - chi) % (p * p)) // p % p
    assert T % p == (fp - chi * w - Ap) % p, p
    assert T % p == (chi * (q2 * pow(2, -1, p) - w) - Ap) % p, p
    if p in table307:
        assert (fp, w, Ap, T % p) == table307[p], (p, fp, w, Ap, T % p)
print("row 5: T(n) = A000568 for n <= 14 (own Burnside); for every prime 3 <= p <= 61 with A_p summed directly over "
      "the non-trivial odd partitions: T(p) = f_p - (2/p) w_p - A_p = chi_8(p)(q_p(2)/2 - w_p) - A_p (mod p); "
      "THM-307 table rows p = 3..13 reproduced", flush=True)

# ---------------- row 6: anchor row 2^a - 3
row = [(a, 2 ** a - 3) for a in range(3, 45) if isprime(2 ** a - 3)]
print("row 6: primes 2^a - 3, 3 <= a <= 44:", row)
wil = [(a, p) for a, p in row if a <= 24 and wq(p) == 0]
print(f"       Wilson among a <= 24: {wil}; number of primes 2^a - 3 with a <= 24 other than 5, 13: "
      f"{sum(1 for a, p in row if a <= 24) - 2}; with a <= 44: {len(row) - 2}; 2^44 - 3 < 2e13: {2**44 - 3 < 2 * 10**13}"
      f"  [{time.time()-t0:.0f}s]", flush=True)

# ---------------- row 7: 2 primitive root
def is_prim2(p):
    for q in factor(p - 1):
        if pow(2, (p - 1) // q, p) == 1:
            return False
    return True


assert all(is_prim2(p) for p in (5, 13, 563))
ps = [p for p in primes if p >= 5]
fr = sum(is_prim2(p) for p in ps) / len(ps)
cls = [p for p in ps if p % 8 in (3, 5)]
frc = sum(is_prim2(p) for p in cls) / len(cls)
phi = lambda n: n * math.prod(1 - Fr(1, q) for q in factor(n))
print(f"row 7: 2 primitive mod 5, 13, 563; frequency among primes 5 <= p < 1e6: {fr:.4f} (cube {fr**3:.4f}); "
      f"among p = 3, 5 mod 8: {frc:.4f} (cube {frc**3:.4f}); phi(p-1)/(p-1) product for 5, 13, 563: "
      f"{float(phi(4)/4 * phi(12)/12 * phi(562)/562):.4f}", flush=True)


# ---------------- row 8: p - 2 prime or Carmichael
def korselt(n):
    if n < 3:
        return False
    f = factor(n)
    return len(f) >= 2 and all(e == 1 for e in f.values()) and all((n - 1) % (q - 1) == 0 for q in f)


absF = lambda n: (spf[n] == n) or korselt(n)
cnt = sum(1 for p in ps if absF(p - 2))
car = [c for c in range(3, LIM) if korselt(c)]
carp = [c for c in car if c + 2 < LIM and spf[c + 2] == c + 2]
print(f"row 8: p - 2 prime or Carmichael for primes 5 <= p < 1e6: {cnt}/{len(ps)} = {cnt/len(ps):.4f}; Carmichael "
      f"numbers < 1e6: {len(car)}, those with c + 2 prime: {carp}", flush=True)
for (lo, hi) in ((5, 14), (5, 50), (400, 800), (300, 1000), (10 ** 5, 10 ** 6)):
    sub = [p for p in ps if lo <= p <= hi]
    print(f"       window [{lo}, {hi}]: frequency of 'p - 2 prime or Carmichael' {sum(absF(p - 2) for p in sub)}/{len(sub)}"
          f" = {sum(absF(p - 2) for p in sub)/len(sub):.3f}")
hl = [2 * 0.6601618 / math.log(p) for p in (5, 13, 563)]
print(f"       Hardy-Littlewood local rate 2 C_2 / ln p at 5, 13, 563: {[round(x, 3) for x in hl]}, product "
      f"{math.prod(hl):.3f} (vs the global-average baseline 0.104^3 = {0.104**3:.4f})", flush=True)


# ---------------- row 9: representations at -1 (complete by the positivity lemma)
def canon_reps(K, emax):
    out = []

    def back(K, ehi, tail):
        if K == -1:
            out.append(tuple(reversed(tail)))
            return
        for e in range(2, ehi + 1):
            r = K - (1 << e)
            if r % 3:
                continue
            P = r // 3
            if P == -1 or 1 <= P <= 3 ** (e - 1) - 2 ** e:
                back(P, e - 1, tail + [e])

    if K == -1:
        return [()]
    if K < 1:
        return []
    back(K, (K + 3).bit_length() - 1, [])
    return out


for K in (5, 13, 563, 561, 1093, 3511):
    print(f"row 9: canonical representations of {K}: {sorted(canon_reps(K, 64))}")

# ---------------- row 10
N = 684106401
fN = factor(N) if N <= LIM else None
fN = {}
n = N
for q in range(2, 4000):
    while n % q == 0:
        fN[q] = fN.get(q, 0) + 1; n //= q
assert n == 1
lam = 1
for q in fN:
    lam = lam * (q - 1) // math.gcd(lam, q - 1)
kors = all(e == 1 for e in fN.values()) and all((N - 1) % (q - 1) == 0 for q in fN)
print(f"row 10: {N} = {fN}, Korselt {kors}, lambda = {lam}, lambda + 1 = {lam + 1} = 27*563: {lam + 1 == 27 * 563}")

# ---------------- row 11
x, y = 563, 2e13
ll = math.log(math.log(y)) - math.log(math.log(x))
print(f"row 11: loglog(2e13) - loglog(563) = {ll:.3f}; exp(-that) = {math.exp(-ll):.3f}")


# ---------------- row 2b: unit-denominator cycles, and an own cycle census
def cyc_point(word):
    k = len(word); A = sum(word); b = 0; S = 0
    for j, a in enumerate(word):
        b += 3 ** (k - 1 - j) * 2 ** S; S += a
    return Fr(b, 2 ** A - 3 ** k), b


for (A, k) in ((1, 1), (2, 1), (3, 2)):
    pts = sorted({int(cyc_point(w)[0]) for w in itertools.product(range(1, A + 1), repeat=k) if sum(w) == A})
    print(f"row 2b: (A, k) = ({A}, {k}): cycle points {pts}")
w17 = None


def U(n):
    n = 3 * n + 1
    while n % 2 == 0:
        n //= 2
    return n


cycles = {}
for n0 in range(-100001, 100002, 2):
    n = n0; seen = {}
    for step in range(2000):
        if n in seen:
            break
        seen[n] = step; n = U(n)
    else:
        raise SystemExit(f"no cycle within 2000 U-steps from {n0}")
    c = [n]; m = U(n)
    while m != n:
        c.append(m); m = U(m)
    key = min(c, key=abs)
    if key not in cycles:
        word = []
        for v in c:
            t = 3 * v + 1; a = 0
            while t % 2 == 0:
                t //= 2; a += 1
            word.append(a)
        cycles[key] = (len(c), sum(word), tuple(word))
print(f"row 2b: own census, odd |n| <= 100001, every orbit closed within 2000 U-steps; cycles (min |.|: (k, A, word)): "
      f"{cycles}")
k17, A17, wd = cycles[-17]
bvals = [cyc_point(wd[i:] + wd[:i]) for i in range(k17)]
print(f"        -17 cycle: k = {k17}, A = {A17}, 2^A - 3^k = {2**A17 - 3**k17}; rotations' numerators b all divisible "
      f"by 139: {all(b % 139 == 0 for _, b in bvals)}; points {sorted(int(x) for x, _ in bvals)}  [{time.time()-t0:.0f}s]")
