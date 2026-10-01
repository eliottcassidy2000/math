#!/usr/bin/env python3
"""Tournament Clock Prime Collatz (TCPC): residue clocks with one-way, doubled, missing and looped pairs
(thread session for the owner's TCPC prompt, 2026-10-01).

A *clock* is a circulant digraph Cay(Z/n, D): x -> y iff y - x is in D.  Every unordered pair {x, y} has
one of four types read off the difference d = y - x:
    one-way (d in D xor -d in D), doubled (both), missing (neither); a loop at every vertex iff 0 in D.
A *power clock* uses D = D_k(n) = {x^k mod n : x in Z/n}.

Checks (dependency-free; every load-bearing check raises through FAILS):
  A. combination law under CRT: the pair type is the AND of the factor types in the monoid
     {doubled/loop = 11, forward = 10, backward = 01, missing = 00}; |D| and |D cap -D| are multiplicative;
     the power clocks that are tournaments are exactly the Paley ones (n = p = 3 mod 4, gcd(k, p-1) = 2);
     character (Jacobi) clocks multiply orientations instead: two Paley tournaments give a graph
  B. fractal law: squares mod 3^m = alternating lexicographic tower (C3, empty K3, C3, ...); cubes mod 3^m
     need two-digit blocks (coset partitions by 3^j are modules iff j != 1 mod 3); general tower law (type
     of p^v w read from v and w mod p^min(m-v, 1+v_p(k))); digit conservation (input period / image
     modulus); prime atoms: every circulant on an odd prime number of vertices is a prime digraph unless
     D minus {0} is empty or everything; a clock on n with >= 2 distinct primes is prime iff no prime-power factor has
     a clique module (finite-exact); CRT products with a tower factor are prime
  C. Collatz mod 9: 3n+1 lands in the unit squares {1,4,7}; T(n) = (-1)^v 4^(n+v) (mod 9), so T(n)^3 =
     (-1)^v and T(n) is a square mod 9 iff v is even; forward law 1/9 on squares, 2/9 on non-squares; orbit
     law (8,16,11,4,2,22)/63 on (1,2,4,5,7,8) from the second iterate on; the 3-adic digits of an orbit
     are a function of the last k halving counts; the duality <2> mod 3^m (doubled, complete 3-partite)
     versus <3> mod 2^m (oriented, bipartite tournament); 2 is a primitive root mod 9 iff 2 is neither a
     square nor a cube mod 9; for odd primes a, 2 is a primitive root mod a^2 iff 2 avoids the q-th power
     clock mod a^2 for every prime q | a(a-1) (the condition of THM-4523's Theorem R_a for n/2, an+1)
  D. triplets: Pythagorean triples are transitive triangles of the square clock; loop placement law
     (4 | even leg, 3 | a leg never the hypotenuse, 5 | any of the three, hypotenuse primes = 1 mod 4);
     unit Pythagorean triangles mod p exist iff p >= 7; leg law: the legs are exchangeable mod p iff
     p != 5 mod 8, and the even leg's factor 2 shows up as the Legendre symbol (2/q) at q = 5 mod 8
  E. F = S + U and p^2qr: F = S + U iff |Q| = |U| (Q = square-containing proper divisors); the two 3-4-3
     splits of p^2qr; the five complementary pairs; Q = lcm(U, p^2); complementation = translation by the
     squareclass of N in F_2^3; profile (2,1^(r-1)) layers (2^(r-1)-1, 2^(r-1), 2^(r-1)-1)
  F. Fermat case I as a clock triangle census: unit triangles of the p-th power clock mod p^2 <-> roots
     t != 0, -1 of ((1+t)^p - 1 - t^p)/p mod p; orbit law #roots = 2[p = 1 mod 3] + 3[Wieferich] + 6j;
     census p < P_FERMAT; 2 (resp. 3) lies in the p-th power clock mod p^2 iff p is Wieferich base 2 (3)
  G. local Waring numbers: with loops the CRT product is a strong product, so distances combine by max;
     squares 4 (at 8), cubes 4 (at 9, residues 4 and 5), fourth powers 15 (at 16)
  H. doubling clocks Cay(Z/p, <2>): complete / Paley tournament / symmetric / oriented-with-missing, read
     off the actual pair types for p < 3000 and from ord_p(2) for p < P_DOUBLING; Paley doubling primes are
     7 mod 8; the Paley share is compared with A/2 and the symmetric share with 17/24

Reproduce: python3 04-computation/experiments/tournament_clock_prime_collatz_20261001.py   (about 1 minute)
"""
import math
import sys
from collections import Counter
from itertools import combinations

P_FERMAT = 30000
P_DOUBLING = 2_000_000
FAILS = []


def check(cond, msg):
    print(("PASS " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def info(msg):
    print("   " + msg)


# ---------------------------------------------------------------- utilities
def spf_sieve(n):
    spf = list(range(n + 1))
    for i in range(2, int(n ** 0.5) + 1):
        if spf[i] == i:
            for j in range(i * i, n + 1, i):
                if spf[j] == j:
                    spf[j] = i
    return spf


SPF = spf_sieve(P_DOUBLING + 10)


def factor(n):
    f = {}
    while n > 1:
        if n < len(SPF):
            p = SPF[n]
        else:  # trial division beyond the sieve (not reached by the checks below)
            p = next((d for d in range(2, math.isqrt(n) + 1) if n % d == 0), n)
        while n % p == 0:
            f[p] = f.get(p, 0) + 1
            n //= p
    return f


def legendre(d, p):
    s = pow(d % p, (p - 1) // 2, p)
    return -1 if s == p - 1 else s


def is_prime(n):
    return n >= 2 and n < len(SPF) and SPF[n] == n


def primes_below(n):
    return [p for p in range(2, n) if SPF[p] == p]


def mult_order(a, p):
    """order of a mod prime p (a not divisible by p)"""
    o = p - 1
    for q in factor(p - 1):
        while o % q == 0 and pow(a, o // q, p) == 1:
            o //= q
    return o


def primitive_root(p):
    qs = list(factor(p - 1))
    for g in range(2, p):
        if all(pow(g, (p - 1) // q, p) != 1 for q in qs):
            return g
    return 1


def divisors(n):
    ds = [1]
    for p, e in factor(n).items():
        ds = [d * p ** i for d in ds for i in range(e + 1)]
    return sorted(ds)


def powset(n, k):
    return frozenset(pow(x, k, n) for x in range(n))


def typ(n, D, d):
    """type of the difference d: (d in D, -d in D) as a 2-bit tuple"""
    return (d % n in D, (-d) % n in D)


def counts(n, D):
    s = len(D)
    ss = sum(1 for d in D if (-d) % n in D)
    loop = 0 in D
    return dict(size=s, sym=ss, oneway=s - ss, doubled=ss - (1 if loop else 0), missing=n - 1 - 2 * s + ss + (1 if loop else 0), loop=loop)


def rel(n, D, z, m):
    return ((m - z) % n in D, (z - m) % n in D)


def closure(n, D, S):
    S = set(S)
    changed = True
    while changed:
        changed = False
        for z in range(n):
            if z in S:
                continue
            it = iter(S)
            t0 = rel(n, D, z, next(it))
            for m in it:
                if rel(n, D, z, m) != t0:
                    S.add(z)
                    changed = True
                    break
    return frozenset(S)


def is_module(n, D, M):
    M = set(M)
    return all(len({rel(n, D, z, m) for m in M}) == 1 for z in range(n) if z not in M)


def is_prime_digraph(n, D):
    return all(len(closure(n, D, {0, y})) == n for y in range(1, n))


def has_clique_module(n, D):
    for y in range(1, n):
        if y in D and (-y) % n in D:
            c = closure(n, D, {0, y})
            if all((b - a) % n in D and (a - b) % n in D for a, b in combinations(c, 2)):
                return True
    return False


def diam(n, D):
    reach, s = {0}, 0
    while len(reach) < n:
        new = {(a + d) % n for a in reach for d in D}
        if new == reach:
            return None
        reach, s = new, s + 1
    return s


def v2(x):
    c = 0
    while x % 2 == 0:
        x //= 2
        c += 1
    return c


def syr(n):
    m = 3 * n + 1
    v = v2(m)
    return m >> v, v


# ---------------------------------------------------------------- A
print("=== A. combination law under CRT, multiplicativity, which power clocks are tournaments ===")
AND = lambda s, t: (s[0] and t[0], s[1] and t[1])
ok = True
for a, b in [(9, 7), (8, 9), (5, 7), (16, 5), (27, 4), (11, 13)]:
    for k in (2, 3, 4):
        n = a * b
        D, Da, Db = powset(n, k), powset(a, k), powset(b, k)
        for d in range(n):
            if typ(n, D, d) != AND(typ(a, Da, d), typ(b, Db, d)):
                ok = False
check(ok, "pair type of a CRT product = AND of factor types (moduli 63, 72, 35, 80, 108, 143; k = 2,3,4)")
info("monoid table: 11 (doubled or loop) is the identity, 00 (missing) absorbs, 10 AND 01 = 00 (opposite arcs annihilate)")
bad = 0
for k in range(2, 7):
    for n in range(2, 400):
        c = counts(n, powset(n, k))
        ps = pss = 1
        for p, e in factor(n).items():
            cq = counts(p ** e, powset(p ** e, k))
            ps *= cq["size"]
            pss *= cq["sym"]
        bad += (c["size"], c["sym"]) != (ps, pss)
check(bad == 0, "|D_k(n)| and |D_k(n) cap -D_k(n)| are multiplicative (n < 400, k = 2..6); one-way, doubled, missing counts follow")
mism = []
for k in range(2, 9):
    for n in range(2, 3000):
        c = counts(n, powset(n, k))
        is_t = c["doubled"] == 0 and c["missing"] == 0
        pred = is_prime(n) and n % 4 == 3 and math.gcd(k, n - 1) == 2
        if is_t != pred:
            mism.append((n, k))
check(not mism, "power clock (nonzero differences) is a tournament iff n = p = 3 mod 4 prime and gcd(k, p-1) = 2: the Paley tournaments (n < 3000, k = 2..8)")
rng = __import__("random").Random(20261001)
ok = True
for _ in range(3000):
    n = rng.randrange(2, 40)
    D = frozenset(d for d in range(n) if rng.random() < rng.random())
    c = counts(n, D)
    ts = Counter(typ(n, D, d) for d in range(1, n))
    ok &= (c["oneway"], c["doubled"], c["missing"]) == (ts[(True, False)], ts[(True, True)], ts[(False, False)])
c7 = counts(7, frozenset({1, 2, 4}))
check(ok and (c7["oneway"], c7["doubled"], c7["missing"]) == (3, 0, 0),
      "per-vertex counts for any D (3000 random D, with and without 0): one-way |D| - |D cap -D|, doubled |D cap -D| - [0 in D], "
      "missing n - 1 - 2|D| + |D cap -D| + [0 in D]; loopless P7 = Cay(Z/7, {1,2,4}) has 0 missing")
info("per vertex: one-way out-arcs, doubled neighbours, missing neighbours (2*one-way + doubled + missing = n - 1)")
for n, k in [(9, 2), (9, 3), (7, 2), (7, 3), (5, 2), (8, 2), (27, 2), (63, 2)]:
    c = counts(n, powset(n, k))
    info("n=%-3d k=%d  D=%s  one-way=%d doubled=%d missing=%d loop=%s" % (n, k, sorted(powset(n, k)) if n < 30 else "...", c["oneway"], c["doubled"], c["missing"], c["loop"]))
# character clocks multiply instead: the Jacobi clock J(pq) = {d : (d/pq) = +1} has orientation (-1/p)(-1/q)
ok = True
ps60 = primes_below(60)[1:]
for p_, q_ in combinations(ps60, 2):
    n = p_ * q_
    J = {d for d in range(n) if legendre(d, p_) * legendre(d, q_) == 1}
    eps = legendre(-1, p_) * legendre(-1, q_)
    units = [d for d in range(1, n) if d % p_ and d % q_]
    if eps == 1:
        ok &= J == {(-d) % n for d in J}
    else:
        ok &= all((d in J) != ((-d) % n in J) for d in units)
check(ok, "Jacobi clocks (distinct odd primes p, q < 60): symmetric iff (-1/p)(-1/q) = +1, else a tournament on unit differences: two Paley tournaments (p, q = 3 mod 4) multiply to an undirected graph")

# ---------------------------------------------------------------- B
print("=== B. fractal towers at prime powers, digit conservation, prime atoms ===")
ok = True
for m in range(2, 6):
    n = 3 ** m
    D = powset(n, 2)
    ok &= all(is_module(n, D, set(range(r, n, 3 ** j))) for j in range(1, m) for r in range(3 ** j))
    ok &= all(typ(n, D, 3 ** j) == ((True, False) if j % 2 == 0 else (False, False)) for j in range(m))
check(ok, "squares mod 3^m (m = 2..5): every 3^j-coset partition is a module partition and level j is C3 (j even) / empty (j odd)")
ok = True
for m in range(2, 7):
    n = 3 ** m
    D = powset(n, 3)
    for j in range(1, m):
        mod = all(is_module(n, D, set(range(r, n, 3 ** j))) for r in range(3 ** j))
        pred = (j % 3 != 1)
        ok &= (mod == pred)
check(ok, "cubes mod 3^m (m = 2..6): the 3^j-coset partition is a module partition iff j != 1 (mod 3): two-digit C9 blocks then one empty digit")
check(sorted(powset(9, 3)) == [0, 1, 8] and is_prime_digraph(9, powset(9, 3)), "cube clock mod 9 = 9-cycle with loops (D = {0, 1, 8}), a prime digraph")
# digit conservation: input period and image-decision modulus
per = lambda k, n: min(t for t in range(1, n + 1) if n % t == 0 and all(pow(x + t, k, n) == pow(x, k, n) for x in range(n)))
check(per(2, 9) == 9 and per(3, 9) == 3, "x^2 mod 9 has period 9, x^3 mod 9 has period 3 (the owner's two cycles)")
ok = True
for m in range(2, 10):
    n = 3 ** m
    U = [x for x in range(n) if x % 3]
    ok &= {x * x % n for x in U} == {x for x in U if x % 3 == 1}
    ok &= {pow(x, 3, n) for x in U} == {x for x in U if x % 9 in (1, 8)}
check(ok, "unit squares mod 3^m = {1 mod 3} (decided mod 3); unit cubes mod 3^m = {+-1 mod 9} (decided mod 9), m = 2..9")
info("conservation: squaring loses no 3-adic input digit and its image is read mod 3; cubing loses one input digit and its image needs one extra digit (mod 9)")
bad = []
for p in primes_below(200)[1:]:
    for k in range(2, 9):
        D = powset(p, k)
        if len(D) == p:
            continue
        if not is_prime_digraph(p, D):
            bad.append((p, k))
check(not bad, "prime-level power clocks Cay(Z/p, D_k(p)) are prime (indecomposable) digraphs unless complete (3 <= p < 200, k = 2..8)")
bad = []
for p in primes_below(60)[1:]:  # p = 2 is degenerate: two vertices have no non-trivial module
    for _ in range(25):
        D = frozenset(d for d in range(p) if rng.random() < 0.5)
        trivial = len(D - {0}) in (0, p - 1)
        if is_prime_digraph(p, D) == trivial:
            bad.append((p, sorted(D)))
check(not bad, "every circulant Cay(Z/p, D) with p an odd prime is a prime digraph iff D minus {0} is neither empty nor everything (3 <= p < 60, 25 random D each)")
cache = {}
for k in range(2, 7):
    for q in range(2, 301):
        f = factor(q)
        if len(f) == 1:
            cache[(q, k)] = has_clique_module(q, powset(q, k))
mism, tested = [], 0
for k in range(2, 7):
    for n in range(6, 301):
        f = factor(n)
        if len(f) < 2:
            continue
        tested += 1
        pred = not any(cache[(p ** e, k)] for p, e in f.items())
        if is_prime_digraph(n, powset(n, k)) != pred:
            mism.append((n, k))
check(not mism, "FINITE-EXACT: a power clock on n with >= 2 distinct primes is prime iff no prime-power factor has a clique module (%d cases, n <= 300, k = 2..6)" % tested)
odd_rule = all(cache[(q, k)] == ((m := sum(factor(q).values())) and (m - 1) % k == 0 and math.gcd(k, p - 1) == 1)
               for (q, k) in cache for p in factor(q) if p != 2)
check(odd_rule, "odd p: Cl(p^e, D_k(p^e)) has a clique module iff k | e-1 and gcd(k, p-1) = 1 (prime powers <= 300, k = 2..6); for e = 1 it is the whole complete clock")
ok = True
for p_ in (3, 5, 7, 11):
    for k in range(2, 7):
        a = 0
        while k % p_ ** (a + 1) == 0:
            a += 1
        for m in range(1, 8):
            n = p_ ** m
            if n > 2500:
                break
            D = powset(n, k)
            seen = {}
            for d in range(1, n):
                v, w = 0, d
                while w % p_ == 0:
                    v, w = v + 1, w // p_
                key = (v, w % p_ ** min(m - v, 1 + a))
                t = typ(n, D, d)
                ok &= seen.setdefault(key, t) == t
                ok &= (v % k == 0) or t == (False, False)
check(ok, "tower law (p = 3,5,7,11; k = 2..6; p^m <= 2500): the type of d = p^v w depends only on v and w mod p^min(m-v, 1+v_p(k)), and is missing unless k | v")
def is_tower(q, k):
    p_ = min(factor(q))
    e = factor(q)[p_]
    D = powset(q, k)
    return any(all(is_module(q, D, set(range(r, q, p_ ** j))) for r in range(p_ ** j)) for j in range(1, e))


ok = all(is_prime_digraph(n, powset(n, k)) for n, k in [(63, 2), (45, 2), (225, 2), (175, 3), (189, 3), (225, 3)])
ok &= is_tower(9, 2) and is_tower(25, 2) and is_tower(25, 3) and is_tower(27, 3) and not is_tower(9, 3)
check(ok, "squares mod 63, 45, 225 and cubes mod 175, 189, 225 are prime, although each has a tower factor (squares mod 9 and 25, "
          "cubes mod 25 and 27): CRT dissolves towers; cubes mod 9 is not a tower (it is the prime 9-cycle)")
ok = True
for p_, k in [(3, 2), (3, 3), (3, 6), (3, 9), (3, 18), (5, 2), (5, 5), (5, 4), (5, 20), (7, 7), (7, 3)]:
    a = 0
    while k % p_ ** (a + 1) == 0:
        a += 1
    for m in range(1, 7):
        n = p_ ** m
        if n > 20000:
            break
        U = [x for x in range(1, n) if x % p_]
        img = {x: pow(x, k, n) for x in U}
        j = min(j for j in range(m + 1) if len({(x % p_ ** j, img[x]) for x in U}) == len({x % p_ ** j for x in U}))
        pred = m - a if m >= a + 2 else (0 if k % (p_ - 1) == 0 else 1)
        ok &= j == pred
check(ok, "digit conservation (Prop 6): on units mod p^m, x^k depends exactly on x mod p^(m-a), a = v_p(k), when m >= a+2; "
          "for m <= a+1 only on x mod p (on nothing when (p-1) | k) (p = 3, 5, 7, several k, p^m <= 20000)")

# ---------------------------------------------------------------- C
print("=== C. Collatz on the mod-9 clocks ===")
ok = True
dist = Counter()
for n in range(1, 2 * 10 ** 6, 2):
    t, v = syr(n)
    ok &= (3 * n + 1) % 9 in (1, 4, 7)
    ok &= pow(t, 3, 9) == (1 if v % 2 == 0 else 8)
    ok &= ((t % 9) in (1, 4, 7)) == (v % 2 == 0)
    ok &= t % 9 == (-1) ** v * pow(4, (n + v) % 3, 9) % 9
    dist[t % 9] += 1
check(ok, "odd n < 2e6: 3n+1 is a unit square mod 9 (digital root 1,4,7); T(n)^3 = (-1)^v (mod 9); T(n) is a square mod 9 iff v is even; "
          "exactly T(n) = (-1)^v * 4^(n+v) (mod 9): sign in the cube clock {1,8}, rotation n+v mod 3 in the square clock {1,4,7}")
N = sum(dist.values())
frac = {r: dist[r] / N for r in sorted(dist)}
check(all(abs(frac[r] - (1 / 9 if r in (1, 4, 7) else 2 / 9)) < 1e-3 for r in frac), "forward Syracuse images mod 9: 1/9 on each square, 2/9 on each non-square (observed %s)" % {r: round(x, 4) for r, x in frac.items()})
# orbit law: from the second iterate on, n_(i+1) mod 9 = (-1)^v_i 4^(r + v_i) with r = n_i mod 3 = 1 + [v_(i-1) odd]
from fractions import Fraction
law = Counter()
for r, pr in ((1, Fraction(1, 3)), (2, Fraction(2, 3))):
    for j in range(1, 7):  # P(v = j mod 6) = 2^(6-j)/63
        law[(-1) ** j * pow(4, (r + j) % 3, 9) % 9] += pr * Fraction(2 ** (6 - j), 63)
check({x: law[x] * 63 for x in sorted(law)} == {1: 8, 2: 16, 4: 11, 5: 4, 7: 2, 8: 22},
      "orbit law mod 9 (exact, two independent geometric halving counts): (1,2,4,5,7,8) get (8,16,11,4,2,22)/63")
dist2 = Counter()
for n in range(1, 2 * 10 ** 6, 2):
    dist2[syr(syr(n)[0])[0] % 9] += 1
N2 = sum(dist2.values())
check(all(abs(dist2[x] / N2 - float(law[x])) < 1e-3 for x in law),
      "second Syracuse iterates of odd n < 2e6 follow the orbit law to 1e-3 (observed x63: %s): 7 mod 9 is rare, 8 mod 9 common"
      % {x: round(63 * dist2[x] / N2, 3) for x in sorted(dist2)})


def digits_from_valuations(k, starts=range(1, 100001, 2), steps=60):
    L = [2 * 3 ** j for j in range(k)]
    table = {}
    for n0 in starts:
        n, vs = n0, []
        for s in range(steps):
            if n == 1 and s > 0:
                break
            t, v = syr(n)
            vs.append(v)
            if len(vs) >= k:
                key = tuple(vs[-k + j] % L[j] for j in range(k))
                if table.setdefault(key, t % 3 ** k) != t % 3 ** k:
                    return False
            n = t
    return True


check(all(digits_from_valuations(k) for k in range(1, 6)),
      "n_(i+1) mod 3^k is a function of (v_(i-k+1) mod 2, ..., v_i mod 2*3^(k-1)) for k = 1..5 (odd starts < 1e5): the 3-adic clock forgets the start")
ok = True
for m in range(2, 13):
    S2, x = {1}, 2
    while x != 1:
        S2.add(x)
        x = 2 * x % 3 ** m
    ok &= len(S2) == 2 * 3 ** (m - 1) and (3 ** m - 1) in S2
    S3, x = {1}, 3
    while x != 1:
        S3.add(x)
        x = 3 * x % 2 ** (m + 1)
    ok &= (2 ** (m + 1) - 1) not in S3 and S3 == {y for y in range(1, 2 ** (m + 1)) if y % 8 in (1, 3)}
ok &= (4 - 1) in {pow(3, i, 4) for i in range(2)}
check(ok, "duality (m = 2..12): <2> mod 3^m = all units, contains -1 (doubled, complete 3-partite clock); <3> mod 2^(m+1) = {1,3 mod 8}, misses -1 "
          "(bipartite tournament); at m = 1, <3> mod 4 = {1, 3} contains -1")
sq9, cu9 = {x * x % 9 for x in range(1, 9) if x % 3}, {pow(x, 3, 9) for x in range(1, 9) if x % 3}
gens9 = {g for g in range(1, 9) if g % 3 and len({pow(g, i, 9) for i in range(6)}) == 6}
check(gens9 == {g for g in range(1, 9) if g % 3 and g not in sq9 and g not in cu9} and 2 in gens9,
      "units mod 9: generators = units that are neither squares {1,4,7} nor cubes {1,8}; 2 is one (so 2 generates every (Z/3^k)^x)")
ok, rows = True, []
for a in primes_below(400)[1:]:
    a2 = a * a
    prim = len({pow(2, i, a2) for i in range(a * (a - 1))}) == a * (a - 1)
    qs = set(factor(a * (a - 1)))
    avoid = all(pow(2, a * (a - 1) // q, a2) != 1 for q in qs)
    ok &= prim == avoid
    if a < 40:
        rows.append((a, prim, sorted(q for q in qs if pow(2, a * (a - 1) // q, a2) == 1)))
check(ok, "odd primes a < 400: 2 is a primitive root mod a^2 iff 2 lies in no q-th power clock mod a^2, q | a(a-1) prime")
info("(a, 2 primitive mod a^2, power clocks containing 2): %s" % rows)
info("this is exactly the condition of THM-4523's Theorem R_a (an iff) for n/2, an+1; at a = 3 it is 'squares mod 9' and 'cubes mod 9'")
ok, qrows = True, []
for a in (3, 5, 7, 11, 13, 17, 1093, 3511):
    a2 = a * a
    w_inv = pow(pow(2, a, a2), -1, a2)          # omega(2)^(-1), omega(2) = 2^a mod a^2 (Teichmuller lift)
    qa = (pow(2, a - 1, a2) - 1) // a % a        # Fermat quotient q_a(2) mod a
    qrows.append((a, qa))
    for n in range(1, 200001, 2):
        m = a * n + 1
        v = v2(m)
        ok &= (m >> v) % a2 == pow(w_inv, v, a2) * pow(1 + a, (n + qa * v) % a, a2) % a2
check(ok, "an+1 maps (a = 3,5,7,11,13,17,1093,3511; odd n < 2e5): (an+1)/2^v = omega(2)^(-v) (1+a)^(n + q_a(2) v) (mod a^2); "
          "each halving turns the mu_(a-1) clock by omega(2)^(-1) and the 1+aZ clock by the Fermat quotient q_a(2), frozen iff a is Wieferich")
info("(a, q_a(2) mod a): %s" % qrows)

# ---------------------------------------------------------------- D
print("=== D. triplets: Pythagorean triples as transitive triangles of the square clock ===")


def ppts(C):
    out, m = [], 2
    while m * m < C:
        for n in range(1, m):
            if (m - n) % 2 and math.gcd(m, n) == 1:
                a, b, c = m * m - n * n, 2 * m * n, m * m + n * n
                if c <= C:
                    out.append((a, b, c))
        m += 1
    return out


T = ppts(200000)
ok, five = True, Counter()
for a, b, c in T:
    ok &= b % 4 == 0 and a % 2 == 1 and (a * b) % 3 == 0 and c % 3 != 0 and (a * b * c) % 5 == 0
    ok &= all(p % 4 == 1 for p in factor(c))
    five["odd leg" if a % 5 == 0 else ("even leg" if b % 5 == 0 else "hypotenuse")] += 1
check(ok, "all %d PPTs with c <= 2e5: 4 | even leg, 3 | a leg and never c, 5 | abc, every prime of c is 1 mod 4" % len(T))
info("which member carries the 5: %s (the mod-5 clock is symmetric, so the hypotenuse may take the loop)" % dict(five))
free = []
for p in primes_below(400):
    Q = {t * t % p for t in range(1, p)}
    if not any((x + y) % p in Q for x in Q for y in Q):
        free.append(p)
check(free == [2, 3, 5], "the square clock mod p has a unit transitive triangle (x + y = z, all nonzero squares) iff p >= 7 (p < 400)")
hyp = {p for a, b, c in T for p in factor(c)}
ok = all((p in hyp) == (powset(p, 2) == frozenset((-d) % p for d in powset(p, 2))) == (p % 4 == 1) for p in primes_below(400)[1:])
ok &= 2 not in hyp and powset(2, 2) == frozenset({0, 1})
check(ok, "odd p < 400: p divides some hypotenuse (c <= 2e5) iff the square clock mod p is symmetric iff p = 1 mod 4; at p = 2 the clock is "
          "symmetric but every hypotenuse is odd")
ok = True
for p_ in primes_below(200)[1:]:
    Nxy = Counter(((m * m - n * n) % p_, 2 * m * n % p_) for m in range(p_) for n in range(p_))
    ok &= all(Nxy[(x, y)] == Nxy[(y, x)] for x in range(p_) for y in range(p_)) == (p_ % 8 != 5)
check(ok, "leg exchange: #{(m, n) mod p : (m^2 - n^2, 2mn) = (x, y)} is symmetric in (x, y) iff p != 5 (mod 8), odd p < 200 (i is a square in F_p[i] iff p != 5 mod 8)")
jac = lambda x, c: math.prod(legendre(x, q) ** e for q, e in factor(c).items())
ok = all(jac(b, c) == 1 and jac(a, c) == jac(2, c) == jac(2, a) for a, b, c in T)
incid = 0
for a, b, c in T:
    for q_ in set(factor(a)) | set(factor(b)) | set(factor(c)):
        if q_ % 4 != 1:
            continue
        incid += 1
        if a % q_ == 0:
            ok &= legendre(b, q_) == legendre(2, q_)
        if b % q_ == 0:
            ok &= legendre(a, q_) == 1
        if c % q_ == 0:
            ok &= legendre(a, q_) == legendre(2, q_) and legendre(b, q_) == 1
check(ok, "PPTs c <= 2e5: Jacobi (b/c) = +1 and (a/c) = (2/c) = (2/a); for every prime q = 1 mod 4 dividing abc (%d incidences): "
          "q | a => (b/q) = (2/q), q | b => (a/q) = +1, q | c => (a/q) = (2/q), (b/q) = +1: the factor 2 of the even leg is visible "
          "exactly where (2/q) = -1, i.e. q = 5 mod 8" % incid)
check(not any(((x + y) % 9) in (1, 4, 7) for x in (1, 4, 7) for y in (1, 4, 7)), "square clock mod 9 (C3 blown up) has no unit transitive triangle: every PPT has 3 | a leg")
check(not any(((x + y) % 9) in (1, 8) for x in (1, 8) for y in (1, 8)), "cube clock mod 9 (9-cycle) has no unit triangle: x^3 + y^3 = z^3 forces 3 | xyz (Fermat case I for p = 3)")

# ---------------------------------------------------------------- E
print("=== E. F = S + U, p^2qr and the 3-4-3 split ===")


def fsuq(N):
    F = [d for d in divisors(N) if 1 < d < N]
    U = [d for d in F if is_prime(d)]
    S = [d for d in F if all(e == 1 for e in factor(d).values())]
    Q = [d for d in F if d not in S]
    return F, S, U, [d for d in S if d not in U], Q


check(all((len(F) == len(S) + len(U)) == (len(Q) == len(U)) for N in range(2, 20000) for F, S, U, M, Q in [fsuq(N)]),
      "F = S + U iff |Q| = |U| (square-containing proper divisors = prime divisors), N < 2e4")
p, q, r = 2, 3, 5
N = p * p * q * r
F, S, U, M, Q = fsuq(N)
lab = {**{d: "U" for d in U}, **{d: "M" for d in M}, **{d: "Q" for d in Q}}
pairs = sorted({tuple(sorted((d, N // d))) for d in F})
check((len(U), len(M), len(Q)) == (3, 4, 3), "N = p^2qr = 60: classes U = %s, M = %s, Q = %s (3-4-3)" % (U, M, Q))
mg = Counter("".join(sorted(lab[a] + lab[b])) for a, b in pairs)
check(mg == Counter({"QU": 2, "MU": 1, "MQ": 1, "MM": 1}), "the 5 complementary pairs %s give the class multigraph U=Q doubled, M looped, U-M, M-Q" % pairs)
lay = {d: factor(d).get(p, 0) for d in F}
lg = Counter(tuple(sorted((lay[a], lay[b]))) for a, b in pairs)
check(lg == Counter({(0, 2): 3, (1, 1): 2}), "the p-valuation layering L0 = %s, L1 = %s, L2 = %s is also 3-4-3 and complementation swaps L0 <-> L2, fixes L1"
      % ([d for d in F if lay[d] == 0], [d for d in F if lay[d] == 1], [d for d in F if lay[d] == 2]))
check(sorted(math.lcm(u, p * p) for u in U) == sorted(Q), "Q = lcm(U, p^2): doubling the exponent of the repeated prime carries U onto Q (the p^3 move)")
sqc = lambda d: tuple(factor(d).get(x, 0) % 2 for x in (p, q, r))
check(all(tuple((x + y) % 2 for x, y in zip(sqc(a), sqc(b))) == sqc(N) for a, b in pairs),
      "in squareclass coordinates F_2^3 complementation is translation by [N] = [qr]; its orbits are the 3 Fano lines through [qr] and {0, [qr]}")
F3, _, U3, M3, Q3 = fsuq(27)
check((U3, M3, Q3) == ([3], [], [9]) and 27 // 3 == 9, "N = p^3: U = {p}, Q = {p^2}, and complementation itself swaps U and Q")
ok, rows = True, []
from functools import reduce
for rr in range(1, 7):
    ps = primes_below(30)[:rr]
    Nn = ps[0] ** 2 * reduce(lambda x, y: x * y, ps[1:], 1)
    F_, S_, U_, M_, Q_ = fsuq(Nn)
    L = Counter(factor(d).get(ps[0], 0) for d in F_)
    ok &= [L[i] for i in range(3)] == [2 ** (rr - 1) - 1, 2 ** (rr - 1), 2 ** (rr - 1) - 1]
    ok &= (len(U_), len(M_), len(Q_)) == (rr, 2 ** rr - 1 - rr, 2 ** (rr - 1) - 1)
    rows.append((rr, (len(U_), len(M_), len(Q_)), [L[i] for i in range(3)]))
check(ok, "profile (2,1^(r-1)): classes (r, 2^r-1-r, 2^(r-1)-1), layers (2^(r-1)-1, 2^(r-1), 2^(r-1)-1); they agree iff r = 3: %s" % rows)

# ---------------------------------------------------------------- F
print("=== F. Fermat case I as a triangle census of the p-th power clock mod p^2 (Teichmuller clock) ===")
check(sorted({pow(x, 3, 9) for x in range(9)}) == [0, 1, 8], "p = 3: the p-th power clock mod p^2 is the owner's cube clock mod 9")
stats = {1: Counter(), 2: Counter()}
orbit_ok, isos, direct_ok, nroot2, tot2 = True, [], True, 0, 0
for p in primes_below(P_FERMAT)[2:]:
    g, p2 = primitive_root(p), p * p
    W = pow(g, p, p2)
    om = [0] * p
    x, w = 1, 1
    for _ in range(p - 1):
        om[x] = w
        x, w = x * g % p, w * W % p2
    roots = [t for t in range(1, p - 1) if (om[t + 1] - 1 - om[t]) % p2 == 0]
    if p < 400:
        Tset = set(om[1:])
        tri = any((a + b) % p2 in Tset for a in om[1:] for b in om[1:])
        direct = [t for t in range(1, p - 1) if (pow(1 + t, p, p2) - 1 - pow(t, p, p2)) % p2 == 0]
        direct_ok &= (direct == roots) and (tri == bool(roots))
    cls = 1 if p % 3 == 1 else 2
    R = set(roots)
    orbit_ok &= R == {pow(t, -1, p) for t in R} and R == {(-1 - t) % p for t in R}
    if cls == 1:
        w3 = pow(g, (p - 1) // 3, p)
        orbit_ok &= w3 in R and w3 * w3 % p in R
    orbit_ok &= (1 in R) == (pow(2, p - 1, p2) == 1)
    stats[cls][len(roots)] += 1
    if cls == 2:
        tot2 += 1
        nroot2 += not roots
    if 1 in roots:
        isos.append(p)
    base = (2 if cls == 1 else 0) + (3 if 1 in roots else 0)
    orbit_ok &= len(roots) >= base and (len(roots) - base) % 6 == 0
check(direct_ok, "p < 400: unit triangles x^p + y^p = z^p (mod p^2) exist iff ((1+t)^p - 1 - t^p)/p has a root t != 0, -1 mod p (direct pow cross-check)")
check(orbit_ok, "orbit law for 5 <= p < %d: the root set is invariant under t -> 1/t and t -> -1-t, contains the Eisenstein pair "
                "(1 + w = -w^2) iff p = 1 mod 3, contains t = 1 iff p is Wieferich, and #roots = 2[p = 1 mod 3] + 3[Wieferich] + 6j" % P_FERMAT)
wief = [p for p in primes_below(P_FERMAT)[1:] if pow(2, p - 1, p * p) == 1]
check(isos == wief == [1093, 3511], "the isosceles (doubling) triangle x + x = 2x lies in the clock exactly at the Wieferich primes %s" % isos)
info("p = 1 mod 3 root counts: %s" % sorted(stats[1].items()))
info("p = 2 mod 3 root counts: %s" % sorted(stats[2].items()))
info("p = 2 mod 3 with NO unit triangle mod p^2 (case I of FLT then follows from this congruence alone): %d / %d = %.4f; Poisson(1/6) heuristic exp(-1/6) = %.4f"
     % (nroot2, tot2, nroot2 / tot2, math.exp(-1 / 6)))
check(nroot2 / tot2 > 0.8, "a large majority of p = 2 mod 3 have a unit-triangle-free Teichmuller clock (first ones: %s)"
      % [p for p in primes_below(200)[2:] if p % 3 == 2 and not any((pow(1 + t, p, p * p) - 1 - pow(t, p, p * p)) % (p * p) == 0 for t in range(1, p - 1))])
ok = all((pow(2, p - 1, p * p) == 1) == (2 in {pow(x, p, p * p) for x in range(1, p)}) and
         (pow(3, p - 1, p * p) == 1) == (3 in {pow(x, p, p * p) for x in range(1, p)}) for p in primes_below(1200)[2:])
check(ok, "2 (resp. 3) lies in the p-th power clock mod p^2 iff p is a base-2 (resp. base-3) Wieferich prime (5 <= p < 1200; 1093 and 11 occur)")

# ---------------------------------------------------------------- G
print("=== G. local Waring numbers: loops make CRT products strong products, so distances take the max ===")
for k, top_expect, where in [(2, 4, 8), (3, 4, 9), (4, 15, 16)]:
    badmax, best = 0, (0, None)
    for n in range(2, 600 if k < 4 else 400):
        g = diam(n, powset(n, k))
        gp = max(diam(p ** e, powset(p ** e, k)) for p, e in factor(n).items())
        badmax += g != gp
        if g > best[0]:
            best = (g, n)
    check(badmax == 0 and best == (top_expect, where), "k = %d: local Waring number of n = max over its prime powers; maximum %d first at n = %d" % (k, best[0], best[1]))
D9 = powset(9, 3)
check(sorted(set(range(9)) - {(a + b + c) % 9 for a in D9 for b in D9 for c in D9}) == [4, 5], "residues mod 9 that are not sums of three cubes: 4, 5 (the antipodes of the 9-cycle at distance 4)")

# ---------------------------------------------------------------- H
print("=== H. doubling clocks Cay(Z/p, <2>) ===")
cnt, paley = Counter(), []
for p in primes_below(P_DOUBLING)[1:]:
    o = mult_order(2, p)
    if o == p - 1:
        cnt["complete"] += 1
    elif o % 2 == 0:
        cnt["symmetric"] += 1
    elif o == (p - 1) // 2:
        cnt["paley"] += 1
        paley.append(p)
    else:
        cnt["oriented+missing"] += 1
tot = sum(cnt.values())
ok = True
for p in primes_below(3000)[1:]:
    S2, x = {1}, 2
    while x != 1:
        S2.add(x)
        x = 2 * x % p
    ts = Counter(typ(p, S2, d) for d in range(1, p))
    if ts[(True, True)] == p - 1:
        kind = "complete"
    elif ts[(True, False)] + ts[(False, True)] == p - 1:
        kind = "paley"
    elif ts[(True, False)] == 0:
        kind = "symmetric"
    else:
        kind = "oriented+missing"
    o = mult_order(2, p)
    pred = "complete" if o == p - 1 else "symmetric" if o % 2 == 0 else "paley" if o == (p - 1) // 2 else "oriented+missing"
    ok &= kind == pred and (kind != "paley" or S2 == {t * t % p for t in range(1, p)})
check(ok, "doubling clocks read off their actual pair types (odd p < 3000) agree with the ord_p(2) rule; the tournaments are the Paley tournaments (<2> = squares)")
check(all(p % 8 == 7 for p in paley) and paley[:8] == [7, 23, 47, 71, 79, 103, 167, 191],
      "doubling clock mod p is a (Paley) tournament iff p = 7 mod 8 and ord_p(2) = (p-1)/2; first: %s" % paley[:12])
info("odd primes < %d: complete %.5f, Paley %.5f, symmetric %.5f, oriented-with-missing %.5f"
     % (P_DOUBLING, cnt["complete"] / tot, cnt["paley"] / tot, cnt["symmetric"] / tot, cnt["oriented+missing"] / tot))
A = 1.0
for ell in primes_below(100000):
    A *= 1 - 1 / (ell * (ell - 1))
info("Artin's constant A = %.6f (product over l < 1e5); A/2 = %.6f" % (A, A / 2))
check(abs(cnt["complete"] / tot - A) < 0.002 and abs(cnt["paley"] / tot - A / 2) < 0.002,
      "observed shares of complete and Paley doubling clocks agree with A and A/2 to 0.002 (heuristic comparison, not a theorem)")
sym_share = (cnt["complete"] + cnt["symmetric"]) / tot
check(abs(sym_share - 17 / 24) < 0.002, "share of symmetric doubling clocks (ord_p(2) even, i.e. -1 in <2>) = %.5f vs 17/24 = %.5f (unconditional density theorem for base 2)" % (sym_share, 17 / 24))
check(3511 in paley and 1093 not in paley, "the Wieferich prime 3511 is a Paley doubling prime (ord = 1755); 1093 is not (ord = 364): NUMEROLOGY, recorded as data")

print()
if FAILS:
    print("FAILED CHECKS: %d" % len(FAILS))
    for f in FAILS:
        print("  - " + f)
    sys.exit(1)
print("ALL CHECKS PASSED")
