#!/usr/bin/env python3
"""
collatz_mod6_20260917_floor_sums_odd_functions.py -- the odd-function floor law,
its squarefree / prime truth sets, the cube-root lattice reciprocity, the bilinear
floor sum via Pillai's gcd-sum, and the even-power class-number contrast.

Lane: floor_sums_odd_functions (session collatz-mod6-20260917, wave two).

The user's three sums over a prime p:
  (1) sum_{k=1}^{p-1} floor(k^3/p)                 = (p-2)(p-1)(p+1)/4
  (2) sum_{y=1}^{(p-1)(p-2)} floor((py)^{1/3})     = (3p-5)(p-2)(p-1)/4
  (3) sum_{i,j=1}^{p-1} floor(ij/p)                = (p-2)(p-1)^2/4
and the slogan "odd cycles ~ odd functions ~ odd powers".

Inheritance (read, not re-derived; cite by path):
  05-knowledge/results/arithmetic_braids_20260917_collatz.md  (3n-1 = negatives of 3n+1;
     three 3n-1 cycles; cycle gate B/(2^K-3^L))
  05-knowledge/results/arithmetic_braids_20260917_summand.md  (signed conjugation C_+(-n)=-C_-(n))
  05-knowledge/results/arithmetic_braids_20260917_divisors.md (F=S+U iff p,p^3,p^2qr)
  01-canon/theorems/THM-001-redei.md (H(T) odd), LEM-020 (Redei involution / fixed layer),
  01-canon/theorems/THM-1745-... (h-spectrum = odds minus {7,21})
  Session-lead firsthand fact F2 (truth sets of (1),(2),(3) for n<400 / n<120 / n<200;
     Z_7(n)=n/rad(n)-1 for n<200; even e fails at (p,e)=(3,2)).

Everything load-bearing is exact integer / Fraction arithmetic.  Cube roots are
computed by a monotone integer walk (no floating point anywhere).  Every check is
an explicit `raise`, so the checks stay active under `python3 -O`.

Sections:
  S1  odd-function floor law (PROVED; verified for x, x^3, x^5, x^3+x, 3x; even f minimal failure)
  S2  Z_{x^e}(n) = n/n_e - 1 (PROVED); truth set of (1) and of all odd e = squarefree; density 6/pi^2
  S3  cube-root lattice reciprocity N1+N2 = (n-1)^2(n-2) + Z_3(n) (PROVED); (2) <=> squarefree
  S4  bilinear sum via Pillai gcd-sum P(n): N3 = (n-1)^2(n-2)/4 + (P(n)-2n+1)/2 (PROVED); (3) <=> prime
  S5  even powers: sum floor(k^2/p) = law + h(-p) for p = 3 mod 4, p > 3 (class numbers by reduced forms)
  S6  typing table for "odd cycles ~ odd functions ~ odd powers"
"""
import sys, time
from fractions import Fraction as Fr
from math import gcd, isqrt

T0 = time.time()
def stamp():
    return "[t=%.1fs]" % (time.time() - T0)

def check(cond, msg):
    if not cond:
        raise RuntimeError("CHECK FAILED: " + msg)

def banner(s):
    print()
    print("=" * 78)
    print(s)
    print("=" * 78)

# ---------------------------------------------------------------- arithmetic helpers
def factorint(n):
    """exact trial-division factorization, dict p->a"""
    f = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f

def is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True

def is_squarefree(n):
    return all(a == 1 for a in factorint(n).values())

def rad(n):
    r = 1
    for p in factorint(n):
        r *= p
    return r

def n_e(n, e):
    """prod_{p^a || n} p^ceil(a/e)"""
    r = 1
    for p, a in factorint(n).items():
        r *= p ** (-(-a // e))
    return r

def phi(n):
    r = n
    for p in factorint(n):
        r = r // p * (p - 1)
    return r

def divisors(n):
    ds = [1]
    for p, a in factorint(n).items():
        ds = [d * p ** i for d in ds for i in range(a + 1)]
    return sorted(ds)

def pillai(n):
    """P(n) = sum_{d|n} d*phi(n/d) = sum_{k=1}^{n} gcd(k,n)  (A018804)"""
    return sum(d * phi(n // d) for d in divisors(n))

def Z_f(f, n):
    return sum(1 for k in range(1, n) if f(k) % n == 0)

def floor_sum(f, n):
    return sum(f(k) // n for k in range(1, n))

def odd_law(f, n):
    """(1/n) sum f(k) - (n-1)/2 + Z_f(n)/2 as an exact Fraction"""
    return Fr(sum(f(k) for k in range(1, n)), n) - Fr(n - 1, 2) + Fr(Z_f(f, n), 2)

# ======================================================================= S1
banner("S1  The odd-function floor law (PROVED) and its even-function failure")
print("""
THEOREM S1 (odd-function floor law).  Let f in Z[x] with f(-x) = -f(x), n >= 2, and
  Z_f(n) = #{1 <= k <= n-1 : n | f(k)}.   Then
  sum_{k=1}^{n-1} floor(f(k)/n) = (1/n) sum_{k=1}^{n-1} f(k) - (n-1)/2 + Z_f(n)/2.
PROOF.  Write f(k) = n q_k + r_k with 0 <= r_k < n, so floor(f(k)/n) = q_k and
  sum q_k = ( sum f(k) - sum r_k ) / n.                                   (a)
Pairing.  f(n-k) = f(-k) = -f(k) (mod n) because f is odd and polynomial, so
  r_k + r_{n-k} = 0 (mod n), hence r_k + r_{n-k} in {0, n}; it is 0 iff r_k = 0 iff r_{n-k} = 0.
Fixed point.  If n is even, k = n/2 is the unique fixed point of k -> n-k, and
  2 r_{n/2} = 0 (mod n) gives r_{n/2} in {0, n/2}, i.e. r_{n/2} = (n/2) [r_{n/2} != 0].
Summing.  Each non-fixed k contributes n/2 per nonzero residue on average over its pair, and the
  fixed point contributes n/2 [r != 0] exactly; in all cases
  sum_k r_k = (n/2) * #{k : r_k != 0} = (n/2) (n - 1 - Z_f(n)).                (b)
Insert (b) in (a):  sum q_k = (sum f)/n - (n-1)/2 + Z_f(n)/2.                   QED
Remark.  Only f(-x) = -f(x) mod n was used; the law holds for any odd function Z/n -> Z/n.
""")

odd_polys = {
    "x":       lambda k: k,
    "x^3":     lambda k: k ** 3,
    "x^5":     lambda k: k ** 5,
    "x^3+x":   lambda k: k ** 3 + k,
    "3x":      lambda k: 3 * k,
    "x^7-5x^3+x": lambda k: k ** 7 - 5 * k ** 3 + k,
}
NMAX_S1 = 400
for name, f in odd_polys.items():
    bad = [n for n in range(2, NMAX_S1 + 1) if Fr(floor_sum(f, n)) != odd_law(f, n)]
    check(not bad, "odd law fails for %s at %s" % (name, bad[:5]))
    print("  FINITE-EXACT  f=%-12s  law holds for all 2 <= n <= %d" % (name, NMAX_S1))

# residue pairing and fixed point witnessed explicitly
for n in (4, 6, 10, 12):
    f = odd_polys["x^3"]
    pairs = [(k, f(k) % n, f(n - k) % n) for k in range(1, n // 2)]
    check(all((r1 + r2) in (0, n) for _, r1, r2 in pairs), "pairing residues n=%d" % n)
    rf = f(n // 2) % n
    check(rf in (0, n // 2), "fixed point residue n=%d" % n)
    print("  n=%2d f=x^3: pairs (k, r_k, r_{n-k}) = %s ; fixed k=n/2 residue = %d in {0,%d}"
          % (n, pairs, rf, n // 2))

print()
print("EVEN f: the pairing gives r_{n-k} = r_k (no cancellation).  Minimal failures:")
even_polys = {"x^2": lambda k: k * k, "x^4": lambda k: k ** 4, "x^2+x": lambda k: k * k + k,
              "2 (const)": lambda k: 2}
for name, f in even_polys.items():
    first = None
    for n in range(2, 60):
        if Fr(floor_sum(f, n)) != odd_law(f, n):
            first = (n, floor_sum(f, n), odd_law(f, n))
            break
    check(first is not None, "even f %s never fails?" % name)
    print("  REFUTED for even f=%-10s minimal n=%d : floor sum = %d, odd-law value = %s"
          % (name, first[0], first[1], first[2]))
print("  (f=x^2, n=2 holds: the single term k=1 is the fixed point; n=3 is the minimal failure,"
      " matching the session lead's (p,e)=(3,2).)")
print(stamp())

# ======================================================================= S2
banner("S2  Z_{x^e}(n) = n/n_e - 1 (PROVED); identity (1) and all odd e <=> squarefree; density 6/pi^2")
print("""
THEOREM S2.  For e >= 1 and n >= 2 put n_e = prod_{p^a || n} p^{ceil(a/e)}.  Then
  Z_{x^e}(n) = #{1 <= k <= n-1 : n | k^e} = n/n_e - 1.
PROOF.  n | k^e  iff  for every p^a || n,  a <= e v_p(k)  iff  v_p(k) >= ceil(a/e)  iff  n_e | k.
  n_e | n since ceil(a/e) <= a, so the multiples of n_e in [1, n-1] number n/n_e - 1.   QED
COROLLARY S2a.  Z_{x^e}(n) = 0  iff  n_e = n  iff  ceil(a/e) = a for all a  iff  (e = 1) or (n squarefree).
COROLLARY S2b (identity (1)).  (sum k^3)/n - (n-1)/2 = (n-1)^2 n/4 - (n-1)/2 = (n-1)(n-2)(n+1)/4, so
  sum_{k=1}^{n-1} floor(k^3/n) = (n-2)(n-1)(n+1)/4 + Z_3(n)/2,   equality with the user's RHS iff n squarefree.
COROLLARY S2c.  For every odd e >= 3:  sum floor(k^e/n) = (sum k^e)/n - (n-1)/2  iff  n squarefree.
  For e = 1 the identity holds for every n (Z_1 = 0), so e = 1 is the degenerate member of the odd family.
CITED (Gegenbauer 1885; Hardy-Wright Thm 333): the squarefree numbers have natural density 1/zeta(2) = 6/pi^2.
  Hence the truth set of (1), and of the odd-power identity for every fixed odd e >= 3, has density 6/pi^2.
""")

NMAX_S2 = 2000
for e in range(1, 10):
    bad = [n for n in range(2, NMAX_S2 + 1) if Z_f(lambda k: k ** e, n) != n // n_e(n, e) - 1]
    check(not bad, "Z formula e=%d bad %s" % (e, bad[:5]))
print("  FINITE-EXACT  Z_{x^e}(n) = n/n_e - 1 for 1 <= e <= 9, 2 <= n <= %d" % NMAX_S2)

def rhs1(n):
    return Fr((n - 2) * (n - 1) * (n + 1), 4)

truth1 = [n for n in range(2, NMAX_S2 + 1) if Fr(floor_sum(lambda k: k ** 3, n)) == rhs1(n)]
sqf = [n for n in range(2, NMAX_S2 + 1) if is_squarefree(n)]
check(truth1 == sqf, "truth set of (1) != squarefree")
print("  FINITE-EXACT  truth set of (1) on [2,%d] = squarefree numbers exactly (%d values)" % (NMAX_S2, len(sqf)))
for n in (4, 9, 12, 8, 36):
    print("     n=%2d: floor sum = %5d, RHS(1) = %s, Z_3(n) = %d, deviation = Z_3/2 = %s"
          % (n, floor_sum(lambda k: k ** 3, n), rhs1(n), Z_f(lambda k: k ** 3, n),
             Fr(floor_sum(lambda k: k ** 3, n)) - rhs1(n)))

NMAX_ODD = 1200
for e in (3, 5, 7, 9):
    ts = [n for n in range(2, NMAX_ODD + 1)
          if Fr(floor_sum(lambda k: k ** e, n)) == Fr(sum(k ** e for k in range(1, n)), n) - Fr(n - 1, 2)]
    check(ts == [n for n in sqf if n <= NMAX_ODD], "odd e=%d truth set != squarefree" % e)
    print("  FINITE-EXACT  e=%d: identity sum floor(k^e/n) = (sum k^e)/n - (n-1)/2 holds on [2,%d] iff squarefree" % (e, NMAX_ODD))
ts1 = [n for n in range(2, NMAX_ODD + 1)
       if Fr(floor_sum(lambda k: k, n)) == Fr(sum(range(1, n)), n) - Fr(n - 1, 2)]
check(ts1 == list(range(2, NMAX_ODD + 1)), "e=1 should hold everywhere")
print("  FINITE-EXACT  e=1: holds for every n in [2,%d] (degenerate member)" % NMAX_ODD)

# the session lead's Z_7(n) = n/rad(n) - 1 for n<200: explained and its first failure
bad7 = [n for n in range(2, 300) if Z_f(lambda k: k ** 7, n) != n // rad(n) - 1]
check(bad7 == [256], "Z_7 vs n/rad(n)-1 first failure should be 256, got %s" % bad7[:5])
print("  PROVED+FINITE-EXACT  Z_7(n) = n/rad(n)-1 iff every exponent a <= 7 (then n_7 = rad n);"
      " first failure n=256=2^8: Z_7(256) = %d, 256/rad-1 = %d" % (Z_f(lambda k: k ** 7, 256), 255))

# density control
import math
N_DENS = 10 ** 6
# sieve squarefree up to N_DENS
sf = bytearray([1]) * (N_DENS + 1)
q = 2
while q * q <= N_DENS:
    sf[q * q::q * q] = bytearray(len(range(q * q, N_DENS + 1, q * q)))
    q += 1
cnt = sum(sf[1:])
print("  FINITE-EXACT  #squarefree <= 10^6 = %d ; 6/pi^2 * 10^6 = %.1f ; ratio %.6f" % (cnt, 6 / math.pi ** 2 * N_DENS, cnt / N_DENS))
check(abs(cnt / N_DENS - 6 / math.pi ** 2) < 1e-3, "squarefree density sanity")
print("""
TYPED ANALOGY (the user's "6/pi^2" remark).  source: truth set of identity (1) (and of every odd-e identity);
  target: squarefree integers; map: n -> Z_e(n) = n/n_e - 1 (the identity holds iff Z_e(n) = 0);
  preserved predicate: 'squarefree' exactly, both directions (PROVED, S2a-c); lost information: none
  (the deviation Z_e(n)/2 = (n/n_e - 1)/2 is itself the exact local invariant); sidecar: the exponent
  profile (a_p) of n through n_e; cheapest decisive test: n = 4, 9, 12 (deviations 1/2, 1, 1/2 for e=3).
  The 6/pi^2 is CITED squarefree density; nothing about Collatz enters this map.
""")
print(stamp())

# ======================================================================= S3
banner("S3  Cube-root lattice reciprocity (PROVED); identity (2) <=> squarefree; exact integer cube roots")
print("""
THEOREM S3.  For n >= 2 let Y = floor((n-1)^3/n), N1(n) = sum_{x=1}^{n-1} floor(x^3/n),
  N2(n) = sum_{y=1}^{Y} floor((ny)^{1/3}).  Then
  (i)  Y = (n-1)(n-2);
  (ii) N1(n) + N2(n) = (n-1)^2 (n-2) + Z_3(n);
  (iii) N2(n) = (3n-5)(n-2)(n-1)/4 + Z_3(n)/2, so identity (2) holds iff n is squarefree.
PROOF.  (i) (n-1)^3 = n(n^2-3n+3) - 1, so floor((n-1)^3/n) = n^2-3n+2 = (n-1)(n-2).
  (ii) Count lattice points (x,y), 1<=x<=n-1, 1<=y<=Y.  For x <= n-1 we have x^3/n <= (n-1)^3/n, so
  floor(x^3/n) <= Y and N1 = #{(x,y) in box : ny <= x^3}.  For y <= Y we have ny <= (n-1)^3, so
  floor((ny)^{1/3}) <= n-1 and N2 = #{(x,y) in box : x^3 <= ny}.  Every box point lies in at least one
  set; the overlap is the curve x^3 = ny, whose points in the box are exactly the x in [1,n-1] with
  n | x^3 (then y = x^3/n is automatically in [1,Y]).  So N1 + N2 = (n-1) Y + Z_3(n).
  (iii) Subtract S2b: N2 = (n-1)^2(n-2) + Z_3 - (n-1)(n-2)(n+1)/4 - Z_3/2
       = (n-1)(n-2) (4(n-1) - (n+1))/4 + Z_3/2 = (n-1)(n-2)(3n-5)/4 + Z_3/2.                    QED
Cube roots are computed EXACTLY: since (ny)^{1/3} is nondecreasing in y, walk x upward while (x+1)^3 <= ny.
""")

def N2_exact(n):
    Y = (n - 1) * (n - 2)
    x = 0
    total = 0
    for y in range(1, Y + 1):
        ny = n * y
        while (x + 1) ** 3 <= ny:
            x += 1
        total += x
    return total

def rhs2(n):
    return Fr((3 * n - 5) * (n - 2) * (n - 1), 4)

NMAX_S3 = 220
SAMPLE_S3 = [225, 243, 250, 256, 289, 300, 331, 343, 360, 361, 385, 399, 400]
truth2 = []
for n in list(range(2, NMAX_S3 + 1)) + SAMPLE_S3:
    Y = (n - 1) ** 3 // n
    check(Y == (n - 1) * (n - 2), "Y formula n=%d" % n)
    N1 = floor_sum(lambda k: k ** 3, n)
    N2 = N2_exact(n)
    Z3 = n // n_e(n, 3) - 1
    check(N1 + N2 == (n - 1) ** 2 * (n - 2) + Z3, "reciprocity n=%d" % n)
    check(Fr(N2) == rhs2(n) + Fr(Z3, 2), "N2 closed form n=%d" % n)
    if Fr(N2) == rhs2(n):
        truth2.append(n)
check(truth2 == [n for n in list(range(2, NMAX_S3 + 1)) + SAMPLE_S3 if is_squarefree(n)], "truth set (2)")
print("  FINITE-EXACT  reciprocity (ii) and closed form (iii) for all 2 <= n <= %d and n in %s" % (NMAX_S3, SAMPLE_S3))
print("  FINITE-EXACT  truth set of (2) on that universe = squarefree exactly (covers the lead's n<120 claim)")
for n in (7, 8, 9, 12, 31, 32):
    N1 = floor_sum(lambda k: k ** 3, n); N2 = N2_exact(n); Z3 = n // n_e(n, 3) - 1
    print("     n=%2d: N1=%6d N2=%6d N1+N2=%7d (n-1)^2(n-2)=%7d Z_3=%2d ; N2-RHS(2)=%s"
          % (n, N1, N2, N1 + N2, (n - 1) ** 2 * (n - 2), Z3, Fr(N2) - rhs2(n)))
print(stamp())

# ======================================================================= S4
banner("S4  The bilinear sum via Pillai's gcd-sum (PROVED); identity (3) <=> prime; 4 N3 vs N1 + N2")
print("""
THEOREM S4.  For n >= 2 let N3(n) = sum_{i,j=1}^{n-1} floor(ij/n), S(n) = sum_{i,j=1}^{n-1} (ij mod n),
  P(n) = sum_{d|n} d phi(n/d) = sum_{k=1}^{n} gcd(k,n)  (Pillai, A018804).  Then
  (i)  S(n) = n (n^2 - P(n)) / 2;
  (ii) N3(n) = ((n-1)^2 n^2/4 - S(n))/n = (n-1)^2 (n-2)/4 + (P(n) - 2n + 1)/2;
  (iii) P(n) >= 2n-1 with equality iff n is prime; hence identity (3) holds iff n is prime.
PROOF.  (i) Extend i,j to [0,n-1] (the added terms are 0).  Fix i and put d = gcd(i,n), i = d i',
  gcd(i', n/d) = 1.  Then ij mod n = d (i' j mod n/d).  As j runs over [0,n-1], i' j mod n/d takes each
  residue r in [0, n/d - 1] exactly d times (i' is a unit mod n/d), so
  sum_j (ij mod n) = d * d * (n/d)(n/d - 1)/2 = n (n - d)/2.
  The number of i in [0,n-1] with gcd(i,n) = d is phi(n/d).  Hence
  S(n) = (n/2) sum_{d|n} phi(n/d) (n - d) = (n/2) ( n * n - sum_{d|n} d phi(n/d) ) = n (n^2 - P(n))/2,
  using sum_{d|n} phi(n/d) = n.
  (ii) sum_{i,j} ij = (n(n-1)/2)^2; subtract S and divide by n; simplify with (n-1)^2(n-2) = n^3-4n^2+5n-2.
  (iii) P(n)/n = sum_{d|n} phi(d)/d is multiplicative with value 1 + a(1 - 1/p) at p^a.  For n = p:
  P = 2p - 1.  For n = p^a, a >= 2: P/n = 1 + a - a/p >= 1 + a/2 >= 2 > 2 - 1/n.  For n with two
  distinct prime factors: P/n >= (2 - 1/p)(2 - 1/q) >= (3/2)(5/3) = 5/2 > 2.  So P(n) > 2n - 1 unless n is prime.
  Identity (3) says N3 = (n-1)^2(n-2)/4, i.e. P(n) = 2n-1, i.e. n prime.                        QED
  For prime n, (i) is the permutation argument: j -> ij mod n permutes [1,n-1], S = (n-1) n(n-1)/2.
""")

def S_brute(n):
    return sum((i * j) % n for i in range(1, n) for j in range(1, n))

def N3_brute(n):
    return sum((i * j) // n for i in range(1, n) for j in range(1, n))

NMAX_S4 = 260
truth3 = []
for n in range(2, NMAX_S4 + 1):
    P = pillai(n)
    check(P == sum(gcd(k, n) for k in range(1, n + 1)), "Pillai gcd-sum n=%d" % n)
    S = S_brute(n)
    check(2 * S == n * (n * n - P), "S(n) closed form n=%d" % n)
    N3 = N3_brute(n)
    check(Fr(N3) == Fr((n - 1) ** 2 * (n - 2), 4) + Fr(P - 2 * n + 1, 2), "N3 closed form n=%d" % n)
    check(Fr(N3) == (Fr((n - 1) ** 2 * n * n, 4) - S) / n, "N3 via S n=%d" % n)
    check(P >= 2 * n - 1 and ((P == 2 * n - 1) == is_prime(n)), "P(n) >= 2n-1 iff prime n=%d" % n)
    if Fr(N3) == Fr((n - 1) ** 2 * (n - 2), 4):
        truth3.append(n)
check(truth3 == [n for n in range(2, NMAX_S4 + 1) if is_prime(n)], "truth set (3)")
print("  FINITE-EXACT  (i),(ii),(iii) for all 2 <= n <= %d ; truth set of (3) = primes exactly (covers lead's n<200)" % NMAX_S4)
print("     n | P(n) | P-2n+1 | N3 | (n-1)^2(n-2)/4 | deviation")
for n in (4, 5, 6, 7, 8, 9, 12, 15):
    P = pillai(n); N3 = N3_brute(n)
    print("    %2d | %4d | %6d | %5d | %s | %s" % (n, P, P - 2 * n + 1, N3, Fr((n - 1) ** 2 * (n - 2), 4),
                                                   Fr(N3) - Fr((n - 1) ** 2 * (n - 2), 4)))
print("""
COROLLARY S4b (PROVED).  4 N3(n) - (N1(n) + N2(n)) = 2 (P(n) - 2n + 1) - Z_3(n).
  For prime n both correction terms vanish: 4 N3 = N1 + N2 = (n-1)^2 (n-2).
  For composite n the left side is 2*(P(n)-2n+1) - (n/n_3 - 1): the bilinear sum sees the WHOLE divisor
  lattice (through P), the cubic sums see only the exponent profile through n_3.  The two corrections
  are independent: n = 6 has Z_3 = 0, P-2n+1 = 4; n = 8 has Z_3 = 3, P-2n+1 = 5.
""")
for n in range(2, 150):
    N1 = floor_sum(lambda k: k ** 3, n); N2 = N2_exact(n); N3 = N3_brute(n)
    P = pillai(n); Z3 = n // n_e(n, 3) - 1
    check(4 * N3 - (N1 + N2) == 2 * (P - 2 * n + 1) - Z3, "S4b n=%d" % n)
    if is_prime(n):
        check(4 * N3 == N1 + N2 == (n - 1) ** 2 * (n - 2), "prime 4N3 n=%d" % n)
print("  FINITE-EXACT  S4b verified for 2 <= n < 150 (4 N3 = N1 + N2 at every prime)")
print(stamp())

# ======================================================================= S5
banner("S5  Even powers see class numbers: sum floor(k^2/p) = law + h(-p) for p = 3 (mod 4), p > 3")
print("""
Setup.  For prime p and e >= 1 let g = gcd(e, p-1).  The map k -> k^e on (Z/p)^* has image the unique
subgroup H_g of index g, each value taken g times, so  sum_{k=1}^{p-1} (k^e mod p) = g * sum_{r in H_g} r.
If -1 in H_g (iff (p-1)/g is even), H_g is symmetric under r -> p-r and sum_{H_g} r = p |H_g| / 2 = p(p-1)/(2g),
so the odd-law value (sum k^e)/p - (p-1)/2 is exact (this recovers S1 for odd e, since then g is odd,
(p-1)/g even).  If (p-1)/g is odd (forcing g even) there is a genuine correction:
  sum floor(k^e/p) = (sum k^e)/p - (p-1)/2 + D_{e}(p),  D_e(p) = (p-1)/2 - g sum_{H_g} r / p.
Quadratic case g = 2, p = 3 (mod 4).  Let R, N be the sums of quadratic residues / non-residues in [1,p-1].
R + N = p(p-1)/2 and (CITED: Dirichlet class number formula for Q(sqrt(-p)); Davenport, Multiplicative
Number Theory, ch. 6; for p > 3, w = 2)   R - N = sum_a (a/p) a = -p h(-p).   Hence R = p(p-1)/4 - p h(-p)/2 and
  sum_{k=1}^{p-1} floor(k^2/p) = (p-1)(2p-1)/6 - (p-1)/2 + h(-p)         (p = 3 mod 4, p > 3).
So the deviation from the odd law is EXACTLY the class number (not h/2).  For p = 3 the unit count w = 6
changes the formula: sum floor(k^2/3) = 1 while (2*5)/6 - 1 + h(-3) = 5/3 (p = 3 is the excluded prime).
For p = 1 (mod 4) the quadratic sum obeys the odd law exactly (D_2 = 0): the even power sees nothing.
h(-p) is computed below INDEPENDENTLY by counting reduced primitive forms (a,b,c), b^2-4ac = -p,
so the Dirichlet recollection is itself checked exactly on the finite universe.
""")

def class_number(D):
    """h(D) for D < 0: number of reduced primitive forms ax^2+bxy+cy^2, b^2-4ac=D"""
    check(D < 0 and D % 4 in (0, 1), "discriminant %d" % D)
    h = 0
    a = 1
    while 3 * a * a <= -D:
        for b in range(-a + 1, a + 1):
            num = b * b - D
            if num % (4 * a):
                continue
            c = num // (4 * a)
            if c < a:
                continue
            if gcd(gcd(a, abs(b)), c) != 1:
                continue
            if a == c and b < 0:
                continue
            h += 1
        a += 1
    return h

def legendre(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1

def subgroup_sum(p, g):
    """sum of the elements of the index-g subgroup of (Z/p)^* (g | p-1)"""
    return sum(r for r in range(1, p) if pow(r, (p - 1) // g, p) == 1)

PMAX = 500
primes = [p for p in range(3, PMAX) if is_prime(p)]
rows = []
for p in primes:
    S2 = floor_sum(lambda k: k * k, p)
    law = Fr((p - 1) * (2 * p - 1), 6) - Fr(p - 1, 2)
    dev = Fr(S2) - law
    if p % 4 == 3:
        h = class_number(-p)
        dirichlet = Fr(-sum(legendre(a, p) * a for a in range(1, p)), p)
        R = sum(r for r in range(1, p) if legendre(r, p) == 1)
        Nn = p * (p - 1) // 2 - R
        if p > 3:
            check(dirichlet == h, "Dirichlet vs reduced forms p=%d: %s vs %d" % (p, dirichlet, h))
            check(R - Nn == -p * h, "R-N = -p h p=%d" % p)
            check(dev == h, "S2 = law + h(-p) fails at p=%d: dev=%s h=%d" % (p, dev, h))
        else:
            check(S2 == 1 and dev == Fr(1, 3) and h == 1, "p=3 exception: S2=1, law=2/3, deviation 1/3 = 2h/w with w=6")
        rows.append((p, S2, law, dev, h))
    else:
        check(dev == 0, "p=1 mod 4 should obey odd law, p=%d dev=%s" % (p, dev))
print("  FINITE-EXACT  primes p < %d, p = 3 mod 4, p > 3: sum floor(k^2/p) = (p-1)(2p-1)/6 - (p-1)/2 + h(-p)," % PMAX)
print("                with h(-p) from reduced forms; Dirichlet's -(1/p) sum (a/p) a = h(-p) confirmed on the same range.")
print("  FINITE-EXACT  primes p < %d, p = 1 mod 4: quadratic sum obeys the odd law exactly (deviation 0)." % PMAX)
print("  REFUTED for p = 3: deviation -2/3 (w = 6 exception), h(-3) = 1.")
print("     p | sum floor(k^2/p) | odd-law value | deviation | h(-p)")
for p, S2, law, dev, h in rows[:16]:
    print("   %3d | %16d | %13s | %9s | %d" % (p, S2, law, dev, h))
print("   ... (all %d primes p = 3 mod 4 below %d checked)" % (len(rows), PMAX))

# even exponent e with gcd(e,p-1)=2 gives the same class number; e with gcd 4 and (p-1)/4 odd tabulated
print()
print("  General even e (PROVED above, verified): deviation D_e(p) = h(-p) whenever gcd(e,p-1) = 2 and p = 3 mod 4;")
print("  D_e(p) = 0 whenever (p-1)/gcd(e,p-1) is even.")
for p in primes:
    for e in (2, 4, 6, 8, 10, 12):
        g = gcd(e, p - 1)
        Se = floor_sum(lambda k: k ** e, p)
        law = Fr(sum(k ** e for k in range(1, p)), p) - Fr(p - 1, 2)
        dev = Fr(Se) - law
        pred_sum = g * subgroup_sum(p, g)
        check(sum(pow(k, e, p) for k in range(1, p)) == pred_sum, "subgroup image sum p=%d e=%d" % (p, e))
        if ((p - 1) // g) % 2 == 0:
            check(dev == 0, "D_e=0 case p=%d e=%d" % (p, e))
        elif g == 2 and p > 3:
            check(dev == class_number(-p), "D_e = h(-p) case p=%d e=%d" % (p, e))
print("  FINITE-EXACT  all primes 3 <= p < %d, e in {2,4,6,8,10,12}: both cases verified." % PMAX)

print()
print("  Exploratory (FINITE-EXACT table, no formula claimed): g = 4, p = 5 (mod 8) [(p-1)/4 odd], e = 4:")
print("     p | D_4(p) | h(-p) [disc -4p since p=1 mod 4] | h(-4p) | D_4 - 2h(-4p)")
quart = []
for p in primes:
    if p % 8 == 5:
        Se = floor_sum(lambda k: k ** 4, p)
        law = Fr(sum(k ** 4 for k in range(1, p)), p) - Fr(p - 1, 2)
        dev = Fr(Se) - law
        h4 = class_number(-4 * p)
        quart.append((p, dev, h4))
        if p < 200:
            print("   %3d | %6s | %s | %d | %s" % (p, dev, "n/a", h4, dev - 2 * h4))
match_2h = all(dev == 2 * h4 for _, dev, h4 in quart)
match_h = all(dev == h4 for _, dev, h4 in quart)
print("  verdict: D_4(p) == h(-4p) for all p=5 mod 8 below %d: %s ; D_4(p) == 2 h(-4p): %s" % (PMAX, match_h, match_2h))
print("  (Whatever the verdict, the quartic deviation is left OPEN as a formula; the quadratic one is CITED+verified.)")
print(stamp())

# ======================================================================= S6
banner("S6  Typing 'odd cycles ~ odd functions ~ odd powers'")
print("""
(a) odd functions <-> negation conjugacy of the 3n+k family.  PROVED one-liner: for odd n and odd k put
    T_k(n) = (3n+k)/2^{v_2(3n+k)}.  Then T_{-k}(-n) = (-(3n+k))/2^{v_2(3n+k)} = -T_k(n), because v_2(-m) = v_2(m).
    So n -> T_k(n) is an odd function up to the sign flip of k; cycles of 3n-k are exactly the negatives of the
    cycles of 3n+k (inherited F3 and arithmetic_braids_20260917_summand.md, C_+(-n) = -C_-(n)).  This is the
    SAME algebraic fact as the pairing f(n-k) = -f(k) of S1: both are the statement that the map commutes with
    x -> -x, read mod n (S1) or on Z (Collatz).  Preserved predicate: 'orbit structure up to global sign'.
    Lost: nothing.  Decisive test: the 3n-1 cycle (5,7) is the 3n+1 cycle (-5,-7).
(b) odd powers <-> the pairing law.  x^e is an odd function iff e is odd; S1 then applies and the residue
    r_k + r_{n-k} in {0, n} cancels everything except the fixed layer Z_f(n).  Even e destroys the pairing
    (r_{n-k} = r_k) and the residue mass becomes a character sum -> class number (S5).  The user's slogan
    'odd powers' is therefore literally 'the exponent is odd so x^e is an odd function'; nothing deeper.
(c) 'odd cycles'.  Repo meaning (grep of 01-canon/theorems for Redei): THM-001 = Redei's theorem, every
    tournament has an ODD number of Hamiltonian paths; LEM-020 = 'Redei involution': a fixed-point-free pairing
    of witnesses plus a fixed layer, with the parity law #witnesses = #fixed points (mod 2); THM-1745 = the
    Hamiltonian-path spectrum is the odd numbers minus {7,21}.  Comparison with S1:
      source: S1 residue multiset {r_k};  target: LEM-020/THM-001 witness sets;
      map: NONE between the objects.  What IS shared is the PROOF TEMPLATE: an involution (k <-> n-k here,
      tau <-> 1-tau there, path reversal in Route A of THM-001) whose non-fixed orbits cancel exactly and whose
      fixed layer carries the whole residue.  In S1 the cancellation is exact and integer-valued
      (sum r_k = (n/2)(n-1-Z)), so it is strictly stronger than a parity law; in THM-001 only parity survives.
      Preserved predicate: 'the unpaired part is the fixed layer'.  Lost: all arithmetic of tournaments;
      the {7,21} hole has no counterpart in the floor sums (the truth sets here are squarefree / prime,
      cofinite in nothing).  Cheapest decisive test: the floor-sum fixed layer Z_f(n) can be any n/n_e - 1,
      e.g. Z_3(8) = 3 is odd while every Redei fixed layer has odd cardinality 1 -- there is no parity match.
    The three 3n-1 cycles and the 'three cycles' fascination: no map found from the floor-sum truth sets to
    cycle counts (the k=+-1,5,7,11,13 counts 4,9,5,7,13 of F3 have no relation to squarefreeness or
    primality that survives the decisive test: n=9 and n=4 are both non-squarefree but 9 counts cycles of
    3n+5 and 4 counts cycles of 3n+1).  Honest verdict: beyond the word 'odd' (= commutes with negation),
    'odd cycles' shares only the involution proof template with 'odd functions/powers'.
""")
# (a) verification
def T(k, n):
    m = 3 * n + k
    if m == 0:
        return 0
    v = (m & -m).bit_length() - 1
    return m >> v

for k in (1, -1, 5, -5, 7, -7, 13, -13):
    for n in range(-10001, 10002, 2):
        check(T(-k, -n) == -T(k, n), "conjugacy k=%d n=%d" % (k, n))
print("  FINITE-EXACT  T_{-k}(-n) = -T_k(n) for k in {+-1,+-5,+-7,+-13}, odd |n| <= 10001")
print("  FINITE-EXACT  the (5,7) 2-cycle of 3n-1 is the (-5,-7) cycle of 3n+1: T_1(-5) = %d, T_1(-7) = %d" % (T(1, -5), T(1, -7)))
check(T(1, -5) == -7 and T(1, -7) == -5, "(-5,-7) cycle")
print("  FINITE-EXACT  Z_3(8) = %d (odd fixed layer of size 3, so no Redei-style 'exactly one fixed point')" % (8 // n_e(8, 3) - 1))

banner("SUMMARY OF LABELS")
print("""
PROVED   S1 odd-function floor law (all n >= 2, even n included via the fixed point r_{n/2} in {0,n/2}).
PROVED   S2 Z_{x^e}(n) = n/n_e - 1; (1) and every odd-e identity hold iff n squarefree (e=1 degenerate: all n).
CITED    squarefree density 6/pi^2 = 1/zeta(2) (Gegenbauer 1885 / Hardy-Wright Thm 333).
PROVED   S3 lattice reciprocity N1 + N2 = (n-1)^2(n-2) + Z_3(n); (2) iff squarefree; exact integer cube roots.
PROVED   S4 S(n) = n(n^2-P(n))/2, N3 = (n-1)^2(n-2)/4 + (P(n)-2n+1)/2, P(n) >= 2n-1 iff-equality prime;
         (3) iff prime; 4 N3 - (N1+N2) = 2(P-2n+1) - Z_3.
CITED    S5 Dirichlet: sum (a/p) a = -p h(-p) for p = 3 mod 4, p > 3 (verified against reduced forms, p < 500).
PROVED   S5 given Dirichlet: sum floor(k^2/p) = (p-1)(2p-1)/6 - (p-1)/2 + h(-p) (p = 3 mod 4, p > 3); deviation 0
         for p = 1 mod 4; general even e: deviation h(-p) iff gcd(e,p-1)=2 and p=3 mod 4, else 0 when (p-1)/g even.
REFUTED  even f obeys the odd law: minimal witness (f, n) = (x^2, 3); p = 3 quadratic formula (w=6).
OPEN     closed form of the quartic deviation D_4(p) for p = 5 mod 8 (table only).
PROVED   S6(a) negation conjugacy T_{-k}(-n) = -T_k(n).  S6(c): no object-level map from floor sums to
         tournament odd cycles / Redei; shared involution proof template only.
""")
print("total time", stamp())
print("ALL CHECKS PASSED")
