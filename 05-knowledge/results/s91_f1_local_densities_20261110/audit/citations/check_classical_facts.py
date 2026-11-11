#!/usr/bin/env python3
"""Auditor 'citations' (S91): arithmetic behind the classical / literature statements quoted in
README section 3, THM-5128 (D), (E), (F), THM-5129 (A2) and HYP-9661 (D).  Exact integers, sympy, mpmath (50 digits).
No code of the audited session is imported."""
import itertools
import sys
from math import gcd, isqrt

import mpmath as mp
import sympy

mp.mp.dps = 50
OK = True


def check(name, cond, detail=""):
    global OK
    OK = OK and bool(cond)
    print(("PASS  " if cond else "FAIL  ") + name + ("   [" + str(detail) + "]" if detail != "" else ""))


# ---------------------------------------------------------------- README section 3: Euler-Sombor abstract
N = lambda i, j: i * i + i * j + j * j
check("norms i^2+ij+j^2 at (1,1),(1,2),(1,3),(2,3),(3,3) = 3,7,13,19,27; (2,2) gives 12",
      [N(*p) for p in [(1, 1), (1, 2), (1, 3), (2, 3), (3, 3)]] == [3, 7, 13, 19, 27] and N(2, 2) == 12)
s = lambda n: mp.sqrt(n)
eps = 2 * s(19) - 5 * s(3)
delta = s(13) + 2 * s(3) - s(7) - s(19)
check("eps = 2 sqrt19 - 5 sqrt3 = sqrt76 - sqrt75 > 0 and eps < delta (paper: 0 < eps < delta)",
      abs(eps - (s(76) - s(75))) < mp.mpf(10) ** -40 and 0 < eps < delta, "eps=%s delta=%s" % (mp.nstr(eps, 8), mp.nstr(delta, 8)))
check("paper's penalties: a = sqrt7+sqrt19-4sqrt3 = 0.076447, b = sqrt13-2sqrt3 = 0.141450, delta = b - a",
      abs((s(7) + s(19) - 4 * s(3)) - mp.mpf("0.076447")) < 1e-6 and abs((s(13) - 2 * s(3)) - mp.mpf("0.141450")) < 1e-6
      and abs(delta - ((s(13) - 2 * s(3)) - (s(7) + s(19) - 4 * s(3)))) < mp.mpf(10) ** -40)
check("1039 = 2+5+13+33+89+240+657 (unicyclic graphs on 4..10 vertices, OEIS A001429)", 2 + 5 + 13 + 33 + 89 + 240 + 657 == 1039)
check("204990 = 11+23+47+106+235+551+1301+3159+7741+19320+48629+123867 (trees on 7..18 vertices, OEIS A000055)",
      sum([11, 23, 47, 106, 235, 551, 1301, 3159, 7741, 19320, 48629, 123867]) == 204990)

# ---------------------------------------------------------------- README section 3: the 2-clock
def order(a, m):
    k, x = 1, a % m
    while x != 1:
        x = x * a % m
        k += 1
    return k


ords = [order(2, p) for p in (3, 5, 7, 13, 19, 27, 1039)]
check("ord_p(2) = 2, 4, 3, 12, 18, 18, 519 for p = 3, 5, 7, 13, 19, 27, 1039", ords == [2, 4, 3, 12, 18, 18, 519], ords)
check("1039 is prime; 1039 = 7 mod 8; 1039 + 1 = 2^4 * 5 * 13", sympy.isprime(1039) and 1039 % 8 == 7 and 1040 == 16 * 5 * 13)
check("-1 is a power of 2 modulo 3, 5, 13, 19, 27 and not modulo 7, 1039",
      all(any(pow(2, k, p) == p - 1 for k in range(1, 600)) for p in (3, 5, 13, 19, 27))
      and not any(pow(2, k, p) == p - 1 for p in (7, 1039) for k in range(1, 600)))
check("5*13 = 2^6+1, 27*19 = 2^9+1, 7 = 2^3-1, 21 = (2^6-1)/3 = 111 in base 4", 65 == 2 ** 6 + 1 and 513 == 2 ** 9 + 1 and 7 == 2 ** 3 - 1 and 21 == 16 + 4 + 1)
check("x^2 + 2^s y^2 - [s=0] with x, y odd is 1,3,5,1,1,.. mod 8 (never 7)",
      all(((x * x + 2 ** t * y * y - (t == 0)) % 8) == [1, 3, 5, 1, 1, 1][t] for t in range(6) for x in (1, 3, 5, 7) for y in (1, 3, 5, 7)))
X, Y, x, y = sympy.symbols("X Y x y")
check("3x^2+9xy+7y^2 = X^2+XY+Y^2 at (X,Y) = (x+y, x+2y); discriminants -3; x^2+xy+2y^2 has discriminant -7",
      sympy.expand((x + y) ** 2 + (x + y) * (x + 2 * y) + (x + 2 * y) ** 2 - (3 * x * x + 9 * x * y + 7 * y * y)) == 0
      and 81 - 84 == -3 and 1 - 8 == -7)
print("      discriminants of x^2 + 2^s y^2, s = 0..5:", [-4 * 2 ** t for t in range(6)],
      " (only s = 0, 1 give the fundamental discriminants -4, -8; s >= 2 are non-maximal orders of conductor 2^floor(s/2))")

# ---------------------------------------------------------------- THM-5129 (A2): Busch's formula against the stated f
def busch(n):
    r = n % 3
    if r == 1:
        return 5 ** ((n - 1) // 3)
    if r == 0:
        return 3 * 5 ** ((n - 3) // 3)
    return 9 * 5 ** ((n - 5) // 3)


check("Busch/Moon formula gives f(3..10) = 3, 5, 9, 15, 25, 45, 75, 125 (the values used in THM-5129 and README 4.5)",
      [busch(n) for n in range(3, 11)] == [3, 5, 9, 15, 25, 45, 75, 125], [busch(n) for n in range(3, 11)])
print("      formula at n = 2 (no strong tournament of order 2): 9 * 5^(-1) =", 9 / 5, " -> the range n >= 3 (or n = 1) must be stated")
check("the formula is strictly increasing for n >= 3 and f(n) >= f(n-1) + 2 there", all(busch(n) >= busch(n - 1) + 2 for n in range(4, 60)))
# exhaustive f(n) for n = 3..6 by brute force (strong tournaments, Hamiltonian paths by DP)
def ham_paths(n, adj):
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for v in range(n):
            c = dp[mask][v]
            if c:
                for w in range(n):
                    if not mask >> w & 1 and adj[v] >> w & 1:
                        dp[mask | 1 << w][w] += c
    return sum(dp[full])


def strong(n, adj):
    def reach(a):
        seen = 1
        stack = [0]
        while stack:
            v = stack.pop()
            m = a[v] & ~seen
            while m:
                w = (m & -m).bit_length() - 1
                m &= m - 1
                seen |= 1 << w
                stack.append(w)
        return seen == (1 << n) - 1
    radj = [0] * n
    for v in range(n):
        for w in range(n):
            if adj[v] >> w & 1:
                radj[w] |= 1 << v
    return reach(adj) and reach(radj)


fmin = {}
for n in range(3, 7):
    pairs = list(itertools.combinations(range(n), 2))
    best = None
    for bits in range(1 << len(pairs)):
        adj = [0] * n
        for k, (i, j) in enumerate(pairs):
            if bits >> k & 1:
                adj[i] |= 1 << j
            else:
                adj[j] |= 1 << i
        if strong(n, adj):
            h = ham_paths(n, adj)
            if best is None or h < best:
                best = h
    fmin[n] = best
check("exhaustive: least H over strong tournaments of order 3..6 is 3, 5, 9, 15", [fmin[n] for n in range(3, 7)] == [3, 5, 9, 15], fmin)
# nearly transitive tournament (one reversed arc of span n-1): H = 2^(n-2) + 1  (Moon 1972, quoted by Busch 2006)
def nearly_transitive(n):
    adj = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if (i, j) == (0, n - 1):
                adj[j] |= 1 << i
            else:
                adj[i] |= 1 << j
    return adj


check("transitive tournament with the arc v_1 v_n reversed has 2^(n-2) + 1 Hamiltonian paths, n = 3..12",
      all(ham_paths(n, nearly_transitive(n)) == 2 ** (n - 2) + 1 for n in range(3, 13)))

# ---------------------------------------------------------------- THM-5128 (D): dilogarithm constants
Li2 = lambda z: mp.polylog(2, z)
l = lambda j: -Li2(-mp.mpf(2) ** (-j))
check("Li_2(-1) = -pi^2/12, so l_0 = pi^2/12", abs(l(0) - mp.pi ** 2 / 12) < mp.mpf(10) ** -40)
k0 = l(0) - 2 * l(1) + l(2)
k1 = l(1) - l(2)
check("kappa_0 = pi^2/12 + 2 Li2(-1/2) - Li2(-1/4) = 0.1615389173", abs(k0 - (mp.pi ** 2 / 12 + 2 * Li2(-0.5) - Li2(-0.25))) < mp.mpf(10) ** -40
      and abs(k0 - mp.mpf("0.1615389173")) < 5e-11, mp.nstr(k0, 12))
check("kappa_1 = Li2(-1/4) - Li2(-1/2) = 0.2125139092", abs(k1 - mp.mpf("0.2125139092")) < 5e-11, mp.nstr(k1, 12))
g = 2 * (k0 + k1)
check("g_bar = pi^2/6 + 2 Li2(-1/2) = 0.7481056530", abs(g - (mp.pi ** 2 / 6 + 2 * Li2(-0.5))) < mp.mpf(10) ** -40 and abs(g - mp.mpf("0.7481056530")) < 5e-11, mp.nstr(g, 12))
check("triple shares w_0 = 0.43186124, w_1 = 0.56813876", abs(k0 / (k0 + k1) - mp.mpf("0.43186124")) < 5e-9 and abs(k1 / (k0 + k1) - mp.mpf("0.56813876")) < 5e-9,
      (mp.nstr(k0 / (k0 + k1), 10), mp.nstr(k1 / (k0 + k1), 10)))
tot = mp.nsum(lambda t: (t + 1) * (l(t) - 2 * l(t + 1) + l(t + 2)), [0, mp.inf])
check("sum_{s>=0} (s+1)(l_s - 2 l_{s+1} + l_{s+2}) = pi^2/12", abs(tot - mp.pi ** 2 / 12) < mp.mpf(10) ** -25, mp.nstr(tot, 15))
check("word-pair shares w_0 = kappa_0/(pi^2/12) = 0.19640777, w_1 = 0.80359223",
      abs(k0 / (mp.pi ** 2 / 12) - mp.mpf("0.19640777")) < 5e-9 and abs(1 - k0 / (mp.pi ** 2 / 12) - mp.mpf("0.80359223")) < 5e-9,
      mp.nstr(k0 / (mp.pi ** 2 / 12), 10))
# the integral identity of (D1)
for t in range(0, 4):
    F = lambda r, t=t: mp.log((r + 2 ** t) / (r / 2 + 2 ** t)) / r
    I = mp.quad(F, [mp.mpf(1) / 2, 1])
    if abs(I - (l(t) - 2 * l(t + 1) + l(t + 2))) > mp.mpf(10) ** -30:
        check("integral identity int_{1/2}^1 F_s = l_s - 2 l_{s+1} + l_{s+2} at s = %d" % t, False)
        break
else:
    check("integral identity int_{1/2}^1 (1/r) ln((r+2^s)/(r/2+2^s)) dr = l_s - 2 l_{s+1} + l_{s+2}, s = 0..3", True)
check("prod over odd p of (1 - 1/p^2) = 8/pi^2 = 0.810569...", abs(8 / mp.pi ** 2 - mp.mpf("0.8105694691")) < 1e-9, mp.nstr(8 / mp.pi ** 2, 12))
check("(1 - 1/p^2 + 1/p)/(1 - 1/p^2) = 1 + p/(p^2-1) and orthogonality probability 1/(p+1), p = 3..97",
      all(sympy.Rational(p * p - 1 + p, p * p) / sympy.Rational(p * p - 1, p * p) == 1 + sympy.Rational(p, p * p - 1)
          and sympy.Rational(p - 1, p * p - 1) == sympy.Rational(1, p + 1) for p in sympy.primerange(3, 98)))

# ---------------------------------------------------------------- THM-5128 (E): the symbol over F_2
z = sympy.symbols("z")
P = sympy.Poly((1 + z) ** 2 * (1 + z + z ** 2), z, modulus=2)
check("(1+z)^2 (1+z+z^2) = 1 + z + z^3 + z^4 over F_2 (= z^2 (z + z^2 + z^-1 + z^-2)) and divides z^6 - 1",
      P == sympy.Poly(1 + z + z ** 3 + z ** 4, z, modulus=2) and sympy.Poly(z ** 6 - 1, z, modulus=2).rem(P).is_zero)
# radical inverse on B_n: each value (2i+1)/2^n once
def phi_num(h, n):  # phi(h) * 2^n with phi(h) = sum h_i 2^(-i-1)
    return sum(((h >> i) & 1) << (n - i - 1) for i in range(n))


check("binary radical inverse on B_n takes each value (2i+1)/2^n exactly once, n = 1..14",
      all(sorted(phi_num(h, n) for h in range(1 << (n - 1), 1 << n)) == list(range(1, 1 << n, 2)) for n in range(1, 15)))
# q(h) recursion and its parity form (Rudin-Shapiro type)
def q(h):
    if h == 0:
        return 0
    if h % 2 == 0:
        return q(h // 2)
    k = h // 2
    return q(k) + (k + 1) // 2


def qpar(h):
    b = [(h >> i) & 1 for i in range(h.bit_length() + 3)]
    return sum(b[i] * (b[i + 1] + b[i + 2]) for i in range(h.bit_length())) % 2


check("q(h) = sum_i h_i (h_{i+1} + h_{i+2}) mod 2 for h < 2^14; 0 <= q(h) <= h/2", all(q(h) % 2 == qpar(h) and 0 <= 2 * q(h) <= h for h in range(1 << 14)))
xs = lambda h: 2 * h + 1 - 2 * q(h)
check("x(h) = y(h) iff h = 0 or h = 2^j, and then x = y = 2^(j+1) + 1 (h < 2^14)",
      all((xs(h) == 2 * h + 1) == (h == 0 or h & (h - 1) == 0) for h in range(1 << 14)) and all(xs(1 << j) == 2 ** (j + 1) + 1 for j in range(13)))
first = {}
for p in (3, 5, 7, 11, 13, 17, 19, 23, 29, 37, 41, 43):
    h = 1
    while not (xs(h) % p == 0 and (2 * h + 1) % p == 0):
        h += 1
    first[p] = h
check("first p-divisible states: 3:1, 5:2, 7:38, 11:16, 13:32, 17:8, 19:218, 23:103, 29:1435, 37:1350, 41:266, 43:64",
      first == {3: 1, 5: 2, 7: 38, 11: 16, 13: 32, 17: 8, 19: 218, 23: 103, 29: 1435, 37: 1350, 41: 266, 43: 64}, first)
print("ALL CHECKS PASSED" if OK else "SOME CHECKS FAILED")
sys.exit(0 if OK else 1)
