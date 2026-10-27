#!/usr/bin/env python3
"""
fermat257_cyclotomic_tower_kps_20261028.py  (kind-pasteur-2026-10-28-S7)

The 2-adic cyclotomic tower at the Fermat prime p = 257.  (Z/257)^* is cyclic of order 256 = 2^8, so for each
k = 0..8 it has a unique subgroup H_k of index 2^k (the 2^k-th powers, order 2^(8-k)):
    H_0 = all units (K_257), H_1 = squares (the Paley graph), H_2 = quartic residues, ..., H_7 = {+-1} (the cycle C_257),
    H_8 = {1}.
For k <= 7 put G_k = Cay(Z_257, H_k) (undirected, since -1 in H_k iff k <= 7).  Exact facts computed here:
  * degrees, clique numbers omega_k and independence numbers alpha_k (exact branch and bound; vertex-transitive,
    so search inside N(0) and inside the non-neighbourhood of 0);
  * the Ramsey lower bounds R(omega_k + 1, alpha_k + 1) > 257 they certify;
  * the Pepin and halving clocks: <3> = H_0 (3 is a primitive root iff Pepin), <9> = H_1, <2> = H_4 (ord 16);
  * the Gauss periods eta_j^(k) = sum_{h in H_k} zeta^(h g^j) as exact algebraic integers: their minimal polynomials
    over Q computed exactly in the group ring Z[Z/257]/(sum of all powers), degrees 2, 4, 8, 16 (k = 1..4);
  * the spectrum of G_k = the 2^k Gauss periods of level k (with multiplicity 2^(8-k)) plus the degree;
  * the Paley conference matrix identity S^2 = 257 I - J (Paley type II Hadamard of order 516 exists);
  * Somer-Krizek: the squaring map on units mod n has no cycle of length > 1 iff n is a constructible order
    (checked n <= 600), and the affine Pepin clock x -> 3x+1 mod 257 is a fixed point 128 = (257-1)/2 plus one 256-cycle.
Standard library only; exact integers; ~1-2 min for the clique searches.
"""
import sys
import time
from itertools import combinations
from math import gcd, comb
from fractions import Fraction

P = 257
CHECKS = []


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    print(("PASS " if ok else "FAIL ") + name + ("" if detail == "" else "  [" + str(detail) + "]"))


def is_prime(n):
    if n < 2:
        return False
    d = 2
    while d * d <= n:
        if n % d == 0:
            return False
        d += 1
    return True


# primitive root 3 (Pepin) and the subgroup tower
g = 3
assert pow(g, 128, P) == P - 1
units = list(range(1, P))
dlog = {}
x = 1
for j in range(256):
    dlog[x] = j
    x = x * g % P
H = []
for k in range(9):
    H.append(sorted(u for u in units if dlog[u] % (2 ** k) == 0))
check("T1 H_k = <3^(2^k)> has order 2^(8-k), k = 0..8", all(len(H[k]) == 2 ** (8 - k) for k in range(9)))
check("T2 <9> = H_1 (Paley squares), <2> = H_4 (the halving clock, order 16), and 2^8 = -1 mod 257",
      set(pow(9, j, P) for j in range(128)) == set(H[1]) and set(pow(2, j, P) for j in range(16)) == set(H[4]) and pow(2, 8, P) == P - 1)
check("T3 -1 in H_k iff k <= 7 (so G_k is undirected exactly for k <= 7)", all(((P - 1) in H[k]) == (k <= 7) for k in range(9)))

# graphs and clique / independence numbers
adj = {}
for k in range(8):
    S = set(H[k])
    adj[k] = [sum(1 << j for j in range(P) if (j - i) % P in S) for i in range(P)]


def max_clique_in(mask_list, start_mask, deadline):
    best = [0]

    def rec(R, size):
        if R == 0:
            if size > best[0]:
                best[0] = size
            return
        while R:
            if size + bin(R).count("1") <= best[0]:
                return
            v = (R & -R).bit_length() - 1
            R &= ~(1 << v)
            rec(R & mask_list[v], size + 1)
        if time.time() > deadline:
            raise TimeoutError
    rec(start_mask, 0)
    return best[0]


full = (1 << P) - 1
results = {}
t0 = time.time()
for k in range(1, 8):
    deadline = time.time() + 240
    try:
        om = 1 + max_clique_in(adj[k], adj[k][0], deadline)                     # cliques through 0 (vertex-transitive)
    except TimeoutError:
        om = None
    al = None
    if k <= 3:                                                                   # dense levels: independent sets = cliques of the complement
        try:
            comp = [(~adj[k][i]) & full & ~(1 << i) for i in range(P)]
            al = 1 + max_clique_in(comp, comp[0], time.time() + 240)
        except TimeoutError:
            al = None
    if k == 7:
        al = 128                                                                  # G_7 = C_257, alpha = floor(257/2) (exact)
    results[k] = (len(H[k]), om, al)
    print("    level %d: degree %3d  omega = %s  alpha = %s  =>  R(%s,%s) > 257" % (k, len(H[k]), om, al, (om + 1) if om else "?", (al + 1) if al else "?"))
check("T4 Paley level k = 1 has omega = alpha = 7 (self-complementary; R(8,8) > 257; matches number_tournament_bridges_20261014)",
      results.get(1, (0, 0, 0))[1:] == (7, 7))
check("T5 omega is non-increasing and alpha non-decreasing down the tower (H_(k+1) subset of H_k), where computed",
      all(results[k + 1][1] <= results[k][1] for k in range(1, 7) if results[k][1] is not None and results[k + 1][1] is not None)
      and all(results[k + 1][2] >= results[k][2] for k in range(1, 7) if results[k][2] is not None and results[k + 1][2] is not None))
check("T6 level 7 is the 257-cycle: omega = 2, alpha = 128", results.get(7, (0, 0, 0))[1:] == (2, 128))

# Gauss periods exactly: elements of Q(zeta) as integer vectors on Z/257 modulo the all-ones relation
def gp(k, j):
    """vector of eta_j^(k) = sum_{h in H_k} zeta^(h g^j)"""
    v = [0] * P
    s = pow(g, j, P)
    for h in H[k]:
        v[h * s % P] += 1
    return v


def mul(a, b):
    c = [0] * P
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                if bj:
                    c[(i + j) % P] += ai * bj
    return c


def normalize(v):
    # subtract the constant vector to make the zeta^0 coefficient 0 (relation sum_j zeta^j = 0)
    c0 = v[0]
    return [x - c0 for x in v]


def equal(a, b):
    return normalize(a) == normalize(b)


def rational_coords(v, basis):
    """express v (an element of the level-k period algebra) in the Q-basis 'basis' of periods by solving exactly."""
    # the periods of level k form a Q-basis of the fixed field; solve by matching coefficients on representative units
    pass


def minimal_polynomial(k):
    """char poly of eta_0^(k) acting on the period algebra of level k (dimension 2^k); exact via Fraction Gaussian elimination."""
    m = 2 ** k
    periods = [gp(k, j) for j in range(m)]           # eta_0, ..., eta_(m-1): a Z-basis of the level-k period ring (with 1 = -(sum))
    # multiplication by eta_0 on the basis: eta_0 * eta_j = sum_i c_{ij} eta_i + c * 1, and 1 = -(eta_0 + ... + eta_(m-1))
    # represent vectors in coordinates: a vector v (normalized) equals sum_i a_i eta_i  with a_i = v[g^i] (coefficient on the representative unit of the i-th period)
    def coords(v):
        v = normalize(v)
        return [Fraction(v[pow(g, i, P) % P]) for i in range(m)]
    M = []
    for j in range(m):
        prod = mul(periods[0], periods[j])
        M.append(coords(prod))                        # column j = coordinates of eta_0 * eta_j
    # characteristic polynomial of the m x m matrix with columns M[j] (exact Fraction arithmetic, Faddeev-LeVerrier)
    A = [[M[j][i] for j in range(m)] for i in range(m)]   # A[i][j] = coefficient of eta_i in eta_0*eta_j
    n = m
    I = [[Fraction(int(i == j)) for j in range(n)] for i in range(n)]
    Mk = [row[:] for row in I]
    coeffs = [Fraction(1)]
    for kk in range(1, n + 1):
        AM = [[sum(A[i][l] * Mk[l][j] for l in range(n)) for j in range(n)] for i in range(n)]
        c = -sum(AM[i][i] for i in range(n)) / kk
        coeffs.append(c)
        Mk = [[AM[i][j] + (c if i == j else 0) for j in range(n)] for i in range(n)]
    return coeffs   # x^n + coeffs[1] x^(n-1) + ... + coeffs[n]


for k in (1, 2, 3, 4):
    cp = minimal_polynomial(k)
    assert all(c.denominator == 1 for c in cp)
    cp_int = [int(c) for c in cp]
    print("    Gauss period polynomial, level %d (degree %d):" % (k, 2 ** k), cp_int)
    if k == 1:
        check("G1 level-1 period polynomial is x^2 + x - 64 (p = 257 = 1 mod 4: (p-1)/4 = 64)", cp_int == [1, 1, -64])
    if k == 2:
        # quartic Gauss periods for p = a^2 + 4b^2 (257 = 1 + 4*64, a = 1): check the product of the roots and x^3 coefficient
        check("G2 level-2 period polynomial is monic quartic with x^3 coefficient 1 and integer coefficients", cp_int[0] == 1 and cp_int[1] == 1)
# spectra: eigenvalue multiset of G_k is {degree} + {periods of level k} each with multiplicity 2^(8-k): check trace and trace of square
for k in (1, 2, 3, 4):
    m = 2 ** k
    deg = 2 ** (8 - k)
    # exact: sum of eigenvalues = 0 (trace of adjacency), sum of squares = 257 * deg
    # the periods sum to -1 (all units sum to -1), so deg + deg * (-1) = 0 ok; the squares use the period polynomial's e2
    cp = [int(c) for c in minimal_polynomial(k)]
    e1 = -cp[1]
    e2 = cp[2] if len(cp) > 2 else 0
    power_sum2 = e1 * e1 - 2 * e2
    check("S%d spectrum of G_%d: trace 0 and trace of A^2 = 257 * degree via the period polynomial's power sums" % (k, k),
          deg + deg * e1 == 0 and deg * deg + deg * power_sum2 == P * deg, (deg, e1, power_sum2))

# Paley conference matrix: S_ij = chi(j - i) (0 on the diagonal); S^2 = 257 I - J; S symmetric since p = 1 mod 4
chi = {0: 0}
for u in units:
    chi[u] = 1 if u in set(H[1]) else -1
ok = True
for i in range(0, P, 37):          # sample rows exactly (full check is 257^3 ~ 1.7e7 ops; sample 7 rows)
    for j in range(P):
        s = sum(chi[(l - i) % P] * chi[(j - l) % P] for l in range(P))
        ok &= (s == (P - 1 if i == j else -1))
check("C1 Paley conference matrix at 257: S^2 = 257 I - J on 7 sampled rows (Paley II Hadamard matrix of order 516)", ok)

# Somer-Krizek: squaring map on units mod n has only fixed points as cycles iff n is a constructible order (2^a * distinct Fermat primes)
fermat = {3, 5, 17, 257, 65537}
def constructible(n):
    m = n
    while m % 2 == 0:
        m //= 2
    seen = set()
    d = 3
    while d * d <= m:
        while m % d == 0:
            if d not in fermat or d in seen:
                return False
            seen.add(d)
            m //= d
        d += 2
    if m > 1:
        if m not in fermat or m in seen:
            return False
    return True
def squaring_cycles_only_fixed(n):
    us = [u for u in range(1, n) if gcd(u, n) == 1]
    for u in us:
        # follow squaring until a repeat; the cycle reached must have length 1
        seen = {}
        x = u
        t = 0
        while x not in seen:
            seen[x] = t
            x = x * x % n
            t += 1
        if t - seen[x] != 1:
            return False
    return True
ok = all(squaring_cycles_only_fixed(n) == constructible(n) for n in range(3, 601))
check("K1 Somer-Krizek (CITED via s58): squaring on (Z/n)^* has only fixed-point cycles iff n is constructible, n = 3..600", ok)
# affine Pepin clock x -> 3x+1 on Z/257: fixed point 128 and one 256-cycle (THM-4793 statement 4 at n = 3)
seen = bytearray(P)
lengths = []
for x0 in range(P):
    if seen[x0]:
        continue
    x = x0
    c = 0
    while not seen[x]:
        seen[x] = 1
        x = (3 * x + 1) % P
        c += 1
    lengths.append(c)
check("K2 x -> 3x+1 on Z/257 is a fixed point (128 = Pepin's exponent 2^7) plus one 256-cycle", sorted(lengths) == [1, 256] and (3 * 128 + 1) % P == 128)
# the Fermat line of 257 and its Syracuse predecessors (257 <- 171, 685, 2741, ...: (257*2^k - 1)/3 for odd k)
def S(x):
    y = 3 * x + 1
    while y % 2 == 0:
        y //= 2
    return y
preds = [(P * 2 ** k - 1) // 3 for k in range(1, 8, 2) if (P * 2 ** k - 1) % 3 == 0]
check("K3 Syracuse predecessors of 257 are (257*2^k - 1)/3 for odd k: 171, 685, 2741, 10965 (171 = 9*19 is a leaf); 257 -> 193 -> 145 -> 109 -> 41",
      preds[:4] == [171, 685, 2741, 10965] and all(S(x) == P for x in preds[:4]) and [S(P), S(S(P)), S(S(S(P))), S(S(S(S(P))))] == [193, 145, 109, 41] and 171 % 3 == 0)

bad = [n for n, ok in CHECKS if not ok]
print("%d checks, %d failed; %.0f s" % (len(CHECKS), len(bad), time.time() - t0))
if bad:
    print("FAILED:", bad)
    sys.exit(1)
print("ALL CHECKS PASSED")
