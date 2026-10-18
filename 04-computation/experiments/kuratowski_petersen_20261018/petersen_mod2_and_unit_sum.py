#!/usr/bin/env python3
"""(a) Reduction of the golden Petersen realization at the prime of Q(zeta20) above 2 (residue field F_16):
    zeta5 -> zeta (a root of Phi_5 over F_2), i -> 1, zeta20 = i zeta5^4 -> zeta^4, sqrt5 = 2(zeta5 + zeta5^4) + 1 -> 1,
    so R = (zeta20 + zeta20^-1)/sqrt5 -> zeta^4 + zeta = w (a primitive cube root of unity), phi -> 1 + w = w^2,
    r = R/phi -> w^2. Outer o_k -> w zeta^k, inner i_k -> w^2 zeta^k. Check: these ten points are exactly the
    non-neighbours of 0 in G_4 = Cay(F_16, mu_5) (= F_16^* minus mu_5), and the fifteen Petersen edges map to edges.
(b) Characteristic zero: on the ten two-sums s_ab = zeta5^a + zeta5^b of the regular pentagon, s_ab + s_cd is a unit
    vector iff {a,b} and {c,d} are disjoint (then the sum is -zeta5^e). So the Petersen graph is the unit-SUM graph of
    the two-sums; mod 2, sum = difference and it becomes the unit-distance second subconstituent of (a)."""
import itertools, cmath, math
# F_16 = F_2[x]/(x^4 + x^3 + x^2 + x + 1)? Phi_5 is irreducible over F_2; use it, so x = zeta has order 5.
MOD = 0b11111
def mul(a, b):
    r = 0
    for i in range(4):
        if b >> i & 1: r ^= a << i
    for d in range(7, 3, -1):
        if r >> d & 1: r ^= MOD << (d - 4)
    return r
def pw(a, n):
    r = 1
    while n:
        if n & 1: r = mul(r, a)
        a = mul(a, a); n >>= 1
    return r
zeta = 0b10
mu5 = [pw(zeta, k) for k in range(5)]
assert len(set(mu5)) == 5 and pw(zeta, 5) == 1
w = mu5[1] ^ mu5[4]                        # zeta + zeta^4
assert pw(w, 3) == 1 and w != 1
phi = 1 ^ w
assert phi == mul(w, w)
outer = [mul(w, z) for z in mu5]; inner = [mul(mul(w, w), z) for z in mu5]
pts = outer + inner
E = [(k, (k + 1) % 5) for k in range(5)] + [(k, k + 5) for k in range(5)] + [(5 + k, 5 + (k + 2) % 5) for k in range(5)]
N = lambda x: pw(x, 5)                     # norm F_16 -> F_4
nonnbr = [x for x in range(1, 16) if N(x) != 1]
print('(a) ten distinct reduced points: %s ; equal to the non-neighbours of 0 in G_4: %s ; all 15 edges map to unit (norm-1) differences: %s'
      % (len(set(pts)) == 10, set(pts) == set(nonnbr), all(N(pts[a] ^ pts[b]) == 1 for a, b in E)))
# (b)
z5 = [cmath.exp(2j * math.pi * k / 5) for k in range(5)]
pairs = list(itertools.combinations(range(5), 2))
ok = True
for P, Q in itertools.combinations(pairs, 2):
    s = z5[P[0]] + z5[P[1]] + z5[Q[0]] + z5[Q[1]]
    unit = abs(abs(s) - 1) < 1e-12
    ok &= unit == (not set(P) & set(Q))
print('(b) for all 45 pairs of two-sums: s + s\' is a unit vector exactly when the index pairs are disjoint:', ok)
