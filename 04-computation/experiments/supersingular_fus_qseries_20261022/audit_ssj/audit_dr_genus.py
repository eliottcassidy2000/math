#!/usr/bin/env python3
"""Independent audit of section 4 (Deligne-Rapoport genus identities, Eichler's H(l, M), cusp
multiplicities) with an own implementation of the genus formula of X_0(N) (no PARI mfdim). (audit_ssj)
Checks, including the primes l = 2, 3 (Kronecker-symbol conventions (-4/2) = 0, (-3/3) = 0):
  (G1) g(X_0(pq)) = 2 g(X_0(q)) + H(p,q) - 1 = 2 g(X_0(p)) + H(q,p) - 1   for all prime pairs p != q < 200
  (G2) g(X_0(p^2 q)) = 2 g(X_0(p^2)) + H(q, p^2) - 1                       for p^2 q <= 60000
  (G3) the same DR identity for X_0(l M) at every prime l not dividing M, M <= 300 (general form)
  (G4) the closed forms for H(p,q), H(q,p^2) and the example rows N = 55, 275, 637
  (G5) Katz-Mazur multiplicities phi(l^min(i,a-i)) equal the cusp weights phi(gcd(d, l^a/d))
Run: python3 -B audit_dr_genus.py > audit_dr_genus.out
"""
from fractions import Fraction as Fr
from math import gcd
from sympy import factorint, primerange, divisors, totient, jacobi_symbol

bad = 0


def fail(m):
    global bad
    bad += 1
    print("FAIL:", m)


def kron(a, n):
    """Kronecker symbol (a/n) for n a positive prime (own; only a = -4, -3, -p needed)"""
    if n == 2:
        if a % 2 == 0:
            return 0
        return 1 if a % 8 in (1, 7) else -1
    if a % n == 0:
        return 0
    return 1 if pow(a % n, (n - 1) // 2, n) == 1 else -1


def psi(N):
    r = Fr(N)
    for p in factorint(N):
        r *= Fr(p + 1, p)
    return r


def nu2(N):
    if N % 4 == 0:
        return 0
    r = 1
    for p in factorint(N):
        r *= 1 + kron(-4, p)
    return r


def nu3(N):
    if N % 9 == 0:
        return 0
    r = 1
    for p in factorint(N):
        r *= 1 + kron(-3, p)
    return r


def cusps(N):
    return sum(int(totient(gcd(d, N // d))) for d in divisors(N))


def genus(N):
    g = 1 + psi(N) / 12 - Fr(nu2(N), 4) - Fr(nu3(N), 3) - Fr(cusps(N), 2)
    assert g.denominator == 1
    return int(g)


def H(l, M):
    return (l - 1) * psi(M) / 12 + Fr(nu2(M), 4) * (1 - kron(-4, l)) + Fr(nu3(M), 3) * (1 - kron(-3, l))


# sanity: known genera
known = {11: 1, 23: 2, 37: 2, 55: 5, 49: 1, 25: 0, 27: 1, 14: 1, 64: 3, 100: 7, 637: 57, 275: 25, 46: 5, 22: 2, 28: 2}
for N, g in known.items():
    if genus(N) != g:
        fail(f"genus({N}) = {genus(N)} != {g}")
P = list(primerange(2, 200))
n1 = 0
for p in P:
    for q in P:
        if p == q:
            continue
        n1 += 1
        hpq = H(p, q)
        if hpq.denominator != 1:
            fail(f"H({p},{q}) not integral")
        if genus(p * q) != 2 * genus(q) + hpq - 1:
            fail(f"G1 {p},{q}")
n2 = 0
for p in P:
    for q in P:
        if p != q and p * p * q <= 60000:
            n2 += 1
            if genus(p * p * q) != 2 * genus(p * p) + H(q, p * p) - 1:
                fail(f"G2 {p},{q}")
n3 = 0
for M in range(1, 301):
    for l in primerange(2, 60):
        if M % l == 0:
            continue
        n3 += 1
        g1 = genus(M) if M > 1 else 0
        if genus(l * M) != 2 * g1 + H(l, M) - 1:
            fail(f"G3 l={l} M={M}")
print(f"(G1) {n1} ordered prime pairs < 200 (incl. 2, 3); (G2) {n2} pairs p^2 q <= 60000; (G3) {n3} pairs (l, M), M <= 300, l < 60")
# H at l = 2, 3 (supersingular j = 0 = 1728 there): a few values
print("H(l, M) at l = 2, 3:", {(l, M): str(H(l, M)) for l in (2, 3) for M in (1, 2, 3, 4, 5, 7, 9, 13, 25, 49) if M % l})


# (G4) closed forms
def Hpq_closed(p, q):
    return Fr((p - 1) * (q + 1), 12) + Fr(1, 4) * (1 - kron(-4, p)) * (1 + kron(-4, q)) + Fr(1, 3) * (1 - kron(-3, p)) * (1 + kron(-3, q))


def Hqp2_closed(q, p):
    n2 = 0 if p == 2 else 1 + kron(-4, p)
    n3 = 0 if p == 3 else 1 + kron(-3, p)
    return Fr((q - 1) * p * (p + 1), 12) + Fr(1, 4) * (1 - kron(-4, q)) * n2 + Fr(1, 3) * (1 - kron(-3, q)) * n3


for p in P[:30]:
    for q in P[:30]:
        if p != q:
            if Hpq_closed(p, q) != H(p, q):
                fail(f"closed H(p,q) {p},{q}")
            if Hqp2_closed(q, p) != H(q, p * p):
                fail(f"closed H(q,p^2) {q},{p}")
            if 2 * H(q, 1) + H(p, q) != 2 * H(p, 1) + H(q, p):
                fail(f"symmetric form {p},{q}")
rows = [(55, H(5, 11), H(11, 5), genus(55)), (275, H(11, 25), genus(275)), (637, H(13, 49), genus(637))]
print("(G4) example rows: N=55: H(5,11) = %s, H(11,5) = %s, g = %d;  N=275: H(11,25) = %s, g = %d;  N=637: H(13,49) = %s, g = %d"
      % (rows[0][1], rows[0][2], rows[0][3], rows[1][1], rows[1][2], rows[2][1], rows[2][2]))
if (rows[0][1:], rows[1][1:], rows[2][1:]) != ((4, 6, 5), (26, 25), (56, 57)):
    fail("example rows")
# (G5) multiplicities
for l in (2, 3, 5, 7):
    for a in range(1, 6):
        km = [int(totient(l ** min(i, a - i))) for i in range(a + 1)]
        cw = [int(totient(gcd(l ** i, l ** (a - i)))) for i in range(a + 1)]
        if km != cw:
            fail("G5")
print("(G5) phi(l^min(i,a-i)) = phi(gcd(l^i, l^(a-i))) for l in {2,3,5,7}, a <= 5; at a = 2: multiplicities",
      {l: [int(totient(l ** min(i, 2 - i))) for i in range(3)] for l in (2, 3, 5, 7)},
      "(the middle component is reduced at l = 2)")
print("\nfailures:", bad)
print("ALL AUDIT CHECKS PASSED" if bad == 0 else "AUDIT FOUND FAILURES")
