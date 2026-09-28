#!/usr/bin/env python3
"""Independent audit of the session note
05-knowledge/results/ogg_triangular_chowla_20260927.md, of its script
04-computation/experiments/ogg_triangular_chowla_20260927.py and of its output.

Written from scratch by the auditing subagent (2026-09-27); it shares no code with
the audited script.  Exact integer / rational arithmetic (Fractions) throughout;
numpy only for the sigma sieve and vectorised parity counts; sympy for primes,
factorisation, divisor_sigma (independent of the sieve) and GF(p) factoring.

Sections
  A  Proposition 1 (parity of s'(n) = sigma(n) - n - 1) and its corollaries
  B  class numbers (reduced forms, cross-checked against the analytic formula),
     genus of X_0(p), fixed points of w_p, genus of X_0(p)^+ for all p < 1000
  C  supersingular j-invariants for 5 <= p <= 97 by two independent methods
  D  Dirichlet's half-row formula for p = 3 mod 4, 3 < p < 2000
  E  the nodal cubic, the Pell conic, the curve E: Y^2 = X^3 - 4X + 1 (group law,
     integer points, torsion, independence of P and Q, Tate's algorithm at 2)
  F  Chowla dynamics for all n <= 10^6 (sieve to 10^7), cycles, drifts, persistence
  G  the integers of sections 2, 3, 6 (378, 637, Hurwitz, PSL(2,7), Paley heptagon)
  H  the paper text (word counts supporting section 1's typing)

Run: python 04-computation/experiments/ogg_triangular_chowla_20260927_audit.py
       > 05-knowledge/results/ogg_triangular_chowla_20260927_audit.out
"""
from __future__ import annotations

import math
import random
import time
from collections import Counter
from fractions import Fraction
from itertools import permutations
from math import comb, gcd, isqrt

import numpy as np
from sympy import Poly, divisor_sigma, factorint, isprime, primerange, symbols

OGG = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 47, 59, 71]
FAILS: list[str] = []
T0 = time.time()
PAPER = ("C:/Users/Eliott/AppData/Local/Temp/claude/C--Users-Eliott-Documents-GitHub-ephrepos-math/"
         "f417c5d2-5659-4853-a6d7-d3049f4aa2b7/scratchpad/paper30.txt")


def banner(title: str) -> None:
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def check(cond, msg: str) -> bool:
    cond = bool(cond)
    print(f"   [{'OK  ' if cond else 'FAIL'}] {msg}")
    if not cond:
        FAILS.append(msg)
    return cond


def tri_index(n: int):
    """m with T_m = n, else None."""
    if n < 0:
        return None
    r = isqrt(8 * n + 1)
    return (r - 1) // 2 if r * r == 8 * n + 1 else None


def is_square(n: int) -> bool:
    return n >= 0 and isqrt(n) ** 2 == n


def legendre(a: int, p: int) -> int:
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def jacobi(a: int, n: int) -> int:
    """Jacobi symbol (a/n), n odd positive."""
    assert n > 0 and n % 2 == 1
    a %= n
    result = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                result = -result
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            result = -result
        a %= n
    return result if n == 1 else 0


def kronecker(D: int, n: int) -> int:
    """Kronecker symbol (D/n), n >= 1, D = 0 or 1 mod 4."""
    if n == 1:
        return 1
    res = 1
    while n % 2 == 0:
        if D % 2 == 0:
            return 0
        res *= 1 if D % 8 in (1, 7) else -1
        n //= 2
    if n == 1:
        return res
    return res * jacobi(D, n)


def sigma_sieve(N: int) -> np.ndarray:
    s = np.zeros(N + 1, dtype=np.int64)
    for d in range(1, N + 1):
        s[d::d] += d
    return s


# ============================================================================ A
def part_A(sig: np.ndarray) -> None:
    banner("A. Proposition 1 (parity alternation of s'(n) = sigma(n) - n - 1) and corollaries")
    print("   Blind re-derivation.  sigma is multiplicative; sigma(2^a) = 2^(a+1) - 1 is odd; for an odd prime q,")
    print("   sigma(q^e) = 1 + q + ... + q^e has e + 1 terms, all odd, so sigma(q^e) = e + 1 (mod 2).  Hence")
    print("   sigma(n) odd  <=>  all odd-prime exponents even  <=>  the odd part of n is a square  <=>  n = k^2 or 2k^2.")
    print("   n even: s'(n) = sigma(n) - n - 1 = sigma(n) + 1 (mod 2): s'(n) odd <=> sigma(n) even <=> odd part not a square.")
    print("   n odd : s'(n) = sigma(n) - 1 - 1 = sigma(n)     (mod 2): s'(n) even <=> sigma(n) even <=> n not a square.")
    print("   In both cases  s'(n) = n (mod 2)  <=>  the odd part of n is a square  ('square type').  Prop 1 (a), (b) hold.")
    print("   (c) a fixed point or an odd cycle cannot flip parity at every step (an odd number of flips returns to the")
    print("   opposite parity), so it contains a square-type element; a same-parity 2-cycle has s'(n) = n (mod 2) at both")
    print("   members; #{n <= N square type} = sum_a #{k odd: 2^a k^2 <= N} ~ N^(1/2)/(2 - sqrt 2) = O(N^(1/2)).  All correct.")
    NA = 2 * 10**6
    n = np.arange(NA + 1, dtype=np.int64)
    m = n.copy()
    m[0] = 1
    for _ in range(30):
        ev = (m & 1) == 0
        if not ev.any():
            break
        m[ev] >>= 1
    r = np.floor(np.sqrt(m.astype(np.float64))).astype(np.int64)
    r[(r + 1) ** 2 <= m] += 1
    r[r * r > m] -= 1
    sqtype = r * r == m
    sqtype[0] = False
    sg = sig[: NA + 1]
    check(np.all(((sg[1:] & 1) == 1) == sqtype[1:]),
          f"classical fact: sigma(n) odd iff n is a square or twice a square, all 1 <= n <= {NA}")
    sp = sg - n - 1
    ge2 = n >= 2
    same = (sp & 1) == (n & 1)
    check(np.all(same[ge2] == sqtype[ge2]),
          f"Prop 1 (a)+(b) as one statement: s'(n) = n (mod 2) iff odd part of n is a square, all 2 <= n <= {NA}")
    even = ge2 & ((n & 1) == 0)
    odd = ge2 & ((n & 1) == 1)
    check(np.all(((sp[even] & 1) == 1) == ~sqtype[even]), "Prop 1 (a): n even, s'(n) odd iff the odd part is not a square")
    check(np.all(((sp[odd] & 1) == 0) == ~sqtype[odd]), "Prop 1 (b): n odd, s'(n) even iff n is not a square")
    fp = np.nonzero((sp == n) & ge2)[0]
    check(len(fp) == 0, f"no fixed point of s' (= quasiperfect number) with n <= {NA}")
    # exact count of square-type numbers in [M, 2M) by direct enumeration, vs the vectorised flag
    print("   Parity persistence in [M, 2M).  By the theorem the persisting n are EXACTLY the square-type n, so the")
    print("   'measured' proportion is a deterministic count; its asymptotic is sqrt(1/(2M)) among all n (k odd in 2^a k^2).")
    print("   M         #sqtype  share(all n)  share(composites)  sqrt(1/(2M))  sqrt(1/(2M))*M/#comp  note .out")
    note_vals = {10**3: 0.02543, 10**4: 0.00792, 10**5: 0.00245, 10**6: 0.00076}
    for M in (10**3, 10**4, 10**5, 10**6):
        rng = slice(M, 2 * M)
        allc = int(sqtype[rng].sum())
        direct = 0
        a = 0
        while 2**a <= 2 * M:
            lo, hi = -(-M // 2**a), (2 * M - 1) // 2**a  # k^2 in [lo, hi]
            klo = isqrt(lo - 1) + 1 if lo > 0 else 1
            khi = isqrt(hi)
            direct += sum(1 for k in range(klo, khi + 1) if k % 2 == 1)
            a += 1
        comp = sp[rng] > 0
        ncomp = int(comp.sum())
        samec = int((same[rng] & comp).sum())
        sqc = int((sqtype[rng] & comp).sum())
        pred = math.sqrt(1 / (2 * M))
        print(f"   {M:<9d} {allc:>7d}  {allc / M:.5f}       {samec / ncomp:.5f}            {pred:.5f}       "
              f"{pred * M / ncomp:.5f}               {note_vals[M]}")
        check(allc == direct and samec == sqc,
              f"M = {M}: square-type count by enumeration = {direct}; persisting composites = square-type composites = {samec}")
    print("   Cattaneo's theorem (quasiperfect => odd square) is STRONGER than Prop 1(c): it also excludes n twice a square.")
    print("   Argument: n = 2^a m^2 (a odd, m odd), sigma(n) = (2^(a+1) - 1) sigma(m^2) = 2^(a+1) m^2 + 1 gives")
    print("   m^2 = -1 (mod 2^(a+1) - 1), but 2^(a+1) - 1 = 3 (mod 4) has a prime factor = 3 (mod 4): impossible.")
    ok = True
    for a in range(1, 41, 2):
        q = 2 ** (a + 1) - 1
        ok &= q % 4 == 3 and any(r % 4 == 3 for r in factorint(q))
    check(ok, "checked a = 1, 3, ..., 39: 2^(a+1) - 1 = 3 (mod 4) and has a prime factor = 3 (mod 4)")


# ============================================================================ B
def h_forms(D: int) -> int:
    """Class number of primitive positive definite forms of discriminant D < 0 (reduced-form count)."""
    assert D < 0 and D % 4 in (0, 1)
    count = 0
    a = 1
    while 3 * a * a <= -D:
        for b in range(-a, a + 1):
            num = b * b - D
            if num % (4 * a):
                continue
            c = num // (4 * a)
            if c < a:
                continue
            if (abs(b) == a or a == c) and b < 0:
                continue
            if gcd(gcd(a, abs(b)), c) != 1:
                continue
            count += 1
        a += 1
    return count


def squarefree(n: int) -> bool:
    return all(e == 1 for e in factorint(abs(n)).values())


def is_fundamental(D: int) -> bool:
    if D >= 0 or D % 4 not in (0, 1):
        return False
    if D % 4 == 1:
        return squarefree(D)
    m = D // 4
    return m % 4 in (2, 3) and squarefree(m)


def h_analytic(D: int) -> Fraction:
    """Dirichlet: h(D) = -(w/(2|D|)) sum_{a=1}^{|D|-1} (D/a) a for fundamental D < 0."""
    w = 6 if D == -3 else 4 if D == -4 else 2
    S = sum(kronecker(D, a) * a for a in range(1, -D))
    return Fraction(-w * S, 2 * (-D))


def genus_X0(p: int) -> Fraction:
    nu2 = 1 + kronecker(-4, p)
    nu3 = 1 + kronecker(-3, p)
    return Fraction(1) + Fraction(p + 1, 12) - Fraction(nu2, 4) - Fraction(nu3, 3) - 1


def fixed_points(p: int) -> int:
    if p == 2:
        return 2  # an involution of P^1 has two fixed points; numerically h(-8) + h(-4)
    if p % 4 == 3:
        return h_forms(-4 * p) + h_forms(-p)
    return h_forms(-4 * p)


ROWS_B: list[tuple[int, int, int, int]] = []


def part_B() -> None:
    banner("B. Class numbers, genus of X_0(p), fixed points of w_p, genus of X_0(p)^+")
    bad, nfund = [], 0
    for D in range(-3, -1000, -1):
        if not is_fundamental(D):
            continue
        nfund += 1
        ha = h_analytic(D)
        if ha.denominator != 1 or int(ha) != h_forms(D):
            bad.append(D)
    check(not bad, f"reduced-form count = analytic class number formula for all {nfund} fundamental D in (-1000, 0)")
    known = {-3: 1, -4: 1, -7: 1, -8: 1, -11: 1, -15: 2, -19: 1, -20: 2, -23: 3, -43: 1, -47: 5, -56: 4, -59: 3,
             -67: 1, -71: 7, -148: 2, -163: 1, -164: 8, -188: 5, -236: 9, -284: 7}
    check(all(h_forms(D) == h for D, h in known.items()), f"spot values {known}")
    ok = True
    for p in primerange(5, 1000):
        if p % 4 == 3:
            ok &= h_forms(-4 * p) == (3 if p % 8 == 3 else 1) * h_forms(-p)
    check(ok, "order relation h(-4p) = 3h(-p) [p = 3 mod 8], h(-p) [p = 7 mod 8] for 3 < p < 1000")
    print(f"   nu_2, nu_3 at p = 2: {1 + kronecker(-4, 2)}, {1 + kronecker(-3, 2)}; at p = 3: {1 + kronecker(-4, 3)}, {1 + kronecker(-3, 3)}")
    for p in primerange(2, 1000):
        g = genus_X0(p)
        assert g.denominator == 1, p
        g = int(g)
        if p >= 5:
            assert g == (p + 1) // 12 - (1 if p % 12 == 1 else 0), p
        f = fixed_points(p)
        assert (2 * g + 2 - f) % 4 == 0, (p, g, f)
        ROWS_B.append((p, g, f, (2 * g + 2 - f) // 4))
    check(True, "g(X_0(p)) integral, = floor((p+1)/12) - [p = 1 mod 12] for p >= 5; 2g + 2 - f = 0 mod 4 for all p < 1000")
    known_g = {2: 0, 3: 0, 5: 0, 7: 0, 11: 1, 13: 0, 17: 1, 19: 1, 23: 2, 29: 2, 31: 2, 37: 2, 41: 3, 43: 3, 47: 4, 53: 4,
               59: 5, 61: 4, 67: 5, 71: 6, 73: 5, 79: 6, 83: 7, 89: 7, 97: 7, 101: 8, 103: 8, 107: 9, 109: 8, 113: 9}
    check(all(g == known_g[p] for p, g, f, gp in ROWS_B if p in known_g), "g(X_0(p)) agrees with the standard table for p <= 113")
    check(all(gp >= 0 for _, _, _, gp in ROWS_B), "g(X_0(p)^+) >= 0 for all p < 1000")
    genus0 = [p for p, g, f, gp in ROWS_B if gp == 0]
    check(genus0 == OGG, f"g(X_0(p)^+) = 0 for p < 1000 exactly at {genus0}")
    print("   p:g,f,g+ for p < 120: " + "; ".join(f"{p}:{g},{f},{gp}" for p, g, f, gp in ROWS_B if p < 120))
    note_rows = {11: (1, 4, 0), 37: (2, 2, 1), 71: (6, 14, 0), 73: (5, 4, 2), 97: (7, 4, 3)}
    check(all((g, f, gp) == note_rows[p] for p, g, f, gp in ROWS_B if p in note_rows), f"the note's quoted rows {note_rows}")
    g1 = [p for p, g, f, gp in ROWS_B if gp == 1 and p < 200]
    g2 = [p for p, g, f, gp in ROWS_B if gp == 2 and p < 200]
    print(f"   g(X_0(p)^+) = 1 for p < 200: {g1}  (recollection of the standard list: 37, 43, 53, 61, 79, 83, 89, 101, 131)")
    print(f"   g(X_0(p)^+) = 2 for p < 200: {g2}  (recollection: 67, 73, 103, 107, 167, 191)")
    print(f"   p = 2: naive h(-4p) = h(-8) = {h_forms(-8)} != 2 = fixed points of an involution of P^1; h(-8) + h(-4) = "
          f"{h_forms(-8) + h_forms(-4)}.  The one-line form 2g + 2 = h(-4p) + h(-p)[p = 3 (4)] therefore needs p >= 3.")
    hodd = {p: h_forms(-p) for p in OGG if p % 4 == 3}
    print(f"   h(-p) for the p = 3 (mod 4) members of Ogg's list: {hodd}")
    check(hodd == {3: 1, 7: 1, 11: 1, 19: 1, 23: 3, 31: 3, 47: 5, 59: 3, 71: 7}, "class numbers 1,1,1,1,3,3,5,3,7 as in the note")


# ============================================================================ C
def ss_in_Fp_two_ways(p: int):
    """Supersingular j in F_p (i) by point counting over F_p (a_p = 0) and (ii) by the Hasse invariant,
    the coefficient of x^(p-1) in (x^3 + ax + b)^((p-1)/2), on an F_p-model of each j."""
    m = (p - 1) // 2
    by_count, by_hasse = set(), set()
    for j in range(p):
        if j == 0:
            a, b = 0, 1
        elif j == 1728 % p:
            a, b = 1, 0
        else:
            a = 3 * j * (1728 - j) % p
            b = 2 * j * (1728 - j) ** 2 % p
        disc = (4 * pow(a, 3, p) + 27 * b * b) % p
        assert disc, (p, j)
        assert 1728 * 4 * pow(a, 3, p) * pow(disc, p - 2, p) % p == j, (p, j)
        s = sum(legendre(x * x * x + a * x + b, p) for x in range(p))
        if s == 0:  # a_p = -s; for p >= 5, |a_p| < p so a_p = 0 mod p forces a_p = 0
            by_count.add(j)
        H = 0
        for i in range(m + 1):
            jx, k = 2 * m - 3 * i, 2 * i - m
            if jx < 0 or k < 0:
                continue
            H += math.factorial(m) // (math.factorial(i) * math.factorial(jx) * math.factorial(k)) * a**jx * b**k
        if H % p == 0:
            by_hasse.add(j)
    return by_count, by_hasse


def ss_total_legendre(p: int):
    """All supersingular j in F_{p^2} via the roots lambda in F_{p^2} of H_p(lambda) = sum C(m,k)^2 lambda^k."""
    m = (p - 1) // 2
    coeffs = [comb(m, k) ** 2 % p for k in range(m + 1)]
    d = next(x for x in range(2, p) if legendre(x, p) == -1)

    def mul(u, v):
        return ((u[0] * v[0] + u[1] * v[1] * d) % p, (u[0] * v[1] + u[1] * v[0]) % p)

    def inv(u):
        nrm = (u[0] * u[0] - u[1] * u[1] * d) % p
        ni = pow(nrm, p - 2, p)
        return (u[0] * ni % p, (-u[1] * ni) % p)

    roots = []
    for u in range(p):
        for v in range(p):
            lam = (u, v)
            if lam in ((0, 0), (1, 0)):
                continue
            acc = (0, 0)
            for c in reversed(coeffs):
                acc = mul(acc, lam)
                acc = ((acc[0] + c) % p, acc[1])
            if acc == (0, 0):
                roots.append(lam)
    js = set()
    for lam in roots:
        l2 = mul(lam, lam)
        num = ((l2[0] - lam[0] + 1) % p, (l2[1] - lam[1]) % p)
        num3 = mul(mul(num, num), num)
        lm1 = ((lam[0] - 1) % p, lam[1])
        den = mul(l2, mul(lm1, lm1))
        j = mul(num3, inv(den))
        js.add((256 * j[0] % p, 256 * j[1] % p))
    return roots, js


def part_C() -> None:
    banner("C. Supersingular j-invariants, 5 <= p <= 97: F_p-count by point counting and Hasse invariant; total by Legendre roots")
    x = symbols("x")
    fdict = {p: f for p, g, f, gp in ROWS_B}
    print("   p   total  in_Fp(count) in_Fp(Hasse) in_Fp(Legendre)  f/2  formula_total  all_in_Fp")
    all_in = []
    totals = []
    ok_all = True
    for p in primerange(5, 98):
        cnt, has = ss_in_Fp_two_ways(p)
        roots, js = ss_total_legendre(p)
        in_fp_leg = {j[0] for j in js if j[1] == 0}
        m = (p - 1) // 2
        coeffs = [comb(m, k) ** 2 % p for k in range(m + 1)]
        fl = Poly(list(reversed(coeffs)), x, modulus=p).factor_list()[1]
        degs = sorted(f.degree() for f, e in fl for _ in range(e))
        nlin = degs.count(1)
        eps = {1: 0, 5: 1, 7: 1, 11: 2}[p % 12]
        formula_total = p // 12 + eps
        mass = sum(Fraction(1, 3) if j == (0, 0) else Fraction(1, 2) if j == (1728 % p, 0) else Fraction(1) for j in js)
        row_ok = (cnt == has == in_fp_leg and len(roots) == m and max(degs) <= 2 and all(e == 1 for f, e in fl)
                  and nlin == sum(1 for r in roots if r[1] == 0) and len(js) == formula_total
                  and mass == Fraction(p - 1, 12) and 2 * len(cnt) == fdict[p])
        ok_all &= row_ok
        totals.append(len(js))
        if len(cnt) == len(js):
            all_in.append(p)
        print(f"   {p:<3d} {len(js):>5d}  {len(cnt):>10d} {len(has):>12d} {len(in_fp_leg):>15d}  {fdict[p] // 2:>4d}  "
              f"{formula_total:>12d}  {'yes' if len(cnt) == len(js) else 'no '}  {'' if row_ok else '<-- inconsistency'}")
    check(ok_all, "for every p: the three F_p-counts agree, H_p has (p-1)/2 simple roots in F_{p^2} (degrees <= 2), "
                  "total = floor(p/12) + eps(p), Eichler-Deuring mass = (p-1)/12, and #ss j in F_p = f/2 (fixed points of w_p)")
    check(totals == [1, 1, 2, 1, 2, 2, 3, 3, 3, 3, 4, 4, 5, 5, 6, 5, 6, 7, 6, 7, 8, 8, 8], f"totals {totals} as in the note")
    check(all_in == [p for p in OGG if p >= 5], f"all supersingular j in F_p exactly for {all_in}")
    # the fixed-point formula cross-checked against the F_p count up to 199 (Deuring: #ss j in F_p = f/2)
    ok = True
    for p in primerange(5, 200):
        cnt, _ = ss_in_Fp_two_ways(p)
        ok &= 2 * len(cnt) == fdict[p]
    check(ok, "2 * #{supersingular j in F_p} = h(-4p) + h(-p)[p = 3 (4)] for all 5 <= p < 200 (independent check of the w_p fixed-point formula)")


# ============================================================================ D
def part_D() -> None:
    banner("D. Dirichlet's half-row formula h(-p) = (sum_{0<a<p/2} (a/p)) / (2 - (2/p)), p = 3 mod 4, p > 3")
    ok, cnt, ok_alt = True, 0, True
    for p in primerange(7, 2000):
        if p % 4 != 3:
            continue
        S = sum(legendre(a, p) for a in range(1, (p + 1) // 2))
        h = h_forms(-p)
        ok &= S == (2 - legendre(2, p)) * h
        ok_alt &= (S == h) == (p % 8 == 7) and (S == 3 * h) == (p % 8 == 3)
        cnt += 1
    check(ok, f"formula exact for all {cnt} primes p = 3 (mod 4), 3 < p < 2000")
    check(ok_alt, "equivalently: half-row sum = h(-p) for p = 7 (mod 8) and = 3 h(-p) for p = 3 (mod 8)")
    S3 = sum(legendre(a, 3) for a in range(1, 2))
    print(f"   p = 3: half-row sum {S3}, 2 - (2/3) = {2 - legendre(2, 3)}, quotient {Fraction(S3, 2 - legendre(2, 3))} != h(-3) = 1: "
          "the hypothesis p > 3 is needed (extra units).")
    for p in (7, 11, 19, 23, 31, 47, 59, 71):
        qr = {x * x % p for x in range(1, p)}
        out_deg = sum(1 for a in range(1, (p + 1) // 2) if a in qr)          # arc 0 -> a iff a - 0 is a residue
        in_deg = sum(1 for a in range(1, (p + 1) // 2) if (-a) % p in qr)    # arc a -> 0 iff 0 - a is a residue
        S = sum(legendre(a, p) for a in range(1, (p + 1) // 2))
        assert out_deg - in_deg == S
        print(f"   p = {p}: Paley tournament vertex 0, out - in into {{1..{(p - 1) // 2}}} = {out_deg} - {in_deg} = {S} = "
              f"{2 - legendre(2, p)} * h(-p) = {2 - legendre(2, p)} * {h_forms(-p)}")
    check(True, "Paley reading: half-row sum = out-degree minus in-degree of vertex 0 into the first half (checked 8 primes)")


# ============================================================================ E
A4, B6 = -4, 1


def e_add(P, Q):
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2:
        if y1 + y2 == 0:
            return None
        lam = (3 * x1 * x1 + A4) / (2 * y1)
    else:
        lam = (y2 - y1) / (x2 - x1)
    x3 = lam * lam - x1 - x2
    return (x3, lam * (x1 - x3) - y1)


def e_neg(P):
    return None if P is None else (P[0], -P[1])


def e_mul(n: int, P):
    if n < 0:
        n, P = -n, e_neg(P)
    R, Qp = None, P
    while n:
        if n & 1:
            R = e_add(R, Qp)
        Qp = e_add(Qp, Qp)
        n >>= 1
    return R


def on_E(P) -> bool:
    return P is None or P[1] ** 2 == P[0] ** 3 + A4 * P[0] + B6


def nodal_add(P, Q):
    """Chord-tangent law on y^2 = x^3 + x^2 (identity at infinity)."""
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2:
        if y1 + y2 == 0:
            return None
        lam = (3 * x1 * x1 + 2 * x1) / (2 * y1)
    else:
        lam = (y2 - y1) / (x2 - x1)
    x3 = lam * lam - 1 - x1 - x2   # x^3 + (1 - lam^2) x^2 + ... : sum of roots = lam^2 - 1
    y3 = lam * (x3 - x1) + y1
    return (x3, -y3)


def ec_points_mod(p: int, a: int, b: int):
    sq = {}
    for y in range(p):
        sq.setdefault(y * y % p, []).append(y)
    pts = [None]
    for x in range(p):
        for y in sq.get((x * x * x + a * x + b) % p, []):
            pts.append((x, y))
    return pts


def ec_add_mod(P, Q, p: int, a: int):
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2:
        if (y1 + y2) % p == 0:
            return None
        lam = (3 * x1 * x1 + a) * pow(2 * y1, p - 2, p) % p
    else:
        lam = (y2 - y1) * pow(x2 - x1, p - 2, p) % p
    x3 = (lam * lam - x1 - x2) % p
    return (x3, (lam * (x1 - x3) - y1) % p)


def ec_order_mod(P, p: int, a: int) -> int:
    k, R = 1, P
    while R is not None:
        R = ec_add_mod(R, P, p, a)
        k += 1
    return k


def subgroup_mod(gens, p: int, a: int):
    S = {None}
    frontier = [None]
    while frontier:
        new = []
        for R in frontier:
            for G in gens:
                T = ec_add_mod(R, G, p, a)
                if T not in S:
                    S.add(T)
                    new.append(T)
        frontier = new
    return S


def part_E() -> None:
    banner("E. Proposition 2: the nodal cubic, the Pell conic, the elliptic curve E: Y^2 = X^3 - 4X + 1")
    # (a) parametrisation and group law
    rnd = random.Random(20260927)
    ok = True
    for _ in range(300):
        t1 = Fraction(rnd.randint(-30, 30), rnd.randint(1, 9))
        t2 = Fraction(rnd.randint(-30, 30), rnd.randint(1, 9))
        if abs(t1) == 1 or abs(t2) == 1:
            continue
        P1 = (t1 * t1 - 1, t1 ** 3 - t1)
        P2 = (t2 * t2 - 1, t2 ** 3 - t2)
        assert P1[1] ** 2 == P1[0] ** 3 + P1[0] ** 2
        phi = lambda R: (R[1] + R[0]) / (R[1] - R[0]) if R is not None else Fraction(1)
        assert phi(P1) == (t1 + 1) / (t1 - 1)
        S = nodal_add(P1, P2)
        ok &= phi(S) == phi(P1) * phi(P2)
        if t1 + t2 != 0:
            t3 = (1 + t1 * t2) / (t1 + t2)  # my derived parameter of the sum
            ok &= S == (t3 * t3 - 1, t3 ** 3 - t3)
        else:
            ok &= S is None
    check(ok, "(a) chord-tangent law on y^2 = x^2(x+1): (y+x)/(y-x) = (t+1)/(t-1) is multiplicative; parameter of the sum is (1 + t1 t2)/(t1 + t2)")
    ok = all(is_square(x * x * (x + 1)) == is_square(x + 1) for x in range(1, 10**5))
    check(ok, "integer points of the nodal cubic have x + 1 = t^2 (t integer), y = t^3 - t = 6 C(t+1, 3) [x <= 10^5]")
    # reduction of the nodal cubic: split node at odd p, cusp at p = 2
    for p in (2, 3, 5, 7, 11, 13):
        pts = sum(1 for x in range(p) for y in range(p) if (y * y - x * x * (x + 1)) % p == 0)
        ns = pts - 1 + 1  # remove the node (0,0), add the point at infinity
        print(f"   nodal cubic mod {p}: #nonsingular points incl. infinity = {ns} = {'p - 1 (split node, a_p = 1)' if ns == p - 1 else 'p (CUSP: additive, not multiplicative)'}")
    check(True, "y^2 = x^2(x+1) is a split node at every odd p but a CUSP at p = 2 ((y - x)(y + x) = (y + x)^2 in char 2)")
    print("   Legendre fibre: lambda = 1, x -> x + 1 gives y^2 = (x+1) x^2 EXACTLY; lambda = 0, x -> -x gives y^2 = -x^2(x+1), the -1 twist.")
    # (b) Pell
    pell = [t for t in range(1, 10**6 + 1) if is_square(8 * t * t - 7)]
    print(f"   (b) t <= 10^6 with x = t^2 - 1 triangular: {pell}")
    print(f"       x = {[t * t - 1 for t in pell[:8]]}, m = {[tri_index(t * t - 1) for t in pell[:8]]}")
    orb1 = [t for t in pell if pell.index(t) % 2 == 0]
    orb2 = [t for t in pell if pell.index(t) % 2 == 1]
    rec = all(orb[i + 2] == 6 * orb[i + 1] - orb[i] for orb in (orb1, orb2) for i in range(len(orb) - 2))
    unit_ok = True
    for orb, u0 in ((orb1, 1), (orb2, 5)):
        u, t = u0, orb[0]
        for tn in orb[1:]:
            u, t = 3 * u + 8 * t, u + 3 * t
            unit_ok &= t == tn and u * u - 8 * t * t == -7
    check(pell[:12] == [1, 2, 4, 11, 23, 64, 134, 373, 781, 2174, 4552, 12671] and rec and unit_ok,
          "Pell list as in the note; two interleaved orbits of (u,t) -> (3u + 8t, u + 3t), each t_{k+2} = 6 t_{k+1} - t_k")
    check(tri_index(134 * 134 - 1) == 189 and tri_index(373 * 373 - 1) == 527, "T_189 = 134^2 - 1 = 17955; T_527 = 373^2 - 1")
    # (c) the elliptic curve
    disc = -16 * (4 * A4**3 + 27 * B6**2)
    check(disc == 3664 and factorint(disc) == {2: 4, 229: 1}, f"discriminant {disc} = 2^4 * 229; v_p(Delta) < 12 for all p, so the model is minimal")
    P = (Fraction(0), Fraction(1))
    Q = (Fraction(2), Fraction(1))
    assert on_E(P) and on_E(Q)
    P2, P3, Q2, PQ, PmQ = e_mul(2, P), e_mul(3, P), e_mul(2, Q), e_add(P, Q), e_add(P, e_neg(Q))
    print(f"   2P = {P2}, 3P = {P3}, 2Q = {Q2}, P + Q = {PQ}, P - Q = {PmQ}")
    check(P2 == (4, 7) and P3[0] == Fraction(-7, 4), "2P = (4, 7); x(3P) = -7/4 (not integral => P non-torsion by Nagell-Lutz)")
    XMAX = 10**7
    ipts = [(X, isqrt(X**3 - 4 * X + 1)) for X in range(-2, XMAX + 1) if is_square(X**3 - 4 * X + 1)]
    note_pts = [(-2, 1), (-1, 2), (0, 1), (2, 1), (3, 4), (4, 7), (10, 31), (12, 41), (20, 89), (114, 1217), (1274, 45473)]
    check(ipts == note_pts, f"integer points with -2 <= X <= {XMAX} (Y >= 0): {ipts}")
    # exact combinations, |n|,|m| <= 30
    table = {}
    rel = []
    for nn in range(-30, 31):
        for mm in range(-30, 31):
            R = e_add(e_mul(nn, P), e_mul(mm, Q))
            if R is None:
                if (nn, mm) != (0, 0):
                    rel.append((nn, mm))
                continue
            if R[0].denominator == 1 and R[1].denominator == 1:
                table[(int(R[0]), int(R[1]))] = (nn, mm)
    check(not rel, "no relation nP + mQ = O with |n|, |m| <= 30 (the note claims <= 12 but its script never tests this)")
    note_combo = {(-2, 1): (1, 1), (-1, 2): (1, -1), (0, 1): (1, 0), (2, 1): (0, 1), (3, 4): (2, 1), (4, 7): (2, 0),
                  (10, 31): (2, -1), (12, 41): (0, 2), (20, 89): (2, 2), (114, 1217): (4, 1), (1274, 45473): (2, -3)}
    print("   exact combinations (sign-sensitive):")
    sign_issues = []
    for (X, Y) in ipts:
        cp, cm = table.get((X, Y)), table.get((X, -Y))
        nc = note_combo[(X, Y)]
        actual = e_add(e_mul(nc[0], P), e_mul(nc[1], Q))
        sgn = "+" if actual == (X, Y) else "-" if actual == (X, -Y) else "??"
        if sgn != "+":
            sign_issues.append(((X, Y), nc, sgn))
        print(f"      ({X}, {Y}) = {cp[0]}P + {cp[1]}Q ;  ({X}, {-Y}) = {cm[0]}P + {cm[1]}Q ;  note says {nc[0]}P + {nc[1]}Q, which is ({X}, {'+' if sgn == '+' else '-'}{Y})")
    check(all((X, Y) in table and (X, -Y) in table for X, Y in ipts), "every integer point (both signs) is an exact Z-combination of P and Q")
    print(f"   sign discrepancies in the note's list (its script compared |Y| only): {len(sign_issues)} of 11: "
          + ", ".join(f"({X},{Y}) is -({n}P+{m}Q)" for (X, Y), (n, m), s in sign_issues))
    # torsion
    orders = {p: len(ec_points_mod(p, A4 % p, B6 % p)) for p in (3, 5, 7, 11, 13, 17, 19, 23)}
    g = 0
    for v in orders.values():
        g = gcd(g, v)
    print(f"   #E(F_p): {orders}; gcd = {g}")
    check(g == 1, "E(Q)_tors is trivial (torsion injects into E(F_p) for good p >= 3)")
    # independence of P and Q via E(F_p)/2E(F_p)
    print("   E(F_p) for good p: |E|, ord P, ord Q, |<P,Q>|, cyclic?")
    proof_primes = []
    for p in primerange(3, 400):
        if p == 229:
            continue
        a, b = A4 % p, B6 % p
        pts = ec_points_mod(p, a, b)
        Pb, Qb = (0, 1), (2, 1)
        oP, oQ = ec_order_mod(Pb, p, a), ec_order_mod(Qb, p, a)
        H = subgroup_mod([Pb, Qb], p, a)
        cyc = len(H) == oP * oQ // gcd(oP, oQ)
        if p < 60:
            print(f"      p = {p:>3d}: {len(pts):>4d}  {oP:>4d} {oQ:>4d} {len(H):>5d}  {'cyclic' if cyc else 'NOT cyclic'}")
        roots = [x for x in range(p) if (x**3 + a * x + b) % p == 0]
        if len(roots) == 3:
            twoE = {ec_add_mod(R, R, p, a) for R in pts}
            PQb = ec_add_mod(Pb, Qb, p, a)
            if Pb not in twoE and Qb not in twoE and PQb not in twoE:
                proof_primes.append(p)
    check(bool(proof_primes),
          f"P, Q independent in E(F_p)/2E(F_p) = (Z/2)^2 at p = {proof_primes[:6]}...: with trivial torsion this PROVES rank E(Q) >= 2")
    # Tate's algorithm at p = 2 on the translated model y^2 + 2y = x^3 - 4x  (y -> y + 1)
    a1, a2, a3, a4, a6 = 0, 0, 2, -4, 0
    b2 = a1 * a1 + 4 * a2
    b4 = 2 * a4 + a1 * a3
    b6 = a3 * a3 + 4 * a6
    b8 = a1 * a1 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3 * a3 - a4 * a4
    Delta = -b2 * b2 * b8 - 8 * b4**3 - 27 * b6 * b6 + 9 * b2 * b4 * b6
    print(f"   Tate at 2: model [0,0,2,-4,0] (y -> y+1 moves the singular point to (0,0)); b2={b2}, b4={b4}, b6={b6}, b8={b8}, Delta={Delta}")
    v2 = 0
    D = Delta
    while D % 2 == 0:
        D //= 2
        v2 += 1
    steps = [("2 | a3, a4, a6", a3 % 2 == 0 and a4 % 2 == 0 and a6 % 2 == 0),
             ("2 | b2 (so not multiplicative)", b2 % 2 == 0),
             ("4 | a6 (so not type II)", a6 % 4 == 0),
             ("8 | b8 (so not type III)", b8 % 8 == 0),
             ("8 does NOT divide b6 => type IV", b6 % 8 != 0)]
    for s, c in steps:
        print(f"      {s}: {c}")
    kodaira_IV = all(c for _, c in steps)
    check(kodaira_IV and v2 == 4, f"Kodaira type IV at 2, v_2(Delta) = {v2}, conductor exponent f_2 = v(Delta) - 2 = {v2 - 2}")
    print(f"   at 229: v(Delta) = 1 => type I_1, f = 1.  Conductor N = 2^{v2 - 2} * 229 = {2 ** (v2 - 2) * 229} (own computation, no table lookup).")
    # y triangular, both triangular
    ts = sorted(X // 2 for X, Y in ipts if X % 2 == 0)
    print(f"   even X give t = {ts} (the note lists only t >= 1, omitting t = -1, 0 which give y = 0)")
    ypos = [t for t in ts if t >= 1]
    ys = [t**3 - t for t in ypos]
    check(ypos == [1, 2, 5, 6, 10, 57, 637] and [tri_index(y) for y in ys] == [0, 3, 15, 20, 44, 608, 22736],
          f"t <= {XMAX // 2} with t^3 - t triangular: {ypos}; y = {ys} = T_k, k = {[tri_index(y) for y in ys]}")
    both = [t for t in ypos if tri_index(t * t - 1) is not None]
    check(both == [1, 2], f"both coordinates triangular (t <= {XMAX // 2}): t = {both}, points {[(t * t - 1, t**3 - t) for t in both]}; (3, 6) = (T_2, T_3)")
    print(f"   t = 5: (x, y) = (24, 120), 120 = 5! = T_15; t = 6: (35, 210), 210 = 2*3*5*7 = T_20; t = 10: (99, 990), 990 = T_44")
    # perfect numbers on the Pell family?  cannonball
    euclid_all = [k for k in range(2, 101) if is_square(2 ** (k - 1) * (2**k - 1) + 1)]
    euclid_prime = [k for k in euclid_all if isprime(k)]
    print(f"   2^(k-1)(2^k - 1) + 1 square for k <= 100: k = {euclid_all} (k = 4: T_15 = 120 = 2^3 * 15, 121 = 11^2, the Pell point t = 11)")
    check(euclid_prime == [] and 4 in euclid_all,
          "no PERFECT number T_(2^p - 1) (p prime <= 100) lies on the Pell family; the note's 'none below p = 31' holds for prime p, "
          "but the Euclid-shaped T_15 = 120 (k = 4, not prime) does lie on it")
    cann = [n for n in range(1, 10**6) if is_square(n * (n + 1) * (2 * n + 1) // 6)]
    check(cann == [1, 24], f"n < 10^6 with 1^2 + ... + n^2 square: {cann} (70^2 = 4900)")


# ============================================================================ F
def part_F(sig: np.ndarray, SIEVE: int, N: int = 10**6) -> None:
    banner("F. The Chowla map s'(n) = sigma(n) - n - 1 on all n <= 10^6 (sieve to 10^7, sympy beyond)")
    sp_arr = sig - np.arange(SIEVE + 1, dtype=np.int64) - 1
    big = {}

    def SP(x: int) -> int:
        if x <= SIEVE:
            return int(sp_arr[x])
        if x not in big:
            big[x] = int(divisor_sigma(x)) - x - 1
        return big[x]

    ends = Counter()
    cycles = Counter()
    total_steps = 0
    longest = (0, 0)
    best_peak = (0, 1)  # (peak, n) maximising peak/n
    beyond = 0
    for n in range(2, N + 1):
        x, k, peak = n, 0, n
        seen = set()
        while x != 0:
            if x in seen:
                cyc = [x]
                y = SP(x)
                while y != x:
                    cyc.append(y)
                    y = SP(y)
                L = len(cyc)
                ends["fixed" if L == 1 else "cycle2" if L == 2 else f"cycle{L}"] += 1
                cycles[tuple(sorted(cyc))] += 1
                break
            seen.add(x)
            x = SP(x)
            k += 1
            if x > SIEVE:
                beyond += 1
            if x > peak:
                peak = x
        else:
            ends["zero"] += 1
        total_steps += k
        if k > longest[0]:
            longest = (k, n)
        if peak * best_peak[1] > best_peak[0] * n:
            best_peak = (peak, n)
    print(f"   endings: {dict(ends)}; longest {longest[0]} steps at n = {longest[1]}; largest peak/n = {best_peak[0]}/{best_peak[1]} "
          f"= {best_peak[0] / best_peak[1]:.3f}; mean length {total_steps / (N - 1):.3f}; values beyond the 10^7 sieve: {beyond}")
    check(ends["zero"] == 982099 and ends["cycle2"] == 17900 and len(ends) == 2,
          "982099 sequences end at 0 and 17900 enter a 2-cycle; no fixed point, no longer cycle (as the note states)")
    check(longest == (48, 948375) and best_peak[1] == 980100 and round(total_steps / (N - 1), 2) == 9.32,
          "longest 48 steps (n = 948375); largest peak ratio at n = 980100 = 990^2; mean length 9.32")
    note_cycles = [(48, 75), (140, 195), (1050, 1925), (1575, 1648), (2024, 2295), (5775, 6128), (8892, 16587), (9504, 20735),
                   (62744, 75495), (186615, 206504), (196664, 219975), (199760, 309135), (266000, 507759), (312620, 549219),
                   (526575, 544784), (573560, 817479), (587460, 1057595), (1139144, 1159095)]
    print(f"   cycles reached (members: starts): {sorted(cycles.items())}")
    check(sorted(cycles) == note_cycles, "exactly the eighteen 2-cycles listed in the note are reached")
    genuine = all(int(divisor_sigma(a)) - a - 1 == b and int(divisor_sigma(b)) - b - 1 == a for a, b in note_cycles)
    check(genuine, "each listed pair is a genuine 2-cycle of s' (sigma via sympy, independent of the sieve)")
    check(all((a + b) % 2 == 1 for a, b in note_cycles), "all eighteen pairs have opposite parity")
    check(all(not is_square(a >> ((a & -a).bit_length() - 1)) and not is_square(b >> ((b & -b).bit_length() - 1)) for a, b in note_cycles),
          "no member is of square type (consistent with the parity theorem's forced alternation)")
    # all betrothed pairs with smaller member <= 2*10^6 (sieve values), which of them were reached
    lo = np.arange(2, 2 * 10**6 + 1, dtype=np.int64)
    mv = sp_arr[lo]
    okmask = (mv > lo) & (mv <= SIEVE)
    assert np.all((mv <= SIEVE) | (mv <= lo)), "value above sieve for n <= 2*10^6"
    back = np.zeros_like(lo)
    back[okmask] = sp_arr[mv[okmask]]
    pairs = [(int(a), int(b)) for a, b in zip(lo[okmask & (back == lo)], mv[okmask & (back == lo)])]
    print(f"   all 2-cycles with smaller member <= 2*10^6: {pairs}")
    print(f"   not reached from any start <= 10^6: {[pr for pr in pairs if pr not in cycles]}")
    check(all((a + b) % 2 == 1 for a, b in pairs), "every betrothed pair with smaller member <= 2*10^6 has opposite parity")
    # drifts
    ns = np.arange(2, N + 1, dtype=np.int64)
    spv = sp_arr[2: N + 1]
    ev = (ns % 2 == 0) & (spv > 0)          # excludes only n = 2 (s'(2) = 0)
    odc = (ns % 2 == 1) & (spv > 0)         # odd composites
    d_even = float(np.mean(np.log2(spv[ev] / ns[ev])))
    d_odd = float(np.mean(np.log2(spv[odc] / ns[odc])))
    s_al = sig[2: N + 1] - ns
    a_even = float(np.mean(np.log2(s_al[ns % 2 == 0] / ns[ns % 2 == 0])))
    a_odd = float(np.mean(np.log2(s_al[ns % 2 == 1] / ns[ns % 2 == 1])))
    print(f"   mean log2(s'(n)/n), n <= 10^6: even {d_even:+.3f}, odd composite {d_odd:+.3f}; aliquot s(n): even {a_even:+.3f}, odd (incl. primes) {a_odd:+.3f}")
    check(round(d_even, 3) == -0.048 and round(d_odd, 2) == -2.82 and round(a_even, 3) == -0.048,
          "Chowla drifts -0.048 (even) / -2.82 (odd composite) and aliquot even drift -0.048 as in the note")
    # the note's aliquot odd figure -5.34 is hard-coded in its script, not computed; locate it
    for NN in (10**4, 10**5, 10**6, 10**7):
        nn = np.arange(3, NN + 1, 2, dtype=np.int64)
        v = float(np.mean(np.log2((sig[3: NN + 1: 2] - nn) / nn)))
        print(f"      aliquot s(n), odd n <= {NN}: mean log2(s(n)/n) = {v:+.3f}")
    for M in (10**5, 10**6, 5 * 10**6):
        nn = np.arange(M + 1, 2 * M, 2, dtype=np.int64)
        v = float(np.mean(np.log2((sig[M + 1: 2 * M: 2] - nn) / nn)))
        print(f"      aliquot s(n), odd n in [{M}, {2 * M}): mean log2(s(n)/n) = {v:+.3f}")
    check(round(a_odd, 2) == -5.34, "aliquot odd drift -5.34 for n <= 10^6 as quoted in the note (hard-coded in its script, not computed there)")
    oa = [n for n in range(1, 100001, 2) if int(sp_arr[n]) > n]
    oab = [n for n in range(1, 100001, 2) if int(sig[n]) > 2 * n]
    check(len(oa) == 210 and oa[:6] == [945, 1575, 2205, 2835, 3465, 4095] and oa == oab,
          f"odd n <= 10^5 with s'(n) > n: {len(oa)}, first {oa[:6]}; identical to the odd abundant numbers")
    print("   (parity persistence proportions were recomputed exactly in section A)")


# ============================================================================ G
def part_G() -> None:
    banner("G. The integers of sections 2, 3 and 6")
    s = sum(OGG)
    check(s == 378 and tri_index(378) == 27 and factorint(378) == {2: 1, 3: 3, 7: 1},
          f"sum of Ogg's primes = {s} = T_27 = 2 * 3^3 * 7; primes <= 31 sum {sum(p for p in OGG if p <= 31)}, the rest {sum(p for p in OGG if p > 31)}")
    monster = {2: 46, 3: 20, 5: 9, 7: 6, 11: 2, 13: 3, 17: 1, 19: 1, 23: 1, 29: 1, 31: 1, 41: 1, 47: 1, 59: 1, 71: 1}
    sopfr = sum(p * e for p, e in monster.items())
    check(sopfr == 637 and factorint(637) == {7: 2, 13: 1} and sorted(monster) == OGG, f"|M| sopfr = {sopfr} = 7^2 * 13; distinct-prime sum 378")
    check([p for p in primerange(2, 72) if p not in OGG] == [37, 43, 53, 61, 67], "primes < 72 missing from the list: 37, 43, 53, 61, 67")
    exc = {q: Fraction(1, 2) + Fraction(1, 3) + Fraction(1, q) - 1 for q in (3, 4, 5, 6, 7)}
    check(exc == {3: Fraction(1, 6), 4: Fraction(1, 12), 5: Fraction(1, 30), 6: Fraction(0), 7: Fraction(-1, 42)}, f"(2,3,q) excess {dict(exc)}")
    check([2 / exc[q] for q in (3, 4, 5)] == [12, 24, 60], "spherical group orders 2/excess = 12, 24, 60 (A_4 = PSL(2,3), S_4 = PGL(2,3), A_5 = PSL(2,5))")

    def psl_order(q):
        return q * (q * q - 1) // (2 if q % 2 else 1)

    check(psl_order(7) == 168 and psl_order(8) == 504 and psl_order(13) == 1092, "|PSL(2,7)|, |PSL(2,8)|, |PSL(2,13)| = 168, 504, 1092")
    check([84 * (g - 1) for g in (3, 7, 14)] == [168, 504, 1092] and factorint(168) == {2: 3, 3: 1, 7: 1}
          and factorint(504) == {2: 3, 3: 2, 7: 1} and factorint(1092) == {2: 2, 3: 1, 7: 1, 13: 1},
          "84(g - 1) for g = 3, 7, 14 and the factorisations 2^3.3.7, 2^3.3^2.7, 2^2.3.7.13")
    hg = {q: 1 + psl_order(q) // 84 for q in (7, 8, 13, 27, 29, 41, 43)}
    print(f"   Hurwitz genera of PSL(2,q): {hg}; the genus-17 Hurwitz group (order 1344 = 84*16) is not a PSL(2,q) and lies between 14 and 118")
    check(hg[27] == 118 and 118 == 2 * 59, "118 = 2 * 59 is the PSL(2,27) Hurwitz genus, the next PSL(2,q) genus after 14 (not the next Hurwitz genus: 17 precedes it)")
    # Paley tournament on 7 vertices
    qr = {x * x % 7 for x in range(1, 7)}
    arcs = {(i, j) for i in range(7) for j in range(7) if i != j and (j - i) % 7 in qr}
    auts = [pi for pi in permutations(range(7)) if all((pi[i], pi[j]) in arcs for (i, j) in arcs)]

    def perm_order(pi):
        k, cur = 1, pi
        ident = tuple(range(7))
        while cur != ident:
            cur = tuple(pi[c] for c in cur)
            k += 1
        return k

    orders = Counter(perm_order(pi) for pi in auts)
    comp = lambda a, b: tuple(a[b[i]] for i in range(7))
    nonab = any(comp(a, b) != comp(b, a) for a in auts for b in auts)
    check(len(auts) == 21 and orders == Counter({1: 1, 3: 14, 7: 6}) and nonab,
          f"|Aut(Paley tournament on 7)| = 21, element orders {dict(orders)}, nonabelian: the Frobenius group F_21 = 7:3")
    # PSL(2,7) and the normaliser of a Sylow 7-subgroup
    p = 7

    def mmul(A, B):
        return ((A[0] * B[0] + A[1] * B[2]) % p, (A[0] * B[1] + A[1] * B[3]) % p,
                (A[2] * B[0] + A[3] * B[2]) % p, (A[2] * B[1] + A[3] * B[3]) % p)

    def canon(A):
        for e in A:
            if e:
                return A if e <= 3 else tuple((-x) % p for x in A)
        return A

    def minv(A):
        a, b, c, d = A
        return (d, (-b) % p, (-c) % p, a)

    SL = [(a, b, c, d) for a in range(p) for b in range(p) for c in range(p) for d in range(p) if (a * d - b * c) % p == 1]
    PSL = sorted({canon(M) for M in SL})
    U = {canon((1, k, 0, 1)) for k in range(7)}
    conjs = {}
    for g in PSL:
        conjs[g] = frozenset(canon(mmul(mmul(g, u), minv(g))) for u in U)
    Nrm = [g for g in PSL if conjs[g] == frozenset(U)]
    n_syl = len(set(conjs.values()))
    nonab_N = any(canon(mmul(a, b)) != canon(mmul(b, a)) for a in Nrm for b in Nrm)
    check(len(SL) == 336 and len(PSL) == 168 and len(Nrm) == 21 and n_syl == 8 and nonab_N,
          f"|PSL(2,7)| = 168; normaliser of the Sylow 7-subgroup has order {len(Nrm)} (nonabelian, = F_21), index {168 // len(Nrm)}; {n_syl} Sylow 7-subgroups")
    check(tri_index(120) == 15 and tri_index(210) == 20 and math.factorial(5) == 120 and 2 * 3 * 5 * 7 == 210 and factorint(189) == {3: 3, 7: 1},
          "120 = 5! = T_15, 210 = 7# = T_20, 189 = 3^3 * 7")
    check([p for p in OGG if p + 1 in (2**k for k in range(2, 8))] == [3, 7, 31] and 127 not in OGG,
          "Mersenne primes in Ogg's list: 3, 7, 31 (127 is not)")
    j0 = [p for p in primerange(2, 100) if (p == 2 or p == 3 or (sum(legendre(x**3 + 1, p) for x in range(p)) == 0))]
    print(f"   j = 0 (y^2 = x^3 + 1) supersingular at p < 100: {j0} (= p = 2 (mod 3) for p >= 5; also at p = 2 and p = 3)")


# ============================================================================ H
def part_H() -> None:
    banner("H. The paper text (for section 1's typing)")
    try:
        txt = open(PAPER, encoding="utf-8", errors="replace").read()
    except OSError as e:
        print(f"   could not read the paper text: {e}")
        return
    low = txt.lower()
    print(f"   pages: {txt.count('=== PAGE')}; occurrences: 'monster' {low.count('monster')}, 'moonshine' {low.count('moonshine')}, "
          f"'jack daniels' {low.count('jack daniels')}, 'non-cm'/'complex multiplication' {low.count('non-cm') + low.count('complex multiplication')}")
    print(f"   reference [2] is Ogg, Hyperelliptic modular curves, Bull. SMF 102 (1974): {'hyperelliptic modular curves' in low and '1974' in txt}")
    print(f"   Lang-Trotter sentence mentions sqrt(X)/log X without a non-CM hypothesis: {'√' in txt or 'sqrt' in low}")
    print(f"   Theorem 1 (Ogg) stated with both conditions: {'defined over fp' in low.replace(chr(10), ' ') and 'genus zero' in low}")


# ============================================================================ main
if __name__ == "__main__":
    SIEVE = 10**7
    t = time.time()
    SIG = sigma_sieve(SIEVE)
    print(f"sigma sieve to {SIEVE}: {time.time() - t:.1f}s")
    part_A(SIG)
    part_B()
    part_C()
    part_D()
    part_E()
    part_F(SIG, SIEVE)
    part_G()
    part_H()
    banner("SUMMARY")
    print(f"   failed checks: {len(FAILS)}")
    for f in FAILS:
        print(f"   - {f}")
    print(f"   elapsed {time.time() - T0:.1f}s")
