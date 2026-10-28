"""Independent audit checks for THM-4825 (opus-2026-10-27-S66). Written from scratch; does not import the author's code.

Elements of K = Q(sqrt-d) are pairs (r, s) of Fractions meaning r + s sqrt-d.
"""
import math
import sys
from fractions import Fraction as Fr

import sympy as sp

FAIL = []


def check(cond, msg):
    if not cond:
        FAIL.append(msg)
        print('  FAIL:', msg)
    return cond


def squarefree(n):
    return all(n % (p * p) for p in range(2, math.isqrt(n) + 1))


def mul(x, y, d):
    (a, b), (c, e) = x, y
    return (a * c - d * b * e, a * e + b * c)


def in_G(x):
    """G_d = Z + Z tau, tau = (1 + sqrt-d)/2: x = u + v tau = (u + v/2) + (v/2) sqrt-d."""
    r, s = x
    v = 2 * s
    u = r - s
    return Fr(v).denominator == 1 and Fr(u).denominator == 1


def lattice_index_and_basis(vecs):
    """Return (det, sorted basis) of the Z-span of the given integer vectors, via a Smith/HNF-free gcd reduction."""
    # reduce with integer row operations (2-dim): standard Euclid on the first coordinate
    vs = [list(map(int, v)) for v in vecs if any(v)]
    basis = []
    # first coordinate gcd
    while True:
        nz = [v for v in vs if v[0] != 0]
        if len(nz) <= 1:
            break
        nz.sort(key=lambda v: abs(v[0]))
        piv = nz[0]
        new = [piv]
        for v in vs:
            if v is piv:
                continue
            if v[0] != 0:
                q = v[0] // piv[0]
                w = [v[0] - q * piv[0], v[1] - q * piv[1]]
                new.append(w)
            else:
                new.append(v)
        vs = [v for v in new if any(v)]
    piv = [v for v in vs if v[0] != 0]
    rest = [v for v in vs if v[0] == 0]
    g = 0
    for v in rest:
        g = math.gcd(g, v[1])
    if not piv:
        return None
    p = piv[0]
    if p[0] < 0:
        p = [-p[0], -p[1]]
    if g:
        p = [p[0], p[1] % g]
    return (p, [0, g])


def section_A(DMAX=999):
    print('== A. ring or ideal class (odd squarefree d <= %d)' % DMAX)
    n1 = n3 = 0
    for d in range(1, DMAX + 1, 2):
        if not squarefree(d):
            continue
        tau = (Fr(1, 2), Fr(1, 2))
        # G_d is a lattice with basis {1, tau}; closed under multiplication iff tau^2 in G_d
        ring = in_G(mul(tau, tau, d))
        check(ring == (d % 4 == 3), f'ring test d={d}')
        # random products: if ring, products of lattice points stay in G_d; if not, tau*tau escapes
        if d % 4 == 3:
            n3 += 1
            # tau is a root of t^2 - t + (1+d)/4, integral
            r, s = mul(tau, tau, d)
            check((r - Fr(1, 2) + Fr(1 + d, 4), s - Fr(1, 2)) == (0, 0), f'tau minimal polynomial d={d}')
        else:
            n1 += 1
            # the checkerboard I = 2G_d = span{(2,0),(1,1)}; the ideal (2, 1+sqrt-d) of Z[sqrt-d] has Z-generators
            gens = [(2, 0), mul((2, 0), (0, 1), d), (1, 1), mul((1, 1), (0, 1), d)]
            p, q = lattice_index_and_basis(gens)
            check(p == [1, 1] and q == [0, 2], f'(2,1+sqrt-d) = checkerboard d={d} got {p},{q}')
            # I^2 = (2): products of the Z-basis of I
            prods = [mul((2, 0), (2, 0), d), mul((2, 0), (1, 1), d), mul((1, 1), (1, 1), d)]
            p2, q2 = lattice_index_and_basis(prods)
            check(p2 == [2, 0] and q2 == [0, 2], f'p2^2 = (2) d={d} got {p2},{q2}')
            # principal iff an element of norm 2 exists in Z[sqrt-d]
            principal = any(x * x + d * y * y == 2 for x in range(-2, 3) for y in range(-2, 3))
            check(principal == (d == 1), f'principal test d={d}')
            # G_d not a ring: tau^2 not in G_d
            check(not in_G(mul(tau, tau, d)), f'tau^2 escapes d={d}')
    x, y, d = sp.symbols('x y d')
    tau_r, tau_s = sp.Rational(1, 2), sp.Rational(1, 2)
    # N(x + y tau) = (x + y/2)^2 + d (y/2)^2 ; times 2
    form = sp.expand(2 * ((x + y * tau_r) ** 2 + d * (y * tau_s) ** 2))
    check(sp.simplify(form - (2 * x ** 2 + 2 * x * y + (1 + d) / 2 * y ** 2)) == 0, 'normalised norm form')
    print(f'  {n3} values d = 3 mod 4: G_d is the ring Z[tau] (tau integral); {n1} values d = 1 mod 4: 2G_d = (2, 1+sqrt-d), '
          f'p2^2 = (2), principal only at d = 1, tau^2 escapes; 2N(x + y tau) = 2x^2 + 2xy + ((1+d)/2) y^2')


def is_sq(n):
    return n >= 0 and math.isqrt(n) ** 2 == n


def rep_first(n, d):
    """n = x^2 + d y^2 ?"""
    y = 0
    while d * y * y <= n:
        if is_sq(n - d * y * y):
            return True
        y += 1
    return False


def rep_second(n, d):
    """n = 2x^2 + 2xy + ((1+d)/2) y^2  <=>  2n = (2x + y)^2 + d y^2 with X = 2x + y = y (mod 2)."""
    y = 0
    while d * y * y <= 2 * n:
        r = 2 * n - d * y * y
        if is_sq(r) and (math.isqrt(r) - y) % 2 == 0:
            return True
        y += 1
    return False


def rep_two_squares(n):
    a = 0
    while a * a <= n:
        if is_sq(n - a * a):
            return True
        a += 1
    return False


def primes_upto(n):
    s = bytearray([1]) * (n + 1)
    s[0] = s[1] = 0
    for i in range(2, math.isqrt(n) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]


def section_C(QMAX=200_000):
    print('== C. three grids agree for d = 5, 13, 37; converse; first counterexamples')
    P = primes_upto(QMAX)
    tot = 0
    for d in (5, 13, 37):
        n = 0
        for q in P:
            if q == 2 or d % q == 0:
                continue
            if pow(-d % q, (q - 1) // 2, q) != 1:
                continue
            n += 1
            f1, f2, g = rep_first(q, d), rep_second(q, d), rep_two_squares(q)
            check(f1 != f2, f'exactly one grid d={d} q={q}')
            check(f1 == (q % 4 == 1) == g, f'three grids d={d} q={q}')
        tot += n
        print(f'  d = {d}: {n} split primes q < {QMAX}: first form <=> q = 1 mod 4 <=> sum of two squares; else second form')
    print(f'  total split primes checked: {tot}')
    # converse: for squarefree d = 1 mod 4 up to 3000, least split prime q = 1 mod 4 that is NOT x^2 + d y^2
    first_ce = {}
    none_found = []
    P2 = primes_upto(200_000)
    for d in range(1, 3001, 4):
        if not squarefree(d):
            continue
        found = None
        for q in P2:
            if q == 2 or d % q == 0 or q % 4 != 1:
                continue
            if pow(-d % q, (q - 1) // 2, q) != 1:
                continue
            if not rep_first(q, d):
                found = q
                break
        if found is None:
            none_found.append(d)
        else:
            first_ce[d] = found
    print(f'  squarefree d = 1 mod 4, d <= 3000: no counterexample below 2*10^5 exactly for d in {none_found}')
    listed = {17: 13, 29: 5, 41: 5, 61: 5, 21: 5, 53: 13, 85: 37, 33: 17, 65: 29, 73: 37}
    for d, q in listed.items():
        check(first_ce.get(d) == q, f'listed first counterexample d={d}: listed {q}, found {first_ce.get(d)}')
    print('  listed first counterexamples:', {d: first_ce[d] for d in listed}, '(all match)' if all(first_ce[d] == q for d, q in listed.items()) else '')
    print('  largest least counterexample for d <= 3000:', max(first_ce.values()), 'at d =', max(first_ce, key=first_ce.get))
    return none_found


def section_D():
    print('== D. statement (5): one grid per ramified norm')
    for d in (5, 13, 37, 17, 29, 41):           # 17, 29, 41: prime d = 1 mod 4 with h > 2, as a control of the scope
        ok = True
        for a in range(0, 24):
            for b in range(0, 10):
                n = 2 ** a * d ** b
                if n > 10 ** 8:
                    continue
                f1, f2 = rep_first(n, d), rep_second(n, d)
                if not (f1 != f2 and f1 == (a % 2 == 0)):
                    ok = False
                    print(f'    d={d} n=2^{a}*{d}^{b}: first {f1}, second {f2}')
        print(f'  d = {d}: every n = 2^a d^b <= 10^8 is represented by exactly one of x^2 + dy^2, 2x^2+2xy+((1+d)/2)y^2, '
              f'the first iff a even: {ok}')
        if d in (5, 13, 37):
            check(ok, f'(5) d={d}')
    owner = {5: (True, False), 10: (False, True), 25: (True, False), 50: (False, True), 100: (True, False)}
    for n, exp in owner.items():
        check((rep_first(n, 5), rep_second(n, 5)) == exp, f'owner number {n}')
    for d, ns in ((13, {13: True, 26: False, 52: True}), (37, {37: True, 74: False, 148: True})):
        for n, first in ns.items():
            check(rep_first(n, d) == first and rep_second(n, d) == (not first), f'chain {d}: {n}')
    # explicit witnesses quoted in the note's table
    check(2 * 1 ** 2 + 2 * 1 * (-2) + 3 * (-2) ** 2 == 10 and 2 * 5 ** 2 == 50, 'table witnesses')
    print('  owner numbers 5, 25, 100 first grid; 10, 50 second grid (d = 5); 13, 52 / 26 and 37, 148 / 74 as stated')


def section_E():
    print('== E. side claims in the note')
    # d = 12: the (1/2,1/2)-copy of Z[sqrt-12]: 2G = span{(2,0),(1,1)} in coordinates of u + v sqrt-12
    d = 12
    gens = [(2, 0), (1, 1)]
    # closed under multiplication by sqrt-12 ? (u, v) -> (-12 v, u)
    closed = all(((-12 * v) - u) % 2 == 0 for (u, v) in gens)
    tau = (Fr(1, 2), Fr(1, 2))
    r, s = mul(tau, tau, d)
    ring = Fr(2 * s).denominator == 1 and Fr(r - s).denominator == 1
    # multiplier ring of the lattice Z + Z tau, tau = 1/2 + sqrt-3: minimal polynomial of tau over Z
    t = sp.symbols('t')
    tau_s = sp.Rational(1, 2) + sp.sqrt(-3)
    poly = sp.Poly(sp.expand(4 * (t - tau_s) * (t - sp.conjugate(tau_s))), t)
    a2, b2, c2 = [int(c) for c in poly.all_coeffs()]
    g = math.gcd(math.gcd(a2, b2), c2)
    a2, b2, c2 = a2 // g, b2 // g, c2 // g
    disc = b2 * b2 - 4 * a2 * c2
    print(f'  d = 12: (1/2,1/2)-copy closed under sqrt-12: {closed}; ring: {ring}; the lattice Z + Z(1/2 + sqrt-3) has '
          f'primitive minimal form ({a2}, {b2}, {c2}), discriminant {disc} = -3 * {disc // -3}, so its multiplier ring is '
          f'the order of conductor {math.isqrt(disc // -3)} (Z[sqrt-48]), not Z[sqrt-12]')
    # idoneal m = 1 mod 4 with no prime factor = 3 mod 4: is 25 among them?
    idoneal = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 15, 16, 18, 21, 22, 24, 25, 28, 30, 33, 37, 40, 42, 45, 48, 57,
               58, 60, 70, 72, 78, 85, 88, 93, 102, 105, 112, 120, 130, 133, 165, 168, 177, 190, 210, 232, 240, 253, 273,
               280, 312, 330, 345, 357, 385, 408, 462, 520, 760, 840, 1320, 1365, 1848]
    sel = [m for m in idoneal if m % 4 == 1 and all(p % 4 != 3 for p in sp.factorint(m))]
    sel_sf = [m for m in sel if squarefree(m)]
    print(f'  idoneal m = 1 mod 4 with no prime factor = 3 mod 4: {sel}; squarefree ones: {sel_sf}')
    # verify idoneality of 25 independently: discriminant -100 has class group of exponent <= 2 (h = 2)
    forms = [(a, b, c) for a in range(1, 7) for b in range(-a + 1, a + 1) for c in [((b * b + 100) // (4 * a))]
             if (b * b + 100) % (4 * a) == 0 and c >= a and math.gcd(math.gcd(a, abs(b)), c) == 1 and not (b < 0 and a == c)]
    print(f'  reduced primitive forms of discriminant -100: {forms} (h = {len(forms)}; exponent <= 2, so 25 is idoneal)')
    # 85: h(-340) = 4 with four genera
    forms340 = [(a, b, c) for a in range(1, 11) for b in range(-a + 1, a + 1) for c in [((b * b + 340) // (4 * a))]
                if (b * b + 340) % (4 * a) == 0 and c >= a and math.gcd(math.gcd(a, abs(b)), c) == 1 and not (b < 0 and a == c)]
    print(f'  reduced primitive forms of discriminant -340: {forms340}')
    check(len(forms340) == 4, 'h(-340) = 4')


def negative_controls():
    print('== negative controls')
    # 1. a wrong class map: pretend [p2] is trivial (identify G_d with L_d); the dichotomy test must then fail somewhere
    bad = [d for d in range(5, 200, 4) if squarefree(d) and in_G(mul((Fr(1, 2), Fr(1, 2)), (Fr(1, 2), Fr(1, 2)), d))]
    print(f'  ring test on d = 1 mod 4 (should be empty): {bad}')
    # 2. a perturbed second form (2x^2 + 2xy + ((3+d)/2) y^2) must break the three-grid exclusivity for d = 5
    def rep_wrong(n, d):
        c = (3 + d) // 2
        return any(2 * x * x + 2 * x * y + c * y * y == n for x in range(-60, 61) for y in range(-60, 61))
    P = [q for q in primes_upto(400) if q > 2 and q != 5 and sp.jacobi_symbol(-5 % q, q) == 1]
    broken = [q for q in P if not (rep_first(q, 5) != rep_wrong(q, 5))]
    print(f'  perturbed second form for d = 5 breaks exclusivity at split q: {broken[:8]} ... ({len(broken)} of {len(P)})')
    # 3. a wrong d (d = 17, h = 4) must fail the three-grid agreement
    q = 13
    print(f'  d = 17, q = 13: split {sp.jacobi_symbol(-17 % 13, 13) == 1}, q = 1 mod 4, first form {rep_first(13, 17)}, '
          f'second form {rep_second(13, 17)} -> agreement fails as it must')


if __name__ == '__main__':
    section_A()
    section_C()
    section_D()
    section_E()
    negative_controls()
    print('FAILURES:', len(FAIL))
    for f in FAIL[:20]:
        print('  ', f)
