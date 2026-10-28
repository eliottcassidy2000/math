"""Independent audit checks for THM-4826 (opus-2026-10-27-S66). Written from scratch; does not import the author's code.

Eisenstein integers are pairs (a, b) = a + b*w, w^2 + w + 1 = 0; product (a,b)(c,d) = (ac - bd, ad + bc - bd).
"""
import math
import random
import sys
from fractions import Fraction as Fr
from itertools import product

import numpy as np
import sympy as sp

FAIL = []


def check(cond, msg):
    if not cond:
        FAIL.append(msg)
        print('  FAIL:', msg)
    return cond


def emul(x, y):
    a, b = x
    c, d = y
    return (a * c - b * d, a * d + b * c - b * d)


def enorm(x):
    a, b = x
    return a * a - a * b + b * b


def primes_upto(n):
    s = bytearray([1]) * (n + 1)
    s[0] = s[1] = 0
    for i in range(2, math.isqrt(n) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]


def is_sq(n):
    return n >= 0 and math.isqrt(n) ** 2 == n


# ------------------------------------------------------------------------------------------------ (1)
def part1():
    print('== (1) the 1:2 split: Z[w] (n <= 20000) and the Hurwitz order (n <= 400)')
    N = 20000
    R = 2 * math.isqrt(N // 3 + 1) + 2
    first = np.zeros(N + 1, dtype=np.int64)
    cent = np.zeros(N + 1, dtype=np.int64)
    a = np.arange(-R, R + 1)
    A, B = np.meshgrid(a, a, indexing='ij')
    Nm = A * A - A * B + B * B
    mask = Nm <= N
    nf = Nm[mask & (B % 2 == 0)]
    nc = Nm[mask & (B % 2 != 0)]
    np.add.at(first, nf, 1)
    np.add.at(cent, nc, 1)
    # r(n) = 6 sum_{e | n} chi_{-3}(e), by a divisor sieve
    chi = np.array([0, 1, -1])
    S = np.zeros(N + 1, dtype=np.int64)
    for e in range(1, N + 1):
        c = chi[e % 3]
        if c:
            S[e::e] += c
    tot = first + cent
    ok_tot = bool(np.all(tot[1:] == 6 * S[1:]))
    odd = np.arange(1, N + 1, 2)
    even = np.arange(2, N + 1, 2)
    ok_odd = bool(np.all(cent[odd] == 2 * first[odd])) and bool(np.all(first[odd] == 2 * S[odd]))
    ok_even = bool(np.all(cent[even] == 0))
    check(ok_tot and ok_odd and ok_even, 'Z[w] shells')
    print(f'  Z[w]: r(n) = 6 sum chi_-3 for all n <= {N}: {ok_tot}; odd n: first = 2 sum chi_-3 = r/3, centred = 2 first: '
          f'{ok_odd}; even n: centred = 0: {ok_even}')
    # Hurwitz: Lipschitz r4(n) and half-integer counts via theta-series convolution
    M = 400
    r1 = np.zeros(M + 1, dtype=np.int64)
    for x in range(-math.isqrt(M), math.isqrt(M) + 1):
        r1[x * x] += 1
    r2 = np.convolve(r1, r1)[:M + 1]
    r4 = np.convolve(r2, r2)[:M + 1]
    # half-integer coordinates: X = 2x odd, sum X^2 = 4n
    M4 = 4 * M
    o1 = np.zeros(M4 + 1, dtype=np.int64)
    for X in range(-math.isqrt(M4), math.isqrt(M4) + 1):
        if X % 2:
            o1[X * X] += 1
    o2 = np.convolve(o1, o1)[:M4 + 1]
    o4 = np.convolve(o2, o2)[:M4 + 1]
    half = np.array([o4[4 * n] for n in range(M + 1)])
    sig = [0] + [sum(d for d in range(1, n + 1) if n % d == 0) for n in range(1, M + 1)]
    ok = True
    for n in range(1, M + 1):
        m = n
        while m % 2 == 0:
            m //= 2
        if n % 2:
            ok &= (r4[n] == 8 * sig[n]) and (half[n] == 16 * sig[n])
        else:
            ok &= (r4[n] == 24 * sig[m]) and (half[n] == 0)
    check(ok, 'Hurwitz shells')
    print(f'  Hurwitz, n <= {M}: odd n: Lipschitz 8 sigma(n), centred 16 sigma(n); even n: Lipschitz 24 sigma(odd part), '
          f'centred 0: {ok}')
    # the mechanism: unit of order 3 is not = 1 mod P and cycles F4*
    w = (0, 1)
    check(enorm((w[0] - 1, w[1])) == 3, 'N(w - 1) = 3')
    zeta = (Fr(-1, 2), Fr(1, 2), Fr(1, 2), Fr(1, 2))
    d1 = tuple(z - (1 if i == 0 else 0) for i, z in enumerate(zeta))
    check(sum(c * c for c in d1) == 3, 'N(zeta - 1) = 3')
    # (1+i)H = D4: index 4 in H and contains only even-norm Lipschitz vectors
    print('  units of order 3: N(w - 1) = 3 and N(zeta - 1) = 3 (odd), so both act on F4* by a 3-cycle')


# ------------------------------------------------------------------------------------------------ (2)
def elements_norm(d, n):
    """Elements x + y tau of Z[(1+sqrt-d)/2] of norm n: (2x + y)^2 + d y^2 = 4n.  Returns list of (x, y)."""
    out = []
    y = 0
    while d * y * y <= 4 * n:
        r = 4 * n - d * y * y
        if is_sq(r):
            X = math.isqrt(r)
            for XX in {X, -X}:
                for yy in {y, -y}:
                    if (XX - yy) % 2 == 0:
                        out.append(((XX - yy) // 2, yy))
        y += 1
    return out


def part2():
    print('== (2) the family Z[(1+sqrt-d)/2], d = 3 mod 4')
    sev = [d for d in range(7, 200, 8) if all(d % (p * p) for p in range(2, 15))]
    for d in sev[:12]:
        bad = 0
        for n in range(1, 3001, 2):
            for (x, y) in elements_norm(d, n):
                if y % 2:
                    bad += 1
        check(bad == 0, f'd={d} = 7 mod 8 odd shells')
    print(f'  d = 7 mod 8, d in {sev[:12]}: no element of odd norm n <= 3000 on the centred grid')
    P = primes_upto(300_000)
    res = {}
    for d in (11, 19, 43, 67, 163, 35, 51, 91, 115, 123, 187, 235, 267, 403, 427):
        lim = 300_000 if d in (11, 19, 43, 67, 163) else 100_000
        nsplit = nonempty = nfirst = 0
        for p in P:
            if p > lim:
                break
            if p == 2 or d % p == 0 or pow(-d % p, (p - 1) // 2, p) != 1:
                continue
            nsplit += 1
            els = elements_norm(d, p)
            if not els:
                continue
            nonempty += 1
            grids = {y % 2 == 0 for (x, y) in els}
            check(len(grids) == 1 and len(els) == 4, f'one grid d={d} p={p}')
            f = grids.pop()
            # first iff p = x^2 + d y^2
            fy = any(is_sq(p - d * yy * yy) for yy in range(0, math.isqrt(p // d) + 1))
            check(f == fy, f'first iff x^2+dy^2 d={d} p={p}')
            nfirst += f
        res[d] = (nsplit, nonempty, nfirst)
        print(f'  d = {d:3d}: split primes < {lim}: {nsplit}; non-empty prime shells {nonempty}; all four elements on one grid; '
              f'first grid for {nfirst} = {nfirst / nonempty:.4f} of the non-empty shells ({nfirst / nsplit:.4f} of all split)')
    return res


# ------------------------------------------------------------------------------------------------ (3)
def eround(z):
    """Nearest Eisenstein integer to a complex number given exactly as a pair of Fractions (u, v) = u + v w."""
    u, v = z
    best = None
    for a in (math.floor(u), math.floor(u) + 1):
        for b in (math.floor(v), math.floor(v) + 1):
            dd = (u - a, v - b)
            nn = dd[0] * dd[0] - dd[0] * dd[1] + dd[1] * dd[1]
            if best is None or nn < best[0]:
                best = (nn, (a, b))
    return best[1]


def ediv(x, y):
    """x / y as exact Fractions in the basis (1, w): x * conj(y) / N(y); conj(c + d w) = (c - d) - d w."""
    c, d = y
    cy = (c - d, -d)
    num = emul(x, cy)
    n = enorm(y)
    return (Fr(num[0], n), Fr(num[1], n))


def egcd(x, y):
    while y != (0, 0):
        q = eround(ediv(x, y))
        qy = emul(q, y)
        r = (x[0] - qy[0], x[1] - qy[1])
        x, y = y, r
    return x


UNITS = [(1, 0), (0, 1), (-1, -1), (-1, 0), (0, -1), (1, 1)]


def primary_prime(p):
    """A primary prime a + b w (a = 2 mod 3, b = 0 mod 3) of norm p = 1 mod 3, by a Euclidean gcd with r - w."""
    while True:
        x = random.randrange(2, p - 1)
        r = pow(x, (p - 1) // 3, p)
        if r != 1:
            break
    assert (r * r + r + 1) % p == 0
    g = egcd((p, 0), (r, -1))
    assert enorm(g) == p, (p, g)
    prim = [emul(u, g) for u in UNITS]
    prim = [z for z in prim if z[0] % 3 == 2 and z[1] % 3 == 0]
    assert len(prim) == 1
    return prim[0]


def part3():
    print('== (3) the cubic character of 2 (all p = 1 mod 3 below 10^6)')
    random.seed(1)
    P = [p for p in primes_upto(1_000_000) if p % 3 == 1]
    mism = nres = 0
    sym_ok = True
    conj_ok = True
    for p in P:
        a, b = primary_prime(p)
        grid = b % 2 == 0
        cub = pow(2, (p - 1) // 3, p) == 1
        f27 = any(is_sq(p - 27 * y * y) for y in range(0, math.isqrt(p // 27) + 1))
        if not (grid == cub == f27):
            mism += 1
        nres += cub
        # supplementary law with values: (2/pi)_3 = w^k with 2^((p-1)/3) = r^k mod p, r = image of w mod pi = -a/b
        r = (-a * pow(b, -1, p)) % p
        s = pow(2, (p - 1) // 3, p)
        k = [0, 1, 2][[1, r, r * r % p].index(s)]
        cls = {(1, 0): 0, (0, 1): 1, (1, 1): 2}[(a % 2, b % 2)]       # pi mod 2 = 1, w, w^2 = 1 + w
        sym_ok &= (k == cls)
        # the conjugate primary prime: (a - b) - b w has the same grid
        ca, cb = a - b, -b
        conj_ok &= (ca % 3 == 2 and cb % 3 == 0 and (cb % 2 == 0) == grid)
    check(mism == 0 and sym_ok and conj_ok, 'cubic character')
    print(f'  {len(P)} primes; mismatches of (b even) <=> (2 cubic residue) <=> (p = x^2 + 27y^2): {mism}; 2 a cubic residue '
          f'for {nres} ({nres / len(P):.4f}); the symbol (2/pi)_3 equals pi mod 2 in {{1, w, w^2}} for every p: {sym_ok}; '
          f'the conjugate primary prime has the same grid: {conj_ok}')
    small = [(p, primary_prime(p)) for p in P if p < 50]
    print('  p < 50: ' + ', '.join(f'{p}: pi = {a}{b:+d}w ({"first" if b % 2 == 0 else "centred"})' for p, (a, b) in small))
    # Mersenne primes
    listed = {5: (2, 1), 7: (10, 1), 13: (46, 15), 17: (362, 1), 19: (298, 127), 31: (46162, 783),
              61: (1505304098, 38443119), 89: (17376907720394, 3426547218625),
              107: (9286834445316902, 1677895317602657), 127: (9328321181472828398, 1754608824145338993)}
    ok = True
    for q, (x, y) in listed.items():
        p = 2 ** q - 1
        ok &= (x * x + 27 * y * y == p)
        m = pow(3, -1, q)
        ok &= pow(2, 3 * m, p) == 2                     # 2 is a cube
        ok &= pow(2, (p - 1) // 3, p) == 1
        if q <= 61:
            a, b = primary_prime(p)
            ok &= (b % 2 == 0) and abs(a - b // 2) == x and abs(b // 6) == y
    # q = 3: 7 is not x^2 + 27 y^2 and 2 is not a cube mod 7
    ok &= not any(is_sq(7 - 27 * y * y) for y in range(0, 1)) and pow(2, 2, 7) != 1
    check(ok, 'Mersenne table')
    print(f'  Mersenne: the listed x^2 + 27y^2 representations of 2^q - 1 (q = 5..127) are correct, 2 is a cube via 3m = 1 '
          f'(mod q), and (q <= 61) the primary prime recomputed by a Euclidean gcd gives the same (x, y): {ok}; q = 3 excluded')
    # is every Mersenne prime > 7 known to be x^2 + 27 y^2 (checked for the 10 listed)?


# ------------------------------------------------------------------------------------------------ (4)
def gegen_ok(X, n, t):
    """t-design test by Gegenbauer sums: sum_{x,y in X} G_k(<x,y>) = 0 for 1 <= k <= t (unit vectors in R^n)."""
    z = sp.symbols('z')
    lam = sp.Rational(n - 2, 2)
    res = []
    for k in range(1, t + 1):
        if n == 2:
            G = sp.chebyshevt(k, z)
        else:
            G = sp.gegenbauer(k, lam, z)
        s = 0
        for x in X:
            for y in X:
                s += G.subs(z, sp.nsimplify(sum(a * b for a, b in zip(x, y))))
        res.append(sp.simplify(s) == 0)
    return res


def moment_ok(X, n, t):
    ys = sp.symbols('y0:%d' % n)
    r2 = sum(v * v for v in ys)
    res = []
    for k in range(1, t + 1):
        s = sp.expand(sum(sum(a * b for a, b in zip(x, ys)) ** k for x in X))
        if k % 2:
            res.append(s == 0)
        else:
            # constant c with s = c |y|^k: compare coefficient-wise
            c = sp.Rational(1)
            for j in range(1, k, 2):
                c *= sp.Rational(j)
            den = sp.Rational(1)
            for j in range(0, k // 2):
                den *= (n + 2 * j)
            c = c / den * len(X)
            res.append(sp.expand(s - c * r2 ** (k // 2)) == 0)
    return res


def strength(flags):
    t = 0
    for f in flags:
        if f:
            t += 1
        else:
            break
    return t


def part4():
    print('== (4) designs and identities')
    r3 = sp.sqrt(3)
    hexa = [(sp.cos(sp.pi * k / 3), sp.sin(sp.pi * k / 3)) for k in range(6)]
    hexa = [tuple(sp.nsimplify(c) for c in p) for p in hexa]
    first2 = [hexa[0], hexa[3]]
    cent4 = [hexa[1], hexa[2], hexa[4], hexa[5]]
    lips = []
    for i in range(4):
        for s in (1, -1):
            v = [0] * 4
            v[i] = s
            lips.append(tuple(sp.Integer(c) for c in v))
    tess = [tuple(sp.Rational(s, 2) for s in sg) for sg in product((1, -1), repeat=4)]
    cell24 = lips + tess
    sq = [(1, 0), (0, 1), (-1, 0), (0, -1)]
    octa = sq + [(1 / sp.sqrt(2), 1 / sp.sqrt(2)), (-1 / sp.sqrt(2), 1 / sp.sqrt(2)), (1 / sp.sqrt(2), -1 / sp.sqrt(2)),
                 (-1 / sp.sqrt(2), -1 / sp.sqrt(2))]
    sets = [('{+-1}', first2, 2), ('{+-w, +-w^2}', cent4, 2), ('hexagon', hexa, 2), ('16-cell', lips, 4),
            ('tesseract', tess, 4), ('24-cell', cell24, 4), ('square', sq, 2), ('octagon (d = 1 two-grid units)', octa, 2)]
    table = {}
    for name, X, n in sets:
        m = moment_ok(X, n, 7)
        g = gegen_ok(X, n, 7)
        check(strength(m) == strength(g), f'two design tests agree on the strength of {name}')
        table[name] = strength(m)
        print(f'  {name:32s}: design strength {strength(m)} (moment test) / {strength(g)} (Gegenbauer test); '
              f'2-design {m[1] if len(m) > 1 else None}, 4-design {m[3]}, 6-design {m[5]}')
    expected = {'{+-1}': 1, '{+-w, +-w^2}': 1, 'hexagon': 5, '16-cell': 3, 'tesseract': 3, '24-cell': 5, 'square': 3,
                'octagon (d = 1 two-grid units)': 7}
    check(table == expected, f'design strengths {table}')
    a, b = sp.symbols('a b')
    x = sp.symbols('x0:4')
    cand = sp.expand(a ** 4 + b ** 4 + (a + b) ** 4 - 2 * (a * a + a * b + b * b) ** 2) == 0
    liou = sp.expand(sum((x[i] + x[j]) ** 4 + (x[i] - x[j]) ** 4 for i in range(4) for j in range(i + 1, 4))
                     - 6 * sum(v * v for v in x) ** 2) == 0
    hex4 = sp.expand(sum((p[0] * a + p[1] * b) ** 4 for p in hexa) - sp.Rational(9, 4) * (a * a + b * b) ** 2) == 0
    c24 = sp.expand(sum(sum(c * v for c, v in zip(u, x)) ** 4 for u in cell24) - 3 * sum(v * v for v in x) ** 2) == 0
    octid = sp.expand(4 * a ** 4 + 4 * b ** 4 + (a + b) ** 4 + (a - b) ** 4 - 6 * (a * a + b * b) ** 2) == 0
    check(cand and liou and hex4 and c24 and octid, 'identities')
    print(f'  Candido: {cand}; Liouville: {liou}; sum over hexagon <u,x>^4 = (9/4)|x|^4: {hex4}; sum over 24-cell = 3|x|^4: '
          f'{c24}; octagon identity 6(a^2+b^2)^2 = 4a^4 + 4b^4 + (a+b)^4 + (a-b)^4: {octid}')

    def qmul(p, q):          # Hamilton product, (1, i, j, k)
        a1, b1, c1, d1 = p
        a2, b2, c2, d2 = q
        return (a1 * a2 - b1 * b2 - c1 * c2 - d1 * d2,
                a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2,
                a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2,
                a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2)
    # sanity: i j = k, j i = -k
    check(qmul((0, 1, 0, 0), (0, 0, 1, 0)) == (0, 0, 0, 1) and qmul((0, 0, 1, 0), (0, 1, 0, 0)) == (0, 0, 0, -1), 'qmul')
    units = [tuple(Fr(int(2 * c), 2) for c in u) for u in cell24]
    roots = set()
    for i in range(4):
        for j in range(i + 1, 4):
            for si in (1, -1):
                for sj in (1, -1):
                    v = [0, 0, 0, 0]
                    v[i], v[j] = si, sj
                    roots.add(tuple(Fr(c) for c in v))
    left = {qmul((1, 1, 0, 0), u) for u in units}
    right = {qmul(u, (1, 1, 0, 0)) for u in units}
    check(left == roots and right == roots and len(roots) == 24, '(1+i) H* = D4 roots')
    print(f'  (1+i) * Hurwitz units = the 24 roots of D4 (left: {left == roots}, right: {right == roots})')
    # negative controls
    bad1 = sp.expand(a ** 4 + b ** 4 + (a + b) ** 4 - 3 * (a * a + a * b + b * b) ** 2) == 0
    bad2 = sp.expand(sum((x[i] + x[j]) ** 4 + (x[i] - x[j]) ** 4 for i in range(4) for j in range(i + 1, 4))
                     - 5 * sum(v * v for v in x) ** 2) == 0
    bad3 = strength(moment_ok(cell24[:-1], 4, 7))
    print(f'  negative controls: Candido with 3 instead of 2: {bad1}; Liouville with 5 instead of 6: {bad2}; '
          f'24-cell minus one point has design strength {bad3}')
    check(not bad1 and not bad2 and bad3 == 0, 'negative controls (4)')


# ------------------------------------------------------------------------------------------------ (5), (6)
def part56():
    print('== (5) Eisenstein fourth powers mod 16, (6) sums of four')

    def p4(a, b, m):
        s = emul((a, b), (a, b))
        f = emul(s, s)
        return (f[0] % m, f[1] % m)
    res = sorted({p4(a, b, 16) for a in range(16) for b in range(16)})
    print('  fourth powers mod 16:', res)
    mu3 = {(1, 0), (0, 1), (15, 15)}
    nine = {(9 * x % 16, 9 * y % 16) for (x, y) in mu3}
    check(set(res) == {(0, 0)} | mu3 | nine, 'fourth powers mod 16 structure')
    ok = True
    for a in range(16):
        for b in range(16):
            f = p4(a, b, 16)
            if a % 2 == 0 and b % 2 == 0:
                ok &= f == (0, 0)
            else:
                u = (f[0] % 2, f[1] % 2)          # alpha^4 mod 2
                uu = {(1, 0): (1, 0), (0, 1): (0, 1), (1, 1): (15, 15)}[u]
                ok &= f in (uu, (9 * uu[0] % 16, 9 * uu[1] % 16))
    check(ok, '(5) u or 9u with u = alpha^4 mod 2')
    m6 = sorted({emul(emul(u, u), emul(u, u)) for u in UNITS})
    check(set(m6) == {(1, 0), (0, 1), (-1, -1)}, 'mu6^4 = mu3')
    print(f'  alpha in 2Z[w] -> 0; else u or 9u with u = alpha^4 mod 2: {ok}; mu6^4 = {m6} = {{1, w, w^2}}')

    def reach(resset, m, terms=4):
        cur = {(0, 0)}
        for _ in range(terms):
            cur = {((x + r[0]) % m, (y + r[1]) % m) for (x, y) in cur for r in resset}
        return sorted(x for (x, y) in cur if y == 0)
    rZw = reach(set(res), 16)
    rZ = sorted({(a + b + c + d) % 16 for a in (0, 1) for b in (0, 1) for c in (0, 1) for d in (0, 1)})
    first_grid = sorted({p4(a, b, 16) for a in range(16) for b in range(0, 16, 2)})
    rL = reach(set(first_grid), 16)
    print(f'  rational residues mod 16 of four fourth powers: over Z {rZ}; over the first grid Z[sqrt-3] {rL}; '
          f'over Z[w] {rZw}')
    check(rZw == [r for r in range(16) if r % 8 != 5], '(6) reachable mod 16')
    # negative control: a wrong residue set (add 5) removes the obstruction
    rbad = reach(set(res) | {(5, 0)}, 16)
    print(f'  negative control: adding the false residue 5 to the fourth-power set makes {sorted(set(rbad) - set(rZw))} reachable')


if __name__ == '__main__':
    part1()
    part2()
    part3()
    part4()
    part56()
    print('FAILURES:', len(FAIL))
    for f in FAIL[:20]:
        print('  ', f)
