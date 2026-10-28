"""The integer grid Z[sqrt-d] plus its (1/2,1/2)-centred copy: ring or ideal class; which grid carries a norm shell.

Companion to 05-knowledge/results/two_grids_ring_or_ideal_class_13_37_20261027.md (opus-2026-10-27-S66).

Notation.  L_d = Z[sqrt-d] (the first grid, coordinates (u, v) of u + v sqrt-d), G_d = L_d u ((1+sqrt-d)/2 + L_d) (the
owner's grid with its (1/2,1/2)-centred copy), so 2 G_d = {u + v sqrt-d : u = v mod 2} (the checkerboard of L_d).
Eisenstein elements are pairs (a, b) for a + b omega, omega^2 + omega + 1 = 0, norm a^2 - ab + b^2.

Sections:
  A  ring or ideal class (odd squarefree d <= 199): G_d is a ring iff d = 3 mod 4; for d = 1 mod 4 the checkerboard
     2G_d is an ideal of L_d with square (2), principal only for d = 1
  B  the two grids are the whole class group: h(-4d) = 2 for squarefree d = 1 mod 4 exactly at d = 5, 13, 37 (exact
     reduced-form search to D_MAX; complete by Baker-Stark); the reduced forms for those d are x^2 + dy^2 and
     2x^2 + 2xy + ((1+d)/2)y^2
  C  three grids agree: for d in {5, 13, 37} a split prime q is a first-grid norm iff q = 1 mod 4 iff q is a Gaussian
     norm (the centred form takes the rest); first counterexamples for other d = 1 mod 4; the owner's 5, 10, 25, 50, 100
     alternate between the two grids of Z[sqrt-5]
  D  which grid carries a shell: Z[omega] (1/3 : 2/3 on odd shells, all on the first grid on even shells), the Hurwitz
     order (Lipschitz 8 sigma(n) of 24 sigma(n)), d = 7 mod 8 (odd shells on the first grid), d = 3 mod 8, d > 3 (prime
     shells entirely on one grid; first-grid share of split primes -> 1/3); the order class numbers h(-4d)/h(-d)
  E  the grid of the primary Eisenstein prime is the cubic character of 2 (p = x^2 + 27y^2), all p = 1 mod 3 below P_MAX;
     Mersenne primes 2^q - 1 (q > 3) lie on the first grid
  F  the centred units complete the quartic designs: hexagon = {+-1} u {+-omega, +-omega^2}, 24-cell = 16-cell u
     tesseract; Candido's and Liouville's identities; (1+i) maps the Hurwitz units onto the D4 roots; the square grid's
     units are not a 4-design; Eisenstein fourth powers mod 16; mu_6^4 = mu_3 (the merge charges)
"""
import math
import sys
import time
from fractions import Fraction

D_MAX = 200_000      # section B: squarefree d = 1 mod 4 below D_MAX searched for h(-4d) = 2
Q_MAX = 200_000      # section C: split primes q below Q_MAX
P_MAX = 1_000_000    # section E: primes p = 1 mod 3 below P_MAX


def primes_upto(n):
    s = bytearray([1]) * (n + 1)
    s[0] = s[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]


def squarefree(n):
    i = 2
    while i * i <= n:
        if n % (i * i) == 0:
            return False
        i += 1
    return True


def is_square(n):
    return n >= 0 and math.isqrt(n) ** 2 == n


def rep(n, f):
    """Is n represented by the positive form f = (a, b, c) (a x^2 + b xy + c y^2)?  Exhaustive, exact."""
    a, b, c = f
    D = b * b - 4 * a * c
    ymax = math.isqrt(4 * a * n // (-D)) + 1
    for y in range(-ymax, ymax + 1):
        # a x^2 + b y x + c y^2 - n = 0
        disc = (b * y) ** 2 - 4 * a * (c * y * y - n)
        if disc < 0 or not is_square(disc):
            continue
        r = math.isqrt(disc)
        for s in (r, -r):
            if (-b * y + s) % (2 * a) == 0:
                return True
    return False


def reduced_forms(D):
    """All reduced primitive positive definite forms of discriminant D < 0."""
    out = []
    a = 1
    while 3 * a * a <= -D:
        for b in range(-a + 1, a + 1):
            if (b * b - D) % (4 * a):
                continue
            c = (b * b - D) // (4 * a)
            if c < a or math.gcd(math.gcd(a, abs(b)), c) != 1:
                continue
            if b < 0 and a == c:
                continue
            out.append((a, b, c))
        a += 1
    return out


def h(D):
    return len(reduced_forms(D))


def has_third_reduced_form(D, exclude):
    """True iff a reduced primitive form of discriminant D exists outside `exclude` (early exit)."""
    a = 1
    while 3 * a * a <= -D:
        for b in range(-a + 1, a + 1):
            if (b * b - D) % (4 * a):
                continue
            c = (b * b - D) // (4 * a)
            if c < a or math.gcd(math.gcd(a, abs(b)), c) != 1:
                continue
            if b < 0 and a == c:
                continue
            if (a, b, c) not in exclude:
                return True
        a += 1
    return False


def legendre(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def sqrt_mod(a, p):
    """Tonelli-Shanks; p an odd prime, a a nonzero quadratic residue."""
    a %= p
    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while legendre(z, p) != -1:
        z += 1
    m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, tt = 0, t
        while tt != 1:
            tt = tt * tt % p
            i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, b * b % p, t * b * b % p, r * b % p
    return r


def cornacchia(n, p):
    """A solution of x^2 + n y^2 = p (p prime) or None."""
    if legendre(-n, p) != 1:
        return None
    r0 = sqrt_mod(-n, p)
    if 2 * r0 < p:
        r0 = p - r0
    a, b = p, r0
    lim = math.isqrt(p)
    while b > lim:
        a, b = b, a % b
    rest = p - b * b
    if rest % n or not is_square(rest // n):
        return None
    return b, math.isqrt(rest // n)


# ---------------------------------------------------------------- Eisenstein arithmetic, (a, b) = a + b omega
def emul(x, y):
    a, b = x
    c, d = y
    return (a * c - b * d, a * d + b * c - b * d)


def enorm(x):
    a, b = x
    return a * a - a * b + b * b


UNITS6 = [(1, 0), (0, 1), (-1, -1), (-1, 0), (0, -1), (1, 1)]   # 1, omega, omega^2, -1, -omega, -omega^2


def first_grid(x):
    """a + b omega lies in Z[sqrt-3] iff b is even (a + b omega = (a - b/2) + (b/2) sqrt-3)."""
    return x[1] % 2 == 0


def primary(x):
    """Ireland-Rosen: a + b omega is primary iff a = 2 mod 3 and b = 0 mod 3."""
    return x[0] % 3 == 2 and x[1] % 3 == 0


def eisenstein_prime_above(p):
    """An Eisenstein prime of norm p (p = 1 mod 3), from 4p = U^2 + 3V^2: pi = (U + V sqrt-3)/2 = ((U+V)/2) + V omega."""
    s = cornacchia(3, p)                     # p = x^2 + 3y^2 always holds for p = 1 mod 3
    x, y = s
    pi = (x + y, 2 * y)                      # x + y sqrt-3 = (x + y) + 2y omega
    assert enorm(pi) == p
    return pi


# ================================================================ A
def section_A():
    rows = []
    for d in range(1, 200, 2):
        if not squarefree(d):
            continue
        # ring test: c = (1 + sqrt-d)/2; c^2 = (1-d)/4 + (1/2) sqrt-d lies in G_d iff (1-d)/4 - 1/2 is an integer
        c2_real = Fraction(1 - d, 4)
        ring = (c2_real - Fraction(1, 2)).denominator == 1
        assert ring == (d % 4 == 3), d
        if d % 4 == 1:
            # checkerboard 2G_d = span{(2,0), (1,1)}; closed under multiplication by sqrt-d: (u, v) -> (-d v, u)
            for (u, v) in ((2, 0), (1, 1)):
                uu, vv = -d * v, u
                assert (uu - vv) % 2 == 0, d
            # p2^2 = span of products {4, 2(1+sqrt-d), (1+sqrt-d)^2 = (1-d) + 2 sqrt-d} = 2 L_d exactly
            gens = [(4, 0), (2, 2), (1 - d, 2)]
            # 2 L_d = {u, v even}; every generator is in it, and (2,0), (0,2) are integer combinations of them
            assert all(g[0] % 2 == 0 and g[1] % 2 == 0 for g in gens)
            k = (1 - d) // 4
            assert (1 - d) % 4 == 0 and (gens[2][0] - k * 4, gens[2][1]) == (0, 2)
            principal = any(x * x + d * y * y == 2 for x in range(2) for y in range(2))
            assert principal == (d == 1)
            rows.append(d)
    return rows


# ================================================================ B
def section_B():
    hits = []
    count = 0
    for d in range(1, D_MAX, 4):
        if not squarefree(d):
            continue
        count += 1
        D = -4 * d
        if d == 1:
            continue
        two = {(1, 0, d), (2, 2, (1 + d) // 2)}
        if not has_third_reduced_form(D, two):
            hits.append(d)
    exact = {d: h(-4 * d) for d in (5, 13, 37, 17, 29, 41, 53, 61, 85, 21, 33, 65)}
    forms = {d: reduced_forms(-4 * d) for d in (5, 13, 37)}
    return count, hits, exact, forms


# ================================================================ C
def section_C():
    primes = primes_upto(Q_MAX)
    out = {}
    for d in (5, 13, 37):
        f0, f1 = (1, 0, d), (2, 2, (1 + d) // 2)
        n_split = 0
        for q in primes:
            if q == 2 or d % q == 0:
                continue
            if legendre(-d, q) != 1:
                continue
            n_split += 1
            first = rep(q, f0)
            second = rep(q, f1)
            gauss = rep(q, (1, 0, 1))
            assert first != second, (d, q)
            assert first == (q % 4 == 1) == gauss, (d, q)
            # genus characters: q = 1 mod 4 and (q/d) = 1 for the principal form
            assert first == (q % 4 == 1 and legendre(q, d) == 1)
        out[d] = n_split
    counter = {}
    for d in (17, 29, 41, 53, 61, 85, 21, 33, 65, 73):
        for q in primes:
            if q % 4 == 1 and d % q and legendre(-d, q) == 1 and not rep(q, (1, 0, d)):
                counter[d] = q
                break
    owner = {}
    f0, f1 = (1, 0, 5), (2, 2, 3)
    for n in (2, 5, 10, 25, 50, 100, 13, 26, 52, 37, 74, 148):
        d = 5
        if n in (13, 26, 52):
            d = 13
        if n in (37, 74, 148):
            d = 37
        g0, g1 = (1, 0, d), (2, 2, (1 + d) // 2)
        owner[n] = (d, rep(n, g0), rep(n, g1))
    return out, counter, owner


# ================================================================ D
def chi3(n):
    return (0, 1, -1)[n % 3]


def section_D():
    N = 20000
    # Eisenstein shells by grid
    first = [0] * (N + 1)
    cent = [0] * (N + 1)
    v = 0
    while 3 * v * v <= N:
        u = 0
        while u * u + 3 * v * v <= N:
            n = u * u + 3 * v * v
            first[n] += (2 if u else 1) * (2 if v else 1)
            u += 1
        v += 1
    V = 1
    while 3 * V * V <= 4 * N:
        U = 1
        while U * U + 3 * V * V <= 4 * N:
            n4 = U * U + 3 * V * V
            if n4 % 4 == 0:
                cent[n4 // 4] += 4
            U += 2
        V += 2
    for n in range(1, N + 1):
        r = 6 * sum(chi3(e) for e in range(1, n + 1) if n % e == 0) if n <= 3000 else None
        if r is not None:
            assert first[n] + cent[n] == r, n
        if n % 2:
            assert 2 * first[n] == cent[n], n
        else:
            assert cent[n] == 0, n
    eis_ok = N
    # Hurwitz: Lipschitz and half-integer counts by norm
    M = 300
    lip = [0] * (M + 1)
    half = [0] * (M + 1)
    R = math.isqrt(M)
    rng = range(-R, R + 1)
    sq = {x: x * x for x in rng}
    for a in rng:
        for b in rng:
            s2 = sq[a] + sq[b]
            if s2 > M:
                continue
            for c in rng:
                s3 = s2 + sq[c]
                if s3 > M:
                    continue
                for e in rng:
                    s4 = s3 + sq[e]
                    if s4 <= M:
                        lip[s4] += 1
    odd = [x for x in range(-2 * R - 1, 2 * R + 2) if x % 2]
    for a in odd:
        for b in odd:
            s2 = a * a + b * b
            if s2 > 4 * M:
                continue
            for c in odd:
                s3 = s2 + c * c
                if s3 > 4 * M:
                    continue
                for e in odd:
                    s4 = s3 + e * e
                    if s4 <= 4 * M and s4 % 4 == 0:
                        half[s4 // 4] += 1
    sigma = lambda n: sum(e for e in range(1, n + 1) if n % e == 0)
    for n in range(1, M + 1):
        if n % 2:
            assert lip[n] == 8 * sigma(n) and half[n] == 16 * sigma(n), n
        else:
            m = n
            while m % 2 == 0:
                m //= 2
            assert lip[n] == 24 * sigma(m) and half[n] == 0, n
    # d = 7 mod 8: odd shells on the first grid only; d = 3 mod 8, d > 3: prime shells on one grid
    def grid_counts(d, n):
        f = c = 0
        y = 0
        while d * y * y <= n:
            if is_square(n - d * y * y):
                x = math.isqrt(n - d * y * y)
                f += (2 if x else 1) * (2 if y else 1)
            y += 1
        Y = 1
        while d * Y * Y <= 4 * n:
            rest = 4 * n - d * Y * Y
            if is_square(rest) and math.isqrt(rest) % 2 == 1:
                c += 4
            Y += 2
        return f, c
    seven = {}
    for d in (7, 15, 23, 31, 39, 47):
        assert all(grid_counts(d, n)[1] == 0 for n in range(1, 3001, 2)), d
        seven[d] = True
    three = {}
    pr = primes_upto(300_000)
    for d in (11, 19, 43, 67):
        n_split = n_first = 0
        for p in pr:
            if p == 2 or d % p == 0 or legendre(-d, p) != 1:
                continue
            f, c = grid_counts(d, p)
            assert (f == 0) != (c == 0) and f + c == 4, (d, p)   # +-pi, +-pibar, all on one grid
            n_split += 1
            n_first += f > 0
        three[d] = (n_split, n_first)
    # order class numbers h(-4d) / h(-d)
    ratio = {}
    for d in range(3, 200, 4):
        if not squarefree(d):
            continue
        r = Fraction(h(-4 * d), h(-d))
        expect = 1 if d == 3 else (3 if d % 8 == 3 else 1)
        assert r == expect, (d, r)
        ratio[d] = r
    return eis_ok, M, seven, three, ratio


# ================================================================ E
def section_E():
    primes = primes_upto(P_MAX)
    n = mismatches = n_res = 0
    for p in primes:
        if p % 3 != 1:
            continue
        n += 1
        pi = eisenstein_prime_above(p)
        prim = [emul(u, pi) for u in UNITS6]
        prim = [x for x in prim if primary(x)]
        assert len(prim) == 1, p
        pi = prim[0]
        grid = first_grid(pi)
        cubic = pow(2, (p - 1) // 3, p) == 1
        form27 = any(is_square(p - 27 * y * y) for y in range(0, math.isqrt(p // 27) + 1))
        mismatches += not (grid == cubic == form27)
        n_res += cubic
    mers = []
    for q in (3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127):
        p = (1 << q) - 1
        sol = cornacchia(27, p)
        cubic = pow(2, (p - 1) // 3, p) == 1
        if sol:
            x, y = sol
            assert x * x + 27 * y * y == p
            pi = (x + 3 * y, 6 * y)            # x + 3y sqrt-3 = (x + 3y) + 6y omega, a prime of norm p
            prim = [emul(u, pi) for u in UNITS6 if primary(emul(u, pi))]
            assert len(prim) == 1
            mers.append((q, x, y, cubic, first_grid(prim[0])))
        else:
            mers.append((q, None, None, cubic, False))
    return n, n_res, mismatches, mers


# ================================================================ F
def section_F():
    import sympy as sp
    x, y, z, w, a, b = sp.symbols('x y z w a b', real=True)
    r3 = sp.sqrt(3)
    hexa = [(sp.cos(sp.pi * k / 3), sp.sin(sp.pi * k / 3)) for k in range(6)]
    first2 = [hexa[0], hexa[3]]
    cent2 = [hexa[1], hexa[2], hexa[4], hexa[5]]

    def power_sum(pts, k, vars_):
        return sp.expand(sum((sum(c * v for c, v in zip(p, vars_))) ** k for p in pts))

    def is_design_k(pts, k, vars_):
        s = power_sum(pts, k, vars_)
        r2 = sum(v * v for v in vars_)
        ratio = sp.simplify(s / r2 ** (k // 2))
        return ratio.free_symbols == set(), ratio
    res = {}
    for name, pts in (('hexagon first {+-1}', first2), ('hexagon centred {+-w,+-w^2}', cent2), ('hexagon (union)', hexa)):
        res[name] = [is_design_k(pts, k, (x, y))[0] for k in (2, 4, 6)]
    e = [tuple(1 if i == j else 0 for j in range(4)) for i in range(4)]
    lips = [tuple(s * c for c in v) for v in e for s in (1, -1)]
    tess = [tuple(sp.Rational(s1, 2) * 1 for s1 in signs) for signs in
            [(s0, s1, s2, s3) for s0 in (1, -1) for s1 in (1, -1) for s2 in (1, -1) for s3 in (1, -1)]]
    cell24 = lips + tess
    for name, pts in (('Lipschitz units (16-cell)', lips), ('centred units (tesseract)', tess), ('Hurwitz units (24-cell)', cell24)):
        res[name] = [is_design_k(pts, k, (x, y, z, w))[0] for k in (2, 4, 6)]
    sq4 = [(1, 0), (0, 1), (-1, 0), (0, -1)]
    res['Gaussian units (square)'] = [is_design_k(sq4, k, (x, y))[0] for k in (2, 4, 6)]
    # the Gaussian two-grid set G_1 = (1/2)(1+i)Z[i]: the square with the normalised centred points (+-1 +- i)/sqrt2 is the
    # regular octagon (a 7-design); its quartic identity is Liouville's at x3 = x4 = 0 (found by the S66 audit)
    r2 = sp.sqrt(2)
    octa = sq4 + [(s1 / r2, s2 / r2) for s1 in (1, -1) for s2 in (1, -1)]
    res['Gaussian square + normalised centred points (octagon)'] = [is_design_k(octa, k, (x, y))[0] for k in (2, 4, 6)]
    assert sp.expand(6 * (x * x + y * y) ** 2 - (4 * x ** 4 + 4 * y ** 4 + (x + y) ** 4 + (x - y) ** 4)) == 0
    candido = sp.expand(a ** 4 + b ** 4 + (a + b) ** 4 - 2 * (a * a + a * b + b * b) ** 2) == 0
    X = (x, y, z, w)
    liou = sp.expand(sum((X[i] + X[j]) ** 4 + (X[i] - X[j]) ** 4 for i in range(4) for j in range(i + 1, 4))
                     - 6 * (x * x + y * y + z * z + w * w) ** 2) == 0
    hurw4 = sp.expand(power_sum(cell24, 4, X) - 3 * (x * x + y * y + z * z + w * w) ** 2) == 0

    # quaternion product, (1+i) * Hurwitz units = the 24 D4 roots +-e_i +- e_j
    def qmul(p, q):
        a1, b1, c1, d1 = p
        a2, b2, c2, d2 = q
        return (a1 * a2 - b1 * b2 - c1 * c2 - d1 * d2, a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2,
                a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2, a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2)
    img = {tuple(sp.nsimplify(c) for c in qmul((1, 1, 0, 0), u)) for u in cell24}
    roots = set()
    for i in range(4):
        for j in range(i + 1, 4):
            for si in (1, -1):
                for sj in (1, -1):
                    v = [0, 0, 0, 0]
                    v[i], v[j] = si, sj
                    roots.add(tuple(sp.Integer(c) for c in v))
    d4 = img == roots and len(roots) == 24
    # Eisenstein fourth powers mod 16, and mu_6^4 = mu_3
    def epow(x_, k):
        r = (1, 0)
        for _ in range(k):
            r = emul(r, x_)
        return r
    fourth_mod16 = sorted({tuple(c % 16 for c in epow((s, t), 4)) for s in range(16) for t in range(16)})
    mu3 = sorted({epow(u, 4) for u in UNITS6})

    # rational residues mod 16 of sums of four fourth powers: bases in Z, in the first grid Z[sqrt-3] (b even), in Z[omega]
    def rational_residues(bases):
        S = {tuple(c % 16 for c in epow(x_, 4)) for x_ in bases}
        R = {(0, 0)}
        for _ in range(4):
            R = R | {((r0 + s0) % 16, (r1 + s1) % 16) for (r0, r1) in R for (s0, s1) in S}
        return sorted(n for n in range(16) if (n, 0) in R)
    box = [(s, t) for s in range(16) for t in range(16)]
    repair = {'Z': rational_residues([(s, 0) for s in range(16)]),
              'Z[sqrt-3]': rational_residues([x_ for x_ in box if x_[1] % 2 == 0]),
              'Z[omega]': rational_residues(box)}
    assert repair['Z'] == [0, 1, 2, 3, 4]
    assert repair['Z[sqrt-3]'] == [0, 1, 2, 3, 4, 9, 10, 11, 12]
    assert repair['Z[omega]'] == [n for n in range(16) if n not in (5, 13)]
    return res, candido, liou, hurw4, d4, fourth_mod16, mu3, repair


def main():
    t0 = time.time()
    print('== A. ring or ideal class: G_d = L_d u ((1+sqrt-d)/2 + L_d), odd squarefree d <= 199')
    rows = section_A()
    print(f'  G_d is a ring exactly for d = 3 mod 4 (checked every odd squarefree d <= 199); for each of the {len(rows)} '
          f'd = 1 mod 4 the checkerboard 2G_d = (2, 1+sqrt-d) is an ideal of L_d with square (2), principal only at d = 1')
    print('== B. the two grids are the whole class group')
    count, hits, exact, forms = section_B()
    print(f'  squarefree d = 1 mod 4 below {D_MAX} ({count} values): h(-4d) = 2 exactly at d = {hits}')
    print(f'  exact class numbers h(-4d): {exact}')
    for d, fs in forms.items():
        print(f'  reduced forms of discriminant {-4 * d}: {fs}')
    print('== C. three grids agree for d = 5, 13, 37 (split primes q below %d)' % Q_MAX)
    out, counter, owner = section_C()
    for d, n_split in out.items():
        print(f'  d = {d}: all {n_split} split primes q: q = x^2 + {d}y^2  <=>  q = 1 mod 4  <=>  q = a^2 + b^2; '
              f'otherwise q = 2x^2 + 2xy + {(1 + d) // 2}y^2')
    print(f'  first split q = 1 mod 4 that is NOT x^2 + dy^2, for other d = 1 mod 4 (h > 2): {counter}')
    print('  the owner\'s numbers in Z[sqrt-d] (d, first-grid form, centred form):')
    for n, (d, r0, r1) in owner.items():
        print(f'    n = {n:3d} in Q(sqrt-{d}): x^2+{d}y^2: {r0}; 2x^2+2xy+{(1 + d) // 2}y^2: {r1}')
    print('== D. which grid carries a norm shell')
    eis_ok, M, seven, three, ratio = section_D()
    print(f'  Z[omega], n <= {eis_ok}: odd shells split first : centred = 1 : 2 exactly; even shells lie on the first grid '
          f'(r = 6 sum chi_-3 checked to 3000)')
    print(f'  Hurwitz order, n <= {M}: odd n: Lipschitz 8 sigma(n), centred 16 sigma(n); even n: Lipschitz 24 sigma(odd part), '
          f'centred 0')
    print(f'  d = 7 mod 8 ({sorted(seven)}): every odd shell n <= 3000 lies on the first grid')
    for d, (ns, nf) in three.items():
        print(f'  d = {d} (= 3 mod 8): every split-prime shell p < 3*10^5 lies on one grid; first grid for {nf} of {ns} '
              f'= {nf / ns:.4f} (Chebotarev in the cubic ring class field: 1/3)')
    print(f'  h(-4d)/h(-d) for squarefree d = 3 mod 4 below 200: 1 at d = 3, 3 for d = 3 mod 8 (d > 3), 1 for d = 7 mod 8 '
          f'({len(ratio)} values checked)')
    print('== E. the grid of the primary Eisenstein prime is the cubic character of 2')
    n, n_res, mism, mers = section_E()
    print(f'  primes p = 1 mod 3 below {P_MAX}: {n}; primary pi on the first grid <=> 2 a cubic residue <=> p = x^2 + 27y^2: '
          f'{mism} mismatches; 2 is a cubic residue for {n_res} ({n_res / n:.4f})')
    for q, x, y, cubic, fg in mers:
        if x is None:
            print(f'  Mersenne 2^{q} - 1: not x^2 + 27y^2; 2 cubic residue: {cubic}')
        else:
            print(f'  Mersenne 2^{q} - 1 = {x}^2 + 27*{y}^2; 2 cubic residue: {cubic}; primary prime on the first grid: {fg}')
    print('== F. the centred units complete the quartic designs')
    res, candido, liou, hurw4, d4, f16, mu3, repair = section_F()
    for name, flags in res.items():
        print(f'  {name}: 2-design {flags[0]}, 4-design {flags[1]}, 6-design {flags[2]}')
    print(f'  Candido a^4 + b^4 + (a+b)^4 = 2(a^2+ab+b^2)^2: {candido}; Liouville 6(sum x_i^2)^2 = sum_(i<j) (x_i +- x_j)^4: '
          f'{liou}; sum over the 24 Hurwitz units of <u,x>^4 = 3|x|^4: {hurw4}; (1+i) maps the Hurwitz units onto the 24 '
          f'D4 roots: {d4}')
    print(f'  Eisenstein fourth powers mod 16 (pairs (a, b) for a + b omega): {f16}')
    print(f'  fourth powers of the six units: {mu3} = 1, omega, omega^2 (the merge charges of THM-4820)')
    for ring, rr in repair.items():
        print(f'  rational residues mod 16 of four fourth powers with bases in {ring}: {rr}')
    print('  -> the first grid adds 9..12 to Z\'s 0..4; the centred grid adds 6, 7, 8, 14, 15; only 5 and 13 (n = 5 mod 8) '
          'stay out')
    print(f'ALL CHECKS PASSED  [{time.time() - t0:.0f} s]')


if __name__ == '__main__':
    main()
