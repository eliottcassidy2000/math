#!/usr/bin/env python3
"""Audit F2 (independent referee code): exact arithmetic in F = Q(zeta15)(j)(w), j^2 = -11, w^3 + 2w^2 - 4w - 4 = 0.

Written from scratch for the audit; it does NOT import the repo's Field / Geo / make helpers.
Basis: zeta^a j^b w^c (a < 8, b < 2, c < 3), flattened index a*6 + b*3 + c (the same ordering as the
repo's Field class, so that Python's sort of coordinate tuples reproduces the repo's sorted direction lists;
coordinates in a basis are canonical, so equal elements give equal tuples).
Elements are stored as (tuple of 48 ints, positive int denominator), normalised by gcd.
Multiplication: schoolbook product in (zeta, j, w) followed by explicit reduction with
  w^3 = -2 w^2 + 4 w + 4,  j^2 = -11,  zeta^8 = zeta^7 - zeta^5 + zeta^4 - zeta^3 + zeta - 1  (Phi_15).
Inverse: exact linear solve of the multiplication matrix with python-flint fmpq_mat (fallback: Fractions).
"""
from fractions import Fraction
from math import gcd
from functools import reduce

N = 48
IDX = {(a, b, c): a * 6 + b * 3 + c for a in range(8) for b in range(2) for c in range(3)}
EXP = {v: k for k, v in IDX.items()}


def _norm(num, den):
    if den < 0:
        num = [-x for x in num]; den = -den
    g = reduce(gcd, num, den)
    if g > 1:
        num = [x // g for x in num]; den //= g
    return (tuple(num), den)


class E:
    __slots__ = ('n', 'd')

    def __init__(self, num, den=1):
        self.n, self.d = _norm(list(num), den)

    # ---- constructors
    @staticmethod
    def zero():
        return E([0] * N)

    @staticmethod
    def scalar(q):
        q = Fraction(q)
        v = [0] * N; v[0] = q.numerator
        return E(v, q.denominator)

    @staticmethod
    def mono(a, b, c, coeff=1):
        q = Fraction(coeff)
        v = [0] * N; v[IDX[(a, b, c)]] = q.numerator
        return E(v, q.denominator)

    # ---- basic ops
    def __add__(self, o):
        if not isinstance(o, E): o = E.scalar(o)
        L = self.d * o.d // gcd(self.d, o.d)
        f1, f2 = L // self.d, L // o.d
        return E([x * f1 + y * f2 for x, y in zip(self.n, o.n)], L)
    __radd__ = __add__

    def __neg__(self):
        return E([-x for x in self.n], self.d)

    def __sub__(self, o):
        if not isinstance(o, E): o = E.scalar(o)
        return self + (-o)

    def __rsub__(self, o):
        return E.scalar(o) - self

    def scale(self, q):
        q = Fraction(q)
        return E([x * q.numerator for x in self.n], self.d * q.denominator)

    def __mul__(self, o):
        if not isinstance(o, E):
            return self.scale(o)
        P = [[[0] * 5 for _ in range(3)] for _ in range(15)]
        A = [(EXP[k], x) for k, x in enumerate(self.n) if x]
        B = [(EXP[k], y) for k, y in enumerate(o.n) if y]
        for (a1, b1, c1), x in A:
            for (a2, b2, c2), y in B:
                P[a1 + a2][b1 + b2][c1 + c2] += x * y
        # reduce w: w^c for c = 4, 3 -> w^(c-3) * (-2 w^2 + 4 w + 4)
        for a in range(15):
            for b in range(3):
                row = P[a][b]
                for c in (4, 3):
                    t = row[c]
                    if t:
                        row[c] = 0; row[c - 1] += -2 * t; row[c - 2] += 4 * t; row[c - 3] += 4 * t
        # reduce j: j^2 = -11
        for a in range(15):
            for c in range(3):
                t = P[a][2][c]
                if t:
                    P[a][2][c] = 0; P[a][0][c] += -11 * t
        # reduce zeta from the top: zeta^a = zeta^(a-8) (zeta^7 - zeta^5 + zeta^4 - zeta^3 + zeta - 1)
        for a in range(14, 7, -1):
            for b in range(2):
                for c in range(3):
                    t = P[a][b][c]
                    if t:
                        P[a][b][c] = 0
                        P[a - 1][b][c] += t; P[a - 3][b][c] -= t; P[a - 4][b][c] += t
                        P[a - 5][b][c] -= t; P[a - 7][b][c] += t; P[a - 8][b][c] -= t
        v = [0] * N
        for (a, b, c), k in IDX.items():
            v[k] = P[a][b][c]
        return E(v, self.d * o.d)
    __rmul__ = __mul__

    def __pow__(self, k):
        if k < 0:
            return self.inv() ** (-k)
        r = E.scalar(1); x = self
        while k:
            if k & 1: r = r * x
            x = x * x; k >>= 1
        return r

    def __eq__(self, o):
        if not isinstance(o, E): o = E.scalar(o)
        return self.n == o.n and self.d == o.d

    def __hash__(self):
        return hash((self.n, self.d))

    def is_zero(self):
        return not any(self.n)

    def fracs(self):
        return tuple(Fraction(x, self.d) for x in self.n)

    def den(self):
        """lcm of the coordinate denominators."""
        return self.d  # normalised: gcd(all numerators, d) = 1, so d is the lcm of coordinate denominators

    def mat(self):
        """48x48 rational matrix of multiplication by self (column k = self * basis_k), as (int rows, den)."""
        cols = []
        for k in range(N):
            e = [0] * N; e[k] = 1
            cols.append(self * E(e))
        L = 1
        for c in cols: L = L * c.d // gcd(L, c.d)
        rows = [[cols[k].n[i] * (L // cols[k].d) for k in range(N)] for i in range(N)]
        return rows, L

    def inv(self):
        rows, L = self.mat()
        try:
            import flint
            M = flint.fmpq_mat(N, N, [flint.fmpq(x, L) for r in rows for x in r])
            rhs = flint.fmpq_mat(N, 1, [flint.fmpq(1)] + [flint.fmpq(0)] * (N - 1))
            sol = M.solve(rhs)
            fr = [Fraction(int(sol[i, 0].p), int(sol[i, 0].q)) for i in range(N)]
        except ImportError:
            fr = _solve_fractions([[Fraction(x, L) for x in r] for r in rows], [Fraction(1)] + [Fraction(0)] * (N - 1))
        D = 1
        for q in fr: D = D * q.denominator // gcd(D, q.denominator)
        out = E([int(q * D) for q in fr], D)
        assert (self * out) == E.scalar(1), 'exact inverse check'
        return out

    def __truediv__(self, o):
        if not isinstance(o, E):
            return self.scale(Fraction(1) / Fraction(o))
        return self * o.inv()


def _solve_fractions(A, b):
    n = len(A)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for i in range(n):
        piv = next(r for r in range(i, n) if M[r][i] != 0)
        M[i], M[piv] = M[piv], M[i]
        inv = 1 / M[i][i]
        M[i] = [x * inv for x in M[i]]
        for r in range(n):
            if r != i and M[r][i] != 0:
                f = M[r][i]
                M[r] = [x - f * y for x, y in zip(M[r], M[i])]
    return [M[i][n] for i in range(n)]


# ---------------------------------------------------------------- generators and automorphisms
ONE = E.scalar(1)
ZETA = E.mono(1, 0, 0)
J = E.mono(0, 1, 0)
W = E.mono(0, 0, 1)
ZPOW = [ONE]
for _ in range(1, 15):
    ZPOW.append(ZPOW[-1] * ZETA)
assert ZPOW[-1] * ZETA == ONE, 'zeta^15 = 1'
JPOW = [ONE, J]
WPOW = [ONE, W, W * W]


def substitute(x, zeta_exp, jsign):
    """Apply the automorphism zeta -> zeta^zeta_exp, j -> jsign*j, w -> w."""
    out = E.zero()
    for k, v in enumerate(x.n):
        if v:
            a, b, c = EXP[k]
            term = ZPOW[(a * zeta_exp) % 15] * (JPOW[b] if jsign > 0 or b == 0 else -J) * WPOW[c]
            out = out + term.scale(Fraction(v, x.d))
    return out


def conj(x):
    """Complex conjugation: zeta -> zeta^14, j -> -j, w -> w (w is real)."""
    return substitute(x, 14, -1)


def norm_cc(x):
    return x * conj(x)


# ---------------------------------------------------------------- the objects of THM-4690(5)
def build():
    O = {}
    alpha = ZPOW[3]                       # zeta_5
    phi = ONE + alpha + alpha ** 4        # golden ratio 1 + 2cos(2pi/5)
    root5 = phi.scale(2) - ONE
    eta = ONE + ZPOW[5]                   # 1 + zeta_3 = -zeta_3^2, a primitive 6th root of unity
    rho = (J + 5).scale(Fraction(1, 6))   # (5 + sqrt(-11))/6, |rho| = 1
    baseD = ONE + eta
    M = [E.zero(), ONE, eta, baseD, rho, rho * eta, rho * baseD]
    ports = [E.zero(), ONE, ONE + alpha, ONE + alpha + alpha ** 2, -(alpha ** 4),
             (E.scalar(3) - J * root5).scale(Fraction(1, 8)), E.scalar(Fraction(3, 4)),
             (E.scalar(3) + J * root5).scale(Fraction(1, 8))]

    def uniq(lst):
        out, seen = [], set()
        for x in lst:
            if x not in seen:
                seen.add(x); out.append(x)
        return out
    C = uniq(M + ports[:5])
    A = uniq(ports + [p + m for p in ports for m in M])
    d = ONE + alpha
    a = conj(d)
    h = W + 3
    u = (h - J) / (h + J)
    u2 = (h + phi - J) / (h + phi + J)
    m2 = (J.scale(11) + 20).scale(Fraction(-1, 19))
    m1 = (E.scalar(20) - J.scale(11)).scale(Fraction(1, 19))
    c = (ONE + a * m2 - a * a * m1 + a ** 3) / (d * a - ONE)
    b = c * d + a * m1 - a * a
    Z = d + u
    Y = -(Z / conj(Z))
    O.update(alpha=alpha, phi=phi, root5=root5, eta=eta, rho=rho, M=M, ports=ports, C=C, A=A, d=d, a=a, h=h,
             u=u, u2=u2, m1=m1, m2=m2, b=b, c=c, Z=Z, Y=Y)
    return O


def roots_of_unity():
    rs = ZPOW + [-x for x in ZPOW]
    assert len(set(rs)) == 30
    return rs


if __name__ == '__main__':
    import time
    t0 = time.time()
    O = build()
    u, b, c, a, d, Y = O['u'], O['b'], O['c'], O['a'], O['d'], O['Y']
    print('sizes |M| = %d, |C| = %d, |A| = %d' % (len(O['M']), len(O['C']), len(O['A'])))
    print('u conj(u) = 1:', norm_cc(u) == ONE, '; u2 conj(u2) = 1:', norm_cc(O['u2']) == ONE)
    print('u^3 = m2 u^2 + m1 u + 1:', u ** 3 == O['m2'] * u * u + O['m1'] * u + ONE)
    print('u^2 = a + b u + c Y (the u-square circuit):', u * u == a + b * u + c * Y)
    print('|d|^2 - 1 = phi:', d * a - ONE == O['phi'], '; phi^2 = phi + 1:', O['phi'] ** 2 == O['phi'] + ONE)
    print('b, c in K = Q(zeta15, j) (no w):', all(O['b'].n[k] == 0 and O['c'].n[k] == 0 for k in range(N) if EXP[k][2] > 0))
    print('denominators: b %d, c %d' % (b.d, c.d))
    print('spindle unit edges:', sum(1 for i in range(7) for k in range(i) if norm_cc(O['M'][i] - O['M'][k]) == ONE))
    print('[%.1f s]' % (time.time() - t0))
