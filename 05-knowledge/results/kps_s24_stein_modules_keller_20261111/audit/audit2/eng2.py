"""eng2.py - engine of the second audit (written from scratch; not the author's stein.py, not the first audit's tw.py).

Two coordinate rings, elements are dicts {key: coefficient}, coefficients are integers modulo PRIME:

  ring 'L' : Laurent polynomials C[x, 1/x, y, 1/y]; key (k, l) is the monomial x^k y^l.
             (C[x, y], C[x, 1/x, y] and C[x, 1/x, y, 1/y] are the subrings picked by the caller.)
  ring 'H' : C[x, y, 1/q], q = 1 + x y (the plane minus the hyperbola x y = -1);
             key (c, n): x^c q^n if c >= 0, y^(-c) q^n if c < 0  (c = charge deg x - deg y, n in Z).
             This is a basis because the charge-c part of the ring is x^c C[q, 1/q] (resp. y^(-c) C[q, 1/q]).

Linear algebra: sparse elimination modulo PRIME (rank, membership), and python-flint for larger dense ranks.
A full rank modulo a prime is a full rank over Q (rank mod p <= rank over Q), so injectivity statements are exact.
"""
from math import comb

PRIME = 2305843009213693951          # 2^61 - 1 (the author uses 2147483629, the first audit 2147483647 and 1000000007)


def set_prime(p):
    global PRIME
    PRIME = p


def inv(a):
    return pow(a % PRIME, PRIME - 2, PRIME)


# ---------------------------------------------------------------- generic dict arithmetic
def add(a, b, s=1):
    out = dict(a)
    for m, c in b.items():
        v = (out.get(m, 0) + s * c) % PRIME
        if v:
            out[m] = v
        else:
            out.pop(m, None)
    return out


def scal(a, s):
    s %= PRIME
    if not s:
        return {}
    return {m: c * s % PRIME for m, c in a.items()}


def lincomb(pairs):
    out = {}
    for s, a in pairs:
        s %= PRIME
        if not s:
            continue
        for m, c in a.items():
            v = (out.get(m, 0) + s * c) % PRIME
            if v:
                out[m] = v
            else:
                out.pop(m, None)
    return out


class RingL:
    """Laurent polynomials in x, y."""
    name = "L"

    def mono(self, k, l, c=1):
        return {(k, l): c % PRIME} if c % PRIME else {}

    def one(self):
        return {(0, 0): 1}

    def mul(self, a, b):
        out = {}
        for (k, l), c in a.items():
            for (k2, l2), c2 in b.items():
                m = (k + k2, l + l2)
                v = (out.get(m, 0) + c * c2) % PRIME
                if v:
                    out[m] = v
                else:
                    out.pop(m, None)
        return out

    def dx(self, a):
        out = {}
        for (k, l), c in a.items():
            if k:
                out[(k - 1, l)] = c * k % PRIME
        return {m: c for m, c in out.items() if c}

    def dy(self, a):
        out = {}
        for (k, l), c in a.items():
            if l:
                out[(k, l - 1)] = c * l % PRIME
        return {m: c for m, c in out.items() if c}


class RingH:
    """C[x, y, 1/q], q = 1 + xy, in the charge basis."""
    name = "H"

    def mono(self, c, n, coef=1):
        return {(c, n): coef % PRIME} if coef % PRIME else {}

    def one(self):
        return {(0, 0): 1}

    def x(self):
        return {(1, 0): 1}

    def y(self):
        return {(-1, 0): 1}

    def q(self, n=1):
        return {(0, n): 1}

    def mul(self, a, b):
        out = {}
        for (c1, n1), v1 in a.items():
            for (c2, n2), v2 in b.items():
                vv = v1 * v2 % PRIME
                if c1 * c2 >= 0:
                    m = (c1 + c2, n1 + n2)
                    v = (out.get(m, 0) + vv) % PRIME
                    if v:
                        out[m] = v
                    else:
                        out.pop(m, None)
                else:
                    mm = min(abs(c1), abs(c2))        # x^mm y^mm = (q - 1)^mm
                    for j in range(mm + 1):
                        m = (c1 + c2, n1 + n2 + j)
                        v = (out.get(m, 0) + vv * comb(mm, j) * (-1) ** (mm - j)) % PRIME
                        if v:
                            out[m] = v
                        else:
                            out.pop(m, None)
        return out

    def dx(self, a):
        terms = []
        for (c, n), v in a.items():
            if c >= 1:          # x^c q^n -> (c + n) x^(c-1) q^n - n x^(c-1) q^(n-1)
                terms.append(((c + n) * v, {(c - 1, n): 1}))
                terms.append((-n * v, {(c - 1, n - 1): 1}))
            else:               # y^m q^n -> n y^(m+1) q^(n-1)
                terms.append((n * v, {(c - 1, n - 1): 1}))
        return lincomb(terms)

    def dy(self, a):
        terms = []
        for (c, n), v in a.items():
            if c <= -1:
                m = -c
                terms.append(((m + n) * v, {(c + 1, n): 1}))
                terms.append((-n * v, {(c + 1, n - 1): 1}))
            else:
                terms.append((n * v, {(c + 1, n - 1): 1}))
        return lincomb(terms)


def power(R, a, n):
    r = R.one()
    for _ in range(n):
        r = R.mul(r, a)
    return r


def dual_fields(R, u, v, Jinv):
    """d_u, d_v for the pair (u, v) with Jacobian J = u_x v_y - u_y v_x and Jinv = 1/J in the ring."""
    ux, uy, vx, vy = R.dx(u), R.dy(u), R.dx(v), R.dy(v)
    J = add(R.mul(ux, vy), R.mul(uy, vx), -1)
    assert R.mul(J, Jinv) == R.one(), "Jinv is not the inverse of the Jacobian"

    def du(h):
        return R.mul(add(R.mul(R.dx(h), vy), R.mul(R.dy(h), vx), -1), Jinv)

    def dv(h):
        return R.mul(add(R.mul(ux, R.dy(h)), R.mul(uy, R.dx(h)), -1), Jinv)

    assert du(u) == R.one() and dv(v) == R.one() and not du(v) and not dv(u)
    return du, dv, J


# ---------------------------------------------------------------- sparse linear algebra modulo PRIME
class Echelon:
    """Row space with pivots at the largest key (in the order given by keyfun)."""

    def __init__(self, keyfun=None):
        self.piv = {}
        self.keyfun = keyfun if keyfun is not None else (lambda m: m)

    def reduce(self, row):
        row = dict(row)
        kf = self.keyfun
        while row:
            lead = max(row, key=kf)
            prow = self.piv.get(lead)
            if prow is None:
                return row, lead
            f = row[lead]                      # pivot rows are normalised to leading coefficient 1
            for m, c in prow.items():
                v = (row.get(m, 0) - f * c) % PRIME
                if v:
                    row[m] = v
                else:
                    row.pop(m, None)
        return row, None

    def insert(self, row):
        row, lead = self.reduce(row)
        if lead is None:
            return False
        s = inv(row[lead])
        self.piv[lead] = {m: c * s % PRIME for m, c in row.items()}
        return True

    def contains(self, row):
        r, lead = self.reduce(row)
        return lead is None

    def rank(self):
        return len(self.piv)


def rank_rows(rows, keyfun=None):
    E = Echelon(keyfun)
    for r in rows:
        E.insert(r)
    return E.rank()


def flint_rank(rows, cols=None):
    """rank of a list of sparse rows with python-flint (dense, modulo PRIME if it fits in 63 bits)."""
    from flint import nmod_mat
    if cols is None:
        cols = sorted({m for r in rows for m in r})
    idx = {m: i for i, m in enumerate(cols)}
    nr, nc = len(rows), len(cols)
    if nr == 0 or nc == 0:
        return 0
    ent = [0] * (nr * nc)
    for i, r in enumerate(rows):
        base = i * nc
        for m, c in r.items():
            j = idx.get(m)
            if j is not None:
                ent[base + j] = c % PRIME
    return nmod_mat(nr, nc, ent, PRIME).rank()


# ---------------------------------------------------------------- self-test against sympy
if __name__ == "__main__":
    import random
    import sympy as sp
    random.seed(7)
    X, Y = sp.symbols("x y")
    Q = 1 + X * Y

    def toH(expr_num, npow):
        """element (polynomial in x, y) / q^npow of ring H from a sympy polynomial: reduce with x y = q - 1."""
        R = RingH()
        out = {}
        for (i, j), c in sp.Poly(expr_num, X, Y).terms():
            m = min(i, j)
            base = R.mono(i - j, -npow, int(c))
            t = base
            if m:
                t = R.mul(base, power(R, add(R.q(), R.one(), -1), m))
            out = add(out, t)
        return out

    def toS(a):
        e = 0
        for (c, n), v in a.items():
            vv = v if v < PRIME // 2 else v - PRIME
            e += vv * (X ** c if c >= 0 else Y ** (-c)) * Q ** n
        return e

    R = RingH()
    bad = 0
    for trial in range(40):
        p1 = sum(random.randint(-3, 3) * X ** random.randint(0, 3) * Y ** random.randint(0, 3) for _ in range(4))
        p2 = sum(random.randint(-3, 3) * X ** random.randint(0, 3) * Y ** random.randint(0, 3) for _ in range(4))
        n1, n2 = random.randint(0, 3), random.randint(0, 3)
        a, b = toH(p1, n1), toH(p2, n2)
        e1, e2 = p1 / Q ** n1, p2 / Q ** n2
        chk = [(toS(a) - e1), (toS(R.mul(a, b)) - e1 * e2), (toS(R.dx(a)) - sp.diff(e1, X)), (toS(R.dy(a)) - sp.diff(e1, Y))]
        for cc in chk:
            if sp.simplify(sp.together(cc)) != 0:
                bad += 1
    print("ring H self-test against sympy (representation, product, d/dx, d/dy on 40 random pairs): %d failures" % bad)
    L = RingL()
    a = {(2, 1): 3, (-1, 2): 5}
    assert L.dx(a) == {(1, 1): 6, (-2, 2): PRIME - 5} and L.dy(a) == {(2, 0): 3, (-1, 1): 10}
    rows = [{1: 1, 2: 2}, {1: 2, 2: 4}, {3: 1}]
    assert rank_rows(rows) == 2 and flint_rank(rows) == 2
    print("ring L and linear algebra self-test: ok")
