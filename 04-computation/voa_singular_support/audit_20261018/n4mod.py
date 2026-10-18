"""Independent implementation (audit) of the small N=4 superconformal vertex algebra (NS sector),
universal vacuum module, with super-PBW monomials and the super-polynomial associated graded.

Generators: L (wt 2, even), J+, J3, J- (wt 1, even), G1, G2, Gb1, Gb2 (wt 3/2, odd).
Small N=4 relations (Pauli-matrix form, e.g. Eguchi-Taormina), rewritten in the J^+-,J^3 basis:
  [L_m,L_n] = (m-n)L_{m+n} + c/12 (m^3-m) d,   [L_m,J_n] = -n J_{m+n},   [L_m,G_r] = (m/2-r) G_{m+r}
  [J3_m,J3_n] = c/12 m d,  [J3_m,J+-_n] = +-J+-_{m+n},  [J+_m,J-_n] = 2 J3_{m+n} + c/6 m d
  [J3,G1] = -1/2 G1, [J3,G2] = +1/2 G2, [J+,G1] = -G2, [J-,G2] = -G1
  [J3,Gb1] = +1/2 Gb1, [J3,Gb2] = -1/2 Gb2, [J+,Gb2] = +Gb1, [J-,Gb1] = +Gb2
  {G1_r,Gb1_s} = 2L - 2(r-s) J3 + c/3 (r^2-1/4) d,  {G2_r,Gb2_s} = 2L + 2(r-s) J3 + c/3 (r^2-1/4) d
  {G1_r,Gb2_s} = -2(r-s) J-,   {G2_r,Gb1_s} = -2(r-s) J+,   {G,G} = {Gb,Gb} = 0.
(The super-Jacobi identity is checked numerically in n4_audit.py.)
Modes are (type, n2) with n2 = 2 * mode index.
"""
import sys
from fractions import Fraction

sys.setrecursionlimit(100000)

TYPES = ['L', 'J+', 'J3', 'J-', 'G1', 'G2', 'Gb1', 'Gb2']
TIDX = {t: i for i, t in enumerate(TYPES)}
DELTA2 = {'L': 4, 'J+': 2, 'J3': 2, 'J-': 2, 'G1': 3, 'G2': 3, 'Gb1': 3, 'Gb2': 3}   # 2*conformal weight
ODD = {'G1', 'G2', 'Gb1', 'Gb2'}
CH2 = {'L': 0, 'J+': 2, 'J3': 0, 'J-': -2, 'G1': -1, 'G2': 1, 'Gb1': 1, 'Gb2': -1}   # 2*J3 eigenvalue


def parity(t):
    return 1 if t in ODD else 0


def bracket_raw(X, Y, c, one):
    """Super bracket [X, Y] of modes X=(t,m2), Y=(u,n2): returns (dict mode->coef, central)."""
    t, m2 = X
    u, n2 = Y
    s2 = m2 + n2
    out = {}
    cen = 0 * one
    half = one / 2

    def add(tt, coef):
        if coef != 0:
            out[(tt, s2)] = out.get((tt, s2), 0) + coef

    m = Fraction(m2, 2)
    n = Fraction(n2, 2)

    def F(x):  # Fraction -> field
        return one * x.numerator / x.denominator

    if t == 'L':
        if u == 'L':
            add('L', F(m - n))
            if s2 == 0:
                cen = c * F(m ** 3 - m) / 12
        elif u in ('J+', 'J3', 'J-'):
            add(u, F(-n))
        else:
            add(u, F(m / 2 - n))
        return out, cen
    if u == 'L':
        o, ce = bracket_raw(Y, X, c, one)
        sign = -1  # even with anything: [X,Y] = -[Y,X] when one of them is even
        return {k: sign * v for k, v in o.items()}, sign * ce
    if t in ('J+', 'J3', 'J-') and u in ('J+', 'J3', 'J-'):
        if t == 'J3' and u == 'J3':
            if s2 == 0:
                cen = c * F(m) / 12
        elif t == 'J3':
            add(u, one if u == 'J+' else -one)
        elif u == 'J3':
            add(t, -one if t == 'J+' else one)
        elif t == 'J+' and u == 'J-':
            add('J3', 2 * one)
            if s2 == 0:
                cen = c * F(m) / 6
        elif t == 'J-' and u == 'J+':
            add('J3', -2 * one)
            if s2 == 0:
                cen = -c * F(n) / 6
        return out, cen
    if t in ('J+', 'J3', 'J-') and u in ODD:
        tab = {('J3', 'G1'): ('G1', -half), ('J3', 'G2'): ('G2', half),
               ('J+', 'G1'): ('G2', -one), ('J-', 'G2'): ('G1', -one),
               ('J3', 'Gb1'): ('Gb1', half), ('J3', 'Gb2'): ('Gb2', -half),
               ('J+', 'Gb2'): ('Gb1', one), ('J-', 'Gb1'): ('Gb2', one)}
        if (t, u) in tab:
            tt, cf = tab[(t, u)]
            add(tt, cf)
        return out, cen
    if t in ODD and u in ('J+', 'J3', 'J-'):
        o, ce = bracket_raw(Y, X, c, one)
        return {k: -v for k, v in o.items()}, -ce
    # odd-odd
    if t in ODD and u in ODD:
        pair = (t, u)
        if pair in (('Gb1', 'G1'), ('Gb2', 'G2'), ('Gb2', 'G1'), ('Gb1', 'G2')):
            return bracket_raw(Y, X, c, one)  # symmetric
        r, s = m, n
        if pair == ('G1', 'Gb1'):
            add('L', 2 * one)
            add('J3', F(-2 * (r - s)))
            if s2 == 0:
                cen = c * F(r * r - Fraction(1, 4)) / 3
        elif pair == ('G2', 'Gb2'):
            add('L', 2 * one)
            add('J3', F(2 * (r - s)))
            if s2 == 0:
                cen = c * F(r * r - Fraction(1, 4)) / 3
        elif pair == ('G1', 'Gb2'):
            add('J-', F(-2 * (r - s)))
        elif pair == ('G2', 'Gb1'):
            add('J+', F(-2 * (r - s)))
        return {k: v for k, v in out.items() if v != 0}, cen
    raise ValueError((X, Y))


def is_creation(X):
    t, n2 = X
    return n2 <= -DELTA2[t]


def modekey(X):
    return (X[1], TIDX[X[0]])   # most negative mode first


class N4Vac:
    def __init__(self, c, one):
        self.c = c * one
        self.one = one
        self.cache = {}
        self.bcache = {}

    def br(self, X, Y):
        key = (X, Y)
        r = self.bcache.get(key)
        if r is None:
            o, ce = bracket_raw(X, Y, self.c, self.one)
            r = ({k: v for k, v in o.items() if v != 0}, ce)
            self.bcache[key] = r
        return r

    def act(self, X, mono):
        key = (X, mono)
        r = self.cache.get(key)
        if r is not None:
            return r
        one = self.one
        res = {}

        def acc(m, cf):
            v = res.get(m, 0) + cf
            if v == 0:
                res.pop(m, None)
            else:
                res[m] = v
        if not mono:
            if is_creation(X):
                res = {(X,): one}
        else:
            Y = mono[0]
            rest = mono[1:]
            if is_creation(X) and (modekey(X) < modekey(Y) or (X == Y and X[0] not in ODD)):
                res = {(X,) + mono: one}
            elif X == Y and X[0] in ODD:
                # X X rest = 1/2 {X,X} rest
                o, ce = self.br(X, X)
                for Z, cz in o.items():
                    for m3, c3 in self.act(Z, rest).items():
                        acc(m3, cz * c3 / 2)
                if ce != 0:
                    acc(rest, ce / 2)
            else:
                sgn = -1 if (X[0] in ODD and Y[0] in ODD) else 1
                for m2, c2 in self.act(X, rest).items():
                    for m3, c3 in self.act(Y, m2).items():
                        acc(m3, sgn * c2 * c3)
                o, ce = self.br(X, Y)
                for Z, cz in o.items():
                    for m3, c3 in self.act(Z, rest).items():
                        acc(m3, cz * c3)
                if ce != 0:
                    acc(rest, ce)
        self.cache[key] = res
        return res

    def actv(self, X, vec):
        out = {}
        for mono, cf in vec.items():
            for m2, c2 in self.act(X, mono).items():
                v = out.get(m2, 0) + cf * c2
                if v == 0:
                    out.pop(m2, None)
                else:
                    out[m2] = v
        return out


# ---------------- super polynomial ring (associated graded) ----------------
def sort_with_sign(lst):
    """Sort modes by modekey; return (tuple, sign) with sign from permuting odd elements; None if odd repeat."""
    arr = list(lst)
    sign = 1
    # insertion sort tracking odd transpositions
    for i in range(1, len(arr)):
        j = i
        while j > 0 and modekey(arr[j - 1]) > modekey(arr[j]):
            if arr[j - 1][0] in ODD and arr[j][0] in ODD:
                sign = -sign
            arr[j - 1], arr[j] = arr[j], arr[j - 1]
            j -= 1
    for i in range(1, len(arr)):
        if arr[i] == arr[i - 1] and arr[i][0] in ODD:
            return None, 0
    return tuple(arr), sign


def smul(p1, p2):
    out = {}
    for a, ca in p1.items():
        for b, cb in p2.items():
            m, s = sort_with_sign(a + b)
            if m is None:
                continue
            v = out.get(m, 0) + s * ca * cb
            if v == 0:
                out.pop(m, None)
            else:
                out[m] = v
    return out


def sd(p):
    """Even derivation d x_{u,j} = (j+1) x_{u,j+1}; variable (t, n2): j = (-n2 - DELTA2[t])/2."""
    out = {}
    for m, cf in p.items():
        for i, X in enumerate(m):
            t, n2 = X
            j = (-n2 - DELTA2[t]) // 2
            new = list(m)
            new[i] = (t, n2 - 2)
            mm, s = sort_with_sign(new)
            if mm is None:
                continue
            v = out.get(mm, 0) + s * (j + 1) * cf
            if v == 0:
                out.pop(mm, None)
            else:
                out[mm] = v
    return out


def wt2(mono):
    return -sum(n2 for _, n2 in mono)


def ch2(mono):
    return sum(CH2[t] for t, _ in mono)


def li(mono):
    return sum((-n2 - DELTA2[t]) // 2 for t, n2 in mono)


def nodd(mono):
    return sum(1 for t, _ in mono if t in ODD)


def creation_modes(maxw2):
    out = []
    for t in TYPES:
        n2 = -DELTA2[t]
        while -n2 <= maxw2:
            out.append((t, n2))
            n2 -= 2
    return sorted(out, key=modekey)


def pbw_basis(maxw2):
    """All super-PBW monomials (sorted by modekey, odd modes not repeated) of doubled weight <= maxw2."""
    modes = creation_modes(maxw2)
    res = {}

    def rec(start, cur, w):
        res.setdefault(w, []).append(tuple(cur))
        for i in range(start, len(modes)):
            X = modes[i]
            w2 = -X[1]
            if w + w2 > maxw2:
                continue
            nxt = i + 1 if X[0] in ODD else i
            cur.append(X)
            rec(nxt, cur, w + w2)
            cur.pop()
    rec(0, [], 0)
    return res
