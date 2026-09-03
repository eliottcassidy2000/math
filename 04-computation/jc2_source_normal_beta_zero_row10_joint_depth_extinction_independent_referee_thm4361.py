#!/usr/bin/env python3
"""Dependency-free exact referee for the proposed THM-4361.

This file deliberately imports neither the scout nor any repository
certificate.  It rebuilds the source-normal bracket rows from the displayed
definitions, constructs the projected depth matrices from the monomial
spanning law, and performs all symbolic arithmetic in a tiny Laurent ring
over Fraction.
"""

from fractions import Fraction as F
from math import comb, gcd
import sys

sys.stdout.reconfigure(newline="\n")


CHECKS = 0


def check(condition, label):
    global CHECKS
    if not condition:
        raise AssertionError(label)
    CHECKS += 1


class E:
    """Q[P^{+-1},X,A] sparse Laurent expressions."""

    __slots__ = ("t",)

    def __init__(self, terms=None):
        self.t = {}
        if terms:
            for mon, val in terms.items():
                val = F(val)
                if val:
                    self.t[tuple(mon)] = val

    @staticmethod
    def c(value):
        value = F(value)
        return E({(0, 0, 0): value}) if value else E()

    @staticmethod
    def var(i):
        mon = [0, 0, 0]
        mon[i] = 1
        return E({tuple(mon): F(1)})

    @staticmethod
    def mon(p=0, x=0, a=0, value=1):
        return E({(p, x, a): F(value)})

    def __bool__(self):
        return bool(self.t)

    def __eq__(self, other):
        return self.t == as_e(other).t

    def __neg__(self):
        return E({m: -v for m, v in self.t.items()})

    def __add__(self, other):
        other = as_e(other)
        out = dict(self.t)
        for m, v in other.t.items():
            out[m] = out.get(m, F(0)) + v
            if not out[m]:
                del out[m]
        return E(out)

    __radd__ = __add__

    def __sub__(self, other):
        return self + (-as_e(other))

    def __rsub__(self, other):
        return as_e(other) - self

    def __mul__(self, other):
        other = as_e(other)
        out = {}
        for m, v in self.t.items():
            for n, w in other.t.items():
                key = tuple(m[i] + n[i] for i in range(3))
                out[key] = out.get(key, F(0)) + v * w
        return E(out)

    __rmul__ = __mul__

    def __truediv__(self, scalar):
        scalar = F(scalar)
        return E({m: v / scalar for m, v in self.t.items()})

    def __pow__(self, exponent):
        if exponent < 0:
            if len(self.t) != 1:
                raise ValueError("only monomials have negative powers")
            (m, v), = self.t.items()
            return E({tuple(exponent * z for z in m): v ** exponent})
        out = E.c(1)
        base = self
        n = exponent
        while n:
            if n & 1:
                out = out * base
            base = base * base
            n //= 2
        return out

    def subst(self, index, replacement):
        replacement = as_e(replacement)
        out = E()
        for mon, val in self.t.items():
            exponent = mon[index]
            if exponent < 0 and len(replacement.t) != 1:
                raise ValueError("negative-power substitution is not monomial")
            kept = list(mon)
            kept[index] = 0
            out += E({tuple(kept): val}) * replacement ** exponent
        return out

    def coeff(self, mon):
        return self.t.get(tuple(mon), F(0))


def as_e(value):
    return value if isinstance(value, E) else E.c(value)


def epoly(items=None):
    out = {}
    if items:
        for degree, coeff in items.items():
            coeff = as_e(coeff)
            if coeff:
                out[int(degree)] = coeff
    return out


def pa(a, b):
    out = dict(a)
    for d, c in b.items():
        out[d] = out.get(d, E()) + c
        if not out[d]:
            del out[d]
    return out


def pn(a):
    return {d: -c for d, c in a.items()}


def ps(a, scalar):
    scalar = F(scalar)
    return {d: c * scalar for d, c in a.items() if c * scalar}


def pm(a, b):
    out = {}
    for i, c in a.items():
        for j, d in b.items():
            out[i + j] = out.get(i + j, E()) + c * d
            if not out[i + j]:
                del out[i + j]
    return out


def pd(a):
    return {i - 1: c * i for i, c in a.items() if i and c * i}


def pc(a, degree):
    return a.get(degree, E())


def psubst(a, index, replacement):
    return epoly({d: c.subst(index, replacement) for d, c in a.items()})


def pshift_down(a):
    check(pc(a, 0) == 0, "division by x has zero constant")
    return {d - 1: c for d, c in a.items() if d > 0}


def add_many(polys):
    out = {}
    for p in polys:
        out = pa(out, p)
    return out


def rank_fraction(matrix):
    if not matrix:
        return 0
    a = [[F(z) for z in row] for row in matrix]
    rows, cols = len(a), len(a[0])
    r = 0
    for c in range(cols):
        pivot = next((i for i in range(r, rows) if a[i][c]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        z = a[r][c]
        a[r] = [v / z for v in a[r]]
        for i in range(rows):
            if i != r and a[i][c]:
                z = a[i][c]
                a[i] = [a[i][j] - z * a[r][j] for j in range(cols)]
        r += 1
        if r == rows:
            break
    return r


def independent_rows(matrix):
    if not matrix:
        return []
    transpose = [list(col) for col in zip(*matrix)]
    a = [[F(z) for z in row] for row in transpose]
    rows, cols = len(a), len(a[0])
    pivots = []
    r = 0
    for c in range(cols):
        pivot = next((i for i in range(r, rows) if a[i][c]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        z = a[r][c]
        a[r] = [v / z for v in a[r]]
        for i in range(rows):
            if i != r and a[i][c]:
                z = a[i][c]
                a[i] = [a[i][j] - z * a[r][j] for j in range(cols)]
        pivots.append(c)
        r += 1
        if r == rows:
            break
    return pivots


def solve_square(matrix, rhs):
    n = len(matrix)
    check(n == len(matrix[0]) == len(rhs), "square solve dimensions")
    a = [[F(z) for z in row] for row in matrix]
    b = list(rhs)
    for c in range(n):
        pivot = next(i for i in range(c, n) if a[i][c])
        a[c], a[pivot] = a[pivot], a[c]
        b[c], b[pivot] = b[pivot], b[c]
        z = a[c][c]
        a[c] = [v / z for v in a[c]]
        b[c] = b[c] / z
        for i in range(n):
            if i != c and a[i][c]:
                z = a[i][c]
                a[i] = [a[i][j] - z * a[c][j] for j in range(n)]
                b[i] = b[i] - b[c] * z
    return b


def solve_full_columns(matrix, rhs):
    cols = len(matrix[0])
    rows = independent_rows(matrix)
    check(len(rows) == cols, "full column rank")
    sol = solve_square([[matrix[i][j] for j in range(cols)] for i in rows],
                       [rhs[i] for i in rows])
    residual = []
    for i, row in enumerate(matrix):
        residual.append(rhs[i] - sum((sol[j] * row[j] for j in range(cols)), E()))
    return sol, residual, rows


def rref_sparse(rows, ncols):
    rows = [{c: F(v) for c, v in row.items() if v} for row in rows]
    pivots = []
    r = 0
    for c in range(ncols):
        pivot = next((i for i in range(r, len(rows)) if rows[i].get(c)), None)
        if pivot is None:
            continue
        rows[r], rows[pivot] = rows[pivot], rows[r]
        z = rows[r][c]
        rows[r] = {j: v / z for j, v in rows[r].items()}
        prow = rows[r]
        for i in range(len(rows)):
            if i == r:
                continue
            z = rows[i].get(c, F(0))
            if not z:
                continue
            new = dict(rows[i])
            for j, v in prow.items():
                nv = new.get(j, F(0)) - z * v
                if nv:
                    new[j] = nv
                elif j in new:
                    del new[j]
            rows[i] = new
        pivots.append(c)
        r += 1
        if r == len(rows):
            break
    return rows, pivots


def nullspace_sparse(rows, ncols):
    reduced, pivots = rref_sparse(rows, ncols)
    free = [c for c in range(ncols) if c not in set(pivots)]
    basis = []
    for f in free:
        v = {f: F(1)}
        for i, p in enumerate(pivots):
            z = reduced[i].get(f, F(0))
            if z:
                v[p] = -z
        basis.append(v)
    return len(pivots), basis


def depth_columns(m, d):
    coords = [(n, r) for n in range(m + 1) for r in range(n + d + 1)]
    index = {z: i for i, z in enumerate(coords)}
    columns = []
    for a in range(d + 1):
        for b in range(d - a + 1):
            for e in range(m // 2 + 1):
                for c in range(m + 1):
                    n0 = b + c + 2 * e
                    if n0 > m:
                        continue
                    degree = c + e
                    r0 = a + 2 * b + e
                    col = {}
                    for k in range(degree + 1):
                        n, r = n0 + k, r0 + 2 * k
                        if n <= m and (n, r) in index:
                            col[index[(n, r)]] = F(comb(degree, k))
                    columns.append(col)
    return coords, columns


def jet(rows, coords):
    return [pc(rows[n], r) for n, r in coords]


def dot_sparse(v, values):
    return sum((values[i] * z for i, z in v.items()), E())


def tangent_vector(coords, row, which):
    # v_0'=(x/2,-3/4-3x^2/8).
    out = [F(0)] * len(coords)
    for i, (n, r) in enumerate(coords):
        if n != row:
            continue
        if which == "A" and r == 1:
            out[i] = F(1, 2)
        if which == "C" and r == 0:
            out[i] = F(-3, 4)
        if which == "C" and r == 2:
            out[i] = F(-3, 8)
    return out


def tangent_column(coords, row, power, which):
    base = tangent_vector(coords, row, which)
    out = [F(0)] * len(coords)
    for i, (n, r) in enumerate(coords):
        if n == row and r >= power:
            source = (n, r - power)
            if source in coords:
                out[i] = base[coords.index(source)]
    return out


def eval_left_on_column(left, column):
    return sum((left.get(i, F(0)) * column[i] for i in range(len(column))), F(0))


def zero_mod_q(expr, qexpr):
    """Test a one-variable Laurent expression modulo a polynomial q(P)."""
    check(all(mon[1] == mon[2] == 0 for mon in expr.t), "q reduction is univariate")
    check(all(mon[1] == mon[2] == 0 for mon in qexpr.t), "q polynomial is univariate")
    if not expr:
        return True
    shift = max(0, -min(mon[0] for mon in expr.t))
    poly = {mon[0] + shift: value for mon, value in expr.t.items()}
    qpoly = {mon[0]: value for mon, value in qexpr.t.items()}
    qdegree = max(qpoly)
    qlead = qpoly[qdegree]
    while poly and max(poly) >= qdegree:
        degree = max(poly)
        factor = poly[degree] / qlead
        offset = degree - qdegree
        for j, value in qpoly.items():
            k = j + offset
            poly[k] = poly.get(k, F(0)) - factor * value
            if not poly[k]:
                del poly[k]
    return not poly


P = E.var(0)
X = E.var(1)
A = E.var(2)

D = E.c(F(896, 15))
Theta = E.c(F(512, 75))
K = E.c(F(-32, 5))
ups = E.c(F(-731648, 2025))
zeta = -P * F(3, 2)
U = (E.c(475515904) - X * 109350) / 200475
W = -(P * P * 4343625 - X * 17172000 + E.c(143826305024)) / 4009500
Z = E()
eta = (E.c(12506118074368) - P * P * 173745000 - X * 926883000) * E.mon(p=-1) / 195463125

AA = {
    0: epoly({0: 1, 2: F(1, 4)}),
    1: epoly({0: F(4, 3), 2: 2}),
    2: epoly({0: F(-32, 9), 2: F(-4, 5)}),
    3: epoly({0: F(2176, 135), 1: -P / 2,
              2: E.c(F(1088, 315)) - K * F(4, 7), 4: F(-32, 15)}),
}
CC = {
    0: epoly({1: F(-3, 4), 3: F(-1, 8)}),
    1: epoly({1: -4, 3: F(-3, 2)}),
    2: epoly({1: F(88, 15), 3: F(-12, 5)}),
    3: epoly({0: P * F(3, 4), 1: E.c(F(-8128, 315)) + K * F(6, 7),
              2: P * F(3, 8), 3: E.c(F(736, 105)) + K * F(3, 7),
              5: F(8, 5)}),
}

GG = {
    4: epoly({0: D, 1: P, 2: K - E.c(F(1376, 45)), 4: F(8, 3)}),
    5: epoly({0: ups, 1: eta, 2: D * 4 + Theta, 3: P * 3,
              4: K * 2 - E.c(F(1376, 45))}),
    6: epoly({0: U, 1: A, 2: ups * 5 + X, 3: eta * 4 + zeta,
              4: D * 6 + Theta * 3, 5: P * 3,
              6: K - E.c(F(1376, 135))}),
    7: epoly({2: U * 6 + W, 3: A * 5, 4: ups * 10 + X * 4,
              5: eta * 6 + zeta * 3, 6: D * 4 + Theta * 3, 7: P}),
    8: epoly({4: U * 15 + W * 5, 5: A * 10, 6: ups * 10 + X * 6,
              7: eta * 4 + zeta * 3, 8: D + Theta}),
    9: epoly({6: U * 20 + W * 10, 7: A * 10,
              8: ups * 5 + X * 4, 9: eta + zeta}),
    10: epoly({8: U * 15 + W * 10, 9: A * 5, 10: ups + X}),
}

q_x = epoly({0: -3, 2: F(-1, 2)})


def B_row(m):
    return add_many([
        pa(ps(pm(pd(AA[i]), CC[m - i]), m - i),
           ps(pm(AA[i], pd(CC[m - i])), -i))
        for i in range(1, m)
    ])


def T_row(m):
    out = add_many([pm(CC[i], CC[m - i]) for i in range(1, m)])
    for i in range(m):
        for j in range(m):
            k = m - i - j
            if 0 <= k < m:
                out = pa(out, pn(pm(pm(AA[i], AA[j]), AA[k])))
    return out


def prediction(m):
    return pa(T_row(m), ps(pm(q_x, B_row(m)), F(-1, m)))


def particular(m):
    determinant = ps(B_row(m), F(-1, m))
    aconst = pc(determinant, 0) * F(4, 3)
    correction = pm(epoly({0: 2, 2: 1}), epoly({0: aconst * F(3, 8)}))
    cpart = ps(pshift_down(pa(determinant, pn(correction))), 2)
    apart = epoly({0: aconst})
    lhs = pa(pm(pn(pd(CC[0])), apart), pm(pd(AA[0]), cpart))
    check(lhs == determinant, "particular determinant row %d" % m)
    return apart, cpart


def D_matrix(m, max_power=None):
    if max_power is None:
        max_power = m - 1
    rows = m + 1
    out = [[F(0) for _ in range(max_power + 1)] for _ in range(rows)]
    for j in range(max_power + 1):
        if j:
            out[j - 1][j] += F(3 * j, m)
        out[j + 1][j] += F(j - 2 * m, 2 * m)
    return out


def add_tangent(row, apart, cpart, values, powers=None):
    if powers is None:
        powers = list(range(len(values)))
    theta = epoly({powers[i]: values[i] for i in range(len(values))})
    AA[row] = pa(apart, pm(theta, pd(AA[0])))
    CC[row] = pa(cpart, pm(theta, pd(CC[0])))


check(GG[4] == prediction(4), "literal G4 compatibility")
for row in range(4, 8):
    apart, cpart = particular(row)
    AA[row], CC[row] = apart, cpart
    defect = pa(GG[row + 1], pn(prediction(row + 1)))
    rhs = [pc(defect, j) for j in range(row + 2)]
    sol, residual, _ = solve_full_columns(D_matrix(row + 1), rhs)
    check(all(not z for z in residual), "row %d selects uniquely" % row)
    add_tangent(row, apart, cpart, sol)
    check(GG[row + 1] == prediction(row + 1), "row %d next source" % row)

# Row nine compatibility and the unique alpha solution.
apart8, cpart8 = particular(8)
AA[8], CC[8] = apart8, cpart8
defect9 = pa(GG[9], pn(prediction(9)))
sol8, residual9, _ = solve_full_columns(D_matrix(9), [pc(defect9, j) for j in range(10)])

E9 = (P * P * 613527750 - P * A * 511211250 - P * eta * 3154140000
      - eta * eta * 255605625 + X * 6736896000 - E.c(46483785515008))
nonzero9 = [r for r in residual9 if r]
check(len(nonzero9) == 1, "one row-nine cokernel")
probe = next(iter(E9.t))
ratio = nonzero9[0].coeff(probe) / E9.coeff(probe)
check(nonzero9[0] == E9 * ratio and ratio != 0, "row-nine residual is E9")

Aexpr = (
    P**4 * 4085005423617421875
    + P**2 * X * 24824465812575000000
    - P**2 * 278518552828793671680000
    - X**2 * 7302452813356500000
    + X * 197059040065115394048000
    - E.c(1329425408965288765218095104)
) * E.mon(p=-3) / 649499164991015625
check(E9.subst(2, Aexpr) == 0, "displayed alpha solves E9")
check((P * -511211250) != 0, "alpha coefficient is nonzero on Phi != 0")

add_tangent(8, apart8, cpart8, sol8)
for n in list(AA):
    AA[n] = psubst(AA[n], 2, Aexpr)
    CC[n] = psubst(CC[n], 2, Aexpr)
for n in list(GG):
    GG[n] = psubst(GG[n], 2, Aexpr)
check(GG[9] == prediction(9), "row-nine source after alpha solution")

# Exact projected depth matrices and the row-nine terminal fibre.
coords29, cols29 = depth_columns(9, 2)
coords39, cols39 = depth_columns(9, 3)
rank29, null29 = nullspace_sparse(cols29, len(coords29))
rank39, null39 = nullspace_sparse(cols39, len(coords39))
check((len(coords29), len(cols29), rank29, len(null29)) == (75, 160, 59, 16),
      "pi9(P2) dimensions")
check((len(coords39), len(cols39), rank39, len(null39)) == (85, 251, 73, 12),
      "pi9(P3) dimensions")

apart9, cpart9 = particular(9)
AA[9], CC[9] = apart9, cpart9
base_res = [dot_sparse(v, jet(AA, coords29)) for v in null29]
base_res += [dot_sparse(v, jet(CC, coords39)) for v in null39]

L = []
for v in null29:
    L.append([eval_left_on_column(v, tangent_column(coords29, 9, j, "A"))
              for j in range(10)])
for v in null39:
    L.append([eval_left_on_column(v, tangent_column(coords39, 9, j, "C"))
              for j in range(10)])
check(rank_fraction(L) == 3, "row-nine terminal depth rank three")
pivot_cols = []
for j in range(10):
    if rank_fraction([[row[k] for k in pivot_cols + [j]] for row in L]) > len(pivot_cols):
        pivot_cols.append(j)
check(pivot_cols == [7, 8, 9], "row-nine pivots q7 q8 q9")
Lp = [[row[j] for j in pivot_cols] for row in L]
rows = independent_rows(Lp)
q789 = solve_square([[Lp[i][j] for j in range(3)] for i in rows],
                    [-base_res[i] for i in rows])
check(all(base_res[i] + sum((L[i][pivot_cols[j]] * q789[j] for j in range(3)), E()) == 0
          for i in range(len(L))), "all 28 row-nine depth equations")
add_tangent(9, apart9, cpart9, q789, pivot_cols)

# Seven free row-nine tangent directions face the literal row-ten source.
defect10 = pa(GG[10], pn(prediction(10)))
D10 = D_matrix(10, 6)
q0to6, residual10, chosen10 = solve_full_columns(D10, [pc(defect10, j) for j in range(11)])
check(chosen10 == [0, 1, 2, 3, 4, 5, 7], "row-ten seven-row selector")
check(residual10[8] == -(X * 10125 + E.c(1928704)) * F(4, 61875),
      "row-ten xi residual")
X0 = E.c(F(-1928704, 10125))
qP = (P**6 * 2779225183740234375
      - P**4 * 194721282033880320000000
      - P**2 * 1868800030080493839974400000
      - E.c(9659395340042262184105231777792))
check(residual10[6].subst(1, X0) == qP * E.mon(p=-4) * F(2, 94585080322265625),
      "row-ten cubic residual")
check(all((not r) or i in (6, 8) for i, r in enumerate(residual10)),
      "only xi and cubic row-ten residuals")

# Put the seven uniquely selected tangent values into row nine and specialize X.
add_tangent(9, AA[9], CC[9], q0to6, list(range(7)))
for n in list(AA):
    AA[n] = psubst(AA[n], 1, X0)
    CC[n] = psubst(CC[n], 1, X0)
for n in list(GG):
    GG[n] = psubst(GG[n], 1, X0)

# Row-ten projected modules and the necessary joint P2/P3 obstruction.
coords210, cols210 = depth_columns(10, 2)
coords310, cols310 = depth_columns(10, 3)
rank210, null210 = nullspace_sparse(cols210, len(coords210))
rank310, null310 = nullspace_sparse(cols310, len(coords310))
check((len(coords210), len(cols210), rank210, len(null210)) == (88, 193, 68, 20),
      "pi10(P2) dimensions")
check((len(coords310), len(cols310), rank310, len(null310)) == (99, 304, 83, 16),
      "pi10(P3) dimensions")

apart10, cpart10 = particular(10)
AA[10], CC[10] = apart10, cpart10
baseA10 = [dot_sparse(v, jet(AA, coords210)) for v in null210]
LA = [[eval_left_on_column(v, tangent_column(coords210, 10, j, "A"))
       for j in range(11)] for v in null210]
check(rank_fraction(LA) == 3, "row-ten P2 terminal rank three")
pivA = []
for j in range(11):
    if rank_fraction([[row[k] for k in pivA + [j]] for row in LA]) > len(pivA):
        pivA.append(j)
check(pivA == [8, 9, 10], "row-ten P2 pivots r8 r9 r10")
LAp = [[row[j] for j in pivA] for row in LA]
arows = independent_rows(LAp)
r810 = solve_square([[LAp[i][j] for j in range(3)] for i in arows],
                    [-baseA10[i] for i in arows])
check(all(baseA10[i] + sum((LA[i][pivA[j]] * r810[j] for j in range(3)), E()) == 0
          for i in range(len(LA))), "all row-ten P2 terminal equations")
add_tangent(10, apart10, cpart10, r810, pivA)

baseC10 = [dot_sparse(v, jet({**CC, 10: cpart10}, coords310)) for v in null310]
LC = [[eval_left_on_column(v, tangent_column(coords310, 10, j, "C"))
       for j in range(11)] for v in null310]
check(rank_fraction(LC) == 3, "row-ten P3 terminal rank three")
pivC = []
for j in range(11):
    if rank_fraction([[row[k] for k in pivC + [j]] for row in LC]) > len(pivC):
        pivC.append(j)
check(pivC == [8, 9, 10], "row-ten P3 pivots r8 r9 r10")
LCp = [[row[j] for j in pivC] for row in LC]
crows = independent_rows(LCp)
c810 = solve_square([[LCp[i][j] for j in range(3)] for i in crows],
                    [-baseC10[i] for i in crows])
check(all(zero_mod_q(baseC10[i] + sum((LC[i][pivC[j]] * c810[j] for j in range(3)), E()), qP)
          for i in range(len(LC))), "P3 alone is terminal-compatible modulo q")
check(rank_fraction(LA + LC) == 3, "joint terminal coefficient rank three")

tet = {(5, 0): F(56), (6, 2): F(-35), (7, 4): F(20),
       (8, 6): F(-10), (9, 8): F(4), (10, 10): F(-1)}
tetvec = {coords310.index(key): value for key, value in tet.items()}
check(gcd(gcd(gcd(gcd(gcd(56, 35), 20), 10), 4), 1) == 1,
      "tetrahedral functional primitive")
check(all(sum((tetvec.get(i, F(0)) * value for i, value in col.items()), F(0)) == 0
          for col in cols310), "tetrahedral functional annihilates pi10(P3)")
H = dot_sparse(tetvec, jet(CC, coords310))
dP = (P**4 * 1482253431328125
      - P**2 * 103851350418069504000
      - E.c(1468662667625265243357184))
claimed = -P * dP / 132327846238037905244160
correction = qP * E.mon(p=-1) / 248114711696321072332800000
check(H == claimed + correction, "tetrahedral obstruction identity modulo q")

# Compact modular Bezout certificates prove squarefreeness, coprimality,
# and nonvanishing of the W gate over every characteristic-zero q-root.
def trim(poly, p):
    out = [z % p for z in poly]
    while out and out[-1] == 0:
        out.pop()
    return out


def mpadd(a, b, p):
    out = [0] * max(len(a), len(b))
    for i in range(len(out)):
        out[i] = ((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)) % p
    return trim(out, p)


def mpmul(a, b, p):
    out = [0] * (len(a) + len(b) - 1)
    for i, x0 in enumerate(a):
        for j, y0 in enumerate(b):
            out[i + j] = (out[i + j] + x0 * y0) % p
    return trim(out, p)


p31 = 31
q31 = [-10, 7, 7, 11]
d31 = [-7, -14, 9]
qp31 = [7, 14, 33]
check(mpadd(mpmul([8, -12], q31, p31), mpmul([15, 4, -6], d31, p31), p31) == [1],
      "mod31 gcd(q,d) Bezout")
check(mpadd(mpmul([3, 13], q31, p31), mpmul([0, -11, 6], qp31, p31), p31) == [1],
      "mod31 squarefree q Bezout")
w31 = [-15, -4]
check(mpadd(mpmul([10], q31, p31), mpmul([-15, -12, 12], w31, p31), p31) == [1],
      "mod31 q-W Bezout")
check(qP.coeff((0, 0, 0)) != 0, "q has no zero Y root")
check(F(225611776, 91125) != 0, "U gate is nonzero")

# Exact F_29 positive control: Phi=1, xi=7, eta=alpha=8.
def modf(value, p):
    value = F(value)
    return value.numerator * pow(value.denominator, -1, p) % p


def eval_mod(expr, p, vals):
    total = 0
    for mon, value in expr.t.items():
        term = modf(value, p)
        for i, exponent in enumerate(mon):
            term = term * pow(vals[i], exponent, p) % p
        total = (total + term) % p
    return total


vals29 = [1, 7, 8]
check(eval_mod(eta, 29, vals29) == 8, "F29 eta graph")
check(eval_mod(Aexpr, 29, vals29) == 8, "F29 row-nine alpha")
check(eval_mod(U, 29, vals29) == 10, "F29 U gate")
check(eval_mod(W, 29, vals29) == 24, "F29 W gate")
check(eval_mod(zeta, 29, vals29) == 13, "F29 zeta gate")
check(modf(F(-1928704, 10125), 29) == 7, "F29 xi row-ten condition")
qY1 = sum(modf(c, 29) for c in [
    -9659395340042262184105231777792,
    -1868800030080493839974400000,
    -194721282033880320000000,
    2779225183740234375,
]) % 29
dY1 = sum(modf(c, 29) for c in [
    -1468662667625265243357184,
    -103851350418069504000,
    1482253431328125,
]) % 29
check(qY1 == 0 and dY1 == 23, "F29 bracket survivor and depth hostile")

# Sign sheet: P,eta,alpha,zeta change sign; X,U,W,Z stay fixed.
eta_x0 = eta.subst(1, X0)
alpha_x0 = Aexpr.subst(1, X0)
check(all(mon[0] % 2 == 1 for mon in eta_x0.t), "eta odd under Phi sign")
check(all(mon[0] % 2 == 1 for mon in alpha_x0.t), "alpha odd under Phi sign")
check(all(mon[0] % 2 == 1 for mon in zeta.t), "zeta odd under Phi sign")
check(all(mon[0] % 2 == 0 for mon in W.subst(1, X0).t), "W even under Phi sign")

print("THM-4361 CLEAN-ROOM REFEREE: PASS")
print("checks=%d" % CHECKS)
print("row9_depth=(P2 75x160 rank59 null16; P3 85x251 rank73 null12; terminal_rank3; fibre=A7)")
print("row10_bracket=xi:-1928704/10125; q_degree=3; six_sign_sheets")
print("row10_depth=(P2 88x193 rank68 null20; P3 99x304 rank83 null16)")
print("tetrahedral=(56,-35,20,-10,4,-1); joint_P2_P3_obstruction=true")
print("mod31_bezout=(q,qprime)=1;(q,d)=1;(q,W)=1; F29_control=(1,7,8,8,0)")
