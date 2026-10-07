"""Audit A misc checks:
 (1) note §7.1: rational unit vectors (a/c, b/c) reduce mod 7 into the circle x^2+y^2=1 over F_7 and (1+2i)*circle = knight set.
 (2) THM-4559: separability-idempotent certificate for Palvolgyi's heptagon with ALGEBRAIC r = 3 (F = Q(sqrt5, sqrt2)),
     checked exactly in F (x) F = Q[a,b,A,B]/(a^2-5, b^2-2, A^2-5, B^2-2); and the same without e is nonzero."""
import itertools, random
from fractions import Fraction as Fr
from math import gcd
import sympy as sp

# (1)
knight = {(1, 2), (2, 1), (-1, 2), (-2, 1), (1, -2), (2, -1), (-1, -2), (-2, -1)}
knight = {(x % 7, y % 7) for x, y in knight}
circle = {(x, y) for x in range(7) for y in range(7) if (x * x + y * y) % 7 == 1}
img = set()
ok = True
for m in range(1, 60):
    for n in range(1, m):
        if gcd(m, n) != 1 or (m - n) % 2 == 0: continue
        a, b, c = m * m - n * n, 2 * m * n, m * m + n * n
        for sa, sb, sw in itertools.product((1, -1), (1, -1), (0, 1)):
            x, y = (sa * a, sb * b) if not sw else (sb * b, sa * a)
            ok &= (c % 7 != 0)
            ci = pow(c, -1, 7)
            r = ((x * ci) % 7, (y * ci) % 7)
            ok &= r in circle
            # multiply by 1+2i: (x + i y)(1 + 2i) = (x - 2y) + i(2x + y)
            img.add(((r[0] - 2 * r[1]) % 7, (2 * r[0] + r[1]) % 7))
for r in [(1, 0), (0, 1), (6, 0), (0, 6)]:
    img.add(((r[0] - 2 * r[1]) % 7, (2 * r[0] + r[1]) % 7))
print("(1) all primitive Pythagorean unit vectors (m<60): 7 not | c, residue on circle:", ok, "; |circle| =", len(circle),
      "; image under (1+2i) = knight set:", img == knight)

# (2)
a, b, A, B = sp.symbols('a b A B')
rel = {a**2: 5, b**2: 2, A**2: 5, B**2: 2}
def red(expr):
    expr = sp.expand(expr)
    # reduce powers: repeatedly substitute squares
    p = sp.Poly(expr, a, b, A, B)
    out = 0
    for (ea, eb, eA, eB), coef in p.terms():
        t = coef * 5 ** (ea // 2) * 2 ** (eb // 2) * 5 ** (eA // 2) * 2 ** (eB // 2)
        out += t * a ** (ea % 2) * b ** (eb % 2) * A ** (eA % 2) * B ** (eB % 2)
    return sp.expand(out)
# heptagon with r = 3: (0,0), (j, +-sqrt(j(6-j))) j = 1,3,4 -> sqrt5, 3, sqrt8 = 2 sqrt2
pts_left = [(0, 0), (1, a), (1, -a), (3, 3), (3, -3), (4, 2 * b), (4, -2 * b)]
pts_right = [(0, 0), (1, A), (1, -A), (3, 3), (3, -3), (4, 2 * B), (4, -2 * B)]
H = sp.Matrix([[0, -3, 0], [-3, 1, 0], [0, 0, 1]])       # x^2 + y^2 - 6x = 0
e = sp.Rational(1, 4) * (1 + a * A / 5) * (1 + b * B / 2)
assert red(e * e - e) == 0                                 # idempotent
# m(e) = 1: substitute A -> a, B -> b then reduce
assert red(red(e).subs({A: a, B: b})) == 1
ok_e, nonzero_without = True, 0
for (x1, y1), (x2, y2) in zip(pts_left, pts_right):
    pl = sp.Matrix([1, x1, y1]); pr = sp.Matrix([1, x2, y2])
    E0 = red((pl.T * H * pr)[0, 0])
    Ee = red(e * E0)
    ok_e &= (Ee == 0)
    nonzero_without += (E0 != 0)
    assert red(E0.subs({A: a, B: b})) == 0                 # multiplied version = sphere equation
print("(2) heptagon r=3: e idempotent, m(e)=1; tensor evaluations vanish with e:", ok_e,
      "; evaluations nonzero without e:", nonzero_without, "of 7")
