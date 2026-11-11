"""a9_rank_one.py - HYP-9668 (E) ('if H^2(C^2, d - dW ^) is one-dimensional, the kernel of its functional is a
Mathieu-Zhao subspace') is FALSE for the note's own example W = x^2 y.

A Mathieu-Zhao subspace M of a commutative algebra that contains 1 is the whole algebra (take a = 1: a^m = 1 lies in
M for all m, so b = b a^m lies in M for every b).  For W = x^2 y:

      1 = 2 (d/dy - W_y)(y) - (d/dx - W_x)(x)                    (exact identity, checked below with sympy)

so 1 lies in M = Im(d/dx - W_x) + Im(d/dy - W_y), while M has codimension one (the functional
x^k y^l -> l! [k = 2l + 1] vanishes on M and equals 1 on x).  Hence M is not a Mathieu-Zhao subspace.

The same happens for every W = x^p y^q with gcd(p, q) = 1 and (p, q) != (1, 1) (the class of 1 is zero by weight).
For Broughton's polynomial x + x^2 y the class of 1 is NOT zero (functional value -1), so (E) is not refuted there.
The same example is a counterexample to the Image Conjecture of Zhao (J. Algebra 324 (2010), Conjecture 1.3) in its
original generality (A = C, z = (x, y), C = {d/dx - 2xy, d/dy - x^2}).
"""
import sympy as sp
from math import factorial
from tw import RingA, add, mul, scal, power, build, codim, P1, P2

x, y = sp.symbols("x y")
W = x ** 2 * y
Dx = lambda h: sp.expand(sp.diff(h, x) - sp.diff(W, x) * h)
Dy = lambda h: sp.expand(sp.diff(h, y) - sp.diff(W, y) * h)
print("identity 1 = 2 D_y(y) - D_x(x) for W = x^2 y:", sp.expand(2 * Dy(y) - Dx(x)) == 1)
lam = lambda k, l: factorial(l) if k == 2 * l + 1 else 0


def lam_poly(h):
    return sum(c * lam(*m) for m, c in sp.Poly(h, x, y).terms())


bad = 0
for k in range(0, 12):
    for l in range(0, 8):
        bad += (lam_poly(Dx(x ** k * y ** l)) != 0) + (lam_poly(Dy(x ** k * y ** l)) != 0)
print("functional x^k y^l -> l! [k = 2l+1] kills D_x, D_y of all monomials with k < 12, l < 8:", bad == 0, "; value on x:", lam(1, 0), "; value on 1:", lam(0, 0))

# machine confirmation with the independent engine: normal forms of 1 and of x
A = RingA()
X, Y, one = {(1, 0): 1}, {(0, 1): 1}, {(0, 0): 1}
pw = lambda a, k: power(A, a, k, one)
for name, Wd in (("x^2 y", mul(A, pw(X, 2), Y)), ("x^3 y^2", mul(A, pw(X, 3), pw(Y, 2))), ("x y^2", mul(A, X, pw(Y, 2))),
                 ("x y  (Gaussian)", mul(A, X, Y)), ("x + x^2 y  (Broughton)", add(X, mul(A, pw(X, 2), Y)))):
    for p in (P1, P2):
        E = build(A, Wd, 20, p)
        print("W = %-24s prime %d: dim = %d; class of 1: %-7s non-pivot monomials of degree <= 8: %s" % (
            name, p, codim(A, E, 8), "ZERO" if not E.nf(one) else "nonzero", E.nonpivots_upto(8)))
