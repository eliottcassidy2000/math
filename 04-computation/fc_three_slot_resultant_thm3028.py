"""THREE-SLOT extension of THM-3022: the resultant invariant R_w(a,b,c).

At three slots f = c1 s^a + c2 s^b + c3 s^c, work projectively in P^2:
  L(f)=0    is a LINE,
  L(f^2)=0  is a CONIC,
  L(f^3)=0  is a CUBIC.
Line meets conic in 2 points (Bezout).  The three-slot threshold is 3 exactly
when NEITHER of those 2 points lies on the cubic -- one condition on the
WEIGHT, not on f.  So the natural invariant is the resultant of the binary
forms obtained by restricting L(f^2), L(f^3) to the line.

This generalises THM-3022: there the 'line' was a point and the invariant was
Q_w; here it is Res(deg 2, deg 3).
"""
import sympy as sp
from math import factorial, comb
from fractions import Fraction

u, v = sp.symbols('u v')

def Lpow(w, exps, coeffs, m):
    """L(f^m) for f = sum coeffs_j s^{exps_j}, weight w."""
    n=len(exps); tot=0
    def rec(j, rem, mono, coef):
        nonlocal tot
        if j==n-1:
            k=rem
            tot += coef*sp.binomial(m, k)*0 if False else 0
        return
    # direct multinomial expansion
    from itertools import product as iproduct
    def compositions(m, n):
        if n==1: yield (m,); return
        for i in range(m+1):
            for rest in compositions(m-i, n-1): yield (i,)+rest
    for k in compositions(m, n):
        mult = factorial(m)
        for kk in k: mult //= factorial(kk)
        term = mult
        for j in range(n): term *= coeffs[j]**k[j]
        e = sum(k[j]*exps[j] for j in range(n))
        tot += term*w(e)
    return sp.expand(tot)

def R_w(w, a, b, c):
    """Resultant invariant.  Parametrise the line L(f)=0 by (u,v)."""
    A,B,C = w(a), w(b), w(c)
    # c1 = u*B, c2 = -u*A + v*C, c3 = -v*B   satisfies c1 A + c2 B + c3 C = 0
    coeffs = [u*B, -u*A + v*C, -v*B]
    F2 = sp.Poly(Lpow(w, [a,b,c], coeffs, 2), u, v)
    F3 = sp.Poly(Lpow(w, [a,b,c], coeffs, 3), u, v)
    if F2.total_degree()==0 or F3.total_degree()==0:
        return sp.Integer(0) if (F2.total_degree()==0 and F2.as_expr()==0) else None
    return sp.resultant(F2.as_expr(), F3.as_expr(), v)

WS = {
 "factorial j!":      lambda j: sp.Integer(factorial(j)),
 "[0,1] Lebesgue 1/(j+1)": lambda j: sp.Rational(1, j+1),
 "central binom":     lambda j: sp.Integer(comb(2*j,j)),
 "geometric 3^j":     lambda j: sp.Integer(3**j),
 "Fibonacci":         lambda j: sp.Integer([0,1,1,2,3,5,8,13,21,34,55,89,144,233,377,610,987,1597,2584,4181,6765,10946,17711,28657][j]),
}
print("R_w(a,b,c): resultant of the restricted quadric and cubic on the line L(f)=0")
print("R_w != 0  =>  three moments force f = 0 (three-slot threshold 3)")
print()
for name,w in WS.items():
    print(f"  {name}:")
    for (a,b,c) in ((0,1,2),(0,1,3),(1,2,3),(0,2,4)):
        try:
            r = R_w(w,a,b,c)
            r = sp.simplify(r)
            print(f"     (a,b,c)=({a},{b},{c}):  R = {r}   -> {'threshold 3' if r!=0 else 'VANISHES'}")
        except Exception as e:
            print(f"     ({a},{b},{c}): error {type(e).__name__}")
