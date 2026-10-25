#!/usr/bin/env python3
"""Exact identity and hostile-control replay for the nonzero-defect (4,6) lane.

The all-degree conclusion is proved in anchors_quartic_sextic.md; this is
supplemental symbolic/finite verification, not an exhaustive Keller census.
Requires SymPy. All guards execute under python -O as well.
"""
from __future__ import annotations

import json
import sympy as s

CHECKS = 0


def check(condition, label):
    global CHECKS
    CHECKS += 1
    if not bool(condition):
        raise RuntimeError(label)


def zero(expr, label):
    check(s.cancel(expr) == 0, label)


def dinf(expr, x):
    num, den = s.fraction(s.cancel(expr))
    if num == 0:
        return None
    return s.degree(num, x) - s.degree(den, x)


x, y, z = s.symbols("x y z")
a, b, d, kap, E, H, L, M = s.symbols("a b d kappa E H L M")
da, db, dd = s.symbols("da db dd")
q = s.Rational


def coeffs(a, b, d, kap, E, H):
    h = q(3, 4)*a*b + q(3, 8)*kap*a + H
    f = q(3, 4)*a*d + q(3, 8)*b*b + q(3, 8)*kap*b - a**3/16 + E*a/2
    P = y**4 + a*y*y + b*y + d
    Q = y**6 + 3*a*y**4/2 + (3*b+kap)*y**3/2 + (3*d/2+3*a*a/8+E)*y*y+h*y+f
    I = 3*b*d + 3*kap*d/2 - 3*a*a*b/4 - 3*kap*a*a/16 + 2*E*b + H*a
    K = 3*d*d/2 - 3*a*a*d/4 + 2*E*d - 3*a*b*b/4 - 3*kap*a*b/8 + H*b + 3*a**4/32 - E*a*a/2
    return P, Q, h, f, I, K


P, Q, h, f, I, K = coeffs(a, b, d, kap, E, H)
D = lambda v: s.diff(v, a)*da + s.diff(v, b)*db + s.diff(v, d)*dd
J = s.Poly(s.expand(D(P)*s.diff(Q,y)-s.diff(P,y)*D(Q)), y)
for j in range(3, 9):
    zero(J.nth(j), f"positive top row {j}")
zero(J.nth(2)-D(I), "y2 derivative")
zero(J.nth(1)-D(K), "y1 derivative")
zero(J.nth(0)-(h*dd-b*D(f)), "actual scalar J")

# Independently retain the constant target shear and its top-row effect.
C, F = s.symbols("C F")
Q_unreduced = Q+C*P+F
zero(D(P)*s.diff(Q_unreduced,y)-s.diff(P,y)*D(Q_unreduced)-J.as_expr(), "constant target shear")

r, U = 2*b+kap, L+E*kap
N = 9*a*(kap*kap*a**3-16*b*r**3)+96*kap*a*a*(H*a-U)+256*(H*a-U)**2+r*r*(384*H*b-256*E*E-384*M)
res = s.resultant(I-L, K-M, d)
zero(res-q(3,512)*N, "independent Sylvester resultant")
d_from_I = a*a*(4*b+kap)/(8*r)-2*E/3-2*H*a/(3*r)+2*U/(3*r)
zero((I-L).subs(d,d_from_I), "exact d chart")
zero(s.cancel((K-M).subs(d,d_from_I))-N/(384*r*r), "K chart normalization")

# The finite control family, not a formal Keller source.
b0 = kap*kap*z**3/16
a0 = z*(2*b0+kap)
d0 = z*z*(2*b0+kap)*(4*b0+kap)/8
P0,Q0,h0,f0,I0,K0 = coeffs(a0,b0,d0,kap,s.Integer(0),s.Integer(0))
zero(I0, "comparison I")
zero(K0, "comparison K")
J0 = s.diff(P0,z)*s.diff(Q0,y)-s.diff(P0,y)*s.diff(Q0,z)
J0_formula = 3*kap**4*z*z*(kap*z**3+8)*(10*kap*kap*z**6+61*kap*z**3+64)/16384
zero(J0-J0_formula, "full two-variable comparison determinant")
check(s.degree(J0,z)==11, "generic comparison degree")
zero(s.Poly(J0,z).LC()-q(15,8192)*kap**7, "generic leading coefficient")

# Polynomial controls include multiple degrees and signs of the nonzero defect.
polynomial_rows=[]
for kval in (s.Integer(1),s.Integer(2),s.Integer(-3)):
    for zz in (x,x*x+x+1,x**3-2*x):
        m=int(s.degree(zz,x))
        pp,qq=[s.expand(v.subs({z:zz,kap:kval})) for v in (P0,Q0)]
        jj=s.expand(s.diff(pp,x)*s.diff(qq,y)-s.diff(pp,y)*s.diff(qq,x))
        expected=s.expand(J0_formula.subs({z:zz,kap:kval})*s.diff(zz,x))
        zero(jj-expected,"polynomial direct replay")
        check(s.degree(jj,y)==0,"polynomial y-independent")
        check(s.degree(jj,x)==12*m-1,"sharp polynomial x-degree")
        zero(s.Poly(jj,x).LC()-q(15,8192)*kval**7*m*s.Poly(zz,x).LC()**12,"sharp polynomial leading term")
        check(s.degree(jj,x)>0,"first integrals do not imply constant Jacobian")
        polynomial_rows.append({"kappa":int(kval),"z":str(zz),"degree_J":int(s.degree(jj,x))})

# Rational charts with finite poles cannot be silently called polynomial changes.
rational_rows=[]
for m in (1,2,3):
    zz=(x**(m+1)+1)/(x-1)
    pp,qq=[v.subs({z:zz,kap:1}) for v in (P0,Q0)]
    jj=s.cancel(s.diff(pp,x)*s.diff(qq,y)-s.diff(pp,y)*s.diff(qq,x))
    expected=s.cancel(J0_formula.subs({z:zz,kap:1})*s.diff(zz,x))
    zero(jj-expected,"rational direct replay")
    check(dinf(zz,x)==m,"rational z degree infinity")
    check(dinf(jj,x)==12*m-1,"rational scalar degree infinity")
    check(s.denom(s.cancel(zz)).subs(x,1)==0,"finite pole hostile retained")
    rational_rows.append({"m":m,"degree_infinity_J":int(dinf(jj,x)),"finite_pole":1})

# Constant-a/nonconstant-b branch; a=H=0 and d=-2E/3 are retained exactly.
p1,q1,_,_,i1,k1=coeffs(s.Integer(0),x,-2*E/3,kap,E,s.Integer(0))
j1=s.expand(s.diff(p1,x)*s.diff(q1,y)-s.diff(p1,y)*s.diff(q1,x))
zero(j1+q(3,8)*x*(2*x+kap),"constant-a hostile branch")
zero(s.diff(i1,x),"constant-a first integral I")
zero(s.diff(k1,x),"constant-a first integral K")
check(s.degree(j1,x)==2,"constant-a scalar degree")

# Degree comparison is checked as an integer inequality for a closed positive
# degree universe; the proof gives its unrestricted two-term argument.
for A in range(1,25):
    for B in range(1,25):
        top=max(4*A,A+4*B)
        other=max(3*A,2*A,A,3*B,2*B,B,0)
        check(other<top,"two-term degree dominance")
        if 4*A==A+4*B:
            check(A%4==0 and B%3==0 and A//4==B//3,"forced degree ratio")

receipt={
    "status":"PASS",
    "scope":"Exact symbolic identities plus a closed finite control universe; all-degree proof is in the companion note.",
    "checks":CHECKS,
    "polynomial_controls":polynomial_rows,
    "rational_pole_controls":rational_rows,
    "degree_universe":{"A":[1,24],"B":[1,24]},
    "identity_J0":str(s.factor(J0_formula)),
    "normalization":"Res_d(I-L,K-M)=3N/512; K(d_from_I)-M=N/(384r^2)",
    "open":"Other genus-six terminal degree pairs, planar JC, and Collatz remain open."
}
print(json.dumps(receipt,indent=2,sort_keys=True))
