#!/usr/bin/env python3
"""Exact audit of the refined affine-tail complete-turn reduction.

This proves the 1<=T<4 bands for raw p>=264, the 4<=T<5 band for
p>=277, and the |c|<=13 reduction when T<1.  It intentionally leaves the
finite bridge p=264..276 and the T<1 affine tails open.
"""
from fractions import Fraction as F
import sympy as sp

DMAX=F(186_636_088_362,11_773_143_757_375);TARGET=DMAX/5

def require(x,msg):
 if not x:raise RuntimeError(msg)

def B(m,p,L):
 # Universal specialization d=8,r=7,z<=Lp,w>=L(p+1)-14.
 w=L*(p+1)-14
 return (F(m*(p-7),49*(m+1)*p)
         -F(8*L,14*w)-F((m+1)**2*64*L,2*w*(p-7)))

def derivative_audit(m,p0):
 p,L,x,y,q=sp.symbols('p L x y q',positive=True)
 w=L*(p+1)-14
 fun=(sp.Rational(m,49*(m+1))*(p-7)/p
      -8*L/(14*w)-sp.Rational((m+1)**2*64,2)*L/(w*(p-7)))
 # p monotonicity: positive coefficients after the hostile corner shift.
 np,dp=sp.fraction(sp.together(sp.diff(fun,p)))
 P=sp.Poly(sp.expand(np.subs({p:x+p0,L:y+168})),x,y)
 require(all(c>0 for c in P.coeffs()),('p derivative',m,P))
 # Both losses decrease as L increases, so B increases in L.  Certify by
 # the same positive-coefficient corner shift.
 nL,dL=sp.fraction(sp.together(sp.diff(fun,L)))
 PL=sp.Poly(sp.expand(nL.subs({p:x+p0,L:y+168})),x,y)
 require(all(c>0 for c in PL.coeffs()),('L derivative',m,PL))
 return len(P.terms()),len(PL.terms())

def main():
 corners=[]
 for m in range(1,5):
  p0=264 if m<4 else 277
  audit=derivative_audit(m,p0)
  value=B(m,p0,168)
  gap=value-TARGET
  require(gap>0,(m,value,gap))
  corners.append((m,p0,value,gap,audit))
 cap=F(8*264,264-7)+F(37,7)
 require(cap==F(24293,1799)<14,cap)
 print('LRC14 refined affine turn-band reduction exact audit')
 print('bands=',corners)
 print('Tlt1_c_strict_cap=',cap,';integer_abs_c_cap=13',sep='')
 print('status=PROVED reductions; OPEN p264..276 bridge and T<1 affine tails')

if __name__=='__main__':main()
