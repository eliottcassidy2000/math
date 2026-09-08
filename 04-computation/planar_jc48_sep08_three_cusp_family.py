#!/usr/bin/env python3
"""Exact three-ordinary-cusp (4,6) family controls; topology has a separate proof.
No parameter census and no new braid-word inference occurs in this source.
"""
import json
import sympy as S
s,a,q,t,z,u,v=S.symbols('s a q t z u v')
GATES=0

def check(value,name):
 global GATES
 GATES+=1
 if not value:raise RuntimeError(name)
def eq(x,y,name):check(S.cancel(x-y)==0,name)
def monic(x):return S.Poly(x,t).monic().as_expr()
def lead(expr,k,value,name):
 numerator,denominator=S.cancel(expr).as_numer_denom();p=S.Poly(numerator,z)
 check(all(j[0]>=k for j,c in p.terms()),name+' lower powers')
 check(denominator.subs(z,0)!=0,name+' denominator unit')
 eq(p.coeff_monomial(z**k),value*denominator.subs(z,0),name+' leading coefficient')

A=(t*t-1)*(t-s)
U=t**4-4*s*t**3/3-2*t*t+4*s*t
V=t**6+6*(a-s)*t**5/5-S.Rational(3,2)*t**4+2*(-a*s*s-a+s)*t**3+6*a*s*s*t


def family_algebra():
 eq(S.diff(U,t),4*A,'full quartic derivative')
 eq(S.diff(V,t),6*A*(t*t+a*t+a*s),'full sextic derivative in the literal b0 chart')
 Vq=S.integrate(6*A*(t*t+a*t+q),t)
 eq(Vq+S.Rational(3,2)*(a*s-q)*U,V,'degree-preserving linear target shear')
 eq(S.discriminant(A,t),4*(s*s-1)**2,'only forbidden collisions of the three cusp parameters')
 Rminus=45*a*a*s*s-162*a*a*s-27*a*a+80*a*s**3-216*a*s*s-216*a*s+128*s**3-432*s*s
 Rplus=45*a*a*s*s+162*a*a*s-27*a*a+80*a*s**3+216*a*s*s-216*a*s-128*s**3-432*s*s
 Rthird=36*a*a*s*s-216*a*a+28*a*s**3-108*a*s+11*s**4-126*s*s+675
 for e,expected in [(-1,-4*(s+1)**2*Rminus/675),(1,-4*(s-1)**2*Rplus/675),(s,(s*s-1)**2*Rthird/675)]:
  du=S.cancel((U-U.subs(t,e))/(t-e)**2);dv=S.cancel((V-V.subs(t,e))/(t-e)**2)
  check(S.Poly(du,t).LC()==1 and S.degree(du,t)==2,'finite cusp-extra-preimage incidence')
  eq(S.resultant(du,dv,t),expected,'explicit no-extra-cusp-preimage polynomial')

 jet=S.diff(U,t,2)*S.diff(V,t,3)-S.diff(V,t,2)*S.diff(U,t,3)
 for e in [-1,1,s]:
  eq(jet.subs(t,e),48*S.diff(A,t).subs(t,e)**2*(2*e+a),'exact ordinary cusp jet')
 eq(U.subs(t,1),-1+8*s/3,'positive cusp U value')
 eq(U.subs(t,-1),-1-8*s/3,'negative cusp U value')
 eq(U.subs(t,s),2*s*s-s**4/3,'third cusp U value')
 eq(U.subs(t,s)-U.subs(t,1),-(s-1)**3*(s+3)/3,'harmless projection coincidence at s minus3')
 eq(U.subs(t,s)-U.subs(t,-1),-(s+1)**3*(s-3)/3,'harmless projection coincidence at s plus3')
 X=S.cancel(U.subs(t,1/z)/V.subs(t,1/z));Z=S.cancel(1/V.subs(t,1/z))
 lead(X,2,1,'infinity coordinate X')
 lead(Z,6,1,'contact with the marked infinity line')
 h7=4*(3*a+2*s)/5
 lead(Z-X**3,7,h7,'first odd infinity term')
 X9=X.subs(a,-2*s/3);Z9=Z.subs(a,-2*s/3);k4=3-4*s*s/3;h9=16*s*(s*s-9)/27
 lead(Z9-X9**3,8,k4,'removable even infinity term')
 lead(Z9-X9**3-k4*X9**4,9,h9,'next odd infinity term')
 ordinary=(s*s-1)*(a*a-4)*(a+2*s)
 for ss in [0,-3,3]:
  eq(ordinary.subs({s:ss,a:-2*S.Rational(ss,3)}),0,'higher infinity boundary is nonordinary')
 for m,N in [(7,4),(9,3)]:check(3+(m-1)//2+N==10,'complete rational sextic genus budget')
 eq(k4.subs(s,S.Rational(3,2)),0,'vanishing even coefficient retained')
 for poly in [U,V]:eq(poly.subs({s:0,a:0}).subs(t,-t),poly.subs({s:0,a:0}),'excluded common quadratic-cover control')
 eq(jet.subs({s:0,a:0,t:0}),0,'quadratic-cover fixed point has zero cusp jet')
 return dict(h7=str(h7),h9=str(h9),ordinary=str(S.factor(ordinary)))


def control(name,ss,aa,H,m,expected_constant):
 UU=S.expand(U.subs({s:ss,a:aa}));VV=S.expand(V.subs({s:ss,a:aa}));es=[S.Integer(-1),S.Integer(1),ss]
 F=S.resultant(UU-u,VV-v,t)
 check(S.Poly(F,v).LC()==1 and S.degree(F,v)==4,'named monic quartic')
 check(S.Poly(F,u,v).total_degree()==6,'named image total degree six')
 eq(F.subs({u:UU,v:VV}),0,'named literal image equation')
 eq(monic(S.gcd(S.diff(UU,t),S.diff(VV,t))),S.prod(t-e for e in es),'only prescribed cusp critical parameters')
 for e in es:
  eq(monic(S.gcd(UU-UU.subs(t,e),VV-VV.subs(t,e))),(t-e)**2,'each cusp has one image preimage')
  check((S.diff(UU,t,2)*S.diff(VV,t,3)-S.diff(VV,t,2)*S.diff(UU,t,3)).subs(t,e)!=0,'named ordinary cusp jet')
 cusp=S.prod(u-UU.subs(t,e) for e in es)
 disc=S.discriminant(F,v)
 eq(disc,expected_constant*cusp**3*H**2,'complete named discriminant with coprojections')
 check(S.discriminant(H,u)!=0,'all residual node projection values simple')
 check(S.resultant(H,cusp,u)!=0,'residual nodes avoid cusp projection values')
 check(S.resultant(H,S.discriminant(UU-u,t),u)!=0,'source parameters simple at residual node values')
 N=S.degree(H,u);check(3+(m-1)//2+N==10,'named genus inventory')
 return dict(name=name,s=str(ss),a=str(aa),infinity=m,nodes=int(N),distinct_cusp_U_values=len(set(UU.subs(t,e) for e in es)),discriminant=str(S.factor(disc)))


def main():
 formulas=family_algebra();rows=[]
 rows.append(control('s0 harmless double cusp projection',S.Integer(0),S.Integer(1),(9*u+5)*(625*u**3+699*u*u+663*u+289),7,-S.Rational(65536,244140625)))
 rows.append(control('s3 harmless double cusp projection',S.Integer(3),S.Integer(1),(9*u+17)*(625*u**3-10881*u*u+49707*u+171549),7,-S.Rational(5308416,244140625)))
 rows.append(control('zero intermediate even coefficient',S.Rational(3,2),S.Integer(-1),16*u**3+36*u*u+27*u-765,9,-S.Integer(1296)))
 rows.append(control('literal certified braid representative',S.Integer(2),-S.Rational(4,3),729*u**3+3483*u*u+5547*u-78359,9,-S.Rational(167772160000,282429536481)))
 print('FINITE-EXACT THREE-CUSP FAMILY ALGEBRA PASS; topology and marked braid transfer proved separately')
 print(json.dumps(dict(formulas=formulas,controls=rows),sort_keys=True,indent=2))
 print('Always-active gates:',GATES)
if __name__=='__main__':main()
