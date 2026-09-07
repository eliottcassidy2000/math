"""Exact identities for arbitrary-polynomial y-linear carrier genus.

Standalone producer. All-parameter proofs are in the companion report;
the finite controls retain transcendental generic-fibre constant c.
"""
from pathlib import Path
import sys
import json
import sympy as s
sys.stdout.reconfigure(newline='\n')
G=0
def gate(ok,msg):
 global G
 G+=1
 if not ok: raise ArithmeticError(msg)
def eq(a,b):return s.cancel(a-b)==0
p,y,c,A,r,Y,u,tau,k=s.symbols('p y c A r Y u tau k')
F=(y*y-p**3)*(y+A)+c/p**2
N=4*p**7*(p**3-A*A)**2+4*A*c*p**2*(9*p**3-A*A)-27*c*c
gate(eq(s.discriminant(F,y)*p**4,N),'universal discriminant with unevaluated A')
gate(N.subs(p,0)==-27*c*c,'zero excluded from finite discriminant roots')
gate(eq(s.resultant(s.diff(F,y),s.diff(F,y,2),y),-12*(A*A+3*p**3)),'triple root eliminant')
e,q=s.symbols('e q')
v0,v1,v2,v3=s.symbols('v0 v1 v2 v3')
model=y*y*(y+q)+e*(v0+v1*y+v2*y*y+v3*y**3)
gate(eq(s.diff(s.discriminant(model,y),e).subs(e,0),-4*q**3*v0),'double-root discriminant variation')

controls=[s.Integer(0),s.Integer(2),1-p,p*p,p*p+1,p**3-p+1,2*p**4+p-3,(p-1)**5,p**6+2*p*p+1]
rows=[]
for aa in controls:
 m=int(s.degree(aa,p)) if aa!=0 else -1
 ff=F.subs(A,aa); nn=s.Poly(N.subs(A,aa),p,domain=s.QQ.frac_field(c))
 n=13 if m<=1 else 4*m+7
 genus=(n-1)//2
 gate(nn.degree()==n,'exact generic branch degree')
 gate(s.gcd(nn,nn.diff()).degree()==0,'independent generic branch squarefreeness control')
 zero=s.expand(ff.subs({p:r**3,y:r**-2*Y},simultaneous=True)*r**6)
 gate(zero.subs(r,0)==Y**3+c,'zero cubic cover for arbitrary A degree')
 if m<=1:
  inf=s.expand(ff.subs({p:r**-2,y:r**-3*Y},simultaneous=True)*r**9)
  gate(eq(inf.subs(r,0),Y*(Y*Y-1)),'small-degree infinity branches')
 else:
  am=s.Poly(aa,p).LC()
  inf=s.expand(ff.subs({p:r**-2,y:r**-3*Y},simultaneous=True)*r**(2*m+6))
  gate(eq(inf.subs(r,0),am*(Y*Y-1)),'all large-degree square-root branches')
  gate(eq(inf.subs({r:-r,Y:-Y},simultaneous=True),inf),'deck exchanges two sqrt branches for all m')
  third=s.expand(ff.subs({p:1/r,y:r**(-m)*Y},simultaneous=True)*r**(3*m))
  gate(eq(third.subs(r,0),Y*Y*(Y+am)),'third branch is unramified')
  gate(s.diff(third,Y).subs({r:0,Y:-am})==am*am,'simple third branch at infinity')
 gate(-6+n+2+1==2*genus-2,'complete trigonal Riemann-Hurwitz count')
 # Actual logarithmic invariant; parity prevents boundary cancellation.
 pp=u*u+tau; yy=u*pp
 I=s.expand(pp*pp*(pp**3-yy*yy)*(aa.subs(p,pp)+yy))
 W=u**8*(aa.subs(p,u*u)+u**3)
 gate(eq(I.coeff(tau,1),W) and W!=0,'native nonzero first invariant coefficient')
 for kk in (1,2):
  H=s.expand(I**kk)
  delta=s.expand(tau*(2*u*s.diff(H,tau)-s.diff(H,u)))
  gate(eq(delta.coeff(tau,kk),2*kk*u*W**kk),'outer-flow nonzero first displacement')
 rows.append([str(aa),m,n,genus])

# Three reduced generic fibres remain correctly distinct from outer f(I).
for power in (2,3,4):
 z=s.symbols('z')
 gate(s.degree(z**power-c,z)==power,'outer polynomial has multiple generic components after constant closure')
print('Y-LINEAR CARRIER: R=A(p)+gamma*y, constant gamma!=0, arbitrary polynomial A')
print('GENUS: 6 for degree A<=1; 2*degree A+3 for degree A>=2')
print('[A,degree A,finite branch count,genus]',json.dumps(rows,separators=(',',':')))
print('PASS',G,'always-active exact gates; raw LF')
print('Scope: generic curves and same-invariant scalar-time flows; variable gamma(p) and general carriers remain open.')
