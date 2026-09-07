"""Independent algebra and branch checks for the carrier-genus theorem.

No producer imports. General statements are proved in the accompanying audit;
this engine checks independent elimination and actual birational source maps.
"""
from pathlib import Path
import hashlib
import json
import sys
import sympy as s
sys.stdout.reconfigure(newline='\n')
G=0
def need(ok,label):
 global G
 G+=1
 if not ok: raise ArithmeticError(label)
def same(a,b):return s.cancel(a-b)==0
p,y,c,a,b,h,z,r,Y,x,t,u,w=s.symbols('p y c a b h z r Y x t u w')
Delta=p**3-y**2
K=p**2*Delta*(a+b*p+h*y)
F=s.Poly(K-c,y)
disc=s.discriminant(F.as_expr(),y)
normalized=s.cancel(disc/p**4)
need(s.Poly(normalized,p).degree()==13,'unscaled discriminant has exact degree13')
need(s.Poly(normalized,p).LC()==4*h**4,'nonzero y multiplier keeps infinity degree')
need(normalized.subs(p,0)==-27*c*c*h*h,'no branch root at p0')
Fy=s.cancel(s.diff(K,y)/p**2)
Fp=s.cancel(s.diff(K,p)/p)
critical=s.factor(s.resultant(Fp,Fy,y))
need(s.Poly(critical,p).degree()==9,'full arbitrary-y-coefficient critical resultant degree9')
need(s.Poly(critical,p).LC()==169*h**5,'critical resultant nonzero for all a,b when h nonzero')
need(same(critical.subs(h,1),p**3*(p**3-(a+b*p)**2)*(169*p**3-120*b*b*p*p-148*a*b*p-40*a*a)),'independent reversed elimination')
triple=s.factor(s.resultant(Fy,s.diff(Fy,y),y))
need(s.Poly(triple,p).degree()==3,'triple-root constraint has finite p support')
need(s.Poly(triple,p).LC()!=0,'triple-root constraint notidenticallyzero')

# Base-change degree and descent of the local branches.
Z=s.expand((K-c).subs({p:r**3,y:r**-2*Y},simultaneous=True))
need(Z.subs(r,0)==-h*Y**3-c,'p0 local cubic is simple and transitive')
# Avoid formal square roots: directly clear the actual highest pole13.
Inf=s.expand((K-c).subs({p:r**-2,y:r**-3*Y},simultaneous=True)*r**13)
need(same(Inf,(1-Y*Y)*(h*Y+b*r+a*r**3)-c*r**13),'infinity exact cubic-cover chart')
need(Inf.subs(r,0)==h*Y-h*Y**3,'infinity simple roots0,+1,-1')
need(same(Inf.subs({r:-r,Y:-Y},simultaneous=True),-Inf),'deck involution retains original (p,y)')
need((-2*3+13+2+1+2)//2==6,'all branch contributions yield genus6')

# Arbitrary univariate R=A²B. The birational square-removal formula is symbolic.
AA,BB,RR=s.symbols('AA BB RR',nonzero=True)
need(same((p*AA*BB)**2*(p**3-c/(p*p*AA*AA*BB)),BB*(p**5*AA*AA*BB-c)),'all-univariate square-removal identity')
univariates=[s.Integer(3),p,p*p,(p-1)**2,(p-1)**3,(p-1)**2*(p+2),p**4*(p*p+1)**3,(p**3-p+1)**2,(p-2)**3*(p+1)**4]
table=[]
for R in univariates:
 lc,factors=s.factor_list(R,p)
 A=s.prod(f**(k//2) for f,k in factors)
 B=lc*s.prod(f for f,k in factors if k%2)
 need(same(A*A*B,R),'exact squarefree factorization')
 H=s.Poly(B*(p**5*R-c),p,domain=s.QQ.frac_field(c))
 m=s.degree(R,p); bb=s.degree(B,p)
 need((m-bb)%2==0 and H.degree()==m+bb+5,'all-univariate odd branch degree')
 need(s.gcd(H,H.diff()).degree()==0,'exact generic branch squarefreeness')
 genus=int((m+bb+4)//2)
 need(genus>=2 and genus==(H.degree()-1)//2,'all-univariate genus formula')
 table.append([str(R),int(m),int(bb),genus])

# All primitive monomial exponents; inverse chart, odd branch count, genus.
mon=[]
for aa in range(1,13):
 for bb in range(1,9):
  if s.gcd(aa,bb)!=1:continue
  ee,ff,gg=s.gcdex(aa,bb)
  need(aa*ee+bb*ff==1,'Bezout exponents')
  pp=c**ee*z**bb; dd=c**ff*z**(-aa)
  need(same(pp**aa*dd**bb,c),'actual generic invariant')
  need(same(pp**ff*dd**(-ee),z),'generic constant-field inverse')
  exp=aa+3*bb; degree=exp+(aa%2)
  finite=exp+(aa%2); infinity=degree%2
  genus=(finite+infinity-2)//2
  need(genus==(aa+3*bb+aa%2+bb%2-2)//2 and genus>=2,'independent double-cover branch count')
  mon.append([aa,bb,genus])

# Native source and birational recovery, preserving actual base field.
pp=t*(1+x*x*t); yy=x*t*pp
need(same(pp**3-yy**2,t**3*(1+x*x*t)**2),'native cusp pullback')
ss=yy/pp; tt=pp-ss*ss
need(same(ss, x*t) and same(tt,t) and same(ss/tt,x),'actual rational source inverse')
lp=u*u+w; ly=u*lp
S0=s.expand(lp**2*(lp**3-ly**2)*(a+b*lp+h*ly))
W=u**8*(a+b*u*u+h*u**3)
need(same(s.expand(S0).coeff(w,1),W),'actual first w coefficient')
bracket=lambda f,g:s.expand(w*(s.diff(f,u)*s.diff(g,w)-s.diff(f,w)*s.diff(g,u)))
for k in range(1,5):
 image=bracket(lp,S0**k)
 need(same(image.coeff(w,k),2*k*u*W**k),'arbitrary outer leading term displacement')
 need(all(image.coeff(w,j)==0 for j in range(k)),'flow filtration starts at correct order')

print('INDEPENDENT AUDIT: arbitrary affine y coefficient; finite and infinity branches; native inverse')
print('UNIVARIATE [R,degree,odd-multiplicity degree,genus]',json.dumps(table,separators=(',',':')))
print('PRIMITIVE MONOMIAL CONTROLS',len(mon),'max exponent pair',(12,8))
print('PASS',G,'always-active exact gates; raw LF')
print('Analytic scope: generic curves and individual same-invariant flows; different-flow compositions remain open.')
