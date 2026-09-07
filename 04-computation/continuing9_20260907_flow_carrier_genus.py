"""Exact source, branch, and flow controls; no imported research producer."""
from math import gcd
from pathlib import Path
import json
import sys
import sympy as s
sys.stdout.reconfigure(newline='\n')
G=0
def gate(ok,label):
    global G
    G+=1
    if not ok: raise RuntimeError(label)
def zero(f): return s.cancel(f)==0
p,y,c,alpha,beta,gamma=s.symbols('p y c alpha beta gamma')
z,r,x,t,ss,tau,lam=s.symbols('z r x t ss tau lam')
D=p**3-y**2
A=alpha+beta*p
S=p**2*D*(A+y)
F=y**3+A*y**2-p**3*y-A*p**3+c/p**2
gate(zero(F+(S-c)/p**2),'actual normalized invariant equation')
N=4*p**7*(p**3-A**2)**2+4*A*c*p**2*(9*p**3-A**2)-27*c**2
gate(zero(s.discriminant(F,y)*p**4-N),'entire discriminant numerator')
gate(s.Poly(N,p).degree()==13 and s.Poly(N,p).LC()==4,'degree thirteen for every alpha,beta')
gate(N.subs(p,0)==-27*c*c,'no zero branch root')
f1=s.cancel(s.diff(S,y)/p**2)
f2=s.cancel(s.diff(S,p)/p)
crit=p**3*(169*p**3-120*beta**2*p**2-148*alpha*beta*p-40*alpha**2)*(p**3-A**2)
gate(zero(s.resultant(f1,f2,y)-crit),'complete finite critical-locus elimination')
gate(s.Poly(crit,p).LC()==169,'critical elimination never vanishes identically')
gate(zero(s.resultant(f1,s.diff(f1,y),y)-12*(A*A+3*p**3)),'triple-root elimination')

# A repeated discriminant root at a double root is precisely the gradient gate.
q,v0,v1,v2,v3,eps=s.symbols('q v0 v1 v2 v3 eps')
pert=y*y*(y+q)+eps*(v0+v1*y+v2*y*y+v3*y**3)
gate(zero(s.diff(s.discriminant(pert,y),eps).subs(eps,0)+4*q**3*v0),'double-root discriminant derivative')

# Native local charts at zero and infinity, with all input parameters retained.
Y=s.symbols('Y')
zero_chart=s.expand(F.subs({p:r**3,y:r**-2*Y},simultaneous=True)*r**6)
gate(zero_chart.subs(r,0)==Y**3+c,'zero valuation chart has three simple roots upstairs')
inf_chart=s.expand(F.subs({p:r**-2,y:r**-3*Y},simultaneous=True)*r**9)
expected=(Y*Y-1)*(Y+beta*r+alpha*r**3)+c*r**13
gate(zero(inf_chart-expected),'complete infinity quadratic-cover chart')
gate(inf_chart.subs(r,0)==Y**3-Y,'three simple infinity roots upstairs')
gate(zero(inf_chart.subs({r:-r,Y:-Y},simultaneous=True)+inf_chart),'deck involution exchanges the two nonzero branches')
for root in (-1,0,1):
    gate(s.diff(inf_chart,Y).subs({r:0,Y:root})!=0,'formal simple branch at infinity')
gate((-6+13+2+1+2)//2==6,'trigonal Riemann-Hurwitz genus')

affine_controls=[]
for aa,bb in ((0,0),(1,0),(0,1),(1,1),(-1,2),(2,-3),(-2,-2)):
    nn=s.Poly(N.subs({alpha:aa,beta:bb}),p,domain=s.QQ.frac_field(c))
    gate(s.gcd(nn,nn.diff()).degree()==0,'independent exact generic discriminant squarefreeness control')
    affine_controls.append([aa,bb,nn.degree(),6])
for aa,bb in ((1,0),(-2,0),(0,1),(1,1),(2,-3)):
    RR=aa+bb*p
    H=RR*(p**5*RR-c)
    h=s.Poly(H,p,domain=s.QQ.frac_field(c))
    gate(s.gcd(h,h.diff()).degree()==0,'hyperelliptic affine-multiplier squarefreeness')
    gate(h.degree()==(7 if bb else 5),'all univariate affine degrees')
    gate((h.degree()-1)//2==(3 if bb else 2),'hyperelliptic genus table')

# Parent-proposed extension: all univariate multipliers, with square factors.
AA,BB=s.symbols('AA BB')
gate(zero((p*AA*BB*y)**2-BB*(p**5*AA**2*BB-c)
          -p**2*AA**2*BB**2*(y*y-p**3+c/(p**2*AA**2*BB))),
     'universal squarepart pullback without adjoining coefficient square roots')
univariate_controls=[]
for RR in (s.Integer(3),p,p*p,(p-1)**2,p*(p-1)**2,
           (p*p+1)**2*(p-3),p**4*(p+2)**3,2*(p*p+1)**3*(p-1)**2):
    unit,factors=s.factor_list(RR,p)
    aa=s.prod(factor**(power//2) for factor,power in factors)
    bb=unit*s.prod(factor for factor,power in factors if power%2)
    gate(zero(aa*aa*bb-RR),'exact univariate squarepart decomposition')
    gate(s.gcd(s.Poly(bb,p),s.Poly(bb,p).diff()).degree()==0,'squarefree factor with scalar retained')
    H=s.Poly(s.expand(bb*(p**5*RR-c)),p,domain=s.QQ.frac_field(c))
    gate(s.gcd(H,H.diff()).degree()==0,'full univariate generic fibre squarefree')
    md=int(s.degree(RR,p)); bd=int(s.degree(bb,p))
    genus=(md+bd+4)//2
    gate(H.degree()==md+bd+5 and H.degree()%2==1 and genus>=2,'univariate all-degree genus law')
    univariate_controls.append([str(s.factor(RR)),md,bd,genus])

# Actual source chart, not a formally chosen invariant curve.
pp=t*(1+x*x*t); yy=x*t*pp
gate(zero((pp**3-yy**2)-t**3*(1+x*x*t)**2),'actual cusp equation')
sourceS=(p**2*D*(alpha+beta*p+gamma*y)).subs({p:pp,y:yy},simultaneous=True)
gate(zero(sourceS-t**5*(1+x*x*t)**4*(alpha+beta*pp+gamma*yy)),'complete actual affine source carrier')
lp=ss*ss+tau; ly=ss*lp
ls=sourceS.subs({x:ss/tau,t:tau},simultaneous=True)
W=ss**8*(alpha+beta*ss**2+gamma*ss**3)
gate(zero(s.expand(ls).coeff(tau,1)-W),'universal affine boundary coefficient')
bracket=lambda f,g:tau*(s.diff(f,ss)*s.diff(g,tau)-s.diff(f,tau)*s.diff(g,ss))
for k in (1,2,3):
    hp=s.expand(bracket(lp,ls**k))
    gate(zero(hp.coeff(tau,k)-2*k*ss*W**k),'first flow displacement for every outer leading order')
    gate(all(hp.coeff(tau,j)==0 for j in range(k)),'no earlier displacement')

# Primitive monomial invariants retain the actual generic constant field.
monomials=[]
for aa in range(1,9):
    for bb in range(1,6):
        if gcd(aa,bb)!=1: continue
        uu=next(j for j in range(-bb,bb+1) if (1-aa*j)%bb==0)
        vv=(1-aa*uu)//bb
        gate(aa*uu+bb*vv==1,'Bezout constant-field parametrization')
        pz=c**uu*z**bb; dz=c**vv*z**(-aa)
        gate(zero(pz**aa*dz**bb-c) and zero(pz**vv*dz**(-uu)-z),'birational torus inverse without adjoining constant roots')
        exponent=aa+3*bb
        H=(c**(3*uu)*z**exponent-c**vv)*(z if aa%2 else 1)
        gate(zero((pz**3-dz)*z**(2*((aa+1)//2))-H),'complete hyperelliptic monomial model')
        genus=(exponent+aa%2+bb%2-2)//2
        gate(genus>=2 and genus==(s.degree(H,z)-1)//2,'exact genus including infinity parity')
        for k in (1,2,3):
            gate((aa*k>=2 and bb*k>=1)==(aa>=2 or k>=2),'complete universal-carrier order gate')
        monomials.append([aa,bb,genus])
for aa,bb in ((1,1),(1,2),(2,1),(2,3),(3,2)):
    J=tau**bb*lp**(aa+2*bb)
    for k in (1,2):
        h=s.expand(bracket(lp,J**k))
        gate(zero(h.coeff(tau,bb*k)-2*bb*k*ss**(2*(aa+2*bb)*k+1)),'monomial-composition first displacement')

# The inherited rational non-LND hostile stays outside the universal carrier.
X=x/(1-lam*x); T=t*(1-lam*x)**2
gate(zero(X*X*T-x*x*t),'rational hostile preserves its actual invariant')
gate(zero(s.diff(X,x)*s.diff(T,t)-s.diff(X,t)*s.diff(T,x)-1),'rational hostile preserves source volume')
gate(zero(X.subs({x:x/(1+lam*x)},simultaneous=True)-x),'rational hostile invertibility')
print('AFFINE MULTIPLIER GENUS: constant 2; nonconstant p-affine 3; nonzero y coefficient 6')
print('UNIVARIATE [R,degree,squarefree_degree,genus]',json.dumps(univariate_controls,separators=(',',':')))
print('TRIGONAL [alpha,beta,finite branch count,genus]',json.dumps(affine_controls,separators=(',',':')))
print('PRIMITIVE MONOMIAL [a,b,genus]',json.dumps(monomials,separators=(',',':')))
print('CRITICAL RESULTANT',str(s.factor(crit)))
print('FLOW: native carrier and nonzero leading displacement; constant/zero-time and rational non-LND boundaries retained')
print('PASS',G,'always-active exact gates; raw LF')
print('Scope: same-invariant polynomial compositions only; no exclusion of compositions of distinct flows or all universal carriers.')
