"""Independent rational-polynomial audit; no producer or symbolic-library imports."""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import hashlib
import sys
sys.stdout.reconfigure(newline='\n')
COUNT=0

def check(ok,msg):
    global COUNT
    COUNT+=1
    if not ok:raise ArithmeticError(msg)

class Poly:
    # Laurent polynomials over Q in (z,e3,e4,e5), immutable-by-convention.
    def __init__(self,p=0):
        if isinstance(p,Poly):self.d=dict(p.d)
        elif isinstance(p,dict):self.d={k:F(v) for k,v in p.items() if v}
        else:self.d={} if not p else {(0,0,0,0):F(p)}
    def __add__(self,o):
        d=dict(self.d)
        for k,v in Poly(o).d.items():d[k]=d.get(k,F(0))+v
        return Poly(d)
    __radd__=__add__
    def __neg__(self):return Poly({k:-v for k,v in self.d.items()})
    def __sub__(self,o):return self+-Poly(o)
    def __rsub__(self,o):return Poly(o)+-self
    def __mul__(self,o):
        d={}
        for a,u in self.d.items():
            for b,v in Poly(o).d.items():
                k=tuple(x+y for x,y in zip(a,b));d[k]=d.get(k,F(0))+u*v
        return Poly(d)
    __rmul__=__mul__
    def __pow__(self,n):
        if n<0:
            if len(self.d)!=1:raise ValueError('negative power only for monomial')
            k,v=next(iter(self.d.items()));return Poly({tuple(n*x for x in k):v**n})
        q=Poly(1)
        for _ in range(n):q=q*self
        return q
    def __eq__(self,o):return self.d==Poly(o).d
    def coeff(self,n):return Poly({(0,)+k[1:]:v for k,v in self.d.items() if k[0]==n})
    def eval(self,x):return sum(v*__import__('math').prod(t**a for t,a in zip(x,k)) for k,v in self.d.items())

def variable(i):return Poly({tuple(int(j==i) for j in range(4)):1})
z,e3,e4,e5=[variable(i) for i in range(4)]
B=z**5-13*z**4+55*z**3-e3*z*z+e4*z-e5
C=z**4-12*z**3+45*z*z-F(2,3)*e3*z+F(3,7)*e4
D=z**4-11*z**3+36*z*z-F(5,12)*e3*z+F(1,7)*e4
check(z*C-F(2,3)*B==F(1,3)*z**3*(z-5)**2-F(5,21)*e4*z+F(2,3)*e5,'formal C elimination')
check(z*D-F(5,12)*B==F(1,12)*z**3*(7*z*z-67*z+157)-F(23,84)*e4*z+F(5,12)*e5,'formal D elimination')
# Separate use of variables for remaining-root quadratic: z=t,e3=delta,e4=ab.
t=z;delta=e3;ab=e4;d=5+delta
sumce=13-t-d
prodce=55-13*t+t*t-ab-d*sumce
check(25-5*sumce+prodce==-3*t+2*delta+t*t+delta*t+delta*delta-ab,'formal remaining roots identity')
check(7*(5+z)**2-67*(5+z)+157==-3+3*z+7*z*z,'formal endpoint quadratic')

# Every bound needed for the universal root-location contradiction.
check(F(118)<121 and 81>59,'third/fourth-root estimates')
eps=F(1,100);tb=3*eps;db=F(1,10)
check(F(5,7)*eps<db*db,'fourth-root closeness')
pairbound=3*tb+2*db+tb*tb+db*tb+db*db+tb*tb/4
check(pairbound==F(2433,8000)<F(9,25),'uniform pair product bound')
check(F(9,25)/F(9,5)==F(1,5),'largest root location')
fbound=F(3,2)*eps*eps
quad=-3+3*F(1,5)+7*F(1,5)**2
check(quad==-F(53,25),'endpoint quadratic sign')
Dupper=F(24,5)**2*quad/12+5*fbound/(12*F(24,5))
check(Dupper==-F(2544,625)+F(1,76800)<0,'weak-largest-root sign contradiction')

# Independently build BOTH ordinary E/O factors; derive, do not assume, their
# Vandermonde multiplier and the exact Laurent support at every exponent.
O=sum((comb(14,2*j+1)*z**j for j in range(7)),Poly())
E=sum((comb(14,2*j)*z**j for j in range(8)),Poly())
mult=O*O+z**(-1)*E*E
for j in range(-2,15):
    expected=comb(28,2*j+2) if 0<=2*j+2<=28 else 0
    check(mult.coeff(j)==expected,'complete E/O multiplier coefficient')
GB=1+13*z+55*z*z+e3*z**3+e4*z**4+e5*z**5
GC=1+12*z+45*z*z+F(2,3)*e3*z**3+F(3,7)*e4*z**4
GD=1+11*z+36*z*z+F(5,12)*e3*z**3+F(1,7)*e4*z**4
beta=z**(-1)*GB;cr=z**(-1)*GC;dr=z**(-1)*GD
mix=beta*beta+2*z*cr*dr
# Coefficientwise product in z only, retaining formal coefficient variables.
Q=sum((mult.coeff(j)*mix.coeff(j)*z**j for j in range(-2,15)),Poly())
P=sum((O.coeff(j)*beta.coeff(j)*z**j for j in range(-1,8)),Poly())
check(P==182+20020*z+2002*e3*z*z+3432*e4*z**3+2002*e5*z**4,'original first-row compiler')
check(Q.coeff(8)==comb(28,18)*e5*e5,'top carried coefficient')
check(Q.coeff(7)==comb(28,16)*(2*e4*e5+F(6,49)*e4*e4),'next carried coefficient')
check(Q.coeff(-1)==28 and mix.coeff(-2)==1,'lower carry retained, erased coefficient accounted')
check(min(k[0] for k in Q.d)==-1 and max(k[0] for k in Q.d)==8,'complete Laurent response range')
for k,v in Q.d.items():check(v>=0,'coefficient polynomial nonnegative')
for cc in (GC,GD):
    check(all(v>=0 for v in (GB-cc).d.values()),'formal contiguous coefficient domination')
check(all(comb(28,j)<=comb(28,14) for j in range(29)),'binomial mass maximum')
# Univariate coefficient mass after substituting z=1 preserves formal e variables.
getmass=lambda poly:sum((poly.coeff(j) for j in range(-2,10)),Poly())
check(getmass(mix)==getmass(GB)**2+2*getmass(GC)*getmass(GD),'full untruncated coefficient mass')
# Reconstruct P(-s)/(2002s^3) as a formal Laurent relation.
Pminus=sum(((-1)**j*P.coeff(j)*z**j for j in range(5)),Poly())
check(F(1,2002)*z**(-3)*Pminus==e5*z-F(12,7)*e4+e3*z**(-1)-10*z**(-2)+F(1,11)*z**(-3),'original-root relation')
check(F(12,7)*comb(28,18)-2*comb(28,16)==-38346750,'leading cross cancellation')
lowmass=3*comb(28,14)*F(18,5)**10
fifthmass=10*comb(28,18)*F(13,5)**5
H=lowmass+fifthmass;floor=F(6,49)*comb(28,16)*eps**2
check(H==F(17194291313831660508,390625) and floor==F(2607579,7000),'exact tail constants')
S=118163898523
check(H/floor<S,'stated integer strictly above displayed coarse ratio')
check(-floor+H/S<0,'strict tail at stated endpoint')
for sv in (F(1),F(3,2),F(2),F(S),F(S+1)):
    for j in range(-1,7):check(sv**(j-7)<=1/sv,'all lower Laurent terms bounded only for s>=1')
check(F(1,2)**(-8)>1/F(1,2),'s>=1 restriction has a literal hostile')

# Nonvacuity from a direct product and exact ordered signs, not a matrix oracle.
roots=tuple(F(n,5) for n in (1,3,9,22,30))
poly=Poly(1)
for root in roots:poly=poly*(z-root)
params=tuple(poly.coeff(j).eval((F(1),F(0),F(0),F(0))) for j in (2,1,0))
u3,u4,u5=-params[0],params[1],-params[2]
check((u3,u4,u5)==(F(2127,25),F(27144,625),F(3564,625)),'literal model coefficients')
check(sum(roots)==13 and sum(x*x for x in roots)==59,'literal anchors')
for pp in (C,D):
    for i,root in enumerate(roots):check((-1)**i*pp.eval((root,u3,u4,u5))>0,'strict ordered interlacing')
check(u4>eps,'strict model floor control')
for root in (F(0),F(3),F(5)):
    check(B.eval((root,F(75),F(0),F(0)))==0,'rejected boundary roots')
check(D.eval((F(5),F(75),F(0),F(0)))==-F(25,4),'second interlacer needed at boundary')

print('PASS: explicit e4>1/100 within both-anchor, both-weak-interlacer model.')
print('Independent raw E/O convolution and Laurent support: -1 through8; q[-1]=28.')
print('D endpoint upper bound',Dupper)
print('Tail ratio',H/floor,'; strict integer endpoint',S)
print('No finite-phase or general actual Laurent noncancellation assertion.')
print('PASS',COUNT,'always-active independent exact gates')
