"""Independent exact identities for the complete near-minimizer classification.

No producer imports. Universal compactness, uniform errors and analytic-product
arguments are audited in the companion report; finite checks do not replace them.
"""
import sympy as S
from fractions import Fraction as F
from math import factorial,comb
COUNT=0
def gate(ok,label):
 global COUNT
 COUNT+=1
 if not ok:raise RuntimeError(label)
def zero(x,label):gate(S.simplify(x)==0,label)
u=S.sqrt(2);v=S.sqrt(3);z=1/v
cstar=(13-8*u)/3
K3=4*v/(3*(1+u)*(1+v))
Kone=u-S.Rational(2,3)
Ktwo=(64-44*u)/3
Kdust=(28*u-32)/3
zero((1-cstar)/(2-u)-Kone,'one-atom quotient limit')
gate(F(7,6)**2<2,'one-atom constant exceeds one half')
gate(F(125,88)**2>2,'two-atom constant exceeds one half')
gate(F(4,3)**2<2,'dust constant exceeds two-atom constant')
gate(F(4,9)<F(1,2),'inherited K3 upper bound leaves strict boundary gap')
zero(((3-4*v/3)-cstar)/(2-2*u/v)-K3,'three-atom quotient limit')

# Independent directional differentiation, uniform over q+v0 scale directions.
eta,A,B=S.symbols('eta A B',real=True)
q=eta*A;w2=eta*B;t=S.sqrt(2-2*q-w2)
p3top=t*(1-q+w2)/2
p4top=((1-q)**2+(2-2*q-w2)*w2)/2
g=u-1-p3top+(2-u)*p4top
zero(g.subs(eta,0),'two-atom objective zero')
linear=S.diff(g,eta).subs(eta,0)
zero(linear-((7*u/4-2)*A+(2-11*u/8)*B),'uniform directional linear term')
zero(S.diff(2-u*t,eta).subs(eta,0)-(A+B/2),'distance linear term')
zero((1-p4top).subs(eta,0)/2-S.Rational(1,4),'two-atom energy')
zero(S.Rational(16,3)*(7*u/4-2)-Kdust,'dust slope')
zero(S.Rational(32,3)*(2-11*u/8)-Ktwo,'imbalance slope')
zero(S.Rational(16,3)*linear-(Kdust*A+Ktwo*B/2),'complete local slope mixture')

# Recover the original one-atom quotient without assuming any tail l1 bound.
Q=S.symbols('Q',nonnegative=True)
aa=S.sqrt(1-Q)
num=5-8*aa**3+3*aa**4
den=3*(1-aa**4)
zero(num.subs(Q,0),'one-atom numerator zero')
zero(S.diff(num,Q).subs(Q,0)-6,'one-atom numerator first term')
zero(S.diff(den,Q).subs(Q,0)-6,'one-atom denominator first term')

# Exact secant-defect terms, including their sign distinction.
r,b,C,w=S.symbols('r b C w')
f=lambda x:x-C*x*x
zero(r*r*(f(b)-f(r))-r*r*(b-r)*(1-C*(b+r)),'positive-root secant loss')
zero(w*w*(f(b)-f(-w))-w*w*f(b)-(w**3+C*w**4),'negative-root secant loss')
C3=(3-v)/2
zero(f(b).subs({b:z,C:C3})-(v-1)/2,'limiting positive square-mass tax')
zero((1-2*C*b).subs({b:z,C:C3})-(2-v),'limiting derivative tax')

# Integer multiplicity is retained in the quantitative square-measure bridge.
gate(F(1,3)/64==F(1,192),'outside-band square-distance gap')
gate(4*F(49,192)==F(49,48)>1,'at most three band roots')
gate(2*F(81,192)==F(27,32),'at most two band square masses')
gate(192*F(5,6144)==F(5,32),'strict small-M threshold')
gate(F(27,32)+F(5,32)==1,'at least three roots required')
gate(F(192)+F(192,49)==F(9600,49),'three-atom distance bound')
zero(z*z-S.Rational(1,3),'moment-measure normalization')

# Entire limit coefficients, reconstructed directly rather than through moments.
lam=1-v
co=[]
for n in range(9):
 co.append(S.simplify(sum(S.binomial(3,j)*z**j*lam**(n-j)/S.factorial(n-j) for j in range(min(3,n)+1))))
mom=[S.Integer(0)]
for n in range(1,9):
 pn=(-1)**(n-1)*n*co[n]
 for j in range(1,n):pn+=(-1)**(j+1)*co[j]*mom[n-j]
 pn=S.simplify(pn);mom.append(pn)
 zero(pn-(1 if n==1 else 3*z**n),'entire limit Newton moment '+str(n))
zero(co[2],'second elementary coefficient remains zero')
zero(mom[4]-2*z*mom[3]+S.Rational(1,3),'moment localization limit')

# Independent algebraic normalization of the full mixed-dust family.
L,c,W=S.symbols('L c W')
T=3+c
Delta=3*L**2+L*(T*T-3+c*c)-c*c
scale=(W-T)/(L-1)
d=T-scale
normdef=S.together(scale**2-3-(c*c+d*d)/L)
nn=S.Poly(S.together(normdef).as_numer_denom()[0],W)
zero(S.rem(nn,S.Poly(W*W-Delta,W)).as_expr(),'all-parameter normalization identity')
zero(Delta-T*T-(L-1)*(3*L+T*T+c*c),'strict positive scale')
zero(L*(T*T-3)-c*c-((L-1)*c*c+6*L*c+6*L),'negative-dust numerator positive')
zero(scale+d-T,'exact total first moment')
zero(c/scale-d/scale-(1-3/scale),'signed net dust identity')
zero(c*z-((c+3)*z-1)-(1-v),'limiting net dust independent of cancelling mass')
zero((4*z-1)-(v-1)-z,'c=1 hostile excess negative mass')

# Finite controls of the independent radical identity, no approximate roots.
for ll,cc in ((2,0),(3,1),(8,2),(16,2),(81,3),(256,4),(625,5)):
 DD=S.Integer(Delta.subs({L:ll,c:cc}));TT=3+cc
 gate(DD>TT*TT and ll>=2 and cc>=0,'eligible exact mixed-dust parameters')
 ss=(S.sqrt(DD)-TT)/(ll-1)
 dd=TT-ss
 zero(ss**2-3-(cc*cc+dd*dd)/ll,'literal normalized square sum')
 zero((3+cc-dd)/ss-1,'literal normalized signed first sum')
 gate(ll*(TT*TT-3)-cc*cc>0,'literal negative dust present')
print('Independent analytic audit: full five-way near-minimizer iff, both quotient degeneracies, and signed-net dust only.')
print('Exact directional expansions, secant losses, three-atom localization and entire-product Newton coefficients PASS.')
print('Exact arbitrary-cancelling-dust normalization and seven literal radical controls PASS.')
print('PASS',COUNT,'always-active exact gates; no producer or repository mathematical imports.')
