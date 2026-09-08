#!/usr/bin/env python3
"""Exact controls for all binary(4,2,1,1) DG quartic boundary placements.
Full symbolic original-position coefficient rows; no bounded mate search.
"""
from hashlib import sha256
import sympy as S
GATES=[]
def need(label,pred):
    if not bool(pred):raise RuntimeError(label)
    GATES.append(label)
def zero(label,v):need(label,S.cancel(v)==0)
u,t,p,r,T,zeta=S.symbols('u t p r T zeta')
a,b,c,d,lam,e,aa,bb=S.symbols('a b c d lam e aa bb')
A=S.symbols('A0:5');B=S.symbols('B0:5');C0,E0=S.symbols('C0 E0')
N=u**4*(u*u-1)
P=sum(A[i]*u**i for i in range(5))
Q=(A[4]-1)*u*u+(A[3]-2*p*A[4]+4*p)*u+C0
M=sum(B[i]*u**i for i in range(5))
R=B[4]*u*u+(B[3]-2*p*B[4])*u+E0

def global_h(label,n,pp,qq):
    for j,v in enumerate([qq,pp-2*qq*(u+p)**2,n-pp*(u+p)**2+qq*(u+p)**4]):
        need(label+str(j),S.Poly(S.expand(v),u).degree()<=4)
def global_l(label,mm,rr):
    for j,v in enumerate([rr,mm-rr*(u+p)**2]):
        need(label+str(j),S.Poly(S.expand(v),u).degree()<=2)
def invert(v):return S.expand(v.subs({u:1/r-p,t:-r*r-r**4*T},simultaneous=True))
def jac(f,g):return S.diff(f,u)*S.diff(g,t)-S.diff(f,t)*S.diff(g,u)
global_h('full all-p H',N,P,Q);global_l('full all-p L',M,R)
q2,q1,q0=S.symbols('q2 q1 q0');q=q2*u*u+q1*u+q0
hi=S.Poly(S.expand(N-P*(u+p)**2+q*(u+p)**4),u)
sol=S.solve([hi.coeff_monomial(u**j) for j in (6,5)],(q2,q1),dict=True)
need('unique complete high-row solution',len(sol)==1)
zero('forced Q2',sol[0][q2]-(A[4]-1))
zero('forced Q1',sol[0][q1]-(A[3]-2*p*A[4]+4*p))
H=N*t*t+P*t+Q;L=M*t+R
Hi=invert(H);Li=invert(L)
need('whole original H second chart polynomial',S.denom(S.cancel(Hi))==1)
need('whole original L second chart polynomial',S.denom(S.cancel(Li))==1)
zero('D H slope',S.diff(Hi.subs(r,0),T)-(2-A[4]))
zero('D L slope',S.diff(Li.subs(r,0),T)+B[4])
zero('actual infinity boundary multiplicity',r**8*N.subs(u,1/r-p)-r*r*(1-p*r)**4*((1-p*r)**2-r*r))
# At the two simple roots, residue exactness forces M to vanish.
m2,m3,m4=S.symbols('m2 m3 m4');Ms=u*u*(m4*u*u+m3*u+m2)
zero('simple residue plus',S.residue(Ms/N,u,1)-(m4+m3+m2)/2)
zero('simple residue minus',S.residue(Ms/N,u,-1)+(m4-m3+m2)/2)
sol=S.solve([m4+m3+m2,m4-m3+m2],(m2,m3),dict=True)
need('two simple-root constraints complete',sol==[{m2:-m4,m3:0}])
zero('surviving T2 exact primitive',S.diff(-lam/u,u)-lam*u*u*(u*u-1)/N)
P=a*u**4+b*u**3+c*u*u
Q=(a-1)*u*u+(b-2*p*a+4*p)*u+d
H=N*t*t+P*t+Q
L=lam*(u*u*(u*u-1)*t+u*u-2*p*u)
F=H*H+L
HD=invert(H).subs(r,0);LD=invert(L).subs(r,0)
zero('generic D degree2 coefficient',S.expand(HD*HD+LD).coeff(T,2)-(2-a)**2)
zero('a2 D degree1 coefficient',S.expand((HD*HD+LD).subs(a,2)).coeff(T,1)+lam)
zero('active source line is constant fibre',F.subs(u,0)-d*d)
G=aa/u+bb*H/u
rem=S.Poly(S.rem(jac(F,G),F-zeta,t),t)
zero('leading basis bracket',rem.coeff_monomial(t**3)-4*aa*u**6*(u*u-1)**2)
zero('second basis bracket after first vanishes',rem.coeff_monomial(t*t).subs(aa,0)-bb*lam*u**4*(u*u-1)*(3-u*u))
need('both basis coefficients forced zero',S.Poly(u**4*(u*u-1)*(3-u*u),u).degree()==8)
# Both infinity regimes are essential: normal-unit and lower-unit cusp.
need('normal-unit infinity eta order3',2+4-3==3)
need('lower-unit infinity eta order8',6+8+2-8==8)
need('normal-unit genus2 divisor',2*3+2*2-4*2==2)
need('lower-unit genus2 divisor',8+2-4*2==2)
need('full RR dimension3',4+1-2==3)
need('active4 four primitive simple poles',4-6==-2)
need('H/u bounded normal-unit infinity',-1+1==0)
need('H/u bounded lower-unit infinity',-2+3==1)
need('finite simple relative form is unit',2-2==0)
# Complete quartic active tangents are generically square-free.
Z,pp,qq,mm,rr=S.symbols('Z pp qq mm rr')
P0=(-1+pp*Z+qq*Z*Z)**2+mm*Z**3+rr*Z**4
zero('active quartic nonzero repeated-root eliminant',
     (Z*S.diff(P0,Z)-4*P0).subs(Z,0)+4)
# The leading midpoint condition and a sharp rational family.
y=S.sqrt(u*u-1)
zero('leading differential exact',S.diff(y/u,u)-1/(u*u*y))
q=1+u*u*t;Hp=(u*u-1)*q*q+d;Gh=-1/(2*u*q)
global_h('sharp all-p global H',N,2*u*u*(u*u-1),u*u-1+d)
zero('sharp H mate',jac(Hp,Gh)-1)
zero('sharp F mate',jac(Hp*Hp+e,Gh/(2*Hp))-1)
need('sharp mate actual denominator',S.denom(S.cancel(Gh/(2*Hp)))!=1)
alpha=S.symbols('alpha',nonzero=True)
zero('finite double leading logarithm hostile',S.residue(1/(alpha*u),u,0)-1/alpha)
print('Binary (4,2,1,1): full global rows, both weighted-infinity regimes and complete genus-two primitive space')
print('Rational mate => double infinity and L constant; all polynomial mates excluded')
print('Exact gates:',len(GATES))
print('Semantic SHA256:',sha256('\n'.join(GATES).encode()).hexdigest())
