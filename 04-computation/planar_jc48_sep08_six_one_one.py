#!/usr/bin/env python3
"""Exact unsaturated inverse coefficients and the complete 6+1+1 reduction.

RESERVED pending independent audit.  No mate-degree bound or coefficient
census is used.  The final coprimality has a separate Fraction Euclidean
replay with a literal primitive remainder sequence.
"""
from fractions import Fraction as QF
from hashlib import sha256
from math import factorial,gcd,lcm
import json
import sympy as S

gates=[]
def need(label,value):
    if not bool(value):raise RuntimeError(label)
    gates.append(label)
def zero(label,value):need(label,S.cancel(S.together(value))==0)
def jac(F,G,x,t):return S.diff(F,x)*S.diff(G,t)-S.diff(F,t)*S.diff(G,x)

# Universal inverse coefficients, BEFORE any source-family specialization.
A=S.Symbol('A',nonzero=True)
D,M,Z,P,T=S.symbols('D M Z P T')
N=A*A
K=-D/(4*N)
expected={
 1:D/(8*A**3),
 2:-M/(4*A**2),
 3:-D*D/(128*A**5)-Z/(4*A),
 5:D**3/(1024*A**7)+(D*Z-M*M)/(32*A**3),
 7:-5*D**4/(32768*A**9)-3*D*D*Z/(512*A**5)
   -3*Z*Z/(32*A)-3*D*M*M/(256*A**5),
}
# Literal multinomial expansion of [y^-1]F^(j/4), including every term
# 2i+3j+4k=n+1.  Only the fixed positive integer n is divided here.
for n,wanted in expected.items():
    coefficient=0
    for i in range((n+1)//2+1):
        for j in range((n+1)//3+1):
            for k in range((n+1)//4+1):
                if 2*i+3*j+4*k!=n+1:continue
                count=i+j+k
                coefficient+=S.binomial(S.Rational(n,4),count)*S.Rational(
                    factorial(count),factorial(i)*factorial(j)*factorial(k))*\
                    (2*K/N)**i*(M/N**2)**j*((K*K+Z)/N**2)**k
    zero(f'unsaturated Lagrange inverse T{n}',-A**n*coefficient/n-wanted)

# Independent coefficient-by-coefficient formal inverse, with all four
# independent coefficients A,D,M,Z retained.  Coefficients not displayed
# above are also solved; no frozen prefix or specialization is imported.
def conv(a,b,n):
    return [S.expand(sum(a[j]*b[k-j] for j in range(k+1)
                         if j<len(a) and k-j<len(b))) for k in range(n+1)]
U=[1/A]
for n in range(1,9):
    trial=U+[0]
    square=conv(trial,trial,n)
    hh=[A*A*c for c in square]
    if n>=2:hh[2]+=K
    ee=conv(hh,hh,n)
    if n>=3:ee[n]+=M*trial[n-3]
    if n==4:ee[n]+=Z
    U.append(S.factor(-ee[n]/(4*A)))
    if n-1 in expected:
        zero(f'independent recursive inverse T{n-1}',U[n]-expected[n-1])
zero('centered inverse constant term',U[1])
q0,r0=S.symbols('q0 r0')
original=(N*T*T+P*T+q0)**2+M*T+r0
centered=(N*(T+P/(2*N))**2+q0-P*P/(4*N))**2+M*(T+P/(2*N))+r0-M*P/(2*N)
zero('unrestricted centering identity',original-centered)

# Full shifted global sections after necessary value/first-jet constraints.
x,u,t,p,a,b,c,k,kap,ell,w,f=S.symbols('x u t p a b c k kap ell w f')
C=u*u-2*u+3
qq=1+u*u*t-2*p/u
H=u*u*C*qq*qq+(b*u*u+c*u+a)*qq+k+2*p*a/u
L=kap*C*qq+ell+6*kap*p/u
H=S.cancel(H);L=S.cancel(L)
NN=u**6*C
PP=2*u**6-(4+4*p)*u**5+(6+8*p+b)*u**4+(c-12*p)*u**3+a*u*u
QQ=u**4-(2+4*p)*u**3+(3+8*p+4*p*p+b)*u*u+\
    (c-12*p-2*p*b-8*p*p)*u+a+k-2*p*c+12*p*p
MM=kap*u*u*C
RR=kap*u*u-kap*(2+2*p)*u+kap*(3+4*p)+ell
zero('entire shifted H source polynomial',H-(NN*t*t+PP*t+QQ))
zero('entire shifted L source polynomial',L-(MM*t+RR))
Nx=S.Poly(S.expand(NN.subs(u,x-p)),x)
Px=S.Poly(S.expand(PP.subs(u,x-p)),x)
Qx=S.Poly(S.expand(QQ.subs(u,x-p)),x)
need('full numerator degree eight',Nx.degree()==8)
zero('global source P6',Px.nth(6)-2*Nx.nth(8))
zero('global source P5',Px.nth(5)-2*Nx.nth(7))
zero('global source Q4',Qx.nth(4)-Nx.nth(8))
zero('global source Q3',Qx.nth(3)-Nx.nth(7))
zero('global source Q2',Qx.nth(2)-Px.nth(4)+Nx.nth(6))
zero('global source Q1',Qx.nth(1)-Px.nth(3)+Nx.nth(5))
Mx=S.Poly(S.expand(MM.subs(u,x-p)),x)
Rx=S.Poly(S.expand(RR.subs(u,x-p)),x)
zero('global source R2',Rx.nth(2)-Mx.nth(4))
zero('global source R1',Rx.nth(1)-Mx.nth(3))
# The remaining free source rows P2,P3,P4,Q0 have determinant one in
# (a,c,b,k); the constant L row is free independently of kap.
need('complete remaining H coefficient map determinant',
     S.Matrix([PP.coeff(u,2),PP.coeff(u,3),PP.coeff(u,4),QQ.coeff(u,0)]).jacobian([a,c,b,k]).det()==1)

# Explicit local logarithmic supplier at a=P2!=0,kap!=0.  First centre
# at the moving zero of the full numerator.  Its derivative in Z at
# u=0 is a, so analytic implicit inversion pays all higher units.
v,y,n6,n7,n8,p3,p4,m2,m3,Q0=S.symbols('v y n6 n7 n8 p3 p4 m2 m3 Q0')
z0=-n6/a
z1=-(n7+p3*z0)/a
ss=u**4*(z0+u*(z1+y))
local_s=S.Symbol('local_s')
local_n=n6*u**6+n7*u**7+n8*u**8+(a*u*u+p3*u**3+p4*u**4)*local_s+Q0*local_s**2
local_E=local_n**2+(m2*u*u+m3*u**3)*local_s**3-f*local_s**4
E=S.expand(local_E.subs(local_s,ss))
zero('sixfold shared logarithmic split',
     E.coeff(u,14)-(a*a*y*y+m2*z0**3))
zero('moving numerator centre retains arbitrary first units',
     S.expand(local_n.subs(local_s,ss)).coeff(u,7)-a*y)
zero('local logarithmic relative denominator leading coefficient',
     S.expand(S.diff(local_E,local_s).subs(local_s,ss)).coeff(u,9)-2*a*a*y)

# Set a=0 only AFTER the separate local argument.  The fibre remains a
# genuine quadratic FIELD extension over C(F,H), even if a geometric
# generic fibre later has a special singular model.
z=f-ell-w*w
EE=(u*z-6*kap*p)**2+kap*(b*u+c)*(u*z-6*kap*p)-kap*kap*C*(w-k)
AA=z*z+kap*b*z-kap*kap*(w-k)
BB=kap*(c-12*p)*z+kap*kap*(2*(w-k)-6*p*b)
CC=kap*kap*(-3*(w-k)+36*p*p-6*p*c)
zero('complete quadratic fibre elimination',EE-(AA*u*u+BB*u+CC))
q_sub=(z-6*kap*p/u)/(kap*C)
q_symbol=S.Symbol('q_symbol')
Hq=u*u*C*q_symbol*q_symbol+(b*u*u+c*u)*q_symbol+k
Lq=kap*C*q_symbol+ell+6*kap*p/u
zero('actual fibre elimination transport',kap*kap*C*(Hq.subs(q_symbol,q_sub)-w)-EE)
Fq=Hq*Hq+Lq
zero('actual relative form factor',
     jac(Fq,Hq,u,q_symbol).subs(q_symbol,q_sub)+S.diff(EE,u)/kap-S.diff(C,u)*EE/(kap*C))
disc=S.expand(BB*BB-4*AA*CC)
zero('quadratic field nonsquare leading coefficient',
     disc.coeff(f,2)-kap*kap*((c-12*p)**2+12*(w-k)-144*p*p+24*p*c))
zero('trace double-pole residue',
     S.diff(BB,w).subs(w,k+12*p*p-2*p*c)-
     (-2*kap*(c-12*p)*(k+12*p*p-2*p*c)+2*kap*kap))

# Odd inverse residues on the two unramified points above u=0.  Coefficients
# are scaled by sqrt(3); all vanish iff these rational polynomials vanish.
def series(exponent,n):
    return S.series((1-S.Rational(2,3)*u+u*u/3)**exponent,u,0,n+1).removeO()
def coeff(poly,power,n):return S.expand(poly*series(power,n)).coeff(u,n)
D2=(b*u+c)**2-4*k*C
zero('full discriminant after a zero',
     (PP*PP-4*NN*QQ).subs(a,0)-u**6*D2)
zero('full centered correction row',
     (RR-MM*PP/(2*NN)).subs(a,0)-(ell-kap*b/2-kap*(c-12*p)/(2*u)))
# Translate F by a constant to take ell=kap*b/2; the constant part of
# Z0 vanishes, while its displayed polar part remains.
r1=coeff(D2,-S.Rational(3,2),2)/24
r3=-coeff(D2*D2,-S.Rational(5,2),2)/1152+\
    kap*(c-12*p)*series(-S.Rational(1,2),3).coeff(u,3)/8
r5=coeff(D2**3,-S.Rational(7,2),2)/27648-\
    kap*(c-12*p)*coeff(D2,-S.Rational(3,2),3)/192-\
    3*kap*kap*series(S.Rational(1,2),4).coeff(u,4)/32
r7=-5*coeff(D2**4,-S.Rational(9,2),2)/(32768*81)+\
    3*kap*(c-12*p)*coeff(D2**2,-S.Rational(5,2),3)/(1024*9)-\
    3*kap*kap*(c-12*p)**2*series(-S.Rational(1,2),4).coeff(u,4)/128-\
    3*kap*kap*coeff(D2,-S.Rational(1,2),4)/256
# Independently obtain the same family residues directly from each full
# universal formula.  A formal root W has W^2=C; expand Laurent monomials
# in u and W, then extract their ordinary residue.  This does not use
# the four hand-collected expressions r1,r3,r5,r7 above.
W=S.Symbol('W',nonzero=True)
family_sub={A:u**3*W,D:u**6*D2,M:kap*u*u*C,Z:-kap*(c-12*p)/(2*u)}
zero('universal even row specializes before tracing',
     expected[2].subs(family_sub,simultaneous=True).subs(W**2,C)+kap/(4*u**4))
for idx,wanted in [(1,r1),(3,r3),(5,r5),(7,r7)]:
    raw=S.expand(expected[idx].subs(family_sub,simultaneous=True))
    got=0
    for term in S.Add.make_args(raw):
        powers=term.as_powers_dict()
        eu=int(powers.get(u,0));ew=int(powers.get(W,0));order=-1-eu
        if order<0:continue
        coefficient=S.cancel(term/(u**eu*W**ew))
        got+=coefficient*3**S.Rational(ew+1,2)*series(S.Rational(ew,2),order).coeff(u,order)
    zero(f'universal T{idx} actual family residue reconstruction',got-wanted)
quad=b*b+2*b*c+c*c/3
zero('first residue gives homogeneous conic',r1-quad/24)
red3=S.rem(S.together(r3).as_numer_denom()[0],quad,b)
zero('T3 reduced relation',red3-4*(6*b*c**3+c**4-12*c*kap+144*kap*p))
scale=S.Symbol('scale',nonzero=True)
for name,eq,weight in [('T1',r1,2),('T3',r3,4),('T5',r5,6),('T7',r7,8)]:
    zero(f'exact coefficient homogeneity {name}',eq.subs(
        {b:scale*b,c:scale*c,p:scale*p,k:scale**2*k,kap:scale**3*kap},simultaneous=True)-scale**weight*eq)
zero('trace coefficient homogeneity',
     ((c-12*p)*(k+12*p*p-2*p*c)-kap).subs(
       {c:scale*c,p:scale*p,k:scale**2*k,kap:scale**3*kap},simultaneous=True)-
     scale**3*((c-12*p)*(k+12*p*p-2*p*c)-kap))

# Only now use the nonzero c ratio and d=c-12p.  All denominators here are
# displayed and known nonzero from the prior trace/T3 branches.
d,V=S.symbols('d V',nonzero=True)
kv=(6*b+1)/(12*d*d)+(1-d*d)/12
lv=(6*b+1)/(12*d)
Qb=3*b*b+6*b+1
E5=(96*b+16)*V*V+(1296*b+240)*V+900*b+165
E7=(96*b+16)*V**4+(1296*b+240)*V**3+(180*b+33)*V*V+(14952*b+2744)*V+10098*b+1853
sub={c:1,p:(1-d)/12,k:kv,kap:lv}
reduced=[]
for name,eq,target in [('T5',r5,E5),('T7',r7,E7)]:
    numerator=S.together(eq.subs(sub)).as_numer_denom()[0]
    rem=S.rem(numerator,Qb,b)
    ratio=S.factor(rem/target.subs(V,d*d))
    need(f'{name} exact reduced factor independent of b',not ratio.has(b))
    # Require a nonzero rational constant times a power of d, with no
    # hidden parameter factor removed.  The actual ratio is printed below.
    num,den=S.fraction(ratio)
    need(f'{name} only declared d factor divided',
         len(S.Poly(num,d).terms())==1 and len(S.Poly(den,d).terms())==1)
    zero(f'{name} full reduced identity',rem-ratio*target.subs(V,d*d))
    reduced.append((name,str(ratio)))

R5=256*V**4+3072*V**3-2208*V*V-2880*V+225
R57=192*V**4-4320*V**3+2656*V*V+3252*V-255
zero('literal first resultant orientation',S.resultant(E5,Qb,b)-3*R5)
zero('literal second resultant orientation',S.resultant(E5,E7,b)+6*R57)

# Independent standard-library polynomial Euclidean replay.  Primitive
# rescaling divides only nonzero rational NUMBERS, never a parameter.
sequence=[
 [225,-2880,-2208,3072,256],
 [-255,3252,2656,-4320,192],
 [1695,-21648,-17248,26496],
 [305385,-4290180,1875232],
 [-27815,379612],
 [1],
]
def trim(poly):
    while len(poly)>1 and not poly[-1]:poly.pop()
    return poly
def remainder(A,B):
    A=list(map(QF,A));B=list(map(QF,B))
    while len(A)>=len(B) and any(A):
        j=len(A)-len(B);coef=A[-1]/B[-1]
        for k0 in range(len(B)):A[k0+j]-=coef*B[k0]
        trim(A)
    return A
def primitive(poly):
    den=1
    for q in poly:den=lcm(den,q.denominator)
    ints=[int(q*den) for q in poly];content=0
    for z0 in ints:content=gcd(content,abs(z0))
    if not content:raise RuntimeError('zero Euclidean remainder')
    if ints[-1]<0:content=-content
    return [z0//content for z0 in ints]
for j in range(2,len(sequence)):
    need(f'Fraction primitive Euclidean remainder {j}',
         primitive(remainder(sequence[j-2],sequence[j-1]))==sequence[j])
need('literal Euclidean final constant',sequence[-1]==[1])

# Sharp same-partition rational controls, at EVERY finite p.  They have
# L constant, exactly the branch left by the theorem.
vh=u+u**3*t-2*p
Hh=C*vh*vh+k
Gh=(u+1)/(12*u*u*vh)
zero('all-finite-position H rational hostile',jac(Hh,Gh,u,t)-1)
zero('same-partition F rational hostile',jac(Hh*Hh+ell,Gh/(2*Hh),u,t)-1)
rr,bb=S.symbols('rr bb')
chart=S.cancel(S.together(Hh.subs({u:1/rr-p,t:-rr*rr-rr**4*bb},simultaneous=True)))
need('rational hostile H actual whole boundary chart is polynomial',S.denom(chart)==1)
need('actual simple roots are distinct and away from sixfold point',
     S.discriminant(C,u)==-8 and C.subs(u,0)==3)

print('six_one_one: PASS')
print('scope: all finite sixfold positions; rational mate forces L constant; polynomial mate impossible')
print('universal inverse rows: T1,T2,T3,T5,T7; all coefficients retained before specialization')
print('formal path: independent recursive inverse and full multinomial residue formula')
print('normalization: homogeneous coefficient ratios only; c, kappa, c-12p branches explicit')
print('final reduced factors:',reduced)
print('primitive Fraction Euclidean sequence:',sequence)
print('hostile: actual global all-p constant-L family has a rational mate')
print('gates:',len(gates))
print('semantic sha256:',sha256(json.dumps([gates,sequence,reduced],separators=(',',':')).encode()).hexdigest())
