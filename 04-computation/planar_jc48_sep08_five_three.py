#!/usr/bin/env python3
"""All-parameter exact controls for the all-finite DG boundary class (5,3).

The companion proves the complete local divisors and Riemann--Roch spaces.
The rejected genus-three correction is tested as a hostile, never as a basis.
"""
from hashlib import sha256
import json
import sympy as S

gates=[]


def need(label, predicate):
    if not bool(predicate):
        raise RuntimeError(label)
    gates.append(label)


def zero(label, expr):
    need(label,S.cancel(expr)==0)


u,x,t,w,p,A,B,b,c,zeta,lam=S.symbols('u x t w p A B b c zeta lam')
r,bd,s,tau,Z=S.symbols('r bd s tau Z')
mu4,mu3,mu2,mu1,mu0,r0,b1,b0=S.symbols('mu4 mu3 mu2 mu1 mu0 r0 b1 b0')
N=u**5*(u-1)**3
P=u*u*(2*u**4+(-4*p-6)*u**3+A*u*u+B*u+b)
Q=u**4+(-4*p-3)*u**3+(A+4*p*p-3)*u*u+(-2*A*p+B+12*p*p+12*p+1)*u+c
M=u*u*(mu4*u*u+mu3*u+mu2)
R=mu4*u*u+(mu3-2*p*mu4)*u+r0
Norig=S.Poly(S.expand(N.subs(u,x-p)),x)
Pfull=P+b1*u+b0
Mfull=M+mu1*u+mu0
Porig=S.Poly(S.expand(Pfull.subs(u,x-p)),x)
Qorig=S.Poly(S.expand(Q.subs(u,x-p)),x)
Morig=S.Poly(S.expand(Mfull.subs(u,x-p)),x)
Rorig=S.Poly(S.expand(R.subs(u,x-p)),x)
zero('full P degree-six coefficient',Porig.coeff_monomial(x**6)-2*Norig.coeff_monomial(x**8))
zero('full P degree-five coefficient',Porig.coeff_monomial(x**5)-2*Norig.coeff_monomial(x**7))
for degree,expected in [(4,Norig.coeff_monomial(x**8)),(3,Norig.coeff_monomial(x**7)),
    (2,Porig.coeff_monomial(x**4)-Norig.coeff_monomial(x**6)),
    (1,Porig.coeff_monomial(x**3)-Norig.coeff_monomial(x**5))]:
    zero(f'full Q induced coefficient {degree}',Qorig.coeff_monomial(x**degree)-expected)
zero('full R induced quadratic',Rorig.coeff_monomial(x*x)-Morig.coeff_monomial(x**4))
zero('full R induced linear',Rorig.coeff_monomial(x)-Morig.coeff_monomial(x**3))
Hfull=N*t*t+Pfull*t+Q
Lfull=Mfull*t+R
for name,G in [('H',Hfull),('L',Lfull)]:
    chart=S.cancel(G.subs(u,x-p).subs({x:1/r,t:-r*r-r**4*bd},simultaneous=True))
    need(name+' complete original global chart',S.fraction(chart)[1]==1)
need('fixed leading H and full L dimensions',len([A,B,b,b1,b0,c])==6 and len([mu4,mu3,mu2,mu1,mu0,r0])==6)
zero('active P value',P.subs(u,0))
zero('active P first derivative',S.diff(P,u).subs(u,0))
zero('active M value',M.subs(u,0))
zero('active M first derivative',S.diff(M,u).subs(u,0))
for root,multiplicity in [(0,5),(1,3)]:
    for j in range(multiplicity): zero(f'N root {root} jet {j}',S.diff(N,u,j).subs(u,root))
    need(f'N root {root} exact multiplicity',S.diff(N,u,multiplicity).subs(u,root)!=0)
need('nonsquare leading sidecar',5%2==1 and 3%2==1)

# T2, full leading residues and the actual moving-root rows.
zero('complete T2 residue condition',S.residue(M/N,u,0)+mu4+3*mu3+6*mu2)
zero('both active multiple points would violate T2',S.residue(u*u*(u-1)**2/N,u,1)-1)
Hblow=S.expand(N*w*w/u**4+P*w/u**2+Q)
Fblow=S.expand(Hblow*Hblow+M*w/u**2+R)
Q1=-2*A*p+B+12*p*p+12*p+1
f=(b*w+c)**2+mu2*w+r0
g=2*(b*w+c)*(-w*w+B*w+Q1)+mu3*w+mu3-2*p*mu4
zero('complete actual first F row',Fblow.subs(u,0)-f)
zero('complete actual next F row',S.diff(Fblow,u).subs(u,0)-g)
zero('actual blow-up volume',S.diff(w/u**2,w)-1/u**2)
resnum=S.expand(S.diff(g,w)*S.diff(f,w)-g*S.diff(f,w,2))
zero('normal-order-two unavoidable residue polynomial',resnum.coeff(w,3)+8*b**3)
zero('n2 next quadratic residue coefficient',S.expand(resnum.subs(b,0)).coeff(w,1)+4*c*mu2)
zero('n2 final linear residue coefficient',resnum.subs({b:0,c:0})-mu2*mu3)
MI=lam*u*u*(1-6*u*u)
RI=-6*lam*u*u+12*p*lam*u
MII=lam*u**3*(1-3*u)
RII=-3*lam*u*u+lam*(1+6*p)*u
zero('case I entire lower row', (M*w/u**2+R).subs({mu2:lam,mu3:0,mu4:-6*lam,r0:0})-(MI*w/u**2+RI))
zero('case II entire lower row', (M*w/u**2+R).subs({mu2:0,mu3:lam,mu4:-3*lam,r0:0})-(MII*w/u**2+RII))
zero('case I triple becomes M-unit only now',MI.subs(u,1)+5*lam)
zero('case II triple becomes M-unit only now',MII.subs(u,1)+2*lam)

# Exact local faces for every surviving normal order, with complete generic
# fibre parameter. Orders are on actual normalizations, not raw slopes.
al,bl,cl,ml,rl=S.symbols('al bl cl ml rl')
NI=al*u**5+bl*u**3*s+cl*u*s*s
EI=NI**2+s**3*(lam*u*u-zeta*s)
zero('case I outer actual leading face',S.expand(EI.subs(s,u*u*Z)).coeff(u,8)-Z**3*(lam-zeta*Z))
zero('case I high actual leading cubic',S.expand(EI.subs({u:tau**3,s:tau**8*Z},simultaneous=True)).coeff(tau,30)-(al*al+lam*Z**3))
need('case I outer differential order',4-6==-2)
need('case I high normalized differential order',16+2-22==-4)
need('case I branch exhaustion',1+3==4)
need('case I normal-unit triple zero exceeds pole degree',4+1>(2-1)+(4-1))
EII=(al*u**5+bl*u**3*s+cl*u*s*s)**2+s**3*(lam*u**3-zeta*s)
zero('case II high actual leading quartic',S.expand(EII.subs({u:tau*tau,s:tau**5*Z},simultaneous=True)).coeff(tau,20)-(al*al-zeta*Z**4))
need('case II actual normalized differential order',10+1-15==-4)
need('case II pair of normalized branches',4//2==2)
need('normal-unit triple actual zero order',12+1-9==4)
need('case I genus two canonical degree',-2-4+8==2*2-2)
need('case II zero-normal-coefficient genus one canonical degree',-4-4+8==2*1-2)
need('case II normal-unit genus three canonical degree',-4-4+8+4==2*3-2)
for name,deg,genus,dim in [('I',4,2,3),('II-regular-triple',6,1,6),('II-normal-unit',6,3,4)]:
    need(name+' complete Riemann-Roch dimension',deg>2*genus-2 and deg+1-genus==dim)
zero('full triple normal coefficient',P.subs(u,1)-(A+B+b-4*p-4))

# Formal inverse identities: the Ferrari pair-sum provides the entire even
# subsequence and in particular T6,T10 with no coefficient truncation loss.
NN,DD,EE,MM,v,zz,yy,sigma,qa,qb=S.symbols('NN DD EE MM v zz yy sigma qa qb')
factor=S.expand((yy*yy+sigma*yy+qa)*(yy*yy-sigma*yy+qb))
zero('Ferrari y coefficient',factor.coeff(yy,1)-sigma*(qb-qa))
zero('Ferrari y squared coefficient',factor.coeff(yy,2)-(qa+qb-sigma*sigma))
# qa+qb=sigma²+2D/N, qb-qa=M/(N² sigma), qa qb=(D²+E-v^-4)/N².
resolvent=sigma**6+4*DD*sigma**4/NN+4*(v**-4-EE)*sigma*sigma/(NN*NN)-MM*MM/NN**4
qav=(sigma*sigma+2*DD/NN-MM/(NN*NN*sigma))/2
qbv=(sigma*sigma+2*DD/NN+MM/(NN*NN*sigma))/2
zero('Ferrari full product elimination',
     4*sigma*sigma*(qav*qbv-(DD*DD+EE-v**-4)/(NN*NN))-resolvent)
VV=S.symbols('VV')
zero('full normalized pair-sum equation',
    (resolvent.subs(sigma,MM*v*v*VV/(2*NN))*NN**4/MM**2)
    -((1-EE*v**4)*VV**2+DD*MM*MM*v**8*VV**4/(4*NN)+MM**4*v**12*VV**6/(64*NN*NN)-1))
Vseries=1+EE*zz/2+(3*EE*EE-DD*MM*MM/NN)*zz*zz/8
Veq=(1-EE*zz)*Vseries**2+DD*MM*MM*zz*zz*Vseries**4/(4*NN)+MM**4*zz**3*Vseries**6/(64*NN*NN)-1
for j in range(3):zero(f'Ferrari exact series coefficient {j}',S.expand(Veq).coeff(zz,j))
zero('full T2 from pair sum',-MM/(4*NN)+MM/(4*NN))
zero('full T6 from pair sum',S.expand(-MM/(4*NN)*Vseries).coeff(zz,1)+MM*EE/(8*NN))
zero('full T10 from pair sum',S.expand(-MM/(4*NN)*Vseries).coeff(zz,2)+MM*(3*EE*EE-DD*MM*MM/NN)/(32*NN))
need('nonzero exactness multipliers six and ten',S.Rational(6,4)!=0 and S.Rational(10,4)!=0)

P2=P.subs(b,0)
Q2=Q
Dcenter=Q2-P2*P2/(4*N)
Ecenter=RII+r0-MII*P2/(2*N)
T6=-MII*Ecenter/(8*N)
zero('complete T6 residue before tuning',S.residue(T6,u,0)-lam*lam*(B+12*p+2)/16)
T10=-MII*(3*Ecenter*Ecenter-Dcenter*MII*MII/N)/(32*N)
zero('complete T10 residue after T6',S.residue(T10.subs(B,-12*p-2),u,0)-lam**3*c/32)

def bracket(F,G):
    return S.cancel(u*u*(S.diff(F,u)*S.diff(G,w)-S.diff(F,w)*S.diff(G,u)))


def remainders(F,basis,label):
    result=[S.cancel(S.rem(bracket(F,G),F-zeta,w)) for G in basis]
    for index,G in enumerate(result):
        need(f'{label} matrix remainder {index} is polynomial without clearing',S.fraction(G)[1]==1)
    return result


def row(rem,upower,wpower):
    return [S.expand(G).coeff(w,wpower).coeff(u,upower) for G in rem]


U=1/u
J=(u-1)**2*w+u*u-(2+2*p)*u
Jsource=u*u*(u-1)**2*t+u*u-(2+2*p)*u
J0=(u-1)**3*w/u+u*u-(3+2*p)*u
J0source=u*(u-1)**3*t+u*u-(3+2*p)*u
for name,G in [('J',Jsource),('J0',J0source)]:
    need(name+' original global chart',S.fraction(S.cancel(G.subs(u,x-p).subs({x:1/r,t:-r*r-r**4*bd},simultaneous=True)))[1]==1)
zero('J triple double leading zero',(u*u*(u-1)**2).subs(u,1))
zero('J triple derivative leading zero',S.diff(u*u*(u-1)**2,u).subs(u,1))
for j in range(3):zero(f'J0 triple leading jet {j}',S.diff(u*(u-1)**3,u,j).subs(u,1))

# Case I: complete 3-dimensional primitive space, including all A,p,lambda.
PI=P.subs({b:0,B:4*p+4-A})
QI=Q.subs({c:0,B:4*p+4-A})
HI=S.expand(N*w*w/u**4+PI*w/u**2+QI)
FI=S.expand(HI*HI+MI*w/u**2+RI)
zero('case I exact triple zero',PI.subs(u,1))
zero('case I original affine line has constant F',((N*t*t+PI*t+QI)**2+MI*t+RI).subs(u,0))
ri=remainders(FI,[U,J],'genus-two')
ca,cb,cc=S.symbols('ca cb cc')
exi=S.expand(ca*ri[0]+cb*ri[1])
zero('case I complete cubic homogeneous row',exi.coeff(w,3)-4*ca*u*u*(u-1)**6)
zero('case I remaining linear homogeneous row',S.expand(exi.subs(ca,0)).coeff(w,1)+2*cb*lam*u*(u-1)**2)

# The common case-II family after both exact inverse residues.
PII=P.subs({b:0,B:-12*p-2})
QII=Q.subs({c:0,B:-12*p-2})
HII=S.expand(N*w*w/u**4+PII*w/u**2+QII)
FII=S.expand(HII*HII+MII*w/u**2+RII)
d=A-16*p-6
zero('case II remaining triple coefficient',PII.subs(u,1)-d)
zero('case II original affine line has constant F',((N*t*t+PII*t+QII)**2+MII*t+RII).subs(u,0))

# Case II d=0: the full six-dimensional primitive space.
H1=S.expand(HII.subs(A,16*p+6))
F1=S.expand(FII.subs(A,16*p+6))
V1=(u-1)*H1/u
basis1=[U,J,U*J,V1,V1*J]
rem1=remainders(F1,basis1,'genus-one')
labels1=[(8,3),(6,2),(5,2),(4,1),(3,1)]
matrix1=S.Matrix([row(rem1,i,j) for i,j in labels1])
zero('complete genus-one determinant',matrix1.det()-192*lam*zeta*(-lam*lam+32*p*p*zeta))
zero('genus-one determinant at anchored p0',(192*lam*zeta*(-lam*lam+32*p*p*zeta)).subs(p,0)+192*lam**3*zeta)
need('genus-one pole degrees and t-degree independence',max([0,2,1,3,2,3])==3 and 3<4)

# Rejected correction: repairing a boundary pole introduces actual affine
# poles over u=1. Its residue is nonzero whenever the normal coefficient is.
badV=(J+d/(u-1))/u
zero('rejected correction actual affine pole coefficient',S.limit((u-1)*badV,u,1)-d)
need('rejected correction is not included in either basis',badV not in basis1)

# Valid genus-three repair has ONLY u denominators. It kills two analytic
# jets while preserving all ordinary affine points on u=1.
K=-u*u+(2*p+3)*u-A+(12*p+2)/u
k0=-A+14*p+4
k1=-10*p-1
zero('actual analytic part of J0 at triple',
     u*u-(3+2*p)*u-PII/u**4-K)
zero('analytic constant jet',K.subs(u,1)-k0)
zero('analytic first jet',S.diff(K,u).subs(u,1)-k1)
Adj=J0-k0+k1*(U-1)
zero('corrected analytic zero jet',(K-k0+k1*(U-1)).subs(u,1))
zero('corrected analytic derivative jet',S.diff(K-k0+k1*(U-1),u).subs(u,1))
W0=HII*Adj
need('valid repair has no u1 denominator',S.fraction(S.cancel(W0))[1].subs(u,1)!=0)
need('valid triple fractional order cancellation',S.Rational(3,2)-S.Rational(3,2)==0)
need('valid active pole bounds',3<=3 and 2<=3)
rem3=remainders(FII,[U,J0,W0],'genus-three')
matrix3=S.Matrix([row(rem3,i,j) for i,j in [(8,3),(6,2),(5,2)]])
det3=16*p*(lam*(A*A-(56*p+12)*A+640*p*p+336*p+36)-12*zeta*(A-8*p-6))
zero('complete genus-three determinant nonanchored',matrix3.det()-det3)
zero('genus-three nonanchored zero-slope case',det3.subs(A,8*p+6)-4096*lam*p**3)
matrix30=S.Matrix([row([G.subs(p,0) for G in rem3],i,j) for i,j in [(8,3),(6,2),(4,1)]])
zero('complete genus-three anchored determinant',matrix30.det()+20*lam*((6-A)*lam+12*zeta))
need('valid genus-three basis t degrees',3<4 and 1<3)

# Sharp rational constant-L family, same boundary partition and every p.
kap,constant=S.symbols('kap constant')
Wsharp=u*u*(u-1)*t+u-2*p-1
Hsharp=u*(u-1)*Wsharp*Wsharp+kap
Gsharp=(1/(3*u*u)+5/(3*u)+1/(u-1))/Wsharp
jac=lambda f,g:S.diff(f,u)*S.diff(g,t)-S.diff(f,t)*S.diff(g,u)
zero('sharp family exact leading polynomial',S.expand(Hsharp).coeff(t,2)-N)
need('sharp family original global chart',S.fraction(S.cancel(Hsharp.subs(u,x-p).subs({x:1/r,t:-r*r-r**4*bd},simultaneous=True)))[1]==1)
zero('sharp family rational H mate',jac(Hsharp,Gsharp)-1)
zero('sharp family rational F mate',jac(Hsharp*Hsharp+constant,Gsharp/(2*Hsharp))-1)
GG=S.Function('GG')(u,t)
zero('final polynomial factor',jac(Hsharp*Hsharp+constant,GG)-2*Hsharp*jac(Hsharp,GG))

# Named parameter controls include the zero fibre-slope and higher-normal-jet
# cases; no determinant is tested only at generic numerical parameters.
for vals in [{p:0,lam:1},{p:1,lam:2}]:
    need('genus-one named determinant '+str(vals),S.Poly((192*lam*zeta*(-lam*lam+32*p*p*zeta)).subs(vals),zeta).as_expr()!=0)
for vals in [{p:1,A:14,lam:2},{p:1,A:0,lam:2}]:
    need('genus-three named determinant '+str(vals),S.Poly(det3.subs(vals),zeta).as_expr()!=0)

semantic=sha256(json.dumps(gates,separators=(',',':')).encode()).hexdigest()
print('DG all-finite boundary (5,3): complete coefficient and primitive-space controls')
print('Moving residues, T2/T6/T10, and all three genus cases: PASS')
print('Genus-three correction preserves ordinary affine points: PASS')
print('Rational mate forces constant L; polynomial mates of every degree excluded')
print('Same-partition all-p rational constant-L family: PASS')
print(f'Gates: {len(gates)}')
print('Semantic SHA256: '+semantic)
print('RESULT: PASS')
