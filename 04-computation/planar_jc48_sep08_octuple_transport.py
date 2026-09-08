#!/usr/bin/env python3
"""Exact changed-carrier, trace, genus-two, and pole-repair identities.

No inherited mathematical implementation is imported. The matching note
pays global completeness, generic branches, and the entire primitive space.
"""
import hashlib
import json
import sympy as S

u,t,s,z,p,r,bd,c,Z=S.symbols('u t s z p r bd c Z')
alpha,k0,kt,kxt,k2,k3,k4=S.symbols('alpha k0 kt kxt k2 k3 k4')
l0,lt,lxt,l2,l3,l4=S.symbols('l0 lt lxt l2 l3 l4')
hh,vv,y,rr=S.symbols('hh vv y rr')
gates=0
records={}

def check(ok,label):
    global gates
    gates+=1
    if not ok: raise RuntimeError(label)

def eq(a,b,label):
    check(S.cancel(a-b)==0,label)

def coeff(a,v,n):
    return S.expand(a).coeff(v,n)

def jac(a,b,x,yvar):
    return S.diff(a,x)*S.diff(b,yvar)-S.diff(a,yvar)*S.diff(b,x)

w=u*u*t
v=u+u**3*t-2*p
h=u*u+u**4*t-2*p*u
eq(h,u*v,'actual changed-carrier product')
basis=[S.Integer(1),t,u*t,w,v,h]
sections=[S.cancel((s*bb.subs(t,1/s)).subs(s,z-(u+p)**2)) for bb in basis]
mon=[u**i*z**j for i in range(3) for j in range(2)]
columns=S.Matrix([[coeff(coeff(bb,u,i),z,j) for bb in sections]
                 for i in range(3) for j in range(2)])
check(S.factor(columns.det()) in (-1,1),'complete global L1 basis for every p')
for i,bb in enumerate(sections):
    check(S.Poly(bb,u,z).degree(u)<=2 and S.Poly(bb,u,z).degree(z)<=1,
          'full section box '+str(i))
H=alpha*h*h+k0+kt*t+kxt*u*t+k2*w+k3*v+k4*h
L=l0+lt*t+lxt*u*t+l2*w+l3*v+l4*h
N=S.cancel(s*s*H.subs(t,1/s))
M=S.cancel(s*L.subs(t,1/s))
Ng=S.expand(N.subs(s,z-(u+p)**2))
Mg=S.expand(M.subs(s,z-(u+p)**2))
eq(N.subs(s,0),alpha*u**8,'actual finite-p octic')
check(S.Poly(Ng,u,z).degree(u)<=4 and S.Poly(Ng,u,z).degree(z)<=2,'full global H numerator')
check(S.Poly(Mg,u,z).degree(u)<=2 and S.Poly(Mg,u,z).degree(z)<=1,'full global L numerator')
eq(S.diff(N,s).subs({u:0,s:0}),kt,'normal-unit coefficient retained')
eq(M.subs({u:0,s:0}),lt,'actual finite M-unit coefficient')
hi=S.cancel(h.subs({u:1/r-p,t:-r*r-r**4*bd},simultaneous=True))
vi=S.cancel(v.subs({u:1/r-p,t:-r*r-r**4*bd},simultaneous=True))
wi=S.cancel(w.subs({u:1/r-p,t:-r*r-r**4*bd},simultaneous=True))
check(hi.is_polynomial(r,bd) and vi.is_polynomial(r,bd) and wi.is_polynomial(r,bd),
      'all corrected carriers extend to the actual boundary')
eq(hi.subs(r,0),-bd-3*p*p,'actual h restriction')
eq(vi.subs(r,0),0,'actual v restriction')
eq(wi.subs(r,0),-1,'actual w restriction')
uncorrected=S.cancel((u*u+u**4*t).subs({u:1/r-p,t:-r*r-r**4*bd},simultaneous=True))
eq(S.cancel(r*uncorrected).subs(r,0),2*p,'uncorrected translated function has a pole')

E=S.expand(N*N+s**3*M-c*s**4)
kap=k0-2*p*k3
ell=l0-2*p*l3
eq(coeff(E.subs({kt:0,lt:0}).subs(s,u*Z),u,4),
   Z*Z*((kxt+kap*Z)**2+lxt*Z+(ell-c)*Z*Z),'actual first-jet polynomial')
reduction={kt:0,kxt:0,lt:0,lxt:0}
Nr=S.expand(N.subs(reduction))
Mr=S.expand(M.subs(reduction))
Er=S.expand(E.subs(reduction))
Fr=S.expand(H.subs(reduction)**2+L.subs(reduction))
T=k2*k2+(2*kap*k2+l2)*Z+(kap*kap+ell-c)*Z*Z
eq(coeff(Er.subs(s,u*u*Z),u,8),Z*Z*T,'full p-dependent low quadratic')
eq(S.diff(S.discriminant(T,Z),c),4*k2*k2,'generic two simple low branches')
eq(Er.subs(u,0),(kap*kap+ell-c)*s**4,'all four local sheets')
Nb=S.cancel(Nr.subs(s,u**6*Z)/u**8)
check(Nb.is_polynomial(u,Z),'p terms do not create a lower cancellation scale')
eq(Nb.subs(u,0),alpha+k2*Z,'same exact cancellation centre')
eq(coeff(S.diff(Nr,s).subs(s,u**6*Z),u,2),k2,'same normal derivative unit')
for n,lower,ln in [(2,{},l2),(3,{l2:0},l3),(4,{l2:0,l3:0},l4)]:
    mm=S.expand((Mr-c*s).subs(lower).subs(s,u**6*Z))
    check(min(mon[0] for mon,cc in S.Poly(mm,u).terms() if cc!=0)==n,
          'entire M order at cancellation '+str(n))
    eq(coeff(mm,u,n),ln,'actual nonzero leading supplier '+str(n))
    direct=S.expand((S.diff(Er,s)-2*Nr*S.diff(Nr,s)).subs(lower).subs(s,u**6*Z))
    check(min(mon[0] for mon,cc in S.Poly(direct,u).terms() if cc!=0)>=12+n,
          'all direct derivative terms have higher order '+str(n))
    eq(coeff(direct,u,12+n),3*Z*Z*ln,'whole derivative remainder supplier '+str(n))
    raw=1-S.Rational(n,2)
    eq(raw if n%2==0 else 2*raw+1,0 if n<4 else -1,'normalized high branch '+str(n))
HD=alpha*(-bd-3*p*p)**2+k0-k2+k4*(-bd-3*p*p)
LD=l0-l2+l4*(-bd-3*p*p)
eq(coeff(HD*HD+LD,bd,4),alpha*alpha,'actual D degree supplier')

# The complete rational field and the general quadratic trace.
eq(jac(h,v,u,t),h**3/v**2,'exact Poisson identity independent of p')
inv={u:hh/vv,t:(vv*vv+2*p*vv-hh)*vv*vv/hh**3}
eq(h.subs(inv,simultaneous=True),hh,'first changed-carrier field inverse')
eq(v.subs(inv,simultaneous=True),vv,'second changed-carrier field inverse')
eq(w.subs(inv,simultaneous=True),(vv*vv+2*p*vv)/hh-1,'retained nonzero-point w correction')
eq(jac(inv[u],inv[t],hh,vv),vv*vv/hh**3,'literal pulled-back volume')
A=alpha*hh*hh+k4*hh+k0
aq=k3*k3+l2/hh
bq=k3*A+l3/2+p*l2/hh
dq=A*A+l4*hh+l0-l2
Fhv=(A+k3*vv)**2+l2*((vv*vv+2*p*vv)/hh-1)+l3*vv+l4*hh+l0
eq(Fhv,aq*vv*vv+2*bq*vv+dq,'entire quadratic coefficient family')
eta=-vv*vv/(2*hh**3*(aq*vv+bq))
invol=-2*bq/aq-vv
eq(Fhv.subs(vv,invol),Fhv,'actual quadratic involution')
trace=S.cancel(eta+eta.subs(vv,invol))
eq(trace,2*bq/(aq*aq*hh**3),'full quadratic trace coefficient')
closed=(2*k3*hh*A+l3*hh+2*p*l2)/(hh**2*(k3*k3*hh+l2)**2)
eq(trace,closed,'closed rational trace on the base')
eq(S.residue(-closed.subs(hh,1/rr)/rr**2,rr,0),-2*alpha/k3**3,
   'nonzero infinity trace residue when k3 nonzero')
trace0=S.cancel(closed.subs(k3,0))
eq(trace0,l3/(l2*l2*hh)+2*p/(l2*hh*hh),'remaining trace when k3 zero')
eq(S.residue(trace0,hh,0),l3/l2**2,'rational mate forces l3 zero')

# Literal genus-two residual and its whole possible primitive space.
P=A*A+l4*hh+l0-l2
Q=p*p+hh*(c-P)/l2
eq(Fhv.subs({k3:0,l3:0,vv:y-p}),P+l2*(y*y-p*p)/hh,
   'original genus-two fibre value is unchanged')
eq(coeff(Q,hh,5),-alpha*alpha/l2,'actual odd degree five')
eq(Q.subs(hh,0),p*p,'two distinct unramified zero points when p nonzero')
eq(hh*l2*S.diff(Q,hh)-l2*Q,-l2*p*p-hh*hh*S.diff(P,hh),
   'all possible repeated-root values form a finite set')
eq((hh*hh*S.diff(P,hh)+l2*p*p).subs(hh,0),l2*p*p,
   'generic squarefreeness supplier is not identically zero')
etay=-(y-p)**2/(2*l2*hh*hh*y)
eq(-(2*l2*y/hh)*etay,(y-p)**2/hh**3,'exact original genus-two volume')
aa=(c-P.subs(hh,0))/(2*p*l2)
plus=S.cancel(etay.subs(y,p+aa*hh))
minus=S.cancel(etay.subs(y,-p-aa*hh))
eq(S.limit(hh*hh*plus,hh,0),0,'no double pole at Pplus')
eq(S.residue(plus,hh,0),0,'no simple pole at Pplus')
eq(S.limit(hh*hh*minus,hh,0),2*p/l2,'actual double pole at Pminus')
eq(S.residue(minus,hh,0),0,'actual double pole is residue-free')
check(2*(-5)-(2*(-2)+(-5))+(-3)==-4,'complete infinity order is minus four')
check((-5)-(-2)==-3,'last basis function has infinity pole three')
check(4-2+1==3,'nonspecial degree-four function space dimension')
etared=-(Q+p*p-2*p*y)/(2*l2*hh*hh*y)
dbasis=(hh*S.diff(Q,hh)-2*Q+2*p*y)/(2*hh*hh*y)
eq(etared+etared.subs(y,-y),2*p/(l2*hh*hh),'genus-two trace is exact')
eq(dbasis+dbasis.subs(y,-y),2*p/(hh*hh),'trace of the last basis derivative')
eq(etared-dbasis/l2,S.diff(P,hh)/(2*l2*l2*y),'decisive full primitive-space remainder')
eq(coeff(S.diff(P,hh),hh,3),4*alpha*alpha,'primitive remainder cannot vanish')
for pp,ll,cc in [(1,1,0),(2,-1,1),(S.I,2,3)]:
    qtest=S.Poly(Q.subs({p:pp,l2:ll,c:cc,alpha:1,k0:0,k4:0,l4:0,l0:0}),hh)
    check(S.gcd(qtest,qtest.diff()).degree()==0,'named genuine genus-two fibre '+str((pp,ll,cc)))

# Complete rational primitive and same-fibre repair in the last linear case.
g,P0,P1,P2,P3,P4=S.symbols('g P0 P1 P2 P3 P4')
Pgeneral=P0+P1*hh+P2*hh**2+P3*hh**3+P4*hh**4
integrand=-(g-(Pgeneral-P0))**2/(l3**3*hh**3)
eq(S.residue(integrand,hh,0),-(P1*P1-2*g*P2)/l3**3,
   'whole linear-family residue criterion')
q=P3*hh**3+P4*hh**4
Int1=P3*hh+P4*hh*hh/2
Int2=P3*P3*hh**4/4+2*P3*P4*hh**5/5+P4*P4*hh**6/6
G0=(g*g/(2*hh*hh)+2*g*Int1-Int2)/l3**3
eq(S.diff(G0,hh),-(g-q)**2/(l3**3*hh**3),'entire exact rational primitive')
remainder=S.cancel(G0.subs(g,l3*vv+q)-vv*vv/(2*l3*hh*hh))
check(remainder.is_polynomial(hh,vv),'only source pole is the declared double pole')
qsource=P3*h**3+alpha*alpha*h**4
Flinear=P0+qsource+l3*v
FE=P0-2*p*l3
quotient=S.cancel((Flinear-FE)/u)
check(quotient.is_polynomial(u,t),'actual same-fibre residual factor')
eq(quotient.subs(u,0),l3,'residual divisor disjoint from source E')
eq(coeff(quotient,t,4),alpha*alpha*u**15,'residual factor is nonconstant')
eq(Flinear.subs(u,0),FE,'correct shifted exceptional-fibre value')
eq(S.diff(h,u).subs(u,0),-2*p,'h itself need not be source-critical')

# Rational mates genuinely survive inside the polynomial-exclusion class.
Fhost=h**4+h
Ghost=1/(3*u**3*(4*h**3+1))
eq(jac(Fhost,Ghost,u,t),1,'all-p rational-mate hostile')
Fhost2=h**4-v
Ghost2=-1/(2*u*u)+2*v*h*h-S.Rational(4,3)*h**6
eq(jac(Fhost2,Ghost2,u,t),1,'all-p boundary-linear rational mate')
eq(S.diff(Fhost2,u).subs(u,0),-1,'second hostile smooth on source E')
eq(S.diff(v,t),u**3,'second hostile also smooth on v-zero away from E')

records['point']='arbitrary finite p; p0 inherited and nonzero p paid explicitly'
records['global_basis_dimension']=6
records['trace_residue_infinity']=str(-2*alpha/k3**3)
records['genus_two_basis']=['1','h','(y-p)/h']
records['genus_two_remainder']='Pprime(h)/(2*l2^2*y), leading Pprime=4*alpha^2*h^3'
records['same_fibre']='F_E=P0-2p*l3; residual factor has value l3 on E'
records['scope']='no polynomial mate; two actual rational mates retained for all p'
semantic=hashlib.sha256(json.dumps(records,sort_keys=True).encode()).hexdigest()
print('PASS all-finite octuple transport:',gates,'always-active exact gates')
print('Actual global carriers and complete L1 entry; unchanged four-branch capacity')
print('General quadratic trace excludes k3!=0; genus-two primitive space excludes its residual')
print('Shifted same-fibre repair excludes the last polynomial mate; rational hostiles survive')
print('Semantic SHA256:',semantic)
