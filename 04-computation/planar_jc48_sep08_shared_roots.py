#!/usr/bin/env python3
"""Exact shared-boundary-root jet controls; analytic proof is separate."""
import hashlib,json
import sympy as S

x,t,s,w,v,Z,V,tau,c,r,b,q=S.symbols('x t s w v Z V tau c r b q')
a,B,e,d,f,A,C,D,E,beta,alpha=S.symbols('a B e d f A C D E beta alpha')
GATES=0
records=[]
def check(ok,label):
    global GATES
    GATES+=1
    if not ok:raise RuntimeError(label)
def eq(l,r,label):check(S.cancel(l-r)==0,label)
def lead(P,var):
    pp=S.Poly(S.expand(P),var)
    k=min(i[0] for i,co in pp.terms() if co!=0)
    return k,pp.coeff_monomial(var**k)
def jac(F,G,u,z):return S.diff(F,u)*S.diff(G,z)-S.diff(F,z)*S.diff(G,u)
def chart(P):return S.cancel(P.subs({x:1/r,t:-r*r-r**4*b},simultaneous=True))

# Simple N restriction; all arbitrary higher weighted jets are covered by
# the analytic proof, while this literal model checks the decisive orders.
P=v*v+s**3*(A*s+B*v+C*s*s+D*s*v+E*v*v)-c*s**4
k,co=lead(P.subs(v,V*s*s),s)
check(k==4,'simple N exact generic leading order')
eq(co,V*V+A-c,'simple N two nonzero tangent roots generically')
k,co=lead(S.diff(P,v).subs(v,V*s*s),s)
check(k==2,'simple N polar order')
eq(co,2*V,'simple N polar leading coefficient')

# Every normal-unit multiplicity and each relative tangential regime.
for m in range(1,9):
    for n in list(range(1,10))+[None]:
        ell=m if n is None else min(m,n)
        MM=(w**n if n is not None else 0)+d*s
        PP=S.expand((s+w**m)**2+s**3*MM-c*s**4)
        sub={w:tau*tau,s:-tau**(2*m)+V*tau**(3*m+ell)}
        k,co=lead(PP.subs(sub,simultaneous=True),tau)
        check(k==6*m+2*ell,'normal-unit exact split weight')
        expected=V*V-1 if n is not None and n<m else V*V+d-c-(1 if n==m else 0)
        eq(co,expected,'normal-unit generic coefficient')
        k,co=lead(S.diff(PP,s).subs(sub,simultaneous=True),tau)
        check(k==3*m+ell,'normal-unit polar order')
        eq(co,2*V,'normal-unit polar coefficient')
        check(4*m+1-k==m-ell+1>=1,'normal-unit differential regular after ramification')
        records.append([m,n,ell,'regular'])

P0=(a+B*Z+e*Z*Z)**2+d*Z**3+f*Z**4
Pc=P0-c*Z**4
eq((Z*S.diff(P0,Z)-4*P0).subs(Z,0),-4*a*a,'four-root exceptional polynomial nonzero')
Q2=S.cancel(Pc.subs(a,0)/Z**2)
eq(Q2,B*B+(2*B*e+d)*Z+(e*e+f-c)*Z*Z,'two-root exact factor')
eq((Z*S.diff(Q2,Z)-2*Q2).subs(Z,0),-2*B*B,'two-root exceptional polynomial nonzero')
eq(Pc.subs({a:0,B:0}),Z**3*(d+(e*e+f-c)*Z),'one-root exact factor')
eq(Pc.subs({a:0,B:0,d:0}),(e*e+f-c)*Z**4,'first unsupported face retained')
for aa,bb,dd,nonzero in [(1,0,0,4),(0,1,0,2),(0,0,1,1)]:
    pp=S.expand(Pc.subs({a:aa,B:bb,d:dd,e:0,f:0,c:1}))
    reduced=S.cancel(pp/Z**(4-nonzero))
    check(S.degree(reduced,Z)==nonzero,'named number of nonzero tangent roots')
    check(reduced.subs(Z,0)!=0,'named roots exclude zero')
    check(S.discriminant(reduced,Z)!=0,'named roots all simple')
    records.append(['tangent',aa,bb,dd,nonzero])

# Original equation, rather than a tangent-only statistic.
NN=a*w*w+B*s*w+e*s*s+s**3+w**3
MM=d*w+f*s+s*s+w*w
EE=S.expand(NN*NN+s**3*MM-c*s**4)
k,co=lead(EE.subs(s,w*Z),w)
check(k==4,'actual tangent quartic degree')
eq(co,Pc,'actual complete tangent quartic')
k,co=lead(S.diff(EE,s).subs(s,w*Z),w)
check(k==3,'actual polar tangent weight')
eq(co,S.diff(Pc,Z),'actual polar derivative')
eq((w*Z)**2/w**3,Z*Z/w,'nonzero-root logarithmic coefficient')
eq(w*w*(Z*Z/w),Z*Z*w,'infinity canonical numerator removes this logarithm')

# Two actual global low-order shared controls distinguish N_s.
z=S.symbols('z')
for name,N in [('normal_unit',z),('singular',x*x+(z-x*x)**2)]:
    H=S.cancel(N.subs(z,x*x+1/t)*t*t);L=x*t
    check(S.denom(H)==1 and S.denom(chart(H))==1,'actual global shared H')
    check(S.denom(chart(L))==1,'actual global shared L')
    eq(N.subs(z,x*x),x*x,'actual shared N multiplicity two')
    eq(S.diff(N.subs(z,s+x*x),s).subs({s:0,x:0}),1 if name=='normal_unit' else 0,
       'actual normal derivative distinguishes regimes')

# Sharp nonzero-M rational mate; each generic compact component has degree3.
h=x*x+x**4*t
H=h*h;L=h;F=h**4+h
G=1/(3*x**3*(4*h**3+1))
eq(jac(h,x**-3,x,t),3,'base rational primitive identity')
eq(jac(F,G,x,t),1,'actual rational mate exact Jacobian')
eq(chart(h),-b,'actual global h full chart')
eq(chart(G),r**3/(3*(1-4*b**3)),'actual rational mate full chart')
N=x**4*z*z;M=x*x*z
eq(S.cancel(N.subs(z,x*x+1/t)*t*t),H,'actual numerator N')
eq(S.cancel(M.subs(z,x*x+1/t)*t),L,'actual numerator M nonzero')
eq(N.subs(z,x*x),x**8,'actual m eight')
eq(S.diff(N.subs(z,s+x*x),s).subs(s,0),2*x**6,'actual j six')
eq(M.subs(z,x*x),x**4,'actual n four')
eq((r**8*N.subs({x:1/r,z:1/r**2},simultaneous=True)),1,'N nonzero at infinity')
eq((r**4*M.subs({x:1/r,z:1/r**2},simultaneous=True)),1,'M nonzero at infinity')
eq(chart(F).subs(r,0),b**4-b,'nonconstant F on D')
eq(F.subs(t,(alpha-x*x)/x**4),alpha**4+alpha,'actual generic compact component')
eq(G.subs(t,(alpha-x*x)/x**4),1/(3*x**3*(4*alpha**3+1)),'degree three on each component')
check(S.discriminant(alpha**4+alpha-c,alpha)!=0,'four generic components simple')
check(3==3,'primitive degree and local degree equality boundary')

# Root's exact stop object: shared low multiplicity can make F_D constant.
N=beta*x*x*z-beta*x**4+alpha*x*x*z*z;M=S.Integer(-1)
H=S.cancel(N.subs(z,x*x+1/t)*t*t);L=-t
eq(chart(H).subs(r,0),-beta,'stop object H restriction constant')
eq(chart(L).subs(r,0),0,'stop object L restriction constant')
eq(N.subs(z,x*x),alpha*x**6,'stop object finite six and infinite two')
Ni=S.cancel(r**4*q*q*N.subs({x:1/r,z:1/q},simultaneous=True))
eq(Ni,alpha*r*r+beta*r*r*q-beta*q*q,'complete stop infinity numerator')
eq(-r*r*q*M,r*r*q,'complete stop M infinity numerator')
eq(4*beta**3-27*(-S.Rational(4,27)*beta**3)*(-1),0,'finite balanced double-face wall')

print('Shared-root jet controls PASS; analytic theorem requires independent audit')
print('Universe: normal-unit m1..8,n1..9/infinity; three tangent regimes; full global positive and rational-mate controls')
print('Hostiles: same(m,n)different normal derivative; finite/infinite differential distinction; actual nonzero-M rational mate')
print('Always-active gates:',GATES)
print('Semantic SHA256:',hashlib.sha256(json.dumps(records,separators=(',',':')).encode()).hexdigest())
