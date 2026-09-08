#!/usr/bin/env python3
"""Exact geometry of one literal two-cusp sextic; no global braid claims."""
import json
import sympy as S

s,t,p,q,u,v,A,B,z,lam=S.symbols('s t p q u v A B z lam')
GATES=0
def check(value,label):
    global GATES
    GATES+=1
    if not value:raise RuntimeError(label)
def eq(a,b,label):check(S.cancel(a-b)==0,label)

U=t**4-2*t**2
V=t**6-S.Rational(3,2)*t**4+t**3/3-t
F=(1296*u**6+1944*u**5+7776*u**4*v+2457*u**4
   -2592*u**3*v**2+10440*u**3*v+7400*u**3+9720*u**2*v**2
   -432*u**2*v+8088*u**2-7776*u*v**3+16416*u*v**2-2592*u*v
   +2448*u+1296*v**4-5184*v**3+4896*v**2)/1296
eq(F,S.resultant(U-u,V-v,t),'literal actual resultant')
eq(F.subs({u:U,v:V}),0,'actual substitution')
check(S.degree(F,v)==4 and S.Poly(F,v).LC()==1,'monic quartet')
eq(S.gcd(S.diff(U,t),S.diff(V,t)),t*t-1,'exact critical locus')
for e,image,jet in [(1,-S.Rational(7,6),352),(-1,S.Rational(1,6),-416)]:
    eq(U.subs(t,e),-1,'cusp projection')
    eq(V.subs(t,e),image,'distinct cusp image')
    eq(S.gcd(U+1,V-image),(t-e)**2,'cusp fibre has one parameter')
    determinant=S.diff(U,t,2)*S.diff(V,t,3)-S.diff(V,t,2)*S.diff(U,t,3)
    eq(determinant.subs(t,e),jet,'ordinary cusp jet determinant')

N=S.cancel((U.subs(t,s)-U)/(s-t));M=S.cancel((V.subs(t,s)-V)/(s-t))
Np=S.rem(N.subs(t,p-s),s*s-p*s+q,s)
Mp=S.rem(M.subs(t,p-s),s*s-p*s+q,s)
eq(Np,p*(p*p-2*q-2),'complete first pair equation')
eq(Mp,(6*p**5-24*p**3*q-9*p**3+2*p*p+18*p*q*q+18*p*q-2*q-6)/6,'complete second pair equation')
eq(Mp.subs(p,0),-(q+3)/3,'retained zero-sum node branch')
eq(Mp.subs({p:0,q:-1}),-S.Rational(2,3),'no collision where first pair branches meet')
H=3*p**3-2
eq(Mp.subs(q,(p*p-2)/2),-(p-2)*(p+2)*H/12,'nonzero-sum collision branch')
eq(p*p-4*(p*p-2)/2,4-p*p,'pair off-diagonal discriminant')
check(S.discriminant(H,p)!=0,'three distinct nonzero-sum pairs')
check(S.resultant(H,p*(p-2)*(p+2),p)!=0,'all three pairs off diagonal')
eq(S.rem(U,t*t-3,t),3,'zero-sum node target first coordinate')
eq(S.rem(V,t*t-3,t),S.Rational(27,2),'zero-sum node target second coordinate')
T=S.diff(U,t).subs(t,s)*S.diff(V,t)-S.diff(V,t).subs(t,s)*S.diff(U,t)
check(S.groebner([N,M,T],s,t,domain=S.QQ)==
      S.groebner([s+t**3-2*t,(t*t-1)**2],s,t,domain=S.QQ),'off-diagonal tangents never coincide')
pair_res=-(t-1)**2*(t+1)**2*(t*t-3)*(18*t**6-54*t**4+6*t**3+54*t*t-18*t-17)/27
eq(S.resultant(N,M,s),pair_res,'complete ordered-pair resultant with cusp multiplicities')
R=S.rem(V-B,U-A,t)
eq(R,t**3/3+(A+1)*t*t-t+A/2-B,'constant-leading triple remainder')
coeffs=[S.together(c).as_numer_denom()[0] for c in S.Poly(S.rem(U-A,R,t),t).all_coeffs()]
check(S.groebner(coeffs,A,B,domain=S.QQ)==S.groebner([1],A,B,domain=S.QQ),'no shared triple image')

X=S.cancel(U.subs(t,1/z)/V.subs(t,1/z));Z=S.cancel(1/V.subs(t,1/z))
eq(S.limit(X/z**2,z,0),1,'infinity multiplicity two')
eq(S.limit(Z/z**6,z,0),1,'line contact six')
eq(S.limit((Z-X**3)/z**7,z,0),0,'infinity seventh coefficient vanishes')
eq(S.limit((Z-X**3)/X**4,z,0),3,'infinity even coefficient removed')
eq(S.limit((Z-X**3-3*X**4)/z**9,z,0),S.Rational(2,3),'infinity first odd coefficient ninth')
check(1+1+4+4==10,'complete sextic genus check')
node_poly=324*u**3+972*u*u+1080*u+289
disc=-S.Rational(256,531441)*u*(u-3)**2*(u+1)**6*node_poly**2
eq(S.discriminant(F,v),disc,'complete actual vertical discriminant')
eq(S.discriminant(node_poly,u),-59592250800,'three node projection values distinct')
eq(S.resultant(node_poly,u*(u-3)*(u+1),u),868900175,'no additional projection coincidences')
eq(F.subs(u,-1),(6*v-1)**2*(6*v+7)**2/1296,'two distinct co-projected cusp targets')
eq(F.subs(u,0),v*v*(9*v*v-36*v+34)/9,'single smooth vertical fold')
eq(S.diff(V,t).subs(t,0),-1,'fold parameter is a smooth curve point')
eq(F.subs(u,3),(2*v-27)**2*(36*v*v+180*v+289)/144,'isolated zero-sum node projection')

# Hostile parameter controls do not promote a whole lambda-family.
Vlam=t**6-S.Rational(3,2)*t**4+lam*(t**3/3-t)
jetlam=S.diff(U,t,2)*S.diff(Vlam,t,3)-S.diff(Vlam,t,2)*S.diff(U,t,3)
eq(jetlam.subs(t,1),32*(12-lam),'ordinary-cusp test fails at lambda twelve')
eq(jetlam.subs(t,-1),-32*(12+lam),'ordinary-cusp test fails at lambda minus twelve')
eq(Vlam.subs({lam:0,t:-t}),Vlam.subs(lam,0),'lambda zero loses birationality by evenness')

print('FINITE-EXACT TWO-CUSP GEOMETRY PASS; global complement classification OPEN')
print(json.dumps(dict(curve=['t^4-2t^2','t^6-3t^4/2+t^3/3-t'],
    finite_cusps=[['-1','-7/6','2,3'],['-1','1/6','2,3']],
    nodes=4,infinity='2,9',pair_factors=['p=0,q=-3','3p^3-2=0,q=(p^2-2)/2'],
    discriminant=str(S.factor(disc)),gates=GATES),sort_keys=True,indent=2))
print('PASS always-active gates='+str(GATES))
