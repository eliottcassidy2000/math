#!/usr/bin/env python3
"""Exact two-ordinary-cusp (4,6) class controls; no finite parameter census.

Default is exact algebra. --scout is a labelled numerical infinity11 sidecar;
its words are not an actual braid certificate.
"""
import json,sys
import sympy as S

t,x,z,s,a,b,c,p,q,A,B=S.symbols('t x z s a b c p q A B')
GATES=0
def check(value,label):
    global GATES
    GATES+=1
    if not value:raise RuntimeError(label)
def eq(lhs,rhs,label):check(S.cancel(lhs-rhs)==0,label)
def monic_gcd(lhs,rhs):return S.Poly(S.gcd(lhs,rhs),t).monic().as_expr()

U=t**4-4*s*t**3/3-2*t*t+4*s*t
Vb=t**6+6*a*t**5/5+S.Rational(3,2)*(b-1)*t**4+2*(c-a)*t**3-3*b*t*t-6*c*t
V=Vb.subs(b,0)
P=t**3+a*t*t+c

def leading(expr,n,coefficient,label):
    numerator,denominator=S.together(expr).as_numer_denom()
    polynomial=S.Poly(S.expand(numerator),z)
    check(all(k[0]>=n for k,val in polynomial.terms()),label+' lower coefficients')
    d0=S.expand(denominator).subs(z,0)
    check(d0!=0,label+' denominator at zero')
    eq(polynomial.coeff_monomial(z**n),coefficient*d0,label+' coefficient')

def classification():
    eq(S.diff(U,t),4*(t*t-1)*(t-s),'complete quartic derivative')
    eq(S.diff(Vb,t),6*(t*t-1)*(t**3+a*t*t+b*t+c),'complete sextic derivative')
    eq(Vb-3*b*U/2,V.subs(c,c+b*s),'degree-preserving target shear')
    eq(S.rem(P,t-s,t),P.subs(t,s),'exact extra critical condition')
    jet=S.diff(U,t,2)*S.diff(V,t,3)-S.diff(V,t,2)*S.diff(U,t,3)
    for e in [-1,1]:
        eq(jet.subs(t,e),192*((e-s)*S.diff(P,t).subs(t,e)-P.subs(t,e)),
           'ordinary-cusp determinant including repeated projection critical point')
    eq(S.discriminant(S.diff(U,t),t),1024*(s*s-1)**2,'only double projection-critical boundaries')
    for e in [-1,1]:
        eq(monic_gcd(S.diff(U,t).subs(s,e),S.diff(U,t,2).subs(s,e)),t-e,'no triple projection-critical root')
    # Local hostile: s=1 still has an ordinary cusp at t=1.
    sub={s:1,a:0,c:S.Rational(1,6)}
    eq(S.diff(U,t,2).subs(sub).subs(t,1),0,'order-two U coefficient may vanish')
    eq(S.diff(V,t,2).subs(sub).subs(t,1),14,'V supplies the order-two cusp coordinate')
    eq(jet.subs(sub).subs(t,1),-224,'ordinary cusp survives the repeated U critical root')
    eq(monic_gcd(S.diff(U,t).subs(sub),S.diff(V,t).subs(sub)),t*t-1,'repeated U root does not add source ramification')

    X=S.cancel(U.subs(t,1/z)/V.subs(t,1/z));Z=S.cancel(1/V.subs(t,1/z))
    leading(X,2,1,'infinity first coordinate')
    leading(Z,6,1,'line contact six')
    leading(Z-X**3,7,4*(3*a+5*s)/5,'first infinity odd coefficient')
    a9=-5*s/3
    X9=X.subs(a,a9);Z9=Z.subs(a,a9)
    k4=3-4*s*s/3
    h9=4*(c-4*s/3-14*s**3/27)
    leading(Z9-X9**3,8,k4,'infinity ninth even removal')
    leading(Z9-X9**3-k4*X9**4,9,h9,'infinity ninth odd coefficient')
    c11=4*s/3+14*s**3/27
    X11=X9.subs(c,c11);Z11=Z9.subs(c,c11)
    k5=(128*s**4-216*s*s+1053)/108
    h11=-16*s*(s-3)**2*(s+3)**2/81
    leading(Z11-X11**3-k4*X11**4,10,k5,'infinity eleventh even removal')
    leading(Z11-X11**3-k4*X11**4-k5*X11**5,11,h11,'infinity eleventh odd coefficient')
    ps=S.factor(P.subs({a:a9,c:c11}).subs(t,s))
    eq(ps,-4*s*(s-3)*(s+3)/27,'higher infinity extra critical factor')
    eq(h11,4*(s*s-9)*ps/3,'infinity thirteen forces an excluded critical condition')
    for e in [-1,1]:
        local=((e-s)*S.diff(P,t).subs(t,e)-P.subs(t,e)).subs({a:a9,c:c11})
        eq(local,-2*(s-3*e)**2*(7*s-3*e)/27,'ordinary cusp factors on infinity11 stratum')
    for m,n in [(7,5),(9,4),(11,3)]:check(2+(m-1)//2+n==10,'complete two-cusp sextic genus')
    even={s:0,a:0,c:0}
    eq(U.subs(even).subs(t,-t),U.subs(even),'even quadratic-cover hostile U')
    eq(V.subs(even).subs(t,-t),V.subs(even),'even quadratic-cover hostile V')
    eq(monic_gcd(S.diff(U,t).subs(even),S.diff(V,t).subs(even)),t*(t*t-1),'excluded degree-two cover has odd derivative gcd degree')
    eq(k4.subs(s,S.Rational(3,2)),0,'vanishing intermediate even coefficient retained')
    return dict(infinity7=str(4*(3*a+5*s)/5),infinity9=str(h9),infinity11=str(h11))

def geometry_control(name,sub,nodes):
    u,v=S.symbols('u v');u0=S.expand(U.subs(sub));v0=S.expand(V.subs(sub))
    eq(monic_gcd(S.diff(u0,t),S.diff(v0,t)),t*t-1,name+' exact critical locus')
    for e in [-1,1]:
        eq(monic_gcd(u0-u0.subs(t,e),v0-v0.subs(t,e)),(t-e)**2,name+' single cusp preimage')
        check((S.diff(u0,t,2)*S.diff(v0,t,3)-S.diff(v0,t,2)*S.diff(u0,t,3)).subs(t,e)!=0,name+' ordinary cusp')
    N=S.cancel((u0.subs(t,x)-u0)/(x-t));M=S.cancel((v0.subs(t,x)-v0)/(x-t))
    Np=S.rem(N.subs(t,p-x),x*x-p*x+q,x)
    Mp=S.rem(M.subs(t,p-x),x*x-p*x+q,x)
    if name=='infinity7':
        eq(Np,p*(p*p-2*q-2),'split first pair equation')
        Hq=18*q*q+25*q-15;Hp=15*p**3+18*p*p-22
        eq(Mp.subs(p,0),Hq/15,'both zero-sum node pairs retained')
        eq(Mp.subs(q,(p*p-2)/2),-(p*p-4)*Hp/60,'three nonzero-sum node pairs')
        check(S.discriminant(Hq,q)!=0 and S.resultant(Hq,q*(q+1),q)!=0,'zero-sum pairs distinct and off diagonal')
        check(S.discriminant(Hp,p)!=0 and S.resultant(Hp,p*(p*p-4),p)!=0,'other pairs distinct and off diagonal')
        check(S.degree(Hq,q)+S.degree(Hp,p)==nodes,'five-pair count')
    else:
        qp=(3*p**3-8*p*p-6*p+24)/(2*(3*p-4));Hp=9*p**3-60*p*p+112*p-144
        eq(Np,(3*p**3-8*p*p-6*p+24-2*(3*p-4)*q)/3,'residual first pair equation')
        eq(Np.subs(p,S.Rational(4,3)),S.Rational(80,27),'residual pair denominator cannot vanish')
        eq(Mp.subs(q,qp),-(p*p-4)*Hp/36,'three residual node pairs')
        eq(p*p-4*qp,-3*(p-4)*(p*p-4)/(3*p-4),'residual pair discriminant')
        check(S.discriminant(Hp,p)!=0 and S.resultant(Hp,(3*p-4)*(p-4)*(p*p-4),p)!=0,'residual pairs distinct and off diagonal')
        check(S.degree(Hp,p)==nodes,'three-pair count')
    tangent=S.diff(u0,t).subs(t,x)*S.diff(v0,t)-S.diff(v0,t).subs(t,x)*S.diff(u0,t)
    check(S.groebner([N,M,tangent],x,t,domain=S.QQ)==S.groebner([x+t**3-2*t,(t*t-1)**2],x,t,domain=S.QQ),name+' only cusp diagonal tangencies')
    R=S.rem(v0-B,u0-A,t)
    check(S.degree(R,t)==3 and not S.Poly(R,t).LC().free_symbols,name+' constant cubic leading coefficient')
    coeffs=[S.together(k).as_numer_denom()[0] for k in S.Poly(S.rem(u0-A,R,t),t).all_coeffs()]
    check(S.groebner(coeffs,A,B,domain=S.QQ)==S.groebner([1],A,B,domain=S.QQ),name+' no triple image')
    F=S.resultant(u0-u,v0-v,t)
    eq(F.subs({u:u0,v:v0}),0,name+' actual resultant substitution')
    check(S.Poly(F,v).LC()==1 and S.degree(F,v)==4,name+' monic actual quartet')
    return dict(name=name,U=str(u0),V=str(v0),nodes=nodes,discriminant=str(S.factor(S.discriminant(F,v))))

def scout():
    import numpy as np,itertools
    choices=list(itertools.permutations(range(4)));mult=1+.25j;base=4+3j;r=1/32
    centres=[('smooth',8/3+0j),('cusp_plus',13/3+0j),('cusp_minus',-19/3+0j),
      ('node_real',-8.42611211+0j),('node_upper',-8.71081226+3.75245026j),
      ('node_lower',-8.71081226-3.75245026j)]
    def roots(u):
        ts=np.roots([1,-8/3,-2,8,-u])
        return [t**6-4*t**5-1.5*t**4+548*t**3/27-368*t/9 for t in ts]
    for steps in [256,512]:
        for name,centre in centres:
            vv=[base,centre+r,centre+1j*r,centre-r,centre-1j*r,centre+r,base]
            old=sorted(roots(base),key=lambda z:(z*mult).real);order=list(range(4));word=[]
            for start,end in zip(vv,vv[1:]):
                for j in range(1,steps+1):
                    raw=roots(start+(end-start)*j/steps)
                    perm=min(choices,key=lambda pp:sum(abs(raw[pp[k]]-old[k])**2 for k in range(4)))
                    new=[raw[perm[k]] for k in range(4)];z0=[z*mult for z in old];z1=[z*mult for z in new];events=[]
                    for i in range(4):
                        for k in range(i):
                            aa=(z0[i]-z0[k]).real;bb=(z1[i]-z1[k]).real
                            if aa*bb<0:events.append((aa/(aa-bb),i,k))
                    for h,i,k in sorted(events):
                        ii,jj=order.index(i),order.index(k)
                        if abs(ii-jj)!=1:raise RuntimeError('heuristic nonadjacent crossing')
                        pos=min(ii,jj);left,right=order[pos:pos+2]
                        imag=((1-h)*(z0[left]-z0[right])+h*(z1[left]-z1[right])).imag
                        if imag==0:raise RuntimeError('heuristic zero crossing separation')
                        letter=(pos+1)*(1 if imag<0 else-1)
                        if word and word[-1]==-letter:word.pop()
                        else:word.append(letter)
                        order[pos:pos+2]=[right,left]
                    old=new
            print('HEURISTIC ONLY',steps,name,word,flush=True)

def main():
    if '--scout' in sys.argv:scout();return
    coefficients=classification()
    controls=[geometry_control('infinity7',{s:0,a:1,c:S.Rational(1,6)},5),
      geometry_control('infinity11',{s:2,a:-S.Rational(10,3),c:S.Rational(184,27)},3)]
    print('FINITE-EXACT TWO-CUSP FAMILY ALGEBRA PASS; topology proved separately; infinity11 group OPEN')
    print(json.dumps(dict(coefficients=coefficients,controls=controls,gates=GATES),sort_keys=True,indent=2))
    print('PASS always-active gates='+str(GATES))
if __name__=='__main__':main()
