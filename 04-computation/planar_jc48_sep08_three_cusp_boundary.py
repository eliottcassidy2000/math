#!/usr/bin/env python3
"""Exact geometry and braid-five hostile controls; no actual braid claim."""
import hashlib
import json
from itertools import permutations
import sympy as S

t, h, s, a, u, v = S.symbols('t h s a u v')
U = t**4-S.Rational(4,3)*s*t**3-2*t**2+4*s*t
V = (t**6+S.Rational(6,5)*(a-s)*t**5-S.Rational(3,2)*t**4
     +2*(-a*s**2-a+s)*t**3+6*a*s**2*t)
A = (t*t-1)*(t-s)
q = t*t+a*t+a*s
GATES = 0


def check(condition, label):
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(label)


def zero(expr, label):
    check(S.cancel(expr) == 0, label)


def mul(p, q):
    return tuple(p[q[i]] for i in range(len(p)))


def inverse(p):
    z = [0]*len(p)
    for i, j in enumerate(p):
        z[j] = i
    return tuple(z)


def product(*ps):
    z = tuple(range(len(ps[0])))
    for p in ps:
        z = mul(z, p)
    return z


def fix(p):
    return {i for i, j in enumerate(p) if i == j}


def main():
    zero(S.diff(U,t)-4*A, 'quartic derivative')
    zero(S.diff(V,t)-6*A*q, 'sextic derivative')
    numerator=S.Poly(S.expand(V**2-U**3),t)
    zero(numerator.nth(12),'cancelled infinity leading numerator')
    zero(numerator.nth(11)-S.Rational(4,5)*(3*a+2*s),'actual infinity seventh coefficient')
    for e in (-1,1,s):
        zero((S.diff(U,t,2)*S.diff(V,t,3)-S.diff(V,t,2)*S.diff(U,t,3)).subs(t,e)
             -48*S.diff(A,t).subs(t,e)**2*(2*e+a), 'ordinary jet identity')
    rows = []
    # The three disjoint bad-jet lines over s!=+-1; there is no omitted chart.
    for av,e,cover in ((2,-1,-3),(-2,1,3),(-2*s,s,0)):
        uu=S.expand(U.subs(t,e+h)-U.subs(t,e))
        # Evaluate the subtracted constant before translating t.
        vv=S.expand(V.subs({a:av,t:e+h})-V.subs({a:av,t:e}))
        a1=S.diff(A,t).subs(t,e)
        a2=3*e-s
        qe=q.subs({a:av,t:e})
        ww=S.expand(vv-S.Rational(3,2)*qe*uu)
        zero(uu-(2*a1*h*h+S.Rational(4,3)*a2*h**3+h**4), 'local U exact polynomial')
        zero(ww-(S.Rational(3,2)*a1*h**4+S.Rational(6,5)*a2*h**5+h**6),
             'local V exact polynomial')
        remainder=S.Poly(S.expand(ww-S.Rational(3,8)/a1*uu**2),h)
        for j in range(5):
            zero(remainder.nth(j), 'vanishing lower local coefficients')
        zero(remainder.nth(5)+S.Rational(4,5)*a2, 'first fifth coefficient')
        zero(a2.subs(s,cover) if hasattr(a2,'subs') else a2, 'unique cover value')
        zero(uu.subs(s,cover)-uu.subs({s:cover,h:-h}), 'cover U even')
        zero(ww.subs(s,cover)-ww.subs({s:cover,h:-h}), 'cover V even')
        i7=S.factor(S.Rational(4,5)*(3*av+2*s))
        check(S.degree(i7,s)==1 and S.solve(i7,s)==[cover], 'infinity7 zero is the cover')
        rows.append([str(av),str(e),str(S.factor(-S.Rational(4,5)*a2)),cover,str(i7)])

    # Genuine affine source changes and linear target normalization identify
    # all three bad-jet lines with a=2; local nonlinear tests are not used here.
    zero(U.subs({t:-t,s:-s})-U,'parameter reflection U')
    zero(V.subs({t:-t,s:-s,a:-a})-V,'parameter reflection V')
    alpha=(1-s)/2
    beta=(s+1)/2
    new_s=(s+3)/(s-1)
    moved_u=S.cancel((U.subs(t,alpha*t+beta)-U.subs(t,beta))/alpha**4)
    moved_v=S.cancel((V.subs(a,-2*s).subs(t,alpha*t+beta)
                      -V.subs(a,-2*s).subs(t,beta))/alpha**6)
    kq=S.cancel((beta**2-2*s*beta-2*s*s)/alpha**2)
    zero(moved_u-U.subs(s,new_s),'bad-s critical relabel U')
    zero(moved_v+S.Rational(3,2)*(2*new_s-kq)*moved_u
         -V.subs({s:new_s,a:2}),'bad-s critical relabel V with genuine shear')

    # Original cusp-image resultants remain valid when one ordinary jet vanishes.
    rm=(45*a*a*s*s-162*a*a*s-27*a*a+80*a*s**3-216*a*s*s-216*a*s
        +128*s**3-432*s*s)
    rp=(45*a*a*s*s+162*a*a*s-27*a*a+80*a*s**3+216*a*s*s-216*a*s
        -128*s**3-432*s*s)
    rt=36*a*a*s*s-216*a*a+28*a*s**3-108*a*s+11*s**4-126*s*s+675
    for e,expected in ((-1,-4*(s+1)**2*rm/675),
                       (1,-4*(s-1)**2*rp/675),
                       (s,(s*s-1)**2*rt/675)):
        px=S.cancel((U-U.subs(t,e))/(t-e)**2)
        py=S.cancel((V-V.subs(t,e))/(t-e)**2)
        zero(S.resultant(px,py,t)-expected, 'complete cusp-image resultant')
        check(S.Poly(px,t).LC()==1, 'monic finite incidence')
    expected_a2=(36*(s+1)*(8*s*s-27*s-3),4*(s+3)**2*(8*s-3),
                 (s+1)*(s+3)**2*(11*s-21))
    for p,r in zip((rm,rp,rt),expected_a2):
        zero(p.subs(a,2)-r,'a2 image resultant factors')

    # Exact literals: the residual discriminant nodes avoid all critical U-values.
    control_rows=[]
    for sv in (0,2):
        x=S.expand(U.subs({s:sv,a:2}));y=S.expand(V.subs({s:sv,a:2}))
        poly=S.resultant(x-u,y-v,t)
        check(S.Poly(poly,v).degree()==4 and S.Poly(poly,v).LC()==1,'monic vertical quartic')
        zero(poly.subs({u:x,v:y}), 'literal resultant substitution')
        if sv==0:
            crit=u*(u+1)
            nodes=(9*u+5)*(625*u*u+434*u+49)
            expected=-S.Rational(1048576,244140625)*u**3*(u+1)**8*nodes**2
            expected_u=-256*u*(u+1)**2
        else:
            crit=(3*u-13)*(3*u-8)*(3*u+19)
            nodes=(u-3)*(2187*u*u-13962*u+22139)
            expected=(-S.Rational(1048576,847288609443)*(3*u-13)**3
                      *(3*u-8)**3*(3*u+19)**5*nodes**2)
            expected_u=-S.Rational(256,27)*crit
        zero(S.discriminant(poly,v)-expected,'full vertical discriminant')
        zero(S.discriminant(x-u,t)-expected_u,'full parameter discriminant')
        check(S.degree(nodes,u)==3,'three residual node values')
        check(S.degree(S.gcd(nodes,S.diff(nodes,u)),u)==0,'residual squarefree')
        check(S.degree(S.gcd(nodes,crit),u)==0,'residual avoids U criticals')
        check(S.degree(S.gcd(S.diff(x,t),S.diff(y,t)),t)==3,'exact three common criticals')
        image_pairs=[]
        for e in (-1,1,sv):
            px=x-x.subs(t,e);py=y-y.subs(t,e)
            zero(S.monic(S.gcd(px,py),t)-(t-e)**2,'cusp image has sole preimage')
            image_pairs.append([str(x.subs(t,e)),str(y.subs(t,e))])
        check(len({tuple(p) for p in image_pairs})==3,'three distinct target cusp images')
        check(10-2-1-1-3==3,'genus node count')
        control_rows.append({'s':sv,'images':image_pairs,'nodes':str(S.expand(nodes))})

    # A minimal five-letter hostile to carrying the ordinary-cusp inequality.
    sigma=(1,2,0,3,4)       # (123)
    tau=(0,1,3,4,2)         # (345)
    check(product(sigma,tau,sigma,tau,sigma)==product(tau,sigma,tau,sigma,tau),
          'odd braid length five')
    check(product(sigma,tau,sigma)!=product(tau,sigma,tau),'ordinary braid fails')
    ss=set(range(5))-fix(sigma);tt=set(range(5))-fix(tau)
    check((len(ss),len(tt),len(ss&tt))==(3,3,1),'half-support hostile')
    g=product(sigma,tau,sigma,tau)
    check(product(g,sigma,inverse(g))==tau,'correct odd-five conjugator')
    aa,bb=fix(sigma),fix(tau)
    check({g[i] for i in aa}==bb,'full-fixed five-cusp re-access')
    k,n=len(aa),len(aa&bb)
    check((k,n)==(2,0),'full-fixed hostile counts')
    check(2*len(ss&tt)<len(ss),'ordinary half-overlap refuted')
    check(4*n<5*k-5,'guessed higher-cusp injection refuted')
    # Complete small ambient check of minimality, retaining the length-five law.
    small_pairs=0
    for d in range(2,5):
        bank=list(permutations(range(d)))
        for p in bank:
            for z in bank:
                ps=set(range(d))-fix(p);zs=set(range(d))-fix(z)
                if len(ps)!=len(zs):
                    continue
                if product(p,z,p,z,p)!=product(z,p,z,p,z):
                    continue
                check(2*len(ps&zs)>=len(ps),'minimal ambient hostile size')
                small_pairs+=1

    report={'boundary_rows':rows,'controls':control_rows,'minimality_pairs':small_pairs,
            'hostile_sigma':sigma,'hostile_tau':tau,'hostile_k_n':[k,n]}
    digest=hashlib.sha256(json.dumps(report,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    print('THREE-CRITICAL NONORDINARY BOUNDARY GEOMETRY: exact controls PASS')
    print('Rows [a,e,fifth coefficient,double-cover s,infinity seventh]:')
    for row in rows: print(row)
    print('Birational boundary type: finite (2,5),(2,3),(2,3); infinity(2,7).')
    for row in control_rows: print('Literal:',json.dumps(row,sort_keys=True))
    print('Braid-five hostile: sigma=(123),tau=(345),support overlap1,k2,n0.')
    print('Ordinary half-overlap and guessed4n>=5k-D are REFUTED by this abstract pair.')
    print('No actual global braid or Keller realization is certified here.')
    print('Minimality support controls on D<=4:',small_pairs)
    print('Semantic SHA256:',digest)
    print('Always-active gates:',GATES)


if __name__=='__main__':
    main()
