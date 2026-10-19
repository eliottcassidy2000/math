#!/usr/bin/env python3
"""Controls for the all-chart logarithmic cusp lower bound.

Exact algebra/rational tests are proof controls. mpmath integrations are
separately labelled diagnostics; they do not certify the theorem.
Explicit checks stay enabled under python -O.
"""
from fractions import Fraction as F
import json
import sys
import sympy as sp
import mpmath as mp

sys.set_int_max_str_digits(0)
CHECKS = 0


def check(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)


def zero(expr, message):
    check(sp.factor(expr) == 0, message)


def spread(u):
    return 4*u*(u*u-1)/(u*u+1)**2


def add(R,S,k):
    if R is None: return S
    if S is None: return R
    x,y=R; v,w=S
    if x==v and y==-w: return None
    slope=(3*x*x-2*k*k*x-1)/(2*y) if R==S else (w-y)/(v-x)
    xx=slope*slope+k*k-x-v
    return xx,-y+slope*(x-xx)


def on_curve(R,k):
    return R is None or R[1]**2==(R[0]-k*k)*(R[0]**2-1)


def inverse(R,k):
    if R is None or R[0] in (0,k*k): return None
    X,Y=R
    u=(k*Y+X-k*k)/(X*(X-k*k))
    z=(2*u*u*X-k*k-2*u-k*k*u*u)/k
    check(z*z==k*k*(u*u+1)**2-4*u*(u*u-1),'exact quartic point')
    return u,z


def exact_first(t,cap=24):
    f0=spread(t)
    check(0<f0<=F(1,1000),'probe lies in explicit cusp range')
    k=(t*t+2*t-1)/(t*t+1)
    Q=((t*t+4*t+1)/(t*t+1),8*t*(t+1)/(t*t+1)**2)
    P=(F(0),k)
    check(on_curve(Q,k),'anchor on cubic')
    R=None
    for q in range(1,cap+1):
        R=add(R,P,k)
        check(on_curve(R,k),'multiple on cubic')
        accepted=[]
        for sign in (1,-1):
            V=R if sign==1 else (R[0],-R[1])
            S=add(Q,V,k)
            check(on_curve(S,k),'translate on cubic')
            pair=inverse(S,k)
            if pair is None: continue
            u,z=pair
            f=spread(u)
            check((z/(u*u+1))**2==1+f0-f,'contact identity')
            if 0<f<F(1,2) and f!=f0:
                check(u>1 or -1<u<0,'all positive-spread charts included')
                accepted.append(sign*q)
        if accepted: return q,accepted
    raise RuntimeError('exact probe exhausted its diagnostic cap')


def numerical_data(delta):
    K=1+delta
    Om=2/mp.sqrt(K+1)*mp.ellipk(2/(K+1))
    lo=1/(4*K)
    a=mp.quad(lambda v:1/(2*mp.sqrt(v*(v+delta)*(v+2+delta))),[lo,1,mp.inf])
    C=mp.pi/4+mp.pi*mp.sqrt(2)/12
    L=mp.quad(lambda th:1/mp.sqrt(K+mp.sin(4*th)),
              [-mp.pi/24,0,mp.pi/4,7*mp.pi/24])
    return Om,a,C,L


def main():
    t,k,v=sp.symbols('t k v',positive=True)
    K=k*k
    zero(((t*t+2*t-1)/(t*t+1))**2-1-spread(t),'K-1=f(t)')
    X=K+1/(4*K)
    Y=(1-4*K*K)/(8*K*k)
    zero(Y*Y-(X-K)*(X*X-1),'2P lies on cubic')
    slope=-1/(2*k)
    zero(slope*slope+K-X,'doubling X')
    zero(-k-slope*X-Y,'doubling Y')
    up=(2*t**3+t*t-2*t+1)/(t*(t*t-4*t+1))
    um=(t*t+4*t+1)/(t**3+2*t*t+t-2)
    zero(up+1-(t-1)*(3*t*t-1)/(t*(t*t-4*t+1)),'odd plus base sign identity')
    zero(1-um-(t+1)*(t*t-3)/(t**3+2*t*t+t-2),'odd minus base sign identity')
    zero(spread(-1/t)-spread(t),'negative accepted chart symmetry')
    zero(spread(t)-(1-(t*t-2*t-1)**2/(t*t+1)**2),'global maximum f<=1')
    check(spread(F(4,3))==F(336,625)>F(1,2),'two anchor ranges separated')
    theta=sp.symbols('theta',real=True)
    zero(sp.trigsimp(spread(sp.tan(theta))+sp.sin(4*theta),method='fu'),'angular spread')
    check(F(565001,562500)==(1+F(4,1000))*(1+F(4,9000)),'explicit denominator endpoint')
    check(F(565001,562500)<F(1003,1000)**2,'explicit denominator bound')
    check(F(3500,3009)>F(341,294),'step strictly larger than cluster')
    check(F(7,5)**2<2<F(10,7)**2,'rational square-root enclosures')
    check(F(341,294)==F(22,7)*(3+F(10,7))/12,'cluster upper constant')
    j=sp.symbols('j',positive=True)
    tl=(j+2)/(j+1); tr=2*j+2
    zero((j+2)**2+(j+1)**2-(2*j*j+6*j+5),'left spine hypotenuse')
    zero((2*j+2)**2+1-(4*j*j+8*j+5),'right spine hypotenuse')
    check(sp.limit(j*spread(tl),j,sp.oo)==2,'left spine spread asymptotic')
    check(sp.limit(j*spread(tr),j,sp.oo)==2,'right spine spread asymptotic')
    check(0<spread(F(5,4))<F(1,2),'left spine admitted at depth three')
    check(0<spread(F(8))<F(1,2),'right spine admitted at depth three')
    exact_algebra_checks=CHECKS

    mp.mp.dps=90
    a0=mp.sqrt(2)*mp.asinh(1)
    beta=1/(2*mp.asinh(1))
    numeric=[]
    for power in (3,6,12,24,40):
        d=mp.mpf(10)**(-power)
        Om,a,C,L=numerical_data(d)
        check(a>C>L,'diagnostic: translating step exceeds actual cluster')
        check(a>=a0/mp.sqrt((1+4*d)*(1+4*d/9)),'diagnostic: integral step lower bound')
        check(abs(a-a0)<d,'diagnostic: alpha cusp approach')
        check(abs(Om-mp.log(32/d)/mp.sqrt(2))<d*mp.log(1/d),'diagnostic: period cusp approach')
        numeric.append({'delta':str(d),'period':mp.nstr(Om,25),
                        'alpha':mp.nstr(a,25),'cluster_actual':mp.nstr(L,25),
                        'all_index_lower_bound':mp.nstr((Om-2*C)/a-1,25)})

    witnesses=[]
    for t in (F(10000),F(10**6),F(10**9),F(10**6+1,10**6)):
        m,signs=exact_first(t)
        d=mp.mpf(spread(t).numerator)/spread(t).denominator
        Om,a,C,L=numerical_data(d)
        bound=(Om-2*C)/a-1
        check(m>=bound,'exact first index above numerical evaluation of lower bound')
        witnesses.append({'anchor':str(t),'exact_first_abs_q':m,'signs':signs,
                          'numerical_lower_bound':mp.nstr(bound,25)})
    return {'exact_algebra_checks':exact_algebra_checks,
            'total_explicit_checks':CHECKS,
            'leading_lower_coefficient':mp.nstr(beta,30),
            'diagnostic_integrals':numeric,'exact_first_success_probes':witnesses,
            'scope':'Quadrature values are diagnostics, not proof certificates.'}


if __name__=='__main__':
    result=main()
    text=json.dumps(result,indent=2)+'\n'
    if len(sys.argv)==2:
        with open(sys.argv[1],'w') as file: file.write(text)
    print(text,end='')
