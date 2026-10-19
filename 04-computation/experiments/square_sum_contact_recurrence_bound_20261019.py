#!/usr/bin/env python3
"""Controls for the logarithmic signed elliptic-contact index bound.

Exact rational/symbolic checks supply algebraic certificates. mpmath values
are explicitly diagnostic, not certified intervals. No census is attempted.
"""
from fractions import Fraction as F
from math import gcd
import json
import sys
import sympy as sp
import mpmath as mp

sys.set_int_max_str_digits(0)
CHECKS = {'exact_rational': 0, 'exact_symbolic': 0, 'numerical': 0}
mp.mp.dps = 80


def check(ok, message, kind='exact_rational'):
    CHECKS[kind] += 1
    if not ok:
        raise RuntimeError(message)


def spread(t):
    return 4*t*(t*t-1)/(t*t+1)**2


def on_curve(R,k):
    return R is None or R[1]**2 == (R[0]-k*k)*(R[0]**2-1)


def add(R,S,k):
    if R is None: return S
    if S is None: return R
    x,y = R
    u,v = S
    if x == u and y == -v: return None
    slope = (3*x*x-2*k*k*x-1)/(2*y) if R == S else (v-y)/(u-x)
    xx = slope*slope+k*k-x-u
    return xx,-y+slope*(x-xx)


def neg(R):
    return None if R is None else (R[0],-R[1])


def mul(n,R,k):
    if n < 0: return mul(-n,neg(R),k)
    out = None
    while n:
        if n & 1: out = add(out,R,k)
        R = add(R,R,k)
        n //= 2
    return out


def forward(u,z,k):
    x = (k*z+k*k+2*u+k*k*u*u)/(2*u*u)
    y = (x-k*k)*(x*u-1)/k
    check(on_curve((x,y),k),'forward point on cubic')
    return x,y


def inverse(R,k):
    if R is None: return None
    x,y = R
    if x == 0 or x == k*k: return None
    u = (k*y+x-k*k)/(x*(x-k*k))
    z = (2*u*u*x-k*k-2*u-k*k*u*u)/k
    check(z*z == k*k*(u*u+1)**2-4*u*(u*u-1), 'inverse quartic')
    return u,z


def normalize(t):
    n,m = t.numerator,t.denominator
    d = 2 if (n+m)%2 == 0 else 1
    a,b,c = (n*n-m*m)//d,2*n*m//d,(n*n+m*m)//d
    if a%2 == 0: a,b = b,a
    check(a>0 and b>0 and a*a+b*b == c*c,'positive PPT')
    check(gcd(a,b)==gcd(b,c)==gcd(a,c)==1,'primitive PPT')
    check(a%2==1 and b%2==0,'odd/even convention')
    check(F(2*a*b,c*c)==spread(t),'spread normalization')
    return (a,b,c),d


def rational_budget(f0,epsilon):
    D = 2+f0 if epsilon == 1 else F(2)
    J = 0
    while 2**J*f0 < D:
        J += 1
    check(2**J*f0 >= D,'exact power-of-two upper bound')
    check(J==0 or 2**(J-1)*f0<D,'least J')
    W = F(5,7)*(4+J)
    budget = 2*((20*W).numerator//(20*W).denominator)
    check(budget == 2*(100*(4+J)//7),'rational stopping formula')
    return J,W,budget


def anchor_data(t,epsilon):
    f0 = spread(t)
    check(t>1 and 0<f0<F(1,2),'admissible anchor')
    k = abs(t*t-1+epsilon*2*t)/(t*t+1)
    check(k*k==1+epsilon*f0,'rational K')
    P = (F(0),k)
    check(on_curve(P,k),'P on cubic')
    Q = forward(epsilon*t,t*t+1,k)
    check(inverse(Q,k)==(epsilon*t,t*t+1),'anchor map round trip')
    A = add(P,P,k)
    K = k*k
    check(A==(K+1/(4*K),(1-4*K*K)/(8*k**3)), 'exact 2P formula')
    check(A[0]>max(K,1) and A[1]<0,'2P identity lower branch')
    check(-1<P[0]<min(K,1),'P nonidentity component')
    return f0,k,P,Q,A


def witness(t,epsilon):
    f0,k,P,Q,A=anchor_data(t,epsilon)
    J,W,budget=rational_budget(f0,epsilon)
    R=None
    for n in range(1,min(budget//2,32)+1):
        R=add(R,A,k)
        check(on_curve(R,k),'even multiple on cubic')
        for sigma in (1,-1):
            S=add(Q,R if sigma==1 else neg(R),k)
            check(on_curve(S,k),'signed translate on cubic')
            pair=inverse(S,k)
            if pair is None: continue
            u,z=pair
            t1=epsilon*u
            if t1<=1 or z<=0: continue
            f1=spread(t1)
            if not (0<f1<F(1,2)) or f1==f0: continue
            root=z/(u*u+1)
            check(root*root==1+epsilon*f0-epsilon*f1,'exact signed contact')
            abc,d=normalize(t1)
            check(abc!=normalize(t)[0],'distinct primitive petal')
            q=2*sigma*n
            check(0<abs(q)<=budget and q%2==0,'exact stopping bound witness')
            return {'anchor_t':str(t),'epsilon':epsilon,'spread':str(f0),
                    'q':q,'direction':'up' if f1>f0 else 'down',
                    'J':J,'exact_rational_budget':budget,
                    'normalization_divisor':d,
                    'neighbor_hypotenuse_digits':len(str(abc[2]))}
    raise RuntimeError('selected exact witness not found in diagnostic cap')


def symbolic_controls():
    u,k,z=sp.symbols('u k z',nonzero=True)
    quartic=k*k*(u*u+1)**2-4*u*(u*u-1)
    X=(k*z+k*k+2*u+k*k*u*u)/(2*u*u)
    Y=(X-k*k)*(X*u-1)/k
    zp=sp.diff(quartic,u)/(2*z)
    Xp=sp.diff(X,u)+sp.diff(X,z)*zp
    num=sp.together(Xp*z+2*Y).as_numer_denom()[0]
    check(sp.rem(num,z*z-quartic,z)==0,'differential pullback', 'exact_symbolic')
    t=sp.symbols('t',positive=True)
    check(sp.factor(sp.diff(spread(t),t)+4*(t**4-6*t*t+1)/(t*t+1)**3)==0,
          'spread derivative','exact_symbolic')
    # sin(4 theta)=2 sin(2 theta) cos(2 theta), tan(theta)=t.
    sin4=2*(2*t/(1+t*t))*((1-t*t)/(1+t*t))
    check(sp.factor(spread(t)+sin4)==0,'angular spread','exact_symbolic')
    check(sp.factor(spread((t+1)/(t-1))-spread(t))==0,'swapped-leg spread','exact_symbolic')
    check(sp.factor(spread(-1/t)-spread(t))==0,'second real chart','exact_symbolic')
    D,m,v=sp.symbols('D m v',nonzero=True)
    x=D/v**2-1
    cubic=(x-(D-1))*(x-(m*D-1))*(x+1)
    ratio=sp.diff(x,v)**2/cubic
    check(sp.factor(ratio-4/(D*(1-v*v)*(1-m*v*v)))==0,
          'period substitution and factor two','exact_symbolic')
    K=sp.symbols('K',positive=True)
    x2=K+1/(4*K)
    check(sp.factor(x2-1-(K-sp.Rational(1,2))**2/K)==0,
          '2P identity component for K below one','exact_symbolic')
    check(sp.factor((1-4*K*K)**2/(64*K**3)-(x2-K)*(x2*x2-1))==0,
          '2P cubic identity','exact_symbolic')


def num(x):
    return mp.mpf(x.numerator)/x.denominator if isinstance(x,F) else mp.mpf(x)


def close(a,b,message,tol=mp.mpf('1e-55')):
    check(abs(a-b)<=tol*max(1,abs(a),abs(b)),message,'numerical')


def values(f0,epsilon):
    D=2+f0 if epsilon==1 else mp.mpf(2)
    r=mp.sqrt(f0/D)
    period=2*mp.ellipk(1-r*r)/mp.sqrt(D) # mpmath uses the parameter.
    integrand=lambda v:1/mp.sqrt((1+v*v)*(1+r*r*v*v))
    period_integral=2/mp.sqrt(D)*mp.quad(integrand,[0,1,1/r,mp.inf])
    close(period,period_integral,'ellipk period vs split integral')
    density=lambda s:1/(4*mp.sqrt(1-s*s)*mp.sqrt(1+epsilon*(f0-s)))
    down=mp.quad(density,[0,f0])
    up=mp.quad(density,[f0,mp.mpf('.5')])
    theta0=mp.pi/4+mp.asin(f0)/4
    theta_density=lambda theta:1/mp.sqrt(1+epsilon*(f0+mp.sin(4*theta)))
    close(down,mp.quad(theta_density,[mp.pi/4,theta0]),'decreasing arc variable change')
    close(up,mp.quad(theta_density,[theta0,7*mp.pi/24]),'increasing arc variable change')
    check(max(up,down)>=mp.pi/(24*mp.sqrt(6)),'global angular arc bound','numerical')
    check(up+down<period/2,'arcs under half period','numerical')
    elementary=(4*mp.asinh(1)+mp.log(D/f0))/mp.sqrt(D)
    global_elementary=(4*mp.asinh(1)+mp.log(mp.mpf('2.5')/f0))/mp.sqrt(2)
    check(period<=elementary<=global_elementary,'elementary period upper bounds','numerical')
    check(down>f0/5 and up>(mp.mpf('.5')-f0)/5,'rational arc lower bounds','numerical')
    return period,up,down


def numerical_controls():
    rows=[]
    for f in [F(1,1000),F(1,10),F(2,5),F(49,100)]:
        for epsilon in (1,-1):
            f0=num(f)
            period,up,down=values(f0,epsilon)
            J,W,budget=rational_budget(f,epsilon)
            check(period<num(W),'rational period upper bound','numerical')
            exact_arc_bound=2*int(mp.floor(period/max(up,down)))
            global_bound=2*int(mp.floor(24*mp.sqrt(6)/mp.pi*period))
            check(exact_arc_bound<=global_bound<=budget,'three stopping budgets ordered','numerical')
            K=1+epsilon*f0
            x2=K+1/(4*K)
            D=2+f0 if epsilon==1 else mp.mpf(2)
            alpha=mp.ellipf(mp.asin(mp.sqrt(D/(x2+1))),1-f0/D)/(mp.sqrt(D)*period)
            check(0<alpha<mp.mpf('.5'),'2P positive half-circle coordinate','numerical')
            N=exact_arc_bound//2
            pts=sorted((mp.frac(j*alpha),j) for j in range(N+1))
            gaps=[(pts[(i+1)%(N+1)][0]-pts[i][0]+(1 if i==N else 0),
                   pts[(i+1)%(N+1)][1]-pts[i][1]) for i in range(N+1)]
            gap,index=min(gaps)
            check(0<abs(index)<=N,'pigeonhole wrap index bound','numerical')
            check(0<gap<=1/mp.mpf(N+1) and gap*period<max(up,down),
                  'pigeonhole return enters longer arc','numerical')
            rows.append({'spread':str(f),'epsilon':epsilon,'period':mp.nstr(period,18),
                         'arc_up':mp.nstr(up,18),'arc_down':mp.nstr(down,18),
                         'longer_arc_bound':exact_arc_bound,'angular_bound':global_bound,
                         'rational_bound':budget,'rotation_2P':mp.nstr(alpha,18)})
    cusp=[]
    for epsilon in (1,-1):
        if epsilon==1:
            C=(mp.atanh(mp.sqrt(3)/2)-mp.atanh(1/mp.sqrt(2)))/(2*mp.sqrt(2))
        else:
            C=(mp.atanh(1/mp.sqrt(2))-mp.atanh(mp.mpf('.5')))/(2*mp.sqrt(2))
        close(C,mp.quad(lambda s:1/(4*mp.sqrt(1-s*s)*mp.sqrt(1-epsilon*s)),[0,mp.mpf('.5')]),
              'closed form limiting upper arc')
        for power in (2,4,8,12,24):
            f0=mp.mpf(10)**(-power)
            period,up,down=values(f0,epsilon)
            main=mp.log(32/f0)/mp.sqrt(2)
            # These are finite asymptotic controls, not the proof of the O terms.
            check(abs(period-main)<2*f0*mp.log(1/f0),'finite cusp period remainder','numerical')
            check(abs(up-C)<f0,'finite increasing arc remainder','numerical')
            check(abs(down/f0-mp.mpf('.25'))<f0,'finite decreasing arc remainder','numerical')
            cusp.append({'epsilon':epsilon,'spread_power10':-power,
                         'period':mp.nstr(period,18),
                         'period_minus_leading':mp.nstr(period-main,10),
                         'arc_up':mp.nstr(up,18),'arc_down_over_f':mp.nstr(down/f0,18),
                         'up_budget':2*int(mp.floor(period/up)),
                         'up_log_coefficient_limit':mp.nstr(mp.sqrt(2)/C,18)})
    return rows,cusp


def run():
    symbolic_controls()
    exact=[]
    for t in [F(9,7),F(8),F(10),F(25),F(100),F(6,5),F(11,10),F(101,100),F(1001,1000)]:
        for epsilon in (1,-1):
            exact.append(witness(t,epsilon))
    check(normalize(F(9,7))[1]==2,'odd/odd normalization control')
    check(normalize(F(9,7))[0]==normalize(F(8))[0],'swapped legs anchor equivalence')
    rows,cusp=numerical_controls()
    return {'status':'PROVED; numeric values are diagnostic, not interval certificates',
            'checks':CHECKS,'exact_selected_witnesses':exact,
            'numerical_period_and_arc_controls':rows,'numerical_cusp_controls':cusp}


if __name__=='__main__':
    result=run()
    print(json.dumps(result,indent=2,sort_keys=True))
