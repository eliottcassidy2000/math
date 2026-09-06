#!/usr/bin/env python3
"""Exact controls for depth-lowering universal Hamiltonian carriers.

No repository producer is imported. Universal conclusions use the proof;
this certificate tests literal source brackets, both depth mechanisms,
core D divisibility, and finite scalar exponentials in the log chart.
"""
import hashlib
import json
from pathlib import Path
import sympy as sp

GATES = 0
x,t,p,y,u,s,tau = sp.symbols('x t p y u s tau')
D = p**3-y**2
source = {p:t*(1+x*x*t),y:x*t*t*(1+x*x*t),u:x*x*t}
log_chart = {p:s*s+tau,y:s*(s*s+tau),u:s*s/tau,x:s/tau}


def check(ok,label):
    global GATES
    if not ok:
        raise RuntimeError(label)
    GATES += 1


def emit(tag,*values):
    print(tag,json.dumps(values,separators=(',',':')),flush=True)


def equal(a,b,label):
    check(sp.expand(a-b) == 0,label)


def src(f):
    return sp.expand(sp.sympify(f).subs(source,simultaneous=True))


def bracket(f,g):
    return sp.expand(sp.diff(f,x)*sp.diff(g,t)-sp.diff(f,t)*sp.diff(g,x))


def core_images(R):
    rp,ry = sp.diff(R,p),sp.diff(R,y)
    AR = (2*D+3*p**3)*R+p*D*rp
    BR = p*(-2*y*R+D*ry)
    CR = 4*y*R+2*p*y*rp+3*p**3*ry
    ER = 10*R+2*p*rp+3*y*ry
    return {
        p:sp.expand(-D*BR),y:sp.expand(D*AR),
        D:sp.expand(-D**2*CR),u:sp.expand(p**3*y*ER),
        x:sp.expand(p*(5*p**3+4*y*y)*R+p*p*(p**3+y*y)*rp
                    +p*y*(2*p**3+y*y)*ry)},(AR,BR,CR,ER)


def core_delta(f,images):
    return sp.expand(sp.diff(f,p)*images[p]+sp.diff(f,y)*images[y])


def d_divisible(f,n,label):
    if n == 0:
        check(sp.Poly(f,p,y).is_multivariate or f == 0 or f.is_number,label)
        return
    remainder = sp.rem(sp.Poly(f,y,p),sp.Poly(D**n,y,p))
    check(remainder.is_zero,label)


def valuation(f,variable):
    f = sp.expand(f)
    if f == 0:
        return None
    return min(int(term.as_powers_dict().get(variable,0)) for term in sp.Add.make_args(f))


controls = [sp.Integer(1),p,y,p*p,p*y,y*y,D,1+p+y+p*y]
for R in controls:
    images,parts = core_images(R)
    S = p*p*D*R
    S_source = src(S)
    for generator in (x,u,p,y,D):
        equal(bracket(src(generator),S_source),src(images[generator]),
              'literal source generator '+str(generator)+' R='+str(R))
    equal(core_delta(D,images),images[D],'derived cusp sign')
    AR,BR,CR,ER = parts
    uniform_x = D*(p*(5+9*u)*R+p*p*(1+2*u)*sp.diff(R,p)
                   +p*y*(2+3*u)*sp.diff(R,y))
    uniform_u = D*(1+u)*y*ER
    equal(src(uniform_x),src(images[x]),'uniform D factor x')
    equal(src(uniform_u),src(images[u]),'uniform D factor u')
    check(valuation(sp.diff(S_source,t),t) >= 4,'source x motion starts at four')
    check(valuation(sp.diff(S_source,x),t) >= 6,'source t motion starts at six')
    equal(sp.diff(sp.diff(S_source,t),x)-sp.diff(sp.diff(S_source,x),t),0,
          'source volume divergence')
    emit('GENERATOR',str(R),'all five source images and uniform D factors PASS')

# A declared full finite depth bank, with literal lowering witnesses.
for R in (sp.Integer(1),p+y):
    images,parts = core_images(R)
    dp_over_D = sp.cancel(images[p]/D)
    dy_over_D = sp.cancel(images[y]/D)
    for a in range(4):
        for b in range(4-a):
            for c in range(3):
                for e in range(3-c):
                    f = x**a*u**b*p**c*y**e
                    response = 0
                    if a:
                        response += a*x**(a-1)*u**b*p**c*y**e*images[x]
                    if b:
                        response += b*x**a*u**(b-1)*p**c*y**e*images[u]
                    lowered_D_monomial = (x**(a-1)*u**b*y*p if a else
                                          x**a*u**(b-1)*y*y if b else D)
                    if c:
                        response += c*lowered_D_monomial*p**(c-1)*y**e*dp_over_D
                    if e:
                        response += e*lowered_D_monomial*p**c*y**(e-1)*dy_over_D
                    actual = bracket(src(f),src(p*p*D*R))
                    equal(src(response),actual,'literal lowered depth witness')
                    response = sp.expand(response)
                    bound = max(0,a+b-1)
                    check(all(sum(monom[:2]) <= bound
                              for monom,_ in sp.Poly(response,x,u,p,y).terms()),
                          'witness depth bound')
                    log_response = sp.expand(response.subs(log_chart,simultaneous=True))
                    val = valuation(log_response,tau)
                    check(val is None or val >= 1-(a+b),'actual Laurent gain')
emit('DEPTH_BANK','R=(1,p+y), a+b<=3, c+e<=2','120 full actual monomials')

# Iterates start in the core after one derivative of x or u. Check the
# stronger core divisibility, not merely divisibility after source pullback.
for R in (sp.Integer(1),p+y):
    images,_ = core_images(R)
    for seed,offset in ((images[x],0),(images[u],0),(p,0),(D,1)):
        f = seed
        for n in range(5):
            d_divisible(f,n+offset,'exact core iterate divisibility')
            f = core_delta(f,images)
emit('CORE_ITERATES','two carriers, four seeds, five derivative levels PASS')

# Independent log-chart scalar-time exponentials: all computations are in
# K[s,tau]/tau^(N+1), so inverse and volume checks retain the whole jet.
N = 6


def cut(f):
    f = sp.Poly(sp.expand(f),tau)
    return sp.expand(sum(coeff*tau**monom[0] for monom,coeff in f.terms() if monom[0] <= N))


def log_delta(f,S):
    return cut(tau*(sp.diff(f,s)*sp.diff(S,tau)-sp.diff(f,tau)*sp.diff(S,s)))


def exponential(f,S,parameter):
    total,iterate = cut(f),cut(f)
    for n in range(1,N+1):
        iterate = log_delta(iterate,S)
        total = cut(total+parameter**n*iterate/sp.factorial(n))
    return total


def compose(f,ss,tt):
    # Horner in tau avoids large irrelevant high-order intermediate powers.
    polynomial = sp.Poly(f,tau)
    powers = [sp.Integer(1)]
    for _ in range(sp.Poly(f,s,tau).degree(s)):
        powers.append(cut(powers[-1]*ss))
    out = 0
    for n in range(polynomial.degree(),-1,-1):
        coefficient = polynomial.nth(n)
        coefficient_image = 0
        for (power_s,),coefficient_s in sp.Poly(coefficient,s).terms():
            coefficient_image = cut(coefficient_image+coefficient_s*powers[power_s])
        out = cut(out*tt+coefficient_image)
    return out


for R in (sp.Integer(1),p):
    S = sp.expand((p*p*D*R).subs(log_chart,simultaneous=True))
    FF = sp.cancel(S/tau)
    equal(log_delta(s,S),cut(tau*(FF+tau*sp.diff(FF,tau))),'log s sign')
    equal(log_delta(tau,S),cut(-tau*tau*sp.diff(FF,s)),'log tau sign')
    sp1,tp1 = exponential(s,S,1),exponential(tau,S,1)
    sm1,tm1 = exponential(s,S,-1),exponential(tau,S,-1)
    equal(compose(sp1,sm1,tm1),s,'scalar inverse s')
    equal(compose(tp1,sm1,tm1),tau,'scalar inverse tau')
    equal(compose(sp1,sp1,tp1),exponential(s,S,2),'scalar group law s')
    equal(compose(tp1,sp1,tp1),exponential(tau,S,2),'scalar group law tau')
    jac = sp.diff(sp1,s)*sp.diff(tp1,tau)-sp.diff(sp1,tau)*sp.diff(tp1,s)
    equal(cut(tau*jac-tp1),0,'entire log-volume identity')
    equal(sp.expand(sp1).coeff(tau,0),s,'residue coordinate fixed')
    equal(sp.expand(tp1).coeff(tau,1),1,'new tau divided by old tau is a unit')
    equal(compose(cut(S),sp1,tp1),cut(S),'Hamiltonian invariant')
    emit('LOG_FLOW',str(R),'cutoff',N,'inverse, group law, invariant, log volume, unit PASS')

# Sharp depth leading coefficient for the simplest actual carrier.
Slog = sp.expand((p*p*D).subs(log_chart,simultaneous=True))
for d in range(1,7):
    f = s**d*tau**(-d)
    response = sp.expand(tau*(sp.diff(f,s)*sp.diff(Slog,tau)-sp.diff(f,tau)*sp.diff(Slog,s)))
    check(valuation(response,tau) == 1-d,'sharp depth drop')
    equal(sp.expand(response*tau**(d-1)).coeff(tau,0),9*d*s**(d+7),'sharp depth symbol')
for n in range(1,8):
    hostile = s**(2*n)*(s*s+tau)**(2*n)*tau**(-n)
    check(valuation(hostile,tau) == -n,'D-adic versus Laurent topology hostile')
check(sp.rem(sp.Poly(2*y*D,p,y),sp.Poly(p,p,y)).as_expr() != 0,
      'fixed-H carrier S=D fails polynomial core image')
emit('HOSTILES','S=D is not universal','D^N*x^(2N) has log valuation -N',
     'non-LND alone does not decide a single polynomial specialization')
emit('PASS',GATES)
emit('SOURCE_SHA256',hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
