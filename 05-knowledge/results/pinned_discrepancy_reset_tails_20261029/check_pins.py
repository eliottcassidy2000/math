"""Numerical proposal search with independent integer Fourier certificates.

SciPy solves only the proposed fractional LP; every reported bound is then
recomputed with rational dual/primal witnesses and outward trigonometric bounds.
This is not an implementation of the Lovett--Meka algorithm.
"""
from fractions import Fraction as F
from math import factorial, gcd
from itertools import product
import json
import numpy as np
from scipy.optimize import linprog

Q=1 << 40
S=1 << 192


def require(ok, msg):
    if not ok:
        raise RuntimeError(msg)


def ceildiv(a,b):
    return -((-a)//b)


def interval_mul(a,b):
    v=[x*y for x in a for y in b]
    return min(v)//S,ceildiv(max(v),S)


def rational_interval(a):
    return a.numerator*S//a.denominator,ceildiv(a.numerator*S,a.denominator)


def sine_centres(n):
    def atan(q):
        value=sum(F((-1)**k,(2*k+1)*q**(2*k+1)) for k in range(20))
        return value,value+F(1,41*q**41)
    a,b=atan(5),atan(239)
    pi_lo,pi_hi=16*a[0]-4*b[1],16*a[1]-4*b[0]
    require(pi_hi < F(22,7), 'pi enclosure')
    tail=ceildiv(((F(22,7)**41)/factorial(41)).numerator*S,
                 ((F(22,7)**41)/factorial(41)).denominator)
    coefficients=[rational_interval(F((-1)**k,factorial(2*k+1))) for k in range(20)]
    centres=[0]*n
    for j in range(1,(n+1)//2):
        xl,xu=2*j*pi_lo/n,2*j*pi_hi/n
        x=(xl.numerator*S//xl.denominator,ceildiv(xu.numerator*S,xu.denominator))
        x2=interval_mul(x,x)
        v=coefficients[-1]
        for c in reversed(coefficients[:-1]):
            v=interval_mul(v,x2)
            v=(v[0]+c[0],v[1]+c[1])
        lo,hi=interval_mul(v,x)
        lo-=tail
        hi+=tail
        centre=(lo+hi)*Q//(2*S)
        require((centre-2)*S <= lo*Q and hi*Q <= (centre+2)*S, 'outward sine error <=2/Q')
        centres[j]=centre
        centres[n-j]=-centre
    return centres


def packets(n,g):
    h=[];v=1
    while v not in h:
        h.append(v);v=v*g % n
    require(v==1 and len(h)%2==1, 'odd subgroup')
    require(all(gcd(v-1,n)==1 for v in h if v!=1), 'free multiplier guard')
    left=set(range(1,n));out=[]
    while left:
        orbit=sorted({min(left)*v % n for v in h})
        require(not set(orbit)&{-v % n for v in orbit}, 'converse guard')
        out.append(orbit)
        left-=set(orbit)|{-v % n for v in orbit}
    return out


def matrix(n,g):
    ps=packets(n,g);s=sine_centres(n)
    return np.array([[2*sum(s[(k*x)%n] for x in orbit) for orbit in ps]
                     for k in range(1,(n+1)//2)],dtype=np.int64)


def sign_interval(C, signs, n):
    peak=int(np.max(np.abs(C@np.array(signs,dtype=np.int64))))
    error=2*(n-1)
    return F(max(0,peak-error),Q),F(peak+error,Q)


def search_signs(A, pins, x0, seed):
    b=A.shape[1]
    free=[i for i in range(b) if i not in pins]
    rng=np.random.default_rng(seed)
    starts=[np.where(x0>=0,1,-1)] + [rng.choice([-1,1],b) for _ in range(5)]
    best=None
    for eps in starts:
        for i,s in pins.items(): eps[i]=s
        v=A@eps
        value=float(np.max(np.abs(v)))
        for _ in range(200):
            if not free: break
            updates=v[:,None]-2*A[:,free]*eps[free][None,:]
            costs=np.max(np.abs(updates),axis=0)
            k=int(np.argmin(costs))
            if costs[k]>=value-1e-10: break
            i=free[k];v=updates[:,k];eps[i]*=-1;value=float(costs[k])
        if best is None or value<best[0]: best=(value,eps.copy())
    return best[1].tolist()


def case(n,g,pin_count,mode):
    C=matrix(n,g);A=C.astype(float)/Q
    M,b=C.shape
    # Fix the largest contributions to the first Fourier row, to stress pin debt.
    order=sorted(range(b),key=lambda i:-abs(int(C[0,i])))
    chosen=order[:pin_count]
    pins={i:(1 if C[0,i]>=0 else -1) for i in chosen}
    if mode=='alternating':
        pins={i:(-1)**k for k,i in enumerate(chosen)}
    free=[i for i in range(b) if i not in pins]
    d=np.array([sum(int(C[k,i])*s for i,s in pins.items()) for k in range(M)],dtype=np.int64)
    B=A[:,free]; df=d.astype(float)/Q
    upper=np.vstack([np.column_stack([B,-np.ones(M)]),np.column_stack([-B,-np.ones(M)])])
    solution=linprog(np.r_[np.zeros(len(free)),1.],A_ub=upper,b_ub=np.r_[-df,df],
                     bounds=[(-1,1)]*len(free)+[(0,None)],method='highs')
    require(solution.success,'fractional LP proposal')
    x=np.zeros(b)
    for i,s in pins.items(): x[i]=s
    x[free]=solution.x[:-1]
    # Rational primal point: feasible in the box; its true norm is an upper bound.
    denom=10**6
    nums=[max(-denom,min(denom,int(round(t*denom)))) for t in x]
    for i,s in pins.items(): nums[i]=denom*s
    maxnum=max(abs(sum(int(c)*v for c,v in zip(row,nums))) for row in C)
    primal=F(maxnum+2*(n-1)*denom,Q*denom)
    # Rational dual point with l1<=1, giving a valid lower bound independently of LP status.
    y=-solution.ineqlin.marginals[:M]+solution.ineqlin.marginals[M:]
    yn=[int(round(float(v)*10**9)) for v in y]
    yd=max(10**9,sum(abs(v) for v in yn))
    obj=sum(a*int(v) for a,v in zip(yn,d))
    obj-=sum(abs(sum(yn[k]*int(C[k,i]) for k in range(M))) for i in free)
    dual=max(F(0),F(obj-2*(n-1)*yd,Q*yd))
    require(dual<=primal,'rational primal/dual sandwich')
    eps=search_signs(A,pins,x,n+pin_count)
    require(all(eps[i]==s for i,s in pins.items()),'retain all prescribed signs')
    low,high=sign_interval(C,eps,n)
    require(high>=dual,'integral completion exceeds dual lower bound')
    optimum=None
    if b<=7:
        lows=[]; highs=[]
        for tail in product((-1,1),repeat=len(free)):
            signs=[pins.get(i,0) for i in range(b)]
            for i,s in zip(free,tail): signs[i]=s
            lo,hi=sign_interval(C,signs,n)
            lows.append(lo);highs.append(hi)
        optimum=[str(min(lows)),str(min(highs))]
        require(dual<=min(highs),'dual survives every integral completion')
    return {'n':n,'g':g,'packets':b,'free_packets':len(free),'pin_mode':mode,
            'pins':[[i,pins[i]] for i in sorted(pins)],'fractional_norm_interval':[str(dual),str(primal)],
            'completion_norm_interval':[str(low),str(high)],'integral_optimum_interval':optimum,
            'completion_upper_over_sqrt_n_diagnostic':float(high)/(n**.5),
            'pin_witness_numerators':nums,'dual_witness_numerators':yn,'dual_denominator':yd,
            'sign_word':''.join('+' if s==1 else '-' for s in eps)}


if __name__=='__main__':
    cases=[(43,36,0,'aligned'),(43,36,2,'aligned'),(43,36,6,'aligned'),
           (133,11,8,'alternating'),(1807,1764,30,'aligned'),(1807,1764,280,'aligned')]
    print(json.dumps({'coefficient_denominator':Q,'sine_error_bound':str(F(2,Q)),
                      'cases':[case(*v) for v in cases]},indent=2,sort_keys=True))
