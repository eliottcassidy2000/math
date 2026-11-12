#!/usr/bin/env python3
"""Exact rational Stein certificates and signed Gaussian brackets; no floating integration."""
import sympy as s
from functools import lru_cache
x,y=s.symbols('x y'); g=s.Rational(1,1000)

def df(n):
    return s.factorial2(n) if n>0 else s.Integer(1)
def expectation(p):
    return s.expand(sum(c*df(i-1)*df(j-1) for (i,j),c in s.Poly(p,x,y).terms() if i%2==0 and j%2==0))
def run(lam, limit):
    W=(x*x+y*y)/2+g*x**(2*lam); R=(lam-1)*g*x**(2*lam)
    def D(h): return s.expand((x*s.diff(h,x)+y*s.diff(h,y))/2+h)
    @lru_cache(None)
    def P(m,h):
        if m==0:return h
        if m==1:return s.expand(D(h)-R*h)
        return s.expand(P(m-1,s.expand(D(h)+(m-1)*h-R*h))+(m-1)*P(m-2,s.expand(R*h)))
    @lru_cache(None)
    def cert(m,h):
        if m==0:return (s.Integer(0),s.Integer(0))
        a=cert(m-1,s.expand(D(h)+(m-1)*h-R*h))
        b=cert(m-2,s.expand(R*h)) if m>=2 else (0,0)
        return tuple(s.expand(-W**(m-1)*h*v+a[k]+(m-1)*b[k]) for k,v in enumerate((x/2,y/2)))
    for m in range(limit+1):
        h=x*x+y*y+1
        p,q=cert(m,h)
        boundary=s.diff(p,x)+s.diff(q,y)-s.diff(W,x)*p-s.diff(W,y)*q
        assert s.expand(W**m*h-P(m,h)-boundary)==0
        if lam==2:
            for h,a in ((s.Integer(1),1),(x*x,2),(x*y,2)):
                closed=sum((-1)**j*s.binomial(m,j)*s.rf(a+j,m-j)*(g*x**4)**j*h for j in range(m+1))
                assert s.expand(P(m,h)-closed)==0
    print('lambda',lam,'exact polynomial certificates m=0..',limit,'PASS')
    return W,P
W,P=run(2,6); run(3,4)
print('Quartic radial corrections P_m(1), z=g*x^4:')
z=s.symbols('z')
for m in range(5):print(m,s.expand(P(m,s.Integer(1))).subs(x**4,z/g).expand())
# For P>=0, Taylor's pointwise alternating bounds are integrated exactly.
def bracket(poly,N):
    a=[(-g)**j*expectation(s.expand(poly*x**(4*j)))/s.factorial(j) for j in range(N+2)]
    u=sum(a[:N+1]);v=u+a[N+1]
    return min(u,v),max(u,v)
print('Normalized Gaussian tilted moments, g=1/1000; rigorous rational intervals:')
for m in (0,1,2,4,6):
    lo,hi=bracket(s.expand(W**m),14)
    print('m',m,'lower',lo,'upper',hi,'width',hi-lo)
    print(' display',s.N(lo,18),s.N(hi,18))
# Pointwise no fixed finite perturbation order is all-order convergence.
print('Gaussian partition Taylor term ratio |a_(j+1)/a_j| = g*(4j+1)*(4j+3)/(j+1); tends to infinity')
print('ALL CHECKS PASS')
def signed_bracket(poly,N):
    lo=hi=0
    for (i,j),c in s.Poly(poly,x,y).terms():
        assert i%2==j%2==0
        a,b=bracket(x**i*y**j,N)
        lo+=c*(a if c>=0 else b); hi+=c*(b if c>=0 else a)
    return lo,hi
print('Exact-representative transport: direct versus Euler-reduced interval width at Taylor order 14')
for m in (1,2,4,6):
    a,b=bracket(s.expand(W**m),14)
    c,d=signed_bracket(P(m,s.Integer(1)),14)
    assert a<=c<=d<=b
    print('m',m,'direct_width',s.N(b-a,8),'reduced_width',s.N(d-c,8),'gain',(s.N((b-a)/(d-c),8)), 'intersection_width',s.N(min(b,d)-max(a,c),8))
