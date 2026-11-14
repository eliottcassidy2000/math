"""Fixed nonnegative boundary profiles, exact inversion of the corridor length."""
from fractions import Fraction as F
from random import Random

ONE=(1,0,0)
RHO=(0,1,0)

def mul(a,b):
    c=[0]*5
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            c[i+j]+=x*y
    for i in (4,3):
        for j in range(i-3,i):
            c[j]+=c[i]
    return tuple(c[:3])

def power(k):
    a=ONE
    b=RHO
    while k:
        if k&1: a=mul(a,b)
        k//=2
        if k: b=mul(b,b)
    return a

def sign(a):
    if not any(a): return 0
    l,r=F(183,100),F(46,25)
    while True:
        lo=sum(c*(l if c>=0 else r)**i for i,c in enumerate(a))
        hi=sum(c*(r if c>=0 else l)**i for i,c in enumerate(a))
        if lo>0: return 1
        if hi<0: return -1
        m=(l+r)/2
        if m*m*m-m*m-m-1<0: l=m
        else: r=m

def eval_count(q,u,v):
    a,b,c=(1,1,1),(1,2,2),(3,3,4)
    for j in range(q):
        a,b,c=b,c,tuple(x+y+z for x,y,z in zip(a,b,c))
    e,f,g=a
    return u[0]*(e*v[0]+f*v[1])+u[1]*(f*v[0]+g*v[1])

def inverse(n,u,v):
    assert n>0 and all(x>=0 for x in u+v) and any(u) and any(v)
    # rho^2 N(1/rho) has entries rho^2+1, rho^2+rho, rho^2+rho+1.
    a=u[0]*v[0]+u[1]*v[1]
    b=u[0]*v[1]+u[1]*v[0]+u[1]*v[1]
    c=sum(u)*sum(v)
    p2=mul((-1,-2,3),(-1,-2,3))
    b2=mul((a,b,c),(a,b,c))
    def hits_upper(q):
        rhs=mul(b2,power(2*q+1))
        return sign(tuple(x-n*n*y for x,y in zip(rhs,p2)))>0
    l=r=4
    while not hits_upper(r):
        r*=2
    while l<r:
        m=(l+r)//2
        if hits_upper(m): r=m
        else: l=m+1
    return tuple(q for q in (*range(4),l) if eval_count(q,u,v)==n)

if __name__=='__main__':
    rng=Random(13)
    for _ in range(440):
        u=(rng.randrange(10**20),rng.randrange(10**20))
        v=(rng.randrange(10**20),rng.randrange(10**20))
        q=rng.randrange(1,130)
        n=eval_count(q,u,v)
        assert q in inverse(n,u,v)
        assert not inverse(n+1,u,v)
    assert inverse(1,(1,0),(1,0))==(0,1)
    assert inverse(25281,(1,2),(1,2))==(13,)
    print('exact cubic-field inverse: 440 generated targets and adjacent nonhits, plus the reset')
