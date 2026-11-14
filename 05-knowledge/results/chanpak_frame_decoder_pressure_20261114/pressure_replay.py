#!/usr/bin/env python3
"""Exact rational certificates for the first-triangle terrace.

No third-party dependencies. The replay is finite arithmetic, not a survey
of large tournaments. Numbers printed as decimal are presentation only.
"""
from fractions import Fraction as F

DEN=10**8
CUT=16
BINS=80

def power_bounds(x,s):
    """Return consecutive multiples of 1/DEN bounding x**(-s), x>=1."""
    x=F(x)
    assert x>=1 and 0<s<1
    n,d=s.numerator,s.denominator
    a=x.numerator**n
    b=x.denominator**n*DEN**d
    lo,hi=0,DEN
    while hi-lo>1:
        k=(hi+lo)//2
        if k**d*a>=b:
            hi=k
        else:
            lo=k
    return F(lo,DEN),F(hi,DEN)

def floor_grid(x): return x.numerator*DEN//x.denominator
def ceil_grid(x): return -floor_grid(-x)

def kernels(k):
    A,B,D=[1,1,3],[1,2,3],[1,2,4]
    for q in range(3,k+1):
        for row in (A,B,D):
            row.append(sum(row[-3:]))
    return tuple(zip(A,B,D,strict=True))

def step(K,r,p):
    A,B,D=K
    y=1<<r
    a=A+(B-A)*p
    b=y*(B+(D-B)*p)
    return a+b,b/(a+b)

def test_upper(s):
    ks=kernels(CUT)
    hi=lambda x: power_bounds(x,s)[1]
    # Missing q and r: assume f<=1 and G >= (1+2^r) A_q.
    R=sum(hi(1+(1<<r)) for r in range(1,CUT+1))
    Q=sum(hi(ks[q][0]) for q in range(1,CUT+1))
    assert all(5*ks[q][0]>=9*ks[q-1][0] for q in (4,5,6))
    t=hi(2)
    R_tail=hi(1<<(CUT+1))/(1-t)
    t=hi(F(9,5))
    Q_tail=hi(ks[CUT][0])*t/(1-t)
    tail_grid=ceil_grid(R_tail*(Q+Q_tail)+R*Q_tail)
    f=lambda p: 1-F(137,1000)*(p-F(2,3))
    worst=F(0)
    for i in range(BINS):
        l=F(2,3)+F(i,3*BINS)
        u=l+F(1,3*BINS)
        grid=tail_grid
        for q in range(1,CUT+1):
            for r in range(1,CUT+1):
                G_l,p_l=step(ks[q],r,l)
                G_u,p_u=step(ks[q],r,u)
                assert G_l<=G_u
                grid+=ceil_grid(hi(G_l)*f(min(p_l,p_u)))
        worst=max(worst,F(grid,DEN)/f(u))
    assert worst<F(999,1000),worst
    return worst,F(tail_grid,DEN)

def test_lower_word(s):
    ks=kernels(CUT)
    lo=lambda x: power_bounds(x,s)[0]
    f=lambda p: 1-F(137,1000)*(p-F(2,3))
    worst=F(2)
    for i in range(BINS):
        l=F(2,3)+F(i,3*BINS)
        u=l+F(1,3*BINS)
        grid=0
        for q in range(1,CUT+1):
            for r in range(1,CUT+1):
                G_l,p_l=step(ks[q],r,l)
                G_u,p_u=step(ks[q],r,u)
                assert G_l<=G_u
                grid+=floor_grid(lo(G_u)*f(max(p_l,p_u)))
        worst=min(worst,F(grid,DEN)/f(l))
    assert worst>F(1003,1000),worst
    return worst

def grow(a,b,q,r):
    if q==1:
        return a+2*b,(1<<(r+1))*(a+b)
    if q==2:
        return 3*(a+b),(1<<r)*(3*a+4*b)
    raise ValueError(q)

def admissible(a,b,q):
    return ((a+2*b)%3!=0) if q==1 else (b%3!=0)

def recover_marked(H,P):
    """Return the guarded q,r history at fixed power of two P, or None."""
    if not (P>=4 and P&(P-1)==0 and H>=1):
        return None
    z,u=divmod(H,P)
    a,b=2*z-u,u-z
    rev=[]
    while (a,b)!=(1,2):
        if a<1 or b<2 or a%2!=1 or b%2!=0:
            return None
        v=(b&-b).bit_length()-1
        if a%3:
            q,r=1,v-1
            z=b>>v
            aa,bb=2*z-a,a-z
        else:
            q,r=2,v
            z=b>>v
            u=a//3
            aa,bb=4*u-z,z-3*u
        if not (1<=r<=20 and 1<=aa and 2<=bb and aa+bb<a+b and admissible(aa,bb,q)):
            return None
        assert grow(aa,bb,q,r)==(a,b)
        rev.append((q,r))
        a,b=aa,bb
    return tuple(rev[::-1])

def test_distinct():
    alpha=F(11,15)
    lo=lambda x:power_bounds(x,alpha)[0]
    def c(q,r):return (2+2**(r+1)) if q==1 else (3+4*2**r)
    even=sum(lo(c(1,r)) for r in range(2,21,2))
    odd=sum(lo(c(1,r)) for r in range(1,21,2))
    both=sum(lo(c(2,r)) for r in range(1,21))
    z=(1000,475,832,469)
    M=((both,even,odd,0),(both,0,0,0),
       (both,0,0,even+odd),(0,odd,even,0))
    v=min(sum(c*w for c,w in zip(row,z,strict=True))/w for row,w in zip(M,z,strict=True))
    assert v>F(1006,1000),v
    states=[(1,2,())]
    for _ in range(6):
        next_states=[]
        for a,b,w in states:
            P=1<<(a+2*b).bit_length()
            H=(a+2*b)+P*(a+b)
            assert recover_marked(H,P)==w
            for q in (1,2):
                if admissible(a,b,q):
                    for r in range(1,4):
                        aa,bb=grow(a,b,q,r)
                        assert aa%2==1 and bb%2==0 and (aa%3 or bb%3)
                        next_states.append((aa,bb,w+((q,r),)))
        states=next_states
    return v

if __name__=='__main__':
    upper,tail=test_upper(F(119,125))
    lower=test_lower_word(F(19,20))
    distinct=test_distinct()
    print(f's=.952 upper weighted row bound {float(upper):.9f}; omitted tail < {float(tail):.9f}')
    print(f's=.95 lower weighted row bound {float(lower):.9f}; counts ordered words only')
    print(f'distinct growth: four-state row lower bound {float(distinct):.9f}; exact marker replay passed')
