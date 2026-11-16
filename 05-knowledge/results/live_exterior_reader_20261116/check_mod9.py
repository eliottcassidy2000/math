#!/usr/bin/env python3
"""Exact finite DP and mod-9 two-sided phase-cell decoder."""
from functools import lru_cache
from fractions import Fraction as F
from random import Random
from check_exteriors import H, rev, model, first

@lru_cache(None)
def table(d):
    pairs={(0,0):0}
    for j in range(14):
        m=1<<j
        pairs={((a+m*e)%d,(2*r+e)%9):b+m*e
           for (a,r),b in pairs.items() for e in (0,1)}
    assert len(pairs)==9*d
    return pairs

def bracket9(N,x0,y0,r):
    n=(N//88).bit_length()-1
    p,q=(n+1)//2,n//2
    k=n//10
    A,B,K,M,R,T,S,err,addr=model(n,k,x0,y0)
    W=1<<16
    w=(N-K+M+1)/(S*S)+(1<<48)
    z=w/(R*T)
    alpha0=85 if z<=165 else F(179,2)
    ell=min((l for l in range(1,84) if l%3),key=lambda l: abs(l-alpha0**2/z))
    c=R//T
    g=c*ell
    sa,sb=g*w,w/g
    u=first(-1,R-1,lambda j:(A+9*j)**2>=sa)
    v=first(0,T,lambda j:(B-j)**2<sb)-1
    cap=min((R-1-u)//g,(T-1-v)//9)
    a=(g&-g).bit_length()-1
    d=g>>a
    block=1<<(14+a)
    fixed=u%(1<<a)
    width=p-k
    base_rev=rev(x0,k)*R+(y0+S*v)*(1<<p)
    target_rev=(-11+82*((2<<n)-base_rev)-r)*pow(82,-1,9)%9
    def good(t):
        # Use one whole aligned block contained in the u-future.
        assert t>=0
        base=((u+g*t+block-1)//block)*block
        assert width>=14+a and base+block<=R
        revbase=rev(base>> (14+a),width-14-a)
        r0=((target_rev-revbase)*pow(pow(2,width-14-a,9),-1,9)-rev(fixed,a)*(1<<14))%9
        u0=((u-base-fixed)*pow(1<<a,-1,d))%d if d>1 else 0
        uu=base+fixed+(table(d)[u0,r0]<<a)
        assert (uu-u)%g==0 and rev(uu,width)%9==target_rev
        tt=(uu-u)//g
        assert t<=tt<=t+W
        return tt
    def val(t):
        uu,vv=u+g*t,v+9*t
        assert 0<=t<=cap and 0<=uu<R and 0<=vv<T
        hh=H(n,addr(uu,vv))
        assert hh==S*S*(A+9*uu)*(B-vv)+K+err(uu,vv)
        assert hh%9==r
        return hh
    lo=good(0)
    assert val(lo)>N
    high_probe=max(1,lo)
    while True:
        hi=good(high_probe)
        if val(hi)<N:
            break
        high_probe*=2
    def lower(which):
        L,U=lo,hi
        while U-L>2*W:
            m=good((L+U)//2)
            assert L<m<U
            if which(val(m)):
                U=m
            else:
                L=m
        return L,U
    t_up=lower(lambda h:h<=N)[0]
    t_dn=lower(lambda h:h<N)[1]
    return n,[(addr(u+g*t,v+9*t),val(t)) for t in [t_dn,t_up]]

if __name__=='__main__':
    for d in range(1,84,2):
        if d%3:
            assert len(table(d))==9*d
    rng=Random(4224)
    for n in [110,111,120,121,180,181]:
        D=1<<n
        S=1<<(n//10)
        for N in [88*D,165*D,176*D-1,rng.randrange(88*D,176*D)]:
            for x0,y0 in [(0,0),(S-1,S-1),(rng.randrange(S),rng.randrange(S))]:
                for r in range(9):
                    nn,[(bb,lo),(cc,hi)]=bracket9(N,x0,y0,r)
                    assert nn==n and lo<N<hi and lo%9==hi%9==r
                    q=n//2
                    for b in [bb,cc]:
                        assert (b>>q)%S==x0 and rev(b%(1<<q),q)%S==y0
    print('The exact mod-9 joint phase-cell DP and two-sided decoder passed.')
