#!/usr/bin/env python3
"""Exact certificates for a target-centred phase-cell walk. Python 3."""
from fractions import Fraction as F
from random import Random

def rev(x, n):
    return int(f'{x:0{n}b}'[::-1], 2) if n else 0

def J(w):
    ans = suffix = prev = 0
    for e in map(int, bin(w)[2:]):
        if e:
            ans += suffix + prev
        suffix = 2*suffix + e
        prev = e
    return ans

def H(n, b):
    return -11 + 82*((2 << n)-rev(b,n)) + 9*(2*b-J(b))

def model(n,k,x0,y0):
    p,q=(n+1)//2,n//2
    P,Q,S=1<<p,1<<q,1<<k
    R,T=P//S,Q//S
    X,Y=rev(x0,k),rev(y0,k)
    a,b=82*P+9*x0,2*Q-y0
    A=F(a,S)+F(9*Y,S*S)
    B=F(b,S)-F(X,S*S)
    K=(-11-82*R*X-9*J(x0)+18*T*Y-9*J(y0)
       -9*(x0%2)*(y0%2)+F(a*X-9*Y*b,S)+F(9*Y*X,S*S))
    assert K==-11+9*(F(x0*X+y0*Y,S)-J(x0)-J(y0)+F(X*Y,S*S)-(x0%2)*(y0%2))
    M=91*R+27*T+18
    def E(u,v):
        return (-82*rev(u,p-k)-9*J(u)+18*rev(v,q-k)-9*J(v)
                -9*((u%2)*(x0>=S//2)+(v%2)*(y0>=S//2)))
    def addr(u,v):
        return ((x0+S*u)<<q)+rev(y0+S*v,q)
    return A,B,K,M,R,T,S,E,addr

def first(a,b,yes):
    """Bisection maintaining a false lower and true upper boundary."""
    assert a<b and not yes(a) and yes(b)
    while b-a>1:
        m=(a+b)//2
        if yes(m):
            b=m
        else:
            a=m
    return b

def bracket(N,x0,y0):
    n=(N//88).bit_length()-1
    k=n//10
    a,b,K,M,R,T,S,E,addr=model(n,k,x0,y0)
    w=(N-K+M+1)/(S*S)
    z=w/(R*T)
    alpha0=85 if z<=165 else 90
    ell=(F(alpha0**2,1)/z+F(1,2)).__floor__()
    c=R//T
    sa,sb=c*ell*w,w/(c*ell)
    u=first(-1,R-1,lambda j:(a+9*j)**2>=sa)
    v=first(0,T,lambda j:(b-j)**2<sb)-1
    dy=c*ell
    def val(t):
        U,V=u+dy*t,v+9*t
        assert 0<=U<R and 0<=V<T
        hh=H(n,addr(U,V))
        assert hh==S*S*(a+9*U)*(b-V)+K+E(U,V)
        assert abs(E(U,V))<=M
        return hh
    cap=min((R-1-u)//dy,(T-1-v)//9)
    assert val(0)>N and cap>0
    hi=1
    while val(hi)>=N:
        if hi==cap:
            raise ValueError('too small for the stated asymptotic walk')
        hi=min(cap,2*hi)
    down=first(0,hi,lambda t:val(t)<N)
    up=first(0,hi,lambda t:val(t)<=N)-1
    return n,[(addr(u+dy*t,v+9*t),val(t)) for t in [down,up]]

def check():
    rng=Random(1849)
    for n in range(4,15):
        p,q=(n+1)//2,n//2
        for k in range(1,q+1):
            S=1<<k
            for _ in range(60):
                x0,y0=rng.randrange(S),rng.randrange(S)
                a,b,K,M,R,T,S,E,addr=model(n,k,x0,y0)
                u,v=rng.randrange(R),rng.randrange(T)
                assert H(n,addr(u,v))==S*S*(a+9*u)*(b-v)+K+E(u,v)
                assert abs(E(u,v))<=M
    worst=0
    for n in range(42,105):
        D=1<<n
        S=1<<(n//10)
        for N in [88*D, 176*D-1,165*D-1,165*D+1]+[rng.randrange(88*D,176*D) for _ in range(6)]:
            for x0,y0 in [(0,0),(S-1,S-1),(S-1,0),(0,S-1)]+[(rng.randrange(S),rng.randrange(S)) for _ in range(4)]:
                nn,[(l,lo),(r,hi)]=bracket(N,x0,y0)
                assert nn==n and lo<N<hi
                q=n//2
                for b in (l,r):
                    assert (b>>q)%S==x0 and rev(b%(1<<q),q)%S==y0
                worst=max(worst,max(N-lo,hi-N)/D**.4)
    print('Phase identity and strict two-sided bracketing verified; maximum tested error/D^0.4:',worst)

if __name__=='__main__':
    check()
