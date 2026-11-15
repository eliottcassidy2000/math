#!/usr/bin/env python3
"""Exact dynamic compiler for a fixed-length complementary modular suffix."""
from math import gcd,isqrt
from window import p,value,K
def phase_size(N):
    return 1<<max(1,(N.bit_length()-1-8)//6)

def extensions(x,w,e,S,m):
    yield (2*x%m,(w+x+e)%m,0),"0"
    yield ((2*x+1)%m,(w+S-x+1-e)%m,1),"1"

def reachable(m,a,b,e,n,word=False):
    # Checksum coordinate is relative to the initial checksum.
    states={(a%m,0,e):''}
    S=(a+b)%m
    for _ in range(n):
        after={}
        for (x,w,t),s in states.items():
            for v,bit in extensions(x,w,t,S,m):
                after.setdefault(v,s+bit if word else "")
        states=after; S=(2*S+1)%m
    return states

def certify(m):
    # Universal exact length; only final x-odd states are needed.
    if m<3 or m%2==0:raise ValueError("odd modulus >=3 required")
    for n in range(1,10*m+1):
        if all(sum(e for x,w,e in reachable(m,a,b,t,n))==m*m
               for a in range(m) for b in range(m) for t in (0,1)):
            return n
    raise ValueError("use the algebraic compiler if short search is inconclusive")

def steer(x,y,N,m,n):
    assert gcd(m,30)==1
    q=1<<n
    assert (x+y)%q==q-1
    a,b=x//q,y//q
    assert (a^b)&1, "need complementary preceding pair"
    S=(x+y)%m
    target=x%m
    other=(S-target)%m
    j=2*N%3
    constant=432*target*other+129*target+258*other+15+K[j]
    w=(N-constant)*pow(45,-1,m)%m
    W=(w-p(a)-p(b))%m
    word=reachable(m,a,b,a&1,n,True)[target,W,1]
    k=int(word,2)
    X,Y=a*q+k,b*q+q-1-k
    assert X+Y==x+y and X>0 and Y>0
    V=tuple(value(X,Y,i) for i in range(3))
    assert V[j]%(3*m)==N%(3*m)
    return X,Y,j,V

def across(lo,hi,condition):
    assert not condition(lo) and condition(hi)
    while hi-lo>1:
        u=(hi+lo)//2
        if condition(u):hi=u
        else:lo=u
    return lo,hi

def bracket(N,m,n):
    g=phase_size(N)
    q=1<<n
    assert g>q,"reserve preceding complementary pair"
    a=isqrt(N//108)
    S0=a+(108*a*a<N)
    S=S0+(g-1-S0)%g
    a=((S+1)//2+q-1)//q
    cap=(S-q)//q
    def upper(u):
        x=u*q;y=S-x
        return 432*x*y+129*x+258*y+45*S+45+K[-1]
    _,b=across(a,cap,lambda u:upper(u)<=N)
    def W(u):
        x=u*q
        return steer(x,S-x,N,m,n)
    i,_=across(a,b,lambda u:W(u)[-1][0]<=N)
    _,j=across(a,b,lambda u:W(u)[-1][-1]<=N)
    lo,hi=W(j),W(i)
    assert lo[0]+lo[1]==hi[0]+hi[1]==S
    assert lo[-1][-1]<=N<hi[-1][0]
    return S,lo,hi

if __name__=='__main__':
    import sys
    m=int(sys.argv[1]);n=certify(m)
    print("Universal exact-length guard",m,n)
    for s in sys.argv[2:]:
        N=int(s);print(N,bracket(N,m,n))
