#!/usr/bin/env python3
"""All odd moduli, with the literal repair and immutable address sum."""
from math import gcd,isqrt
from transfer import steer,append,count,suffix_length

K=(135,137,139)

def choose_sum(N,m,g):
    s=isqrt(N//108)
    s+=108*s*s<N
    s+=(g-1-s)%g
    if m%5==0:
        j=2*N%3
        s+=g*((N-K[j]-3*s)*pow(3*g,-1,5)%5)
    return s

def solve_residue(N,m,S):
    """Literal polynomial/residue solution; no impossibility fallback."""
    j=2*N%3
    d=gcd(45,m)
    C=lambda x:432*x*(S-x)+129*x+258*(S-x)+15+K[j]
    a=next(r for r in range(d) if (N-C(r))%d==0)
    reduced=m//d
    w=((N-C(a))//d*pow(45//d,-1,reduced))%reduced if reduced>1 else 0
    return a,w,j

def program(N,m,x,y,short=None):
    Kbits=suffix_length(m) if short is None else short-1
    assert (x+y)%(1<<(Kbits+1))==(1<<(Kbits+1))-1
    a,b=append(x>>(Kbits+1),y>>(Kbits+1),"1")
    z,w,j=solve_residue(N,m,x+y)
    if short is None:
        tail=steer(a,b,m,z,w)
    else:
        from short_guard import reachable
        W=(w-p_sum(a,b))%m
        tail=reachable(m,a,b,1,Kbits,True)[z,W,1]
    X,Y=append(a,b,tail)
    assert X+Y==x+y and (count(X,Y,j)-N)%m==0
    return X,Y,j

def p_sum(a,b):
    from window import p
    return p(a)+p(b)

def across(lo,hi,P):
    assert not P(lo) and P(hi)
    while hi-lo>1:
        mid=(lo+hi)//2
        if P(mid):hi=mid
        else:lo=mid
    return lo,hi

def bracket(N,m,short=None):
    s=(N.bit_length()-1-8)//6
    Kbits=suffix_length(m)+1 if short is None else short
    assert s>Kbits,"large-target regime for the algebraic guard"
    g,q=1<<s,1<<Kbits
    S=choose_sum(N,m,g)
    a=((S+1)//2+q-1)//q
    cap=(S-q)//q
    def upper(t):
        x=t*q;y=S-x
        return 432*x*y+129*x+258*y+45*S+45+K[-1]
    _,b=across(a,cap,lambda t:upper(t)<=N)
    def V(t):
        x,y,j=program(N,m,t*q,S-t*q,short)
        return x,y,j,[count(x,y,i) for i in range(3)]
    i,_=across(a,b,lambda t:V(t)[-1][0]<=N)
    _,j=across(a,b,lambda t:V(t)[-1][-1]<=N)
    lo,hi=V(j),V(i)
    assert lo[0]+lo[1]==hi[0]+hi[1]==S
    assert lo[-1][-1]<=N<hi[-1][0]
    return lo,hi

def tests():
    from random import Random
    checks=0
    for m in range(3,120,2):
        for n in range(30):
            for s in range(25):
                if m%5==0 and (3*s+K[2*n%3]-n)%5:continue
                a,w,j=solve_residue(n,m,s)
                assert (432*a*(s-a)+129*a+258*(s-a)+15+K[j]+45*w-n)%m==0
                checks+=1
    print(checks,"arithmetic residue controls")
    rng=Random(1511)
    for m in (3,5,9,15,25,27,45):
        k=suffix_length(m)+4
        g=1<<k
        n=rng.randrange(m)
        S=7*g-1
        if m%5==0:
            S+=g*((n-K[2*n%3]-3*S)*pow(3*g,-1,5)%5)
        x=(S+1)//2;y=S-x
        a,b,j=program(n,m,x,y)
        assert (count(a,b,j)-n)%m==0
    print("7 literal all-odd transfers, including powers of 3 and 5")
    from short_guard import certify
    for m in (3,5,9):
        short=certify(m)+1
        N=(1<<128)+131
        a,b=bracket(N,m,short)
        assert a[-1][a[2]]%m==N%m==b[-1][b[2]]%m
        print("same-sum exact residue bracket, modulus",m,"bits",N.bit_length())
if __name__=='__main__':tests()
