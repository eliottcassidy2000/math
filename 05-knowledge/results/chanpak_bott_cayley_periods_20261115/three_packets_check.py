"""Exact compiler and replay for nonlinear literal tournament guards."""
from functools import lru_cache
from random import Random


def p(x):
    w = 0
    while x:
        x, b = divmod(x,2)
        if not b:
            w += x + (x&1)
    return w


def F(x):
    w = 0
    while x:
        x, b = divmod(x,2)
        if b:
            w += (x+(x&1))//2
    return w


def reverse(x,n):
    y=0
    for _ in range(n):
        x,b=divmod(x,2)
        y=2*y+b
    return y


def length(k):
    return 10*k-7 if k>=3 else k+2


@lru_cache(None)
def packet_table(k):
    n=k-2
    q=1<<n
    if n<1:
        raise ValueError("need k>=3")
    fs=[F(x) for x in range(q)]
    rs=[reverse(x,n) for x in range(q)]
    witness={}
    for x in range(q):
        for y in range(q):
            f=fs[x]+fs[y]+rs[x]*y
            word=(y<<(n+1))|x
            if f!=F(word):
                raise AssertionError("bilinear identity")
            witness.setdefault(f%q,word)
    pair={}
    for a,aw in witness.items():
        for b,bw in witness.items():
            pair.setdefault((a+b)%q,(aw,bw))
    triple={}
    for target in range(q):
        for c,cw in witness.items():
            if (target-c)%q in pair:
                triple[target]=(*pair[(target-c)%q],cw)
                break
    if len(triple)!=q:
        raise AssertionError("missing clock residue")
    return triple


def guard(k,correction,low):
    if low not in (0,2):
        raise ValueError("bad terminal state")
    if k==2:
        return low
    word=0
    for w in packet_table(k)[correction%(1<<(k-2))]:
        word=(word<<(3*k-3))|(w<<k)
    return word<<2|low


def H(x,y,j):
    return (432*x*y+129*x+258*y+15*(x&1)+30*(y&1)
            +45*(p(x)+p(y))+135+2*j)


def steer(k,u,v,N):
    m=1<<k
    if (u+v+1)%m or min(u,v)<1 or not N&1:
        raise ValueError("need positive prefixes summing to -1 modulo 2^k and odd target")
    state=p(u)+p(v)-u+(u&1)
    j=2*N%3
    low=(N-state-3-2*j)%4
    X,Y=2*low+1,-2*low-1
    constant=432*X*Y+129*X+258*Y+45+45*low+135+2*j
    want=(N-constant)*pow(45,-1,m)%m
    delta=(state-want)%m
    if delta%4:
        raise AssertionError("mod-4 state")
    word=guard(k,delta//4,low)
    q=1<<length(k)
    if word>=q:
        raise AssertionError("long clock")
    a,b=2*(u*q+word)+1,2*(v*q+q-1-word)+1
    val=tuple(H(a,b,i) for i in range(3))
    if val[1]-val[0]!=2 or val[2]-val[1]!=2 or (val[j]-N)%(3*m):
        raise AssertionError("physical counts")
    if a+b!=2*q*(u+v+1):
        raise AssertionError("address sum")
    return a,b,val


if __name__=="__main__":
    rng=Random(20261115)
    nchecks=0
    for k in range(2,10):
        m=1<<k
        for _ in range(80):
            u,v=(rng.randrange(1,1<<80) for _ in (0,1))
            v+=-(u+v+1)%m
            N=2*rng.randrange(1<<80)+1
            steer(k,u,v,N)
            nchecks+=1
    print(f"{nchecks} exact prefix-to-clock checks; complete packet coverage for k=3..9.")
