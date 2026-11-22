#!/usr/bin/env python3
"""Independent absolute test of complementary suffixes and the fixed-z inverse."""
from math import isqrt
from random import Random

def F(v):
    a=1
    while v:
        p=v>>1
        if v%2==0:
            a+=p+p%2
        v=p
    return a

def J(r):
    z=0
    while r:
        p=r>>1
        z+=(1 if r%2==0 else -1)*(p+p%2)
        r=p
    return z

def rev(r,k):
    return int(f"{r:0{k}b}"[::-1],2)

def T(q):
    return 64*q*q+142*q+81

def direct(G,x,y,r,dest):
    u=G*(16*x +(9 if dest else 5))+r
    w=G*(16*y +(0 if dest else 4))+G-1-r
    return 81+142*(u+w)+256*u*w+5*(F(u)+F(w))

def decode(N,G,z,dest):
    c=16*z+8*dest
    xmin=max(2,z+1)
    xmin+=xmin%2
    S0=G*(32*xmin-16*z+10)-1
    if S0>N:
        return None
    lo,hi=0,(N-S0)//(64*G)+1
    A=256*(abs(c)+1)*G*G
    while lo<hi:
        mid=(lo+hi)//2
        S=S0+64*G*mid
        if T(S)-64*c*c*G*G+A+5*(S+2)<N:
            lo=mid+1
        else:
            hi=mid
    S=S0+64*G*lo
    if S>N or T(S)-64*c*c*G*G-A>N:
        return None
    x=xmin+2*lo
    y=x-z
    a=16*x+(9 if dest else 5)
    b=16*y+(0 if dest else 4)
    C=T(S)+5*(F(a)+F(b)+(G-1)*a+1)
    diff=C-N
    if diff<=0:
        return None
    q=isqrt(diff//64)
    if diff>=16*(2*q+1)**2:
        q+=1
    if c<0:
        q=-q
    r2=q-c*G-1
    if r2%2:
        return None
    r=r2//2
    if r<=0 or r>=G or r%2==0 or direct(G,x,y,r,dest)!=N:
        return None
    return x,y,r

def require(v,msg):
    if not v:
        raise ArithmeticError(msg)

def test():
    for k in range(2,14):
        G=1<<k
        for r in range(1,G,2):
            s=G-1-r
            jj=J(r)
            require(abs(jj)<=r-1,"size of compensated weight")
            require(F(r)+F(s)-rev(r,k)+(r<G//2)==3+jj-2*(r>=G//2),
                    "complement parity")
    rng=Random(44003)
    for i in range(4000):
        k=rng.randrange(2,100)
        G=1<<k
        r=2*rng.randrange(G//2)+1
        x=2*rng.randrange(1,10**20)
        y=(x if i%4==0 else rng.randrange(1,10**20))
        z=x-y
        for dest in (False,True):
            N=direct(G,x,y,r,dest)
            require(decode(N,G,z,dest)==(x,y,r),"lost exact point")
            require(decode(N+2,G,z,dest) is None,"accepted a nonpoint neighbor")
            require(decode(N-2,G,z,dest) is None,"accepted a nonpoint neighbor")
    print("verified all odd complementary suffixes through 13 bits and 8,000 absolute source/destination cases, with 16,000 nearby rejection controls")

if __name__=="__main__":
    test()
