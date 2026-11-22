#!/usr/bin/env python3
"""Exact small-denominator rational tangent corridors at E61, spine four."""
from fractions import Fraction as F
from math import gcd,isqrt
from functools import lru_cache
from argparse import ArgumentParser

H,q,c,t,e,r,k=2465,2293760,1344270,712944,1740,986,423385
N=(1<<60)*((1<<61)-1)-1
def p(x):
    z=0
    while x:
        x,b=divmod(x,2)
        if not b:z+=x+(x&1)
    return z
def height(x,y):
    return k+t*x+c*y+q*x*y+r*(x&1)+e*(y&1)+H*(p(x)+p(y))
def ceil(x):
    return -((-x.numerator)//x.denominator)

def ndata(b0,target):
    Q=10*q; C=c+q*b0; T=10*t; K=k+t*b0+r*(b0&1)
    M=target-K
    return Q,C,T,K,M,M*Q+T*C
def endpoints(nd,a,b,A,B,X):
    Q,C,T,K,M,delta=nd
    def y(t):return F(M-T*t,Q*t+C)
    YA=y(A);YB=y(B)
    if YA<0:return 1,0
    YY=YA.numerator//YA.denominator
    G=F(H*(10*B+1+YY)+e,Q*A+C)
    U=max(b*YA+a*A,b*YB+a*B)
    slope=a-F(b*delta,(Q*X+C)**2)
    L=b*y(X)+a*X+min(slope*(A-X),slope*(B-X))-b*G
    return ceil(L),U.numerator//U.denominator

@lru_cache(None)
def roots(q,A,B,C):
    return tuple(i for i in range(q) if (A*i*i+B*i+C)%q==0)
@lru_cache(None)
def inva(a,b):
    return pow(a,-1,b) if b>1 else 0

def scans(nd,b0,a,b,A,B,X,target,stat):
    Q,C,T,K,M,delta=nd
    low,high=endpoints(nd,a,b,A,B,X)
    stat[0]+=1;stat[1]+=B-A+1
    if low>high:stat[2]+=1
    for z in range(low,high+1):
        stat[3]+=1
        first=A+(z*inva(a,b)-A)%b
        if first>B:continue
        y=(z-a*first)//b
        if y<0:continue
        length=min((B-first)//b,y//a)
        for parity in (0,1):
            if length<parity:continue
            x0=first+b*parity; y0=y-a*parity
            A2=-4*Q*a*b
            B2=2*(Q*(b*y0-a*x0)+T*b-C*a)
            C2=K+Q*x0*y0+T*x0+C*y0+e*(y0&1)-target
            R17=roots(17,A2%17,B2%17,C2%17)
            R29=roots(29,A2%29,B2%29,C2%29)
            top=(length-parity)//2
            for i in R17:
                for j in R29:
                    residue=(i+17*((j-i)*12%29))%493 # 17^{-1} mod 29=12
                    for u in range(residue,top+1,493):
                        x=x0+2*b*u;y=y0-2*a*u
                        stat[4]+=1
                        if height(10*x+b0,y)==target:
                            yield 10*x+b0,y

def directions(nd,L,R,M):
    Q,C,T,K,goal,delta=nd
    l=F(delta,(Q*R+C)**2)
    h=F(delta,(Q*L+C)**2)
    for b in range(1,M+1):
        for a in range(ceil(b*l),(b*h).numerator//(b*h).denominator+1):
            if not a or gcd(a,b)!=1:continue
            X=(isqrt(b*delta//a)-C)//Q
            if L<=X<=R:
                w=isqrt((Q*X+C)**3//(128*b*delta*Q))
                w=min(w,X//8)
                if w:yield max(L,X-w),min(R,X+w),a,b,X

def run(target,lo,hi,M,b0):
    nd=ndata(b0,target)
    stats=[0]*5
    covered=lo-1
    hits=[]
    for A,B,a,b,X in sorted(directions(nd,lo,hi,M)):
        A=max(A,covered+1)
        if A>B:continue
        # A previous interval might clip the original anchor; tangent
        # remains a valid supporting line even when the anchor is outside.
        hits += list(scans(nd,b0,a,b,A,B,X,target,stats))
        covered=B
    return stats,hits

def control():
    # Compare the complete exact strip (not sampled points) with ordinary height.
    x,y=2871,509
    n=height(x,y)
    b0=1
    for a,b in ((1,1),(2,3),(7,5),(23,4)):
        nd=ndata(b0,n);A,B,X=277,297,287
        actual=set(scans(nd,b0,a,b,A,B,X,n,[0]*5))
        expected={(10*z+b0,v)
                  for z in range(A,B+1) for v in range(1+n//(q*(10*z+b0)))
                  if height(10*z+b0,v)==n}
        assert actual==expected and (x,y) in actual
if __name__=='__main__':
    pp=ArgumentParser(description=__doc__)
    pp.add_argument('--den',type=int,default=128)
    pp.add_argument('--lo',type=int,default=115_000_000_000_000)
    pp.add_argument('--hi',type=int,default=200_000_000_000_000)
    o=pp.parse_args()
    control()
    for b0 in (0,1):
        stats,hits=run(N,o.lo,o.hi,o.den,b0)
        print('left residue',b0,'corridors/width/line-free/lines/exact',*stats,'hits',hits)
