"""Replay the signed mirrors, physical base-three probe and short CRT selector."""
from collections import deque
from math import gcd
from check_seams import flipped, p as f, star_rows, strong, subpaths

D=(4,5,24,19,35,15)

def step(x,eps):
    v=3*x+eps
    k=(v&-v).bit_length()-1
    return v>>k,k

def lift(x,eps):
    y,a=step(x,eps);z,b=step(y,eps)
    return z,(a,b)

def verify_mirror():
    for x in range(1,10000,2):
        for eps in (-1,1):
            a,sa=lift(x,eps);b,sb=lift(x+2,eps)
            cond=sum(sa)==sum(sb) and abs(b-a)==2
            assert cond==(x%16==(9 if eps==1 else 5))
            if cond:
                assert b-a==2
                assert (sa,sb)==(((2,1),(1,2)) if eps==1 else ((1,2),(2,1)))
    for ell in (5,7,11,13,1093,1006003):
        q=ell*ell
        a=-pow(3,-1,ell)%ell
        rep=pow(a,ell,q)
        assert ((3*rep+1)%q==0)==(pow(3,ell-1,q)==1)
    print('legal mirror and Teichmuller zero-preimage fixtures pass')

def ordinal(a,b):
    v=len(a)
    return tuple(r|(((1<<len(b))-1)<<v) for r in a)+tuple(r<<v for r in b)

def physical(v,y,d=D):
    word=bin(v)[2:]
    spine=tuple((1<<len(word))-(2<<i) for i in range(len(word)))
    V=star_rows(spine,int(word[::-1],2),0,0)
    tags=1|(1<<(len(D)+len(spine)))
    return star_rows(ordinal(d,V),tags,0,y)

def verify_physical():
    for v,y,M in [(1,y,150+162*y+25*f(y)) for y in range(12)]+[
          (2**k,0,56*3**k+106*2**k-12) for k in range(1,6)]:
        A,B=physical(v,y),physical(v,y,flipped(D,1,3))
        assert flipped(A,1,3)==B and strong(A) and strong(B)
        assert (subpaths(A)[-1],subpaths(B)[-1])==(M-1,M+1),(v,y,M)
    print('fixed physical arc and base-three observer fixtures pass')

def order(m):
    if m==1: return 1
    d,x=1,2%m
    while x!=1:
        d,x=d+1,2*x%m
    return d

def program(r,m,N):
    q2=1<<r
    g=gcd(25,m)
    d=order(m);t=((r+d-1)//d)*d
    y0=q2*(1+((N-149)*pow(162*q2,-1,g)-1)%g)
    s=N-149-162*y0
    b2=(s*pow(25,-1,q2)-f(y0))%q2
    bm=((s//g)*pow(25//g,-1,m//g)-f(y0))%m
    # Any representative mod m/g works. Choose the block count using bm mod m.
    want=(-(b2-2)//2)%(q2//2)
    blocks=[]
    for j in range(want.bit_length()):
        if want>>j&1: blocks.extend((j+1,2))
    blocks.extend([1]*((bm-len(blocks)-1)%m))
    blocks.append(r+(-sum(blocks)-r)%d)
    y=y0
    for a in blocks:
        y=(((y+1)<<a)-1)<<t
    assert y%(q2*m)==y0%(q2*m)
    assert (149+162*y+25*f(y)-N)%(q2*m)==0
    return y

def verify_program():
    for r,m in [(1,1),(2,3),(4,5),(3,9),(3,25),(1,29),(2,125),
                (4,35),(9,1),(40,1),(80,1)]:
        q=(1<<r)*m
        picks=range(1,q,2) if q<1025 else [1,q-1,q//2+1,345%q]
        for N in picks: program(r,m,N)
    print('short mixed-modulus and high dyadic selectors pass')

if __name__=='__main__':
    verify_mirror()
    verify_physical()
    verify_program()
