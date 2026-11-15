#!/usr/bin/env python3
from random import Random
from verify_matching_retention import graph,pi,plug

def mirrored(a,mark,u,v):
    l=len(a); n=2*l+5; apex=n
    m=[0]*(n+1)
    rev=list(range(n-1,-1,-1))
    for i in range(n):
        for j in range(i+1,n):
            if j<l:right=(a[i]>>j)&1
            elif i>=l+5:right=(a[rev[j]]>>rev[i])&1
            else:right=True
            if right:m[i]|=1<<j
            else:m[j]|=1<<i
    center=[u,1-v,0,v,u]
    paint=mark+center+[1-x for x in mark[::-1]]
    for i,b in enumerate(paint):
        if b:m[apex]|=1<<i
        else:m[i]|=1<<apex
    return m

def main():
    rng=Random(592055)
    for trial in range(32):
        l=rng.randrange(1,4)
        a=graph(l,rng.randrange(1<<(l*(l-1)//2)))
        mark=[rng.randrange(2) for _ in range(l)]
        v=rng.randrange(2)
        s,t=mirrored(a,mark,0,v),mirrored(a,mark,1,v)
        for k in (1,2,3,4):
            assert pi(t,k)-pi(s,k)==2*k**4*pi(a,k)**2
        assert pi(plug([2,4,1],t),1)-pi(plug([2,4,1],s),1)==2*pi(a,1)**2+32*pi(a,2)**2
    print('128 exact mirror profile differences and 32 literal strong-reader steps passed')

if __name__=='__main__':main()
