#!/usr/bin/env python3
"""Direct small tournament path counts and old-chart exclusions."""
from functools import lru_cache
from math import isqrt

def digits(n):
    return bin(n)[2:]


def q(n):
    s=0
    while n:
        n,b=divmod(n,2)
        s+=b*((n+1)//2)
    return s


def old_has(t,h,a,b,c,k):
    for v in range(max(0,(t-k)//(2*b-h)+1)):
        d=t-k-2*b*v+2*h*q(v)
        if d<0:
            continue
        L=2*(a*v+c)
        for u in range(-(-d//L),d//(L-h)+1):
            if L*u-2*h*q(u)==d:
                return True
    return False


def graph(word,guard,u,endpoint):
    core_masks=((177,181),(11524,11588))
    out=(19,7)[guard]
    word=digits(word)
    tail=''.join(str(1-int(z)) for z in digits(2*u)[::-1])
    p=len(word)
    n=p+8+len(tail)
    t=n-1
    A=[0]*n
    def add(i,j):
        A[i]|=1<<j
    for i in range(t):
        for j in range(i+1,t):
            if j==p and i<p and word[i]=='1':
                add(j,i)
            else:
                add(i,j)
    m=0
    for i in range(p+1,p+7):
        for j in range(i+1,p+7):
            if core_masks[guard][endpoint]>>m&1:
                A[i]^=1<<j
                add(j,i)
            m+=1
    for i in range(t):
        outer_beats=(i<p or
          (p<i<p+7 and ((out>>(i-p-1))&1)) or
          (i>=p+7 and tail[i-p-7]=='1'))
        add(*( (t,i) if outer_beats else (i,t) ))
    return A


def hp(A):
    @lru_cache(None)
    def go(i,todo):
        if not todo:
            return 1
        out=A[i]&todo
        s=0
        while out:
            bit=out&-out
            out-=bit
            s+=go(bit.bit_length()-1,todo^bit)
        return s
    n=len(A)
    allbits=(1<<n)-1
    return sum(go(i,allbits^(1<<i)) for i in range(n))


if __name__=='__main__':
    for w,c,u,t in [(30,0,0,31151),(20,1,1,176459)]:
        A,B=(graph(w,c,u,e) for e in (0,1))
        assert sum((x^y).bit_count() for x,y in zip(A,B))==2
        assert [hp(A),hp(B)]==[t,t+2]
        for x in [t,t+2]:
            assert all(x%d for d in range(2,isqrt(x)+1))
        assert not any(old_has(x,*params) for x in (t-2,t,t+2)
            for params in ((29,640,191,158,111),(651,17616,6768,2691,2313)))
        print(len(A),t,t+2)
