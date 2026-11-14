#!/usr/bin/env python3
from math import isqrt
def g(q):
    z=1
    while q:
        q,b=divmod(q,2)
        if not b:z+=q+(q&1)
    return z
def height(u,s,q):
    M=1<<s
    return M*M*u*(u+1)+M*M//4-u+2*g(u)+u%2-3+2*g(q)-q*q-2*q+q%2-2*(q>>(s-2))
def inverse(N):
    for s in range(4,(N+2).bit_length()):
        u=isqrt(N)>>s
        lo,hi=0,(1<<(s-1))-1
        while lo<=hi:
            q=(lo+hi)//2
            h=height(u,s,q)
            if h==N:
                yield u,s,q
                break
            if h>N:lo=q+1
            else:hi=q-1
def run():
    for s in range(4,14):
        eq=[]
        for q in range((1<<(s-1))-1):
            k=((q+1)&-(q+1)).bit_length()-1
            v=q>>(k+1)
            gap=6*v+2+2*(v%2)+2*(((q+1)>>(s-2))-(q>>(s-2)))
            assert height(3,s,q)-height(3,s,q+1)==gap
            if gap==2:eq.append(q)
        assert eq==[(1<<j)-1 for j in range(s-2)]
    for r in range(4,15):
        A=1<<r;bad=longest=0;hits=[]
        for N in range(A*A+3,A*A+2*A,2):
            if list(inverse(N)):
                hits.append(N);bad=0
            else:
                bad+=1;longest=max(bad,longest)
        assert len(hits)<=2*(r-1)
        assert longest*(2*r-1)>=A-2*r+1
        print(r,len(hits),longest)
    print('PASS: all exact small-slice gaps and pinned dry bands.')
if __name__=='__main__':run()
