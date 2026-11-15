#!/usr/bin/env python3
"""Complete exact inverse for the specified physical two-address bank."""
from math import isqrt
Q,H,A,B,F,E=432,45,129,258,15,30
K=(135,137,139)

def p(x):
    s=0
    while x:
        if not x&1:
            m=x//2;s+=m+(m&1)
        x//=2
    return s

def value(x,y,j):
    return Q*x*y+A*x+B*y+F*(x&1)+E*(y&1)+H*(p(x)+p(y))+K[j]

def scalar(a,d,x):
    return a*x+H*p(x)+d*(x&1)

def solve_row(a,d,R):
    # Return all nonnegative m with a*m+H*p(m)+d*(m mod 2)=R.
    if R<0:return
    stop=2*isqrt(H*R)+2*H+2*d+4
    depth=0
    while a<<depth<stop:depth+=1
    st=[(a,d,R,0,0)]
    while st:
        u,v,c,b,k=st.pop()
        if c<0:continue
        lo=max(0,(c-v+u+H-1)//(u+H))
        hi=c//u
        if hi<lo:continue
        if lo==hi or k==depth:
            assert lo==hi
            if scalar(u,v,lo)==c:yield (lo<<k)|b
        else:
            st.append((2*u+H,H,c,b,k+1))
            st.append((2*u,0,c-u-v,b|(1<<k),k+1))

def decode(N):
    if N%2==0 or N<min(K):return []
    ans=set()
    L=isqrt(N//Q)
    for x in range(1,L+1):
        tail=A*x+F*(x&1)+H*p(x)
        for j,k in enumerate(K):
            if (N-k)%3:continue
            for y in solve_row(Q*x+B,E,N-tail-k):
                if y>=1:ans.add((x,y,j))
    for y in range(1,L+1):
        tail=B*y+E*(y&1)+H*p(y)
        for j,k in enumerate(K):
            if (N-k)%3:continue
            for x in solve_row(Q*y+A,F,N-tail-k):
                if x>=1:ans.add((x,y,j))
    return sorted(ans)

if __name__=='__main__':
    import sys
    for t in sys.argv[1:]:
        N=int(t);print(N,decode(N))
