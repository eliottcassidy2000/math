"""End-to-end positive-address placement check using the three-packet solver."""
from math import isqrt
from random import Random
from three_packets_check import length, steer

def across(lo,hi,f):
    if f(lo) or not f(hi):
        raise RuntimeError("bracket endpoints not opposite")
    while hi-lo>1:
        mid=(lo+hi)//2
        if f(mid): hi=mid
        else: lo=mid
    return lo,hi

def bracket(N,k):
    B=length(k)
    q=1<<B
    g=1<<max(B+k,(N.bit_length()-1-4*B)//6)
    R=isqrt(N//108)
    R+=108*R*R<N
    S0=(R-1)//2
    S=S0+(g-1-S0)%g
    lo=((S+1)//2+q-1)//q
    hi=S//q-1
    def W(a):
        return steer(k,a,S//q-a,N)
    i,_=across(lo,hi,lambda a:W(a)[-1][0]<=N)
    _,j=across(lo,hi,lambda a:W(a)[-1][-1]<N)
    left,right=W(j),W(i)
    if not(left[2][-1]<N<right[2][0] and sum(left[:2])==sum(right[:2])):
        raise RuntimeError("not literal whole windows")
    return left,right

if __name__=='__main__':
    rng=Random(3411)
    total=0
    for k in range(3,9):
        for _ in range(10):
            h=160*k
            N=(1<<h)+2*rng.getrandbits(h-2)+1
            lo,hi=bracket(N,k)
            eps=max(N-lo[2][0],hi[2][-1]-N)
            if eps**8>1000**8*N**3:
                raise RuntimeError("unexpected placement radius at k="+str(k))
            total+=1
    print(f"{total} two-sided brackets at the growing-precision scale")
