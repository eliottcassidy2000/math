"""Exact antidiagonal inversion of the modulus-dependent cube-root placement."""
from argparse import ArgumentParser
from functools import lru_cache
from itertools import product
from math import isqrt
from pathlib import Path
from random import Random
import sys

# This script sits beside or in the sibling of the inherited packet.
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'chanpak_decimal_kernel_20261122'))
from hyperbola_controls import chart
from fermat_receiver import TARGET, address_height

PRIMES={4:(5,17,29),8:(5,29,257)}

def envelope(data,N,S):
    h,q,c,t,e,r,k=data
    a=q*S+t-c
    D=a*a-4*q*(N-k-c*S)
    H=D+4*q*(h*S+r+e)
    if H<0:return ()
    j=isqrt(H)
    i=isqrt(D-1)+1 if D>0 else 0
    den=2*q
    intervals=[]
    for L,R in ((-j,-max(i,1)),(i,j)):
        lo=max(0,(a+L+den-1)//den);hi=min(S,(a+R)//den)
        if lo<=hi:intervals.append((lo,hi))
    return intervals

def prepare(n,data,N):
    h,q,c,t,e,r,k=data
    assert h==PRIMES[n][0]*PRIMES[n][1]*PRIMES[n][2]
    @lru_cache(None)
    def roots(p,s,px,py):
        return tuple(x for x in range(p)
                     if (k+c*s+(q*s+t-c)*x-q*x*x+r*px+e*py-N)%p==0)
    fac=PRIMES[n]
    @lru_cache(None)
    def crt(Smod,parity):
        rr=[roots(p,Smod%p,parity,(Smod-parity)&1) for p in fac]
        return [((parity*h + sum((2*h//p)*pow(2*h//p,-1,p)*x
                                  for p,x in zip(fac,rrr)))%(2*h))
                for rrr in product(*rr)]
    return crt

def inverse(n,data,N,S,crt):
    h=data[0]
    for lo,hi in envelope(data,N,S):
        if hi-lo<=20:
            _,q,c,t,e,r,k=data
            for x in range(lo,hi+1):
                y=S-x
                if (N-k-t*x-c*y-q*x*y-r*(x&1)-e*(y&1))%(2*h)==0:
                    yield x,y
        else:
            for a in (0,1):
                for u in crt(S%(2*h),a):
                    yield from ((x,S-x) for x in range(lo+(u-lo)%(2*h),hi+1,2*h))

def scan(n,last=None,first=0):
    N=TARGET; data=chart(n)
    h,q,c,t,e,r,k=data;z=isqrt(N//q)
    cap=(N-k)//min(c,t)
    crt=prepare(n,data,N)
    presented=enclosed=survivors=0
    closest=None;hits=[];exponents=[]
    for b in range(first,cap.bit_length() if last is None else min(cap.bit_length(),last+1)):
        B=1<<b
        # Complete superset: replace gcd(h,M) in the placement upper
        # limit by h, and drop all remaining congruences on S.
        lo=2*z+2*B;hi=min(cap,2*z+(h+2)*B-1)
        S0=lo+(-1-lo)%B
        if S0>hi: continue
        exponents.append(b)
        for S in range(S0,hi+1,B):
            presented+=1
            if not envelope(data,N,S):continue
            enclosed+=1
            for x,y in inverse(n,data,N,S,crt):
                survivors+=1
                diff=address_height(data,x,y)-N
                if closest is None or abs(diff)<abs(closest[0]):closest=(diff,x,y,b,S)
                if not diff:hits.append((x,y,b))
    print("spine",n,"z",z,"exponents",exponents[0],exponents[-1],
          "sums",presented,"nonempty envelope",enclosed,
          "phase candidates",survivors,"nearest",closest,"hits",hits,flush=True)
    return presented,enclosed,survivors,closest,hits

def terminal(n,b):
    N=TARGET;data=chart(n);h,q,c,t,e,r,k=data
    cap=(N-k)//min(c,t)
    crt=prepare(n,data,N)
    B=1<<b; seen=passed=0;closest=None;hits=[]
    for S in range(B-1,cap+1,B):
        for x,y in inverse(n,data,N,S,crt):
            passed+=1
            diff=address_height(data,x,y)-N
            if closest is None or abs(diff)<abs(closest[0]):closest=(diff,x,y,S)
            if not diff:hits.append((x,y))
        seen+=1
    print("all aspect ratios: spine",n,"terminal bits",b,"sums",seen,
          "phase candidates",passed,"nearest",closest,"hits",hits,flush=True)

def controls():
    R=Random(61)
    for n in (4,8):
        data=chart(n)
        for _ in range(70):
            S=R.randrange(0,250)
            t=R.randrange(S+1)
            N=address_height(data,t,S-t)+2*R.randrange(-3,4)
            got={(x,y) for x,y in inverse(n,data,N,S,prepare(n,data,N))
                 if address_height(data,x,y)==N}
            want={(x,S-x) for x in range(S+1) if address_height(data,x,S-x)==N}
            assert got==want
    print("literal antidiagonal controls pass",flush=True)

if __name__=='__main__':
    p=ArgumentParser()
    p.add_argument('--spine',type=int,choices=(4,8))
    p.add_argument('--first',type=int,default=0)
    p.add_argument('--last',type=int)
    p.add_argument('--terminal',type=int)
    a=p.parse_args()
    if a.spine and a.terminal is not None:terminal(a.spine,a.terminal)
    elif a.spine:scan(a.spine,a.last,a.first)
    else:controls()
