"""Certified prime-safe phase, cycle-height inverse, and rail constructors."""
from math import gcd,isqrt
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'moonshine_prime_blocks'))
from replay import R,pad,wrap,reads,profile,dp
from clone_replay import P,TT,vectors
from frame_replay import join

C11=(16,17,3,7,12)
D5=(4,1,18,7,11); I5=24
RAIL={1:((8,9,3,52,7,23),48),7:((20,9,10,97,14,23,55,127),96)}

def toggle(a,j):
    if (a,j)==(11,0):return C11,11,100
    if a==3:
        k=j-1;D,I,*_=R[1,1];B=-12+16*2**k+28*3**k
    else:
        assert a==5
        k=j-2;D,I=D5,I5;B=-88+112*2**k+192*3**k
    assert k>=1
    return wrap(D,I,TT(k)),1+a*2**j,B

def phase(p,L):
    """Return a physical source with exact H,B, H congruent p mod 15*2^L."""
    assert L>=1 and gcd(p,30)==1
    T=();h=b=1
    def add(a,j):
        nonlocal T,h,b
        U,s,t=toggle(a,j);T=join(T,U);h*=s;b*=t
    for j in range(1,L):
        if ((h^p)>>j)&1:
            if j==1:add(11,0)
            else:add(5 if j%4==3 else 3,j)
    a=L+(1-L)%4
    if a<2:a+=4
    for _ in range(3):
        if h%5==p%5:break
        add(3,a)
    if h%3!=p%3:
        a=max(3,L)
        add(5,a+(1-a)%2)
    assert (h-p)%(15<<L)==0
    return (T,h,b) if T else ((0,),1,2)

def paired_phase(p,L,scale):
    T,h,b=phase(p,L)
    D,I=RAIL[scale]
    T=pad(T,L)
    return wrap(D,I,T),scale*(60*(b<<L)+h)

def two_valued_rail(p,scale):
    """Inverse on the valuation-separated branch of the long-cycle rail."""
    if p<=1 or p%2==0:return None
    m=((p-1)&-(p-1)).bit_length()-1
    h,b=1+2**m,4*2**m+2*3**m-2
    z,r=divmod(p-h,60*b)
    if r or z<=0 or z&(z-1):return None
    k=z.bit_length()-1
    alpha=60*b
    if k+(alpha&-alpha).bit_length()-1<=m:return None
    D,I=RAIL[scale]
    return wrap(D,I,pad(P(m+2),k)),scale*p

def run():
    for scale,(D,I) in RAIL.items():
        d,g,w=reads(D,I);assert (g,d+w-g)==(60*scale,scale)
    assert profile(C11)==(11,100)
    for a in (3,5):
        D,I=R[1,1][:2] if a==3 else (D5,I5)
        Q=wrap(D,I,TT(1)); t=(len(D),)
        assert dict(vectors(Q,t,1,2))==({(1,):13,(2,):6} if a==3 else {(1,):41,(2,):20})
        assert dict(vectors(Q,t,2,3))==({(1,):104,(2,):100,(3,):28} if a==3 else {(1,):712,(2,):688,(3,):192})
    for a,start in ((3,2),(5,3)):
        for j in range(start,8):
            T,h,b=toggle(a,j);assert profile(T)==(h,b)
    for m in range(1,9):
        assert profile(P(m+2))==(1+2**m,4*2**m+2*3**m-2)
    for L in range(1,10):
        for p in range(1,15<<L,2):
            if gcd(p,30)!=1:continue
            T,h,b=phase(p,L)
            assert (60*(b<<L)+h-p)%(15<<L)==0
    print('Every unit phase through L=9 reproduced')
    for p in (13457,53777,27525137,217252097):
        assert all(p%i for i in range(2,isqrt(p)+1))
        for s in (1,7):
            T,ans=two_valued_rail(p,s)
            assert ans==s*p
            if p<60000:assert sum(dp(T)[-1])==ans
        print('cycle inverse',p,7*p)

if __name__=='__main__':run()
