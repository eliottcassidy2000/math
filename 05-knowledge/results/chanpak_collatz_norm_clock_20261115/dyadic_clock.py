"""Both-1 tournament clock, no third-party dependencies."""
from math import isqrt

def check(x):
    w=0
    while x:
        x,b=divmod(x,2)
        if not b: w+=x+(x&1)
    return w

def read(x,y,j):
    return (432*x*y+129*x+258*y+15*(x&1)+30*(y&1)
            +45*(check(x)+check(y))+135+2*j)

def excursion_guard(k):
    return k+sum(h+1+k for h in range(2,k))+2

def carry(x,y,N,k):
    # x is the controlled pre-address. y gets the complemented word.
    g=1<<k
    if (x+y+1)%g: raise ValueError("bad rail")
    z=check(x)+check(y)-x+(x&1)
    j=2*N%3
    b=(N-z-3-2*j)%4
    assert b in (0,2)
    x0=b
    # final x0 modulo 2^k, for the two final odd addresses.
    # At reset K may have changed by multiple of 4; zeros do not change K.
    # For final t=10 from reset, K adds S+1 = 0 mod 2^k.
    X,Y=2*x0+1,-2*x0-1   # actual final y ≡ −x0−1 before both-1
    constant=(432*X*Y+129*X+258*Y+45+45*x0+135+2*j)
    want=(N-constant)*pow(45,-1,g)%g
    delta=(z-want)%g
    assert delta%4==0, (x,y,j,x0,z,want,constant)
    word=0;length=k
    for h in range(2,k):
        b=h+1+k; word <<= b; length+=b
        if (delta>>h)&1:
            word+=(1<<h | 1)<<k
    word=(word<<2)+x0;length+=2
    B=1<<length
    assert length==excursion_guard(k)
    u,v=B*x+word,B*y+(B-1-word)
    assert u%g==x0
    H=tuple(read(2*u+1,2*v+1,i) for i in range(3))
    assert (H[j]-N)%(3*g)==0
    return u,v,H

def crossing(lo,hi,f):
    assert not f(lo) and f(hi),(lo,hi)
    while hi-lo>1:
        b=(lo+hi)//2
        if f(b):hi=b
        else:lo=b
    return lo,hi

def bracket(N,k):
    B=1<<excursion_guard(k)
    g=1<<max(B.bit_length()-1, (N.bit_length()-9)//6)
    R0=isqrt(N//108); R0+=108*R0**2<N
    S0=(R0+1)//2-1
    S=S0+(g-1-S0)%g
    lo=((S+1)//2+B-1)//B
    hi=(S-B)//B
    def W(i): return carry(i,S//B-i,N,k)
    b0,_=crossing(lo,hi,lambda i:W(i)[2][0]<=N)
    _,b1=crossing(lo,hi,lambda i:W(i)[2][-1]<=N)
    low,high=W(b1),W(b0)
    assert low[0]+low[1]==high[0]+high[1]==S
    assert low[2][-1]<=N<high[2][0]
    return low,high

if __name__=='__main__':
    from random import Random
    rng=Random(20261115)
    for k in range(2,16):
        for _ in range(100):
            u,v=rng.randrange(1,1<<40),rng.randrange(1,1<<40)
            v+=(-(u+v+1))%(1<<k)
            carry(u,v,2*rng.randrange(1<<80)+1,k)
    for k in range(2,8):
        for _ in range(8):
            n=(1<<max(160,24*excursion_guard(k)))+2*rng.randrange(1<<130)+1
            low,high=bracket(n,k)
            radius=max(n-low[2][0],high[2][-1]-n)
            assert radius**8 < 10000**8*n**3
    print("1400 carry controls, 48 literal address brackets")
