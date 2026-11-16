"""An exact-integer implementation of the paired-exterior approximation.
For very small inputs the proof's interior phase margin may be insufficient.
"""
from fractions import Fraction
from itertools import combinations
from math import isqrt

Z, g, I = 82, 9, -11

def rev(x,n):
    return int(f'{x:0{n}b}'[::-1],2)

def count(b,n):
    D=1<<n
    sites=[i for i in range(n) if b>>i&1]
    J=sum((1<<(j-i-1))+(j==i+1) for i,j in combinations(sites,2))
    return I+Z*(2*D-rev(b,n))+g*(2*b-J)

def floor(f): return f.numerator//f.denominator
def ceil(f): return -floor(-f)

def first(lo,hi,pred):
    """Binary search inclusive interval with pred monotone False -> True."""
    if lo>hi or not pred(hi): raise ValueError('insufficient phase margin')
    while lo<hi:
        m=(lo+hi)//2
        if pred(m): hi=m
        else: lo=m+1
    return lo

def live(N, M=1):
    if N<88 or N%2==0 or M<1: raise ValueError('need odd N≥89 and M≥1')
    n=(N//88).bit_length()-1
    m=n//2; k=n-m
    B=1<<m; C=1<<k
    D=B*C
    L=1<<max(0,((B//M).bit_length()-1)//5)
    z=Fraction(N-I,D)
    f=lambda j: Fraction(2*B,L)-Fraction(N-I,L*(C*Z+g*L*j))
    deriv=lambda j: Fraction(g*(N-I),(C*Z+g*L*j)**2)
    # Find an interior point close to the opposite-corners diagonal x=1-u.
    w=4096
    j=first(82,4014,lambda j:(Z+g*Fraction(j,w))*(1+Fraction(j,w))>=z)
    u=Fraction(j,w); radius=Fraction(1,128)
    lo=ceil((u-radius)*C/L); hi=floor((u+radius)*C/L)
    je0=ceil((u-radius/2)*C/L); je1=floor((u+radius/2)*C/L)
    # Use a fixed rational grid for a bounded denominator. Reducing the
    # result improves constants, but is not needed by the proof.
    lower,upper=deriv(je1),deriv(je0)
    for denom in range(1,1<<20):
        tau=Fraction(ceil(denom*lower),denom)
        if tau<=upper: break
    else: raise ValueError('insufficient slope margin')
    r,s=tau.numerator,tau.denominator
    j0=first(je0,je1,lambda j:deriv(j)<=tau)
    hmax=(hi-j0)//s
    F=lambda h:f(j0+s*h)-r*h
    q0=floor(F(0))-1
    out=[]
    for i in range(M):
        q=q0-i
        h=first(0,hmax,lambda h:F(h)<=q)
        U=L*(j0+s*h); X=L*(q+r*h)
        if not (0<=U<C and 0<=X<B):raise ValueError('insufficient exterior margin')
        b=B*U+rev(X,m)
        out.append((n,b,count(b,n),U,X))
    return out

if __name__=='__main__':
    for n in (64,65,96,97,128,129,256,257):
        for t in (88,113,164,175):
            N=(t<<n)+1
            ans=live(N,3)
            assert len({r[1] for r in ans})==3
            B=1<<(n//2)
            # Integer bound for a readable loose constant.
            assert max(abs(r[2]-N)**5 for r in ans) < 10**40*B**4*3
    print('Paired-exterior decoder checked at both parities and 64–257 bits.')
