#!/usr/bin/env python3
"""Deterministic two exact physical triples bracketing a prescribed target."""
from inverse import value, K, p

def first_true(test):
    lo,hi=0,1
    while not test(hi):hi*=2
    while hi-lo>1:
        mid=(hi+lo)//2
        if test(mid):hi=mid
        else:lo=mid
    return hi

def low_signature(g,b):
    s=g.bit_length()-1
    if not s:return 0,0
    bits=[b>>i&1 for i in reversed(range(s))]
    return sum((1-z)<<j for j,z in enumerate(bits)),1-bits[0]

def bracket(N):
    if N<10000: raise ValueError("use exact decoder for small N")
    g=1
    while 432*(2*g)**8<=N:g*=2
    c=first_true(lambda c:432*c*c+387*c+K[0]>=N)
    u,b=divmod(c,g)
    A,z=low_signature(g,b)
    fixed_par=45*(c&1) if g>1 else 0
    baseline=432*c*c+387*c+90*(A*u+p(b))+fixed_par
    def deduct(t):return 432*(g*t)**2+129*g*t
    upper_correction=90*(u+z)+(45 if g==1 else 0)
    t=first_true(lambda t: deduct(t)>baseline+K[0]-N)-1
    T=first_true(lambda t: deduct(t)>=baseline+K[-1]+upper_correction-N)
    shots=[]
    for j in (T,t):
        x,y=c+g*j,c-g*j
        if y<=0:raise ValueError("small scale where the asymptotic guard does not apply")
        shots.append((x,y,[value(x,y,i) for i in range(3)]))
    at,bt=shots
    assert at[2][-1]<=N<=bt[2][0],(N,at,bt)
    return at,bt

if __name__=='__main__':
    import sys
    for s in sys.argv[1:]:
        N=int(s)
        print(N,bracket(N))
