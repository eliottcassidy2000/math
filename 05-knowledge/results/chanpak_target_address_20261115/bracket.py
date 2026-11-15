#!/usr/bin/env python3
"""Deterministic two exact physical triples bracketing a prescribed target."""
from inverse import value, K

def first_true(test):
    lo,hi=0,1
    while not test(hi):hi*=2
    while hi-lo>1:
        mid=(hi+lo)//2
        if test(mid):hi=mid
        else:lo=mid
    return hi

def bracket(N):
    if N<10000: raise ValueError("use exact decoder for small N")
    g=1
    while 432*(2*g)**10<=N:g*=2
    # Least positive source address on the grid with diagonal above N.
    l=first_true(lambda l: 432*(l*g-1)**2+387*(l*g-1)+K[0]>=N)
    c=l*g-1
    def deduct(t):return 432*(g*t)**2+129*g*t
    B=432*c*c+387*c
    # Above: keep lowest branch's affine part >=N. Treat equality
    # explicitly because first_true assumes its inequality fails at zero.
    t=first_true(lambda t: deduct(t)>B+K[0]-N)-1
    # Below: account for every positive checksum/parity term, and
    # the highest branch. Both x and y have log2(g)-bit suffix of ones.
    upper_correction=45*(2*c//g+1)
    T=first_true(lambda t: deduct(t)>=B+K[-1]+upper_correction-N)
    shots=[]
    for z in (T,t):
        x,y=c+g*z,c-g*z
        if y<=0:raise ValueError("small scale where the asymptotic guard does not apply")
        shots.append((x,y,[value(x,y,j) for j in range(3)]))
    a,b=shots
    assert a[2][-1]<=N<=b[2][0],(N,a,b)
    return a,b

if __name__=='__main__':
    import sys
    for s in sys.argv[1:]:
        N=int(s)
        print(N,bracket(N))
