"""Exact forward and inverse tests on the mirror-padding route."""
from fractions import Fraction as F
from itertools import product
import argparse
from prime_support_replay import L,R
from padding_replay import slots, rooted
from frame_replay import join, add, cut
from replay import pad, profile

def inverse(L,R):
    x,z=F(L[1],L[2]),F(L[3],L[2])
    y,w=F(R[1],R[2]),F(R[3],R[2])
    d=0
    while x>y*F(2)**(-d) and z<w*F(3,2)**d:d-=1
    while x<y*F(2)**(-d) and z>w*F(3,2)**d:d+=1
    # Crude exact endpoints, broaden to avoid losing endpoint coincidences.
    lo=hi=d
    while not (x<y*F(2)**(-lo) and z>w*F(3,2)**lo):lo-=1
    while not (x>y*F(2)**(-hi) and z<w*F(3,2)**hi):hi+=1
    for d in range(lo,hi+1):
        A=x-y*F(2)**(-d);B=z-w*F(3,2)**d
        if A*B<=0:continue
        C=A/B;a=max(0,-d)
        while F(3)**a<=C:
            k=C/3**a;m=1
            while F(3,2)**(m-1)<=k:
                if F(2)**(m-1)>=k:
                    for wbits in range(1,1<<((m+1)//2)):
                        I=sum(1<<i for i in range(m) if wbits>>min(i,m-1-i)&1)
                        t=sum(2**i for i in range(m) if I>>i&1)
                        u=sum(3**i*2**(m-1-i) for i in range(m) if I>>i&1)
                        if F(u,t)==k:yield a,a+d,m,I
                m+=1
            a+=1

def ternary_sums(vals):
    out=[0]
    for v in vals:out=[x+z for z in (0,v,-v) for x in out]
    return out

def test_m(m):
    n=m-1
    vals=[2**i*(3**(n-i)-1)-2**(n-i)*(3**i-1) for i in range(m//2)]
    z=len(vals)//2
    a=ternary_sums(vals[:z]);b=ternary_sums(vals[z:])
    return a.count(0)==b.count(0)==1 and set(a).intersection(b)=={0}

def replay(limit):
    a,b=([v[-1] for v in slots(P,3)] for P in (L,R))
    out=list(inverse(a,b))
    assert out==[(0,0,7,20),(1,1,5,10),(2,2,3,5)]
    frames=[]; inner=[]
    for x,y,m,I in out:
        h=pad(L,x);t=pad(R,y)
        U=join(h,join(pad((),m),t));frames.append(add(U,I<<len(h)))
        t=join(pad((),y),R)
        U=join(h,join(pad((),m),t));inner.append(add(U,I<<len(h)))
    assert all(t==inner[0] for t in inner)
    assert len({tuple(sorted(r.bit_count() for r in t)) for t in frames})==3
    print('all decompositions have order 24, identical masks with inner-side padding',out)
    print('outer-side padding gives three distinct degree multisets')
    for m in range(1,8):
        U=pad((),m)
        for I in range(1<<m):
            a=add(U,I);b=cut(a,1<<m)
            l=[int(I>>i&1)-int(I>>(m-1-i)&1) for i in range(m)]
            d=sum(z*2**i for i,z in enumerate(l))
            e=sum(z*2**(m-1-i)*3**i for i,z in enumerate(l))
            h,B=profile(a);hh,BB=profile(b)
            assert BB-B==2*(d+e) and hh-h==d+l[-1]
            E,S=rooted(U,I)[1:3]
            assert profile(a)[1]+2*(S[0]+E[0]-S[1]-E[1])==profile(b)[1]
    assert all(test_m(m) for m in range(1,limit+1))
    print('nonmirror neutral masks excluded through exterior order',limit)

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--limit',type=int,default=48)
    replay(parser.parse_args().limit)
