#!/usr/bin/env python3
from itertools import product
from math import isqrt

def sidecar(w):
    # Direct Hamilton count: one permutation per assignment of chain sides.
    count=0;n=len(w)
    for R in range(1<<n):
        L=((1<<n)-1)^R
        if (not R or w[(R&-R).bit_length()-1]=='O') and (not L or w[L.bit_length()-1]=='I'):
            count+=1
    return count

def word(a,b,c):return 'O'+'I'*a+'O'+'I'*b+'O'+'I'*c
def flip(w,i):return w[:i]+'I'+w[i+1:]
def profile(a,b,c):
    x,y,z=1<<a,1<<b,1<<c
    H=4*x*y*z+(2*y+1)*(z-x)+1-y
    return H,x-2*y*z+y,2*y*x-z+y,4*y*y-1
for a,b in product(range(2,6),range(1,6)):
    w='O'*a+'I'*b; x,y=1<<a,1<<b
    h=sidecar(w);delta=sidecar(flip(w,a-1))-h
    assert h==(x-1)*(y-1)+2 and delta==x//2-y
    d=isqrt((2*delta+1)**2+8*(h-2))
    assert ((2*delta+1+d)//2)+1==x
    assert h%7
hits={}
for a,b,c in product(range(1,4),repeat=3):
    w=word(a,b,c);H,g1,g2,J=profile(a,b,c)
    assert (sidecar(w),sidecar(flip(w,a+1)),sidecar(flip(w,a+b+2)))==(H,H+g1,H+g2)
    if H%7==0:hits[(a%3,b%3,c%3)]=(g1%7,g2%7,J%7)
assert hits=={(2,0,0):(3,1,3),(0,1,0):(6,5,1)}
a,b,c=5,3,6;w=word(a,b,c)
assert (sidecar(w),sidecar(flip(w,a+1)),sidecar(flip(w,a+b+2)))==(66073,65089,66529)
for p in (9439,65089,66529,73):assert all(p%d for d in range(2,isqrt(p)+1))
print("20 two-block inverses; 27 three-run residue classes and exact paths; 18-vertex prime/7p fork")
