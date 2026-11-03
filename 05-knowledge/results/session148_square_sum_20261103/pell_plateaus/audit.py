#!/usr/bin/env python3
"""Independent integer checks for the session-148 Pell plateau claim."""
from math import isqrt

def check(pred, detail):
    if not pred:
        raise RuntimeError(detail)

def square(x):
    return x >= 0 and isqrt(x)**2 == x

def degree(n):
    return isqrt(2*n-1)-isqrt(n)

def degree_direct(n):
    return sum(square(n+j) for j in range(1,n))

def center(n):
    return square(n) and (square(2*n-2) or square(2*n-1))

def above(n):
    return degree(n-1)==degree(n)==degree(n+1)

def threshold(k):
    # Independent block search, rather than the closed formula under test.
    a=0
    while (a+k)**2 > 2*(a+1)**2-3:
        a+=1
    return a,max(a*a,(a+k)**2//2+1)

def closed_threshold(k):
    q=2*(k-1)**2+3
    ceil=isqrt(q)+(isqrt(q)**2 != q)
    a=k-2+ceil
    return a,max(a*a,(a+k)**2//2+1)

for n in range(1,4001):
    check(degree(n)==degree_direct(n), ('direct degree',n))
for k in range(1,1001):
    a,n=threshold(k)
    check((a,n)==closed_threshold(k),('threshold formula',k,a,n))
    check(degree(n)>=k,('threshold entry',k,n))
    check(n==1 or degree(n-1)<k,('previous vertex',k,n))
    # Every earlier block has been checked empty by the block-search loop.
check([threshold(k)[1] for k in range(1,7)]==[3,13,33,61,99,163], 'threshold values')

centers=[]
for n in range(2,1_000_001):
    birth=square(2*n-2) or square(2*n-1)
    death=square(n)
    check(degree(n)-degree(n-1)==int(birth)-int(death),('difference',n))
    check(not (square(2*n-2) and square(2*n-1)),('two births',n))
    if birth and death:
        centers.append(n)
        r=isqrt(n)
        if square(2*n-2):
            q=isqrt(2*n-2)
            check(q%2==0 and r*r-2*(q//2)**2==1,('plus class',n,q))
        else:
            q=isqrt(2*n-1)
            check(q%2==1 and q*q-2*r*r==-1,('minus class',n,q))
        check(above(n),('plateau',n))
check(centers[:6]==[9,25,289,841,9801,28561],('first centers',centers))
check([degree(i) for i in (8,9,10,24,25,26)]==[1,1,1,2,2,2], 'small plateau degrees')

# Multiply the units as integer pairs, checking each entire interval condition.
r,t=3,2
s,v=7,5
for j in range(1,201):
    check(r*r-2*t*t==1 and 2*r*r-2==(2*t)**2 and above(r*r),('plus',j))
    check(s*s-2*v*v==-1 and 2*v*v-1==s*s and above(v*v),('minus',j))
    r,t=3*r+4*t,2*r+3*t
    s,v=3*s+4*v,2*s+3*v

# Verify the exact claimed G_15 path and unique unoriented Hamilton path.
path=(8,1,15,10,6,3,13,12,4,5,11,14,2,7,9)
check(sorted(path)==list(range(1,16)), 'path permutation')
check(all(square(x+y) for x,y in zip(path,path[1:])), 'path edges')
check([path[i]+path[i+1] for i in range(3)]==[9,16,25], 'initial edge labels')
adj=[sum(1<<j for j in range(15) if i!=j and square(i+j+2)) for i in range(15)]
states={(1<<i,i):1 for i in range(15)}
for count in range(1,15):
    nxt={}
    for (mask,end),ways in states.items():
        options=adj[end]&~mask
        while options:
            bit=options&-options; options-=bit
            key=(mask|bit,bit.bit_length()-1)
            nxt[key]=nxt.get(key,0)+ways
    states=nxt
check(sum(states.values())==2,('G15 ordered Hamilton paths',sum(states.values())))

# Hostile plausible misreadings must all fail, with concrete counterexamples.
check(degree(9)==1 and degree(3)==1 and 3<9, 'hostile: center is first degree 1')
check(degree(25)==2 and degree(13)==2 and 13<25, 'hostile: center is first degree 2')
check(square(2*9-2) and not square(2*9-1), 'hostile: 9 birth is odd')
check(square(2*25-1) and not square(2*25-2), 'hostile: 25 birth is even')
check(square(3*3+1)==False and square(2*3*3+1)==False, 'hostile: adjacent square birth')
check(square(4) and square(1) and 4-1==3 and 2*2*2-2>1, 'difference-3 exception')
check(degree(3)==1 and degree(4)==0, 'hostile: degree monotone')
check(9+16==25 and (path[0],path[-1])==(8,9), 'hostile: path endpoints 8 and 10')
print('PASS: direct degrees <=4000; exact threshold blocks k<=1000; recurrence and classification <=1000000; 200 terms per Pell family; G15 uniqueness; hostile cases; no assert')
