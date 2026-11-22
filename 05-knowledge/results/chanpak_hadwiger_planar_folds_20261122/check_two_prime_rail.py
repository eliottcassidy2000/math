"""Exact transfer and simultaneous-phase replay for the 11010 1^k crossed rail."""
from collections import defaultdict
from itertools import product

def stream(tags, swapped, forward_arc):
    d={(True,0,False):1}
    for p,q in tags:
        if swapped:p,q=q,p
        e=defaultdict(int)
        for (left,mid,right),n in d.items():
            e[not p,mid,right]+=n
            if mid or p:e[left,1 if q else 2,right]+=n
            if right or q:e[left,mid,True]+=n
        d=e
    return sum(n for (left,mid,right),n in d.items()
               if left and (mid==2 or (mid==0 and forward_arc)))

def counts(k,arc):
    W=tuple(map(int,"11010"+"1"*k))
    tags=tuple((v,v) for v in W)+((0,0),(1,0),(0,1),(0,0))+tuple((1-v,1-v) for v in reversed(W))
    return stream(tags,False,arc)+stream(tags,True,not arc)

def midpoint(X,Z):
    return 2910-59246*X+18304*Z+459162*X*X-1498068*X*Z+3312738*Z*Z

def phase(s,t,r2,r3):
    assert r2%32==0 and r3%3==1
    B,C=2**s,3**t
    u=0
    for i in range(s-5):
        if ((104186+1660945*u+6625476*u*u-r2//32)>>i)&1:u+=1<<i
    want=1+8*u
    m2=0
    for i in range(3,s):
        if ((pow(9,m2,1<<(i+1))-want)>>i)&1:m2+=1<<(i-3)
    m3=0
    for i in range(1,t):
        for a in range(3):
            candidate=m3+a*3**(i-1)
            x=pow(4,candidate,3**(i+1))
            if midpoint(x,0)%3**(i+1)==r3%3**(i+1):
                m3=candidate
                break
        else: raise AssertionError
    b,c=2**(s-3),3**(t-1)
    m=m2+b*((m3-m2)*pow(b,-1,c)%c) if c>1 else m2
    return 2*(m+b*c*max(s,t))

for k in range(1,10):
    m=midpoint(2**k,3**k)
    assert (counts(k,True),counts(k,False))==(m-1,m+1)
assert all(z%11==0 for z in (-59246,18304,459162,-1498068,3312738))
E=2**60*(2**61-1)
assert 2910%11==6 and E%11==1
for s,t in product((5,6,14,25,42),(1,2,9,17)):
    B,C=2**s,3**t
    for r2,r3 in ((0,1),(64,16),(E%B,E%C)):
        k=phase(s,t,r2,r3)
        for q,R in ((B,r2),(C,r3)):
            assert midpoint(pow(2,k,q),pow(3,k,q))%q==R%q
assert all(midpoint(pow(2,k,15),pow(3,k,15))%15 !=
           (pow(2,p-1,15)*(pow(2,p,15)-1))%15
           for k,p in product(range(1,9),range(2,10)))
print("simultaneous dyadic-ternary phase replay PASS")
