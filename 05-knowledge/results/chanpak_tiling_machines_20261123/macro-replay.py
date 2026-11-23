#!/usr/bin/env python3
"""Finite exhaustive replay for the guarded macro residue automata."""
from collections import deque
from math import gcd

def mul(x,y,q):
    a,b,d,g=x
    A,B,D,G=y
    return ((A+a*D)%q,(b+d*B)%q,d*D%q,(g+G+a*B)%q)

def multiplicative_order_2(q):
    k,x=1,2%q
    while x!=1:
        k,x=k+1,2*x%q
    return k

def minimize(q,tiles,r):
    states=[(a,g) for a in range(q) for g in range(q)]
    cls={x:int(x[1]==(r-1)%q) for x in states}
    while True:
        index={}; nxt={}
        for a,g in states:
            signature=(cls[(a,g)], *(cls[((d*a+A)%q,(g+G+a*b)%q)] for A,b,d,G in tiles))
            nxt[(a,g)]=index.setdefault(signature,len(index))
        if len(index)==len(set(cls.values())):
            return len(index)
        cls=nxt

for q in (3,5,7,9,15,21,25,27,33):
    tiles=tuple(tuple(x%q for x in t) for t in ((2,2,4,2),(6,4,8,4),(4,6,8,4)))
    base=(0,0,1,0)
    distance={base:0}
    queue=deque([base]); period=0
    while queue:
        x=queue.popleft()
        for t in tiles:
            y=mul(x,t,q)
            if y not in distance:
                distance[y]=distance[x]+1
                queue.append(y)
            period=gcd(period,distance[x]+1-distance[y])
    assert period==gcd(q,3)
    assert len(distance)==multiplicative_order_2(q)*q**3//gcd(q,3)
    assert len({(a,g) for a,b,d,g in distance})==q*q
    if q%3==0:
        assert all((a-b)%3==0 if d%3==1 else (a+b)%3==1 for a,b,d,g in distance)
        assert all((g-k-int(a%3==2))%3==0 for (a,b,d,g),k in distance.items())
    for r in range(q):
        assert minimize(q,tiles,r)==q*q
    print(q,len(distance),q*q,period)

# Direct history count versus integer macro transfer, with no signed machinery.
from itertools import product
names=("OI","OOI","OII")
tiles=((2,2,4,2),(6,4,8,4),(4,6,8,4))
checked=0
for k in range(5):
    for word in product(range(3),repeat=k):
        letters=''.join(names[j] for j in word)
        sidecar_count=0
        for side in product((0,1),repeat=len(letters)):
            L=[j for j,x in enumerate(side) if not x]
            R=[j for j,x in enumerate(side) if x]
            sidecar_count+=(not L or letters[L[-1]]=='I') and (not R or letters[R[0]]=='O')
        state=(0,0,1,0)
        for j in word:
            state=mul(state,tiles[j],1<<40)
        assert 1+state[3]==sidecar_count
        checked+=1
print(checked,'exact Hamilton sidecar checks')
