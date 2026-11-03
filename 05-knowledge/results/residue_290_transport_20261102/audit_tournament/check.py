#!/usr/bin/env python3
"""Root's independent coefficient audit: full set partitions, direct permutations."""
from fractions import Fraction as F
from collections import Counter
from itertools import permutations
from math import gcd, factorial
import json
C=0
def check(ok,msg):
    global C
    C+=1
    if not ok: raise RuntimeError(msg)
def parts(n,low=1):
    if n==0: yield ()
    for x in range(low,n+1,2):
        for tail in parts(n-x,x): yield (x,)+tail
def setparts(n):
    blocks=[]
    def rec(i):
        if i==n:
            yield tuple(tuple(b) for b in blocks); return
        for b in blocks:
            b.append(i); yield from rec(i+1); b.pop()
        blocks.append([i]); yield from rec(i+1); blocks.pop()
    yield from rec(0)
def connected(l):
    result=0
    for p in setparts(len(l)):
        e=sum(gcd(l[i],l[j]) for b in p for i in b for j in b if i<j)
        result+=(-1)**(len(p)-1)*factorial(len(p)-1)*2**e
        check(e>=0,'gcd exponent')
    return result

g=[0]
for n in range(1,10):
    coefficient=F(0)
    for l in parts(n):
        z=1
        for j,m in Counter(l).items(): z*=j**m*factorial(m)
        coefficient+=F(2**((n-len(l))//2)*connected(l),z)
    coefficient*=n
    check(coefficient.denominator==1 and coefficient>0,'positive integral coefficient')
    g.append(int(coefficient))
a=[F(1)]
for n in range(1,10):
    a.append(sum(g[k]*a[n-k] for k in range(1,n+1))/n)
check(a==[1,1,1,2,4,12,56,456,6880,191536],'actual sequence prefix')
# Burnside computed by individual permutations, literal signed pair orbits.
for n in range(1,8):
    total=0
    pairs={(i,j) for i in range(n) for j in range(i+1,n)}
    for p in permutations(range(n)):
        unseen=set(pairs); valid=True; orbit_count=0
        while unseen:
            initial=min(unseen); x,y=initial; parity=0; seen=[]
            while (x,y) not in seen:
                seen.append((x,y)); u,v=p[x],p[y]
                parity^=int(u>v); x,y=sorted((u,v))
            check((x,y)==initial,'pair orbit closes')
            unseen.difference_update(seen); orbit_count+=1
            if parity: valid=False
        total+=2**orbit_count if valid else 0
    check(F(total,factorial(n))==a[n],'permutation Burnside agrees')
# First full-weighted-history hostile; direct unit-tail factors.
def fake_history(core_cycles):
    h=[F(2**(r*(r-1)//2),factorial(r)) for r in range(11)]
    for r in range(2):
        h[9+r]+=F(2**(r*(r-1)//2+core_cycles*r),factorial(r))
    return h
x,y=fake_history(1),fake_history(3)
check(x[:10]==y[:10] and y[10]-x[10]==6,'history hostile exact')
# Actual B weights and exact moment annihilator, independent direct profile sums.
profiles=[]
for mass in range(3,25):
    weights={}
    for l in parts(mass,3):
        z=1
        for j,m in Counter(l).items(): z*=j**m*factorial(m)
        e=(mass-len(l))//2+sum(gcd(l[i],l[j]) for i in range(len(l)) for j in range(i+1,len(l)))
        weights[len(l)]=weights.get(len(l),F(0))+F(2**e,z)
    t=sorted(weights)
    check(t==[j for j in range(1,mass//3+1) if j%2==mass%2],'actual support positive')
    d=len(t); q=[1]
    for j in t:
        r=[F(0)]*(len(q)+1)
        for i,v in enumerate(q): r[i]-=2**j*v; r[i+1]+=v
        q=r
    moments=[sum(weights[j]*2**(j*r) for j in t) for r in range(2*d+4)]
    for r in range(d+4): check(sum(q[i]*moments[r+i] for i in range(d+1))==0,'exact annihilation')
    # Gaussian elimination computes rank directly without Vandermonde formula.
    h=[[moments[i+j] for j in range(d)] for i in range(d)]; rank=0
    for col in range(d):
        pivot=next((row for row in range(rank,d) if h[row][col]),None)
        if pivot is None: continue
        h[rank],h[pivot]=h[pivot],h[rank]; v=h[rank][col]
        h[rank]=[z/v for z in h[rank]]
        for row in range(rank+1,d):
            v=h[row][col]; h[row]=[u-v*z for u,z in zip(h[row],h[rank])]
        rank+=1
    check(rank==d,'Hankel rank')
    profiles.append([mass,d])
print(json.dumps({'status':'PASS','gates':C,'g':g[1:],'a':[int(x) for x in a],
 'universes':{'connected_coefficients_through':9,'individual_permutations_through':7,
              'actual_core_mass':[3,24]},'hankel_dimensions':profiles},indent=2,sort_keys=True))
