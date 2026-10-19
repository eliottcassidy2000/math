#!/usr/bin/env python3
"""Targeted exact controls for projective Gold budgets and signed Golay action."""
from itertools import combinations, product
from fractions import Fraction
from collections import Counter
import json
checks=0
def require(ok,msg):
    global checks
    checks+=1
    if not ok: raise RuntimeError(msg)
def mob(x,a,b,c,d,p):
    if x is None: num,den=a%p,c%p
    else: num,den=(a*x+b)%p,(c*x+d)%p
    return None if not den else num*pow(den,-1,p)%p
def g(x,p):return mob(x,-2,1,-1,1,p)
def z(x):return mob(x,1,-5,1,-9,11)
def zi(x):return mob(x,9,-5,1,-1,11)
def orbits(domain,fun):
    todo=set(domain);out=[]
    while todo:
        a=next(iter(todo));v=a;cyc=[]
        while v not in cyc:cyc.append(v);todo.remove(v);v=fun(v)
        require(v==a,'permutation cycle closes')
        out.append(len(cyc))
    return sorted(out)
p1=list(range(11))+[None]
for x in p1:
    require(zi(z(x))==x,'cross ratio inverse')
    require(z(g(x,11))==(None if z(x) is None else 2*z(x)%11),'one step')
    require(z(g(g(x,11),11))==(None if z(x) is None else 4*z(x)%11),'two steps')
    ax=mob(x,2,7,1,6,11)
    require(z(ax)==(None if z(x) is None else (z(x)+1)%11),'translation')
carrier=[x for x in p1 if x!=9]
require(orbits(carrier,lambda x:g(x,11))==[1,10],'one-step cycles')
require(orbits(carrier,lambda x:g(g(x,11),11))==[1,5,5],'two-step cycles')
points=[Fraction(0),*[Fraction(x,13) for x in (1,2,4,8,16)],*[Fraction(x) for x in (-5,-14,-7,-20,-10)]]
labels=[0,3,9,5,4,1,8,10,7,6,2]
kappa=dict(zip(points,labels))
def C(x): return 3*x+1 if x.numerator%2 else x/2
for q in points: require(kappa[C(q)]==4*kappa[q]%11,'standard Collatz carrier')
# The classical cyclic ternary Golay factor and its canonical small errors.
POLY=(2,0,1,2,1,1);ZERO=(0,)*5;ONE=(1,0,0,0,0);A=(0,1,0,0,0)
def add(a,b):return tuple((x+y)%3 for x,y in zip(a,b))
def scale(c,a):return tuple(c*x%3 for x in a)
def mul(a,b):
    v=[0]*9
    for i,x in enumerate(a):
        for j,y in enumerate(b):v[i+j]=(v[i+j]+x*y)%3
    for i in range(8,4,-1):
        c=v[i]
        for j in range(6):v[i-5+j]=(v[i-5+j]-c*POLY[j])%3
    return tuple(v[:5])
def power(a,n):
    ans=ONE
    while n:
        if n&1:ans=mul(ans,a)
        a=mul(a,a);n//=2
    return ans
AP=[power(A,i) for i in range(11)]
def S(v,sign=1):
    ans=ZERO
    for i,c in enumerate(v):ans=add(ans,scale(c,AP[sign*i%11]))
    return ans
E=[(0,)*11]
for w in (1,2):
    for supp in combinations(range(11),w):
        for cs in product((1,2),repeat=w):
            v=[0]*11
            for i,c in zip(supp,cs):v[i]=c
            E.append(tuple(v))
dec={S(v):v for v in E};require(len(dec)==243,'perfect small-error bijection')
def perm(v,k):
    ans=[0]*11
    for i,c in enumerate(v):ans[k*i%11]=c
    return tuple(ans)
def J(s):return S(perm(dec[s],-1))
def G(s):return J(power(s,9))
for s,v in dec.items():
    require(G(s)==S(perm(v,2)),'budget step on decoder')
    require(G(G(s))==power(s,81),'square equals forward golden Frobenius')
    u=s
    for _ in range(5):u=G(u)
    require(u==J(s),'fifth power reflection')
    require(S(perm(v,2))==power(S(v,-1),9),'reciprocal module first coordinate')
    require(S(perm(v,2),-1)==power(S(v),9),'reciprocal module second coordinate')
# Every coordinate two-plane realizes the sharp dimension-two linear patch.
for i,j in combinations(range(11),2):
    basis=(AP[i],AP[j]);image=[]
    for a,b in product(range(3),repeat=2):
        s=add(scale(a,basis[0]),scale(b,basis[1]));image.append(s)
        require(J(s)==add(scale(a,AP[-i%11]),scale(b,AP[-j%11])),'linear two-plane')
    require(len(set(image))==9,'two-plane dimension')
# Prime-family cycle classification, no conjecture census.
prime_cases={}
for p in range(2,102):
    if any(p%d==0 for d in range(2,int(p**.5)+1)):continue
    cyc=orbits(list(range(p))+[None],lambda x:g(g(x,p),p))
    if p==2:require(cyc==[3],'p2')
    elif p==5:require(cyc==[1,5],'p5')
    else:
        split=pow(5,(p-1)//2,p)==1
        if split:
            require(cyc.count(1)==2,'split fixed points')
            non=[x for x in cyc if x>1];require(len(set(non))==1,'split equal lengths')
            require((p-1)//2%non[0]==0,'split square subgroup')
        else:
            require(len(set(cyc))==1 and cyc[0]>1,'inert equal lengths')
            require((p+1)//2%cyc[0]==0,'inert norm-one square subgroup')
    prime_cases[p]=cyc
require(prime_cases[11]==[1,1,5,5],'eleven family')
require(prime_cases[19]==[1,1,9,9],'nineteen family')
require(prime_cases[31]==[1,1,15,15],'thirty-one family')
print(json.dumps({'checks':checks,'eleven_budget_cycles':[1,10], 'eleven_squared_cycles':[1,5,5], 'syndrome_states':243,'linear_coordinate_planes':55,'prime_squared_cycle_partitions':prime_cases},indent=2))
