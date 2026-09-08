#!/usr/bin/env python3
"""Complete type-(5)(2) pair certificate and uniform H4 block argument.

RESERVED pending independent audit.  A pair has support union at most14;
no fourteen-label bound is imposed on a four-generator action.  The
global proof uses the strict majority of the long block, not a tuple census.
"""
from collections import Counter
from hashlib import sha256
from itertools import combinations,permutations
from math import comb
import json

gates=0
def need(value,label):
    global gates
    gates+=1
    if not value:raise RuntimeError(label)
def mul(p,q):return tuple(p[x] for x in q)
def power(p,n):
    q=tuple(range(len(p)))
    for _ in range(n):q=mul(q,p)
    return q
def cycle(p,c):
    q=list(p)
    for x,y in zip(c,c[1:]+c[:1]):q[x]=y
    return tuple(q)
def cycles(p):
    seen=set();result=[]
    for x in range(len(p)):
        if x in seen or p[x]==x:continue
        row=[];y=x
        while y not in seen:row.append(y);seen.add(y);y=p[y]
        result.append(tuple(row))
    return tuple(sorted(result,key=lambda c:(-len(c),c)))
def flags(p,q):
    pq=mul(p,q);qp=mul(q,p)
    return (mul(pq,p)==mul(qp,q),
            mul(mul(pq,pq),p)==mul(mul(qp,qp),q),pq==qp)
def letters(p,q):
    n=len(p)
    ordinary=all(p[q[p[x]]]==q[p[q[x]]] for x in range(n))
    fifth=all(p[q[p[q[p[x]]]]]==q[p[q[p[q[x]]]]] for x in range(n))
    commuting=all(p[q[x]]==q[p[x]] for x in range(n))
    return ordinary,fifth,commuting
def joint_orbits(p,q):
    seen=set();result=[]
    for x in range(len(p)):
        if x in seen:continue
        seen.add(x);todo=[x]
        for y in todo:
            for g in (p,q):
                z=g[y]
                if z not in seen:seen.add(z);todo.append(z)
        if len(todo)>1:result.append(tuple(sorted(todo)))
    return tuple(sorted(result,key=lambda c:(len(c),c)))

D=14;I=tuple(range(D))
sigma=cycle(cycle(I,(0,1,2,3,4)),(5,6))
A5=set(range(5));A2={5,6}
expected=(
 {((3,2,2,0),(7,)):10,
  ((4,0,0,1),(3,6)):420,
  ((4,0,0,2),(2,6)):35,
  ((5,0,0,1),(3,5)):84,
  ((5,0,0,2),(2,5)):6},
 {((2,1,1,0),(10,)):2100,
  ((4,0,0,2),(2,6)):105,
  ((4,1,1,1),(7,)):10,
  ((5,0,0,2),(2,5)):6},
 {((0,0,0,0),(2,2,5,5)):504,
  ((0,0,0,2),(2,5,5)):504,
  ((5,0,0,0),(2,2,5)):84,
  ((5,0,0,2),(2,5)):4},
)
tables=[Counter() for _ in range(3)]
admitted=[[],[],[]]
data={}
raw=sha256();total=0
# Unique unfiltered construction: pick the5-set, fix its least element
# as the cyclic first entry, take all24 cyclic orders, then an unordered
# disjoint2-set.  No support, parity, graph or braid filter is applied.
for five in combinations(range(D),5):
    remaining=tuple(x for x in range(D) if x not in five)
    for order in permutations(five[1:]):
        base=cycle(I,(five[0],)+order)
        for two in combinations(remaining,2):
            tau=cycle(base,two);total+=1
            ff=flags(sigma,tau)
            need(ff==letters(sigma,tau),'independent literal pair relation check')
            # Fixed-width stream: fourteen image bytes, then three flags.
            raw.update(bytes(tau));raw.update(bytes(ff))
            if not any(ff):continue
            cc=cycles(tau)
            need(tuple(map(len,cc))==(5,2),'actual admitted cycle type')
            B5,B2=map(set,cc)
            need(B5==set(five) and B2==set(two),'unique length gauge matches construction')
            cell=(len(A5&B5),len(A5&B2),len(A2&B5),len(A2&B2))
            oo=joint_orbits(sigma,tau);sizes=tuple(map(len,oo))
            data[tau]=(cell,sizes)
            for j,yes in enumerate(ff):
                if not yes:continue
                tables[j][cell,sizes]+=1;admitted[j].append(tau)
                need((cell,sizes) in expected[j],'exact typed matrix and joint orbit coupling')
            if ff[0]:
                need(cell[0]>=3,'ordinary strict majority in the five-block')
                need(len(sizes)==len(set(sizes)),'ordinary nontrivial joint sizes are distinct')
            if ff[1]:
                need(cell[0]>=2,'fifth positive five-block intersection')
            if ff[2]:
                need(not(A5&B2) and not(A2&B5),'commuting cross-length blocks do not mix')
                need(A5==B5 or not(A5&B5),'commuting five-block equal or disjoint')
                need(A2==B2 or not(A2&B2),'commuting two-block equal or disjoint')
                need(tuple(tau[x] for x in sorted(A5)) in
                     [tuple(power(sigma,j)[x] for x in sorted(A5)) for j in range(5)],
                     'commuting restriction is a power on the original five-cycle')

need(total==comb(14,5)*24*comb(9,2)==comb(14,7)*504==1729728,
     'full unfiltered pair universe counted in two ways')
for j in range(3):
    need(dict(tables[j])==expected[j],'complete admitted typed table')
need(tuple(map(len,admitted))==(555,2221,1096),'all admitted totals')
need(1096==comb(7,5)*24+comb(7,5)*24+4*comb(7,2)+4,
     'independent commuting centralizer count')

# Structural inequalities used in the all-degree proof.  These are
# integer checks of the written set argument, not a finite-action proof.
for overlap in (3,4,5):
    need(5-overlap<3,'disjoint from A5 cannot have ordinary majority in B5')
for overlap in (2,4,5):
    need(overlap>0,'fifth overlap defeats the commuting disjoint alternative')
# The pair inventory also retains the unused joint-centralizer sidecar:
# a semiregular nontrivial restriction of a global (5)(2) permutation
# can occupy only a two-point or a five-point transitive orbit.
for size in range(2,15):
    possible=[typ for typ in ((2,),(5,),(5,2))
              if sum(typ)==size and len(set(typ))==1]
    need(bool(possible)==(size in (2,5)),
         'all semiregular restriction sizes for the exact global cycle type')

# After the common five-set has been PROVED invariant, c and d restrict
# to powers of a there.  Every nonidentity power is still a5-cycle.
long=cycle(tuple(range(5)),(0,1,2,3,4))
for i in range(1,5):
    for j in range(1,5):
        pp=power(long,i);qq=power(long,j)
        need(flags(pp,qq)[2],'two centralizer powers commute')
        need(flags(pp,qq)[1]==(i==j),'odd fifth relation on commuting five-cycle powers')

# On the complementary invariant subset only one transposition remains.
# Every pair has support union<=4, so this small pair universe is complete.
I4=tuple(range(4));s2=cycle(I4,(0,1));orders=set()
for two in combinations(range(4),2):
    t2=cycle(I4,two);product=mul(s2,t2)
    order=next(k for k in (1,2,3) if power(product,k)==I4)
    orders.add(order)
    need(flags(s2,t2)[1]==(s2==t2),
         'complete complementary transposition fifth relation')
    need(flags(s2,t2)[1]==(power(product,5)==I4),
         'involutive odd relation and product power agree')
need(orders=={1,2,3},'all transposition pair product orders')

# Minimal seven-label hostiles to premature extraction of the2-cycle.
# These are NOT full H4 quadruples.
hostile3=cycle(cycle(I,(0,5,3,6,1)),(2,4))
hostile5=cycle(cycle(I,(0,3,2,5,1)),(4,6))
hostiles=[]
small_support=set(range(7))
for j,tau in ((0,hostile3),(1,hostile5)):
    need(flags(sigma,tau)[j],'full hostile satisfies the required odd relation')
    extracted=cycle(I,cycles(tau)[1])
    need(not flags(cycle(I,(5,6)),extracted)[j],
         'isolated transpositions lose the odd relation')
    need(set().union(*map(set,cycles(tau)),A5,A2)==small_support,
         'hostile attains minimal ambient support seven')
    hostiles.append((3 if j==0 else 5,cycles(tau),data[tau]))
# Same long block alone does not force equality; cyclic orders are kept.
same_block_fifth=next(tau for tau in admitted[1]
                      if tau!=sigma and data[tau][0]==(5,0,0,2))
need(flags(sigma,same_block_fifth)[1] and not flags(sigma,same_block_fifth)[2],
     'same-five-block fifth hostile needs the global commuting sidecar')
need(flags(sigma,sigma)==(True,True,True),'equal quadruple positive control')
need(tuple(map(len,joint_orbits(sigma,sigma)))==(2,5),
     'equal quadruple already has two distinct nontrivial orbits')

semantic=sha256(json.dumps([sorted(t.items()) for t in tables],
                          separators=(',',':')).encode()).hexdigest()
print('h4_mixed52: PASS')
print('scope: uniform equality for all H4 generators of type(5)(2); no ambient tuple bound')
print('full unfiltered fourteen-label pair universe:',total)
for label,j in (('ordinary',0),('fifth',1),('commuting',2)):
    print(label,'partners:',len(admitted[j]))
    print(label,'typed matrix / joint orbit table:',sorted(tables[j].items()))
print('minimal cycle-extraction hostiles:',hostiles)
print('same-five-block fifth hostile:',cycles(same_block_fifth))
print('gates:',gates)
print('raw fixed-width pair/flags SHA256:',raw.hexdigest())
print('typed table SHA256:',semantic)

