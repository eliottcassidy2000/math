#!/usr/bin/env python3
"""Independent exact controls for induced-book height; no producer imports."""
from math import isqrt, prod
from itertools import combinations
import json
checks=0

def ck(ok,label):
    global checks
    checks+=1
    if not ok:raise RuntimeError(label)

def square(n):return n>=0 and isqrt(n)**2==n

def primes_after_three(count):
    out=[];n=5
    while len(out)<count:
        if all(n%d for d in range(2,isqrt(n)+1)):out.append(n)
        n+=2
    return out

def divisors_squarefree(ps):
    ds=[1]
    for p in ps:ds += [d*p for d in ds]
    return sorted(ds)

ps=primes_after_three(14); records=[];pair_checks=0
for k in range(4,15):
    chosen=ps[:k];Q=prod(chosen);P=max(chosen);ds=divisors_squarefree(chosen)
    ck(Q>=5005 and Q%6 in (1,5),'Q domain')
    h=isqrt(isqrt(k**3))
    while h**4<k**3:h+=1
    ck((h-1)**4<k**3<=h**4,'exact h ceiling')
    ck(2*k<=h*h,'finite Chebyshev fraction at most 1/8')
    central=[a for a in ds if Q<=4*a*a and a*a<=4*Q]
    low=[a for a in ds if 4*a*a<Q]
    high=[a for a in ds if a*a>4*Q]
    ck(len(low)==len(high),'reciprocal low-high pairing')
    ck(4*len(low)>=len(ds),'quarter low divisors')
    rest=divisors_squarefree(chosen[1:])
    for a in rest:
        count=sum(Q<=4*d*d and d*d<=4*Q for d in (a,5*a))
        ck(count<=1,'toggle-five central exclusion')
    ck(len(central)+len(low)+len(high)==2**k,'divisor partition')
    retained=[a for a in low if a*a*P**(2*h)>=Q]
    discarded=[a for a in ds if a*a*P**(2*h)<Q]
    ck(4*h*h*len(discarded)<=k*2**k,'finite exact Chebyshev control')
    ck(len(retained)>=2**(k-3),'retained divisor count')
    S=isqrt(9*Q)+1
    while S%2==0 or S%3==0:S+=1
    ck(S*S>9*Q and (S-4)**2<9*Q,'source square within gap four')
    u=(S*S-9*Q)//2;v=u+9*Q;N=10*Q*P**(2*h)
    ck(u>0 and v<=10*Q,'positive bounded anchors')
    ck(u+v==S*S and u%3==v%3==2,'anchor edge and residues')
    leaves=[]
    for a in retained:
        b=Q//a
        ck(a*b==Q and a%2==b%2==1 and b>a,'factor domain')
        X=9*((b-a)//2)**2-u
        ck(0<X<=N,'positive bounded leaf')
        ck(u+X==(3*(b-a)//2)**2,'first anchor square')
        ck(v+X==(3*(b+a)//2)**2,'second anchor square')
        ck(X%3==1 and X not in (u,v),'leaf collision and nonedge residue')
        leaves.append(X)
    ck(len(set(leaves))==len(leaves),'injective leaves')
    ck(all(x>y for x,y in zip(leaves,leaves[1:])),'strict monotonicity')
    if k<=8:
        for x,y in combinations(leaves,2):
            ck(not square(x+y),'direct leaf nonedge')
            pair_checks+=1
    records.append({'k':k,'low_divisors':len(low),'retained':len(retained),
                    'guaranteed_leaves':2**(k-3),'height_bits':max([u,v]+leaves).bit_length()})

# Direct factorization injection for every bounded pair of anchors and
# every common square-sum neighbor in this finite interval.
neighbor_count=0
for u in range(1,41):
    for v in range(u+1,81):
        D=v-u;seen={}
        for x in range(1,201):
            if square(u+x) and square(v+x):
                p=isqrt(u+x);q=isqrt(v+x);a=q-p;b=q+p
                ck(a*b==D and a>0 and b>a,'neighbor-to-divisor factorization')
                ck(a not in seen,'factorization injectivity')
                seen[a]=x;neighbor_count+=1
        tau=sum(D%d==0 for d in range(1,D+1))
        ck(len(seen)<=tau,'bounded common-neighbor divisor bound')

# Exact parity capacities, independently counted from positive strict factors.
for D in range(1,2001):
    factors=[a for a in range(1,isqrt(D)+1) if D%a==0]
    exact=sum(a<D//a and (a+D//a)%2==0 for a in factors)
    t=D if D%2 else D//4 if D%4==0 else 0
    if t:
        tau=2*sum(t%d==0 for d in range(1,isqrt(t)+1))-int(square(t))
        formula=tau//2
    else:formula=0
    ck(exact==formula,'exact parity capacity')

# Integer ceiling log_2 avoids floating errors at powers of two.
for m in list(range(1,1000))+[2**j-1 for j in range(1,101)]+[2**j for j in range(1,101)]+[2**j+1 for j in range(1,101)]:
    k=max(4,(m-1).bit_length()+3)
    ck(k>=4 and 2**(k-3)>=m,'all-m k selection')

print(json.dumps({'status':'PASS','checks':checks,'prime_products_k_range':[4,14],
                  'exact_full_leaf_pair_checks':pair_checks,
                  'bounded_common_neighbors':neighbor_count,'construction_records':records,
                  'scope':'Finite independent arithmetic controls only; no producer imports. The asymptotic is supplied by the reviewed proof.'},indent=2,sort_keys=True))
