#!/usr/bin/env python3
"""Session69 exact powers/doubling controls; stdlib only, no external data."""
from collections import Counter
from itertools import permutations
from math import gcd, lcm
import json

checks = 0
def require(p, msg=''):
    global checks
    checks += 1
    if not p: raise AssertionError(msg)

def odd_data(n):
    a=0
    while n%2==0: a+=1; n//=2
    u=n
    es=[]; fs=[]; p=3
    while p*p<=n:
        e=0
        while n%p==0: e+=1; n//=p
        if e: es.append(e);fs.append((p,e))
        p+=2
    if n>1: es.append(1);fs.append((n,1))
    if u==1:return a,1,0
    e=0
    for z in es:e=gcd(e,z)
    b=1
    for p,z in fs:b*=p**(z//e)
    return a,b,e

def sq_data(n):
    a,b,e=odd_data(n)
    if b==1:return a,1,0
    t=0
    while e%2==0:t+=1;e//=2
    return a,b**e,t

def row(n,cap,square=False):
    vals=set();v=n
    while v<=cap:
        z=v
        while z<=cap:vals.add(z);z*=2
        if n==1:break
        v=v*v if square else v*n
    return vals

def common(n,m,square=False):
    a,b,e=(sq_data(n) if square else odd_data(n))
    c,d,f=(sq_data(m) if square else odd_data(m))
    if b!=d:return None
    if b==1:return 2**max(a,c)
    if square:
        t=max(e,f)
        return 2**max(a*2**(t-e),c*2**(t-f))*b**(2**t)
    g=gcd(e,f)
    return 2**max(a*f//g,c*e//g)*b**lcm(e,f)

# Independent value enumeration of origins and multiplicities.
N=5000
counts={False:Counter(),True:Counter()}
origins={False:[None]*(N+1),True:[None]*(N+1)}
for square in [False,True]:
    for n in range(1,N+1):
        for x in row(n,N,square):
            if origins[square][x] is None:origins[square][x]=n
        if n==1:continue # base1 would produce infinite exponent multiplicity
        v=n
        while v<=N:
            x=v
            while x<=N:counts[square][x]+=1;x*=2
            v=v*v if square else v*n
    for x in range(1,N+1):
        a,b,e=(sq_data(x) if square else odd_data(x))
        require(origins[square][x]==b,(square,x,'origin',b))
        if b==1:
            predicted=sum(a//(2**j) for j in range(a.bit_length())) if square else sum(a//r for r in range(1,a+1))
        elif square:
            predicted=sum(a//(2**j)+1 for j in range(e+1))
        else:
            predicted=sum(a//r+1 for r in range(1,e+1) if e%r==0)
        require(counts[square][x]==predicted,(square,x,'multiplicity',predicted,counts[square][x]))

# Complete pairwise equality of truncated intersections with the proved row formula.
cap=10**6
for square in [False,True]:
    rows={n:row(n,cap,square) for n in range(1,121)}
    for n in rows:
        for m in range(n,121):
            L=common(n,m,square)
            require(rows[n]&rows[m]==(set() if L is None else row(L,cap,square)),(square,n,m,L))

# Independent square-later search with canonical ODD chain roots.
canonical_later=[];arbitrary_later=[]
for n in range(1,10001):
    roots=[];z=n*n;k=0
    while z%2==0:
        k+=1;z//=2
        if z>n:roots.append((z,k))
    a,b,e=odd_data(n)
    u=n//2**a
    actual=bool(k>0 and z>n) # final z is the unique odd root
    require(actual==(a>=1 and u>2**a),(n,'odd root'))
    if actual:canonical_later.append(n)
    if roots:arbitrary_later.append(n)
# Exact subtraction bridge to primitive abc triples, with no abc assumption.
def rad(n):
    out=1;p=2
    while p*p<=n:
        if n%p==0:
            out*=p
            while n%p==0:n//=p
        p+=1
    if n>1:out*=n
    return out
abc_records=[]
for n in canonical_later:
    a,_,_=odd_data(n);u=n//2**a
    A=u-2**a;B=2**a;C=u
    require(A>0 and A+B==C and gcd(A,B)==gcd(A,C)==gcd(B,C)==1)
    require(B*C==n and u*u-n==u*A)
    if C>rad(A*B*C):abc_records.append({'n':n,'triple':[A,B,C],'radical':rad(A*B*C)})
require(abc_records[0]=={'n':72,'triple':[1,8,9],'radical':6})
require(canonical_later[0]==6)
require(arbitrary_later[0]==4)
require(common(6,9)==36 and common(6,9,True)==36)
require(common(6,27)==216 and common(6,27,True) is None)
require(min(n for n in range(1,N+1) if origins[False][n]!=origins[True][n])==27)

# Semigroup law and noncommutativity; 3 is an odd support control.
for n in range(1,13):
    for k in range(4):
        for r in range(1,5):
            for ell in range(4):
                for s in range(1,5):
                    require(2**k*(2**ell*n**s)**r==2**(k+r*ell)*n**(r*s))
require((2*3)**2!=2*(3**2))

# Paley7 native edge predicate, complete labelled automorphism/path controls.
Q={1,2,4}
def arc(i,j):return (j-i)%7 in Q
require(Q=={pow(x,2,7) for x in range(1,7)})
require({2*x%7 for x in Q}==Q)
require(Counter((x-y)%7 for x in Q for y in Q if x!=y)==Counter({d:1 for d in range(1,7)}))
auts=[];hp=0
for p in permutations(range(7)):
    if all(arc(p[i],p[i+1]) for i in range(6)):hp+=1
    if all(arc(i,j)==arc(p[i],p[j]) for i in range(7) for j in range(i+1,7)):auts.append(p)
expected={tuple((a*i+b)%7 for i in range(7)) for a in Q for b in range(7)}
require(set(auts)==expected)
require(hp==189)
modules=[]
for mask in range(1,127):
    M=[i for i in range(7) if mask>>i&1]
    if 1<len(M)<7 and all(len({arc(x,m) for m in M})==1 for x in range(7) if x not in M):modules.append(M)
require(modules==[])
# F8 = F2[x]/(x^3+x+1), alpha=x; exact polynomial-basis multiplication.
def mult(a,b):
    c=0
    while b:
        if b&1:c^=a
        b>>=1;a<<=1
        if a&8:a^=11
    return c
powers=[1]
for i in range(1,8):powers.append(mult(powers[-1],2))
require(len(set(powers[:7]))==7 and powers[7]==1)
for i in range(7):require(mult(powers[i],powers[i])==powers[(2*i)%7])
require(powers[1]^powers[2]^powers[4]==0)

result={
 'status':'ALL CHECKS PASSED','checks':checks,
 'universe':{'origin_and_occurrence_values':N,'row_intersection_bases':120,'row_value_cap':cap,'square_later_bases':10000},
 'canonical_later_first_30':canonical_later[:30],
 'canonical_later_count_upto_10000':len(canonical_later),
 'arbitrary_later_first':arbitrary_later[0],
 'abc_binary_term_first_quality_above_one':abc_records[0],
 'abc_binary_term_quality_above_one_count_upto_10000':len(abc_records),
 'first_allpower_vs_square_origin_disagreement':27,
 'sample_rows':[{'n':n,'allpower_origin':origins[False][n],'square_origin':origins[True][n],'allpower_occurrences_excluding_base1':counts[False][n],'square_occurrences_excluding_base1':counts[True][n]} for n in [1,2,3,4,6,8,9,12,15,18,25,27,36,49,64,81,216]],
 'paley7':{'automorphisms':len(auts),'hamiltonian_paths':hp,'proper_nonsingleton_modules':modules,'Q':sorted(Q),'F8_alpha_powers':powers[:7]},
}
print(json.dumps(result,indent=2))
