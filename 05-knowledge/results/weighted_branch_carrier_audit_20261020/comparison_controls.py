#!/usr/bin/env python3
"""Small exact controls of the actual and the non-existent proposed bridges."""
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations
from math import factorial, gcd
import json

def check(x,s):
    if not x:raise RuntimeError(s)
def parts(n,lo=1):
    if not n:yield ()
    for k in range(lo,n+1):
        for rest in parts(n-k,k):yield (k,)+rest

def v2(k):
    e=0
    while k%2==0:k//=2;e+=1
    return e

def branch(n):
    odd=lev=Fraction(0)
    for mu in parts(n):
        c=Counter(mu);z=1
        for l,m in c.items():z*=l**m*factorial(m)
        o=sum(l//2 for l in mu)+sum(gcd(a,b) for a,b in combinations(mu,2))
        if all(l%2 for l in mu):odd+=Fraction(2**(o-len(mu)+1),z)
        elif all(l%2==0 for l in mu) and len({v2(l) for l in mu})==1:lev+=Fraction(2**(o-len(mu)),z)
    return str(odd),str(lev)

def carry(w):
    b=0;p=0
    # Recurrence constructs the composed map in chronological order.
    for i,e in enumerate(w):
        if e=='1':b=3*b+2**i;p+=1
    return b,2**len(w)-3**p

def necklaces(q,p):
    seen={}
    for ones in combinations(range(q),p):
        bits=['0']*q
        for i in ones:bits[i]='1'
        w=''.join(bits);rep=min(w[r:]+w[:r] for r in range(q))
        if rep in seen:continue
        b,d=carry(w);seen[rep]=gcd(abs(d),b)
    return Counter(seen.values())

branches={str(n):dict(zip(('N_odd','N_lev'),branch(n))) for n in range(1,17)}
check(branches['2']=={'N_odd':'1/2','N_lev':'1/2'},'n2 boundary')
check([branches[str(n)]['N_lev'] for n in (4,6,8)]==['1','2','15'],'branch census')
# K3 graph switching has triangle parity 0/1 invariant under every permutation.
# The tournament torsor has the same two states, but odd permutations exchange them.
p=list(permutations(range(3)))
parity=lambda x:sum(x[i]>x[j] for i in range(3) for j in range(i+1,3))%2
regular=[x for x in p if parity(x)==0]
linear_total=sum(2 for _ in p)//6;affine_total=sum(2 if parity(x)==0 else 0 for x in p)//6
regular_branch=Fraction(sum(2 for _ in regular),6)
check((linear_total,affine_total,regular_branch)==(2,1,1),'linear affine loss')
w='11101000100';b,d=carry(w);check((b,d,gcd(b,d))==(475,1805,95),'valuation hostile')
counts=necklaces(15,8);check(counts=={1:418,73:10,359:1},'common witness hostile')
print(json.dumps({'partition_branches':branches,'n3_gain_comparison':{'linear_full_orbits':linear_total,'affine_tournament_full_orbits':affine_total,'common_regular_branch':str(regular_branch),'linear_complement':str(linear_total-regular_branch),'tournament_N_lev':str(affine_total-regular_branch)},'collatz_prime_depth_hostile':{'word':w,'carry':b,'clock':d,'gcd':gcd(b,d),'reduced_denominator':d//gcd(b,d),'v19_carry':1,'v19_clock':2},'collatz_same_word_hostile':{'q':15,'ones':8,'clock':2**15-3**8,'necklace_content_counts':dict(sorted(counts.items()))}},indent=2,sort_keys=True))
