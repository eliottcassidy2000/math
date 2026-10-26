# The Legendre relation on primes: for p,q = 3 (mod 4), QR gives (p/q) = -(q/p): a tournament on P3 = {primes = 3 mod 4};
# for p or q = 1 (mod 4) the relation is symmetric: a graph. Census of the tournament on the first N primes of P3.
import sys, json
from itertools import combinations
from collections import Counter
def primes_upto(n):
    s=bytearray([1])*(n+1); s[0]=s[1]=0
    for i in range(2,int(n**0.5)+1):
        if s[i]: s[i*i::i]=bytearray(len(s[i*i::i]))
    return [i for i in range(n+1) if s[i]]
def legendre(a,p):
    a%=p
    if a==0: return 0
    return 1 if pow(a,(p-1)//2,p)==1 else -1
N=int(sys.argv[1]) if len(sys.argv)>1 else 200
P=[p for p in primes_upto(200000) if p%4==3][:N]
# tournament: p -> q iff (p/q) = 1  (p is a square mod q)
arc={}
for p,q in combinations(P,2):
    s=legendre(p,q); t=legendre(q,p)
    assert s==-t, 'QR antisymmetry on P3'
    arc[(p,q)] = (s==1)
def wins(p,q): return arc[(p,q)] if p<q else not arc[(q,p)]
# score sequence, 3-cycle density, 4-vertex census
scores=Counter()
for p in P:
    scores[p]=sum(wins(p,q) for q in P if q!=p)
c3=0; triples=0
for a,b,c in (combinations(P,3) if N>=4 else []):
    triples+=1
    w=[wins(a,b),wins(b,c),wins(c,a)]
    if all(w) or not any(w): c3+=1
def cls4(v):
    # classes by score multiset: TT4 (3,2,1,0); source+C3 (3,1,1,1); sink+C3 (2,2,2,0); strong (2,2,1,1)
    sc=tuple(sorted(sum(wins(x,y) for y in v if y!=x) for x in v))
    return {(0,1,2,3):'TT4',(1,1,1,3):'source+C3',(0,2,2,2):'C3+sink',(1,1,2,2):'strong'}[sc]
census=Counter(cls4(v) for v in (combinations(P,4) if N>=4 else []))
tot=max(1,sum(census.values()))
exp={'TT4':3/8,'source+C3':1/8,'C3+sink':1/8,'strong':3/8}
if N>=4: print(json.dumps({'N':N,'largest_prime':P[-1],'triangle_density':round(c3/max(1,triples),5),'random_expectation':0.25,
 'four_vertex_census':{k:[census[k],round(census[k]/tot,4),exp[k]] for k in exp},
 'all_four_classes_occur':all(census[k]>0 for k in exp),
 'score_spread':[min(scores.values()),max(scores.values()),(N-1)/2]}))

def sampled(N, samples=300000, seed=20261025):
    import random
    random.seed(seed)
    P=[p for p in primes_upto(3000000) if p%4==3][:N]
    idx={p:i for i,p in enumerate(P)}
    cache={}
    def w(p,q):
        key=(p,q) if p<q else (q,p)
        if key not in cache: cache[key]=(legendre(key[0],key[1])==1)
        return cache[key] if p<q else not cache[key]
    census=Counter(); c3=0
    for _ in range(samples):
        v=random.sample(P,4)
        sc=tuple(sorted(sum(w(x,y) for y in v if y!=x) for x in v))
        census[{(0,1,2,3):'TT4',(1,1,1,3):'source+C3',(0,2,2,2):'C3+sink',(1,1,2,2):'strong'}[sc]]+=1
        a,b,c=v[:3]; ww=[w(a,b),w(b,c),w(c,a)]
        c3+= (all(ww) or not any(ww))
    tot=sum(census.values())
    print(json.dumps({'sampled_N':N,'largest_prime':P[-1],'samples':samples,'triangle_density':round(c3/samples,4),
                      'four_vertex_frequencies':{k:round(census[k]/tot,4) for k in exp},'expected':exp}))
if len(sys.argv)>2: sampled(int(sys.argv[2]))
