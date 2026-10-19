#!/usr/bin/env python3
"""Independent square-scale, dyadic/mod-7, and odd-3-power constructions."""
from itertools import combinations
from math import isqrt
import json, random
CHECKS=0

def check(value,msg):
    global CHECKS
    CHECKS+=1
    if not value: raise RuntimeError(msg)

def square(x): return x>=0 and isqrt(x)**2==x

def audit(original, edges, new, inherited=False):
    labels=original+new
    target=set()
    if inherited:
        target.update((i,j) for i,j in combinations(range(len(original)),2) if square(original[i]+original[j]))
    for e,(i,j) in enumerate(edges):
        target.add((i,len(original)+e)); target.add((j,len(original)+e))
    check(len(labels)==len(set(labels)), 'collision')
    for i,j in combinations(range(len(labels)),2):
        check(square(labels[i]+labels[j])==((i,j) in target),('adjacency mismatch',i,j,labels[i],labels[j]))
    if new: check(min(new)>max(original),'disjoint ranges')
    return max(labels)

def dyadic(n,edges):
    power=8
    while power<64*n: power*=2
    lam=7*power
    orig=[lam*(i+1)+3 for i in range(n)]
    new=[(lam*(j-i)//4-1)**2-orig[i] for i,j in edges]
    h=audit(orig,edges,new)
    check(h<=50176*n**4,'dyadic height')
    return h

def triadic(n,edges):
    L=3
    while L<=4*n+2: L*=9
    q=[3*i+1 for i in range(n)]
    orig=[4*L*z for z in q]
    new=[(L*(q[j]-q[i])-1)**2-orig[i] for i,j in edges]
    h=audit(orig,edges,new)
    check(h<=26244*n**4,'triadic height')
    return h

def square_scale(n,edges):
    lam=1
    while lam<16*n: lam*=4
    orig=[lam*(4*i+1) for i in range(n)]
    new=[(lam*(j-i)-1)**2-orig[i] for i,j in edges]
    h=audit(orig,edges,new)
    check(h<=4096*n**4,'square-scale height')
    return h

def seed_scale(seed,edges):
    H=max(seed)
    lam=1
    while lam<32*(H+1): lam*=4
    orig=[4*lam*x for x in seed]
    new=[(lam*(seed[j]-seed[i])-1)**2-orig[i] for i,j in edges]
    h=audit(orig,edges,new,True)
    for i,j in combinations(range(len(seed)),2):
        check(square(seed[i]+seed[j])==square(orig[i]+orig[j]),'seed changed')
    check(h<=65536*H**4,'general seed height')
    return h

families={'dyadic_mod7':dyadic,'odd_power_three':triadic,'square_scale':square_scale}
counts={}
for n in range(1,6):
    pairs=list(combinations(range(n),2))
    for mask in range(1<<len(pairs)):
        edges=[e for z,e in enumerate(pairs) if mask>>z&1]
        for f in families.values(): f(n,edges)
    counts[n]=1<<len(pairs)
rows=[]
for n in range(6,33):
    edges=list(combinations(range(n),2))
    rows.append({'n':n,'m':len(edges),'heights':{key:f(n,edges) for key,f in families.items()}})
rng=random.Random(20261019)
for z in range(100):
    n=rng.randrange(2,81)
    edges=list(combinations(range(n),2));rng.shuffle(edges)
    edges=edges[:rng.randrange(1,min(len(edges),150)+1)]
    for f in families.values():f(n,edges)
seed_counts=0
for H in range(1,31):
    seed=list(range(1,H+1))
    seed_scale(seed,list(combinations(range(H),2)))
    seed_counts+=1
for z in range(100):
    seed=sorted(rng.sample(range(1,10001),rng.randrange(2,31)))
    edges=list(combinations(range(len(seed)),2));rng.shuffle(edges)
    seed_scale(seed,edges[:min(len(edges),100)])
    seed_counts+=1
# Exact all-residue test for the claimed smallest prime and unique p=7 pattern.
residue_rows=[]
for p in [2,3,5,7]:
    squares={a*a%p for a in range(p)}
    admissible=[a for a in range(p) if 2*a%p not in squares]
    triples=[list(t) for t in combinations(admissible,3) if all((a+b)%p in squares for a,b in combinations(t,2))]
    residue_rows.append({'prime':p,'squares':sorted(squares),'admissible':admissible,'triples':triples})
check([r['triples'] for r in residue_rows]==[[],[],[],[[3,5,6]]],'minimal residue pattern')
print(json.dumps({'status':'PASS','always_active_checks':CHECKS,'exhaustive_graphs':counts,'complete_graphs':rows,'random_graphs':100,'scaled_seed_cases':seed_counts,'minimal_prime_pattern':residue_rows},indent=2,sort_keys=True))
