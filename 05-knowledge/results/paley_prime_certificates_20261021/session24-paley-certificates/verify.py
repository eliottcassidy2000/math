"""Exact affine triple tests. Checks remain active under python -O."""
import itertools, json, hashlib
from pathlib import Path

def check(x, msg):
    if not x: raise RuntimeError(msg)
def primes(n):
    return [p for p in range(3,n+1) if all(p%d for d in range(2,int(p**.5)+1))]
def stab(S,p,full=False):
    S=set(S)
    # Independent paths: all p(p-1) maps, or six possible images of 0 and 1.
    maps=itertools.product(range(1,p),range(p)) if full else ((v-u,u) for u,v in itertools.permutations(S,2))
    return {(a%p,b) for a,b in maps if {(a*x+b)%p for x in S}==S}
rows=[]
for p in primes(31):
    if p==3: continue
    counts={1:0,2:0,3:0}
    for S in itertools.combinations(range(p),3):
        H=stab(S,p,True)
        check(len(H) in counts,('stabilizer order',p,S,H))
        counts[len(H)]+=1
    I=int(p%3==1)
    expected={1:p*(p-1)*(p-5-2*I)//6,2:p*(p-1)//2,3:I*p*(p-1)//3}
    check(counts==expected,('counts',p,counts,expected))
    rows.append({'p':p,'stabilizer_histogram':counts,'free_orbits':counts[1]//(p*(p-1))})
norm_cases=0
for p in primes(1009):
    if p==3: continue
    bad={p-1,2,pow(2,-1,p)}|{t for t in range(p) if (t*t-t+1)%p==0}
    for t in range(2,p):
        H=stab({0,1,t},p)
        check((len(H)>1)==(t in bad),('normalization',p,t,H))
        norm_cases+=1
    check((len(stab({0,1,3},p))==1)==(p>=11),('fixed triple',p))
exception_rows=[]
for p in (3,7):
    maps=list(itertools.product(range(1,p),range(p)))
    regular=0
    for mask in range(1<<p):
        S={x for x in range(p) if mask>>x&1}
        regular+=sum({(a*x+b)%p for x in S}==S for a,b in maps)==1
    check(regular==0,('exception',p,regular))
    # Exhaustive permutation check independently confirms full signed group.
    Q={x*x%p for x in range(1,p)}
    signed=[]
    for perm in itertools.permutations(range(p)):
        signs={((j-i)%p in Q)==((perm[j]-perm[i])%p in Q) for i in range(p) for j in range(i+1,p)}
        if len(signs)==1: signed.append(perm)
    affine={tuple((a*x+b)%p for x in range(p)) for a,b in maps}
    check(set(signed)==affine,('full signed group',p))
    exception_rows.append({'p':p,'signed_group_order':len(signed),'regular_subsets':regular})
print(json.dumps({'status':'PASS','exhaustive_affine_triples':rows,'normalized_primes_bound':1009,'normalized_cases':norm_cases,'small_Paley_full_permutation_controls':exception_rows},sort_keys=True,indent=2))
