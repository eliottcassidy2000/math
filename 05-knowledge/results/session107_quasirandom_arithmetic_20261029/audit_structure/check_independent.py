#!/usr/bin/env python3
"""Independent finite guards; no statistical frequency claimed as proof."""
from itertools import permutations, product, combinations
from collections import Counter
from math import gcd,prod
import json
from pathlib import Path

def check(x,msg):
    if not x: raise RuntimeError(msg)
def prime(n): return n>=2 and all(n%d for d in range(2,int(n**.5)+1))
def chi(a,p):
    r=pow(a%p,(p-1)//2,p)
    return -1 if r==p-1 else r

def jacobi(a,n):
    out=1;a%=n
    while a:
        while a%2==0:
            a//=2
            if n%8 in (3,5):out=-out
        a,n=n,a
        if a%4==n%4==3:out=-out
        a%=n
    return out if n==1 else 0
# Odd composite q may be used here: reciprocity evaluates the Dirichlet character.
base_sets=[[-1,2,-2],[2,3,6],[4,2,5,3,7,13],[2,3,7,5,2],[-12,18,8,-1,1,4]]
crt=[]
for bases in base_sets:
    odd=[p for p in range(3,max(map(abs,bases))+1,2) if prime(p) and any(a%p==0 for a in bases)]
    Q=8*prod(odd)
    counts=Counter(tuple(jacobi(a,q) for a in bases) for q in range(3,Q,4) if gcd(q,Q)==1)
    check(len(set(counts.values()))==1,'CRT signatures not uniform')
    crt.append({'bases':bases,'modulus':Q,'attainable':len(counts),'classes_per_signature':next(iter(counts.values()))})
check([r['attainable'] for r in crt]==[2,4,32,16,4],'squareclass ranks')
# Quadratic reciprocity orientation and exclusions.
ps=[p for p in range(3,500) if prime(p)]
for p,q in combinations(ps,2):
    check(chi(p,q)*chi(q,p)==(-1 if p%4==q%4==3 else 1),'orientation')
for p in ps:
    check(chi(4,p)==1,'square port')
    check(chi(-1,p)==(-1 if p%4==3 else 1),'minus-one gauge')
    check(chi(2,p)==(-1 if p%8 in (3,5) else 1),'two supplemental law')
# All maps of Paley7; signed and ordinary action are distinct.
Q={1,2,4}
def edge(x,y):return (y-x)%7 in Q
ordinary=[];reverse=[]
for f in permutations(range(7)):
    same=all(edge(f[x],f[y])==edge(x,y) for x,y in combinations(range(7),2))
    opposite=all(edge(f[x],f[y])!=edge(x,y) for x,y in combinations(range(7),2))
    if same: ordinary.append(f)
    if opposite:reverse.append(f)
check(set(ordinary)=={tuple((a*x+b)%7 for x in range(7)) for a in Q for b in range(7)},'Paley ordinary group')
check(len(reverse)==21,'Paley converse group')
check(not any(all(f[f[x]]==x for x in range(7)) and any(f[x]!=x for x in range(7)) for f in ordinary),'involution obstruction')
# Monster order valuation map and corrected 2.B centralizer normalization.
M={2:46,3:20,5:9,7:6,11:2,13:3,17:1,19:1,23:1,29:1,31:1,41:1,47:1,59:1,71:1}
B={2:41,3:13,5:6,7:2,11:1,13:1,17:1,19:1,23:1,31:1,47:1}
delta={p:e-B.get(p,0) for p,e in M.items() if e>B.get(p,0)}
check([delta[p] for p in [13,2,5,3,7]]==[2,5,3,7,4],'valuation chain')
check(prod(p**e for p,e in delta.items())==194478922284018372000,'order ratio')
check(1+4371+96255+96256==196883,'character degrees')
check(47*59*71==196883 and 3*31*47==4371,'degree factors')
# Sylvester identity and character directions for all prime terms in first 7 terms.
s=[2]
for _ in range(6):s.append(s[-1]**2-s[-1]+1)
for i in range(1,len(s)):check(s[i]==prod(s[:i])+1,'product recurrence')
prime_terms=[x for x in s[:6] if prime(x) and x%4==3]
for p,q in combinations(prime_terms,2):check(chi(q,p)==1 and chi(p,q)==-1,'transitive prime terms')
for p in [3,7,43]:check(chi(13,p)==chi(139,p) and chi(p,13)==-chi(p,139),'factor-sidecar')
# Wieferich lifting does not determine a shared quadratic sign.
for p,expected in [(1093,-1),(3511,1)]:check(pow(2,p-1,p*p)==1 and chi(2,p)==expected,'Wieferich guard')
check(s[5]==3263443 and prime(s[5]),'later prime term')
fac=[547,607,1033,31051]
check(prod(fac)==s[6] and all(prime(p) for p in fac),'later factorization')
for p in [3,7,43,13,139,3263443]:check(prod(chi(q,p) for q in fac)==1,'later zero syndrome')
fermat=[2**(2**j)+1 for j in range(6)]
for j in range(1,6):check(fermat[j]==2+prod(fermat[:j]),'Fermat product')
check(641*6700417==fermat[5] and prime(641) and prime(6700417),'F5 factorization')
for p in fermat[:5]:check(chi(641,p)*chi(6700417,p)==chi(2,p),'Fermat syndrome')
for p,q in combinations(fermat[2:5],2):check(chi(p,q)==chi(q,p)==1,'Fermat residue clique')
result={'status':'PASS','CRT':crt,'odd_primes_checked':len(ps),'Paley_ordinary':len(ordinary),'Paley_reversers':len(reverse),'Sylvester_prime_terms_tested':prime_terms,'proof_scope':'Finite guards only; no infinite Sylvester-prime assertion.'}
print(json.dumps(result,indent=2))
