#!/usr/bin/env python3
"""Independent exact rational and signed-orbit audit; no producer code imported."""
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations, permutations
from math import factorial, gcd, prod
from pathlib import Path
import hashlib, json
HERE=Path(__file__).resolve().parent
SOURCE=Path('/workspace/shared/session29-selfconverse-asymptotics')
def need(x, message):
    if not x: raise RuntimeError(message)
def parts(n, minimum=1):
    if n==0: yield ()
    for i in range(minimum,n+1):
        for t in parts(n-i,i): yield (i,)+t
def z(p): return prod(a**b*factorial(b) for a,b in Counter(p).items())
def representative(p):
    out=[]; start=0
    for r in p:
        out.extend(list(range(start+1,start+r))+[start]);start+=r
    return tuple(out)
def signed_fix(p):
    """Follow ordered endpoints until the unordered pair returns; impose swap iff odd."""
    unseen=set(combinations(range(len(p)),2)); c=0
    while unseen:
        i,j=next(iter(unseen)); a,b=i,j; t=0
        while True:
            unseen.remove(tuple(sorted((a,b))))
            a,b=p[a],p[b];t+=1
            if {a,b}=={i,j}: break
        if bool(t%2)!=(a==j): return 0
        c+=1
    return 1<<c
def admissible(p):
    return p.count(1)<=1 and all(x==1 or x%4==2 for x in p)
def falling(m,r): return prod(range(m-r+1,m+1)) if r<=m else 0
def pow2(e): return Q(2**e) if e>=0 else Q(1,2**(-e))
def oddparts(n, minimum=1):
    if n==0: yield ()
    for i in range(minimum,n+1,2):
        for t in oddparts(n-i,i):yield (i,)+t
def formula(m):
    even=odd=Q(0); buckets={}; counts={}; checked=0
    for p in oddparts(m):
        e=sum(p)+2*sum(gcd(a,b) for a,b in combinations(p,2))
        ev=Q(1<<(e-len(p)),z(p)); od=Q(1<<e,z(p))
        even+=ev;odd+=od;d=m-len(p)
        relative=ev*factorial(m)*pow2(-m*(m-1))
        buckets[d]=buckets.get(d,0)+relative
        counts[d]=counts.get(d,0)+1
        non=[a for a in p if a>1]; k=sum(non); ell=len(non)
        L=m*(m-1)-(e-len(p)); S=sum(gcd(a,b) for a,b in combinations(non,2))
        need(L==2*m*d-(k*k-2*k*ell+2*k-ell+2*S),'exact loss')
        need(2*S<=(ell-1)*k,'gcd sum')
        need(k<=3*d//2 and L>=d*(2*m-k-1) and L>=d*(m-1),'loss bound')
        need(Q(factorial(m),z(p))<=m**k,'denominator bound')
        need(od*factorial(m)*pow2(-m*m)==relative*pow2(-d),'odd transfer')
        checked+=1
    for d,count in counts.items():need(count<=2**(3*d//2+1),'partition count bound')
    A=Q(falling(m,3),3)*pow2(8-4*m)
    C=(Q(falling(m,5),5)*2**24+Q(falling(m,6),18)*2**28)*pow2(-8*m)
    D=(Q(falling(m,7),7)*2**48+Q(falling(m,8),15)*2**48+Q(falling(m,9),162)*2**60)*pow2(-12*m)
    need(buckets.get(2,0)==A,'A coefficient');need(buckets.get(4,0)==C,'C coefficient')
    need(buckets.get(6,0)==D,'defect-six coefficient')
    need(even.denominator==odd.denominator==1,'integral census')
    need(even*factorial(m)*pow2(-m*(m-1))>=1,'mean below factorial')
    return int(even),int(odd),checked

def main():
    results={}; partition_checks=0
    for m in range(36):
        a,b,c=formula(m);results[2*m]=a;results[2*m+1]=b;partition_checks+=c
    prior=json.loads((SOURCE/'receipt.json').read_text())['counts']
    need(all(results[int(n)]==v for n,v in prior.items()),'producer census mismatch')
    typecounts={}; types_checked=0
    for n in range(31):
        total=Q(0)
        for p in parts(n):
            f=signed_fix(representative(p));types_checked+=1
            need(bool(f)==admissible(p),'cycle admissibility')
            if f:
                orbit_exp=sum(x//2 for x in p)+sum(gcd(a,b) for a,b in combinations(p,2))
                need(f==1<<orbit_exp,'pair-orbit exponent')
            total+=Q(f,z(p))
        need(total==results[n],'unrestricted partition census')
        typecounts[n]=int(total)
    direct={}
    for n in range(9):
        raw=sum(signed_fix(p) for p in permutations(range(n)))
        need(raw%factorial(n)==0,'direct integrality')
        direct[n]=raw//factorial(n)
        need(direct[n]==results[n],'all permutation census')
    sources={str(SOURCE/f):hashlib.sha256((SOURCE/f).read_bytes()).hexdigest() for f in ['REPORT.md','check.py','receipt.json']}
    out={'status':'PASS','all_permutations_n0_through8':direct,'unrestricted_cycle_types_n0_through30':typecounts,'cycle_types_checked':types_checked,'odd_partitions_checked_m0_through35':partition_checks,'producer_n0_through61_match':True,'census_n0_through71':results,'source_sha256':sources,'controls':'No producer code imported. Fraction accumulation. Explicit RuntimeError guards survive -O. Edge-orbit endpoint-swap test independent of admissible-cycle filter.'}
    print(json.dumps(out,indent=2,sort_keys=True))
if __name__=='__main__':main()
