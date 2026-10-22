#!/usr/bin/env python3
"""Exact self-converse census and independent anti-permutation Burnside check."""
from collections import Counter
from itertools import permutations
from math import factorial, gcd
from fractions import Fraction
import json
from pathlib import Path

def require(condition):
    if not condition:
        raise RuntimeError("Exact census validation failed")

def oddparts(n, lo=1):
    if not n:
        yield ()
    for d in range(lo,n+1,2):
        for tail in oddparts(n-d,d):
            yield (d,)+tail

def census(m):
    fac=factorial(m); ev=od=0; buckets={}
    for p in oddparts(m):
        z=1
        for d,a in Counter(p).items(): z*=d**a*factorial(a)
        e=m+2*sum(gcd(x,y) for i,x in enumerate(p) for y in p[i+1:])
        weight=fac//z
        require(fac%z==0)
        ev+=weight*2**(e-len(p));od+=weight*2**e
        defect=m-len(p)
        buckets[defect]=buckets.get(defect,0)+weight*2**(e-len(p))
        k=sum(x for x in p if x>1);l=sum(x>1 for x in p)
        loss=m*(m-1)-(e-len(p))
        require(loss>=defect*(2*m-k-1))
        require(k<=3*defect//2)
    require(ev%fac==od%fac==0)
    return ev//fac,od//fac,buckets

def direct(n):
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    idx={p:i for i,p in enumerate(pairs)};total=0
    for perm in permutations(range(n)):
        nxt=[];flip=[]
        for i,j in pairs:
            a,b=perm[i],perm[j]
            nxt.append(idx[tuple(sorted((a,b)))]);flip.append(1^(a>b))
        seen=set();components=0;ok=True
        for e in range(len(pairs)):
            if e in seen:continue
            x=e;parity=0;components+=1
            while x not in seen:
                seen.add(x);parity^=flip[x];x=nxt[x]
            if parity:ok=False;break
        if ok:total+=2**components
    require(total%factorial(n)==0)
    return total//factorial(n)

def falling(m,k):return factorial(m)//factorial(m-k) if m>=k else 0

def power2(k):return Fraction(2**k) if k>=0 else Fraction(1,2**(-k))

def main():
    known=[1,1,2,2,8,12,88,176,2752,8784,279968,1492288,95458560,872687552,111698291584,1787154671104,457509297625088,13013584213369088,6662951988432581120,341143107490935724032,349330527429800077778944,32519496073514216703585280,66599927718689841715666229248]
    results={};details=[]
    for m in range(31):
        ev,od,b=census(m);results[2*m]=ev;results[2*m+1]=od
        if m>=6:
            base=2**(m*(m-1))
            a=Fraction(falling(m,3),3)*power2(8-4*m)
            c=(Fraction(falling(m,5),5)*2**24+Fraction(falling(m,6),18)*2**28)*power2(-8*m)
            require(Fraction(b.get(2,0),base)==a)
            require(Fraction(b.get(4,0),base)==c)
            remainder=Fraction(ev*factorial(m),base)-1-a-c
            require(remainder>=0)
            details.append({'m':m,'even_relative_excess':float(Fraction(ev*factorial(m),base)-1),'after_defect4_scaled_by_m9_2neg12m':float(remainder/power2(-12*m)/m**9)})
    require([results[n] for n in range(1,24)]==known)
    small={n:direct(n) for n in range(1,8)}
    require(all(results[n]==v for n,v in small.items()))
    out={'scope':'Exact partition census n=0..61; direct permutation anti-Burnside n=1..7; OEIS displayed values n=1..23; exponent inequalities each partition m<=30','known_match':True,'direct':small,'counts':results,'diagnostics_not_proof':details}
    Path(__file__).with_name('receipt.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'known_match':True,'direct':small,'computed_through_n':61,'last_diagnostic':details[-1]}))
if __name__=='__main__':main()
