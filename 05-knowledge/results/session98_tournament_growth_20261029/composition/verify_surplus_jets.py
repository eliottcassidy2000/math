#!/usr/bin/env python3
"""Finite-surplus cycle-index jets: exact identities and independent palette orbits."""
from fractions import Fraction as F
from itertools import product
from collections import defaultdict
from math import comb
import json
from inherited_recreate import prime_index, recompose, scalar_product, partitions

def require(ok,msg):
    if not ok: raise RuntimeError(msg)

def jet(index,h):
    out=defaultdict(F)
    for part,c in index.items():
        term={():c}
        for d in part:
            new=defaultdict(F,term)
            for q,v in term.items():
                if sum(q)+d<=h: new[tuple(sorted(q+(d,)))]+=v
            term=new
        for q,v in term.items():out[q]+=v
    return dict(out)

def evaluate(poly,child,h,subtract_one):
    out=[F(0)]*(h+1)
    for part,c in poly.items():
        term=[c]+[F(0)]*h
        for d in part:
            v=[F(0)]*(h+1)
            for j in range(h//d+1):
                if j<len(child):v[d*j]=child[j]
            if subtract_one:v[0]-=1
            term=scalar_product(term,v,h)
        out=[a+b for a,b in zip(out,term)]
    return out

def main():
    bound=15;P=prime_index(bound);a=recompose(P,bound)
    checks=0;table=[]
    for k in range(3,bound+1):
        pk={p:v for p,v in P.items() if sum(p)==k}
        h=min(8,bound-k)
        j=jet(pk,h)
        actual=evaluate(pk,a[1:],h,False)
        require(actual==evaluate(j,a[1:],h,True),'jet mismatch')
        for r in range(h+1):
            require(jet(pk,r)=={p:v for p,v in j.items() if sum(p)<=r},'projection mismatch')
            checks+=1
        require(all(v.denominator==1 for v in actual),'fractional inflation count')
        table.append({'arity':k,'surplus_coefficients':list(map(int,actual)),
                      'jet_coordinates':len(j)})
    # Reconstruct a_n using only the triangular surplus-jet recurrence.
    recovered=[1]+[0]*bound; strong=[0]*(bound+1)
    for n in range(1,bound+1):
        strong[n]=F(n==1)
        for k in range(3,n+1):
            pk={p:v for p,v in P.items() if sum(p)==k}
            strong[n]+=evaluate(jet(pk,n-k),recovered[1:],n-k,True)[n-k]
        recovered[n]=sum(strong[j]*recovered[n-j] for j in range(1,n+1))
    require(recovered==a,'triangular jet recurrence mismatch')
    # Direct cyclic-action orbits on actual rigid TT1,TT2,TT3 palette.
    orbit_checks=[]
    for k in (3,5,7):
        index={(1,)*k:F(1,k),(k,):F(k-1,k)}
        h=2*k
        pred=evaluate(jet(index,h),[1,1,1],h,True)
        seen=set();counts=[0]*(h+1)
        for word in product(range(3),repeat=k):
            canonical=min(word[i:]+word[:i] for i in range(k))
            if canonical not in seen:
                seen.add(canonical);counts[sum(word)]+=1
        require(pred==counts,'direct palette orbit mismatch')
        orbit_checks.append({'root_action':'C'+str(k),'counts':counts})
    # Missing z3 is genuinely wrong at surplus 3 for C3.
    c3={(1,1,1):F(1,3),(3,):F(2,3)}
    j3=jet(c3,3); full=evaluate(j3,a[1:],3,True)
    without=evaluate({p:v for p,v in j3.items() if 3 not in p},a[1:],3,True)
    require(full[3]-without[3]==F(2,3),'cycle onset hostile failed')
    require(full==[1,1,3,9],'C3 formula')
    print(json.dumps({'status':'PASS','bank_bound':bound,'horizon_checks':checks,'triangular_reconstruction':list(map(int,recovered)),
      'dimension_by_h':[sum(1 for n in range(h+1) for p in partitions(n)) for h in range(9)],
      'prime_inflations':table,'independent_palette_orbits':orbit_checks,
      'c3_missing_three_cycle_error':str(full[3]-without[3])},indent=2))
if __name__=='__main__':main()
