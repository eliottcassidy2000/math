#!/usr/bin/env python3
from itertools import combinations_with_replacement, permutations
from collections import Counter
from pathlib import Path
import math,json

def norm(xs):
    return tuple(sorted(y for x in xs for y in ([2,2,3] if x==6 else [x])))
byvalue={};bylength={};checks=0;collisions=0
for k in range(7):
    for xs in combinations_with_replacement(range(2,21),k):
        v=math.prod((1<<x)-1 for x in xs);n=norm(xs)
        if v in byvalue:
            assert byvalue[v]==n,(v,xs,n,byvalue[v]);collisions+=1
        byvalue[v]=n
        key=(k,v)
        assert key not in bylength or bylength[key]==xs,(key,xs,bylength[key])
        bylength[key]=xs;checks+=1
power_counts=[]
for k in range(13):
    v=63**k
    candidates=[n for n in range(2,v.bit_length()+1) if v%((1<<n)-1)==0]
    assert candidates==([2,3,6] if k else [])
    sol=[]
    def rec(i,left,word):
        if i==len(candidates):
            if left==1:sol.append(tuple(word))
            return
        n=candidates[i];m=(1<<n)-1
        rec(i+1,left,word)
        while left%m==0:
            left//=m;word=word+[n];rec(i+1,left,word)
    rec(0,v,[])
    assert len(sol)==k+1
    power_counts.append(len(sol))
matching=[]
for r in range(1,9):
    good=[p for p in permutations(range(r)) if all(i==0 or p[i]==0 or p[i]==i for i in range(r))]
    assert len(good)==r
    matching.append(len(good))
# Boundary: allowing the unit generator destroys fixed-length uniqueness.
assert math.prod([(1<<6)-1,1,1])==math.prod([3,3,7])
receipt={'status':'PASS','normal_form_universe':{'indices':[2,20],'factor_counts':[0,6],'products_checked':checks,'repeated_values':collisions},'powers63_exponents':[0,12],'factorization_counts':power_counts,'star_diagonal_matching_sizes':matching,'unit_hostile':'M6*M1*M1 = M2*M2*M3 = 63, three factors on both sides'}
Path(__file__).with_name('mersenne_independent_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
