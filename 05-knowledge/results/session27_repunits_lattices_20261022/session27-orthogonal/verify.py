#!/usr/bin/env python3
"""Independent small split-orthogonal enumeration; Python standard library only."""
from itertools import product, combinations
from collections import Counter
from fractions import Fraction
import json, random
from pathlib import Path
q=3; n=4; r=2

def basis(rows):
    a=[list(v) for v in rows if any(v)]
    k=0
    for j in range(n):
        p=next((i for i in range(k,len(a)) if a[i][j]%q),None)
        if p is None: continue
        a[k],a[p]=a[p],a[k]
        z=pow(a[k][j]%q,-1,q); a[k]=[(x*z)%q for x in a[k]]
        for i in range(len(a)):
            if i!=k:
                z=a[i][j]%q; a[i]=[(x-z*y)%q for x,y in zip(a[i],a[k])]
        k+=1
        if k==len(a): break
    return tuple(tuple(v) for v in a[:k])

def Q(v): return (v[0]*v[2]+v[1]*v[3])%q

def bil(v,w): return (v[0]*w[2]+v[1]*w[3]+v[2]*w[0]+v[3]*w[1])%q

def isotropic(w): return all(Q(v)==0 for v in w) and all(bil(v,u)==0 for v,u in combinations(w,2))

def subset(u,w): return basis(u+w)==w

lines=sorted({basis([v])[0] for v in product(range(q),repeat=n) if any(v)})
planes=sorted({basis([v,w]) for v,w in combinations(lines,2)})
assert all(len(p)==2 for p in planes)
L=[w for w in planes if isotropic(w)]
B=[v for v in lines if Q(v)==0]
assert len(lines)==40 and len(planes)==130 and len(B)==16 and len(L)==8
masks=[sum(1<<i for i,v in enumerate(B) if subset((v,),w)) for w in L]
assert all(m.bit_count()==4 for m in masks)

# Independently enumerate all spans of subsets of the 16 singular lines.
spans=[()]*(1<<len(B))
weighted=[0]*len(spans)
for s in range(len(spans)):
    if s:
        bit=s&-s; i=bit.bit_length()-1
        spans[s]=basis(spans[s^bit]+(B[i],))
    w=spans[s]
    # Formula prediction used only on this side, not to construct L.
    denom=1
    if isotropic(w):
        for j in range(r-len(w),r): denom*=q**j+1
        predicted=Fraction(len(L),denom)
        assert predicted.denominator==1
        predicted=int(predicted)
    else: predicted=0
    actual=sum(subset(w,l) for l in L)
    assert actual==predicted, (s,w,actual,predicted)
    weighted[s]=(-1)**s.bit_count()*predicted

# Subset zeta transform gives IE prediction for every forbidden subset.
avoid=weighted[:]
for i in range(len(B)):
    for s in range(len(avoid)):
        if s&(1<<i): avoid[s]+=avoid[s^(1<<i)]
for s,predicted in enumerate(avoid):
    actual=sum((s&mask)==0 for mask in masks)
    assert actual==predicted,(s,actual,predicted)

# Separately construct actual span-posets and compute mu recursively.
def mobius_check(s):
    sub=[]; t=s
    while True:
        sub.append(t)
        if not t: break
        t=(t-1)&s
    poset=sorted({spans[t] for t in sub},key=lambda w:(len(w),w))
    mu={}
    for w in poset:
        mu[w]=1 if not w else -sum(mu[u] for u in mu if subset(u,w))
    grouped=Counter()
    for t in sub: grouped[spans[t]]+=(-1)**t.bit_count()
    assert all(grouped[w]==mu[w] for w in poset)
    val=Fraction(0)
    for w in poset:
        if isotropic(w):
            denom=1
            for j in range(r-len(w),r): denom*=q**j+1
            val+=Fraction(mu[w],denom)
    assert val==Fraction(avoid[s],len(L))
    return len(poset)

rng=random.Random(20261022)
cases=[0,(1<<len(B))-1]+[rng.randrange(1<<len(B)) for _ in range(100)]
for s in cases: mobius_check(s)
# Three distinct projective lines in an isotropic plane.
plane=L[0]
indices=[i for i,v in enumerate(B) if subset((v,),plane)]
hostile=sum(1<<i for i in indices[:3])
assert len(spans[hostile])==2
assert Fraction(avoid[hostile],8)==1-3*Fraction(1,4)+2*Fraction(1,8)==Fraction(1,2)
result={
 'field':q,'ambient_dimension':n,'all_projective_lines':len(lines),
 'all_two_subspaces':len(planes),'isotropic_lines':len(B),'maximal_isotropic_planes':len(L),
 'isotropic_line_incidence_degrees':sorted({sum(bool(mask&(1<<i)) for mask in masks) for i in range(len(B))}),
 'forbidden_subsets_exhaustively_checked':len(spans),
 'span_containment_formula_subsets_checked':len(spans),
 'mobius_posets_checked':len(cases),
 'avoidance_count_histogram':dict(sorted(Counter(avoid).items())),
 'hostile_three_lines': [B[i] for i in indices[:3]],
 'hostile_span_rank':2,'hostile_avoiding_planes':avoid[hostile],
 'hostile_avoidance_probability':str(Fraction(avoid[hostile],8)),
 'complete_forbidden_set_avoidance':avoid[-1],
 'status':'all assertions passed'
}
print(json.dumps(result,indent=2))
Path(__file__).with_name('results.json').write_text(json.dumps(result,indent=2)+'\n')
