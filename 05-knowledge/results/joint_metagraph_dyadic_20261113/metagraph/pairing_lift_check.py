#!/usr/bin/env python3
"""Small exact controls for the converse-orbit lift polynomial; no census import."""
from collections import Counter, defaultdict
from itertools import permutations, product
from pathlib import Path
import hashlib, json

checks=0
def need(ok,msg):
    global checks
    checks+=1
    if not ok:raise RuntimeError(msg)

def lift_polynomial(sizes,G):
    P=Counter()
    for g in G:
        if any(g[g[i]]!=i for i in range(len(g))):continue
        fs=sum(sizes[i]==1 and g[i]==i for i in range(len(g)))
        fc=sum(sizes[i]==2 and g[i]==i for i in range(len(g)))
        pc=sum(sizes[i]==2 and g[i]>i for i in range(len(g)))
        Q=Counter({fs:1<<pc})
        for _ in range(fc):
            Z=Counter(Q)
            for k,v in Q.items():Z[k+2]+=v
            Q=Z
        P+=Q
    return P

def group_from_labels(labels):
    buckets=defaultdict(list)
    for i,key in enumerate(labels):buckets[key].append(i)
    result=[]
    for images in product(*(permutations(ids) for ids in buckets.values())):
        g=list(range(len(labels)))
        for ids,im in zip(buckets.values(),images):
            for i,j in zip(ids,im):g[i]=j
        result.append(tuple(g))
    return result

def brute_polynomial(sizes,G):
    X=[(i,b) for i,size in enumerate(sizes) for b in range(size)]
    index={v:i for i,v in enumerate(X)}
    R=[index[i,1-b] if sizes[i]==2 else j for j,(i,b) in enumerate(X)]
    G=set(G);P=Counter()
    # All ordinary permutations, independent of the proposed bit-lift formula.
    for f in permutations(range(len(X))):
        if any(f[f[i]]!=i or f[R[i]]!=R[f[i]] for i in range(len(X))):continue
        g=tuple(X[f[index[i,0]]][0] for i in range(len(sizes)))
        if g not in G:continue
        P[sum(f[i]==i for i in range(len(X)))]+=1
    return P

def main():
    generic=[]
    for s,t in [(s,t) for s in range(3) for t in range(3)]+[(3,0),(0,3)]:
        sizes=[1]*s+[2]*t;G=group_from_labels(sizes)
        P=lift_polynomial(sizes,G);B=brute_polynomial(sizes,G)
        need(P==B,'complete small centralizer involution polynomial')
        generic.append({'singletons':s,'double_fibres':t,'group_order':len(G),'polynomial':dict(sorted(P.items()))})
    # A proper subgroup which ties a singleton transposition to a double-fibre transposition.
    sizes=[1,1,2,2];G=[(0,1,2,3),(1,0,3,2)]
    P=lift_polynomial(sizes,G)
    need(P==brute_polynomial(sizes,G)==Counter({0:2,2:1,4:2,6:1}),'restricted diagonal subgroup')
    need(lift_polynomial([2,2],[(0,1),(1,0)])==Counter({0:3,2:2,4:1}),'double two-cycle has two lifts')
    root=Path(__file__).resolve().parents[1]
    receipt=root/'audit/small.normal.json';data=json.loads(receipt.read_text())
    rows=[]
    for row in data['rows'][:3]:
        meta=row['node_metadata'];codes=[x['code'] for x in meta];ix={c:i for i,c in enumerate(codes)}
        sizes=[1 if x['self_converse'] else 2 for x in meta];H=[x['H'] for x in meta]
        G=[tuple(ix[c] for c in g) for g in row['weighted_vertex_automorphisms']]
        noH=lift_polynomial(sizes,G)
        GH=[g for g in G if all(H[g[i]]==H[i] for i in range(len(g)))]
        withH=lift_polynomial(sizes,GH)
        expected={int(k):v for k,v in row['pairing_guard']['wired_fixed_histogram'].items()}
        need(withH==expected,'weighted graph H-guard comparison')
        abstract=lift_polynomial(sizes,group_from_labels(list(zip(sizes,H))))
        need(abstract=={int(k):v for k,v in row['pairing_guard']['abstract_fixed_histogram'].items()},'abstract H-fibre comparison')
        need(noH==[{0:1,2:1},{2:1,4:1},{8:1,10:2,12:1}][row['n']-3],'unrestricted graph minima')
        rows.append({'n':row['n'],'no_H_guard':dict(sorted(noH.items())), 'with_H_guard':dict(sorted(withH.items())),
                     'abstract_H_preserving':dict(sorted(abstract.items()))})
    print(json.dumps({'status':'PASS','always_active_checks':checks,'synthetic_complete_permutation_controls':generic,
                      'restricted_subgroup_polynomial':dict(sorted(P.items())),'metagraph':rows,
                      'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                      'root_receipt_sha256':hashlib.sha256(receipt.read_bytes()).hexdigest()},sort_keys=True,indent=2))
if __name__=='__main__':main()
