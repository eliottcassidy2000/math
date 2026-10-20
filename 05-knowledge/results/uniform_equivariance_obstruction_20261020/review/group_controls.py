#!/usr/bin/env python3
"""Independent local group controls, without subgroup census."""
from itertools import combinations
from independent_controls import graph_data, cut_space
import json

def compose(a,b): return tuple(a[b[i]] for i in range(len(a)))
def inverse(a):
    b=[0]*len(a)
    for i,x in enumerate(a): b[x]=i
    return tuple(b)
def support(a): return {i for i,x in enumerate(a) if i!=x}
def double_transposition(a):
    return len(support(a))==4 and compose(a,a)==tuple(range(len(a)))

def cycles(n):
    out=[]
    for a,b,c in combinations(range(n),3):
        p=list(range(n)); p[a]=b;p[b]=c;p[c]=a
        out += [tuple(p), inverse(p)]
    return out

def local_checks(n):
    c=cycles(n); counts={1:0,2:0}
    for a in c:
        for b in c:
            overlap=len(support(a)&support(b))
            if overlap not in (1,2): continue
            counts[overlap]+=1
            if overlap==1:
                d=compose(compose(b,a),inverse(b))
                assert len(support(a)&support(d))==2
                b=d
            assert any(double_transposition(compose(a,z)) for z in [b,inverse(b)])
    return {'n':n,'ordered_overlap_one_pairs':counts[1],'ordered_overlap_two_pairs':counts[2]}

def source_checks(n):
    edges,ps,images,invs,affines,gens=graph_data(n)
    u=cut_space(n,edges); size=1<<len(edges)
    canon=[min(t^x for x in u) for t in range(size)]
    identity=tuple(range(n))
    involutions=[i for i,p in enumerate(ps) if p!=identity and compose(p,p)==identity]
    fixedpoint_involutions=[i for i in involutions if len(support(ps[i]))<n]
    dt=[i for i in involutions if double_transposition(ps[i])]
    assert all(row[t]!=t for i in involutions for row in [affines[i]] for t in range(size))
    assert all(canon[affines[i][t]]!=t for i in fixedpoint_involutions for t in set(canon))
    dt_classes=sum(any(canon[affines[i][t]]==t for i in dt) for t in set(canon))
    if n>=5: assert dt_classes==0
    return {'n':n,'nontrivial_involutions':len(involutions),'fixed_point_involutions':len(fixedpoint_involutions),'switching_classes_with_double_transposition':dt_classes}

print(json.dumps({'cycle_pair_controls':[local_checks(n) for n in range(3,8)],'source_involution_controls':[source_checks(n) for n in range(2,6)]},indent=2))
