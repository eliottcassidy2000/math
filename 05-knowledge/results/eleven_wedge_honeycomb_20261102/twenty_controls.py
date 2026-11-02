#!/usr/bin/env python3
"""Exact twenty-object shadow, including reversal and a failed rotation."""
from itertools import combinations
from math import isqrt
from collections import Counter
import json
def need(b,s):
    if not b:raise RuntimeError(s)
def run():
    n=25;adj={v:[] for v in range(1,n+1)}
    for a,b in combinations(adj,2):
        if isqrt(a+b)**2==a+b:adj[a].append(b);adj[b].append(a)
    need([v for v in adj if len(adj[v])==1]==[18],'forced endpoint')
    colors={v:len(adj[v]) for v in adj};refinement=[]
    for _ in range(25):
        signatures={v:(colors[v],tuple(sorted(colors[w] for w in adj[v]))) for v in adj}
        classes=sorted(set(signatures.values()));new={v:classes.index(signatures[v]) for v in adj}
        refinement.append(len(classes))
        if len(set(colors.values()))==len(classes):break
        colors=new
    need(len(set(colors.values()))==25,'rigid guard via invariant refinement')
    paths=[]
    def dfs(p,seen):
        if len(p)==25:paths.append(tuple(p));return
        for v in adj[p[-1]]:
            if v not in seen:dfs(p+(v,),seen|{v})
    dfs((18,),{18});paths=sorted(paths+[p[::-1] for p in paths])
    need(len(paths)==20,'path universe')
    arcs=[set(zip(p,p[1:])) for p in paths];R={i:paths.index(p[::-1]) for i,p in enumerate(paths)};B={}
    for i,j in combinations(range(20),2):
        if not any((b,a) in arcs[j] for a,b in arcs[i]):
            need(i not in B and j not in B,'compatibility not matching');B[i]=j;B[j]=i
    E=list(combinations(range(4),2));Y=list(combinations(range(6),3));full=set(range(6));pi=[1,0,3,2]
    def complement(y):return tuple(sorted(full-set(y)))
    def perm(y):return tuple(sorted(E.index(tuple(sorted((pi[E[e][0]],pi[E[e][1]])))) for e in y))
    def degree(y):return sorted(sum(v in E[e] for e in y) for v in range(4))
    RY={y:complement(perm(y)) for y in Y}
    BY={y:complement(y) for y in Y if degree(y) in [[0,2,2,2],[1,1,1,3]]}
    need(len(B)==len(BY)==8 and all(RY[y]!=y and RY[RY[y]]==y for y in Y),'involutions')
    need(all(R[B[i]]==B[R[i]] for i in B),'commuting source involutions')
    f={}
    for active in [True,False]:
        for i in range(20):
            if i in f or (i in B)!=active:continue
            y=min(y for y in Y if y not in f.values() and (y in BY)==active)
            entries={i:y,R[i]:RY[y]}
            if active:entries.update({B[i]:BY[y],R[B[i]]:RY[BY[y]]})
            need(len(entries)==(4 if active else 2),'orbit shape')
            f.update(entries)
    need(len(f)==len(set(f.values()))==20,'bijection')
    for i in range(20):
        need(f[R[i]]==RY[f[i]],'reversal transport')
        need((i in B)==(f[i] in BY),'compatibility domain')
        if i in B:need(f[B[i]]==BY[f[i]],'compatibility transport')
    rows=[{'index':i,'path':paths[i],'K4_edge_indices':f[i],'K4_edges':[E[e] for e in f[i]],'reverse':R[i],'compatible_partner':B.get(i)} for i in range(20)]
    # Arbitrary rotation of six edge slots preserves weight but not this relation.
    y=(0,1,2);rot=tuple(sorted((e+1)%6 for e in y))
    need(y in BY and rot not in BY,'rotation hostile')
    unseen=set(Y);necklaces=[]
    while unseen:
        a=min(unseen);orb={tuple(sorted((e+k)%6 for e in a)) for k in range(6)};unseen-=orb
        necklaces.append({'representative':a,'orbit_size':len(orb)})
    need(sorted(r['orbit_size'] for r in necklaces)==[2,6,6,6],'six-bit necklace spectrum')
    guard_edges=sum(map(len,adj.values()))//2
    need(guard_edges==32 and all(len(a)==24 for a in arcs),'guard edge clock')
    pair_sizes=[len(arcs[i]|arcs[j]) for i,j in B.items() if i<j]
    need(pair_sizes==[26]*4,'compatible intersection clocks')
    # Compatibility is a matching, so no three path cubes meet. Expand the
    # exact factorial-moment identity, rather than scanning 2^32 orientations.
    first=sum(2**(guard_edges-len(a)) for a in arcs)
    second=sum(2**(guard_edges-s) for s in pair_sizes)
    polynomial={'0':2**guard_edges-first+second,'1':first-2*second,'2':second,
                'outside_free_bits':n*(n-1)//2-guard_edges}
    need(polynomial=={'0':4294962432,'1':4608,'2':256,'outside_free_bits':268},'weighted path polynomial')
    return {'status':'passed','scope':'Finite relation isomorphism; not an isomorphism of integer-labelled source graphs or arbitrary native moves.',
            'guard_color_refinement_class_counts':refinement,
            'mapping':rows,'K4_shape_counts':dict(Counter(str(degree(y)) for y in Y)),
            'reversal_on_reference':'edge complement after vertex permutation (01)(23)',
            'six_slot_rotation_hostile':{'source':y,'source_degrees':degree(y),'target':rot,'target_degrees':degree(rot)},
            'weight_three_necklaces':necklaces,'guarded_polynomial':polynomial}
if __name__=='__main__':print(json.dumps(run(),indent=2,sort_keys=True))
