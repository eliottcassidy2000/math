#!/usr/bin/env python3
"""Deletion refinement and common-inflation controls, using the public H oracle."""
from collections import defaultdict
from itertools import combinations, product
import json
from check_lifting import tournament, triangles, paths

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def code(a):
    return sum(a[i][j]<<b for b,(i,j) in enumerate(combinations(range(len(a)),2)))

def restrict(a,vertices):return [[a[i][j] for j in vertices] for i in vertices]

def inflate(a,children):
    addresses=[(i,v) for i,t in enumerate(children) for v in range(len(t))]
    return [[children[i][v][w] if i==j else a[i][j] for j,w in addresses] for i,v in addresses]

def main():
    states=[];matrices=[];observations=[];lookup={}
    for n in range(6):
        for c in range(1<<(n*(n-1)//2)):
            lookup[n,c]=len(states);states.append((n,c));a=tournament(n,c);matrices.append(a);observations.append((n,paths(a)[0][-1]))
    def ids(signatures):
        mapping={};result=[]
        for sig in signatures:
            if sig not in mapping:mapping[sig]=len(mapping)
            result.append(mapping[sig])
        return result
    actions=[]
    for a in matrices:
        n=len(a);actions.append([lookup[n-1,code(restrict(a,[j for j in range(n) if j!=i]))] for i in range(n)])
    colors=ids(observations);history=[]
    while True:
        history.append([len({colors[j] for j,(size,_) in enumerate(states) if size==n}) for n in range(6)])
        new=ids([(colors[j],tuple(colors[k] for k in actions[j])) for j in range(len(states))])
        if len(set(new))==len(set(colors)):break
        colors=new
    need(history[-1]==[1,1,1,2,11,208],'stable C3 counts')
    for n in range(6):
        c3_to_color={};color_to_c3={}
        for j,(size,_) in enumerate(states):
            if size!=n:continue
            sig=triangles(matrices[j]);col=colors[j]
            need(sig not in c3_to_color or c3_to_color[sig]==col,'C3 implies color')
            need(col not in color_to_c3 or color_to_c3[col]==sig,'color implies C3')
            c3_to_color[sig]=col;color_to_c3[col]=sig
    groups=defaultdict(list)
    for j,state in enumerate(states):groups[observations[j]].append(j)
    hostile=None
    for group in groups.values():
        for u,v in combinations(group,2):
            if colors[u]!=colors[v]:
                hostile=(u,v);break
        if hostile:break
    need(hostile is not None,'scalar hostile exists')
    u,v=hostile;remaining=list(range(states[u][0]));trace=[]
    while observations[u]==observations[v]:
        i=next(i for i,(x,y) in enumerate(zip(actions[u],actions[v])) if colors[x]!=colors[y])
        trace.append(remaining.pop(i));u,v=actions[u][i],actions[v][i]
    witness={'source_states':[states[x] for x in hostile],'deleted_original_labels':trace,'final_observations':[observations[u],observations[v]]}
    basegroups=defaultdict(list)
    for c in range(8):
        a=tournament(3,c);basegroups[triangles(a)].append(a)
    checked=0
    for group in basegroups.values():
        for a in group:
            for codes in product(range(1),range(2),range(8)):
                children=[tournament(i+1,c) for i,c in enumerate(codes)]
                x=inflate(group[0],children);y=inflate(a,children)
                need(triangles(x)==triangles(y),'inflation C3')
                need(paths(x)[0]==paths(y)[0],'inflation all subset H');checked+=1
    print(json.dumps({'status':'PASS','refinement_counts_by_round':history,'states':len(states),'inflation_checks':checked,'scalar_hostile':witness},sort_keys=True))
if __name__=='__main__':main()
