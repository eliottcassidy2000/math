#!/usr/bin/env python3
"""Independent small weighted-quotient audit of the exact six-defect boundary.
No cube, graph-isomorphism, or global-sign engine imports.
"""
from itertools import combinations,permutations
from collections import defaultdict
import json
CHECKS=0

def need(x,msg):
    global CHECKS
    CHECKS+=1
    if not x:raise RuntimeError(msg)

def canon(weights,A):
    k=len(weights);return min((tuple(weights[i] for i in q),tuple(A[q[i]][q[j]] for i,j in combinations(range(k),2))) for q in permutations(range(k)))

def collapse(weights,A):
    vertices=[]
    for i,w in enumerate(weights):vertices += [i]*w
    G=[[int(i!=j and A[vertices[i]][vertices[j]]) for j in range(len(vertices))] for i in range(len(vertices))]
    classes=defaultdict(list)
    for i,row in enumerate(G):classes[tuple(row)].append(i)
    cs=list(classes.values());ws=[len(c) for c in cs];B=[[G[c[0]][d[0]] for d in cs] for c in cs]
    return ws,B

def analyze(weights,A,triples=2):
    weights,A=collapse(weights,A);k=len(weights)
    eligible=[i for i,w in enumerate(weights) if w>=5 and sum((v-5 if j==i else v)//3 for j,v in enumerate(weights))>=triples]
    aut=[]
    for q in permutations(range(k)):
        if all(weights[i]==weights[q[i]] for i in range(k)) and all(A[i][j]==A[q[i]][q[j]] for i,j in combinations(range(k),2)):
            sg=-1 if sum(weights[i]*weights[j]*A[i][j] for i,j in combinations(range(k),2) if q[i]>q[j])%2 else 1
            aut.append((q,sg))
    even=all(sg==1 for _,sg in aut)
    rel=[i for i in eligible if all(sg==1 for q,sg in aut if q[i]==i)]
    seen=set();orbits=[]
    for i in rel:
        if i in seen:continue
        orb=sorted({q[i] for q,_ in aut});seen.update(orb);orbits.append(orb)
    return {'weights':weights,'quotient_edges':[[i,j] for i,j in combinations(range(k),2) if A[i][j]],'eligible_classes':eligible,'globally_even':even,'relative_copy_orbits':orbits,'relative_copy_multiplicity':len(orbits),'fusion_defect':max(0,len(orbits)-1) if even else 0,'escape_orbits':len(orbits) if not even else 0,'key':canon(weights,A)}

def main():
    cases={}
    # Other large-class partitions have at most one eligible class and no
    # equal-weight odd class swap, by the written support-capacity proof.
    for weights in [(5,5,3),(6,6,1),(5,8),(6,7)]:
        pairs=list(combinations(range(len(weights)),2))
        for s in range(1<<len(pairs)):
            A=[[0]*len(weights) for _ in weights]
            for i,(u,v) in enumerate(pairs):A[u][v]=A[v][u]=(s>>i)&1
            r=analyze(weights,A)
            cases[r['key']]=r
    anomalies=[r for r in cases.values() if r['fusion_defect'] or r['escape_orbits']]
    need(sum(r['fusion_defect'] for r in anomalies)==4,'expected four target copy-fusion losses')
    need(sum(r['escape_orbits'] for r in anomalies)==2,'expected two escaping copy orbits')
    need(sum(r['globally_even'] for r in anomalies)==4 and len(anomalies)==6,'expected six anomaly types')
    # Hostile: K55 with two isolated points admits one, not two, triples off K.
    A=[[0,1,0],[1,0,0],[0,0,0]]
    small=analyze([5,5,2],A);need(not small['eligible_classes'],'n12 double-triple hostile incorrectly admitted')
    one=analyze([5,5,2],A,triples=1);need(one['escape_orbits']==1,'single-triple escape control missing')
    boundary=analyze([5,5,3],A);need(boundary['escape_orbits']==1,'n13 double-triple escape missing')
    for r in anomalies:r.pop('key')
    print(json.dumps({'status':'PASS','scope':'weighted twin-quotient anomalies proved exhaustive in SHARP_N13_BOUNDARY.md, not an all-graph enumeration','checks':CHECKS,'distinct_candidate_global_types':len(cases),'target_copy_fusion_defect':4,'escaping_relative_copy_orbits':2,'exact_predicate_deficit':6,'anomalies':sorted(anomalies,key=lambda x:(x['weights'],x['quotient_edges']))},sort_keys=True,indent=2))
if __name__=='__main__':main()
