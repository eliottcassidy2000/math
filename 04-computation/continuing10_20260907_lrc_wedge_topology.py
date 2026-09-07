#!/usr/bin/env python3
"""Actual connected-complement topology consequence of two forbidden wedges."""
from pathlib import Path
from itertools import combinations,product
from hashlib import sha256
import json
import sys
sys.stdout.reconfigure(encoding='utf-8',newline='\n')
HERE=Path(__file__).resolve();ROOT=HERE.parent.parent if HERE.parent.name=='04-computation' else Path('C:/w/s0905')
DEST=ROOT/'05-knowledge/results' if HERE.parent.name=='04-computation' else HERE.parent
G=0
def need(v,s):
    global G
    G+=1
    if not v:raise ArithmeticError(s)
def canon(x):return json.dumps(x,sort_keys=True,separators=(',',':')).encode()
raw=(ROOT/'05-knowledge/results/continuing8_20260906_lrc_minimum_tree_certificate.json').read_bytes()
need(sha256(raw).hexdigest()=='580a7c930103aab3bea867ad463a90b0e0208323a90ee95a685ff811a761582d','complete inherited minimum-tree certificate')
X=json.loads(raw);C=next(r for r in X['clocks'] if r['t']==7200);W={tuple(d):v for d,v in C['weights']}
wpath=DEST/'continuing10_20260907_lrc_composite_wedges_certificate.json'
wr=wpath.read_bytes()
need(sha256(wr).hexdigest()=='890d00d44f0195765d62fe1d40b59ad102311f34a84a42a5fddc229d037209e9','two complete native wedge certificates')
J=json.loads(wr)
need([r['margins'] for r in J['wedges']]==[[4,24,18],[16,24,18]],'precise forbidden wedge scope')
def mst(d,removed):
    edges=sorted((W[tuple(sorted((d[i],d[j])))][0],i,j) for i,j in combinations(range(7),2) if (i,j) not in removed and tuple(sorted((d[i],d[j]))) in W)
    P=list(range(7));selected=[]
    def find(a):
        while P[a]!=a:a=P[a]
        return a
    for cost,i,j in edges:
        a,b=find(i),find(j)
        if a!=b:P[a]=b;selected.append((cost,i,j))
    if len(selected)==6:
        # Independent cut optimality verifies each chosen edge's MST condition.
        for c,i,j in selected:
            seen={i};change=True
            while change:
                change=False
                for cc,u,v in selected:
                    if (cc,u,v)==(c,i,j):continue
                    if (u in seen)!=(v in seen):seen.update([u,v]);change=True
            need(j not in seen,'each selected edge is a true tree bridge')
            need(all(cc>=c for cc,u,v in edges if (u in seen)!=(v in seen)),'literal cut-property minimum tree audit')
    return (sum(e[0] for e in selected) if len(selected)==6 else None),selected,edges
results=[]
for wi,(d,E,M) in enumerate(C['survivors']):
    need(mst(d,set())[0]==M,'unrestricted inherited projected cost recovered')
    vertices=[i for i,v in enumerate(d) if v==24];branches=[]
    for choices in product([0,1],repeat=len(vertices)):
        removed=set()
        for v,choice in zip(vertices,choices):
            for j,x in enumerate(d):
                if (choice==0 and x==18) or (choice==1 and x in [4,16]):removed.add(tuple(sorted((v,j))))
        cost,tree,edges=mst(d,removed)
        branches.append(dict(choices=choices,removed=sorted(removed),cost=cost,tree=tree,allowed_edges=edges,closed=cost is None or cost>E))
    results.append(dict(word_index=wi,word=d,E=E,vertices24=vertices,branches=branches,closed=all(b['closed'] for b in branches)))
deleted=[r['word_index'] for r in results if r['closed']];remaining=[r for r in results if not r['closed']]
need(deleted==[1,2,3,4,5,6,7,8,9,10,13],'complete11-word actual topological deletion')
need([r['word_index'] for r in remaining]==[0,11,12,14],'exact remaining four words')
need(sum(len(r['branches']) for r in results)==29,'complete binary branching including no24 case')
cert=dict(status='FINITE-EXACT topology consumer; dependent native-wedge proof promotion is separate',clock=7200,
    minimum_tree_pin=sha256(raw).hexdigest(),wedges_pin=sha256(wr).hexdigest(),words=results,deleted_indices=deleted,remaining_words=[r['word'] for r in remaining])
data=canon(cert)+b'\n';(DEST/(HERE.stem+'_certificate.json')).write_bytes(data)
print('ACTUAL_GRAPH: forbid simultaneous24--18 and24--(4 or16); complete binary choices per24 vertex')
print('UNIVERSE15 inherited words;29 branches; MINIMUM tree on each permitted graph; every selected edge passes its cut condition')
for r in results:print(r['word_index'],r['word'],'E',r['E'],'branch_costs',[b['cost'] for b in r['branches']],'closed',r['closed'])
print('DELETED',deleted,'REMAINING',[r['word_index'] for r in remaining])
print('SCOPE actual connected seven-complement with selected six gcd7200; no claim on disconnected complements or all7200 closure')
print('CERTIFICATE_SHA256',sha256(data).hexdigest())
print('PASS',G,'always-active exact gates; raw LF')
