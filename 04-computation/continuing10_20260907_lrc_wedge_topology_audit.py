#!/usr/bin/env python3
"""Independent role-cover proof controls and exhaustive labelled-tree referee."""
from pathlib import Path
from itertools import product, combinations
from hashlib import sha256
import json
import sys
sys.stdout.reconfigure(encoding='utf-8',newline='\n')
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent if HERE.name=='04-computation' else Path('C:/w/s0905')
DEST=ROOT/'05-knowledge/results' if HERE.name=='04-computation' else HERE
G=0
def need(v,label):
    global G
    G+=1
    if not v:raise ArithmeticError(label)
def readpin(path,pin=None):
    raw=path.read_bytes()
    if pin is not None:need(sha256(raw).hexdigest()==pin,'input '+str(path))
    return json.loads(raw)
OLD=readpin(ROOT/'05-knowledge/results/continuing8_20260906_lrc_minimum_tree_certificate.json','580a7c930103aab3bea867ad463a90b0e0208323a90ee95a685ff811a761582d')
WEDGE=readpin(DEST/'continuing10_20260907_lrc_composite_wedges_certificate.json','890d00d44f0195765d62fe1d40b59ad102311f34a84a42a5fddc229d037209e9')
J=readpin(DEST/'continuing10_20260907_lrc_wedge_topology_certificate.json','4008a6ae441672637becabb0f0ba17f8e0fe9e5975657fe9dddefe2bb5766e79')
need(sha256((HERE/'continuing10_20260907_lrc_wedge_topology.py').read_bytes()).hexdigest()=='6fc998e7596bfa2c4182bfe41423f3d4cbdd9d9a4f7ce77724afebdfdb653e6b','frozen producer source, never imported')
C=next(r for r in OLD['clocks'] if r['t']==7200)
WEIGHT={tuple(d):v[0] for d,v in C['weights']}
EDGES=list(combinations(range(7),2));EINDEX={e:i for i,e in enumerate(EDGES)}

# Every labelled tree exactly once, by the classical elementary Pruefer
# bijection. This referee calls no spanning-tree optimization routine.
TREES=[]
for word in product(range(7),repeat=5):
    degree=[1]*7
    for a in word:degree[a]+=1
    edges=[]
    for a in word:
        leaf=next(i for i in range(7) if degree[i]==1)
        edges.append(tuple(sorted((a,leaf))));degree[a]-=1;degree[leaf]-=1
    leafs=[i for i in range(7) if degree[i]==1]
    need(len(leafs)==2,'Pruefer final pair')
    edges.append(tuple(leafs));ids=tuple(sorted(EINDEX[e] for e in edges))
    need(len(ids)==len(set(ids))==6,'six distinct tree edges')
    reached={0}
    for _ in range(6):
        reached|={j for i,j in edges if i in reached}|{i for i,j in edges if j in reached}
    need(len(reached)==7,'each decoded graph connected')
    TREES.append(ids)
need(len(TREES)==len(set(TREES))==7**5==16807,'complete labelled-tree universe')

QUERIES=0
def complete_minimum(allowed):
    global QUERIES
    QUERIES+=1
    best=None;attainer=None;count=0
    for tree in TREES:
        if all(e in allowed for e in tree):
            count+=1;cost=sum(allowed[e] for e in tree)
            if best is None or cost<best:best=cost;attainer=tree
    return best,attainer,count

need(len(C['survivors'])==len(J['words'])==15,'all inherited surviving words')
deleted=[];remaining=[];branch_count=0
for wi,((d,E,oldM),rec) in enumerate(zip(C['survivors'],J['words'])):
    need([wi,d,E]==[rec['word_index'],rec['word'],rec['E']],'word indexing and marginal bound')
    base={k:WEIGHT[tuple(sorted((d[i],d[j])))] for k,(i,j) in enumerate(EDGES) if tuple(sorted((d[i],d[j]))) in WEIGHT}
    need(complete_minimum(base)[0]==oldM,'unrestricted inherited minimum from all trees')
    vertices=[i for i,x in enumerate(d) if x==24]
    need(vertices==rec['vertices24'],'all labelled middle vertices')
    # The local role cover is checked on ALL subsets of possible incident
    # graph edges, not merely on the producer's branch extrema.
    for v in vertices:
        neighbors=[j for j in range(7) if j!=v and EINDEX[tuple(sorted((v,j)))] in base]
        for mask in range(1<<len(neighbors)):
            N={j for k,j in enumerate(neighbors) if mask>>k&1}
            has18=any(d[j]==18 for j in N);has416=any(d[j] in [4,16] for j in N)
            forbidden=has18 and has416
            roles=[all(d[j]!=18 for j in N),all(d[j] not in [4,16] for j in N)]
            need((not forbidden)==any(roles),'exact local no-wedge iff binary role cover')
    choices=list(product(range(2),repeat=len(vertices)))
    need(len(choices)==len(rec['branches']),'every product of local roles retained')
    minima=[]
    for role,b in zip(choices,rec['branches']):
        branch_count+=1;removed=set()
        for v,r in zip(vertices,role):
            forbidden_types={18} if r==0 else {4,16}
            removed|={tuple(sorted((v,j))) for j in range(7) if d[j] in forbidden_types}
        allowed={k:cost for k,cost in base.items() if EDGES[k] not in removed}
        need(list(role)==b['choices'] and sorted(map(list,removed))==b['removed'],'exact branch deletions')
        expected=sorted([cost,*EDGES[k]] for k,cost in allowed.items())
        need(expected==b['allowed_edges'],'complete permitted weighted graph')
        best,tree,count=complete_minimum(allowed)
        need(best==b['cost'],'minimum over every labelled spanning tree')
        need((best is None or best>E)==b['closed'],'strict closure versus equality boundary')
        if best is not None:
            ids=tuple(sorted(EINDEX[tuple(sorted((i,j)))] for cost,i,j in b['tree']))
            need(ids in set(TREES) and all(e in allowed for e in ids),'saved attainer is an allowed spanning tree')
            need(sum(allowed[e] for e in ids)==best,'saved attaining cost')
        else:need(count==0,'disconnected branch has no spanning tree')
        minima.append(best)
    closed=all(m is None or m>E for m in minima)
    need(closed==rec['closed'],'whole-word branch universal quantifier')
    if closed:deleted.append(wi)
    else:remaining.append((wi,min(m for m in minima if m is not None),E))
    print('WORD',wi,'E',E,'all_branch_minima',minima,'closed',closed)
need(deleted==J['deleted_indices']==[1,2,3,4,5,6,7,8,9,10,13],'exact eleven closed connected word classes')
need(remaining==[(0,114,116),(11,66,79),(12,102,103),(14,66,87)],'all four unclosed word costs and excesses')
need([C['survivors'][i][0] for i,_,_ in remaining]==J['remaining_words'],'exact remaining words')
need(branch_count==29 and QUERIES==44,'all29 branches plus15 unrestricted controls')
print('UNIVERSE 16807 labelled trees per query;44 complete minima;29 role branches on15 words')
print('DELETED',deleted,'REMAINING',remaining)
print('SCOPE actual connected complements only; no physical realization or unsafe claim for remaining words')
print('PASS',G,'always-active independent exact gates; raw LF')
