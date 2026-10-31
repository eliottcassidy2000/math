#!/usr/bin/env python3
"""Independent exact witness, cut-flow, component, and root-word checks."""
from itertools import product
from pathlib import Path
from math import isqrt
import json
import sys

GATES=0
BLOCK={0:0,1:0,2:1,7:1,3:2,6:2,4:3,5:3}
FLOOR=[4,7,8,7,8,7,8,7]

def need(ok,msg):
    global GATES
    GATES+=1
    if not ok:raise RuntimeError(msg)


def check_cycle(fixture):
    n=fixture['n'];cycle=fixture['cycle']
    need(len(cycle)==n and set(cycle)==set(range(1,n+1)),'Hamiltonian vertex set')
    edges=list(zip(cycle,cycle[1:]+cycle[:1]));roots=[]
    for u,v in edges:
        s=isqrt(u+v);need(s*s==u+v,'nonsquare edge');roots.append(s)
    need(cycle[0] in compiled(roots,n),'long Hamiltonian compiler positive control')
    need(compiled(roots,n)==direct(roots,n),'long compiler vs direct walking')
    cuts=[j for j,s in enumerate(roots) if s%2==0]
    need(len(cuts)==fixture['even_square_edges'] and len(cuts)>=FLOOR[n%8],'seam count and lower bound')
    if fixture.get('sharp',True):need(len(cuts)==FLOOR[n%8],'sharp seam count')
    adjacency={v:[] for v in cycle}
    endpoints=[0]*8
    for (u,v),s in zip(edges,roots):
        if s%2:adjacency[u].append(v);adjacency[v].append(u)
        else:endpoints[u%8]+=1;endpoints[v%8]+=1
    unseen=set(cycle);ks=[0]*4;components=[]
    while unseen:
        root=min(unseen);comp={root};stack=[root];unseen.remove(root)
        while stack:
            for v in adjacency[stack.pop()]:
                if v in unseen:unseen.remove(v);comp.add(v);stack.append(v)
        block=BLOCK[root%8]
        need(all(BLOCK[v%8]==block for v in comp),'odd-root component crosses block')
        ks[block]+=1;components.append(sorted(comp))
    need(sum(ks)==len(cuts) and min(ks)>=1,'cut forest components')
    nr=[sum(v%8==r for v in range(1,n+1)) for r in range(8)]
    deltas=[nr[1]-nr[0],nr[7]-nr[2],nr[3]-nr[6],nr[5]-nr[4]]
    for i,(positive,negative) in enumerate([(1,0),(7,2),(3,6),(5,4)]):
        need(endpoints[positive]==ks[i]+deltas[i] and endpoints[negative]==ks[i]-deltas[i],'endpoint incidence formula')
    x,y=ks[0]+ks[3],ks[1]+ks[2];dl,dr=deltas[0]+deltas[3],deltas[1]+deltas[2]
    need(x-y==dr-dl and (x-dl)%2==0 and (y-dr)%2==0,'aggregate flow congruences')
    # Delete each original edge to audit every Hamiltonian-path endpoint formula.
    for (u,v),s in zip(edges,roots):
        b=[0]*8;b[u%8]+=1;b[v%8]+=1
        pk=ks[:]
        if s%2:pk[BLOCK[u%8]]+=1
        px,py=pk[0]+pk[3],pk[1]+pk[2]
        need(px-py==dr-dl+b[1]+b[5]-b[3]-b[7],'path imbalance')
        need((px-dl-b[0]-b[4])%2==0 and (py-dr-b[2]-b[6])%2==0,'path endpoint parity')
        need(sum(pk)==len(cuts)-(s%2==0)+1,'path component count')
        if BLOCK[u%8]==BLOCK[v%8]:need(pk[BLOCK[u%8]]>=2,'same-block endpoints connectedness')
    return {'n':n,'seams':len(cuts),'path_counts':ks,'endpoint_residues':endpoints,'singleton_paths':sum(len(c)==1 for c in components)}


def compiled(word,n):
    constants=[0]
    for s in word:constants.append(s*s-constants[-1])
    m=len(word)
    if m%2==0 and constants[-1]!=0:return []
    lower,upper=1,n
    for i,c in enumerate(constants[:-1]):
        lo,hi=(1-c,n-c) if i%2==0 else (c-n,c-1)
        lower=max(lower,lo);upper=min(upper,hi)
    if m%2:
        if constants[-1]%2:return []
        lower=max(lower,constants[-1]//2);upper=min(upper,constants[-1]//2)
    forbidden=set()
    for i in range(m):
        for j in range(i+1,m):
            if (i-j)%2==0:
                if constants[i]==constants[j]:return []
            else:
                even,odd=(i,j) if i%2==0 else (j,i)
                d=constants[odd]-constants[even]
                if d%2==0:forbidden.add(d//2)
    return [x for x in range(lower,upper+1) if x not in forbidden]


def direct(word,n):
    starts=[]
    for x in range(1,n+1):
        walk=[x]
        for s in word:walk.append(s*s-walk[-1])
        if walk[-1]==x and len(set(walk[:-1]))==len(word) and all(1<=v<=n for v in walk[:-1]):starts.append(x)
    return starts


def odd_components(n):
    parent=list(range(n+1))
    def root(v):
        while parent[v]!=v:v=parent[v]
        return v
    edges=[]
    for u in range(1,n+1):
        for s in range(3,isqrt(2*n)+1,2):
            v=s*s-u
            if u<v<=n:
                edges.append((u,v));parent[root(u)]=root(v)
    comps={}
    for u in range(1,n+1):comps.setdefault(root(u),[]).append(u)
    return sorted(sorted(c) for c in comps.values()),edges


def main():
    if '--test-failing-gate' in sys.argv:need(False,'intentional failure control')
    data=json.loads((Path(__file__).resolve().parent.parent/'square_sum'/'witnesses.json').read_text())
    rows=[check_cycle(x) for x in data['cycles']]
    need([r['n'] for r in rows]==list(range(80,88)),'eight consecutive residue witnesses')
    control32=check_cycle({'n':32,'even_square_edges':16,'sharp':False,'cycle':[1,8,28,21,4,32,17,19,30,6,3,13,12,24,25,11,5,31,18,7,29,20,16,9,27,22,14,2,23,26,10,15]})
    need(control32['singleton_paths']>0,'singleton endpoint incidence control')
    comps,edges=odd_components(40)
    expected=[]
    for a in range(1,5):expected.extend([[24+a,25-a,a,9-a,16+a,33-a],[32+a,17-a,8+a,41-a]])
    need(comps==sorted(sorted(c) for c in expected),'H40 exact paths')
    need(len(edges)==32,'H40 no unlisted odd edges')
    for n in range(44,257):
        comps,_=odd_components(n)
        desired=sorted(sorted(v for v in range(1,n+1) if BLOCK[v%8]==b) for b in range(4))
        need(comps==desired,'four-component theorem n='+str(n))
    for n in range(13,2001):
        need(any(1<=s*s-n<n for s in range(3,isqrt(2*n)+1,2)),'odd-square arrival neighbor')
    words=0
    for n in range(8,13):
        for m in range(3,7):
            for word in product(range(2,6),repeat=m):
                need(compiled(word,n)==direct(word,n),'compiled root word vs direct walking')
                words+=1
    need(compiled((3,3,3,3),20)==[],'closed-word repeated-vertex hostile')
    print(json.dumps({'status':'PASS','gates':GATES,'sharp_cycles':rows,'nonsharp_singleton_control32':control32,'root_word_cases':words,
        'odd_component_scope':'H40 explicit base; every N=44..256; arrival13..2000'},indent=2,sort_keys=True))

if __name__=='__main__':main()
