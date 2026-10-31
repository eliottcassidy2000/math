#!/usr/bin/env python3
"""Independent literal-permutation OCF and reflection-carry audit."""
from itertools import combinations, permutations, product
from pathlib import Path
import json
import sys

GATES=0

def need(ok,msg):
    global GATES
    GATES+=1
    if not ok:raise RuntimeError(msg)


def analyze(rows,rho):
    n=len(rows)
    arc=lambda u,v:bool(rows[u]&(1<<v))
    need(all(arc(u,v)!=arc(v,u) for u,v in combinations(range(n),2)),'tournament')
    need(all(rho[rho[u]]==u for u in range(n)),'involution')
    need(all(arc(u,v)==arc(rho[v],rho[u]) for u in range(n) for v in range(n)),'anti-automorphism')
    h=sum(all(arc(word[k],word[k+1]) for k in range(n-1)) for word in permutations(range(n)))
    cycles=[]
    for k in range(3,n+1,2):
        for vertices in combinations(range(n),k):
            for tail in permutations(vertices[1:]):
                word=(vertices[0],)+tail
                if all(arc(word[t],word[(t+1)%k]) for t in range(k)):cycles.append(word)
    lookup={c:i for i,c in enumerate(cycles)}
    supports=[set(c) for c in cycles]
    mirror=[]
    for c in cycles:
        d=tuple(rho[v] for v in c[::-1]);j=d.index(min(d));mirror.append(lookup[d[j:]+d[:j]])
    q=sum(mirror[j]==j for j in range(len(cycles)))
    intersect=sum(j<mirror[j] and bool(supports[j]&supports[mirror[j]]) for j in range(len(cycles)))
    disjoint=sum(j<mirror[j] and not(supports[j]&supports[mirror[j]]) for j in range(len(cycles)))
    alpha=[1];fixed=[1]
    for size in range(1,n//3+1):
        packs=[]
        for pack in combinations(range(len(cycles)),size):
            if sum(len(supports[j]) for j in pack)==len(set().union(*(supports[j] for j in pack))):packs.append(pack)
        alpha.append(len(packs));fixed.append(sum(tuple(sorted(mirror[j] for j in pack))==pack for pack in packs))
    alpha.extend([0]*max(0,4-len(alpha)));fixed.extend([0]*max(0,4-len(fixed)))
    need(alpha[1]==q+2*intersect+2*disjoint,'single-cycle orbit types')
    need(fixed[2]==disjoint,'fixed two-packings')
    need((alpha[2]-disjoint)%2==0,'two-packing orbit integrality')
    pairorbits=(alpha[2]-disjoint)//2;t=fixed[3]
    need(h==sum((1<<k)*a for k,a in enumerate(alpha)),'OCF vs literal words')
    need(h%8==(1+2*q+4*intersect)%8,'mod8')
    need(h%16==(1+2*q+4*intersect+8*(disjoint+pairorbits+t))%16,'mod16')
    need((h>>2)&1==(q//2+intersect)%2,'third digit lower carry')
    need((h>>3)&1==((q+2*intersect)//4+disjoint+pairorbits+t)%2,'fourth digit lower carry')
    if n<=8:need(h==1+2*q+4*intersect+8*(disjoint+pairorbits),'small-order exact identity')
    if n%2==0:need(q==t==0 and h%4==1,'even-order fixed odd packing obstruction')
    while len(alpha)>1 and alpha[-1]==0:alpha.pop()
    return {'n':n,'H':h,'alpha':alpha,'q':q,'i':intersect,'j':disjoint,'p':pairorbits,'t':t}


def family(n):
    pairs=list(combinations(range(n),2));index={p:i for i,p in enumerate(pairs)}
    orbits=[];seen=set()
    for j,(u,v) in enumerate(pairs):
        if j not in seen:
            mate=index[(n-1-v,n-1-u)];o={j,mate};seen|=o;orbits.append(o)
    for bits in product(range(2),repeat=len(orbits)):
        rows=[0]*n
        for o,b in zip(orbits,bits):
            for j in o:
                u,v=pairs[j]
                if b:rows[u]|=1<<v
                else:rows[v]|=1<<u
        yield rows


def main():
    if '--test-failing-gate' in sys.argv:need(False,'intentional failure control')
    counts={}
    for n in range(1,7):
        count=0
        for rows in family(n):analyze(rows,tuple(reversed(range(n))));count+=1
        counts[n]=count
    receipt=json.loads((Path(__file__).resolve().parent.parent/'binary_carries'/'normal.json').read_text())
    checked={}
    for key,fixture in receipt['hostiles'].items():
        for d in fixture if isinstance(fixture,list) else [fixture]:
            identity=tuple(d['rows'])
            if identity not in checked:checked[identity]=analyze(d['rows'],d['rho'])
            fresh=checked[identity]
            for field in ['n','H','alpha','q','i','j','p','t']:need(fresh[field]==d[field],'hostile data field '+field)
    print(json.dumps({'status':'PASS','gates':GATES,'exhaustive_standard_reflection_counts':counts,
        'scope':'all 599 reversal-compatible tournaments n<=6, plus all eight distinct named hostile tournaments including order9 triangle stack',
        'hostile_tournaments':list(checked.values())},indent=2,sort_keys=True))

if __name__=='__main__':main()
