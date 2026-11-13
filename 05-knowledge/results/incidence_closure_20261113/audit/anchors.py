#!/usr/bin/env python3
"""Independent exact budget-eight rejection, without producer imports.
Regenerate comparisons from ordered-cover vectors and integer rows.
Search residual column antichains with fresh candidate-neighborhood packing.
"""
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import json
HERE=Path(__file__).resolve().parent
PACKET=HERE.parent
checks=0
def need(ok,message):
    global checks
    checks+=1
    if not ok: raise RuntimeError(message)
def union(xs):
    z=0
    for x in xs:z|=x
    return z

def can_cover(columns,target,budget):
    states=0
    @lru_cache(None)
    def search(left,remaining):
        nonlocal states
        states+=1
        if not left:return True
        if not remaining:return False
        work=sorted({x & left for x in columns}-{0},key=lambda x:(x.bit_count(),x),reverse=True)
        kept=[]
        for x in work:
            if not any(x & y == x for y in kept):kept.append(x)
        if not kept or union(kept)!=left:return False
        if kept[0].bit_count()*remaining<left.bit_count():return False
        if remaining==1:return kept[0]==left
        demands=[];bits=left
        while bits:
            bit=bits & -bits;bits-=bit
            neighbors=tuple(j for j,x in enumerate(kept) if x&bit)
            demands.append((len(neighbors),bit,neighbors))
        demands.sort();occupied=set();packed=0
        for _,_,nb in demands:
            if occupied.isdisjoint(nb):
                occupied.update(nb);packed+=1
                if packed>remaining:return False
        return any(search(left & ~kept[j],remaining-1) for j in demands[0][2])
    value=search(target,budget)
    return value,states

def controls():
    cases=0
    for r in range(6):
        for cols in combinations(range(1,8),r):
            target=union(cols)
            for k in range(4):
                ref=any(union(cc)==target for j in range(min(k,len(cols))+1)
                        for cc in combinations(cols,j))
                need(can_cover(cols,target,k)[0]==ref,'literal small cover reference')
                cases+=1
    return cases

def main():
    datafile=PACKET/'closure/data.json';claimsfile=PACKET/'closure/analysis.json'
    data=json.loads(datafile.read_text());claimed=json.loads(claimsfile.read_text())
    row=next(r for r in data['rows'] if r['n']==7)
    cert=next(r for r in claimed['rows'] if r['n']==7)
    nodes=row['merged_nodes'];codes=[x['code'] for x in nodes];idx={x:i for i,x in enumerate(codes)}
    n=len(nodes);A=[[0]*n for _ in nodes];B=[[0]*n for _ in nodes]
    for a,b,x in row['merged_A']:A[idx[a]][idx[b]]=x
    for a,b,x in row['merged_blue']:B[idx[a]][idx[b]]=x
    degree=list(map(sum,A))
    pairs=[(a,b) for a,b in combinations(range(n),2)
           if nodes[a]['ordered_covers']==nodes[b]['ordered_covers']]
    need(len(pairs)==184,'declared seven-vertex comparison universe')
    results={}
    for name in ('A','P','A_colored','P_colored'):
        norm=name.startswith('P');colored=name.endswith('colored')
        mats=[A,B] if colored else [A]
        columns=[]
        for c in range(n):
            mask=0
            for k,(a,b) in enumerate(pairs):
                da,db=(degree[a],degree[b]) if norm else (1,1)
                if any(M[a][c]*db!=M[b][c]*da for M in mats):mask |= 1<<k
            columns.append(mask)
        target=(1<<len(pairs))-1
        witness=cert['anchor_optimization'][name]['columns']
        need(len(witness)==9 and len(set(witness))==9,'nine distinct anchors')
        need(union(columns[idx[c]] for c in witness)==target,'literal nine-anchor upper bound')
        possible,states=can_cover(columns,target,8)
        need(not possible,'independent exact eight-anchor infeasibility')
        results[name]={'minimum':9,'budget_eight_feasible':possible,
                       'independent_search_states':states,'witness':witness}
    cases=controls()
    return {'status':'PASS; independent exact n7 lower and upper bounds',
            'controls':cases,'results':results,'checks':checks,
            'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
            'data_sha256':sha256(datafile.read_bytes()).hexdigest(),
            'claims_sha256':sha256(claimsfile.read_bytes()).hexdigest()}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
