#!/usr/bin/env python3
"""Find literal three-black certificates in a declared block-word bank."""
from itertools import combinations, permutations
from pathlib import Path
import json
import sys


def arc(T,u,v):return T[u]>>v&1
def directed(T,q):return all(arc(T,u,v) for u,v in zip(q,q[1:]))
def blue(T,q):
    r=dict(zip(q,q[::-1]))
    return all(arc(T,u,v)==arc(T,r[v],r[u]) for u in range(len(T)) for v in range(u+1,len(T)))
def edge(q):return frozenset((min(u,v),max(u,v)) for u,v in zip(q,q[1:]))
def move(T,q):
    E=edge(q);n=len(T)
    return tuple(sum((arc(T,u,v)^int((min(u,v),max(u,v)) not in E))<<v
                     for v in range(n) if v!=u) for u in range(n))
def word_from_edges(n,E):
    if len(E)!=n-1:return None
    adj=[[] for _ in range(n)]
    for u,v in E:adj[u].append(v);adj[v].append(u)
    ends=[i for i,a in enumerate(adj) if len(a)==1]
    if len(ends)!=2 or any(len(a) not in (1,2) for a in adj):return None
    for a in ends:
        q=[a];previous=-1
        while len(q)<n:
            options=[v for v in adj[q[-1]] if v!=previous]
            if len(options)!=1 or options[0] in q:break
            previous=q[-1];q.append(options[0])
        if len(q)==n:return q
    return None
def candidates(n,blocks=3):
    out=set()
    for cuts in combinations(range(1,n),blocks-1):
        cuts=(0,)+cuts+(n,)
        chunks=[tuple(range(a,b)) for a,b in zip(cuts,cuts[1:])]
        for order in permutations(chunks):out.add(sum(order,()))
    return sorted(out)
def find(T,blocks=3):
    n=len(T);P=tuple(range(n));target=move(T,P);EP=edge(P)
    bank=[(q,edge(q)) for q in candidates(n,blocks)]
    tried=0
    for q1,E1 in bank:
        if not directed(T,q1) or blue(T,q1):continue
        T1=move(T,q1)
        for q2,E2 in bank:
            if not directed(T1,q2):continue
            tried+=1
            E3=EP^E1^E2
            # Test degree/connectivity before constructing a whole third source.
            q3=word_from_edges(n,E3)
            if q3 is None or blue(T1,q2):continue
            T2=move(T1,q2)
            if not directed(T2,q3):q3=q3[::-1]
            if not directed(T2,q3) or blue(T2,q3):continue
            if move(T2,q3)!=target:raise RuntimeError('literal target mismatch')
            return {'words':[q1,q2,q3],'candidate_words':len(bank),'legal_black_pairs_tried':tried}
    return {'words':None,'candidate_words':len(bank),'legal_black_pairs_tried':tried}
def main():
    base=Path(__file__).resolve().parent
    rows=json.loads((base.parent/'port_residual_integral_20261025/low_orders.normal.json').read_text())['failed_ledger_hostiles']
    out=[]
    for row in rows:
        if row['n']!=14:continue
        result=find(tuple(row['source_rows']),int(sys.argv[1]) if len(sys.argv)>1 else 3)
        result.update({'n':14,'ports':row['ports'],'source_rows':row['source_rows']})
        out.append(result)
    print(json.dumps({'status':'FINITE-EXACT SEARCH','blocks':int(sys.argv[1]) if len(sys.argv)>1 else 3,'scope':'four literal sources, forward contiguous blocks from the marked path','results':out},indent=2))
if __name__=='__main__':main()
