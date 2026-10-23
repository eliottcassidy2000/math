#!/usr/bin/env python3
"""Direct graph orbit and automorphism-parity test, no tournament enumeration imports."""
import itertools,json
from pathlib import Path
OUT=Path(__file__).resolve().parent

def need(x,why):
    if not x:raise RuntimeError(why)

def main():
    rows=[];expected_even=[0,1,1,2,4,12,56];expected_connected=[0,1,0,1,2,8,43]
    for n in range(1,7):
        pairs=list(itertools.combinations(range(n),2));where={p:i for i,p in enumerate(pairs)}
        transforms=[]
        for p in itertools.permutations(range(n)):
            mapping=[where[tuple(sorted((p[i],p[j])))] for i,j in pairs]
            inverted=sum((p[i]>p[j])<<k for k,(i,j) in enumerate(pairs))
            transforms.append((mapping,inverted))
        def apply(b,mp):return sum(((b>>i)&1)<<j for i,j in enumerate(mp))
        left=set(range(1<<len(pairs)));allcount=evens=conn=0
        while left:
            b=min(left);orb={apply(b,mp) for mp,inv in transforms};left-=orb;allcount+=1
            odd=any(apply(b,mp)==b and (b&inv).bit_count()%2 for mp,inv in transforms)
            if odd:continue
            evens+=1
            adj=[set() for _ in range(n)]
            for k,(i,j) in enumerate(pairs):
                if b>>k&1:adj[i].add(j);adj[j].add(i)
            seen={0};todo=[0]
            while todo:
                for j in adj[todo.pop()]-seen:seen.add(j);todo.append(j)
            conn+=len(seen)==n
        need(evens==expected_even[n],('even graph equality',n,evens))
        need(conn==expected_connected[n],('primitive connected equality',n,conn))
        row=dict(n=n,graph_isotypes=allcount,even_graph_isotypes=evens,connected_even_graph_isotypes=conn)
        rows.append(row);print(row,flush=True)
    report=dict(status='PASS',rows=rows,definition='Every automorphism reverses an even number of edges, not all vertex degrees even.')
    (OUT/'even_graph_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
if __name__=='__main__':main()
