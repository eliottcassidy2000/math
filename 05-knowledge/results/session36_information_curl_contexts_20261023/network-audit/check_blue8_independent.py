#!/usr/bin/env python3
"""Bounded independent blue-only order-eight enumeration, not full G8 census."""
import itertools,json,time
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
from check_metagraph_independent import check,binary_kernel,laplacian,components,determinant
from symmetric_dyadic import valuations
HERE=Path(__file__).resolve().parent

def run():
    start=time.time();n=8
    chords=[(a,b) for a in range(n) for b in range(a+2,n)];remaining=set(chords);orbits=[]
    while remaining:
        e=min(remaining);rev=(n-1-e[1],n-1-e[0]);orb={e,rev};orbits.append(sorted(orb));remaining-=orb
    cache={}
    def canonical(rows):
        key=tuple(rows)
        if key in cache:return cache[key]
        cells=defaultdict(list)
        for i,row in enumerate(rows):cells[row.bit_count()].append(i)
        best=None
        for chunks in itertools.product(*(itertools.permutations(cells[d]) for d in sorted(cells,reverse=True))):
            p=sum(chunks,());word=tuple(sum(((rows[u]>>v)&1)<<j for j,v in enumerate(p)) for u in p)
            if best is None or word<best:best=word
        cache[key]=best;return best
    nodes={};labels=[]
    for state in range(1<<len(orbits)):
        rows=[(1<<(i+1)) if i<n-1 else 0 for i in range(n)]
        for k,orb in enumerate(orbits):
            for a,b in orb:
                u,v=(a,b) if state>>k&1 else (b,a);rows[u]|=1<<v
        check(all(((rows[a]>>b)&1)==((rows[n-1-b]>>(n-1-a))&1) for a in range(n) for b in range(n)),'not literal blue')
        key=canonical(rows)
        if key not in nodes:nodes[key]=len(nodes)
        labels.append(nodes[key])
    edges=Counter();full=len(labels)-1
    for x in range(len(labels)//2):edges[tuple(sorted((labels[x],labels[full-x])))]+=1
    N=len(nodes);L=laplacian(N,edges);cc=components(N,edges);rank,_=binary_kernel(L)
    check(cc==1,'blue8 not connected');A=[r[:-1] for r in L[:-1]];_,basis=binary_kernel(A)
    V=np.array(basis,dtype=np.int64);Q=(V@np.array(A,dtype=np.int64)@V.T//2)%2;brank,_=binary_kernel(Q)
    vals=[v for v in valuations(A) if v];tau=determinant(A)
    check(sum(vals)==(tau&-tau).bit_length()-1,'blue8 determinant valuation mismatch')
    masks=[sum(x<<i for i,x in enumerate(v)) for v in basis];best=None
    for state in range(1<<len(masks)):
        mask=0
        for i,m in enumerate(masks):
            if state>>i&1:mask^=m
        cut=sum(w for (u,v),w in edges.items() if ((mask>>u)^(mask>>v))&1)
        if cut%4==2 and (best is None or cut<best[0]):best=(cut,mask)
    check(best is not None,'missing nonalternating witness')
    out=dict(status='PASS',n=n,blue_presentations=len(labels),selfconverse_nodes=N,components=cc,bicycle_dimension=N-cc-rank,two_primary_exponents=vals,half_form_rank=brank,alternating=not any(Q[i,i] for i in range(len(basis))),minimum_nonalternating_cut=best[0],binary_kernel_states=1<<len(masks),tree_count=str(tau),elapsed_seconds=round(time.time()-start,3))
    source=json.loads((HERE.parent/'metagraph/blue_only_n8.json').read_text())
    for k in ('blue_presentations','selfconverse_nodes','components','bicycle_dimension','two_primary_exponents'):check(out[k]==source[k],f'blue8 mismatch {k}')
    check(brank==source['half_forms'][0]['half_form_rank'],'blue8 form rank mismatch')
    if 'nonalternating_bicycle_witness' in source:check(best[0]==source['nonalternating_bicycle_witness']['crossing_blue_lines'],'minimal witness mismatch')
    (HERE/'blue8_independent_receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2));print('ALL INDEPENDENT BLUE8 CHECKS PASSED')
if __name__=='__main__':run()
