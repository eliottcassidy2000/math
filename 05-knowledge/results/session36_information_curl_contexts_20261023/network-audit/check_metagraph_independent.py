#!/usr/bin/env python3
"""Independent fixed-path enumeration, half-form, and exact determinant audit.
No research-worker code is imported. Python integers throughout determinants.
"""
import itertools, json, time, hashlib
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent

def check(ok,msg):
    if not ok: raise RuntimeError(msg)

def binary_kernel(A):
    n=len(A[0]) if len(A) else 0
    rows=[sum((int(x)&1)<<j for j,x in enumerate(row)) for row in A]
    piv=[];r=0
    for j in range(n):
        q=next((k for k in range(r,len(rows)) if rows[k]>>j&1),None)
        if q is None: continue
        rows[r],rows[q]=rows[q],rows[r]
        for k in range(len(rows)):
            if k!=r and rows[k]>>j&1: rows[k]^=rows[r]
        piv.append(j);r+=1
    basis=[]
    for j in set(range(n))-set(piv):
        v=1<<j
        for k,p in enumerate(piv):
            if rows[k]>>j&1: v|=1<<p
        basis.append([(v>>k)&1 for k in range(n)])
    return r,basis

def laplacian(n,edges):
    L=[[0]*n for _ in range(n)]
    for (a,b),w in edges.items():
        if a==b: continue
        L[a][a]+=w; L[b][b]+=w; L[a][b]-=w; L[b][a]-=w
    return L

def components(n,edges):
    parent=list(range(n))
    def root(x):
        while parent[x]!=x: parent[x]=parent[parent[x]];x=parent[x]
        return x
    for a,b in edges:
        a,b=root(a),root(b)
        if a!=b: parent[b]=a
    return len({root(x) for x in range(n)})

def determinant(A):
    A=[r[:] for r in A];d=1
    for k in range(len(A)-1):
        p=A[k][k];check(p>0,'non-positive principal pivot')
        for i in range(k+1,len(A)):
            x=A[i][k]
            for j in range(i,len(A)):
                z,rem=divmod(p*A[i][j]-x*A[k][j],d)
                check(rem==0,'Bareiss inexact division')
                A[i][j]=z;A[j][i]=z
        for i in range(k+1,len(A)): A[i][k]=A[k][i]=0
        d=p
    return A[-1][-1] if A else 1

def enumerate_graph(n):
    chords=[(i,j) for i in range(n) for j in range(i+2,n)]
    cache={}
    def canonical(rows):
        rows=tuple(rows)
        if rows in cache:return cache[rows]
        # Opposite ordering/encoding to the research census: descending degrees,
        # lexicographic full outgoing adjacency rows, not upper-triangle integers.
        cells=defaultdict(list)
        for i,row in enumerate(rows):cells[row.bit_count()].append(i)
        ordered=[cells[d] for d in sorted(cells,reverse=True)]
        best=None
        for chunks in itertools.product(*(itertools.permutations(c) for c in ordered)):
            p=sum(chunks,())
            word=tuple(sum(((rows[u]>>v)&1)<<j for j,v in enumerate(p)) for u in p)
            if best is None or word<best:best=word
        cache[rows]=best
        return best
    vertex=[];blue=[];verts={}
    for x in range(1<<len(chords)):
        rows=[(1<<(i+1)) if i<n-1 else 0 for i in range(n)]
        for k,(a,b) in enumerate(chords):
            if x>>k&1:rows[a]|=1<<b
            else:rows[b]|=1<<a
        opp=[((1<<n)-1)^r^(1<<i) for i,r in enumerate(rows)]
        key=min(canonical(rows),canonical(opp))
        if key not in verts:verts[key]=len(verts)
        vertex.append(verts[key])
        # Literal reversal-converse equality on every ordered pair.
        blue.append(all(((rows[a]>>b)&1)==((rows[n-1-b]>>(n-1-a))&1) for a in range(n) for b in range(n)))
    B=Counter();K=Counter();full=len(vertex)-1
    for x in range((full+1)//2):
        y=full-x;e=tuple(sorted((vertex[x],vertex[y])))
        check(blue[x]==blue[y],'inconsistent color under complement')
        (B if blue[x] else K)[e]+=1
    check(all(w%2==0 for w in K.values()),'black multiplicity odd')
    return len(verts),B,K

def run():
    expected=json.loads((HERE.parent/'metagraph/dyadic_network_receipt.json').read_text())
    expected_lifts=json.loads((HERE.parent/'metagraph/dyadic_lift_receipt.json').read_text())
    out=[]
    for source,lifts in zip(expected,expected_lifts):
        start=time.time();n=source['n'];N,B,K=enumerate_graph(n)
        Lb=laplacian(N,B);Lbnp=np.array(Lb,dtype=np.int64)
        L=laplacian(N,B+K);A=[r[:-1] for r in L[:-1]]
        Ab=np.array(A,dtype=np.int64)
        rank,vecs=binary_kernel(A);rankB,_=binary_kernel(Lb)
        cb=components(N,B);bicycle=N-cb-rankB;crit=N-1-rank
        if vecs:
            V=np.array(vecs,dtype=np.int64);raw=V@Ab@V.T
            check(np.all(raw%2==0),'half form not integral')
            beta=(raw//2)%2;brank,_=binary_kernel(beta);rad=len(vecs)-brank
            # Alternate integral lifts are indistinguishable. Deterministic probe.
            W=V+2*np.array([[(3*i+5*j+1)%7-3 for j in range(N-1)] for i in range(len(vecs))],dtype=np.int64)
            check(np.array_equal((W@Ab@W.T//2)%2,beta),'half form depends on integral lifts')
        else:brank=rad=0
        tau=determinant(A);v2=(tau&-tau).bit_length()-1
        check(crit==cb-1+bicycle,'rank-bicycle identity')
        check(v2==crit+rad,'unexpected exponent >= 3 in actual graph')
        vals={str(0):rank}
        if crit-rad:vals['1']=crit-rad
        if rad:vals['2']=rad
        row=dict(n=n,nodes=N,blue_lines=sum(B.values()),black_lines=sum(K.values()),blue_components=cb,blue_bicycle_dimension=bicycle,critical_two_rank=crit,half_form_rank=brank,half_form_radical=rad,valuation_histogram=vals,spanning_tree_valuation2=v2,spanning_tree_count=str(tau),elapsed_seconds=round(time.time()-start,3))
        for k in ('nodes','blue_lines','black_lines','blue_components','blue_bicycle_dimension','critical_two_rank'):
            check(row[k]==source[k],f'census mismatch {n}:{k}')
        check(vals==lifts['valuation_histogram'],'Smith census mismatch')
        check(v2==lifts['spanning_tree_valuation2'],'tree valuation mismatch')
        if 'spanning_trees' in source:check(tau==source['spanning_trees'],'tree count mismatch')
        out.append(row); print({k:v for k,v in row.items() if k!='spanning_tree_count'},flush=True)
        (HERE/'metagraph_independent_receipt.json').write_text(json.dumps(dict(status='PASS',rows=out),indent=2)+'\n')
    print('ALL INDEPENDENT METAGRAPH CHECKS PASSED',flush=True)
if __name__=='__main__': run()
