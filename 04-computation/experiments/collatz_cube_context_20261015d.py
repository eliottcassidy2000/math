#!/usr/bin/env python3
"""Exact free-terminal cube optimum and a one-edge repair. No sign oracle.
A frustrated-triangle packing upper bound meets an explicit vertex lower bound.
"""
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
from collatz_gram_slack_20261015b import observe,close_graph,T,need,transport_slack
from collatz_odd_ray_overlap_20261015d import rows

ROOT=Path(__file__).resolve().parents[2]
WITNESS_COLUMNS=[1,54244,61688,61726,65117,70417,70924,71491,71956,78007,79184,186583,189458,196289,199675,203768,204643,206686,207623,215104,220249,222416,227731,228086,231943,232703,233113,237761,510299,537299,563966,570962,575180,580259,594809,606032,622525,630949,644966,667763,678644,688814,1727648,1799714,1885436,1963304,2030629,2144794,6183539,6295643,15756743]
WITNESS_SIGNS=[1,1,-1,-1,1,-1,1,1,1,1,-1,1,-1,-1,-1,1,-1,-1,1,-1,1,-1,-1,-1,-1,1,-1,1,1,1,1,-1,1,1,-1,1,1,-1,-1,-1,-1,1,-1,-1,-1,1,-1,1,1,-1,-1]

def matrix(B):
    need(len(B)*int(np.max(abs(B)))**2<2**63,'int64 Gram entry budget')
    G=B.T@B
    return [[int(x) for x in row] for row in G]

def norm(G):
    return sum(abs(v) for row in G for v in row)

def energy(G,z):
    return sum(z[i]*G[i][j]*z[j] for i in range(len(z)) for j in range(len(z)))

def packing(G):
    n=len(G);C=[[abs(v) for v in row] for row in G];tri=[]
    for i,j,k in itertools.combinations(range(n),3):
        if G[i][j]*G[i][k]*G[j][k]<0:
            tri.append((min(C[i][j],C[i][k],C[j][k]),i,j,k))
    tri.sort(reverse=True);out=[]
    for _,i,j,k in tri:
        p=min(C[i][j],C[i][k],C[j][k])
        if p:
            out.append([i,j,k,p])
            for a,b in ((i,j),(i,k),(j,k)):
                C[a][b]-=p;C[b][a]-=p
    # Independent verification from the exported packing.
    used={}
    for i,j,k,p in out:
        need(p>0 and G[i][j]*G[i][k]*G[j][k]<0,'frustrated triangle')
        for e in ((i,j),(i,k),(j,k)):used[e]=used.get(e,0)+p
    need(all(v<=abs(G[i][j]) for (i,j),v in used.items()),'exact edge capacities')
    return out,norm(G)-4*sum(p for i,j,k,p in out),len(tri)

def controls():
    for G,opt in (([[2,1,1],[1,2,-1],[1,-1,2]],8),([[2,1,1],[1,2,1],[1,1,2]],12)):
        p,upper,_=packing(G)
        need(upper==opt,'balanced and frustrated triangle controls')
        need(max(energy(G,z) for z in itertools.product((-1,1),repeat=3))==upper,'cube vertices')
        need(all(energy(G,z)<=upper for z in itertools.product((-1,0,1),repeat=3)),'interior toy controls')

def main():
    controls()
    k=12;X=1<<k;b=1<<(k//2);N=8*X+1
    edges=observe(N,16);initial_edges=len(edges)
    base_cache,_=close_graph(edges,range(1,N))
    boundary=sorted({base_cache[n][0] for n in range(1,N)})
    for e in boundary:
        x=e
        for _ in range(16):
            if x==1:break
            edges[x]=T(x);x=T(x)
    cache,cycles=close_graph(edges,range(1,N));need(not cycles,'no nonROOT cycle')
    B,cols=rows(k,cache);keep=np.any(B!=0,axis=0);B=B[:,keep]
    cols=[e for e,f in zip(cols,keep) if f]
    need(cols==WITNESS_COLUMNS,'literal terminal identities')
    G=matrix(B);p,upper,tri_count=packing(G)
    z=WITNESS_SIGNS
    need(z[cols.index(1)]==1,'marked ROOT sign')
    lower=energy(G,z)
    # Independent row evaluation of the lower witness.
    need(lower==sum(sum(int(v)*s for v,s in zip(row,z))**2 for row in B),'row and Gram witness')
    need(lower==upper==9716780760,'matching exact cube bounds')
    budget=9*b*X*X
    need(lower>budget,'no box-only certificate can meet this budget')
    old_edges=len(edges)
    need(70924 not in edges and T(70924)==35462,'one genuinely new native edge')
    suffix=[];x=35462
    while x!=1:
        need(x in edges,'entire suffix already observed')
        suffix.append([x,edges[x]]);x=edges[x]
    need(len(suffix)==112,'old suffix integer clock')
    edges[70924]=35462
    cc,cycles=close_graph(edges,range(1,N))
    need(not cycles and cc[70924]==(1,-1,113),'new edge compiles the exact ROOT path')
    controller={e:cc[e][:2] for e in cols}
    tcols,H,slack=transport_slack(G,cols,controller)
    newB,newcols=rows(k,cc);keep=np.any(newB!=0,axis=0);newB=newB[:,keep]
    newcols=[e for e,f in zip(newcols,keep) if f];newG=matrix(newB)
    ix={e:i for i,e in enumerate(tcols)}
    need(all(newG[i][j]==H[ix[e]][ix[f]] for i,e in enumerate(newcols) for j,f in enumerate(newcols)), 'transport agrees with rebuilt fixed-source rows')
    newnorm=norm(newG)
    need(newnorm==norm(H)==9209359132 and newnorm<budget,'one-edge paid certificate')
    need(slack['old_norm']-slack['new_norm']==norm(G)-newnorm,'weighted payment equality')
    report=dict(status='FINITE-EXACT optimum for one fixed controller; one-edge repair; no uniform conclusion',
       k=k,source_count=N-1,initial_depth=16,boundary_depth=16,
       initial_observed_edges=initial_edges,boundary_added_edges=old_edges-initial_edges,
       columns=cols,gram=G,gram_sha256=hashlib.sha256(json.dumps(G,separators=(',',':')).encode()).hexdigest(),
       gram_norm=[norm(G),b],free_cube_optimum=[upper,b],budget=[budget,b],
       witness_signs=z,frustrated_triangles=tri_count,packing=p,
       new_edge=[70924,35462],old_suffix=suffix,new_native_clock=113,
       repaired_gram=[newnorm,b],weighted_payment=slack,
       unresolved_before=sum(cache[n][0]!=1 for n in range(1,N)),
       unresolved_after=sum(cc[n][0]!=1 for n in range(1,N)),
       max_native_clock_after=max(cc[n][2] for n in range(1,N)),
       discovery_scope='Candidate discovery scanned at most 16 native steps at each of 50 active nonROOT terminals; one edge is the exported additional certificate, not the total search cost.')
    path=ROOT/'05-knowledge/results/collatz_cube_context_20261015d.json'
    data=json.dumps(report,indent=2,sort_keys=True)+'\n';path.write_text(data)
    print(json.dumps(dict(output=str(path),sha256=hashlib.sha256(data.encode()).hexdigest(),
                         exact_optimum_over_X2=upper/(b*X*X),repaired_bound_over_X2=newnorm/(b*X*X),packed_triangles=len(p)),sort_keys=True))

if __name__=='__main__':main()
