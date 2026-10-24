#!/usr/bin/env python3
"""Actual chart connectors and exact uniform-marginal transportation tests."""
from collections import defaultdict,deque,Counter
from itertools import combinations,permutations,product
from pathlib import Path
import json,subprocess

CHECKS=0
def need(ok,label):
    global CHECKS
    CHECKS+=1
    if not ok:raise RuntimeError(label)

def canonical(B):
    groups=defaultdict(list)
    for v in range(8):groups[sum(B[v])].append(v)
    best=None
    for pieces in product(*(permutations(groups[d]) for d in sorted(groups))):
        p=sum(pieces,());code=0
        for i,j in combinations(range(8),2):code=(code<<1)|B[p[i]][p[j]]
        if best is None or code<best:best=code
    return best

def transport(left,right,edges):
    # Uniform marginals scaled by |L|*|R|. Integral max flow is exact.
    n=len(left);m=len(right);S=('s',);T=('t',);cap=defaultdict(int);adj=defaultdict(set)
    def edge(u,v,c):cap[u,v]+=c;adj[u].add(v);adj[v].add(u)
    for i in left:edge(S,('l',i),m)
    for j in right:edge(('r',j),T,n)
    for i,j in edges:edge(('l',i),('r',j),n*m)
    flow=0
    while True:
        prev={S:None};todo=deque([S])
        while todo and T not in prev:
            u=todo.popleft()
            for v in sorted(adj[u],key=repr):
                if v not in prev and cap[u,v]>0:prev[v]=u;todo.append(v)
        if T not in prev:break
        amount=n*m;v=T
        while prev[v] is not None:amount=min(amount,cap[prev[v],v]);v=prev[v]
        v=T
        while prev[v] is not None:
            u=prev[v];cap[u,v]-=amount;cap[v,u]+=amount;v=u
        flow+=amount
    obstruction=None
    for bits in range(1,1<<n):
        A={left[k] for k in range(n) if bits>>k&1};N={j for i,j in edges if i in A}
        if len(A)*m>len(N)*n:
            obstruction=dict(left=sorted(A),neighbors=sorted(N),left_numerator=len(A)*m,right_numerator=len(N)*n,denominator=n*m);break
    need((flow==n*m)==(obstruction is None),'weighted Hall subsets agree with exact flow')
    witness=[dict(left=i,right=j,numerator=cap[('r',j),('l',i)],denominator=n*m)
             for i,j in edges if cap[('r',j),('l',i)]] if flow==n*m else []
    return flow==n*m,witness,obstruction

def main():
    root=Path(__file__).resolve().parents[3]
    raw=subprocess.check_output(['git','show','0f9dec9c3e1d:05-knowledge/results/guarded_coverage_20261023/atlas_receipt.json'],cwd=root)
    atlas=json.loads(raw)['charts_data']
    A=[[0]*8 for _ in range(8)]
    for k,(i,j) in enumerate(combinations(range(8),2)):
        A[i][j]=(29558785>>(27-k))&1;A[j][i]=1-A[i][j]
    rows=[];fibres=[]
    for c,chart in enumerate(atlas):
        fibre=[]
        for P in chart['paths']:
            p=tuple(P);pos={v:i for i,v in enumerate(p)}
            K=[[0 if i==j else A[i][j]^(abs(pos[i]-pos[j])!=1) for j in range(8)] for i in range(8)]
            fibre.append(len(rows));rows.append(dict(path=p,chart=c,target=canonical(K)))
        fibres.append(fibre)
    coarse=defaultdict(list);edges=defaultdict(list)
    for i,row in enumerate(rows):coarse[row['path'][0],row['path'][-1],row['target']].append(i)
    for ids in coarse.values():
        for i,j in combinations(ids,2):
            c,d=rows[i]['chart'],rows[j]['chart']
            if c!=d:
                if c>d:i,j,c,d=j,i,d,c
                edges[c,d].append((i,j))
    result=[]
    for (c,d),E in sorted(edges.items()):
        stages=[]
        for width in range(2,7):
            filtered=[]
            for i,j in E:
                support=[k for k in range(8) if rows[i]['path'][k]!=rows[j]['path'][k]]
                if max(support)-min(support)+1<=width:filtered.append((i,j))
            feasible,witness,cut=transport(fibres[c],fibres[d],filtered)
            stages.append(dict(window=width,edges=filtered,balanced=feasible,flow=witness,cut=cut))
        result.append(dict(charts=[c,d],sizes=[len(fibres[c]),len(fibres[d])],edges=E,
                           stages=stages,minimum_balanced_window=next(s['window'] for s in stages if s['balanced'])))
    return dict(status='FINITE-EXACT actual source1238; endpoints and ordinary target class fixed',
                paths=len(rows),pointed_charts=len(fibres),coarse_fibres=len(coarse),
                coarse_chart_multiplicity=dict(Counter(len({rows[i]['chart'] for i in ids}) for ids in coarse.values())),
                cross_chart_pairs=len(result),balanced_pairs=len(result),
                connectors=result,rows=rows,checks=CHECKS)

if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
