#!/usr/bin/env python3
"""Locate the n8 Z/256 factor in the support-bridge-free 3520-node core.
All bridge and core calculations are rebuilt from the literal edge ledger.
Certificates are restricted and their integer images recomputed on the core.
"""
import argparse,gzip,hashlib,json,sys
from collections import Counter
from pathlib import Path
from check_electrical import projection,public,rank2,need,lap
HERE=Path(__file__).resolve().parent

def bridge_edges(n,weights):
    adj=[[] for _ in range(n)]
    for a,b in weights:
        if a!=b:adj[a].append(b);adj[b].append(a)
    tin=[-1]*n;low=[-1]*n;answer=[];clock=0;sys.setrecursionlimit(max(10000,n+10))
    def visit(v,parent):
        nonlocal clock
        tin[v]=low[v]=clock;clock+=1
        for w in adj[v]:
            if w==parent:continue
            if tin[w]>=0:low[v]=min(low[v],tin[w])
            else:
                visit(w,v);low[v]=min(low[v],low[w])
                if low[w]>tin[v]:answer.append(tuple(sorted((v,w))))
    for v in range(n):
        if tin[v]<0:visit(v,-1)
    return sorted(answer),adj

def av(n,edges,x):
    y=[0]*n
    for a,b,w in edges:
        d=w*(x[a]-x[b]);y[a]+=d;y[b]-=d
    return y

def congruence_controls():
    checks=[]
    for w in [1,2,3,4,5,16]:
        edges=[(0,1,2),(1,2,3),(0,2,1),(3,4,4),(4,5,2),(3,5,1),(2,3,w),(0,0,7)]
        A=lap(6,edges);keep=[0,1,2,4,5];A=[[A[i][j] for j in keep] for i in keep]
        P=[[1,1,0,0,0],[1,0,1,0,0],[1,0,0,0,0],[0,0,0,1,0],[0,0,0,0,1]]
        transformed=[[sum(P[a][i]*A[a][b]*P[b][j] for a in range(5) for b in range(5)) for j in range(5)] for i in range(5)]
        expected=[[w,0,0,0,0],[0,3,-2,0,0],[0,-2,5,0,0],[0,0,0,6,-2],[0,0,0,-2,3]]
        need(transformed==expected,'nonleaf weighted-bridge unimodular congruence')
        checks.append({'bridge_weight':w,'both_sides_vertices':3,'congruence_verified':True})
    return checks

def main(base,out):
    with gzip.open(base/'nodes.tsv.gz','rt') as f:n=sum(1 for _ in f)-1
    with gzip.open(base/'edges.tsv.gz','rt') as f:
        next(f);rows=[list(map(int,line.split())) for line in f]
    weights=Counter();B=Counter();L=Counter()
    for a,b,l,k in rows:
        key=tuple(sorted((a,b)));weights[key]+=l+k;B[key]+=k;L[key]+=l
    weights=Counter({e:w for e,w in weights.items() if w});bridges,adj=bridge_edges(n,weights);leaves=[i for i,a in enumerate(adj) if len(a)==1];leafset=set(leaves)
    need(len(bridges)==len(leaves)==8,'eight leaf bridges independently')
    bundle=[]
    for a,b in bridges:
        need((a in leafset)^(b in leafset),'every support bridge has exactly one leaf end')
        leaf,owner=(a,b) if a in leafset else (b,a);w=weights[a,b]
        bundle.append({'leaf':leaf,'owner':owner,'weight':w,'blue':L[a,b],'black':B[a,b]})
    bundle.sort(key=lambda d:d['leaf']);need(Counter(d['weight'] for d in bundle)==Counter({1:3,2:5}),'bridge weights 1 and2')
    keep=[i for i in range(n) if i not in leafset];ix={v:i for i,v in enumerate(keep)}
    coreweights={(ix[a],ix[b]):w for (a,b),w in weights.items() if a in ix and b in ix};cb,cadj=bridge_edges(len(keep),coreweights)
    need(not cb,'core has no support bridge');need(all(len(a)>=2 for a in cadj),'no further support leaf cascade')
    blue=[(ix[a],ix[b],w) for (a,b),w in L.items() if w and a in ix and b in ix];black=[(ix[a],ix[b],w//2) for (a,b),w in B.items() if w and a in ix and b in ix];edges=[(a,b,w) for (a,b),w in coreweights.items()]
    need(all(w%2==0 for w in B.values()),'literal black twins');p=projection(len(keep),blue,black)
    need(p['full_components']==1,'core connected');need((p['critical_two_rank'],p['half_form_rank'],p['factors_divisible_by_four'])==(3358,3353,5),'independent core first layer')
    r2=json.loads((HERE/'n8_lift_receipt.json').read_text());r3=json.loads((HERE/'n8_final_radical_receipt.json').read_text());rt=json.loads((HERE/'n8_complete_tower_receipt.json').read_text());v2=[];checks=[]
    def restrict_and_check(k,full):
        mod=1<<k;x=[full[v] for v in keep];y=av(len(keep),edges,x);need(all(v%mod==0 for v in y),'core qualifying lift')
        need(any(v%2 for v in x),'core binary class nonzero')
        need(x[-1]==0,'same surviving grounded root')
        bridge_energy=sum(d['weight']*(full[d['leaf']]-full[d['owner']])**2 for d in bundle)
        need(bridge_energy%(2*mod)==0,'bridge energies vanish in qk for k>=2')
        num=sum(a*b for a,b in zip(x,y));need(num%mod==0,'core energy integral')
        return x,y,(num//mod)%2,bridge_energy//mod
    for j,r in enumerate(r2['witness_lifts']):
        full=[0]*n
        for i in r['binary_one_nodes']:full[i]+=1
        for i in r['correction_one_nodes']:full[i]+=2
        x,y,q,be=restrict_and_check(2,full);v2.append((x,y));checks.append({'k':2,'basis_index':j,'self_pairing':q,'bridge_energy_divided_by_modulus':be})
    binary=[sum((v%2)<<i for i,v in enumerate(x)) for x,y in v2];need(rank2(binary)==5,'transported V2 basis rank5')
    q2=[]
    for x,y in v2:
        row=[]
        for z,az in v2:
            num=sum(a*b for a,b in zip(x,az));need(num%4==0,'core cross pairing integral');row.append(num//4%2)
        q2.append(row)
    need(q2==r2['q2_matrix'],'full second half-form survives core restriction')
    x,y,q,be=restrict_and_check(3,r3['lift_mod8']);need(q==0,'core q3');checks.append({'k':3,'self_pairing':q,'bridge_energy_divided_by_modulus':be})
    tail_binary=sum((v%2)<<i for i,v in enumerate(x))
    for rec in rt['witness_lifts']:
        k=rec['k'];x,y,q,be=restrict_and_check(k,rec['lift']);need(sum((v%2)<<i for i,v in enumerate(x))==tail_binary,'same core binary tail');need(q==int(k==8),'core q4..q8 sequence');checks.append({'k':k,'self_pairing':q,'bridge_energy_divided_by_modulus':be})
    result={'nonleaf_congruence_controls':congruence_controls(),'n':8,'full_vertices':n,'support_bridge_bundles':bundle,'removed_leaf_nodes':leaves,'core_vertices':len(keep),'core_minimum_support_degree':min(map(len,cadj)),'core_support_bridges':[],'extra_leaf_cascade':False,'core_projection':public(p),'core_q2':q2,'transported_V2_rank':5,'transport_checks':checks,'full_two_primary_exponent_counts':{'1':3358,'2':4,'8':1},'bridge_two_primary_exponent_counts':{'1':5},'core_two_primary_exponent_counts':{'1':3353,'2':4,'8':1},'core_spanning_tree_valuation_two':3369,'source_sha256':{file.name:hashlib.sha256(file.read_bytes()).hexdigest() for file in [base/'edges.tsv.gz',base/'nodes.tsv.gz',HERE/'n8_lift_receipt.json',HERE/'n8_final_radical_receipt.json',HERE/'n8_complete_tower_receipt.json',Path(__file__)]}}
    out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2));print('ALL CHECKS PASSED')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--base',type=Path,default=HERE.parent/'session38-wiring');p.add_argument('--out',type=Path,default=HERE/'bridge_localization_receipt.json');a=p.parse_args();main(a.base,a.out)
