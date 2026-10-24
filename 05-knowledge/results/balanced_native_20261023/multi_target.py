#!/usr/bin/env python3
"""Directed-cycle then zero-sum packet compiler, with physical target rows retained."""
from collections import Counter
from pathlib import Path
from random import Random
import json,subprocess
CHECKS=0
def need(x,m):
    global CHECKS
    CHECKS+=1
    if not x:raise RuntimeError(m)
def compile_edges(edges,t,mods):
    unused=set(range(len(edges)));cycles=[]
    while unused:
        v=edges[min(unused)][0];seen={};walk=[]
        while v not in seen:
            seen[v]=len(walk)
            candidates=[i for i in unused if edges[i][0]==v]
            need(bool(candidates),'balanced remainder supplies outgoing edge')
            i=min(candidates);walk.append(i);unused.remove(i);v=edges[i][1]
        cut=seen[v];unused.update(walk[:cut]);C=walk[cut:]
        need(1<=len(C)<=t,'simple target cycle bound');cycles.append(C)
    need(sorted(i for C in cycles for i in C)==list(range(len(edges))),'actual slot partition')
    def charge(ids):
        z=[[0]*len(mods) for _ in range(t)]
        for i in ids:
            u,v,g=edges[i]
            for j,m in enumerate(mods):z[v][j]=(z[v][j]+g[j])%m;z[u][j]=(z[u][j]-g[j])%m
        need(all(sum(row[j] for row in z)%m==0 for j,m in enumerate(mods)),'full physical row sum retained')
        return tuple(x for row in z[:-1] for x in row)
    M=mods*(t-1);zero=(0,)*len(M);D=1+sum(m-1 for m in M)
    add=lambda a,b:tuple((x+y)%m for x,y,m in zip(a,b,M))
    neg=lambda a:tuple(-x%m for x,m in zip(a,M))
    need(charge(range(len(edges)))==zero,'accepted endpoints have equal full target charge')
    todo=list(range(len(cycles)));packets=[]
    while todo:
        bank={};chosen=None
        for c in todo[:D]:
            g=charge(cycles[c])
            if g==zero:chosen=(c,);break
            if neg(g) in bank:chosen=bank[neg(g)]+(c,);break
            updates={add(h,g):w+(c,) for h,w in bank.items()};updates.setdefault(g,(c,))
            for h,w in updates.items():bank.setdefault(h,w)
        need(chosen is not None,'p-group zero-sum extraction')
        slots=[i for c in chosen for i in cycles[c]]
        need(len(slots)<=t*D and charge(slots)==zero,'bounded full-row-neutral packet')
        packets.append(slots);todo=[c for c in todo if c not in chosen]
    # Recompute the complete physical target data after each atomic event.
    current=[u for u,v,g in edges];counts=Counter(current)
    def current_rows(state):
        z=[[0]*len(mods) for _ in range(t)]
        for v,(_,_,g) in zip(state,edges):
            for j,m in enumerate(mods):z[v][j]=(z[v][j]+g[j])%m
        return z
    baseline=current_rows(current);used=set()
    for packet in packets:
        need(not used&set(packet),'no repeated changed slot');used.update(packet)
        for i in packet:current[i]=edges[i][1]
        need(Counter(current)==counts and current_rows(current)==baseline,'every actual count and target row preserved')
    need(current==[v for u,v,g in edges],'specified endpoint reached')
    return {'slots':len(edges),'cycles':len(cycles),'events':len(packets),'max_support':max(map(len,packets)),'bound':t*D}
def main():
    rng=Random(789264);reports=[]
    for mods in [(2,4),(3,)]:
        for t in range(2,6):
            for trial in range(4):
                edges=[]
                for _ in range(8):
                    vertices=rng.sample(range(t),rng.randrange(2,t+1))
                    cycle=[]
                    for u,v in zip(vertices,vertices[1:]+vertices[:1]):
                        g=tuple(rng.randrange(m) for m in mods);cycle.append((u,v,g))
                    edges.extend(cycle);edges.extend((v,u,g) for u,v,g in cycle)
                rng.shuffle(edges);reports.append({'group':mods,'targets':t,**compile_edges(edges,t,mods)})
    root=Path(__file__).resolve().parents[3]
    def load(rev,p):return json.loads(subprocess.check_output(['git','show',rev+':05-knowledge/results/'+p],cwd=root))
    d=load('0f9dec9c3e1d','metagraph_order8_electrical_20261023/receipt.json')['order8'];X=d['higher_layers'][-1]['lifts'][0];Y=d['radical_lifts_mod4']
    lab=lambda i:tuple(y[i]%2 for y in Y[:4])+(X[i]%128,)
    fibre=[i for i in range(len(X)) if lab(i)==lab(1965)];need(fibre==[570,1965,2234],'complete literal three-target charge fibre')
    domains=load('66db753a59','session46_native_packet_decomposition_20261024/pinned-domains/normal.json')['domains']
    paths=[p for r in domains for p in r['complete_pinned_paths']]
    need(len(paths)==57 and not any(p['target']==570 for p in paths),'third equal-charge target lacks a pinned native option in this full supply')
    return {'status':'PASS conditional compiler controls, not abstract-label native realizations','cases':reports,
      'original_target_fibre':fibre,'current_pinned_paths':len(paths),'paths_to_third_target':0,'checks':CHECKS}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
