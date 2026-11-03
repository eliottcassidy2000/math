#!/usr/bin/env python3
"""Modular reachable/observable ranks, with actual prefix/suffix witnesses.
The source carrier is the pinned prior active-port implementation.
"""
import importlib.util
from pathlib import Path
from collections import deque
import argparse,hashlib,json

SOURCE=Path(__file__).resolve().parents[1]/'adaptive_metagraph_reformulation_20261103'/'active_port_transfer.py'
spec=importlib.util.spec_from_file_location('active_source',SOURCE)
src=importlib.util.module_from_spec(spec);spec.loader.exec_module(src)


def build(d):
    T=src.birth((2,4,1),0);seed=src.degree_control(T,d)
    initial=src.literal_cover_vector(T,src.active(T,d))
    seen={(seed,s) for s in initial};todo=deque(seen);blocks={}
    while todo:
        c,s=todo.popleft()
        for a in range(len(c)):
            dest,terms=src.endpoint_step(s,c,a,d)
            for t,_ in terms:
                pair=(dest,t)
                if pair not in seen:seen.add(pair);todo.append(pair)
    for c,s in seen:blocks.setdefault(c,set()).add(s)
    blocks={c:tuple(sorted(states)) for c,states in sorted(blocks.items())}
    index={c:{s:i for i,s in enumerate(states)} for c,states in blocks.items()}
    transitions={};reverse={c:[] for c in blocks}
    for c,states in blocks.items():
        for a in range(len(c)):
            dest=src.control_step(c,a,d)[0]
            mat=tuple(tuple((index[dest][t],m) for t,m in src.endpoint_step(s,c,a,d)[1]) for s in states)
            transitions[c,a]=(dest,mat);reverse[dest].append((c,a,mat))
    u=[0]*len(blocks[seed])
    for s,w in initial.items():u[index[seed][s]]=w
    final={c:tuple(int(len(s)==1) for s in states) for c,states in blocks.items()}
    return seed,blocks,transitions,reverse,tuple(u),final


class Basis:
    def __init__(self,p):self.p=p;self.rows={};self.original=[];self.words=[]
    def add(self,v,word):
        w=list(v);p=self.p
        for pivot,row in sorted(self.rows.items()):
            a=w[pivot]
            if a:
                for j in range(pivot,len(w)):w[j]=(w[j]-a*row[j])%p
        pivot=next((i for i,x in enumerate(w) if x),None)
        if pivot is None:return False
        scale=pow(w[pivot],-1,p);w=[x*scale%p for x in w]
        self.rows[pivot]=w;self.original.append(tuple(v));self.words.append(tuple(word));return True


def row_apply(row,mat,n,p):
    out=[0]*n
    for i,x in enumerate(row):
        if x:
            for j,m in mat[i]:out[j]=(out[j]+x*m)%p
    return tuple(out)


def col_apply(col,mat,p):return tuple(sum(m*col[j] for j,m in terms)%p for terms in mat)


def probe(d,p):
    seed,blocks,trans,reverse,u,f=build(d)
    reach={c:Basis(p) for c in blocks};obs={c:Basis(p) for c in blocks}
    reach[seed].add(u,());todo=deque([(seed,u,())])
    while todo:
        c,row,w=todo.popleft()
        for a in range(len(c)):
            dest,mat=trans[c,a];v=row_apply(row,mat,len(blocks[dest]),p)
            if reach[dest].add(v,w+(a,)):todo.append((dest,v,w+(a,)))
    todo=deque()
    for c in blocks:
        obs[c].add(f[c],());todo.append((c,f[c],()))
    while todo:
        c,col,w=todo.popleft()
        for prev,a,mat in reverse[c]:
            v=col_apply(col,mat,p)
            if obs[prev].add(v,(a,)+w):todo.append((prev,v,(a,)+w))
    result=[]
    for c in blocks:
        hankel=Basis(p)
        for row,w in zip(reach[c].original,reach[c].words):
            vals=tuple(sum(x*y for x,y in zip(row,col))%p for col in obs[c].original)
            hankel.add(vals,w)
        result.append({'control':c,'carrier':len(blocks[c]),'reachable_mod_rank':len(reach[c].rows),
                       'observable_full_carrier_mod_rank':len(obs[c].rows),'local_hankel_mod_rank':len(hankel.rows),
                       'reachable_prefixes':reach[c].words,'observable_suffixes':obs[c].words})
    return {'degree_cap':d,'prime_modulus':p,'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
            'scope':'modular lower bounds only; no rational upper/minimality claim','blocks':result}


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--cap',type=int,default=2);ap.add_argument('--prime',type=int,default=1000003)
    args=ap.parse_args();print(json.dumps(probe(args.cap,args.prime),indent=2,sort_keys=True))
