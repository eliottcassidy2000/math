#!/usr/bin/env python3
"""Universal exact high-forest compiler for every connected-low finite head."""
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import argparse, os

ROOT = Path(os.environ.get('LRC_REPO', Path(__file__).resolve().parents[1])).resolve()
MASS_PATH = ROOT/'04-computation/lrc_general_reflected_pair_mass_thm3352.py'
TAIL_PATH = ROOT/'04-computation/lrc14_connected_low_uniform_high_forest_tail_thm3350.py'
EXPECTED_MASS_SHA256 = 'afd417297131401254769e1ef172d89c109ad2f9a843ea55e2badc3e7891435b'
EXPECTED_TAIL_SHA256 = '78daaf73966d283c0c0bafa1c0975684e6167d2ef6375a3abeece4e00cdc87f9'
EXPECTED = {
    'contexts': 4044,
    'context_semantic': '9f7c7f24f81c409c09becfa921aa53387e5093710657bd3e5a0e935ecf4ea6c2',
    'heads': 261254,
    'channels': 4148,
    'scale_hist': ((1,204627),(2,45399),(3,8768),(4,1564),(5,518),(6,343),(7,14),(8,7),(9,7),(10,7)),
    'channel_semantic': '40472e196369baf0db300d96c028c2b7f78836bf8a06373b2170d4769ee9eff4',
    'forest_size_hist': ((2,48),(3,576),(4,6110),(5,254520)),
    'outcome_semantic': 'c84b550ace3ce3f0433227ad3587c92d674037a1b5f37531f27075742e038a9d',
    'weakest_margin': F(1017482276391320181594,73150827741949428345875),
}


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def load(name, path):
    s = spec_from_file_location(name, path)
    m = module_from_spec(s)
    s.loader.exec_module(m)
    return m


if sha256(TAIL_PATH.read_bytes().replace(b'\r\n',b'\n')).hexdigest()!=EXPECTED_TAIL_SHA256:
    raise RuntimeError(('tail hash',TAIL_PATH))
T = load('all_heads_tail', TAIL_PATH)
G = None
CONTEXTS = None


def init_worker(contexts):
    global G, CONTEXTS
    G = load('all_heads_mass', MASS_PATH)
    CONTEXTS = contexts


def channel_min(channel):
    scale,d,P,Q = channel
    best = None
    for L,cell,e,f in CONTEXTS:
        value = G.mass(L,cell,e,scale*d*P,f,scale*d*Q)
        record = (value,L,cell,e,f)
        if best is None or record < best:
            best = record
    return channel,best


def feasible_contexts():
    out=set()
    for body,L in T.SEL.MS.body_universe():
        cell,*_ = T.SEL.body_geometry(body,L)
        for e in body:
            for f in body:
                if e != f:
                    out.add((L,cell,e,f))
    return tuple(sorted(out))


def head_bank():
    heads=[];channels=set();scale_hist=Counter()
    for row in T.rows():
        if row in T.DENSE:
            continue
        S,_,_ = T.threshold(row)
        edges=T.high_edges(row)
        for scale in range(1,S):
            heads.append((row,scale,edges))
            scale_hist[scale]+=1
            channels.update((scale,d,P,Q) for _,_,d,P,Q in edges)
    return tuple(heads),tuple(sorted(channels)),scale_hist


def forest(row,scale,edges,mins):
    parent=list(range(6));chosen=[];total=F()
    def root(x):
        while parent[x]!=x:
            parent[x]=parent[parent[x]];x=parent[x]
        return x
    weighted=[]
    for i,j,d,P,Q in edges:
        value=mins[scale,d,P,Q]
        weighted.append((value,(i,j,d,P,Q)))
    for value,edge in sorted(weighted,reverse=True):
        a,b=root(edge[0]),root(edge[1])
        if a!=b:
            parent[a]=b;total+=value;chosen.append((value,edge))
    return total,tuple(chosen)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=10);ap.add_argument('--limit',type=int);a=ap.parse_args()
    require(sha256(MASS_PATH.read_bytes().replace(b'\r\n',b'\n')).hexdigest()==EXPECTED_MASS_SHA256,('mass hash',MASS_PATH))
    contexts=feasible_contexts();heads,channels,scale_hist=head_bank()
    require(len(contexts)==EXPECTED['contexts'],len(contexts));require(len(heads)==EXPECTED['heads'],len(heads))
    context_semantic=sha256()
    for context in contexts:context_semantic.update((repr(context)+'\n').encode())
    require(context_semantic.hexdigest()==EXPECTED['context_semantic'],context_semantic.hexdigest())
    require(tuple(sorted(scale_hist.items()))==EXPECTED['scale_hist'],scale_hist)
    if a.limit is not None:channels=channels[:a.limit]
    else:require(len(channels)==EXPECTED['channels'],len(channels))
    print('inventory contexts',len(contexts),'heads',len(heads),'channels',len(channels),'scale_hist',tuple(sorted(scale_hist.items())),flush=True)
    print('context_semantic',context_semantic.hexdigest(),flush=True)
    mins={};args={};sem=sha256()
    with ProcessPoolExecutor(max_workers=a.workers,initializer=init_worker,initargs=(contexts,)) as pool:
        for index,(channel,best) in enumerate(pool.map(channel_min,channels,chunksize=1),1):
            mins[channel]=best[0];args[channel]=best
            sem.update((repr((channel,best))+'\n').encode())
            if index%100==0:print('channel_progress',index,'of',len(channels),flush=True)
    if a.limit is not None:
        print('partial_channel_semantic',sem.hexdigest());return
    failures=[];weak=None;forest_hist=Counter();outsem=sha256()
    for row,scale,edges in heads:
        credit,chosen=forest(row,scale,edges,mins)
        margin=credit-T.DEBT_MAX/scale
        record=(margin,row,scale,credit,chosen)
        forest_hist[len(chosen)]+=1
        if weak is None or record<weak:weak=record
        if margin<=0:failures.append(record)
        outsem.update((repr(record)+'\n').encode())
    print('channel_semantic',sem.hexdigest())
    actual_forest_hist=tuple(sorted(forest_hist.items()))
    require(sem.hexdigest()==EXPECTED['channel_semantic'],sem.hexdigest())
    require(actual_forest_hist==EXPECTED['forest_size_hist'],actual_forest_hist)
    require(not failures,failures[:30]);require(weak[0]==EXPECTED['weakest_margin'],weak)
    require(outsem.hexdigest()==EXPECTED['outcome_semantic'],outsem.hexdigest())
    print('forest_size_hist',actual_forest_hist)
    print('failures',len(failures),'weakest',weak)
    print('failure_scale_hist',tuple(sorted(Counter(x[2] for x in failures).items())))
    print('first_failures',tuple(sorted(failures)[:30]))
    print('outcome_semantic',outsem.hexdigest())
    print('status=PROVED universal physical lower forest + VERIFIED-EXACT finite atlases')
    print('conclusion=all connected-low primitive common-dilation rays close at every scale on all 649 upper-median bodies and all labellings')


if __name__=='__main__':main()
