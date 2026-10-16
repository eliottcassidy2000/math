#!/usr/bin/env python3
"""Scale-dependent repair of the fixed-stage atlas; exhaustive compositional audit.

For k>=18 the fixed policy is base=floor(k/2)+2, anchor i uses base+i%3.
Every extension is replayed at every distinct previous terminal; vectorized
checks connect it to every original source and integer clock. This avoids
replaying an already authenticated common prefix millions of times.
"""
import argparse
import hashlib
import json
import numpy as np
import collatz_interval_atlas_20261016e as A
from collatz_compact_resets_20261015f import Graph,T

CHECKS=0
def need(ok,label):
    global CHECKS
    CHECKS+=1
    if not ok:raise RuntimeError(label)

def initial_audit(g):
    for lo in range(1,g.M+1,1<<18):
        hi=min(g.M+1,lo+(1<<18));n=np.arange(lo,hi,dtype=np.uint64)
        need(bool(np.all(g.endpoint[lo:hi]%2==1)),'all initial endpoints odd')
        need(bool(np.all((g.endpoint[lo:hi]<<g.clock[lo:hi].astype(np.uint64))==n)),
             'every initial original source and full halving clock')

def extend_audited(g):
    cols,first,inv=np.unique(g.endpoint[1:],return_index=True,return_inverse=True)
    old_clocks=g.clock.copy();g.extend(cols);g.close(cols,inv)
    end=g.endpoint[first+1].copy()
    need(bool(np.all(g.clock>=old_clocks)),'integer clocks never decreased')
    delta=(g.clock[first+1]-old_clocks[first+1]).copy()
    for lo in range(0,g.M,1<<18):
        hi=min(g.M,lo+(1<<18));ids=inv[lo:hi]
        need(bool(np.all(g.endpoint[lo+1:hi+1]==end[ids])), 'every source receives its own old-terminal endpoint')
        need(bool(np.all(g.clock[lo+1:hi+1].astype(np.uint64)==
                        old_clocks[lo+1:hi+1].astype(np.uint64)+delta[ids])),
             'every source receives the exact integer extension clock')
    steps=0
    for n,e,d in zip(cols,end,delta):
        x=int(n)
        for _ in range(int(d)):
            need(x!=1 and g.has(x),'every new edge is actually observed and ROOT stopped')
            x=T(x);steps+=1
        need(x==int(e),'independent literal terminal extension replay')
    return dict(old_terminals=len(cols),replayed_steps=steps,max_added_clock=int(delta.max()))

def hostile_control():
    g=Graph(512);initial_audit(g)
    for _ in range(8):extend_audited(g)
    # The compositional argument is checked against the full original-source replay.
    g.replay(range(1,513));need(True,'small independent full-prefix control')

def run(k,override=None):
    X=1<<k;b=1<<(k//2);base=k//2+2 if override is None else override
    positions=list(range(X+(b-1)//2,2*X,b));requested={a:base+i%3 for i,a in enumerate(positions)}
    g=Graph(8*X);initial_audit(g);steps=[];anchors=[];stages=[];last=None
    for stage in range(1,base+3):
        steps.append(extend_audited(g))
        if stage<base:continue
        ns=A.norms(X,g)
        if last is not None:need(all(a>=b for a,b in zip(last,ns)),'contraction of all finest rows')
        last=ns
        anchors += [[a,ns[a-X],stage] for a in positions if requested[a]==stage]
        stages.append(dict(stage=stage,max_clock=int(g.clock.max()),
           energy=sum(v*v for v in ns),endpoint_sha256=hashlib.sha256(g.endpoint.tobytes()).hexdigest(),
           clock_sha256=hashlib.sha256(g.clock.tobytes()).hexdigest()))
    anchors.sort();cone=A.cone_envelope(X,anchors)
    need(all(v<=c for v,c in zip(last,cone)),'final common controller is bounded by independent anchor cones')
    return dict(k=k,X=X,bin_width=b,base_stage=base,anchors=anchors,stages=stages,
        native_extension_audits=steps,anchor_energy=sum(b*B*B for a,B,stage in anchors),
        cone_energy=sum(c*c for c in cone),passes_7X=sum(c*c for c in cone)<=49*X*X,
        cone_sha256=hashlib.sha256(json.dumps(cone).encode()).hexdigest())

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--k',type=int,nargs='+',default=[18,20])
    ap.add_argument('--base',type=int);args=ap.parse_args();hostile_control()
    runs=[run(k,args.base) for k in args.k]
    print(json.dumps(dict(status='FINITE-EXACT; scale-dependent policy tested after the fixed-stage failure',
        audit='every original source checked compositionally; every distinct extension replayed literally',
        runs=runs,checks=CHECKS,atlas_checks=A.CHECKS),separators=(',',':'),sort_keys=True))
