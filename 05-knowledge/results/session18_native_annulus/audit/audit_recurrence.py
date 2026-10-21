#!/usr/bin/env python3
"""Independent source-resolved recurrence audit; imports no producer code."""
import argparse,json,hashlib
from collections import defaultdict
from itertools import combinations
from pathlib import Path
import independent_audit as a

def charge(M,s):
    x=y=0
    while s<4*M:
        dx,dy=(a.L if s<2*M else a.U)[s%3];x+=dx;y+=dy;s=4*s+1
    return x,y

def changes(M,paths):
    delta=defaultdict(lambda:[0,0])
    # Evaluate the literal weight difference at the only three changed odd m.
    for m in (2*M+1,4*M+1,4*M+3):
        old=(0,0) if m>=4*M else (a.L if m<2*M else a.U)[m%3]
        new=(a.L if m<2*(M+1) else a.U)[m%3]
        s=a.ancestor(m);delta[s][0]+=new[0]-old[0];delta[s][1]+=new[1]-old[1]
    out=[]
    for s,(da,db) in sorted(delta.items()):
        va,vb=charge(M,s);na,nb=charge(M+1,s)
        a.require((va+da,vb+db)==(na,nb),'source-labelled endpoint identity')
        u,t,d=paths.get(s);out.append((s,u,t,d,va,vb,da,db))
    return out

def step(M,state,paths):
    edits=changes(M,paths);nV=len(paths.V);gs=[]
    for h in range(a.H):
        selected=[r for r in edits if r[3]<=h];D=defaultdict(lambda:[0,0]);g=0;masschange=0
        for s,u,t,d,va,vb,da,db in selected:
            sg=(-1)**t;k=h*nV+u
            g+=a.bet(state['a'][k]-sg*va,state['b'][k]-sg*vb,sg*da,sg*db)
            D[u][0]+=sg*da;D[u][1]+=sg*db
            masschange+=a.q(va+da,vb+db)-a.q(va,vb)
        for r,z in combinations(selected,2):
            if r[1]==z[1]:g+=a.bet((-1)**r[2]*r[6],(-1)**r[2]*r[7],(-1)**z[2]*z[6],(-1)**z[2]*z[7])
        grouped=sum(a.q(da,db)+a.bet(state['a'][h*nV+u],state['b'][h*nV+u],da,db) for u,(da,db) in D.items())-masschange
        a.require(g==grouped,'pairwise and grouped local coherence derivative')
        oldF=state['E'][h]-state['M'][h]
        for u,(da,db) in D.items():
            state['a'][h*nV+u]+=da;state['b'][h*nV+u]+=db
        state['E'][h]=sum(a.q(state['a'][h*nV+u],state['b'][h*nV+u]) for u in range(nV))
        state['M'][h]+=masschange
        a.require(g==state['E'][h]-state['M'][h]-oldF,'actual local coherence difference')
        gs.append(g)
    return gs,{'M':M,'edits':[[s,paths.V[u],t,d,va,vb,da,db] for s,u,t,d,va,vb,da,db in edits],'g_sha256':a.sh(a.enc(gs))}

def main():
    p=argparse.ArgumentParser();p.add_argument('--package',type=Path,required=True);p.add_argument('--out',type=Path,required=True);args=p.parse_args()
    fd=json.loads((args.package/'frontier.json').read_text());paths=a.Paths(fd['frontier'],dict(fd['frontier_ranks']),4*262208)
    low=a.direct(4096,paths);high=a.direct(262144,paths)
    G0=[e-m-f+n for e,m,f,n in zip(high['E'],high['M'],low['E'],low['M'])]
    g,record=step(4096,low,paths);steps=[record];total=[0]*a.H
    for M in range(262144,262208):
        z,record=step(M,high,paths);steps.append(record);total=[x+y for x,y in zip(total,z)]
    G1=[e-m-f+n for e,m,f,n in zip(high['E'],high['M'],low['E'],low['M'])]
    a.require([b-c for b,c in zip(G1,G0)]==[b-c for b,c in zip(total,g)],'exact 64-scale gap recurrence H')
    for M,state in ((4097,low),(262208,high)):
        final=a.direct(M,paths)
        for field in ('a','b','E','M'):a.require(final[field]==state[field],'independent complete endpoint universe '+str(M)+' '+field)
    producer=json.loads((args.package/'recurrence_receipt.json').read_text())
    for name,actual in [('G0',G0),('G1',G1),('old_g',g),('large_64step_g',total)]:a.require(actual==producer[name],'producer recurrence witness '+name)
    a.require((G0[-1],G1[-1],G1[-1]-G0[-1])==(-201674,-201616,58),'nonmonotonicity witness')
    # A simultaneous event does not require ordering its old and new arrivals.
    cases=0
    from itertools import product
    for W,H,X,Y in product(product((-1,0,1),repeat=2),repeat=4):
        before=a.q(W[0]+H[0],W[1]+H[1])-a.q(*W)
        after=a.q(W[0]+Y[0]+H[0]+X[0],W[1]+Y[1]+H[1]+X[1])-a.q(W[0]+Y[0],W[1]+Y[1])
        rhs=a.q(*X)+a.bet(W[0]+H[0],W[1]+H[1],*X)+a.bet(*Y,H[0]+X[0],H[1]+X[1])
        a.require(after-before==rhs,'simultaneous event identity C');cases+=1
    result={'status':'PASS; exact local identity, not a closed scalar recurrence or uniform estimate','steps':steps,'horizons':[0,256],'initial_source_scales':[4096,262144],'final_source_scales':[4097,262208],'final_G0_G1': [G0[-1],G1[-1]],'event_identity_cases':cases,'checks':a.CHECKS}
    args.out.write_bytes(a.enc(result));print(json.dumps({'status':'PASS','audit_sha256':a.sh(a.enc(result)),'steps':65,'horizons':257,'checks':a.CHECKS,'final_G0_G1':result['final_G0_G1']},sort_keys=True))
if __name__=='__main__':main()
