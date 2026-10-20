#!/usr/bin/env python3
"""Guarded exact exhaustive circulant certificates; never converts timeout to exclusion.
Run e.g. --max-order 25 --seconds 45. Default starts at17, checks ALL skew sets.
Certificates per order are gzip JSON; mtime=0 for deterministic serialization.
"""
import argparse,gzip,json,time
from pathlib import Path
from probe_cycles import circulant,bad_pair,verify_certificate
ROOT=Path(__file__).resolve().parent
class GuardExceeded(Exception):pass

def bounded_bad_cycle(A,deadline,node_cap,cycle_cap):
    n=len(A); C=[0]; nodes=0; cycles=0
    def dfs(mask):
        nonlocal nodes,cycles
        nodes+=1
        if nodes>node_cap:raise GuardExceeded('per-presentation node guard')
        if time.monotonic()>deadline:raise GuardExceeded('whole-run time guard')
        if len(C)==n:
            if not A[C[-1]][0]:return None
            cycles+=1
            if cycles>cycle_cap:raise GuardExceeded('cycle guard')
            p=[0]*n
            for j in range(n):p[C[j]]=C[(j+1)%n]
            pair=bad_pair(A,p)
            if pair is not None:
                for cut in range(n):
                    P=C[cut:]+C[:cut];rev=[0]*n
                    for j in range(n):rev[P[j]]=P[-j-1]
                    rp=bad_pair(A,rev,anti=True)
                    if rp is not None:return dict(cycle=C.copy(),successor_bad_pair=pair,
                         bad_cut_path=P,reversal_bad_pair=rp,cycles_checked=cycles)
                raise AssertionError('odd cycle reduction failed')
            return None
        for v in range(1,n):
            if not mask>>v&1 and A[C[-1]][v]:
                C.append(v);w=dfs(mask|1<<v);C.pop()
                if w is not None:return w
        return None
    w=dfs(1)
    return w,nodes,cycles

def main():
    p=argparse.ArgumentParser();p.add_argument('--min-order',type=int,default=17)
    p.add_argument('--max-order',type=int,default=25);p.add_argument('--seconds',type=float,default=45)
    p.add_argument('--node-cap',type=int,default=200000);p.add_argument('--cycle-cap',type=int,default=16)
    args=p.parse_args();start=time.monotonic();deadline=start+args.seconds
    for n in range(args.min_order,args.max_order+1,2):
        order_start=time.monotonic();rows=[];mxnodes=mxcycles=0
        total=1<<((n-1)//2)
        try:
            for mask in range(total):
                S,A=circulant(n,mask);w,nodes,cycles=bounded_bad_cycle(A,deadline,args.node_cap,args.cycle_cap)
                if w is None:raise GuardExceeded('all cycles survived: purity not decided')
                verify_certificate(n,S,w)
                rows.append(dict(mask=mask,certificate=w));mxnodes=max(mxnodes,nodes);mxcycles=max(mxcycles,cycles)
        except GuardExceeded as e:
            print(json.dumps(dict(n=n,status='INCOMPLETE',certified=len(rows),universe=total,reason=str(e))),flush=True)
            return
        payload=dict(n=n,universe='all skew subsets of nonzero Z/nZ',presentations=total,
                     excluded=total,max_nodes=mxnodes,max_cycles=mxcycles,rows=rows)
        raw=(json.dumps(payload,separators=(',',':'))+'\n').encode()
        (ROOT/f'circulant_n{n}_certificates.json.gz').write_bytes(gzip.compress(raw,mtime=0))
        print(json.dumps({k:payload[k] for k in ['n','presentations','excluded','max_nodes','max_cycles']}|
                         dict(seconds=round(time.monotonic()-order_start,3))),flush=True)
    print('ALL REQUESTED ORDERS COMPLETE',flush=True)
if __name__=='__main__':main()
