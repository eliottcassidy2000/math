#!/usr/bin/env python3
"""Compile a native route in any conserved fibre of the fixed28-slot inventory."""
from common import *
from collections import defaultdict,deque
from itertools import combinations
import argparse
def compile_route(start,finish):
    R=inventory();ids=[r['source'] for r in R];n=len(ids);full=(1<<n)-1
    need(0<=start<=full and 0<=finish<=full,'states must be28-bit assignments')
    lab=[(sum(b<<j for j,b in enumerate(r['source_label'][:4])),r['source_label'][4]) for r in R]
    def signature(m):
        p=x=0
        for i,(q,y) in enumerate(lab):
            if m>>i&1:p^=q;x+=y
        return m.bit_count(),p,x%128
    need(signature(start)==signature(finish),'different integer margins or original dyadic currents')
    moves=set()
    for k in (2,3):
        bank=defaultdict(list)
        for I in combinations(range(n),k):
            m=sum(1<<i for i in I);bank[signature(m)].append(m)
        for bucket in bank.values():
            for a,b in combinations(bucket,2):
                A=a&~b;B=b&~a;moves.add(tuple(sorted((A,B))))
    extra=json.loads((HERE/'bridges.json').read_text());need(extra['sources']==ids,'bridge owners agree')
    moves.update(tuple(x) for x in extra['moves']);moves=sorted(moves)
    need(len(moves)==2164,'complete sparse alphabet')
    todo=deque([start]);parent={start:None}
    while finish not in parent and todo:
        m=todo.popleft()
        for A,B in moves:
            if m&A==A and not m&B:nxt=(m^A)|B
            elif m&B==B and not m&A:nxt=(m^B)|A
            else:continue
            if nxt not in parent:parent[nxt]=m;todo.append(nxt)
    need(finish in parent,'audited complete-fibre compiler unexpectedly failed')
    states=[];m=finish
    while m is not None:states.append(m);m=parent[m]
    states.reverse();events=[]
    for a,b in zip(states,states[1:]):
        changed=[i for i in range(28) if (a^b)>>i&1];cost=sum(R[i]['cost'] for i in changed)
        need(len(changed)<=8 and cost<=22,'literal event support and clock bounds')
        need(signature(a)==signature(b)==signature(start),'complete conserved source signature at every arrival')
        events.append({'from_mask':a,'to_mask':b,'slots':len(changed),'owner_clocks':cost,
                       'changes':[{'source':ids[i],'pins':R[i]['pins'],
                                   'from_path':R[i]['after' if a>>i&1 else 'before']['path'],
                                   'to_path':R[i]['after' if b>>i&1 else 'before']['path']} for i in changed]})
    need(len(events)<=10,'root-tree upper bound controls the shortest-event route')
    return {'status':'FINITE-EXACT route in original complete28-source supply','sources':ids,'signature':signature(start),
            'start_V':[s for i,s in enumerate(ids) if start>>i&1],'finish_V':[s for i,s in enumerate(ids) if finish>>i&1],
            'states':states,'events':events,'total_owner_clocks':sum(e['owner_clocks'] for e in events),
            'convention':'fixed source/path/pins, original paired multiplicity2, unit capacity with return; disjoint background required',
            'input_revision':REV,'inputs':INPUTS}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--start',type=int,required=True);p.add_argument('--finish',type=int,required=True);a=p.parse_args()
    print(json.dumps(compile_route(a.start,a.finish),indent=2,sort_keys=True))
