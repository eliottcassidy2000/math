"""Hostile probe of moving positive-part collateral; no route beyond the cut."""
from collections import defaultdict
from importlib.util import spec_from_file_location, module_from_spec
from pathlib import Path
import argparse
import json
import sys

HERE = Path(__file__).resolve().parent
spec = spec_from_file_location('prefix_native', HERE.parent/'session18_native_annulus/verify_annulus.py')
native = module_from_spec(spec)
sys.modules[spec.name] = native
spec.loader.exec_module(native)


def delta(v):
    return native.Q(v[2:])+native.beta(v[:2],v[2:])


def add(a,b):
    return tuple(x+y for x,y in zip(a,b))


def run(N,D):
    V=set(json.loads((HERE.parent/'session18_native_annulus/frontier.json').read_text())['frontier'])
    arrived={}; pending={}; cache={}; rows=[]; witnesses=[]
    for s in range(1,256*N,2):
        if not native.source(s): continue
        q=native.charge(64*N,s)
        if not any(q): continue
        v=q+(0,0) if s<4*N else (0,0)+q
        (arrived if s in V else pending)[s]=v
    for d in range(D+1):
        e=sum(delta(v) for v in arrived.values())
        moving=sum(max(0,delta(v)) for v in pending.values())
        rows.append([d,e,moving,e+moving,len(pending)])
        if d==D: break
        nxt=defaultdict(lambda:(0,0,0,0)); buckets=defaultdict(list)
        for u,v in pending.items():
            if u not in cache: cache[u]=native.edge(u)
            z,c=cache[u]; sg=1-2*(c%2); w=tuple(sg*x for x in v)
            nxt[z]=add(nxt[z],w); buckets[z].append([u,c,list(v)])
        growth=0; maximum=None
        for z,v in nxt.items():
            before=sum(max(0,delta(tuple(row[2]))) for row in buckets[z])
            base=arrived.get(z,(0,0,0,0))
            change=(delta(add(base,v))-delta(base) if z in V else max(0,delta(v)))-before
            growth+=change
            if change>0 and (maximum is None or change>maximum['increase']):
                maximum={'depth':d+1,'target':z,'arrived':z in V,'base':base,'children':buckets[z],'increase':change}
            if z in V: arrived[z]=add(base,v)
        if growth>0 and len(witnesses)<3:
            witnesses.append({'depth':d+1,'global_increase':growth,'largest_local':maximum})
        pending={z:v for z,v in nxt.items() if z not in V and any(v)}
    return {'N':N,'D':D,'budget':153*N,'rows':rows,'positive_growth_witnesses':witnesses,
            'maximum_positive_part_potential':max(row[3] for row in rows)}


if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('N',type=int); p.add_argument('--depth',type=int,default=150)
    a=p.parse_args(); result=run(a.N,a.depth)
    (HERE/f'probe_{a.N}.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='rows'}))
