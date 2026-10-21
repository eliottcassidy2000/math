"""Exploration: synchronous prefix blocks, without looking beyond the cut."""
from collections import defaultdict,Counter
from importlib.util import spec_from_file_location,module_from_spec
from pathlib import Path
import argparse,json,sys,time

HERE=Path(__file__).resolve().parent
spec=spec_from_file_location('tail_native',HERE.parent/'session18_native_annulus/verify_annulus.py')
native=module_from_spec(spec);sys.modules[spec.name]=native;spec.loader.exec_module(native)

def Q(v):return v[0]*v[0]-v[0]*v[1]+v[1]*v[1]
def beta(v,w):return 2*v[0]*w[0]-v[0]*w[1]-v[1]*w[0]+2*v[1]*w[1]
def delta(v):return Q(v[2:])+beta(v[:2],v[2:])
def cross(v,w):return beta(v[2:],w[2:])+beta(v[:2],w[2:])+beta(w[:2],v[2:])
def add(a,b):return tuple(x+y for x,y in zip(a,b))

def bound(arrived,pending):
    current=sum(delta(v) for v in arrived.values())
    ps=list(Counter(pending.values()).items());aa=set(arrived.values())
    linear=sum(n*max(0,delta(v)+max([0]+[abs(cross(v,w)) for w in aa])) for v,n in ps)
    pairs=sum(n*(n-1)//2*abs(cross(v,v)) for v,n in ps)
    pairs+=sum(n*m*abs(cross(v,w)) for i,(v,n) in enumerate(ps) for w,m in ps[:i])
    return current,linear,pairs,current+linear+pairs

def run(N,D=160):
    start=time.perf_counter();V=set(json.loads((HERE.parent/'session18_native_annulus/frontier.json').read_text())['frontier'])
    arrived={};pending={};zeros=0;initial=0
    for s in range(1,256*N,2):
        if not native.source(s):continue
        q=native.charge(64*N,s)
        if q==(0,0):zeros+=1;continue
        v=q+(0,0) if s<4*N else (0,0)+q
        initial+=1
        (arrived if s in V else pending)[s]=v
    rows=[];cancelled=0;cache={}
    for d in range(D+1):
        if d%5==0 or not pending:
            b=bound(arrived,pending)
            row={'depth':d,'pending_blocks':len(pending),'arrived_ports':len(arrived),'bound':list(b),'budget':153*N,'cancelled_groups':cancelled}
            rows.append(row);print(json.dumps(row),flush=True)
            if b[-1]<=153*N:break
        if d==D:break
        nxt=defaultdict(lambda:(0,0,0,0))
        for u,v in pending.items():
            if u not in cache:cache[u]=native.edge(u)
            z,c=cache[u];sg=1-2*(c%2);w=tuple(sg*x for x in v)
            if z in V:arrived[z]=add(arrived.get(z,(0,0,0,0)),w)
            else:nxt[z]=add(nxt[z],w)
        pending={u:v for u,v in nxt.items() if any(v)};cancelled+=len(nxt)-len(pending)
    return {'N':N,'initial_nonzero':initial,'initial_zero':zeros,'rows':rows,'seconds':time.perf_counter()-start,'prefix_edges':len(cache)}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('N',type=int);p.add_argument('--depth',type=int,default=160);a=p.parse_args()
    result=run(a.N,a.depth);print(json.dumps({k:v for k,v in result.items() if k!='rows'}))
