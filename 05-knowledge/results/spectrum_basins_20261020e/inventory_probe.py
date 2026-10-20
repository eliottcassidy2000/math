"""Small literal discovery probes; no global connectivity inference."""
from pathlib import Path
import sys, json
from itertools import permutations
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'black_descent_20261020c'))
from scc_probe import paths,blue,reverse_path,word_tournament,need

def join(a,b):
    m=len(a); n=len(b)
    return tuple(x|(((1<<n)-1)<<m) for x in a)+tuple(x<<m for x in b)

def relabel(t,p):
    out=[0]*len(t)
    for u in range(len(t)):
        out[p[u]]=sum(1<<p[v] for v in range(len(t)) if t[u]>>v&1)
    return tuple(out)

def ball2(t,edge=None):
    out={}
    for p in paths(t):
        if blue(t,p) or edge is not None and edge not in {frozenset(e) for e in zip(p,p[1:])}:continue
        u=reverse_path(t,p)
        for q in paths(u):
            if blue(u,q) or edge is not None and edge not in {frozenset(e) for e in zip(q,q[1:])}:continue
            v=reverse_path(u,q)
            out.setdefault(v,(p,q))
    return out

def probe(n):
    source=word_tournament((1,)*(n-3)+(3,))
    target=word_tournament((1,)*(n-6)+(3,3))
    # Fix the prefix pointwise. Seek an even route with a common prefix edge.
    edge=frozenset((0,1)) if n>=8 else None
    forward=ball2(source,edge)
    for p in permutations(range(n-6,n)):
        p=tuple(range(n-6))+p
        t=relabel(target,p)
        back=ball2(t,edge)
        common=forward.keys()&back.keys()
        if common:
            v=min(common)
            route=forward[v]+tuple(tuple(reversed(q)) for q in reversed(back[v]))
            u=source
            for q in route:
                need(not blue(u,q),'blue route')
                u=reverse_path(u,q)
            need(u==t,'wrong target')
            return dict(n=n,source=source,target=t,permutation=p,paths=route)
    return dict(n=n,source=source,found=False)

if __name__=='__main__':
    records=[]
    for n in (7,8,9):
        row=probe(n); records.append(row); print(json.dumps(row),flush=True)
    HERE.joinpath('inventory_probe.json').write_text(json.dumps(records,indent=2)+'\n')
