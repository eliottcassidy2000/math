"""Endpoint-orbit obstruction to the generalized terminal-ear detour."""
from collections import Counter
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
spec = spec_from_file_location('outer_old', HERE.parent / 'scc_basins_20261020d/moves.py')
G = module_from_spec(spec)
spec.loader.exec_module(G)


def fixed_paths(n):
    free = [(u,v) for u,v in combinations(range(n),2) if v != u+1]
    base = [1 << (u+1) if u+1<n else 0 for u in range(n)]
    for mask in range(1 << len(free)):
        a = base[:]
        for bit,(u,v) in enumerate(free):
            x,y=(u,v) if mask>>bit&1 else (v,u)
            a[x] |= 1<<y
        yield tuple(a)


def paths_to(a, z):
    n=len(a)
    reverse=[sum(1<<u for u in range(n) if a[u]>>v&1) for v in range(n)]
    def walk(rev, unused):
        if not unused:
            yield tuple(reversed(rev))
        else:
            choices=reverse[rev[-1]]&unused
            while choices:
                bit=choices&-choices
                choices-=bit
                u=bit.bit_length()-1
                yield from walk(rev+[u], unused^bit)
    yield from walk([z], ((1<<n)-1)^(1<<z))


def strong(a):
    reached=1
    while True:
        new=reached
        for u in range(len(a)):
            if reached>>u&1:
                new |= a[u]
        if new==reached:
            break
        reached=new
    if reached!=(1<<len(a))-1:
        return False
    rev=tuple(sum(1<<v for v in range(len(a)) if a[v]>>u&1) for u in range(len(a)))
    reached=1
    while True:
        new=reached
        for u in range(len(a)):
            if reached>>u&1:
                new |= rev[u]
        if new==reached:
            return reached==(1<<len(a))-1
        reached=new


def endpoint_scan(max_n=7):
    out=[]
    for n in range(4,max_n+1):
        counts=Counter()
        rigid=[]
        for a in fixed_paths(n):
            if not strong(a):
                continue
            counts['strong_presentations']+=1
            different=next((p for p in paths_to(a,n-1) if G.induced(a,p)!=a), None)
            if different is not None:
                counts['endpoint_flexible']+=1
                continue
            counts['endpoint_rigid']+=1
            G.need(a[-1]&1,'rigid terminal path does not close')
            rotation=tuple(range(1,n))+(0,)
            G.need(G.induced(a,rotation)!=a, 'rigid strong nontriangle rotates automorphically')
            chord=bool(a[-1]&2)
            counts['ear_good' if chord else 'ear_bad']+=1
            counts['rigid_degree_'+str(a[-1].bit_count())]+=1
            if not chord:
                rigid.append(dict(rows=a, degrees=[x.bit_count() for x in a], hp_to_end=sum(1 for _ in paths_to(a,n-1))))
        out.append(dict(n=n,counts=dict(counts),hostiles=rigid))
        print(n,dict(counts),flush=True)
    (HERE/'outer_probe.json').write_text(json.dumps(out,indent=2)+'\n')
    return out


if __name__=='__main__':
    endpoint_scan()
