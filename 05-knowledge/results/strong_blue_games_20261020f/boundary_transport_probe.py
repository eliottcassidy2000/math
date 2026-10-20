"""Small complete test of open-path rotation hubs, with literal replays."""
from collections import Counter, defaultdict, deque
from importlib.util import spec_from_file_location, module_from_spec
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
spec = spec_from_file_location('prior_outer', HERE.parent/'spectrum_basins_20261020e/outer_probe.py')
prior = module_from_spec(spec)
spec.loader.exec_module(prior)
G = prior.G


def paths(t):
    n=len(t)
    def visit(p, left):
        if not left:
            yield tuple(p)
            return
        opts=left if not p else left&t[p[-1]]
        while opts:
            bit=opts&-opts; opts-=bit
            yield from visit(p+[bit.bit_length()-1],left^bit)
    yield from visit([], (1<<n)-1)


def edge(u,v): return tuple(sorted((u,v)))


def canon(t):
    return min(G.induced(t,p) for p in paths(t))


def hub_graph(t, marked):
    ps=list(paths(t)); black={p for p in ps if not G.blue(t,p)}
    if not black:return None
    opened=[p for p in ps if t[p[0]]>>p[-1]&1]
    branches=defaultdict(list); incidence=defaultdict(list)
    for p in opened:
        partner=G.complement(t,p)
        for k in range(len(t)-1):
            rotation=p[k+1:]+p[:k+1]
            G.need(G.is_path(partner,rotation),'rotation legality')
            if G.blue(partner,rotation):continue
            key=tuple(sorted((edge(p[0],p[-1]), edge(p[k],p[k+1]))))
            target=list(t)
            for u,v in key:target[u]^=1<<v;target[v]^=1<<u
            G.need(G.complement(partner,rotation)==tuple(target),'hub identity')
            branches[key].append((p,rotation)); incidence[p].append(key)
    seen={marked}; todo=deque([marked]); pred={}; found=None
    while todo:
        p=todo.popleft()
        if p in black:found=p;break
        for hub in incidence[p]:
            for q,rotation in branches[hub]:
                if q not in seen:
                    seen.add(q);pred[q]=(p,hub);todo.append(q)
    witness=[]
    if found:
        at=found
        while at!=marked:
            prev,hub=pred[at];witness.append((prev,hub,at));at=prev
        witness.reverse()
    return dict(success=found is not None, reached=len(seen), paths=len(ps),
                black_paths=len(black), open_paths=len(opened),
                black_open=sum(p in black for p in opened),
                blue_component=list(sorted(seen)) if found is None else [],
                witness=witness)


def main(max_n=7):
    summary=[]; residual=[]
    for n in range(4,max_n+1):
        counts=Counter(); examples=[];p=tuple(range(n))
        for t in G.blue_tilings(n):
            counts['blue_presentations']+=1
            if not t[0]>>(n-1)&1:continue
            if not prior.strong(t) or not prior.strong(G.complement(t,p)):continue
            result=hub_graph(t,p)
            if result is None:continue
            partner=G.complement(t,p)
            if all(G.blue(partner,q) for q in paths(partner)):continue
            counts['open_strong_mixed_pairs']+=1
            ct,cp=canon(t),canon(partner)
            if ct==cp:
                counts['merged_loop']+=1
            elif result['success']:
                counts['hub_connected']+=1
                counts['length_'+str(2*len(result['witness'])+1)]+=1
            else:
                counts['residual']+=1
                residual.append(dict(n=n,rows=t,partner=partner,source_class=ct,
                                     partner_class=cp,degrees=[v.bit_count() for v in t],**result))
        summary.append(dict(n=n,counts=dict(counts)))
        print(n,dict(counts),flush=True)
    result=dict(status='FINITE-EXACT mechanism test, not global bridge test',
                universe='all blue fixed-path presentations orders4..'+str(max_n),
                summary=summary,residual=residual)
    (HERE/'boundary_transport_probe.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
